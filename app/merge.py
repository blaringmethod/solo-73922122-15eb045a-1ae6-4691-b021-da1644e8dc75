"""Amplicon read merging for succeeded pairing tasks.

Given the published R1/R2 FASTQ files of a *succeeded* task, this module
stitches each paired-end read into a single sequence for downstream variant
analysis:

* R2 is reverse-complemented (its quality string is merely reversed);
* an overlap between the R1 suffix and the rc(R2) prefix is searched for;
* a candidate overlap needs >= ``MERGE_MIN_OVERLAP`` aligned bases and
  <= ``MERGE_MAX_MISMATCHES`` mismatching columns;
* only the *longest* candidate is used; when several candidates share the
  longest length the pair is skipped (no arbitrary tie break);
* pairs without any admissible candidate are skipped as well.

Everything is streamed one pair at a time with bounded per-line reads, so a
maximum-size request (5000 pairs of up to 200 bases) keeps a tiny working set
and stays well inside the 256 MiB process memory cap.

Validation is all-or-nothing: any structural or content defect raises
``MergeValidationError`` before the response is built, so a caller can never
observe a partial merge output. The published R1/R2 files are opened read-only
and never modified.
"""

from __future__ import annotations

from pathlib import Path

from . import config
from .fastq import extract_read_id

# bytes.translate table: complement of the ACGT alphabet; anything else maps
# to N so reverse-complement is total even before validation rejects it.
_COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


class MergeValidationError(ValueError):
    """Raised when a merge input violates the merge contract.

    The whole request fails with this error; no partial output is produced.
    """


def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(_COMPLEMENT)[::-1]


def _read_bounded_line(fh, side: str, record_no: int, what: str) -> bytes:
    """Read one logical line, refusing pathological unbounded lines."""
    line = fh.readline(config.MERGE_MAX_LINE_BYTES + 1)
    if line and len(line) > config.MERGE_MAX_LINE_BYTES and not line.endswith(b"\n"):
        raise MergeValidationError(
            f"{side}: record {record_no} {what} line exceeds "
            f"{config.MERGE_MAX_LINE_BYTES} bytes"
        )
    return line


def _read_record(fh, side: str, record_no: int) -> tuple[str, bytes, bytes] | None:
    """Parse one strict four-line FASTQ record from ``fh``.

    Returns ``(normalized_read_id, sequence, quality)`` with line endings
    stripped, or ``None`` at clean EOF. Every structural defect raises
    ``MergeValidationError`` (truncated records, bad marker lines, lengths).
    """
    header = _read_bounded_line(fh, side, record_no, "header")
    if not header:
        return None
    sequence = _read_bounded_line(fh, side, record_no, "sequence")
    plus = _read_bounded_line(fh, side, record_no, "plus")
    quality = _read_bounded_line(fh, side, record_no, "quality")
    if not sequence or not plus or not quality:
        raise MergeValidationError(
            f"{side}: record {record_no} is truncated "
            "(FASTQ records must span exactly four lines)"
        )
    header = header.rstrip(b"\r\n")
    if not header.startswith(b"@"):
        raise MergeValidationError(
            f"{side}: record {record_no} header does not start with '@'"
        )
    try:
        read_id = extract_read_id(header)
    except ValueError as exc:
        raise MergeValidationError(f"{side}: record {record_no}: {exc}") from exc
    if not plus.rstrip(b"\r\n").startswith(b"+"):
        raise MergeValidationError(
            f"{side}: record {record_no} third line does not start with '+'"
        )
    seq = sequence.rstrip(b"\r\n")
    qual = quality.rstrip(b"\r\n")
    if not seq or not qual:
        raise MergeValidationError(
            f"{side}: record {record_no} has an empty sequence or quality string"
        )
    if len(seq) != len(qual):
        raise MergeValidationError(
            f"{side}: record {record_no} sequence length {len(seq)} "
            f"!= quality length {len(qual)}"
        )
    return read_id, seq, qual


def _validate_content(
    side: str, record_no: int, seq: bytes, qual: bytes
) -> None:
    """Enforce the merge-specific alphabet, length and Phred+33 ranges."""
    n = len(seq)
    if not config.MERGE_MIN_READ_LEN <= n <= config.MERGE_MAX_READ_LEN:
        raise MergeValidationError(
            f"{side}: record {record_no} read length {n} outside "
            f"[{config.MERGE_MIN_READ_LEN}, {config.MERGE_MAX_READ_LEN}]"
        )
    # Strict ACGT alphabet: every base byte must be an uppercase ACGT column.
    bad = next((b for b in seq if b not in (0x41, 0x43, 0x47, 0x54)), None)
    if bad is not None:
        raise MergeValidationError(
            f"{side}: record {record_no} sequence contains non-ACGT base "
            f"{bytes([bad])!r}; only uppercase A/C/G/T are accepted"
        )
    qbad = next(
        (b for b in qual if not config.MERGE_QUAL_MIN <= b <= config.MERGE_QUAL_MAX),
        None,
    )
    if qbad is not None:
        raise MergeValidationError(
            f"{side}: record {record_no} quality byte {qbad} outside the "
            f"Phred+33 range [{config.MERGE_QUAL_MIN}, {config.MERGE_QUAL_MAX}]"
        )


def _iter_paired_records(
    r1_path: Path, r2_path: Path
):
    """Yield validated ``(read_id, r1_seq, r1_qual, rc_seq, rev_qual)`` tuples.

    Both files are read strictly in lockstep, which verifies the FASTQ
    four-line structure and the pairing of read ids. The 5000-pair cap is
    enforced while streaming, so an oversized input never accumulates in RAM.
    """
    with r1_path.open("rb") as f1, r2_path.open("rb") as f2:
        record_no = 0
        while True:
            r1 = _read_record(f1, "R1", record_no + 1)
            r2 = _read_record(f2, "R2", record_no + 1)
            if r1 is None or r2 is None:
                if r1 is not None or r2 is not None:
                    raise MergeValidationError(
                        "R1 and R2 contain different numbers of records; "
                        "every read must have its mate"
                    )
                if record_no == 0:
                    raise MergeValidationError("input files contain no records")
                return
            record_no += 1
            if record_no > config.MERGE_MAX_PAIRS:
                raise MergeValidationError(
                    f"merge limit exceeded: at most {config.MERGE_MAX_PAIRS} "
                    "pairs per request"
                )
            id1, seq1, qual1 = r1
            id2, seq2, qual2 = r2
            if id1 != id2:
                raise MergeValidationError(
                    f"record {record_no}: read id mismatch between R1 "
                    f"({id1!r}) and R2 ({id2!r}) after normalization"
                )
            _validate_content("R1", record_no, seq1, qual1)
            _validate_content("R2", record_no, seq2, qual2)
            yield id1, seq1, qual1, reverse_complement(seq2), qual2[::-1]


def _mismatch_count(a: bytes, b: bytes) -> int:
    """Count differing columns with an early exit past the allowed budget."""
    mismatches = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            mismatches += 1
            if mismatches > config.MERGE_MAX_MISMATCHES:
                return mismatches
    return mismatches


def find_overlap(r1: bytes, r2rc: bytes) -> int | None:
    """Return the unique longest admissible overlap length, or ``None``.

    Enumerates every shift in which at least ``MERGE_MIN_OVERLAP`` columns of
    R1 suffix and rc(R2) prefix align. Candidates with at most
    ``MERGE_MAX_MISMATCHES`` mismatches are admissible; the longest
    admissible length is returned only when it occurs for exactly one shift,
    otherwise (no candidate or a longest-length tie) ``None`` is returned.
    """
    min_ov = config.MERGE_MIN_OVERLAP
    n1, n2 = len(r1), len(r2rc)
    max_shift = n1 - min_ov
    if max_shift < 0:
        return None
    best_len = 0
    best_tied = False
    for shift in range(max_shift + 1):
        ov = min(n1 - shift, n2)
        if ov < min_ov:
            continue
        if _mismatch_count(
            r1[shift:shift + ov], r2rc[:ov]
        ) <= config.MERGE_MAX_MISMATCHES:
            if ov > best_len:
                best_len = ov
                best_tied = False
            elif ov == best_len:
                best_tied = True
    if best_len == 0 or best_tied:
        return None
    return best_len


def merge_pair(
    read_id: str,
    r1: bytes,
    q1: bytes,
    r2rc: bytes,
    q2rev: bytes,
) -> bytes | None:
    """Build the merged four-line FASTQ record, or ``None`` when skipped.

    Overlap column policy: equal bases keep the higher quality (ties keep
    R1's); disagreeing bases take the higher-quality side's base and quality;
    an equal-quality disagreement is written as ``N`` with quality ``!``.
    Non-overlapping flanks are concatenated in alignment direction.
    """
    overlap = find_overlap(r1, r2rc)
    if overlap is None:
        return None

    n1, n2 = len(r1), len(r2rc)
    # Place R1 starting at column 0; rc(R2) starts at n1 - overlap. The
    # admissible shift family always keeps rc(R2)'s start inside R1, so
    # r2_start >= 0 and columns split into R1-only / overlap / rc(R2)-only.
    r2_start = n1 - overlap
    merged_len = max(n1, r2_start + n2)

    seq = bytearray(merged_len)
    qual = bytearray(merged_len)

    # R1-only prefix (columns before rc(R2) begins), including any dovetail.
    seq[:n1] = r1
    qual[:n1] = q1

    for j in range(n2):
        col = r2_start + j
        if col >= n1:
            # rc(R2)-only suffix: no R1 base at this column.
            seq[col] = r2rc[j]
            qual[col] = q2rev[j]
            continue
        b1, b2 = r1[col], r2rc[j]
        qa, qb = q1[col], q2rev[j]
        if b1 == b2:
            # Same base: keep the higher quality; ties keep R1's.
            if qb > qa:
                seq[col], qual[col] = b2, qb
            else:
                seq[col], qual[col] = b1, qa
        elif qa > qb:
            seq[col], qual[col] = b1, qa
        elif qb > qa:
            seq[col], qual[col] = b2, qb
        else:
            # Equal-quality disagreement: ambiguous base, zeroed quality.
            seq[col], qual[col] = ord("N"), ord("!")

    header = ("@" + read_id + "\n").encode("latin-1")
    return b"".join((header, bytes(seq), b"\n+\n", bytes(qual), b"\n"))


def merge_task_files(
    r1_path: Path, r2_path: Path
) -> tuple[bytes, int, int]:
    """Merge the published FASTQ pair of one succeeded task.

    Returns ``(merged_fastq_bytes, merged_count, skipped_count)``. Raises
    ``MergeValidationError`` on any defect, in which case nothing partial is
    returned. The input files are only read, never altered.
    """
    if not r1_path.is_file():
        raise MergeValidationError("R1 result file is missing for this task")
    if not r2_path.is_file():
        raise MergeValidationError("R2 result file is missing for this task")

    chunks: list[bytes] = []
    merged_count = 0
    skipped_count = 0
    for read_id, seq1, qual1, seq2rc, qual2rev in _iter_paired_records(
        r1_path, r2_path
    ):
        record = merge_pair(read_id, seq1, qual1, seq2rc, qual2rev)
        if record is None:
            skipped_count += 1
            continue
        chunks.append(record)
        merged_count += 1
    return b"".join(chunks), merged_count, skipped_count
