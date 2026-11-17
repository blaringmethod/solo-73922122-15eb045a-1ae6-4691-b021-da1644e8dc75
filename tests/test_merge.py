"""Functional self-tests for the amplicon merge endpoint and merge algorithm.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_merge.py
Each test run gets an isolated FASTQ_DATA_DIR under tests/output.
"""

from __future__ import annotations

import asyncio
import io
import os
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402
from app.merge import (  # noqa: E402
    find_overlap,
    merge_pair,
    merge_task_files,
    reverse_complement,
)

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes,
                suffix: str = "") -> bytes:
    header = b"@" + rid.encode() + suffix.encode()
    return b"\n".join((header, seq, b"+", quality)) + b"\n"


# Backwards-style short alias used throughout this file.
rec = make_record


def fq(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.02)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload(client, r1: io.BytesIO, r2: io.BytesIO) -> dict:
    r1.seek(0)
    r2.seek(0)
    r = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", r1, "application/octet-stream"),
               "r2": ("r2.fastq", r2, "application/octet-stream")},
    )
    assert r.status_code == 202, r.text
    return await wait_terminal(client, r.json()["id"])


def parse_merged(text: str) -> list[tuple[str, bytes, bytes]]:
    lines = text.splitlines()
    assert len(lines) % 4 == 0
    out = []
    for i in range(0, len(lines), 4):
        assert lines[i].startswith("@")
        assert lines[i + 2] == "+"
        rid = lines[i][1:].split()[0]
        out.append((rid, lines[i + 1].encode(), lines[i + 3].encode()))
    return out


import hashlib

# Deterministic, well-mixed pseudo-random ACGT pattern so shifted windows do
# not accidentally satisfy the <=1-mismatch overlap rule.
def pseudo_seq(n: int, seed: int = 0) -> bytes:
    out = bytearray()
    i = 0
    while len(out) < n:
        out.extend(hashlib.sha256(f"{seed}-{i}".encode()).digest())
        i += 1
    return bytes(b"ACGT"[b % 4] for b in bytes(out[:n]))


def algorithm_unit_checks() -> None:
    # reverse complement
    record("reverse_complement basic",
           reverse_complement(b"ACGTN") == b"NACGT",
           reverse_complement(b"ACGTN"))

    # --- overlap detection: hand-built 60bp pair sharing a 30bp overlap ---
    r1 = b"A" * 30 + b"C" * 30                       # length 60
    r2rc = b"C" * 30 + b"G" * 30                     # length 60, shift=30
    ov = find_overlap(r1, r2rc)
    record("unique longest overlap found (30)", ov == 30, str(ov))

    # No admissible overlap: random 40bp reads, min overlap 20.
    a = pseudo_seq(40, seed=1)
    b = pseudo_seq(40, seed=9)
    record("no overlap returns None", find_overlap(a, b) is None)

    # Tie: r2rc fully contained in a homogeneous R1 -> every shift matching
    # the containment window yields the same maximum overlap length.
    r1_long = b"A" * 100
    r2_short_rc = b"A" * 60
    record("longest-length tie returns None",
           find_overlap(r1_long, r2_short_rc) is None)

    # One mismatch inside an otherwise perfect overlap is accepted. Random
    # flanks make the intended shift the only admissible candidate.
    flank_l = pseudo_seq(25, seed=21)
    flank_r = pseudo_seq(25, seed=22)
    shared = bytearray(pseudo_seq(30, seed=23))
    r1m = flank_l + bytes(shared)
    r2m = bytearray(shared) + flank_r
    r2m[0] = ord("T") if r2m[0] != ord("T") else ord("G")  # one mismatch
    record("single-mismatch overlap accepted",
           find_overlap(r1m, bytes(r2m)) == 30)

    # Two mismatches reject that candidate; no other candidate exists.
    r2m2 = bytearray(shared) + flank_r
    r2m2[0] = ord("T") if r2m2[0] != ord("T") else ord("G")
    r2m2[1] = ord("T") if r2m2[1] != ord("T") else ord("G")
    record("two-mismatch overlap rejected",
           find_overlap(r1m, bytes(r2m2)) is None)

    # --- merge sequence assembly: flank + overlap + flank ---
    q1 = b"E" * 60   # Phred 37
    q2 = b"F" * 60   # Phred 38, higher on every column
    blob = merge_pair("read-x", r1, q1, r2rc, q2[::-1])
    expected_seq = b"A" * 30 + b"C" * 30 + b"G" * 30  # length 90
    expected_qual = b"E" * 30 + b"F" * 60
    want = b"@read-x\n" + expected_seq + b"\n+\n" + expected_qual + b"\n"
    record("merge assembly seq+qual correct", blob == want,
           repr(blob) + " != " + repr(want))

    # Equal bases, differing qualities: higher quality wins, ties keep R1.
    r1e = b"ACGT" * 15
    r2e = b"ACGT" * 15
    q1e = bytes([35] * 60)
    q2e = bytes([40] * 60)
    blob = merge_pair("eq", r1e, q1e, r2e, q2e[::-1])
    _, s, q = parse_merged(blob.decode("latin-1"))[0]
    record("equal bases keep higher quality",
           s == r1e and all(v == 40 for v in q), repr(q))

    # Conflicting bases: higher-quality side's base/quality is taken.
    r1c = b"ACGT" * 15
    r2c = bytearray(b"ACGT" * 15)
    r2c[0] = ord("C")  # overlap alignment starts at shift 0 here? need
    # a deterministic unique shift: use explicit 40bp overlap construction.
    r1c = b"A" * 20 + b"ACGT" * 10          # 60
    r2c_rc = bytearray(b"ACGT" * 10 + b"G" * 20)
    r2c_rc[20] = ord("T")  # mismatch at overlap column 20 (R1 base 'A')
    q1c = bytes([50] * 60)
    q2c = bytes([60] * 60)
    blob = merge_pair("conf", bytes(r1c), q1c, bytes(r2c_rc), q2c[::-1])
    _, s, q = parse_merged(blob.decode("latin-1"))[0]
    # Column 40 = overlap offset 20: R2 side higher quality -> T / 60.
    ok_conf = s[40] == ord("T") and q[40] == 60
    record("conflict resolved toward higher-quality side", ok_conf,
           f"base={chr(s[40])} qual={q[40]}")

    # Equal-quality conflict -> N and '!' (33).
    q_eq = bytes([50] * 60)
    blob = merge_pair("n", bytes(r1c), q_eq, bytes(r2c_rc), q_eq[::-1])
    _, s, q = parse_merged(blob.decode("latin-1"))[0]
    record("equal-quality conflict -> N and !",
           s[40] == ord("N") and q[40] == 33,
           f"base={chr(s[40])} qual={q[40]}")

    # Skipped pair returns None (tie case from above).
    record("merge_pair skips tied pair",
           merge_pair("t", r1_long, b"~" * 100, r2_short_rc,
                      b"~" * 60) is None)

    # Minimum overlap boundary, built from deterministic blocks so adjacent
    # shifts cannot produce accidental <=1-mismatch candidates.
    # 19 shared columns: the alignment at shift 21 is too short and the
    # 20-column window at shift 20 contains 2 mismatches -> no candidate.
    r1_19 = b"A" * 21 + b"T" * 19
    r2rc_19 = b"T" * 19 + b"G" * 21
    r2rc_19 = b"T" * 18 + b"A" + b"G" * 21  # 19-window gets 1 mm, 20-window 2
    record("overlap shorter than 20 rejected",
           find_overlap(r1_19, r2rc_19) is None)
    # 20 perfectly shared columns: shift 20 is the unique admissible candidate.
    r1_20 = b"A" * 20 + b"T" * 20
    r2rc_20 = b"T" * 20 + b"G" * 20
    record("overlap of exactly 20 accepted",
           find_overlap(r1_20, r2rc_20) == 20)


async def api_checks(client) -> None:
    # ---- 1. end-to-end merge via the API: 3 pairs, one skipped ----
    # Pair a: perfect 30bp overlap; pair b: no overlap (skip); pair c: tie.
    a1 = b"A" * 30 + b"C" * 30
    a2 = (b"C" * 30 + b"G" * 30)                       # R2 forward = rc(rc)
    r2_a = reverse_complement(a2)
    b1 = pseudo_seq(40, seed=3)
    b2 = pseudo_seq(40, seed=7)
    c1 = b"A" * 60
    c2 = reverse_complement(b"A" * 40)                 # rc stays all-A
    r1_file = fq([
        rec("pa", a1, b"E" * 60, suffix="/1"),
        rec("pb", b1, b"E" * 40),
        rec("pc", c1, b"E" * 60),
    ])
    r2_file = fq([
        rec("pa", r2_a, b"F" * 60, suffix="/2"),
        rec("pb", b2, b"F" * 40),
        rec("pc", c2, b"F" * 40),
    ])
    body = await upload(client, r1_file, r2_file)
    record("setup task succeeded", body["status"] == "succeeded",
           body.get("message", ""))
    tid = body["id"]

    # Snapshot published bytes; merge must not mutate them.
    p1 = (config.TASKS_DIR / tid / "output" / "R1.fastq").read_bytes()
    p2 = (config.TASKS_DIR / tid / "output" / "R2.fastq").read_bytes()

    mr = await client.post(f"/tasks/{tid}/merge")
    record("merge endpoint returns 200", mr.status_code == 200, mr.text)
    mbody = mr.json()
    record("merge counts (1 merged, 2 skipped)",
           mbody["merged_count"] == 1 and mbody["skipped_count"] == 2
           and mbody["total_pairs"] == 3, str(mbody))
    merged = parse_merged(mbody["merged_fastq"])
    record("merged FASTQ uses normalized read id and R1 order",
           [rid for rid, _, _ in merged] == ["pa"], str(merged))
    _, s, q = merged[0]
    record("merged sequence length 90", s == b"A" * 30 + b"C" * 30 + b"G" * 30)
    record("merged quality follows higher-quality side",
           q == b"E" * 30 + b"F" * 60, repr(q))

    # Original results untouched.
    record("published R1 unchanged by merge",
           (config.TASKS_DIR / tid / "output" / "R1.fastq").read_bytes() == p1)
    record("published R2 unchanged by merge",
           (config.TASKS_DIR / tid / "output" / "R2.fastq").read_bytes() == p2)

    # Re-running merge is deterministic and leaves no extra artifacts.
    before = sorted(p.name for p in
                    (config.TASKS_DIR / tid / "output").iterdir())
    mr2 = await client.post(f"/tasks/{tid}/merge")
    after = sorted(p.name for p in
                   (config.TASKS_DIR / tid / "output").iterdir())
    record("merge is repeatable and writes no side files",
           mr2.status_code == 200 and before == after
           and mr2.json()["merged_fastq"] == mbody["merged_fastq"],
           f"{before} vs {after}")

    # ---- 2. unknown task -> 404 ----
    r = await client.post("/tasks/deadbeef/merge")
    record("merge unknown task -> 404", r.status_code == 404, str(r.status_code))

    # ---- 3. non-succeeded task -> 409 ----
    # A truncated (three-line) record fails the pairing pipeline.
    bad_r1 = fq([rec("x", b"A" * 40, b"E" * 40)])
    bad_r2 = io.BytesIO(b"@x\n" + b"A" * 40 + b"\n+\n")  # missing quality line
    body = await upload(client, bad_r1, bad_r2)
    record("setup malformed task fails", body["status"] == "failed",
           body.get("message", ""))
    r = await client.post(f"/tasks/{body['id']}/merge")
    record("merge non-succeeded task -> 409", r.status_code == 409,
           str(r.status_code))

    # ---- 4. corrupted back segment: truncated R2 mid-record ----
    rec_g = rec("g", b"A" * 40, b"E" * 40)
    rec_h = rec("h", b"C" * 40, b"E" * 40)
    body = await upload(client, fq([rec_g, rec_h]), fq([rec_g, rec_h]))
    tid2 = body["id"]
    r2_path = config.TASKS_DIR / tid2 / "output" / "R2.fastq"
    original = r2_path.read_bytes()
    # Keep record g intact, cut record h halfway through its sequence line.
    r2_path.write_bytes(original[: len(rec_g) + 10])
    r = await client.post(f"/tasks/{tid2}/merge")
    record("truncated back record -> 422 whole-request error",
           r.status_code == 422 and "truncated" in r.json()["detail"],
           r.text)
    r2_path.write_bytes(original)

    # ---- 5. missing end file: the whole request is rejected, no partial ----
    # An orphan pair cannot reach 'succeeded' through the pairing pipeline, so
    # emulate a missing end by deleting one published file of a good task.
    good_body = await upload(client, fq([rec_g, rec_h]), fq([rec_g, rec_h]))
    tid3 = good_body["id"]
    missing = config.TASKS_DIR / tid3 / "output" / "R2.fastq"
    missing.unlink()
    r = await client.post(f"/tasks/{tid3}/merge")
    record("missing end file -> 404, whole request rejected",
           r.status_code == 404, str(r.status_code))

    # ---- 6. read id mismatch between ends -> 422 ----
    # Corrupt the second R2 header of a succeeded task so the normalized ids
    # no longer agree at that position; the whole merge request must fail.
    good_body = await upload(client, fq([
        rec("k1", b"A" * 40, b"E" * 40),
        rec("k2", b"C" * 40, b"E" * 40),
    ]), fq([
        rec("k1", b"A" * 40, b"E" * 40),
        rec("k2", b"C" * 40, b"E" * 40),
    ]))
    tid4 = good_body["id"]
    r2p = config.TASKS_DIR / tid4 / "output" / "R2.fastq"
    r2p.write_bytes(r2p.read_bytes().replace(b"@k2", b"@k9", 1))
    r = await client.post(f"/tasks/{tid4}/merge")
    record("read id mismatch -> 422",
           r.status_code == 422 and "mismatch" in r.json()["detail"], r.text)

    # ---- 7. validation: non-ACGT base, bad length, bad quality -> 422 ----
    async def expect_422(label, s1, q1, s2, q2):
        b = await upload(client, fq([rec("v", s1, q1)]),
                         fq([rec("v", s2, q2)]))
        if b["status"] != "succeeded":
            # The pairing pipeline may reject first; build a succeeded task
            # then patch its output to inject the defect deterministically.
            b = await upload(
                client, fq([rec("v", b"A" * 40, b"E" * 40)]),
                fq([rec("v", b"A" * 40, b"E" * 40)]))
            t = b["id"]
            out = config.TASKS_DIR / t / "output"
            out.joinpath("R1.fastq").write_bytes(rec("v", s1, q1))
            out.joinpath("R2.fastq").write_bytes(rec("v", s2, q2))
        else:
            t = b["id"]
        rr = await client.post(f"/tasks/{t}/merge")
        record(label, rr.status_code == 422, rr.text)

    await expect_422("non-ACGT base -> 422",
                     b"A" * 39 + b"N", b"E" * 40,
                     b"A" * 40, b"E" * 40)
    await expect_422("length below 40 -> 422",
                     b"A" * 39, b"E" * 39, b"A" * 40, b"E" * 40)
    await expect_422("length above 200 -> 422",
                     b"A" * 201, b"E" * 201, b"A" * 40, b"E" * 40)
    await expect_422("quality outside Phred+33 -> 422",
                     b"A" * 40, b"E" * 39 + bytes([32]),
                     b"A" * 40, b"E" * 40)

    # ---- 8. 5000-pair cap: 5001 pairs rejected, nothing partial ----
    big = [rec(f"r{i:05d}", b"A" * 40, b"E" * 40) for i in range(5001)]
    # Uploads themselves are fine (pairing has no such cap); merge must fail.
    # Build in chunks to keep the test's own memory modest.
    buf1 = io.BytesIO(b"".join(big))
    buf2 = io.BytesIO(b"".join(
        rec(f"r{i:05d}", b"A" * 40, b"E" * 40) for i in range(5001)))
    b = await upload(client, buf1, buf2)
    record("5001-pair task pairs successfully", b["status"] == "succeeded",
           b.get("message", ""))
    rr = await client.post(f"/tasks/{b['id']}/merge")
    record("merge enforces 5000-pair cap -> 422",
           rr.status_code == 422 and "5000" in rr.json()["detail"], rr.text)

    # Exactly 5000 pairs merges (all homogeneous all-A pairs tie/skip, but the
    # cap itself must not trigger).
    ok5000 = [rec(f"s{i:05d}", b"A" * 40, b"E" * 40) for i in range(5000)]
    b = await upload(client, io.BytesIO(b"".join(ok5000)),
                     io.BytesIO(b"".join(ok5000)))
    rr = await client.post(f"/tasks/{b['id']}/merge")
    record("5000 pairs accepted by cap", rr.status_code == 200, rr.text)

    # ---- 9. /docs and OpenAPI expose the new endpoint ----
    docs = await client.get("/docs")
    record("/docs reachable", docs.status_code == 200)
    schema = (await client.get("/openapi.json")).json()
    record("merge path present in OpenAPI",
           "/tasks/{task_id}/merge" in schema["paths"])


def direct_file_checks() -> None:
    """Exercise merge_task_files on hand-written files without the API."""
    tmp = Path(tempfile.mkdtemp(dir=str(ROOT / "tests" / "output")))
    try:
        r1 = tmp / "R1.fastq"
        r2 = tmp / "R2.fastq"
        r1.write_bytes(
            rec("x", b"A" * 30 + b"C" * 30, b"E" * 60)
            + rec("y", pseudo_seq(40, 2), b"E" * 40)
        )
        r2.write_bytes(
            rec("x", reverse_complement(b"C" * 30 + b"G" * 30), b"F" * 60)
            + rec("y", pseudo_seq(40, 5), b"F" * 40)
        )
        blob, merged, skipped = merge_task_files(r1, r2)
        record("direct merge counts", merged == 1 and skipped == 1,
               f"{merged}/{skipped}")
        record("direct merge bytes start with normalized header",
               blob.startswith(b"@x\n"))

        # Missing end file raises MergeValidationError.
        r2.unlink()
        try:
            merge_task_files(r1, r2)
            record("missing end raises", False)
        except Exception as exc:
            record("missing end raises MergeValidationError",
                   type(exc).__name__ == "MergeValidationError", str(exc))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


async def main_async() -> int:
    algorithm_unit_checks()
    direct_file_checks()

    tmp_root = Path(tempfile.mkdtemp(prefix="fastq-merge-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    config.PROGRESS_UPDATE_EVERY_BYTES = 64 * 1024

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=120) as client:
        await api_checks(client)

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [n for n, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"self peak RSS: {peak_kib / 1024:.1f} MiB (cap 256 MiB)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
