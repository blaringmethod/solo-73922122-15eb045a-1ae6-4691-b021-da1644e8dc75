"""Runtime configuration for the paired FASTQ service.

All limits can be overridden through environment variables, but the defaults
match the delivery requirement: each uploaded file <= 256 MiB and the whole
process is expected to live inside a 256 MiB memory cap (see
docker-compose.yml).
"""

from __future__ import annotations

import os
from pathlib import Path

MAX_FILE_BYTES = int(os.environ.get("FASTQ_MAX_FILE_BYTES", 256 * 1024 * 1024))
"""Maximum size of a single uploaded FASTQ file (256 MiB)."""

DATA_DIR = Path(os.environ.get("FASTQ_DATA_DIR", "/data"))
"""Root directory for the SQLite database and all task files."""

DB_PATH = DATA_DIR / "tasks.db"
TASKS_DIR = DATA_DIR / "tasks"

# Uploads are spooled to a temporary file once they exceed this size, so a big
# multipart upload never accumulates in RAM.
UPLOAD_SPOOL_THRESHOLD = int(
    os.environ.get("FASTQ_UPLOAD_SPOOL_THRESHOLD", 1 * 1024 * 1024)
)

# Records are inserted in batches; only one batch is held in memory at a time.
INSERT_BATCH_SIZE = int(os.environ.get("FASTQ_INSERT_BATCH_SIZE", 512))

# SQLite page cache kept deliberately small. The heavy data stays on disk.
SQLITE_CACHE_PAGES = int(os.environ.get("FASTQ_SQLITE_CACHE_PAGES", 2048))
SQLITE_BUSY_TIMEOUT_MS = 30_000

# How often the worker refreshes the progress percentage in the database.
PROGRESS_UPDATE_EVERY_BYTES = 8 * 1024 * 1024

STAGE_UPLOAD = "upload"
STAGE_PARSE = "parse"
STAGE_VALIDATE = "validate"
STAGE_PUBLISH = "publish"
STAGE_DONE = "done"

# --- Amplicon merge endpoint limits (all overridable via environment) ---
MERGE_MAX_PAIRS = int(os.environ.get("FASTQ_MERGE_MAX_PAIRS", 5000))
"""Maximum number of R1/R2 pairs accepted by one merge request."""

MERGE_MIN_READ_LEN = int(os.environ.get("FASTQ_MERGE_MIN_READ_LEN", 40))
MERGE_MAX_READ_LEN = int(os.environ.get("FASTQ_MERGE_MAX_READ_LEN", 200))
"""Each end must carry between 40 and 200 ACGT bases."""

MERGE_MIN_OVERLAP = int(os.environ.get("FASTQ_MERGE_MIN_OVERLAP", 20))
"""Shortest admissible R1/rc(R2) overlap, in bases."""

MERGE_MAX_MISMATCHES = int(os.environ.get("FASTQ_MERGE_MAX_MISMATCHES", 1))
"""At most one mismatch is allowed inside an accepted overlap."""

MERGE_QUAL_MIN = 33
MERGE_QUAL_MAX = 126
"""Phred+33 quality bytes span ASCII 33 ('!') through 126 ('~')."""

MERGE_MAX_LINE_BYTES = int(
    os.environ.get("FASTQ_MERGE_MAX_LINE_BYTES", 1024 * 1024)
)
"""Defensive per-line cap while streaming the published FASTQ files."""
