#!/usr/bin/env python3
"""Atomically maintain a one-row progress snapshot for compare.sh."""

from __future__ import annotations

import argparse
import datetime as dt
import math
import os
import re
import signal
import sys
import threading
import time
from pathlib import Path


TERMINAL_STATUSES = {
    "OK",
    "Timeout",
    "LoadTimeout",
    "Killed",
    "RunError",
    "RunnerError",
    "ParseError",
    "ConfigError",
}
INTEGER_RE = re.compile(r"^-?[0-9]+$")
UNSIGNED_RE = re.compile(r"^[0-9]+$")
DECIMAL_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)?$")

HEADER = (
    "updated_at_utc",
    "state",
    "total",
    "completed",
    "succeeded",
    "failed",
    "pending",
    "session_elapsed_seconds",
    "completion_rate_tasks_per_second",
    "eta_seconds",
    "estimated_finish_utc",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tasks", required=True, type=Path)
    parser.add_argument("--status-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stop-file", required=True, type=Path)
    parser.add_argument("--ready-file", required=True, type=Path)
    parser.add_argument("--baseline-completed", required=True, type=int)
    parser.add_argument("--interval", type=float, default=30.0)
    args = parser.parse_args()
    if args.baseline_completed < 0:
        parser.error("--baseline-completed must be nonnegative")
    if args.interval <= 0:
        parser.error("--interval must be positive")
    return args


def read_tasks(path: Path) -> list[tuple[Path, str]]:
    tasks: list[tuple[Path, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"invalid task manifest row {line_number}: expected 9 fields"
                )
            case_id, _, _, algorithm, _, _, _, _, expected_output = fields
            tasks.append((Path(f"{case_id}_{algorithm}.status"), expected_output))
    return tasks


def valid_terminal_status(path: Path, expected_output: str) -> str | None:
    try:
        fields = path.read_text(encoding="utf-8").rstrip("\n").split("\t")
    except (FileNotFoundError, OSError, UnicodeError):
        return None
    if len(fields) != 19 or fields[0] != "SSM_TASK_V2":
        return None
    status = fields[1]
    output = fields[16]
    if status not in TERMINAL_STATUSES or output != expected_output:
        return None
    if not INTEGER_RE.fullmatch(fields[2]):
        return None
    if fields[3] not in {"0", "1"} or fields[4] not in {"0", "1"}:
        return None
    for index in (5, 9, 10, 11, 14):
        if fields[index] != "NA" and not UNSIGNED_RE.fullmatch(fields[index]):
            return None
    for index in (6, 7, 8, 12, 13):
        if fields[index] != "NA" and not DECIMAL_RE.fullmatch(fields[index]):
            return None
    if fields[15] not in {"NA", "load", "algorithm"}:
        return None
    if fields[18] not in {"valid", "invalid"} or not fields[17]:
        return None
    if status == "OK" and (
        fields[3] != "1" or fields[4] != "1" or fields[18] != "valid"
    ):
        return None
    if fields[3] == "1" and (
        not UNSIGNED_RE.fullmatch(fields[5])
        or not DECIMAL_RE.fullmatch(fields[7])
        or not DECIMAL_RE.fullmatch(fields[12])
        or not DECIMAL_RE.fullmatch(fields[13])
    ):
        return None
    if not Path(output).is_file():
        return None
    return status


def scan(
    tasks: list[tuple[Path, str]], status_dir: Path
) -> tuple[list[tuple[Path, str]], int, int]:
    remaining: list[tuple[Path, str]] = []
    succeeded = 0
    failed = 0
    for relative_status, expected_output in tasks:
        status = valid_terminal_status(status_dir / relative_status, expected_output)
        if status == "OK":
            succeeded += 1
        elif status is not None:
            failed += 1
        else:
            remaining.append((relative_status, expected_output))
    return remaining, succeeded, failed


def utc_timestamp(epoch_seconds: float | None = None) -> str:
    if epoch_seconds is None:
        timestamp = dt.datetime.now(dt.timezone.utc)
    else:
        timestamp = dt.datetime.fromtimestamp(epoch_seconds, dt.timezone.utc)
    return timestamp.isoformat(timespec="seconds").replace("+00:00", "Z")


def atomic_write(path: Path, values: tuple[object, ...]) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    payload = "\t".join(HEADER) + "\n" + "\t".join(map(str, values)) + "\n"
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def snapshot(
    args: argparse.Namespace,
    total: int,
    succeeded: int,
    failed: int,
    started: float,
    forced_state: str | None = None,
) -> int:
    completed = succeeded + failed
    pending = max(0, total - completed)
    elapsed = max(0.0, time.monotonic() - started)
    new_completed = max(0, completed - args.baseline_completed)
    rate = new_completed / elapsed if new_completed > 0 and elapsed > 0 else 0.0

    if forced_state is not None:
        state = forced_state
    elif pending == 0:
        state = "completed_with_failures" if failed else "completed"
    else:
        state = "running"

    if pending == 0:
        eta: int | str = 0
        estimated_finish = utc_timestamp()
    elif rate > 0:
        eta = math.ceil(pending / rate)
        estimated_finish = utc_timestamp(time.time() + eta)
    else:
        eta = "NA"
        estimated_finish = "NA"

    values = (
        utc_timestamp(),
        state,
        total,
        completed,
        succeeded,
        failed,
        pending,
        f"{elapsed:.1f}",
        f"{rate:.6f}",
        eta,
        estimated_finish,
    )
    atomic_write(args.output, values)
    return pending


def main() -> int:
    args = parse_args()
    try:
        tasks = read_tasks(args.tasks)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"progress monitor error: {error}", file=sys.stderr)
        return 1

    wake_event = threading.Event()
    stop_requested = False

    def request_stop(_signum: int, _frame: object) -> None:
        nonlocal stop_requested
        stop_requested = True
        wake_event.set()

    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)

    args.ready_file.write_text("ready\n", encoding="utf-8")

    started = time.monotonic()
    total = len(tasks)
    pending_tasks = tasks
    succeeded = 0
    failed = 0
    while True:
        pending_tasks, succeeded_delta, failed_delta = scan(
            pending_tasks, args.status_dir
        )
        succeeded += succeeded_delta
        failed += failed_delta
        pending = snapshot(args, total, succeeded, failed, started)
        if pending == 0:
            return 0
        if stop_requested:
            final_state = "incomplete" if args.stop_file.exists() else "interrupted"
            snapshot(
                args, total, succeeded, failed, started, forced_state=final_state
            )
            return 0
        wake_event.wait(args.interval)
        wake_event.clear()


if __name__ == "__main__":
    raise SystemExit(main())
