#!/usr/bin/env python3
"""Run one benchmark with separate loading and algorithm timeouts.

When --load-timeout is present, the child receives a private READY/ACK pipe
pair.  The loading deadline applies until the child reports that both graphs
are loaded; only then is the algorithm deadline started.  Stdout remains
directly connected to the result file, so checkpoint volume cannot block the
control handshake.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import select
import signal
import subprocess
import sys
import time
from typing import Callable, Optional, Sequence, Tuple

try:
    import resource
except ImportError:  # pragma: no cover - compare.sh targets Unix-like hosts.
    resource = None  # type: ignore[assignment]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="run a command with loading/algorithm timeouts and resource accounting"
    )
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("--load-timeout", type=float)
    parser.add_argument("--kill-after", type=float, default=5.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be greater than zero")
    if args.load_timeout is not None and (
        not math.isfinite(args.load_timeout) or args.load_timeout <= 0
    ):
        parser.error("--load-timeout must be greater than zero")
    if not math.isfinite(args.kill_after) or args.kill_after < 0:
        parser.error("--kill-after must be non-negative")
    return args


def peak_child_rss_kb() -> Optional[int]:
    if resource is None:
        return None
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak = float(usage.ru_maxrss)
    if peak <= 0:
        return None
    if sys.platform == "darwin":
        return int(math.ceil(peak / 1024.0))
    return int(math.ceil(peak))


def signal_process_group(process: subprocess.Popen[bytes], signum: int) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        pass


def format_ms(value: Optional[float]) -> str:
    return "NA" if value is None else f"{value:.4f}"


def write_status(
    path: Path,
    *,
    timed_out: bool,
    timeout_phase: str,
    return_code: int,
    peak_rss_kb: Optional[int],
    ready_seen: bool,
    load_elapsed_ms: Optional[float],
    algorithm_elapsed_ms: Optional[float],
    interrupted_signal: Optional[int] = None,
    runner_error: str = "",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    error_value = runner_error.replace("\t", " ").replace("\r", " ").replace("\n", " ")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"timed_out={int(timed_out)}\n")
        handle.write(f"timeout_phase={timeout_phase}\n")
        handle.write(f"return_code={return_code}\n")
        handle.write(
            f"peak_rss_kb={peak_rss_kb if peak_rss_kb is not None else 'NA'}\n"
        )
        handle.write(f"ready_seen={int(ready_seen)}\n")
        handle.write(f"load_elapsed_ms={format_ms(load_elapsed_ms)}\n")
        handle.write(f"algorithm_elapsed_ms={format_ms(algorithm_elapsed_ms)}\n")
        handle.write(
            "interrupted_signal="
            f"{interrupted_signal if interrupted_signal is not None else 'NA'}\n"
        )
        handle.write(f"runner_error={error_value}\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def shell_return_code(
    return_code: int,
    timed_out: bool,
    interrupted_signal: Optional[int],
    runner_error: str,
) -> int:
    if interrupted_signal is not None:
        return 128 + interrupted_signal
    if timed_out:
        return 124
    if runner_error:
        return 125
    if return_code < 0:
        return 128 + (-return_code)
    return return_code


def wait_until(
    process: subprocess.Popen[bytes],
    timeout: float,
    was_interrupted: Optional[Callable[[], bool]] = None,
) -> Tuple[Optional[int], bool]:
    deadline = time.monotonic() + timeout
    while True:
        return_code = process.poll()
        if return_code is not None:
            return return_code, False
        if was_interrupted is not None and was_interrupted():
            return None, True
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return None, False
        time.sleep(min(0.02, remaining))


def wait_for_ready(
    process: subprocess.Popen[bytes],
    ready_fd: int,
    timeout: float,
    was_interrupted: Callable[[], bool],
) -> Tuple[str, Optional[int]]:
    """Return (ready|exited|interrupted|timeout|closed, return_code)."""
    deadline = time.monotonic() + timeout
    received = bytearray()
    while True:
        return_code = process.poll()
        if return_code is not None:
            return "exited", return_code
        if was_interrupted():
            return "interrupted", None

        remaining = deadline - time.monotonic()
        if remaining <= 0:
            # Give a READY already queued at the boundary priority over timeout.
            poll_timeout = 0.0
        else:
            poll_timeout = min(0.02, remaining)
        readable, _, _ = select.select([ready_fd], [], [], poll_timeout)
        if readable:
            chunk = os.read(ready_fd, 64)
            if not chunk:
                return_code = process.poll()
                if return_code is not None:
                    return "exited", return_code
                return "closed", None
            received.extend(chunk)
            if b"READY\n" in received:
                return "ready", None
        if time.monotonic() >= deadline:
            return "timeout", None


def terminate_process(
    process: subprocess.Popen[bytes], signum: int, kill_after: float
) -> int:
    signal_process_group(process, signum)
    return_code_or_none, _ = wait_until(process, kill_after)
    if return_code_or_none is not None:
        return return_code_or_none
    signal_process_group(process, signal.SIGKILL)
    return process.wait()


def append_termination(
    output_handle: object,
    *,
    status: str,
    timeout_phase: str,
    return_code: int,
    peak_rss_kb: Optional[int],
    load_elapsed_ms: Optional[float],
    algorithm_elapsed_ms: Optional[float],
    runner_error: str,
) -> None:
    error_value = runner_error.replace(" ", "_").replace("\t", "_").replace("\n", "_")
    line = (
        f"\nSSM_TERMINATION status={status} timeout_phase={timeout_phase} "
        f"return_code={return_code} "
        f"peak_rss_kb={peak_rss_kb if peak_rss_kb is not None else 'NA'} "
        f"load_elapsed_ms={format_ms(load_elapsed_ms)} "
        f"algorithm_elapsed_ms={format_ms(algorithm_elapsed_ms)} "
        f"runner_error={error_value or 'NA'}\n"
    ).encode("utf-8", errors="replace")
    handle = output_handle
    handle.seek(0, os.SEEK_END)  # type: ignore[attr-defined]
    handle.write(line)  # type: ignore[attr-defined]
    handle.flush()  # type: ignore[attr-defined]
    os.fsync(handle.fileno())  # type: ignore[attr-defined]


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    process: Optional[subprocess.Popen[bytes]] = None
    timed_out = False
    timeout_phase = "NA"
    return_code = 125
    interrupted_signal: Optional[int] = None
    runner_error = ""
    ready_seen = False
    process_started_at: Optional[float] = None
    ready_at: Optional[float] = None
    load_end_at: Optional[float] = None
    algorithm_end_at: Optional[float] = None
    ended_at: Optional[float] = None
    ready_read_fd: Optional[int] = None
    ready_write_fd: Optional[int] = None
    ack_read_fd: Optional[int] = None
    ack_write_fd: Optional[int] = None

    def forward_signal(signum: int, _frame: object) -> None:
        nonlocal interrupted_signal
        if interrupted_signal is None:
            interrupted_signal = signum
        if process is not None:
            signal_process_group(process, signum)

    signal.signal(signal.SIGINT, forward_signal)
    signal.signal(signal.SIGTERM, forward_signal)

    try:
        with args.output.open("wb") as output_handle:
            environment = None
            pass_fds: Tuple[int, ...] = ()
            if args.load_timeout is not None:
                ready_read_fd, ready_write_fd = os.pipe()
                ack_read_fd, ack_write_fd = os.pipe()
                environment = os.environ.copy()
                environment["SSM_READY_FD"] = str(ready_write_fd)
                environment["SSM_ACK_FD"] = str(ack_read_fd)
                pass_fds = (ready_write_fd, ack_read_fd)

            if interrupted_signal is None:
                process = subprocess.Popen(
                    args.command,
                    stdout=output_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    env=environment,
                    pass_fds=pass_fds,
                )
                process_started_at = time.monotonic()

            if ready_write_fd is not None:
                os.close(ready_write_fd)
                ready_write_fd = None
            if ack_read_fd is not None:
                os.close(ack_read_fd)
                ack_read_fd = None

            if process is None:
                assert interrupted_signal is not None
                return_code = -interrupted_signal
            elif args.load_timeout is not None:
                assert ready_read_fd is not None and ack_write_fd is not None
                state, early_return_code = wait_for_ready(
                    process,
                    ready_read_fd,
                    args.load_timeout,
                    lambda: interrupted_signal is not None,
                )
                if state == "ready":
                    ready_seen = True
                    ready_at = time.monotonic()
                    load_end_at = ready_at
                    if interrupted_signal is None:
                        os.write(ack_write_fd, b"A")
                    os.close(ack_write_fd)
                    ack_write_fd = None

                    if interrupted_signal is not None:
                        algorithm_end_at = time.monotonic()
                        return_code = terminate_process(
                            process, interrupted_signal, args.kill_after
                        )
                    else:
                        result, interrupted = wait_until(
                            process,
                            args.timeout,
                            lambda: interrupted_signal is not None,
                        )
                        algorithm_end_at = time.monotonic()
                        if result is not None:
                            return_code = result
                        elif interrupted:
                            assert interrupted_signal is not None
                            return_code = terminate_process(
                                process, interrupted_signal, args.kill_after
                            )
                        else:
                            timed_out = True
                            timeout_phase = "algorithm"
                            return_code = terminate_process(
                                process, signal.SIGTERM, args.kill_after
                            )
                elif state == "exited":
                    load_end_at = time.monotonic()
                    assert early_return_code is not None
                    return_code = early_return_code
                    if return_code == 0:
                        runner_error = "child exited successfully before READY"
                elif state == "interrupted":
                    load_end_at = time.monotonic()
                    assert interrupted_signal is not None
                    return_code = terminate_process(
                        process, interrupted_signal, args.kill_after
                    )
                elif state == "timeout":
                    load_end_at = time.monotonic()
                    timed_out = True
                    timeout_phase = "load"
                    return_code = terminate_process(
                        process, signal.SIGTERM, args.kill_after
                    )
                else:
                    load_end_at = time.monotonic()
                    runner_error = "ready control pipe closed before READY"
                    return_code = terminate_process(
                        process, signal.SIGTERM, args.kill_after
                    )
            else:
                ready_seen = True
                ready_at = process_started_at
                result, interrupted = wait_until(
                    process, args.timeout, lambda: interrupted_signal is not None
                )
                algorithm_end_at = time.monotonic()
                if result is not None:
                    return_code = result
                elif interrupted:
                    assert interrupted_signal is not None
                    return_code = terminate_process(
                        process, interrupted_signal, args.kill_after
                    )
                else:
                    timed_out = True
                    timeout_phase = "algorithm"
                    return_code = terminate_process(
                        process, signal.SIGTERM, args.kill_after
                    )

            if process is not None and process.poll() is None:
                runner_error = runner_error or "child remained alive after wait"
                return_code = terminate_process(process, signal.SIGTERM, args.kill_after)

            ended_at = time.monotonic()
            if interrupted_signal is not None:
                timed_out = False
                timeout_phase = "NA"

            if process_started_at is None or ended_at is None:
                load_elapsed_ms = None
                algorithm_elapsed_ms = None
            elif args.load_timeout is not None:
                load_end = load_end_at if load_end_at is not None else ended_at
                load_elapsed_ms = (load_end - process_started_at) * 1000.0
                algorithm_elapsed_ms = (
                    ((algorithm_end_at if algorithm_end_at is not None else ended_at) - ready_at)
                    * 1000.0
                    if ready_at is not None
                    else None
                )
            else:
                load_elapsed_ms = None
                algorithm_elapsed_ms = (
                    (algorithm_end_at if algorithm_end_at is not None else ended_at)
                    - process_started_at
                ) * 1000.0

            peak_rss_kb = peak_child_rss_kb()
            if interrupted_signal is not None:
                termination_status = "interrupted"
            elif timed_out:
                termination_status = "timeout"
            elif runner_error:
                termination_status = "runner_error"
            elif return_code == 0:
                termination_status = "completed"
            elif return_code < 0:
                termination_status = "killed"
            else:
                termination_status = "error"

            if termination_status != "completed":
                append_termination(
                    output_handle,
                    status=termination_status,
                    timeout_phase=timeout_phase,
                    return_code=return_code,
                    peak_rss_kb=peak_rss_kb,
                    load_elapsed_ms=load_elapsed_ms,
                    algorithm_elapsed_ms=algorithm_elapsed_ms,
                    runner_error=runner_error,
                )

        write_status(
            args.status,
            timed_out=timed_out,
            timeout_phase=timeout_phase,
            return_code=return_code,
            peak_rss_kb=peak_rss_kb,
            ready_seen=ready_seen,
            load_elapsed_ms=load_elapsed_ms,
            algorithm_elapsed_ms=algorithm_elapsed_ms,
            interrupted_signal=interrupted_signal,
            runner_error=runner_error,
        )
        return shell_return_code(
            return_code, timed_out, interrupted_signal, runner_error
        )
    except Exception as exc:
        runner_error = f"{type(exc).__name__}: {exc}"
        if interrupted_signal is not None:
            timed_out = False
            timeout_phase = "NA"
        if process is not None and process.poll() is None:
            signal_process_group(process, signal.SIGKILL)
            try:
                return_code = process.wait(timeout=1.0)
            except subprocess.TimeoutExpired:
                pass
        ended_at = time.monotonic()
        load_elapsed_ms = None
        algorithm_elapsed_ms = None
        if process_started_at is not None:
            if ready_at is not None:
                load_elapsed_ms = (ready_at - process_started_at) * 1000.0
                algorithm_elapsed_ms = (ended_at - ready_at) * 1000.0
            elif args.load_timeout is not None:
                load_elapsed_ms = (ended_at - process_started_at) * 1000.0
            else:
                algorithm_elapsed_ms = (ended_at - process_started_at) * 1000.0
        peak_rss_kb = peak_child_rss_kb()
        try:
            with args.output.open("ab") as output_handle:
                append_termination(
                    output_handle,
                    status="interrupted" if interrupted_signal is not None else "runner_error",
                    timeout_phase=timeout_phase,
                    return_code=return_code,
                    peak_rss_kb=peak_rss_kb,
                    load_elapsed_ms=load_elapsed_ms,
                    algorithm_elapsed_ms=algorithm_elapsed_ms,
                    runner_error=runner_error,
                )
        except OSError:
            pass
        write_status(
            args.status,
            timed_out=timed_out,
            timeout_phase=timeout_phase,
            return_code=return_code,
            peak_rss_kb=peak_rss_kb,
            ready_seen=ready_seen,
            load_elapsed_ms=load_elapsed_ms,
            algorithm_elapsed_ms=algorithm_elapsed_ms,
            interrupted_signal=interrupted_signal,
            runner_error=runner_error,
        )
        return shell_return_code(
            return_code, timed_out, interrupted_signal, runner_error
        )
    finally:
        for fd in (ready_read_fd, ready_write_fd, ack_read_fd, ack_write_fd):
            if fd is not None:
                try:
                    os.close(fd)
                except OSError:
                    pass


if __name__ == "__main__":
    raise SystemExit(main())
