"""Regressions for transient Windows locks on atomic JSON replacement."""
import errno
import json
import os
from pathlib import Path

import numpy as np
import pytest


from copa2026 import figures as workflow


def windows_error(code):
    error = PermissionError(errno.EACCES, "simulated Windows replacement lock")
    error.winerror = code
    return error


@pytest.mark.parametrize("code", [5, 32, 33])
def test_transient_lock_retries_atomic_replacement_without_rng_draws(tmp_path, monkeypatch, code):
    target = tmp_path / "progress.json"
    target.write_text('{"CompletedJobs":922}', encoding="utf8")
    original_replace = Path.replace
    calls, sleeps = [], []

    def replacement(source, destination):
        calls.append(source)
        # The old published file remains complete while replacement is denied.
        assert json.loads(target.read_text(encoding="utf8")) == {"CompletedJobs": 922}
        assert json.loads(source.read_text(encoding="utf8")) == {"CompletedJobs": 923}
        if len(calls) <= 2:
            raise windows_error(code)
        return original_replace(source, destination)

    monkeypatch.setattr(Path, "replace", replacement)
    monkeypatch.setattr(workflow.time, "sleep", sleeps.append)
    before = np.random.get_state()
    workflow._atomic_json(target, {"CompletedJobs": 923})
    after = np.random.get_state()
    assert json.loads(target.read_text(encoding="utf8")) == {"CompletedJobs": 923}
    assert len(calls) == 3 and len(set(calls)) == 1
    assert sleeps == [0.05, 0.05]
    assert not list(tmp_path.glob("*.tmp"))
    assert before[0] == after[0] and before[2:] == after[2:]
    np.testing.assert_array_equal(before[1], after[1])


def test_persistent_lock_fails_bounded_and_retains_old_json_and_temporary(tmp_path, monkeypatch):
    target = tmp_path / "progress.json"
    target.write_text('{"CompletedJobs":922}', encoding="utf8")
    calls, sleeps = [], []

    def replacement(source, destination):
        calls.append(source)
        raise windows_error(32)

    monkeypatch.setattr(Path, "replace", replacement)
    monkeypatch.setattr(workflow.time, "sleep", sleeps.append)
    with pytest.raises(PermissionError) as caught:
        workflow._atomic_json(target, {"CompletedJobs": 923})
    assert caught.value.winerror == 32
    assert len(calls) == 50 and len(sleeps) == 49
    assert sum(sleeps) < 2.5
    assert json.loads(target.read_text(encoding="utf8")) == {"CompletedJobs": 922}
    leftovers = list(tmp_path.glob("progress.json.*.tmp"))
    assert len(leftovers) == 1
    assert json.loads(leftovers[0].read_text(encoding="utf8")) == {"CompletedJobs": 923}


@pytest.mark.parametrize("error", [
    OSError(errno.ENOSPC, "simulated disk full"),
    PermissionError(errno.EACCES, "not a Windows transient code"),
])
def test_unrelated_errors_propagate_immediately(tmp_path, monkeypatch, error):
    target = tmp_path / "progress.json"
    target.write_text('{"CompletedJobs":922}', encoding="utf8")
    calls, sleeps = [], []

    def replacement(source, destination):
        calls.append(source)
        raise error

    monkeypatch.setattr(Path, "replace", replacement)
    monkeypatch.setattr(workflow.time, "sleep", sleeps.append)
    with pytest.raises(type(error)):
        workflow._atomic_json(target, {"CompletedJobs": 923})
    assert len(calls) == 1 and not sleeps
    assert json.loads(target.read_text(encoding="utf8")) == {"CompletedJobs": 922}


@pytest.mark.skipif(os.name != "nt", reason="Native Windows sharing regression")
def test_real_windows_reader_without_delete_sharing_is_released(tmp_path):
    """Use a temporary test file, never the experiment's live progress file."""
    import ctypes
    from ctypes import wintypes
    import threading

    target = tmp_path / "progress.json"
    target.write_text('{"CompletedJobs":922}', encoding="utf8")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    create = kernel.CreateFileW
    create.argtypes = [wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD,
                       ctypes.c_void_p, wintypes.DWORD, wintypes.DWORD, wintypes.HANDLE]
    create.restype = wintypes.HANDLE
    close = kernel.CloseHandle
    close.argtypes, close.restype = [wintypes.HANDLE], wintypes.BOOL
    # GENERIC_READ, FILE_SHARE_READ|FILE_SHARE_WRITE (no FILE_SHARE_DELETE),
    # OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL.
    handle = create(str(target), 0x80000000, 3, None, 3, 0x80, None)
    if handle == ctypes.c_void_p(-1).value:
        raise ctypes.WinError(ctypes.get_last_error())
    handle_lock, handle_closed = threading.Lock(), [False]

    def close_reader():
        with handle_lock:
            if not handle_closed[0]:
                close(handle)
                handle_closed[0] = True

    release = None
    try:
        probe = tmp_path / "replacement_probe.tmp"
        probe.write_text('{"probe":true}', encoding="utf8")
        # Establish that this actual reader prevents atomic replacement.
        with pytest.raises(PermissionError):
            probe.replace(target)
        probe.unlink()
        release = threading.Timer(0.15, close_reader)
        release.start()
        try:
            workflow._atomic_json(target, {"CompletedJobs": 923})
        finally:
            release.join(timeout=3)
        assert not release.is_alive()
        assert json.loads(target.read_text(encoding="utf8")) == {"CompletedJobs": 923}
    finally:
        if release is not None:
            release.cancel()
            release.join(timeout=3)
        close_reader()
