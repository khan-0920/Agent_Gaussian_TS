"""Own-process lifecycle only; shell=False and private Gaussian scratch directory."""
import os
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from .models import finite


@dataclass
class JobHandle:
    process: subprocess.Popen
    directory: Path
    started: float
    streams: tuple
    cancelled: bool = False
    result: object = None


class LocalRunner:
    synthetic = False

    def __init__(self, command):
        if not isinstance(command, (list, tuple)) or not command or not all(isinstance(x, str) and x for x in command):
            raise ValueError('Gaussian command must be an explicit argv list')
        self.command = list(command)
        self._owned = {}

    def start(self, directory, input_text):
        directory = Path(directory)
        directory.mkdir(parents=False, exist_ok=False)
        (directory / 'scratch').mkdir()
        (directory / 'input.gjf').write_text(input_text, encoding='utf-8')
        source = (directory / 'input.gjf').open('rb')
        output = (directory / 'output.log').open('wb')
        env = dict(os.environ, GAUSS_SCRDIR=str((directory / 'scratch').resolve()))
        try:
            process = subprocess.Popen(self.command, stdin=source, stdout=output, stderr=subprocess.STDOUT,
                                       cwd=str(directory), env=env, shell=False, start_new_session=True)
        except BaseException:
            source.close()
            output.close()
            raise
        handle = JobHandle(process, directory, time.monotonic(), (source, output))
        self._owned[id(handle)] = handle
        return handle

    def get_job_status(self, handle):
        self._check(handle)
        code = handle.process.poll()
        return 'running' if code is None else ('completed' if code == 0 else 'failed')

    def _check(self, handle):
        if self._owned.get(id(handle)) is not handle:
            raise ValueError('job handle is not owned by this runner')

    def cancel(self, handle):
        self._check(handle)
        handle.cancelled = True
        if handle.process.poll() is None:
            try:
                os.killpg(handle.process.pid, signal.SIGTERM)
                handle.process.wait(timeout=2.)
            except subprocess.TimeoutExpired:
                os.killpg(handle.process.pid, signal.SIGKILL)
                handle.process.wait()
            except ProcessLookupError:
                handle.process.wait()

    def wait(self, handle, timeout):
        self._check(handle)
        timeout = finite(timeout, 'job timeout')
        if timeout <= 0:
            raise ValueError('timeout must be positive')
        if handle.result is not None:
            return handle.result
        timed_out = False
        try:
            handle.process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            self.cancel(handle)
        except BaseException:
            self.cancel(handle)
            raise
        finally:
            for stream in handle.streams:
                stream.close()
        status = 'timed_out' if timed_out else ('cancelled' if handle.cancelled else ('completed' if handle.process.returncode == 0 else 'failed'))
        handle.result = {'status': status, 'returncode': handle.process.returncode, 'wall_seconds': time.monotonic()-handle.started,
                         'log_path': str(handle.directory / 'output.log')}
        return handle.result

    def run(self, directory, input_text, timeout):
        return self.wait(self.start(directory, input_text), timeout)
