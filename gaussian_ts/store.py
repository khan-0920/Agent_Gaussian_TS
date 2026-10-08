"""Private run directory, atomic JSON checkpoints and single-writer locking."""
import fcntl
import json
import os
import uuid
from contextlib import contextmanager
from pathlib import Path


class Store:
    def __init__(self, root):
        raw = Path(root)
        if raw.is_symlink():
            raise ValueError('run directory cannot be a symlink')
        self.root = raw.resolve()
        if not self.root.is_dir():
            raise ValueError('run directory does not exist')

    def path(self, relative):
        p = Path(relative)
        if p.is_absolute() or '..' in p.parts or not p.parts:
            raise ValueError('path must be relative and stay inside run directory')
        result = self.root / p
        try:
            result.resolve().relative_to(self.root)
        except ValueError:
            raise ValueError('path escapes run directory') from None
        cursor = self.root
        for part in p.parts:
            cursor = cursor / part
            if cursor.is_symlink():
                raise ValueError('symlinks are not accepted in run artifacts')
        return result

    def read_json(self, relative):
        with self.path(relative).open(encoding='utf-8') as source:
            return json.load(source)

    def write_json(self, relative, data):
        self.write_text(relative, json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

    def write_text(self, relative, text):
        target = self.path(relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        temp = target.with_name('.' + target.name + '.' + uuid.uuid4().hex)
        try:
            with temp.open('x', encoding='utf-8') as output:
                output.write(text)
                output.flush()
                os.fsync(output.fileno())
            os.replace(str(temp), str(target))
        finally:
            if temp.exists():
                temp.unlink()

    @contextmanager
    def lock(self):
        with self.path('.lock').open('a') as lockfile:
            try:
                fcntl.flock(lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError('session is busy; only one process may mutate a run') from None
            try:
                yield
            finally:
                fcntl.flock(lockfile, fcntl.LOCK_UN)
