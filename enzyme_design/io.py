"""Small, dependency-free artifact primitives shared by controller and workers."""
import hashlib
import json
import math
import os
import re
from pathlib import Path


def read_json(path):
    def invalid(value):
        raise ValueError('non-finite JSON number: ' + value)
    with Path(path).open(encoding='utf-8') as stream:
        return json.load(stream, parse_constant=invalid)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp-' + str(os.getpid()))
    try:
        with temporary.open('w', encoding='utf-8') as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def safe_name(value, label='name'):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', value):
        raise ValueError(label + ' must contain only letters, digits, underscore, dot or hyphen')
    return value


def finite(value, label, minimum=None, maximum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(label + ' must be finite')
    if minimum is not None and value < minimum or maximum is not None and value > maximum:
        raise ValueError(label + ' is outside its allowed range')
    return value


def integer(value, label, minimum=1, maximum=None):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(label + ' must be an integer')
    return finite(value, label, minimum, maximum)


def within(root, relative):
    root = Path(root).resolve()
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise ValueError('artifact paths must be relative')
    path = (root / relative).resolve()
    if path == root or root not in path.parents:
        raise ValueError('artifact path escapes its run directory')
    return path


def object_fields(value, required, optional=(), label='object'):
    if not isinstance(value, dict) or set(required) - set(value) or set(value) - set(required) - set(optional):
        raise ValueError(label + ' has missing or unknown fields')
    return value
