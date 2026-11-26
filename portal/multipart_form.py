"""Stdlib-only multipart/form-data and urlencoded POST parser (Python 3.11-3.14+).

Replaces deprecated `cgi.FieldStorage` (removed in Python 3.13) with a strict
size-bounded parser backed by `email.parser.BytesParser` and `urllib.parse`.
"""

from __future__ import annotations

import io
import urllib.parse
from email.parser import BytesParser
from email.policy import default as default_policy


class FormPart:
    """Single parsed form field or uploaded file part."""

    __slots__ = ("name", "filename", "_data")

    def __init__(self, name: str, data: bytes, filename: str = ""):
        self.name = name
        self.filename = filename or ""
        self._data = data

    @property
    def file(self) -> io.BytesIO:
        return io.BytesIO(self._data)

    @property
    def value(self) -> str:
        return self._data.decode("utf-8", errors="replace")

    def read(self) -> bytes:
        return self._data


class FormData:
    """Mapping-like container compatible with the subset of FieldStorage used by SUSS Portal."""

    def __init__(self, parts: list[FormPart]):
        self._by_name: dict[str, list[FormPart]] = {}
        for part in parts:
            self._by_name.setdefault(part.name, []).append(part)

    def __contains__(self, key: object) -> bool:
        return isinstance(key, str) and key in self._by_name

    def __getitem__(self, key: str) -> FormPart:
        items = self._by_name[key]
        return items[0]

    def getfirst(self, key: str, default: str = "") -> str:
        items = self._by_name.get(key)
        if not items:
            return default
        return items[0].value

    def getlist(self, key: str) -> list[str]:
        return [item.value for item in self._by_name.get(key, [])]


def parse_form(fp, headers, *, max_bytes: int) -> FormData:
    """Read and parse an HTTP POST body bounded by `max_bytes`."""
    raw_len = headers.get("Content-Length")
    try:
        length = int(raw_len) if raw_len is not None else 0
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid Content-Length header") from exc
    if length < 0 or length > max_bytes:
        raise ValueError(f"Request body exceeds limit ({max_bytes} bytes)")

    body = fp.read(length) if length > 0 else b""
    if len(body) > max_bytes:
        raise ValueError(f"Request body exceeds limit ({max_bytes} bytes)")

    ctype = str(headers.get("Content-Type", "") or "")
    ctype_lower = ctype.lower()

    if ctype_lower.startswith("multipart/form-data"):
        header_block = (
            f"Content-Type: {ctype}\r\nMIME-Version: 1.0\r\n\r\n"
        ).encode("utf-8")
        message = BytesParser(policy=default_policy).parsebytes(header_block + body)
        parts: list[FormPart] = []
        if message.is_multipart():
            for sub in message.iter_parts():
                name = sub.get_param("name", header="content-disposition")
                if not name:
                    continue
                filename = sub.get_filename() or ""
                payload = sub.get_payload(decode=True)
                if payload is None:
                    raw_payload = sub.get_payload()
                    payload = (
                        raw_payload.encode("utf-8", errors="replace")
                        if isinstance(raw_payload, str)
                        else b""
                    )
                parts.append(FormPart(name=str(name), data=payload, filename=str(filename)))
        return FormData(parts)

    if ctype_lower.startswith("application/x-www-form-urlencoded") or not ctype_lower:
        parsed = urllib.parse.parse_qs(
            body.decode("utf-8", errors="replace"),
            keep_blank_values=True,
        )
        parts = [
            FormPart(name=k, data=v.encode("utf-8"))
            for k, vals in parsed.items()
            for v in vals
        ]
        return FormData(parts)

    raise ValueError(f"Unsupported Content-Type: {ctype}")
