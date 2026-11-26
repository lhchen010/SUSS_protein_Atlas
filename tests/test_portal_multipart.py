"""Tests for stdlib multipart/urlencoded POST parsing and Python 3.13+ portal compatibility."""

import io
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "portal"))

from multipart_form import parse_form


def test_parse_multipart_form_with_text_and_binary_file():
    boundary = "----SUSSTestBoundary7MA4YWxkTrZu0gW"
    binary_payload = b"\x1f\x8b\x08\x00fake-tar-gz-payload\x00\xff"
    body = (
        f"--{boundary}\r\n"
        f"Content-Disposition: form-data; name=\"csrf\"\r\n\r\n"
        f"secret-token-123\r\n"
        f"--{boundary}\r\n"
        f"Content-Disposition: form-data; name=\"code\"\r\n\r\n"
        f"bc_strain\r\n"
        f"--{boundary}\r\n"
        f"Content-Disposition: form-data; name=\"pdb_tar\"; filename=\"structures.tar.gz\"\r\n"
        f"Content-Type: application/gzip\r\n\r\n"
    ).encode("utf-8") + binary_payload + f"\r\n--{boundary}--\r\n".encode("utf-8")

    headers = {
        "Content-Type": f"multipart/form-data; boundary={boundary}",
        "Content-Length": str(len(body)),
    }
    form = parse_form(io.BytesIO(body), headers, max_bytes=1024 * 1024)
    assert form.getfirst("csrf") == "secret-token-123"
    assert form.getfirst("code") == "bc_strain"
    assert form.getfirst("missing", "fallback") == "fallback"
    assert "pdb_tar" in form
    assert form["pdb_tar"].filename == "structures.tar.gz"
    assert form["pdb_tar"].file.read() == binary_payload


def test_parse_urlencoded_and_enforce_size_limit():
    body = b"csrf=tok&id=job-123"
    headers = {
        "Content-Type": "application/x-www-form-urlencoded",
        "Content-Length": str(len(body)),
    }
    form = parse_form(io.BytesIO(body), headers, max_bytes=1024)
    assert form.getfirst("csrf") == "tok"
    assert form.getfirst("id") == "job-123"

    with pytest.raises(ValueError, match="exceeds limit"):
        parse_form(io.BytesIO(body), headers, max_bytes=5)


def test_portal_does_not_import_deprecated_cgi_module():
    portal_src = (ROOT / "portal" / "suss_portal.py").read_text(encoding="utf-8")
    assert "import cgi" not in portal_src
    assert "FieldStorage" not in portal_src
