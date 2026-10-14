"""Tests of `get-genes`: BSD checksum, CHECKSUMS, release lookup, download.

No network is used.
"""
from pathlib import Path
import random
import re
import shutil
import subprocess

import pytest

from dmrdetector.steps import annotate


def test_bsd_sum_known_and_matches_loop(tmp_path):
    f = tmp_path / "a"
    f.write_bytes(b"a")
    assert annotate.bsd_sum(str(f)) == (97, 1) == annotate._bsd_sum_py(
        str(f))  # rotate(0) + 97
    rng = random.Random(1)
    for n in (0, 1, 1023, 1024, 1025, 3 * (1 << 20) + 7):
        f.write_bytes(rng.randbytes(n))
        want = annotate._bsd_sum_py(str(f))
        assert want[1] == -(-n // 1024)
        assert annotate.bsd_sum(str(f)) == want, n
        if shutil.which("sum"):
            got = subprocess.run(["sum", "-r", str(f)],
                                 capture_output=True,
                                 text=True).stdout.split()
            assert (int(got[0]), int(got[1])) == want, n


def test_bsd_sum_falls_back_without_sum_binary(tmp_path, monkeypatch):
    f = tmp_path / "b"
    f.write_bytes(b"abc" * 1000)
    want = annotate.bsd_sum(str(f))
    monkeypatch.setattr(annotate.shutil, "which", lambda _n: None)
    assert annotate.bsd_sum(str(f)) == want


CHECKSUMS = ("12345    138 Homo_sapiens.GRCh38.115.gtf.gz\n"
             "54321   99 Homo_sapiens.GRCh38.115.chr.gtf.gz\n"
             "00042 7 Homo_sapiens.GRCh38.115.abinitio.gtf.gz\n")


def test_parse_checksums_ensembl_format():
    assert annotate.parse_checksums(CHECKSUMS,
                                    "Homo_sapiens.GRCh38.115.gtf.gz") == (12345,
                                                                          138)
    assert annotate.parse_checksums(
        CHECKSUMS, "Homo_sapiens.GRCh38.115.abinitio.gtf.gz") == (42, 7)
    assert annotate.parse_checksums(CHECKSUMS,
                                    "Homo_sapiens.GRCh38.116.gtf.gz") is None


def fake_ftp(monkeypatch, listing: list[int], has_gtf: set[int]):

    def _get(url):
        if url.endswith("/pub/"):
            return "".join(
                f'<a href="release-{n}/">release-{n}/</a>\n' for n in listing)
        n = int(re.search(r"release-(\d+)/", url).group(1))
        if n not in has_gtf and n not in listing:
            raise OSError("404")
        return CHECKSUMS.replace("115",
                                 str(n)) if n in has_gtf else "1 1 other.gz\n"

    monkeypatch.setattr(annotate, "_get", _get)
    monkeypatch.setattr(annotate, "_ssl_fallback", lambda: None)


def test_latest_release_picks_newest_with_gtf_and_falls_back(monkeypatch):
    fake_ftp(monkeypatch, [114, 115, 116, 117], {114, 115, 116, 117})
    assert annotate.latest_release() == 117
    fake_ftp(monkeypatch, [114, 115, 116, 117],
             {115, 116})  # 117 half-populated
    assert annotate.latest_release() == 116
    fake_ftp(monkeypatch, [114, 115, 116, 117],
             {115})  # two newest half-populated
    assert annotate.latest_release() == 115
    fake_ftp(monkeypatch, [114, 115, 116, 117],
             {114})  # more than two older releases are not tried
    with pytest.raises(SystemExit):
        annotate.latest_release()


def test_get_genes_latest_downloads_and_verifies(tmp_path, monkeypatch, capsys):
    gtf = b"#!genome-build GRCh38\n" * 50
    f = tmp_path / "ref.gz"
    f.write_bytes(gtf)
    s, blocks = annotate.bsd_sum(str(f))
    ck = f"{s:05d} {blocks} Homo_sapiens.GRCh38.116.gtf.gz\n"
    monkeypatch.setattr(annotate, "_ssl_fallback", lambda: None)
    monkeypatch.setattr(annotate, "latest_release", lambda: 116)

    def _download(url, dest):
        Path(dest).write_text(ck) if url.endswith("CHECKSUMS") else Path(
            dest).write_bytes(gtf)

    monkeypatch.setattr(annotate, "_download", _download)
    out = annotate.get_genes("latest", str(tmp_path / "ens"))
    assert out == str(
        tmp_path / "ens" /
        "Homo_sapiens.GRCh38.116.gtf.gz") and Path(out).read_bytes() == gtf
    monkeypatch.setattr(
        annotate, "_download", lambda u, d: Path(d).write_text(ck)
        if u.endswith("CHECKSUMS") else Path(d).write_bytes(b"bad"))
    (tmp_path / "ens" / "Homo_sapiens.GRCh38.116.gtf.gz").unlink()
    with pytest.raises(SystemExit, match="checksum mismatch"):
        annotate.get_genes(116, str(tmp_path / "ens"))


def test_get_genes_cli_release_argument():
    p = __import__("dmrdetector.cli", fromlist=["x"]).build_parser()
    assert p.parse_args(["get-genes", "--release", "latest", "--outdir",
                         "d"]).release == "latest"
    assert p.parse_args(["get-genes", "--release", "86", "--outdir",
                         "d"]).release == "86"
    with pytest.raises(SystemExit):
        p.parse_args(["get-genes", "--release", "newest", "--outdir", "d"])
