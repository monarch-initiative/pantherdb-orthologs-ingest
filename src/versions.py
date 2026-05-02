"""Upstream source version fetcher for pantherdb-orthologs-ingest.

Two logical sources:
  - infores:panther — the AllOrthologs.tar.gz from PantherDB. The
    `current_release/` symlink doesn't carry an explicit version, so use
    HTTP Last-Modified.
  - infores:ncbi-gene — gene_info.gz used for ID mapping. HTTP Last-Modified.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from kozahub_metadata_schema import (
    now_iso,
    urls_from_download_yaml,
    version_from_http_last_modified,
)


INGEST_DIR = Path(__file__).resolve().parents[1]
DOWNLOAD_YAML = INGEST_DIR / "download.yaml"


def get_source_versions() -> list[dict[str, Any]]:
    panther_urls = urls_from_download_yaml(DOWNLOAD_YAML, contains=["pantherdb.org"])
    ncbi_urls = urls_from_download_yaml(DOWNLOAD_YAML, contains=["ftp.ncbi.nlm.nih.gov"])
    now = now_iso()

    sources: list[dict[str, Any]] = []

    if panther_urls:
        ver, method = version_from_http_last_modified(panther_urls[0])
        sources.append({
            "id": "infores:panther",
            "name": "PantherDB Orthologs",
            "urls": panther_urls,
            "version": ver,
            "version_method": method,
            "retrieved_at": now,
        })

    if ncbi_urls:
        ver, method = version_from_http_last_modified(ncbi_urls[0])
        sources.append({
            "id": "infores:ncbi-gene",
            "name": "NCBI Gene",
            "urls": ncbi_urls,
            "version": ver,
            "version_method": method,
            "retrieved_at": now,
        })

    return sources
