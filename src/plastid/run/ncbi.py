"""NCBI E-utilities client — the only module in the package that talks to the network.

Reference Panel policy (PROJECT.md §10.8): search and summaries build a *dated cache*; runs
never query NCBI. Without an API key NCBI allows 3 requests/second, so calls are spaced.
Tests inject a transport callable — no test ever touches the network.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

__all__ = ["EUTILS_BASE", "NcbiClient", "NcbiError", "http_transport"]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MIN_INTERVAL_S = 0.34  # 3 requests/second without an API key
USER_AGENT = "plastid-suite/0.1 (academic research; contact: project maintainers)"


class NcbiError(Exception):
    """NCBI returned an error, or something this client cannot parse."""


Transport = Callable[[str], str]


def http_transport(url: str) -> str:
    """Default transport: one HTTPS GET."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310 — https only
        return response.read().decode("utf-8", "replace")


@dataclass
class NcbiClient:
    """Thin, testable wrapper over esearch/esummary/efetch."""

    transport: Transport = http_transport
    api_key: str | None = None
    interval_s: float = MIN_INTERVAL_S
    sleep: Callable[[float], None] = time.sleep
    last_call: list[float] = field(default_factory=lambda: [0.0])

    def _get(self, endpoint: str, params: dict[str, str]) -> str:
        if self.api_key:
            params = {**params, "api_key": self.api_key}
        url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
        wait = self.interval_s - (time.monotonic() - self.last_call[0])
        if wait > 0:
            self.sleep(wait)
        text = self.transport(url)
        self.last_call[0] = time.monotonic()
        return text

    def search_refseq_plastomes(self, taxon: str, *, retmax: int = 200) -> list[str]:
        """UIDs of RefSeq chloroplast genomes for one taxon.

        Two measured traps avoided (2026-10-03, Iridaceae):
        * no title substring filter — "complete genome[title]" drops valid records;
        * **both** organelle filters — ``chloroplast[filter]`` alone returns 71 while
          ``plastid[filter]`` returns 76, because records worded "plastid, complete genome"
          (e.g. NC_024936.1, our own comparison genome) do not match "chloroplast".
        Completeness is enforced structurally, by length, at proposal time.
        """
        term = f"{taxon}[Organism] AND refseq[filter] AND (chloroplast[filter] OR plastid[filter])"
        payload = json.loads(
            self._get(
                "esearch.fcgi",
                {"db": "nuccore", "retmode": "json", "retmax": str(retmax), "term": term},
            )
        )
        ids = payload.get("esearchresult", {}).get("idlist")
        if ids is None:
            raise NcbiError(f"esearch returned no idlist for {taxon!r}: {payload}")
        return [str(uid) for uid in ids]

    def summaries(self, ids: Sequence[str]) -> list[dict]:
        """esummary records (accession, organism, length, title, taxid) for nuccore UIDs."""
        if not ids:
            return []
        payload = json.loads(
            self._get(
                "esummary.fcgi",
                {"db": "nuccore", "retmode": "json", "id": ",".join(ids)},
            )
        )
        result = payload.get("result", {})
        return [result[uid] for uid in result.get("uids", []) if uid in result]

    def lineages(self, taxids: Iterable[str]) -> dict[str, str]:
        """taxid → lineage string ("Eukaryota; …; Iridaceae; Iris")."""
        unique = sorted({str(taxid) for taxid in taxids if str(taxid)})
        if not unique:
            return {}
        xml_text = self._get(
            "efetch.fcgi", {"db": "taxonomy", "retmode": "xml", "id": ",".join(unique)}
        )
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as error:
            raise NcbiError(f"taxonomy efetch returned unparseable XML: {error}") from error
        return {
            taxon.findtext("TaxId") or "": (taxon.findtext("Lineage") or "")
            for taxon in root.findall("Taxon")
        }

    def fetch_genbank(self, accession: str) -> str:
        """One GenBank record with its feature table (gbwithparts)."""
        return self._get(
            "efetch.fcgi",
            {"db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "text"},
        )
