"""NCBI client — URL construction, parsing and rate limiting, all offline (transport injected)."""

from __future__ import annotations

import json
import urllib.parse

import pytest

from plastid.core.panel import entry_from_summary
from plastid.run.ncbi import EUTILS_BASE, NcbiClient, NcbiError

LINEAGE = (
    "cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Embryophyta; "
    "Magnoliopsida; Asparagales; Iridaceae; Iris"
)


def _urls_of(transport_calls: list[str]) -> list[dict[str, list[str]]]:
    return [urllib.parse.parse_qs(urllib.parse.urlparse(url).query) for url in transport_calls]


def test_search_builds_the_query_and_parses_ids() -> None:
    calls: list[str] = []

    def transport(url: str) -> str:
        calls.append(url)
        return json.dumps({"esearchresult": {"idlist": ["111", "222"]}})

    client = NcbiClient(transport=transport, sleep=lambda _: None)
    assert client.search_refseq_plastomes("Iridaceae", retmax=5) == ["111", "222"]

    assert calls[0].startswith(f"{EUTILS_BASE}/esearch.fcgi?")
    params = _urls_of(calls)[0]
    assert params["db"] == ["nuccore"]
    assert params["retmax"] == ["5"]
    term = params["term"][0]
    assert "Iridaceae[Organism]" in term
    assert "refseq[filter]" in term
    # both organelle filters: chloroplast[filter] alone misses "plastid"-worded records
    assert "chloroplast[filter]" in term and "plastid[filter]" in term


def test_search_without_an_idlist_raises() -> None:
    client = NcbiClient(
        transport=lambda url: json.dumps({"esearchresult": {}}), sleep=lambda _: None
    )
    with pytest.raises(NcbiError, match="idlist"):
        client.search_refseq_plastomes("Nope")


def test_summaries_returns_records_in_uid_order() -> None:
    payload = {
        "result": {
            "uids": ["111", "222"],
            "111": {
                "accessionversion": "NC_050833.1",
                "organism": "Iris domestica",
                "slen": 153724,
            },
            "222": {"accessionversion": "NC_024936.1", "organism": "Iris gatesii", "slen": 153441},
        }
    }
    client = NcbiClient(transport=lambda url: json.dumps(payload), sleep=lambda _: None)
    records = client.summaries(["111", "222"])
    assert [record["accessionversion"] for record in records] == ["NC_050833.1", "NC_024936.1"]
    assert client.summaries([]) == []  # no request for an empty id list


def test_lineages_parses_taxonomy_xml() -> None:
    xml_text = (
        "<TaxaSet>"
        "<Taxon><TaxId>12345</TaxId><Lineage>" + LINEAGE + "</Lineage></Taxon>"
        "<Taxon><TaxId>67890</TaxId><Lineage>A; B; Iridaceae; Crocus</Lineage></Taxon>"
        "</TaxaSet>"
    )
    client = NcbiClient(transport=lambda url: xml_text, sleep=lambda _: None)
    lineages = client.lineages(["12345", "67890", "12345"])
    assert lineages["12345"].endswith("Iridaceae; Iris")
    assert lineages["67890"].endswith("Iridaceae; Crocus")


def test_unparseable_taxonomy_xml_raises() -> None:
    client = NcbiClient(transport=lambda url: "<not-xml", sleep=lambda _: None)
    with pytest.raises(NcbiError, match="XML"):
        client.lineages(["1"])


def test_calls_are_spaced_to_respect_the_rate_limit() -> None:
    slept: list[float] = []
    client = NcbiClient(
        transport=lambda url: json.dumps({"esearchresult": {"idlist": []}}),
        sleep=slept.append,
        interval_s=0.34,
    )
    client.search_refseq_plastomes("A")
    client.search_refseq_plastomes("B")
    assert len(slept) == 1  # only the second call waits
    assert 0 < slept[0] <= 0.34


def test_entry_from_summary_maps_family_genus_order() -> None:
    entry = entry_from_summary(
        {"accessionversion": "NC_050833.1", "organism": "Iris domestica", "slen": 153724},
        LINEAGE,
    )
    assert entry.family == "Iridaceae"
    assert entry.order == "Asparagales"
    assert entry.genus == "Iris"
    assert entry.species == "Iris domestica"
    assert entry.length == 153_724
    assert entry.refseq is True


def test_entry_from_summary_needs_an_accession() -> None:
    with pytest.raises(ValueError, match="accession"):
        entry_from_summary({"organism": "Iris"}, LINEAGE)


def test_resolve_reads_one_accession_without_a_search() -> None:
    """A named accession costs one esummary call plus its lineage — never an esearch."""
    calls: list[str] = []

    def transport(url: str) -> str:
        calls.append(url)
        if "esummary.fcgi" in url:
            return json.dumps(
                {
                    "result": {
                        "uids": ["3202415722"],
                        "3202415722": {
                            "accessionversion": "PX923121.1",
                            "organism": "Coleus amboinicus",
                            "taxid": 204180,
                            "slen": 152544,
                        },
                    }
                }
            )
        return (
            "<TaxaSet><Taxon><TaxId>204180</TaxId>"
            "<Lineage>Eukaryota; Lamiaceae; Coleus</Lineage></Taxon></TaxaSet>"
        )

    client = NcbiClient(transport=transport, sleep=lambda _: None)
    summary, lineage = client.resolve("PX923121.1")

    assert summary["accessionversion"] == "PX923121.1"
    assert lineage.endswith("Coleus")
    endpoints = [urllib.parse.urlparse(url).path.rsplit("/", 1)[-1] for url in calls]
    assert endpoints == ["esummary.fcgi", "efetch.fcgi"]
    assert _urls_of(calls)[0]["id"] == ["PX923121.1"]


def test_resolve_names_an_accession_ncbi_does_not_have() -> None:
    client = NcbiClient(
        transport=lambda url: json.dumps({"result": {"uids": []}}), sleep=lambda _: None
    )
    with pytest.raises(NcbiError, match="XX_000000.1"):
        client.resolve("XX_000000.1")


def test_resolve_reports_a_uid_without_a_record_as_ncbi_error() -> None:
    # esummary can list a uid the result block does not carry; that is an NCBI problem,
    # not a KeyError leaking out of the client.
    payload = json.dumps({"result": {"uids": ["123"]}})
    client = NcbiClient(transport=lambda url: payload, sleep=lambda _: None)
    with pytest.raises(NcbiError, match="XX_000000.1"):
        client.resolve("XX_000000.1")
