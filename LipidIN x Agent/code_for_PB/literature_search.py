"""Keyless online literature retrieval for the 'similar reports' verification check.

用于"是否有相似报道"验证的免密钥在线文献检索。

Uses public REST APIs that need no API key: Europe PMC (primary, returns abstracts) and CrossRef
(fallback). Network failures degrade gracefully to an empty list so verification never crashes.

使用无需密钥的公共 REST 接口：Europe PMC（主，返回摘要）与 CrossRef（备用）。网络失败时优雅降级为
空列表，保证验证流程不崩。
"""

import json
import urllib.parse
import urllib.request
import urllib.error


EUROPE_PMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF_URL = "https://api.crossref.org/works"
USER_AGENT = "LipidIN-agent/1.0 (research use)"


def _get_json(url, timeout=25):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def search_europepmc(query, limit=5, timeout=25):
    params = {
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": str(limit),
        "sort": "CITED desc",
    }
    url = f"{EUROPE_PMC_URL}?{urllib.parse.urlencode(params)}"
    data = _get_json(url, timeout=timeout)
    results = []
    for item in data.get("resultList", {}).get("result", []):
        results.append({
            "title": item.get("title", ""),
            "authors": item.get("authorString", ""),
            "year": item.get("pubYear", ""),
            "journal": (item.get("journalInfo", {}) or {}).get("journal", {}).get("title", ""),
            "doi": item.get("doi", ""),
            "abstract": (item.get("abstractText", "") or "")[:600],
            "cited_by": item.get("citedByCount", 0),
            "source": "EuropePMC",
        })
    return results


def search_crossref(query, limit=5, timeout=25):
    params = {"query": query, "rows": str(limit), "sort": "relevance"}
    url = f"{CROSSREF_URL}?{urllib.parse.urlencode(params)}"
    data = _get_json(url, timeout=timeout)
    results = []
    for item in data.get("message", {}).get("items", []):
        title = item.get("title", [""])
        results.append({
            "title": title[0] if title else "",
            "authors": ", ".join(
                f"{a.get('given', '')} {a.get('family', '')}".strip()
                for a in item.get("author", [])[:6]
            ),
            "year": str((item.get("issued", {}).get("date-parts", [[None]]) or [[None]])[0][0] or ""),
            "journal": (item.get("container-title", [""]) or [""])[0],
            "doi": item.get("DOI", ""),
            "abstract": "",
            "cited_by": item.get("is-referenced-by-count", 0),
            "source": "CrossRef",
        })
    return results


def find_similar_reports(query, limit=5):
    """Search Europe PMC first, fall back to CrossRef, and never raise on network failure.

    先查 Europe PMC，失败则退回 CrossRef；网络异常一律不抛错。
    """
    for search in (search_europepmc, search_crossref):
        try:
            results = search(query, limit=limit)
            if results:
                return results
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError, ValueError):
            continue
    return []


def format_reports_for_prompt(reports, max_items=5):
    """Render retrieved references into a compact block for an LLM verification prompt.

    把检索到的文献整理成紧凑文本块，供大模型验证提示使用。
    """
    if not reports:
        return "No references retrieved. / 未检索到文献。"
    lines = []
    for i, ref in enumerate(reports[:max_items], start=1):
        lines.append(
            f"[{i}] {ref.get('title', '')} ({ref.get('year', '')}, {ref.get('journal', '')}; "
            f"cited {ref.get('cited_by', 0)}; DOI {ref.get('doi', '')})"
        )
        abstract = ref.get("abstract", "")
        if abstract:
            lines.append(f"    Abstract: {abstract}")
    return "\n".join(lines)


if __name__ == "__main__":
    hits = find_similar_reports("phosphatidylcholine remodeling in hepatocellular carcinoma", limit=3)
    print(format_reports_for_prompt(hits))
