"""Dermo and MSX by disease bar and year from the Maryland DNR Fall Oyster Survey (PLAN.md P1.7).

Source: Maryland DNR, "The 2024 Fall Survey" report,
https://dnr.maryland.gov/fisheries/Documents/oysters/2024RptFinal.pdf
- Table 3 (pp. 62-68): Perkinsus marinus prevalence (%) and mean intensity
  (scale 0-7) at the 43 disease monitoring bars, 1990-2024.
- Table 4 (pp. 69-72): Haplosporidium nelsoni (MSX) prevalence (%), 1990-2024.
- Table 5 (pp. 73-77): Observed Mortality Index (%), 1985-2024, "an indicator of
  annual natural mortality ... calculated as the total number of dead oysters
  (boxes and gapers) divided by the sum of live and dead oysters" (report p. 13);
  written to mortality_obs with period_days = 365 and method "box_count".

Rows are listed in the same fixed bar order on every page, so the k-th data row
on a page is the k-th bar of BAR_ORDER; each row is checked against that name.
Two bars share the name "Ragged Point" (Little Choptank River, then Potomac
River) and are told apart by position.

Encoding choices (also in PLAN.md):
- NA (too few oysters) and ND (collected, not diagnosed; the report assumes 0)
  are not written: only measured values become rows.
- The tables give no per-bar sample size, so ``n`` is null.
- The tables give only the survey year. The 2024 survey ran 8 October to
  25 November (report p. 5), midpoint 1 November; every row is dated 1 November
  of its year and ``source`` says the date is year-resolution only.
- Bar coordinates are not in the report. They come from Smith, G.F. & Jordan, S.J.
  (1992) Monitoring Maryland's Chesapeake Bay Oysters: A Comprehensive
  Characterization of Modified Fall Survey Results 1990-1991. MD DNR report
  CBRM-OX-93-3, Table 1 (pp. 6-7), which gives the fixed sampling position of each
  of these 43 disease bars (uncorrected LORAN-C, degrees-minutes-seconds;
  https://www.govinfo.gov/content/pkg/CZIC-sh365-m3-s55-1992/pdf/CZIC-sh365-m3-s55-1992.pdf).
  Piney Island East is "Piney Island" (TSPI) in 1992; its position lies nearer the
  Piney Island West polygon of the MD Historic Oyster Bottom layer, recorded in
  its site notes.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from datetime import date
from pathlib import Path

import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table
from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

URL = "https://dnr.maryland.gov/fisheries/Documents/oysters/2024RptFinal.pdf"
RAW_FILE = "mddnr_fall_survey/2024RptFinal.pdf"
DERMO_PAGES = range(62, 69)   # 1-based report pages holding Table 3
MSX_PAGES = range(69, 73)     # Table 4
MORTALITY_PAGES = range(73, 78)   # Table 5
INTENSITY_SCALE = "mddnr_0_7"

BAR_ORDER = [
    "Swan Point", "Hackett Point", "Holland Point (S)", "Stone Rock", "Flag Pond (S)", "Hog Island",
    "Butler", "Buoy Rock", "Old Field (S)", "Bugby", "Parsons Island", "Hollicutt Noose",
    "Bruffs Island (S)", "Turtle Back", "Long Point (S)", "Cook Point (S)", "Royston", "Lighthouse",
    "Sandy Hill (S)", "Oyster Shell Pt. (S)", "Tilghman Wharf", "Deep Neck", "Double Mills (S)",
    "Cason (S)", "Ragged Point", "Norman Addition", "Goose Creek", "Wilson Shoals (S)", "Georges (S)",
    "Holland Straits", "Sharkfin Shoal", "Back Cove", "Piney Island East", "Old Woman's Leg",
    "Marumsco", "Broome Island", "Chicken Cock", "Pagan (S)", "Lancaster", "Mills West",
    "Cornfield Harbor", "Ragged Point", "Lower Cedar Point",
]
SITE_SUFFIX = {24: "little_choptank", 41: "potomac"}   # the two "Ragged Point" bars (0-based index)

BAR_COORDS_DDMMSS = [
    ('BNSP', '390827', '761810'),
    ('UBHA', '385859', '762500'),
    ('WSHP', '384407', '763008'),
    ('MESR', '383920', '762259'),
    ('WSFP', '382606', '762609'),
    ('WSHI', '381854', '762301'),
    ('WSBU', '380632', '761937'),
    ('CHBR', '385938', '761242'),
    ('CHOF', '390448', '760952'),
    ('EBBU', '385255', '761320'),
    ('EBPI', '385420', '761602'),
    ('EBHN', '385114', '762106'),
    ('MRBI', '385129', '761135'),
    ('MRTU', '385119', '761421'),
    ('MRLP', '384613', '761032'),
    ('CRCP', '383909', '761725'),
    ('CRRO', '384115', '761430'),
    ('CRLI', '383927', '761122'),
    ('CRSH', '383539', '760700'),
    ('CROS', '383518', '760001'),
    ('CRTW', '384247', '761915'),
    ('BCDN', '384417', '761433'),
    ('TADM', '384347', '760825'),
    ('LCCA', '383159', '761421'),
    ('LCRP', '383218', '761750'),
    ('HRNO', '381519', '760815'),
    ('FBGC', '381702', '760130'),
    ('NRWS', '381735', '755518'),
    ('MAGE', '380727', '755124'),
    ('HOHO', '380644', '760430'),
    ('TSSS', '381256', '755929'),
    ('TSBC', '380225', '755939'),
    ('TSPI', '380409', '755734'),
    ('TSOW', '375747', '755823'),
    ('PSMA', '375733', '754409'),
    ('PXBI', '382428', '763351'),
    ('SMCC', '380723', '762613'),
    ('SMPA', '381130', '762635'),
    ('WWLA', '381635', '764945'),
    ('WWMW', '382009', '765129'),
    ('PRCH', '380253', '762001'),
    ('PRRP', '380922', '763633'),
    ('PRLC', '381959', '765850'),
]


def _ddmmss(text: str) -> float:
    d, m, sec = int(text[:-4]), int(text[-4:-2]), int(text[-2:])
    return d + m / 60 + sec / 3600


def site_rows() -> list[dict[str, object]]:
    """Registry rows (sites contract) for the 43 disease bars."""
    rows = []
    for k, (code, lat, lon) in enumerate(BAR_COORDS_DDMMSS):
        note = (f"MD DNR Fall Survey disease bar {BAR_ORDER[k]}; Oxford Lab code {code}; position from "
                "Smith & Jordan 1992 (CBRM-OX-93-3) Table 1, uncorrected LORAN-C.")
        if code == "TSPI":
            note += " Listed as Piney Island in 1992; position is nearer the Piney Island West polygon."
        rows.append({"site_id": site_id(k), "name": BAR_ORDER[k].replace("(S)", "(sanctuary)").strip(),
                     "lat": round(_ddmmss(lat), 6), "lon": -round(_ddmmss(lon), 6), "region": "chesapeake",
                     "growing_area": None, "notes": note})
    return rows

_VALUE = r"(?:\d+(?:\.\d+)?|NA|ND)"


def site_id(index: int) -> str:
    name = BAR_ORDER[index].replace("(S)", "").replace("'", "")
    slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if index in SITE_SUFFIX:
        slug += "_" + SITE_SUFFIX[index]
    return "mddnr_" + slug


def fetch(raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> Path:
    body = get(URL)
    if not body.startswith(b"%PDF"):
        raise RuntimeError("Fall Survey download is not a PDF")
    return save_raw(body, raw_dir / RAW_FILE, overwrite=overwrite)


def _normalise(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower().replace("�", "'"))


def parse_page(text: str, values_per_year: int) -> list[tuple[int, list[int], list[str]]]:
    """Return [(bar_index, years, value tokens)] for one table page."""
    years: list[int] = []
    rows: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        tokens = line.split()
        found = [int(t) for t in tokens if re.fullmatch(r"(19|20)\d\d", t)]
        if not years and len(found) >= 1 and len(found) == len([t for t in tokens if re.fullmatch(r"\d{4}", t)]):
            if re.search(r"Oyster Bar|^\s*(19|20)\d\d", line):
                years = found
                continue
        if not years:
            continue
        m = re.match(rf"^(.*?)((?:\s+{_VALUE})+)\s*$", line)
        if not m or not re.search(r"[A-Za-z]", m.group(1)):
            continue
        values = m.group(2).split()
        rows.append((m.group(1).strip(), values))
    out = []
    data_rows = [(label, vals) for label, vals in rows
                 if not re.match(r"(?i)(annual means|average|frequency|bar freq)", label)]
    for k, (label, values) in enumerate(data_rows[:len(BAR_ORDER)]):
        expected = _normalise(BAR_ORDER[k])
        if not _normalise(label).endswith(expected):
            raise ValueError(f"row {k} label {label!r} does not end with bar {BAR_ORDER[k]!r}")
        n_needed = values_per_year * len(years)
        if len(values) < n_needed:
            raise ValueError(f"row {label!r}: {len(values)} values for {len(years)} years")
        out.append((k, years, values[:n_needed]))
    if len(out) != len(BAR_ORDER):
        raise ValueError(f"expected {len(BAR_ORDER)} bar rows, parsed {len(out)}")
    return out


def _num(token: str) -> float | None:
    return None if token in ("NA", "ND") else float(token)


def rows_from_pages(dermo_pages: list[str], msx_pages: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    src = "mddnr_fall_survey:2024RptFinal:{table};date=survey_year_only"
    for text in dermo_pages:
        for k, years, vals in parse_page(text, 2):
            for i, y in enumerate(years):
                prev, inten = _num(vals[2 * i]), _num(vals[2 * i + 1])
                if prev is None:
                    continue
                rows.append({"site_id": site_id(k), "date": date(y, 11, 1), "pathogen": "P_marinus",
                             "prevalence": prev / 100.0, "intensity_mean": inten,
                             "intensity_scale": INTENSITY_SCALE, "n": None, "source": src.format(table="Table3")})
    for text in msx_pages:
        for k, years, vals in parse_page(text, 1):
            for i, y in enumerate(years):
                prev = _num(vals[i])
                if prev is None:
                    continue
                rows.append({"site_id": site_id(k), "date": date(y, 11, 1), "pathogen": "H_nelsoni",
                             "prevalence": prev / 100.0, "intensity_mean": None,
                             "intensity_scale": "none", "n": None, "source": src.format(table="Table4")})
    return rows


def disease_obs(dermo_pages: list[str], msx_pages: list[str]) -> pa.Table:
    rows = rows_from_pages(dermo_pages, msx_pages)
    rows.sort(key=lambda r: (r["pathogen"], r["site_id"], r["date"]))
    table = pa.Table.from_pylist(rows, schema=arrow_schema("disease_obs"))
    validate_table(table, SPECS["disease_obs"])
    return table


def mortality_obs(mortality_pages: list[str]) -> pa.Table:
    """Contract mortality_obs rows from Table 5 page texts (NA rows are not written)."""
    rows = []
    for text in mortality_pages:
        for k, years, vals in parse_page(text, 1):
            for i, y in enumerate(years):
                v = _num(vals[i])
                if v is None:
                    continue
                rows.append({"site_id": site_id(k), "date": date(y, 11, 1), "period_days": 365,
                             "mortality_fraction": v / 100.0, "method": "box_count",
                             "source": "mddnr_fall_survey:2024RptFinal:Table5;date=survey_year_only"})
    rows.sort(key=lambda r: (r["site_id"], r["date"]))
    table = pa.Table.from_pylist(rows, schema=arrow_schema("mortality_obs"))
    validate_table(table, SPECS["mortality_obs"])
    return table


def load_mortality_pages(raw_dir: Path = RAW_DIR) -> list[str]:
    import pdfplumber

    with pdfplumber.open(str(raw_dir / RAW_FILE)) as pdf:
        return [pdf.pages[p - 1].extract_text() or "" for p in MORTALITY_PAGES]


def load_pages(raw_dir: Path = RAW_DIR) -> tuple[list[str], list[str]]:
    import pdfplumber

    with pdfplumber.open(str(raw_dir / RAW_FILE)) as pdf:
        dermo = [pdf.pages[p - 1].extract_text() or "" for p in DERMO_PAGES]
        msx = [pdf.pages[p - 1].extract_text() or "" for p in MSX_PAGES]
    return dermo, msx
