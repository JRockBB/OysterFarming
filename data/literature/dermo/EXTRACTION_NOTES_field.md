# Dermo NE field data search notes (2026-10-14)

## NOAA Mussel Watch (NS&T) oyster histopathology  [DOWNLOADED]
- URL: https://cdn.coastalscience.noaa.gov/datasets/nsandt/dataFiles/mw/eastCoast/MW_EastCoast_Oyster_Histopath.tsv (+ Sites_MusselWatch_East_Coast.tsv, Histopath_Oyster_Data_Dictionary.pdf) -> dl/mw/
- Machine-readable TSV, per-oyster records; columns incl. NST_Site, lat/lon, Fiscal Year, Dermo (Y/N), Dermo Infection Intensity (Negative..Heavy), Dermo Numerical Value 0-5 (Mackin-style 0,0.33,...5 in thirds), Dermo Description.
- East Coast years 1988-2010; 3638 oysters. n=12/site/yr (60 in 1988).
- NORTHERNMOST OYSTER SITES = Delaware Bay (39.38N). NO oyster sites in NJ coast north/NY/CT/RI/MA/NH/ME in this file (New England Mussel Watch sites are Mytilus; mussel histopath file has no Dermo).
- Summary written: mw_oyster_dermo_by_site_year.csv (site-year prevalence, mean intensity).
- No co-located T/S in this file. -> Not useful for cool+salty criterion.

## MAINE
### Marquis et al. 2020 Pathogens 9:256 (Bigelow Lab; Fernandez Robledo) [DOWNLOADED fulltext XML + figs: dl/marquis2020/]
- doi:10.3390/pathogens9040256; PMC7238206. qPCR (P. marinus-specific primers, Marsh et al. 1995), pooled rectum/gill/mantle, NO RFTM/histology.
- 6 commercial leases: Bagaduce R (Brooksville), Weskeag R (Thomaston), Jack's Point (Newcastle, Damariscotta R), Prentiss Island (Bristol? - Damariscotta R. estuary), New Meadows R (Bath), Webhannet R (Wells). No coordinates in text (map only).
  [NB: WebFetch summary wrongly placed Bagaduce/Weskeag in DRE; only Jack's Pt & Prentiss Is. are DRE.]
- Monthly Jun-Oct 2016 & 2017, n=24/site/month, 1440 oysters.
- P. marinus qPCR-prevalence overall 57% assays positive (45% 2016, 68% 2017); site means 7-39% (Webhannet 7%, Jack's Pt 39%, Prentiss 34%, New Meadows 33%, Weskeag 37%, Bagaduce 32%).
- Fig 3 (dl/marquis2020/supp/pathogens-09-00256-g003.jpg): monthly bars, intensity log10 copies/ug DNA ranges ~ -3 to +2 (i.e. 0.001-100 copies/ug) = extremely light; individual months reach 100% prev (Bagaduce 7/2016, Weskeag 7/2016, Jack's Pt 6/2016, 6/2017 ~90%). Data only in figure -> digitization needed; no data table.
- Fig 7: LOBO buoy near Prentiss Is. 2016/2017 daily T: peaks ~22-23 C (2016) and ~21-22 C (2017) late Jul-Aug; T>20 C roughly mid-Jul to early Sep 2016; Salinity 28-32.3 PSU Jun-Oct. => EXACT cool+salty regime (S ~30, DD>20 small).
- No die-offs; no clinical Dermo; authors: parasites "in Maine to stay". Notes Gulf of Maine warming.
- Key Maine refs: Kleinschuster & Parent 1995 JSR 14:489 (first subclinical Perkinsus in ME); Marquis et al 2015 Parasitol Int 64:299 (PCR survey Gulf of Maine); Fernandez Robledo et al 2018 Aquaculture 493:9-17 (historical perspective, 15-65x increase P. marinus/P. chesapeaki by PCR 2002->2014); Pecher et al 2008 JSR (northern range PCR survey, P. chesapeaki 6.9% Walpole ME 2002).

### Fernandez Robledo et al. 2018 Aquaculture 493:9-17 "Pathogens of marine bivalves in Maine: historical perspective" [DOWNLOADED author ms: dl/robledo2018/ms.pdf, NOAA IR noaa/57739]
- 11 refs to Perkinsus in Maine. 1990s "sub-clinical" (Kleinschuster & Parent 1995, RFTM); 2000s low PCR prevalence (Pecher et al 2008 J Parasitol 94:410); Marquis et al 2015 Parasitol Int 64:299: DRE (Jones Cove) 2012 PCR prev P. marinus 65.2%, P. chesapeaki 43.5% = 65x / 15x increase over 12 yrs.
- RFTM cannot separate P. marinus vs P. chesapeaki -> early Maine RFTM records ambiguous.
- "So far, no oyster die-offs have been ascribed to Perkinsus spp. in Maine". No histopathology.
- Penobscot Bay buoy (NERACOOS, since 2001): salinity 17-35 psu; currently <1 day/yr with daily mean T > 20 C; +2C -> ~6 d/yr; +4C -> ~20 d/yr.
- No intensity data (Mackin) anywhere in Maine literature found so far.

## NEW HAMPSHIRE - Great Bay estuary (NH Fish & Game, "Testing of Great Bay Oysters for Two Protozoan Pathogens")  [KEY DATASET, TRANSCRIBED]
- Annual fall (Oct-Dec) samples, natural beds: Adams Point, Nannie Island, Oyster River, Woodman Point, Piscataqua River, Squamscott River, + Salmon Falls, Sturgeon Bed, Seal Rock, Fox Pt, Bellamy, Crommet Cr, Berrys Brook, Lamprey R, Little Bay A/B/C aquaculture (2013-14), experimental reefs.
- Method: RFTM (rectal+mantle) by Rutgers Haskin lab; Mackin categories 0.5,1,2,3,4,5 counts per sample (n=5-40, usually 20 or 10). Histology for MSX.
- 1996-2014 table transcribed from 2015 report (UNH PREP #356, article 1359, via Wayback id_ snapshot 20240226195037) -> greatbay_dermo_1996_2014.csv (100 site-dates; computed prevalence and weighted prevalence = sum(count*cat)/n). One all-zero row 11/04/00 Piscataqua dropped by parser (typo "110/4/00"). Source-table inconsistencies (reported vs count-derived prevalence): 2002-11-20 Salmon Falls native, 2004-11-18 Nannie, 2008-10-22 Nannie, 2008-10-27 Squamscott (counts sum 19 > n=10; n likely 20).
- Also downloaded 2008 report (a1110), 2006 report (a1175), 2018 Review of Oyster Data (a1417, abundance only). Later disease reports (2013 #255, 2019 #415 article 1416, #1490) are behind Cloudflare at scholars.unh.edu and not in Wayback -> NOT retrieved; would extend to ~2019+.
- Trajectory: Dermo first seen 1996 (Spinney Creek, ME side of Piscataqua); 1997-2003 near zero (light only, ~0-20%); rise from 2004; 2005-2014 prevalence 75-100% at most beds, WP ~1-3.3 (2014 Nannie WP 3.3). Little Bay aquaculture oysters (new) 0-10%.
- No co-located T/S in reports. Great Bay NERR SWMP stations (Great Bay buoy GRBGB, Squamscott GRBSQ, Lamprey GRBLR, Oyster River GRBOR) give T/S 1990s-present via CDMO.
- Barber, Langan & Howell 1997 J Parasitol 83:148 (MSX epizootic Piscataqua 1995).
- T/S near Adams Pt/Nannie Is: UNH Great Bay Buoy (ERDDAP erddap.sensors.ioos.us dataset ism-neracoos-unh_gbb, 43.0715N -70.8678W, 2022-05 to present) -> dl/gbb/gbb_raw.csv; scripts/dd20.py
  DD>20 (sum daily mean-20): 2022 191 (71 d>20; Tmax daily 25.4); 2023 229 (gappy; wet yr, JJAS S mean 19.9); 2025 106 (very gappy); 2026 187 (86 d>20).
  JJAS salinity mean 29.2 (2022), 27.8 (2025), 27.0 (2026); p10 ~19-26.
  => Great Bay proper is ~cool+salty: S ~27-29, DD20 ~190-230 => MEETS criterion (S>=25, DD<300) in typical years, and Dermo there reached 80-100% prevalence with WP 2-3.3 (2009-2014). (Caveat: T/S years 2022-26 vs Dermo 1996-2014; need SWMP GRBGBWQ history for overlap.)

## DELAWARE BAY (Rutgers HSRL) [DOWNLOADED dl/saw/]
- SAW reports https://hsrl.rutgers.edu/wp-content/uploads/... (SAW2019,2020,2022,2023 = /2022/12/ or /2023/12/SAWyyyy.pdf; SAW2024 /2024/04/SAW2024.pdf; 2025 /2025/05/Documents_SAW2025.pdf). Downloaded SAW2024, Documents_SAW2025, StockReport1994.
- Content: fall stock survey; Dermo by RFTM, Mackin 0-5, "weighted prevalence (WP) = mean intensity of all oysters, 0-5"; 1.5 WP = mortality threshold. Recent SAWs give WP only by management REGION in figures (VLM, LM, MMT, MMM, SR, HM; "dermo era" since 1990); bed-level WP not tabulated in text of SAW2025 (figures / appendices D-H = abundance). 2024: MMT & MMM WP > 1.5 first time since 2017; SR & HM > 1.5; VLM nearly undetectable.
- 1990s "Stock Report" (StockReport1994/1996/1998.pdf) have bed-level tables: % prevalence and WP for 3 yrs (e.g. Shell Rock 1994/93/92 prev 90/100/98, WP 2.9/3.0/3.2) but scanned OCR -> needs manual transcription.
- Machine-readable alternative: Bushek et al. 2012 J Mar Res 70:225 "Long-term patterns of an estuarine pathogen along a salinity gradient" (Delaware Bay Dermo 1990-2009 by bed). Data probably on request (Rutgers HSRL).
- Regime: Delaware Bay salinity ~7-28 (upper seed beds <15; Shell Rock/Bennies ~15-22; HM/lower ~20-28); summer T 25-28 C => DD>20 likely >400. NOT cool+salty.
- NOAA Mussel Watch histopath also has DE Bay oysters (see above).

## CANADA - P. marinus NO LONGER ABSENT (major update)
- CFIA confirmations: New Brunswick, Spence Cove (with MSX), notified 2024-11-08, confirmed 2024-11-19 (statement 2024-11-21) = first in Canada; Nova Scotia, Merigomish (statement 2024-11-27); Newfoundland & Labrador, Notre Dame Bay (2025-06-05); PEI, Egmont Bay (2025-07-04); Quebec, Chaleur Bay (2025-07-04, with MSX).
- 2025-09-02: all waters of NB, NS, PEI, NL and eastern Quebec declared areas where MSX and Dermo are present or highly likely. Domestic Movement Control Program. Map: inspection.canada.ca MSX-Dermo-DIA.jpg (image only).
- Statements contain NO prevalence, intensity, sample counts, sampling dates, or mortality data. No public dataset. Source URLs: https://inspection.canada.ca/node/665 and canada.ca/en/food-inspection-agency/news/2024/11/..., 2025/06/..., 2025/07/...
- Relevance: these are cool waters (Gulf of St Lawrence / Northumberland Strait shallow bays summer 20-25 C, salinity ~25-31; Notre Dame Bay colder) => presence-only evidence that P. marinus establishes in cool+salty water; intensity unknown. DFO prior page said absent in Canada.
- CFIA "Locations infected" table (https://inspection.canada.ca/en/animal-health/aquatic-animals/diseases/reportable-diseases/perkinsus-marinus/locations-infected, current to 2026-09-30): 19 confirmed detections, province-level only, all C. virginica: 2024 NB x4 (Nov 19, Nov 25, Dec 4, Dec 13), NS Nov 25; 2025 NB Apr 25, NS May 27, NL May 29, QC Jun 27, PEI Jul 2, NB Dec 1, Dec 2, Dec 10 x2; 2026 NS Feb 16, NL May 14, QC Jun 22, NB Jun 26, NS Jul 23. No waterbody, no prevalence/intensity.
- DFO species page (dfo-mpo.gc.ca/science/aah-saa/diseases-maladies/pmdoy-eng.html) updated to note detections.

## CONNECTICUT / Long Island Sound (CT DoAg Bureau of Aquaculture, DABA) [DOWNLOADED dl/ct/]
- Histology surveillance 1997-2016 by state pathologist Inke Sunila: 8-37 sample sets/yr x ~30 oysters, Dermo intensity 0.5-5 (Mackin-like: 0.5 very light ... 5 very heavy); "weighed intensity" = sum scores / n. 2017-18 no data. 2019-2021 triplex PCR by R. Smolowitz (Roger Williams Univ ADL), USDA-APHIS funded.
- 2021 Statewide Shellfish Disease Update PDF (dl/ct/ct2021.pdf): coastline ANNUAL AVERAGE Dermo prevalence & WI 1997-2021 in figures only (prev ~ high since 1997; WI mostly <2; PCR WI higher). No station table. Dermo "takes up to 3 yrs in CT" to reach lethal intensity; mortality only in slow-growth/restoration areas.
- 2024 Fall Survey Report (html, portal.ct.gov/doag/knowledge-base/articles/aquaculture/2024-oyster-disease-report): 14 samples, 418 oysters, 14 sites listed by name/town (Stony Creek Branford, L598 Westport, Mid Quinnipiac, Fence Creek, Hammonasset, Housatonic Track 6, Mayhew Greenwich, L49 Norwalk, Ledyard Ground Thames, Sasco Beach, L-420 Stratford, East River Guilford, Clam Point Stonington, Milford Harbor). Statewide Dermo prev 87.35%, WPrev 1.646 (MSX 12.89%). Pooled only; no per-site values, no T/S.
- Station-level 1997-2016 histology data exist at DABA (not public) -> request (lydia.bienlien@ct.gov per page).
- LIS regime: S ~24-30 (eastern higher), summer T 22-25 => DD20 ~ 200-400 at eastern/central LIS. Borderline-relevant.
- NOAA Milford (NEFSC) monthly triplex qPCR Dermo/MSX/SSO at 4 LIS sites + water quality (project page fisheries.noaa.gov "Understanding the Health of Long Island Sound's Oyster Beds") - data not located.

## MASSACHUSETTS
### Green-Beach 2011 MS thesis Rutgers "Dermo disease on Martha's Vineyard" [DOWNLOADED dl/mv/gb2011.pdf] doi:10.7282/T3251HJ6
- 4 coastal salt ponds (barrier-beach, periodically breached): Edgartown Great Pond (EGP; first Dermo 1994, 85% prev 1995 + mass mortality), Tisbury Great Pond (TGP; first 1999 8%, ~75% mortality 2000), Oyster Pond (OYP; 2004 100%/WP 3.1), Squibnocket Pond (SQP; uninfected; S 10.5-13).
- 2008-2009 summers, ~30 oysters/sample, RFTM, prevalence & WP (0-5); Hobo T loggers 3-h 2008; salinity at each collection. Values in figures (TGP WP peak 1.55 Aug-2008, 2.2 Sep-2009; prev 100% Sep). EGP S>15 2008, <15 2009; OYP S<10 Jul-Oct 2008 yet WP to 1.63.
- Ponds are BRACKISH (mostly 10-25 ppt) -> not salty; data in figures (digitization needed).
### Ford & Smolowitz 2007 Mar Biol 151:119-133 "Infection dynamics of an oyster parasite in its newly expanded range" (paywalled)
- 1996-98 sampling at multiple sites across new range (Delaware Bay to Cape Cod), embayments & estuaries, wild+cultured; high prevalence persists; dynamics resemble southern range; persistence attributed to warmer winters + shallow sites with warm summers. Quantitative tables not retrieved.
### Cape Cod: 2006 press (Cape Cod Coop Extension/SEMAC): Wellfleet 2005 cohort 100% prevalence, intensity 3-4 (0-5). MA DMF requires disease screening for transfers by approved labs (results not public).

## VIRGINIA (VIMS Molluscan Health / Shellfish Pathology)
- Program: Dermo & MSX monitoring since 1960; fall surveys of major VA populations; monthly (later quarterly) on 4 upper James River reefs since 1987. RFTM, Mackin 0-5, weighted prevalence (0-1 minimal, >~1.5 dangerous).
- Public annual reports found: 2002 (monitoring_rept02.pdf), 2003, 2004, 2005, 2009 (covers 2006-2008) at https://www.vims.edu/research/units/labgroups/molluscan_health/_documents/monitoring_rpt_YYYY.pdf ; publications index .../molluscan_health/publications/index.php. Downloaded 2009 -> dl/vims/ (tables 1A-C by site, text PDF, transcribable). Also W&M ScholarWorks reports/602 (2004).
- T from NOAA CO-OPS Sewells Point; salinity 10-25. DD20 high (~600-900). NOT cool. No machine-readable dataset found.

### Pecher 2007 PhD diss. (UMD; = Pecher et al. 2008 J Parasitol 94:410-422) [DOWNLOADED dl/pecher/diss.pdf]
- PCR (genus + species-specific), pooled gut/gill/mantle, monthly Jun-Sep 2002, 18-60 oysters/site/month.
- Table 4 (C. virginica, N / P. marinus n (%)): Walpole ME (Damariscotta R.) 116 / 1 (0.9%) [genus Perkinsus 12 (10.3%); monthly genus prev 19.4, 10.3, 0, 10.7% Jun-Sep]; Martha's Vineyard MA 180 / 173 (96.1%); Narragansett Bay RI 112 / 64 (57.1%); Oyster Bay NY 99 / 18 (18.2%); Delaware Bay NJ 111 / 84 (75.7%); Sandy Point MD 57 / 54 (94.7%).
- => Damariscotta PCR time line: 2002 0.9% (Pecher) -> 2012 Jones Cove 65.2% (Marquis et al 2015) -> 2016-17 site means 34-39% (Marquis et al 2020). All PCR presence, no intensity/histology.

## Great Bay WP (computed from counts) by year, Adams Pt / Nannie Is / Oyster R / Woodman Pt:
1997-2003 ~0-0.3; 2004 1.08/0.68/0.20/-; 2006 1.80/0.90/1.52/1.42; 2009 2.17/0.92/2.42/1.70; 2012 1.30/2.00/2.05/3.10; 2014 2.40/3.33/2.60/2.15. Prevalence 2009-2014 mostly 0.75-1.0.

## REVIEWS (quantitative content verified where noted)
- Ford 1996 J Shellfish Res 15:45-56 (range extension NE US; response to climate change?): 1990-91 outbreaks DE Bay -> Cape Cod (~500 km); repeated introductions + warm winters; EGP Martha's Vineyard 85% prev 1995.
- Cook et al. 1998 ECSS 46:587-597: SST 5 stations; strong winter-SST rise 1986-1991 (smaller rise post-1960 at 4 stations), no summer trend; DE Bay prevalence more strongly associated with winter than summer T. (from Rutgers abstract summary)
- Ford & Smolowitz 2007 Mar Biol 151:119-133: 1996-98 multi-site new range; high prevalence persists; dynamics like south; persistence due to warmer winters & shallow warm summer sites. (abstract summary; tables not retrieved)
- Carnegie et al. 2021 Sci Rep 11:12872 (rapid phenotype change / hypervirulence explains 1980s emergence better than climate alone; range expansion north ~1986+, DE Bay 1990, LIS & S New England 1991-94).
- Fernandez Robledo et al 2018 (Penobscot Bay <1 d/yr above 20C; +2C ~6 d; +4C ~20 d).
- Burge et al. 2014 Annu Rev Mar Sci 6:249-277 - not retrieved.
