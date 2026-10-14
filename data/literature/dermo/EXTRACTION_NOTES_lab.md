# Perkinsus marinus lab/controlled-experiment rates vs T and S -- research notes
Started 2026-10-14. Downloads in dl/ (untrusted). Scripts in scripts/. CSVs in this dir.

## Chu & La Peyre 1993, DAO 16:223-234 (open, int-res d016p223.pdf)  DOI 10.3354/dao016223
- In vivo, 10/15/20/25 C, 17.5 ppt, 1e6 trophozoites shell cavity, 46 d. Prevalence 23/46/91/100 % (n=44). Intensity fig not digitizable. -> chu_lapeyre1993_temperature.csv
- Lysozyme highest at 10 C; CI decreased with T. No mortality reported.

## Chu & Volety 1997, DAO 28:61-68 (open, d028p061.pdf) DOI 10.3354/dao028061
- USES DAMARISCOTTA RIVER (ME) OYSTERS (naive). Dose-response: min infective dose 1e2 meronts/oyster (shell cavity); prev 0/13/50/71% (92 d, 24.5C/13.6ppt) and 0/13/67/80% (60 d, 22C/20.6ppt) for 10,1e2,1e4,1e5 meronts; prezoosporangia less infective.
- T x S x dose: 10,15,25 C x 3,10,20 ppt x 0/2.5e3/2.5e4 meronts, 60 d. Fig 4 digitized by eye (3-D chart, rough). Infections DO establish at 10 C/20 ppt (~27% at 2.5e4) but are light (WP ~0.25 on 0-1-3-5 scale). Temperature most important, then dose, then salinity; T x S and T x dose interactions significant for intensity.
- No mortality during 60-d challenge; acclimation mortality Table 1 (3 ppt) is host stress.
- CSVs: chu_volety1997_dose_response.csv, chu_volety1997_TxS.csv, chu_volety1997_acclimation_mortality.csv

## Ford, Chintala & Bushek 2002 (I), DAO 51:187-201, DOI 10.3354/dao051187 (open)
- Mook Sea Farm (Damariscotta) hatchery oysters, held 10C/25ppt, challenged 26-28C/25ppt with 1e6 cells/g wwt.
- Table 2: log10 burdens of dead oysters (late dead 5.8-6.7), survivors 4.1-5.4. Discussion p199: oysters died with 1e6-1e7 cells/g wwt in most trials, ~1e4 in Expt II. -> ford2002_lethal_burden.csv

## Bushek, Ford & Chintala 2002 (III), DAO 51:217-225, DOI 10.3354/dao051217 (open)
- Fecal shedding log-log regression on body burden; days to death vs shedding. -> bushek2002_fecal_shedding.csv

## La Peyre, Casas & La Peyre 2006, DAO 71:59-74 (open) DOI 10.3354/dao071059
- In vitro, 28 C, 7/11/15/25/35 ppt. P. marinus: no net growth at 7 ppt (viability reduced even after acclimation); slow at 11; fast and similar at 15, 25, 35 ppt (r ~0.8/d in exponential phase after ~6-d lag). -> lapeyre2006_invitro_salinity.csv
- So 35 ppt is NOT inhibitory in vitro at 28 C.

## Volety 1995 PhD dissertation (VIMS, W&M ScholarWorks etd/1539616895; open PDF)
- Ch.2 = Chu & Volety 1997 data (Fig 4 same 3-D chart; Fig 6a/6b marginal means). Fig 6a marginal means by salinity look inconsistent with paper Fig 4 (20 ppt row ~0.8/1.6/2.1 at 10/15/25 C) -> NOT used.
- Ch.4 in vitro: meronts 4e6 cells/ml, 48 h at 4, 12, 20, 28 C (N=6). Cell numbers only visible as clusters in Fig 5 scatter (x 1e6): ~0.3-1 (4C?), ~1.3-1.6 (12C?), ~3.3-4.8 (20C?), ~11-15.6 (28C?). Cluster-to-temperature assignment inferred (monotone), units/volume ambiguous. Text: 'Temperatures above 28C are inhibitory to the growth (Volety & Chu unpubl.)'; 'growth rate of P. marinus meronts in culture is 0.7-1 doubling per day' (p~?? line 5017).
- Osmolality 400/570/840 mOsm (~14/20/28 ppt) at 28C, 48 h: growth highest at 570 (20 ppt); 28 ppt and 14 ppt lower (Fig 6 cluster ~3-4.5 vs ~7-8.6 x1e6).

## Queiroga et al. 2016 Parasitology 143:475-487 (open, PMC4800715) DOI 10.1017/S0031182015001602
- Brazilian (tropical) isolate. 15/25/35 C at 20 psu; 5/20/35 psu at 25 C; 15 d. No net proliferation at 15 C (viability fell to 58-68%) but cells recovered at 25 C. -> queiroga2016_invitro_TS.csv. Caveat: tropical isolate; temperature optimum likely shifted relative to Atlantic isolates.

## Ragone & Burreson 1993 JSR 12:1-7 (open via archive.org JSR v12 OCR)
- Naturally infected oysters, 23 C, 6/9/12/20 ppt, 56 d. Infections persist at 6 ppt (no expulsion) but do not progress at 6-9 ppt; critical range 9-12 ppt. Mortality 9.1/11.6/21.1/27.8%. -> ragone_burreson1993_salinity.csv

## Chu 1996 review JSR 15:57-66 (archive.org JSR v15, pp 57-66, PDF pages 63-72)
- Clearer replots: Fig 4 (Chu & La Peyre 1993 temperature WP), Fig 5 (Chu, La Peyre & Burreson 1993 JIP salinity exp: 3/10/20 ppt, 35 d, 1e6 meronts), Fig 6 (Chu & Volety TxS; N=7-15/group, 3-4 at 25C/3ppt), Fig 7 + Table 1 (Chu, Volety & Constantin 1996: C. virginica vs C. gigas at 10/15/25 C and 3/10/20 ppt; C. virginica deaths after challenge 1/3/10 of 80 at 10/15/25 C).
- Text: meront proliferation in culture at 5, 12, 20, 28 C positively correlated with T (Volety 1995).
- Chu 1999 (Bull EAFP 19:265, open on VIMS ScholarWorks vimsarticles/1724): prezoosporangia survive up to 4 d at 9 C (Chu & Greene 1989); winters <=5 C for 6-8 wk did not lower next summer's infection (Burreson & Ragone Calvo 1996).
- CSVs: chu_lapeyre1993_temperature.csv (now with WP), chu1993_salinity.csv, chu_volety_constantin1996_Cvirginica.csv

## Saunders, Powell & Lewis 1993 JSR 12:229-240 (archive.org v12)
- In vivo (field pond, 28-32 C, 17 ppt) DNA-synthesis doubling times: 1-10 h at low density; density-dependent; eq.19 printed: log10(Td[h]) = 0.014954 exp(0.47364 log10(cells/g dry wt)), 1e4-1e9 cells/g dw, R=0.69. CAUTION: as printed the equation gives Td ~ 1-11 h over the whole range, contradicting the text (1e4 -> ~1 d, 1e7 -> ~1 yr) and Fig 10; digitize Fig 10 (PDF page 250 of archive.org v12) instead. Already used in Hofmann et al. 1995 (in repo).

## La Peyre 1996 review JSR 15:89-101 (archive.org v15, Table 3 p.93) -> lapeyre1996_invitro_optima_summary.csv
- KEY for low T: Dungan & Hamilton 1995 (MTS assay) found proliferation at 10 C, increasing with T 10-35 C; Gauthier & Vasta 1995 (OD600, 10 d) found NO significant growth 4-20 C; La Peyre & Faisal: slow growth at 17 C; optimum 28 C everywhere. Salinity optimum ~24-36 ppt (G&V: 36 not lower than 24-30); D&H near-optimal 17-34 ppt, max 28 ppt.
- Primary D&H 1995 (J Eukaryot Microbiol 42:379, DOI 10.1111/j.1550-7408.1995.tb01598.x) and G&V 1995 (JIP 66:156, DOI 10.1006/jipa.1995.1079) are paywalled; actual relative-proliferation values at 4-20 C not obtained.

## O'Farrell et al. 2000 JSR 19:139-145 (archive.org v19) -> ofarrell2000_invitro_osmolality.csv
- 28 C: doubling 35 h at 6.5-9.7 ppt, 23-26 h at 12.7-27 ppt (acclimated cells).

## Ragone Calvo & Burreson 1994 JSR 13:123-130 (archive.org v13) -- overwintering (field + lab warm-up)
- York River tray (18.8-23 ppt) Nov 1991-May 1992: FTM tissue prevalence 100% (Nov) -> 32% (May); WP 2.0 -> 0.56. Holding samples 30 d at 25 C ALWAYS raised prevalence (e.g. Apr 36->52% hemolymph / 71% tissue; May 36->76% / 100%). Immunoassay: meronts persist at very low intensity in most oysters. Decline earlier/larger at James R. (4-12 ppt). => winter 'regression' is largely to sub-detection levels, not clearance, at high salinity.
## Brousseau & Baglivo 2000 JSR 19:133-138: Long Island Sound field model, WP vs lagged cumulative temperature (53-d lag). Field, not lab.

## Ragone Calvo et al. 2003 DAO 56:75-86 (open) -> ragone_calvo2003_sentinel_transmission.csv
- Naive Damariscotta sentinels at York River, biweekly; acquisition peaks Aug-Sep (max 71% per 15 d in 1994, 29% in 1995) at 24-28 C; sporadic low-burden acquisition at all seasons incl. 16% in Feb 1995 (~4 C) with 2-28 cells/oyster (after 30-d warm hold). Temperatures read from Fig 3A by eye.
- Spearman: T vs water-column cells 0.78; T vs host mortality 0.78; T vs incidence 0.30 (ns).

## Cold survival (in vitro) -> cold_survival_invitro.csv
- La Peyre 2008 Parasitology (abstract via PubMed 18248687): 4 C 30 d: 49% viable, no proliferation; progressive cooling to 2-6 C: 77% (30 d), 54% (60 d). Paywalled.
- La Peyre 2010 JIP (PubMed 20600094): synergy of low T and low S; 0% viability only at 3.5 ppt (2 C by d30; 10 C by d60); <1% at 7 ppt by d90; 15 and 25 ppt not reported in abstract (much higher). Paywalled.
- Ford & Chintala 2006 JEMBE 339:226: no cold adaptation in northern isolates; similar net proliferation 5-20 C. Paywalled; abstract not machine-retrievable (ScienceDirect 403).
- Gignoux-Wolfsohn et al. (bioRxiv 10.1101/2020.11.03.366781; Parasitology 2021): release from live oysters higher at 30-32 vs 23-25 C (76% of oysters), diurnal; salinity (5-25 ppt) no effect on release but higher intensity at 20-25 ppt after 4 months.

## Hewatt & Andrews 1956 Proc NSA 46:129-133 (archive.org) -> hewatt_andrews1956_temperature_mortality.csv
- Experimentally infected oysters: 15 C for 6 wk -> 10% mortality, infections persist (light) but arrested; 28 C closed vat -> 99% dead within 4 wk; pier 26-30 C -> 53%.

## PAYWALLED / NOT OBTAINED (user could supply PDFs)
- Dungan & Hamilton 1995 J Eukaryot Microbiol 42:379-388, DOI 10.1111/j.1550-7408.1995.tb01598.x (relative proliferation 4-40 C, 320-1920 mOsm)
- Gauthier & Vasta 1995 J Invertebr Pathol 66:156-168, DOI 10.1006/jipa.1995.1079 (4-38 C; 12-42 ppt)
- La Peyre, Casas, Villalba & La Peyre 2008 Parasitology 135:505-519, DOI 10.1017/S0031182008004150 (15 and 4 C time courses)
- La Peyre, Casas, Gayle & La Peyre 2010 J Invertebr Pathol 105:176-181, DOI 10.1016/j.jip.2010.06.010 (2/10 C x 7/15/25 ppt, 90 d)
- Ford & Chintala 2006 J Exp Mar Biol Ecol 339:226-235, DOI 10.1016/j.jembe.2006.08.004 (isolate proliferation 5-35 C)
- Chu & Greene 1989 J Invertebr Pathol 53:260-268, DOI 10.1016/0022-2011(89)90016-5 (prezoosporangia development vs T, S)
- Chu, La Peyre & Burreson 1993 J Invertebr Pathol 62:226-232, DOI 10.1006/jipa.1993.1104 (salinity exp; temperature not in open sources)
- Ford & Smolowitz 2007 Mar Biol 151:119-133, DOI 10.1007/s00227-006-0454-6 (field dynamics MA/ME)
- Ray 1954 Rice Inst. Pamphlet special issue (Rice repository hdl 1911/18411; JS-only site, not fetched)
- Chu, Volety & Constantin 1996 JSR 15:375 primary PDF at VIMS ScholarWorks is image-only; numbers taken from Chu 1996 review Fig 7/Table 1.


## User-supplied PDFs (user_pdfs/, 2026-10-14)
Text via pypdfium2 (txt/); scanned PDFs (chu1993, gauthier1995) read from rendered PNGs (img/u/). Scripts: scripts/pdftext.py, rend.py, paths.py, digitize_fs07.py, fs07_tempstats.py, write_user_pdf_csvs.py, make_constraints_v2.py.

### Ford & Smolowitz 2007 Mar Biol 151:119 (ford2006.pdf)
- NO MAINE data: sites Wellfleet, Cotuit, Martha's Vineyard (Butler's Neck, Mashacket Cove) MA; Westport, Milford CT; Oyster Bay NY; Great Sound, Delaware Bay NJ. Jun 1996-Jun 1998, n=25, RFTM 0-5.
- Figs 2-3 are vector graphics: prevalence, WP, cumulative mortality and temperature traces recovered EXACTLY from PDF path coordinates (scripts/digitize_fs07.py). Caveats found: Oyster Bay panels span 24 months (start assumed 1 Jun 1996, +/-1 mo); Fig 2 Milford/Westport panels use an unlabeled 36-month axis (Dec 1995-Dec 1998) -> re-mapped to Fig 3 sampling dates (resid <3 d).
- Derived (ford2007_smolowitz_site_temperature_metrics.csv): 1997 DD>20 C: Wellfleet 227, Cotuit 352, MV 314, Oyster Bay 153, LIS shore station ~115-120 (max 22.9 C); winter 1996-97 min 1.7-3.6 C, 51-106 d <5 C.
- Westport CT reached prev ~100%, WP 2.8-3.1 with only ~115-130 DD>20 (shore-station T); Milford CT with same T record stayed <=51% prev, WP <=0.9 -> thermal index alone does not determine outcome.
- Salinity Table 3: Wellfleet 30-32, Great Sound 28-31, Cotuit 25-27, Oyster Bay 21-27, Delaware Bay 14-26, Edgartown 6-24.
- WP lags monthly T by 4-5 mo in MA (3 mo south). Proliferation window >18-20 C for 10-14 wk suffices.
- CSVs: ford2007_smolowitz_fig2_fig3_digitized.csv, ford2007_smolowitz_fig2_temperature_trace.csv, ford2007_smolowitz_site_temperature_metrics.csv, ford2007_smolowitz_sites_tables.csv

### Ford & Chintala 2006 JEMBE 339:226 ('ford2006 (1).pdf')
- 9 isolates MA/NJ/SC, tetrazolium 72 h at 5-35 C. Ratio to T0 <1 at 5-10 C (decline 25-35%), net growth from 15 C, 15-20 C equal among regions (no cold adaptation); SC higher at 25-35 C. Fig 4 digitized by eye (+/-0.03).
- Table 5: field onset of prevalence rise 22-23 C at MA sites (13-16 C LIS, 15 C Great Sound). Mentions Maine first report 1995 (Kleinschuster & Parent).
- CSVs: ford2006_chintala_invitro_isolates.csv, ford2006_chintala_invivo_onset_temps.csv

### Dungan & Hamilton 1995 J Euk Microbiol 42:379 (dungan1995.pdf)
- MTS 72 h, 4-40 C: no net growth at 4 C; ratio 1.30 at 10 C; 2.98 at 15 C; max 35 C. Osmolality 340-1930 mOsm at 28 C: max 790 (~25 ppt); 30 ppt 85%, 37 ppt 70% of max. Doubling 13 h optimized. Signal = biovolume. -> dungan1995_invitro_TS.csv

### Gauthier & Vasta 1995 JIP 66:156 (gauthier1995.pdf, scanned)
- Temperatures tested only 4, 20, 28, 32, 38 C (OD600, 10 d): 4 = 20 C (no detectable growth), optimum 28-32. Salinity 12-42 ppt at 28 C: 24-36 equal, 42 lower. Maine (Mook Sea Farm) oysters at 26 C/30 ppt: 4/4 died within 5 wk of cultured-cell challenge. -> gauthier1995_invitro_TS.csv

### Chu, La Peyre & Burreson 1993 JIP 62:226 (chu1993.pdf, scanned)
- TEMPERATURE RESOLVED: 21.0+/-1.0 C (methods) / 22.0+/-1.0 C (abstract, Table 1). 35 d after 1e6 cells: prevalence 50/79/82% (3/10/20 ppt), controls 4/19/35% (text; Chu 1996 replot values 7/12 for controls were wrong). WP 0.62/1.76/1.87. -> chu1993_salinity_primary.csv (supersedes chu1993_salinity.csv)

### Chu & Greene 1989 JIP 53:260 (chu1989.pdf)
- Prezoosporangia: no development at 4 C; survive 1-3 d at 4 C (2/3 fail after 4 d); 18-20 C sporulate in 2-8 d vs 2-3 d at 28 C; no sporulation <=4 ppt; 32 ppt as good as 11-23. Zoospores die within 1 d at 4 C. Abstract '64 days' at 4 C not supported by body. -> chu1989_prezoosporangia_TS.csv

### La Peyre et al. 2008 Parasitology 135:505 (peyre2008.pdf)
- 4 C direct: viability 49% at d30, no division. 15 C: lag ~16-20 d, r_net 0.02-0.09/d vs ~0.48/d at 28 C. Gradual cooling to 2-6 C: 91% end of cooling, 78/54/25% after ~30/60/90 d. Direct 28->2 C: 0% by d90. Fig 7 digitized by eye (+/-3%). -> lapeyre2008_invitro_lowT.csv

### La Peyre et al. 2010 JIP 105:176 (lapeyre2010.pdf)
- After 0.5 C/d acclimation: 25 ppt 2 C 95/88/45/21.9% (d1/30/60/90), 10 C 99/89/61/43.2%; 15 ppt similar; 7 ppt and 3.5 ppt lethal at low T. -> lapeyre2010_invitro_TxS.csv

### Constraints
- key_constraints_for_ude_v2.csv written (v1 kept; 4 abstract-level rows replaced; 26 new rows, maine_priority column).
