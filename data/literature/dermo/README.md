# Dermo literature data (extracted, not fetched)

Numbers transcribed or digitized from published papers and reports for the Dermo modelling work
(PLAN.md P4.2c-P4.2d). They were extracted by hand and by figure digitization on 2026-10-14; each CSV
header and the two EXTRACTION_NOTES files give the paper, page, table or figure, and the reading
uncertainty for digitized values. Papers supplied by the user as PDFs are not stored here.

- `EXTRACTION_NOTES_lab.md`: laboratory and controlled-experiment rates vs temperature and salinity
  (Chu & La Peyre 1993; Chu & Volety 1997; Chu, La Peyre & Burreson 1993; Chu & Greene 1989; Dungan &
  Hamilton 1995; Gauthier & Vasta 1995; La Peyre et al. 2006, 2008, 2010; O'Farrell et al. 2000; Ford &
  Chintala 2006; Ford, Chintala & Bushek 2002; Bushek et al. 2002; Ragone & Burreson 1993; Ragone Calvo
  et al. 2003; Hewatt & Andrews 1956; Queiroga et al. 2016) and field series from Ford & Smolowitz 2007.
- `key_constraints_for_ude_v2.csv`: condensed constraints, with a `maine_priority` column.
- `EXTRACTION_NOTES_field.md`: search for Northeast / Maine field data (Great Bay NH, Maine PCR studies,
  Delaware Bay, Connecticut, Canada); `mw_oyster_dermo_by_site_year.csv` summarizes NOAA Mussel Watch
  per-oyster Dermo scores (East Coast oyster sites, 1988-2010).

The Great Bay survey table used in the pipeline is in `data-pipeline/src/oyster_data/sources/nh_greatbay.py`.
