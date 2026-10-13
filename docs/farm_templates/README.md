# Farm record templates

These six spreadsheets let your farm's records improve the growth and risk model for your own site.
Open each `.csv` file in Excel, Google Sheets or LibreOffice, keep the first row (the column names)
exactly as it is, add one row per record, and save as CSV. Leave a cell empty if you did not measure it.

Dates are written `YEAR-MONTH-DAY`, for example `2025-06-15`. A **cohort** is any name you use for one
batch of seed, for example `2025-spring-triploid`. Use the same cohort name in every file.

| File | One row per | What to fill in |
|---|---|---|
| `farm_info.csv` | lease site | Farm name, a short site label, latitude and longitude in decimal degrees (for example 43.95, -69.55), Maine DMR growing area if known, notes. Your farm name is replaced by a code before analysis. |
| `growth_checks.csv` | oyster measured | Date, cohort, oyster number within that check (1, 2, 3, ...), shell height in mm (longest dimension, hinge to bill). Dry tissue weight in g only if you dried meats. Measuring 30 or more oysters per cohort per check is ideal; any number helps. |
| `mortality_counts.csv` | count of a bag, cage or sample | Date, cohort, live and dead oysters counted, and days since that group was last counted (leave empty if unknown). Count gapers and empty paired shells as dead. |
| `seeding.csv` | seed batch put out | Date, cohort, number of seed, average seed shell height in mm, ploidy (`diploid` or `triploid`), gear (for example `floating bags`), seed source (hatchery name). |
| `harvest.csv` | harvest day for a cohort | Date, cohort, number harvested, average shell height in mm if measured, price received per oyster in US dollars if you are willing to share it. |
| `temperature_logger.csv` | logger reading | Date and local time (`2025-06-15 14:30`), water temperature in °C, logger depth in m. Export straight from a HOBO or similar logger and rename the columns. |

Send the files back together. They are converted with `oyster-data farm-ingest` (see `data-pipeline/`).
Your farm and site names never leave the conversion step; analyses use an anonymous site code.
