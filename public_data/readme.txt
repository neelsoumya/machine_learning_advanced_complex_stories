Data

Mortality data
	https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables

	https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregistrationsummarytables/2017

	https://www.gov.uk/government/statistics



https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths

Higham:

     https://www.scotpho.org.uk/media/1403/inphorm-6-final.pdf






Notes from Rudolf

https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathregistrationssummarytablesenglandandwalesreferencetables

... deathsummarytables2017final.xls

and possibly deathsyoacause2017...

... nope! perhaps not.


"full" method = deaths / expected deaths PER YEAR using year-specific, age-specific population tables

... or "quick" method = assume all years are like 2017 (for example)

(1) an SMR table with columns (from ONS data)

- year

- age_band

- age_standardized_mortality_rate_this_year_this_age_range


(2) table with columns:

- patient_id

- year

- proportion_of_year_patient_exposed [i.e. in CPFT, alive]

- did_patient_die_in_this_period [0, 1]

- age_band_this_period (e.g. 50-54)

- age_standardized_mortality_rate_this_year_this_age_range (EXPRESSED AS: deaths per person, not per 1000 or 100,000 people)


derived columns:

- expected_n_deaths := proportion_of_year_patient_exposed * age_standardized_mortality_rate_this_year_this_age_range

- actual_n_deaths := did_patient_die_in_this_period

across all patients:

SMR = sum(actual_n_deaths) / sum(expected_n_deaths)