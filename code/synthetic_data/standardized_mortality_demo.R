#!/usr/bin/env Rscript

# Principles of calculating standardized mortality rates.
# Rudolf Cardinal, 17 Sep 2019 + 3 Feb 2020.

library(data.table)
library(ggplot2)
library(lubridate)
library(tidyr)
set.seed(54321)  # arbitrary number; for consistency across runs


SEP1 <- "===============================================================================\n"
SEP2 <- "-------------------------------------------------------------------------------\n"


# =============================================================================
# Framework data: years -- SPECIMEN
# =============================================================================

make_years <- function(start_year, end_year)
{
    years <- data.table(year = start_year:end_year)
    years[, year_start_date := as.Date(paste0(year, "-01-01"))]
    years[, year_end_date := as.Date(paste0(year, "-12-31"))]
    return(years)
}


# Define years of interest
START_YEAR <- 2013
END_YEAR <- 2018
YEARS <- make_years(start_year = START_YEAR, end_year = END_YEAR)


# =============================================================================
# Framework data: age ranges -- ONS standard
# =============================================================================
# Define (inclusive) age ranges of interest, according to the UK Office of
# National Statistics (ONS) standard.
# The ONS uses 5-year bands, e.g. 10-14, 15-19, etc., except for initial
# categories "under 1" and "1-4". Its final category is "90 and over".

MALE <- "M"
FEMALE <- "F"
SEXES <- data.table(sex = c(MALE, FEMALE))

AGE_RANGES <- data.table(age_range_start = c(0, 1, seq(5, 90, 5)))
AGE_RANGES[, age_range_end := age_range_start + 4]
AGE_RANGES[age_range_start == 0, age_range_end := 0]
AGE_RANGES[age_range_start == 1, age_range_end := 4]
AGE_RANGES[age_range_start == 90, age_range_end := Inf]
AGE_RANGES[, age_range := paste0(age_range_start, "-", age_range_end)]
AGE_RANGES[age_range_start == 90, age_range := "90+"]

get_age_range <- function(age, age_ranges = AGE_RANGES)
{
    # Returns a descriptive age range from an actual age.
    #
    # Args:
    #   age:
    #       Age of a person, or vector of ages
    #   age_ranges:
    #       Data table with columns age_range_start (numeric), age_range_end
    #       (numeric), age_range (character).
    #
    # Returns:
    #   For each age element: the descriptive age range, e.g. "40-45",
    #   or NA_character_.
    #
    # Test code:
    #
    #   get_age_range(70)  # "70-74"
    #   get_age_range(c(1,2,50))  # c("1-4", "1-4", "50-54")
    return(
        sapply(
            as.integer(age),  # vector over which to apply
            function(x) {
                ifelse(x < 0, NA_character_,
                       age_ranges[age_range_start <= x & x <= age_range_end]$age_range)
            }
        )
    )
}


# =============================================================================
# Population data -- FICTIONAL; would plug in ONS mortality data here
# =============================================================================
# e.g.
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathregistrationssummarytablesenglandandwalesreferencetables
# ... deathsummarytables2017final.xls, Table 1
# ... ?use East of England data, or UK data if regional data not available

DAYS_PER_YEAR <- 365  # for simplicity

make_fictional_population_mortality <- function(
    use_fixed_mortality_rate = TRUE,
    fixed_deaths_per_thousand_per_year = 100,
    years = YEARS,
    age_ranges = AGE_RANGES,
    sexes = SEXES
)
{
    # Create fictional population mortality data.
    # Args:
    #   use_fixed_mortality_rate:
    #       Boolean. Use a fixed mortality rate? If not, made-up values
    #       based on age and sex are used.
    #   fixed_deaths_per_thousand_per_year:
    #       Numeric. If use_fixed_mortality_rate is TRUE, use this number of
    #       deaths per thousand people per year.
    #   years:
    #       Data table defining years of interest (see above).
    #
    # Returns:
    #   A data table with columns:
    #   - ... those from year
    #   - ... those from age_ranges
    #   - ... those from sexes
    #   - deaths_per_thousand_per_year

    # Years and age ranges for which we require mortality rates:
    population_mortality <- data.table(tidyr::crossing(years, age_ranges, sexes))

    if (use_fixed_mortality_rate) {
        population_mortality[,
            deaths_per_thousand_per_year := fixed_deaths_per_thousand_per_year
        ]
    } else {
        population_mortality[,
            deaths_per_thousand_per_year := age_range_start ^ 2 / 40 *
                                            ifelse(sex == FEMALE, 0.8, 1)
        ]
    }

    # Some recalculations:
    population_mortality[,
        deaths_per_person_per_year := deaths_per_thousand_per_year / 1000]
    population_mortality[,
        deaths_per_person_per_day := deaths_per_person_per_year / DAYS_PER_YEAR]

    validate_population_mortality_table(population_mortality, years, age_ranges)

    return(population_mortality)
}


validate_population_mortality_table <- function(population_mortality,
                                                years, age_ranges)
{
    # Check we have enough mortality data (irrelevant here, with synthetic data,
    # but for when we scrape ONS data):
    stopifnot(all(years$year %in% population_mortality$year))
    stopifnot(all(age_ranges$age_range %in% population_mortality$age_range))
    stopifnot(all(!is.na(population_mortality$deaths_per_person_per_year)))
}


# Now for an entirely fictional set of mortality rates, in these examples
# static over time:
DEMOFIXED_DEATHS_PER_THOUSAND_PER_YEAR <- 100
DEMOFIXED_DEATHS_PER_PERSON_PER_YEAR <- DEMOFIXED_DEATHS_PER_THOUSAND_PER_YEAR / 1000
DEMOFIXED_DEATHS_PER_PERSON_PER_DAY <- DEMOFIXED_DEATHS_PER_PERSON_PER_YEAR / DAYS_PER_YEAR

DEMO_POPULATION_MORTALITY <- make_fictional_population_mortality(
    fixed_deaths_per_thousand_per_year = DEMOFIXED_DEATHS_PER_THOUSAND_PER_YEAR
)


# =============================================================================
# Simulate lifespan
# =============================================================================
# See https://onlinelibrary.wiley.com/doi/pdf/10.1002/sim.2059

simulate_lifespan_days <- function(n, deaths_per_person_per_day)
{
    # Return simulated lifespans in days.
    #
    # Args:
    #   n:
    #       (integer) Number of lifespans to simulate.
    #   deaths_per_person_per_day:
    #       (numeric) Mean rate of death.
    #
    # Returns:
    #   a vector of length n containing lifespans in days
    #
    # Uses the exponential distribution; see ?rexp
    # The mean of the resulting values is 1/rate.
    return(rexp(n = n, rate = deaths_per_person_per_day))
}


# =============================================================================
# Sample data
# =============================================================================

make_sample_data <- function(n_subjects, p_death_per_day,
                             study_start_date, study_end_date,
                             sexes = SEXES)
{


    # Invent some subjects:
    sample_data <- data.table(subject_id = 1:n_subjects)
    # With different sexes:
    sample_data[, sex := sexes$sex[subject_id %% 2 + 1]]
    # With dates of birth:
    sample_data[, dob := as.Date("1950-01-01") - subject_id * 30]  # try modifying the base year
    # Who enter the study as different times:
    sample_data[, subject_start_date := study_start_date + 5 * subject_id]

    # A subject can't start before the study does:
    sample_data[, subject_start_date := pmax(subject_start_date, study_start_date)]

    # Make them die probabilistically:
    sample_data[, lifespan_days_from_subject_start_date := as.integer(
        simulate_lifespan_days(
            n = nrow(sample_data),
            deaths_per_person_per_day = p_death_per_day)
    )]
    # Now calculate death dates/flags:
    sample_data[, death_date := subject_start_date + lifespan_days_from_subject_start_date]
    sample_data[, subject_end_date := pmin(death_date, study_end_date)]
    sample_data[death_date > subject_end_date, death_date := NA]
    sample_data[, died_during_study := as.integer(!is.na(death_date))]
    sample_data[, death_year := lubridate::year(death_date)]

    # Validity checks:
    stopifnot(all(sample_data$subject_start_date >= study_start_date))
    stopifnot(all(sample_data$subject_subject_end_date <= study_end_date))

    return(sample_data)
}


# =============================================================================
# Calculations
# =============================================================================

calculate_smr <- function(sample_data,
                          population_mortality,
                          years = YEARS, ci_alpha = 0.05) {
    # Calculate standardized mortality rates.
    #
    # Args:
    #   sample_data:
    #       Data table with one row per subject, and columns:
    #       - subject_id
    #       - sex (M, F)
    #       - dob (date)
    #       - subject_start_date (date)
    #       - subject_end_date (date)
    #       - death_date (date or NA)
    #   population_mortality:
    #       Data table with population mortality information, as above.
    #   years:
    #       Data table with columns:
    #       - year (integer)
    #       - year_start_date (date; 1 January for that year)
    #       - year_end_date (date; 31 December for that year)
    #   ci_alpha:
    #       Alpha for confidence interval.
    METHOD <- "

1. For every subject/calendar year combination during which the subject was
   alive and under 'observation' (i.e. when their death would have been be
   recorded, had they died), the subject's age that year was calculated as of 1
   January that year. Their start date during that year was taken as the later
   of 1 January that year and the time the subject entered 'observation'. Their
   end date during that year was taken as the earlier of 31 December and the
   date they died (or otherwise left observation). The year was scored as one
   in which they did or did not die whilst under observation.

2. For every subject/year combination, the subject's 'exposure' was calculated
   (from start date to end date within that year), as the fraction of the whole
   year during which they were under 'observation'.

3. For every subject/year combination, the expected number of deaths for that
   person-year was calculated as the published population mortality figure
   (expressed as deaths per person per year) for that subject's sex and age
   that year, multiplied by the proportion of the year to which the subject was
   'exposed' (as above).

4. The age-standardized mortality ratio (SMR) was then derived as the sum of
   observed deaths divided by the sum of expected deaths, across all
   subject/year combinations.

5. The standard error (SE) of the ASMR was calculated as
        sqrt(n_observed_deaths) / n_expected_deaths
   and confidence intervals as
        SMR ± 1.96 * SE
   according to ERPHO (2005).

   -- Eastern Region Public Health Observatory (erpho). “Standardisation.”
   INphoRM, June 2005. https://www.scotpho.org.uk/media/1403/inphorm-6-final.pdf.

This method was checked for validity by simulating artificial subjects dying
according to a daily hazard rate; the method accurately recovered SMRs from
simulated SMRs.

    "

    # Count subjects
    n_subjects <- nrow(sample_data)

    # Make a table with a row for each subject/year combination:
    subject_mortality <- data.table(tidyr::crossing(sample_data, years))

    # Filter to relevant years only:
    subject_mortality <- subject_mortality[
        year >= lubridate::year(subject_start_date) &
        year <= lubridate::year(subject_end_date)
    ]

    # Calculate a per-subject start date for a given year
    subject_mortality[,
        subject_start_date_this_year := pmax(year_start_date, subject_start_date)]

    # Calculate whether a death occurred in each year
    subject_mortality[, died_this_year := as.integer(
        !is.na(death_date) &
        death_date >= subject_start_date_this_year &
        death_date <= year_end_date
    )]

    # Calculate an age for each year.
    # (Slightly non-obvious because birthday might be midway through the year.)
    subject_mortality[, age_this_year :=
        lubridate::time_length(
            difftime(
                year_start_date,  # later date: start of year for this patient
                # ... ARGUABLE; could use pmax(subject_start_date, year_start_date)
                #     but probably better to be completely consistent across all
                #     patients
                dob  # earlier date: date of birth
            ),
            unit = "years"
        )
    ]

    # Demarcate into the age ranges supplied by the ONS:
    subject_mortality[, age_range_this_year := get_age_range(age_this_year, AGE_RANGES)]

    # Calculate the proportion of the year for which the subject is at risk of
    # death. This should clearly not include time before the subject enters the
    # study. However, should it include time after death? I suspect it should.
    # If someone dies on 1 Jan 2000, that is 1 death in 2000, to be compared to
    # an expected number of deaths in 2000 -- not to a period of time determined
    # by the fact of their death.
    # -- No, I was wrong; that method leads to bias.
    subject_mortality[, exposure_proportion_of_year :=
        # Subject's exposure: as above:
        # ... end date ARGUABLE but the assumption is probably that exposure
        #     lasts beyond death to the end of that year, as per:
        #       pmin(year_end_date, STUDY_END_DATE)  # NO, THIS IS BIASED.
        #     ... NOPE; biased; underestimates high mortality rates (verified by
        #     sim).
        # ... or until death, as per:
        #       pmin(year_end_date, subject_end_date)  # YES.
        lubridate::time_length(
            pmin(year_end_date, subject_end_date) -  # end
            subject_start_date_this_year,  # start
            unit = "days"
        ) /
        # ... divided by total ONS mortality period for this year:
        lubridate::time_length(
            year_end_date -  # end
            year_start_date,  # start
            unit = "days"
        )
    ]

    # Validity checks:
    stopifnot(all(subject_mortality$age_this_year >= 0))
    stopifnot(all(!is.na(subject_mortality$age_range_this_year)))
    stopifnot(all(subject_mortality$exposure_proportion_of_year >= 0))
    stopifnot(all(subject_mortality$exposure_proportion_of_year <= 1))

    # Merge in standardized mortality (we'll create a new table, for clarity):
    subject_and_population_mortality <- merge(
        subject_mortality,
        population_mortality[, c("year", "sex", "age_range", "deaths_per_person_per_year")],
        by.x = c("year", "sex", "age_range_this_year"),
        by.y = c("year", "sex", "age_range")
    )

    # Calculated expected number of deaths, for each person-year contributing:
    subject_and_population_mortality[,
        expected_n_deaths := exposure_proportion_of_year * deaths_per_person_per_year]

    # Calculate overall summaries:
    observed_n_deaths <- sum(subject_and_population_mortality$died_this_year)
    expected_n_deaths <- sum(subject_and_population_mortality$expected_n_deaths)
    age_standardized_mortality_ratio <- observed_n_deaths / expected_n_deaths
    se_asmr <- sqrt(observed_n_deaths) / expected_n_deaths
    asmr_ci_lower <- age_standardized_mortality_ratio + qnorm(ci_alpha / 2) * se_asmr
    asmr_ci_lower <- max(0, asmr_ci_lower)  # can't be negative
    asmr_ci_upper <- age_standardized_mortality_ratio + qnorm(1 - (ci_alpha / 2)) * se_asmr
    #...  as per ERPHO (2005), "Standardisation", Information on Public Health
    #     Observatory recommended methods
    stopifnot(observed_n_deaths <= n_subjects)

    # Report:
    cat(paste0(
        SEP2,
        "Age-standardized mortality rate calculations\n",
        SEP2,
        "n_subjects: ", n_subjects, "\n",
        "observed_n_deaths: ", observed_n_deaths, "\n",
        "expected_n_deaths: ", expected_n_deaths, "\n",
        "age_standardized_mortality_ratio: ", age_standardized_mortality_ratio, "\n",
        "asmr_ci_lower: ", asmr_ci_lower, "\n",
        "asmr_ci_upper: ", asmr_ci_upper, "\n"
    ))

    return(list(
        # input
        sample_data = sample_data,
        population_mortality = population_mortality,
        years = years,
        ci_alpha = ci_alpha,
        # output
        subject_mortality = subject_mortality,
        subject_and_population_mortality = subject_and_population_mortality,
        observed_n_deaths = observed_n_deaths,
        expected_n_deaths = expected_n_deaths,
        age_standardized_mortality_ratio = age_standardized_mortality_ratio,
        se_asmr = se_asmr,
        asmr_ci_lower = asmr_ci_lower,
        asmr_ci_upper = asmr_ci_upper
    ))
}


# =============================================================================
# Validation
# =============================================================================

plot_smr_from_sim <- function(population_mortality,
                              n_subjects = 200,
                              smr_values = c(0.2, 0.5, 1, 2, 5, 10, 20))
{
    simdata <- data.table(expected_smr = smr_values)
    simdata[, observed_smr := NA_real_]
    simdata[, observed_smr_ci_lower := NA_real_]
    simdata[, observed_smr_ci_upper := NA_real_]
    for (rownum in 1:nrow(simdata)) {
        expected_smr <- simdata[rownum, expected_smr]
        sim_p_death_per_day <- DEMOFIXED_DEATHS_PER_PERSON_PER_DAY * expected_smr
        sample_data <- make_sample_data(
            n_subjects = n_subjects,
            p_death_per_day = sim_p_death_per_day,
            study_start_date = as.Date(paste0(START_YEAR, "-01-01")),
            study_end_date = as.Date(paste0(END_YEAR, "-12-31")),
            sexes = SEXES
        )
        result <- calculate_smr(sample_data, population_mortality)
        simdata[rownum, observed_smr := result$age_standardized_mortality_ratio]
        simdata[rownum, observed_smr_ci_lower := result$asmr_ci_lower]
        simdata[rownum, observed_smr_ci_upper := result$asmr_ci_upper]
    }
    p <- (
        ggplot(simdata, aes(x = expected_smr, y = observed_smr))
        + geom_abline(slope = 1, intercept = 0,  # y = x (perfect)
                      linetype = "dashed", colour = "grey")
        + geom_errorbar(aes(ymin = observed_smr_ci_lower,
                            ymax = observed_smr_ci_upper))
        + geom_line()
        + geom_point()
        + scale_x_log10()
        + scale_y_log10()
    )
    return(p)
}


validation_fig1 <- plot_smr_from_sim(DEMO_POPULATION_MORTALITY, n_subjects = 20)
validation_fig2 <- plot_smr_from_sim(DEMO_POPULATION_MORTALITY, n_subjects = 2000)
