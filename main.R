library(data.table)
library(lubridate)
library(ggplot2)
library(table1)
library(car) # Durbin Watson test
library(lmtest) # Scaling std errors
library(sandwich) # Newey West
source("code/functions.R")

# config ------------------------------------------------------------------

relevant_atc <-  c("R03BA", "R03AK")
intervention_date <- as.POSIXct("2021-11-02")
study_start <- as.POSIXct("2020-07-01")
icpc_comorbs <- list(
    kidneydmg = c("U99.01"),
    heartfailure = c("K77", "K77.01", "K77.03", "K77.04"),
    cvd = c("K74", "K74.01", "K74.02", "K75", "K76", 
                  "K76.01", "K77", "K78", "K90", "K92"),
    diabetes = c("T90"),
    copd = c("R95"),
    liver =  c("D97"),
    asthma = c("R96")
)
scaling_coeff <- 50

df_pat <- readRDS(fpath("pat.rds"))
df_eps <- readRDS(fpath("eps_all.rds"))
df_med <- readRDS(fpath("med.rds"))

# determine which episodes had ICS prescription ---------------------------

df_eps[, is_covid := dICPC %in% c("R83", "R83.03") & dBegindatum >= study_start]
get_episodes_w_med <- function(df_eps, df_med, atc) {
    regex_relevant_atc <-
        paste0("^(", paste(atc, collapse = "|"), ")")
    df_med <- df_med[grepl(regex_relevant_atc, dATC),]
    df_med[, med_id := .I]
    df_eps_med <- merge(df_eps, df_med)
    df_eps_med[, `:=`(
        delta_beg_med = difftime(dVoorschrijfdatum, dBegindatum, units = "days"),
        delta_mut_med = difftime(dVoorschrijfdatum, dMutatiedatum, units = "days")
    )]
    df_eps_med <-
        df_eps_med[(delta_beg_med >= 0 & delta_beg_med <= 14) |
                       (delta_mut_med >= 0 &
                            delta_mut_med <= 14),]
    df_eps_med <- unique(df_eps_med, by = "med_id")
    
    return(df_eps_med)
}
df_ics_eps <- get_episodes_w_med(df_eps[is_covid == T, ], df_med, relevant_atc)

# studypop ----------------------------------------------------------------

get_comorb_table <- function(df_eps) {
    df_comorbs <- df_eps[, .(kidneydmg    = any(dICPC %in% icpc_comorbs$kidneydmg),
                             heartfailure = any(dICPC %in% icpc_comorbs$heartfailure),
                             cvd          = any(dICPC %in% icpc_comorbs$cvd),
                             diabetes     = any(dICPC %in% icpc_comorbs$diabetes),
                             copd         = any(dICPC %in% icpc_comorbs$copd),
                             liver        = any(dICPC %in% icpc_comorbs$liver),
                             asthma       = any(dICPC %in% icpc_comorbs$asthma)),
                         by = "PATNR"]
    df_comorbs[, any_comorb := kidneydmg | heartfailure | cvd | diabetes |
                   copd | liver | asthma]
    return(df_comorbs)    
}

get_n_covid_table <- function(df_eps) {
    return(df_eps[is_covid == T, .(n_covid_episodes = .N), by = "PATNR"])
}

df_pat <- merge(df_pat, get_comorb_table(df_eps))
df_pat <- merge(df_pat, get_n_covid_table(df_eps))

df_earliest_covid <- df_eps[is_covid == T, 
                            .(earliest_covid_consult = min(dBegindatum)), 
                            by = PATNR]
df_pat[, earliest_covid := df_earliest_covid[
    match(PATNR, df_earliest_covid$PATNR), "earliest_covid_consult"]]

df_pat[, dob := make_date(year = iGeboortejaar, month = 7, day = 2)]
df_pat[, age := floor(decimal_date(earliest_covid) - decimal_date(dob))]
df_pat[, nhggroup := (age >= 65 | (age >= 50 & any_comorb))]
df_pat[, was_prescribed_ics := PATNR %in% df_ics_eps$PATNR ]
df_pat[, prescribed_before := PATNR %in% df_ics_eps[dVoorschrijfdatum < intervention_date]$PATNR]
df_pat[, prescribed_after := PATNR %in% df_ics_eps[dVoorschrijfdatum >= intervention_date]$PATNR]
df_pat[, is_male := dGeslacht == "M"]

# table1 ------------------------------------------------------------------

strata <- list(before = df_pat[prescribed_before == T, ], 
               after = df_pat[prescribed_after == T, ], 
               overall = df_pat[was_prescribed_ics == T, ],
               overall = df_pat[was_prescribed_ics == F, ])

labels <- list(
    variables=list(
        is_male = "sex (male)", 
        age = "age",
        nhggroup = "nhg group", 
        n_covid_episodes = "n episodes", 
        any_comorb = "any of the following comorbidities",
        asthma = "Asthma", 
        copd = "COPD", 
        cvd = "Cardiovascular disease", 
        diabetes = "Diabetes", 
        kidneydmg = "Chronic kidney damage", 
        liver = "Liver cirrhosis"
    ),
    groups=list("Prescribed", "Not Prescribed"))

table1 <- table1(strata, labels = labels, groupspan = c(3,1))

write.csv2(
    x = as.data.frame(table1),
    file = "results/2024_08_12_table1.csv", 
    row.names = F
)


# add extra info to episode table -----------------------------------------

add_nhg_info <- function(df_eps, df_pat) {
    df_pat <- df_pat[, .(PATNR, dob, dGeslacht, any_comorb)]
    df_eps_all <- merge(df_eps, df_pat)
    df_eps_all[, age := floor(decimal_date(dBegindatum) - decimal_date(dob))]
    df_eps_all[, nhggroup := age >= 65 | (age >= 50 & any_comorb)]
    return(df_eps_all)
}

df_eps <- add_nhg_info(df_eps, df_pat)
df_ics_eps <- add_nhg_info(df_ics_eps, df_pat)

# get weekly counts -------------------------------------------------------

get_weekly_counts <- function(df_eps) {
    df_studyperiod <- data.table(dBegindatum = seq(
        as.Date(study_start), as.Date("2022-08-01"), by = "days"))
    df_daily_eps <- df_eps[, .N, by = dBegindatum]
    df_daily_eps[, dBegindatum := as.Date(dBegindatum)]
    
    df_daily_eps <- merge(df_studyperiod, df_daily_eps, all.x = T)
    df_daily_eps[is.na(N), N := 0]
    df_daily_eps[, week_number := (.I - 1) %/% 7]
    df_weekly_eps <- df_daily_eps[, .(
        n_week = sum(N),
        week_datum = min(dBegindatum)
    ), by = week_number]
    return(df_weekly_eps)
}

make_count_table <- function(df_eps, df_ics_eps) {
    df_eps_counts <- get_weekly_counts(df_eps)
    df_ics_counts <- get_weekly_counts(df_ics_eps)
    setnames(df_ics_counts, "n_week", "n_ics_week")
    df_counts <- merge(df_eps_counts, df_ics_counts)
    df_counts[, `:=`(
        after_intervention = week_datum >= intervention_date,
        Q1 = month(week_datum) <= 3,
        Q2 = month(week_datum) >= 4 & month(week_datum) <= 6,
        Q3 = month(week_datum) >= 7 & month(week_datum) <= 9
    )][, after_intervention_t := pmax(0, cumsum(after_intervention) - 1)]
    return(df_counts)
}

df_covid_eps <- df_eps[is_covid == T]

list_counts <- list(
    df_counts_overall = make_count_table(df_covid_eps, 
                                          df_ics_eps),
    df_counts_nhg     = make_count_table(df_covid_eps[nhggroup == T],
                                          df_ics_eps[nhggroup == T]),
    df_counts_nonnhg  = make_count_table(df_covid_eps[nhggroup == F], 
                                          df_ics_eps[nhggroup == T])
)

# model -------------------------------------------------------------------

extract_results <- function(coefficients) {
    # s <- summary(model)
    # coeff <- s$coefficients
    c <- coefficients[, 1]
    se <- coefficients[, 2]
    LB_95 <- round(exp(c - 1.96 * se), 2) 
    UB_95 <- round(exp(c + 1.96 * se), 2)
    
    
    results <- data.table(
        predictor = rownames(coefficients),
        coefficient = round(c, 2),
        stderror = round(se, 10),
        p = round(coefficients[, 4], 3),
        IRR = round(exp(c), 2),
        CI = paste0(LB_95, ", ", UB_95)
    )
    lapply(results, as.character)
}

fit_poisson <- function(data) {
    model <- glm(
        n_ics_week ~
            week_number +
            after_intervention +
            after_intervention_t +
            Q1 +
            Q2 +
            Q3 + 
            offset(log(n_week)),
        data = data,
        family = quasipoisson)
    return(model)
}

models <- lapply(list_counts, fit_poisson)

lapply(models, function(x) durbinWatsonTest(x, max.lag = 3))
coef <- lapply(models, function(x) coeftest(x, vcov. = NeweyWest(x, lag = 3, prewhite = T)))
results <- lapply(coef, extract_results)


save_results <- function(data, name) {
    filename <- paste0("results/2024_08_12", name, ".csv")
    fwrite(data, file = filename)
}
lapply(
    X = seq_along(results), 
    FUN = function(x) save_results(results[[x]], names(results)[[x]])
)

# sensitivity analysis ----------------------------------------------------

fit_poisson_without_seasons <- function(data) {
    model <- glm(
        n_ics_week ~
            week_number +
            after_intervention +
            after_intervention_t +
            offset(log(n_week)),
        data = data,
        family = quasipoisson)
    return(model)
}
model_no_season <- fit_poisson_without_seasons(list_counts$df_counts_overall)
durbinWatsonTest(model_no_season, max.lag = 3)
coef_no_season <- coeftest(model_no_season, vcov. = NeweyWest(model_no_season, lag = 3, prewhite = T))
results_no_season <- as.data.table(extract_results(coef_no_season))
fwrite(results_no_season, "results/no_season.csv")

# graphs ------------------------------------------------------------------

list_counts$df_counts_overall[, p_ics_week := predict(models$df_counts_overall, 
                                                      type = "response")]
list_counts$df_counts_overall[, p_rate := p_ics_week/n_week]
list_counts$df_counts_overall[, n_ics_week_scaled := as.integer(scaling_coeff*n_ics_week)]
list_counts$df_counts_overall[, rate := n_ics_week/n_week]

CF_data <- copy(list_counts$df_counts_overall)
CF_data[, `:=`(after_intervention = F, after_intervention_t = 0)]
CF_predictions <- predict(models$df_counts_overall, type = "response", newdata = CF_data)
list_counts$df_counts_overall[, p_rate_CF := CF_predictions/n_week]

###


melt(list_counts$df_counts_overall, 
     id.vars = c("week_datum"), 
     measure.vars = c("n_week", "n_ics_week_scaled")) |> 
ggplot(aes(x = week_datum, y = value, color = variable, linetype = variable)) + 
    geom_line() + 
    geom_vline(xintercept = as.Date(intervention_date), color = "purple") +
    scale_y_continuous(
        name = "Consultations, n",
        sec.axis = sec_axis(~./scaling_coeff, name = "ICS prescriptions, n")
    ) + 
    annotate(
        "text",
        x = as.Date(intervention_date)+25,
        y = 8000,
        label = "Guideline revision",
        hjust = "right",
        vjust = "top",
        angle = 270
    ) +
    labs(
        x = "Date"
    ) +
    theme_minimal() +
    theme(
        legend.position = c(0.3, 0.8),
        legend.background = element_rect(fill = "white", color = "black")
    ) +
    scale_color_hue(
        name = "Weekly counts",
        labels = c("n_ics_week_scaled" = "Consultations with ICS prescription", 
                   "n_week" = "Consultations")
    ) +
    scale_linetype(
        name = "Weekly counts",
        labels = c("n_ics_week_scaled" = "Consultations with ICS prescription", 
                   "n_week" = "Consultations")
    )

ggsave(filename = "results/2024_08_12_counts_p_week.png", 
       width = 159, 
       height = 117.7, 
       units = "mm",
       bg = "white")


###
###

make_nice_plot <- function(model, df) {
    
    c <- model$coefficients
    a <- c[2]
    b <- c[1] + 1/4 * (c[5] + c[6] + c[7])
    a2 <- c[4]
    b2 <- c[3]
    
    data <- copy(df)
    
    
    data[, p_ics_week := predict(model, type = "response")]
    data[, p_rate := p_ics_week/n_week]
    data[, n_ics_week_scaled := as.integer(scaling_coeff*n_ics_week)]
    data[, rate := n_ics_week/n_week]
    
    CF_data <- copy(data)
    CF_data[, `:=`(after_intervention = F, after_intervention_t = 0)]
    CF_predictions <- predict(model, type = "response", newdata = CF_data)
    data[, p_rate_CF := CF_predictions/n_week]
    
    data[, trend1_extend := exp(a*week_number + b)]
    data[, trend1 := fifelse(after_intervention == F, exp(a*week_number + b), NA)]
    data[, trend2 := fifelse(after_intervention == T, 
                             exp(a*week_number + a2*after_intervention_t + b2 + b),
                             NA)]
    data[, predictions_pre_int := fifelse(after_intervention == F, p_rate, NA)]
    data[, predictions_post_int := fifelse(after_intervention == T, p_rate, NA)]
    
    labels <- c("1obs" = "Observed prescription rate",
                "2trend" = "Deseasonalised trends",
                "3p" = "Seasonal effect",
                "4cf" = "Counterfactual trend"
    )
    guidelgnd <- guide_legend(title = NULL)
    ggplot(data = data, aes(x = week_datum)) +
        geom_line(aes(y = trend1_extend, color = "4cf", linetype = "4cf", shape = "4cf")) +
        geom_line(aes(y = predictions_pre_int, color = "3p", linetype = "3p", shape = "3p")) +
        geom_line(aes(y = predictions_post_int, color = "3p", linetype = "3p", shape = "3p")) +
        geom_point(aes(y = rate, color = "1obs", linetype = "1obs", shape = "1obs"), size = .5) +
        geom_line(aes(y = trend1, color = "2trend", linetype = "2trend", shape = "2trend")) +
        geom_line(aes(y = trend2, color = "2trend", linetype = "2trend", shape = "2trend")) +
        geom_vline(xintercept = as.Date(intervention_date), color = "purple") +
        scale_color_manual(name = "Legend", labels = labels,
                           values = c("2trend" = "#619CFF", "4cf" = "#619CFF","3p" = "gray", "1obs" = "#F8766D")) + 
        scale_linetype_manual(name = "Legend", labels = labels,
                              values = c("2trend" = "solid", "4cf" = "dashed","3p" = "solid", "1obs" = "blank")) + 
        scale_shape_manual(name = "Legend",  labels = labels,
                           values = c("2trend" = NA, "4cf" = NA,"3p" = NA, "1obs" = 3)) +
        theme_minimal() + 
        theme(
            legend.position = c(0.2,0.85),
            legend.background = element_rect(fill = "white", color = "black")
        ) +
        annotate(
            "text", 
            x = as.Date(intervention_date) + 25, 
            y = 0.030, 
            label = "Guideline revision",
            hjust = "right",
            vjust = "top",
            angle = 270
        ) +
        labs(
            x = "Date",
            y = "Prescription rate",
        ) + 
        ylim(0,0.05) + 
        guides(
            color = guidelgnd,
            shape = guidelgnd,
            linetype = guidelgnd
        )
    
}
make_nice_plot(models$df_counts_overall, list_counts$df_counts_overall)

ggsave(filename = "results/2024_08_12_rate_p_week.png", 
       width = 159, 
       height = 117.7, 
       units = "mm",
       bg = "white")



