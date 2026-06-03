/*******************************************************************************
* Stata Do-File: Labour Market Panel Data Analysis
* Project:       Returns to Education, Experience, and Union Membership
* Author:        [Your Name]
* Date:          2024
* Data:          data/labor_panel.csv (Balanced panel: 20 individuals x 8 years)
* Stata version: 17 or higher recommended
*******************************************************************************/

version 17
clear all
set more off
capture log close

* ── 0. Paths ─────────────────────────────────────────────────────────────────
global root  "."
global data  "$root/data"
global out   "$root/output"
global log   "$root/logs"

cap mkdir "$out"
cap mkdir "$log"

log using "$log/analysis_`c(current_date)'.log", replace text

* ── 1. Load Data ──────────────────────────────────────────────────────────────
import delimited "$data/labor_panel.csv", clear varnames(1)

* ── 2. Describe & Inspect ─────────────────────────────────────────────────────
describe
codebook, compact
list in 1/5

* Check for duplicates in panel identifiers
duplicates report person_id year

* Missing values overview
misstable summarize

* ── 3. Declare Panel Structure ────────────────────────────────────────────────
xtset person_id year, yearly
xtdescribe

* ── 4. Variable Labelling ─────────────────────────────────────────────────────
label variable person_id     "Person identifier"
label variable year          "Survey year"
label variable wage_hourly   "Hourly wage (USD)"
label variable ln_wage       "Log hourly wage"
label variable age           "Age (years)"
label variable education_years "Years of schooling"
label variable experience    "Potential experience (years)"
label variable experience_sq "Potential experience squared"
label variable gender        "Gender"
label variable race          "Race/Ethnicity"
label variable region        "US Census region"
label variable industry      "Industry sector"
label variable occupation    "Occupation group"
label variable employed      "Employed (1=Yes)"
label variable union_member  "Union member (1=Yes)"
label variable firm_size     "Firm size category"
label variable part_time     "Part-time (1=Yes)"
label variable tenure_years  "Job tenure (years)"
label variable hours_week    "Usual hours per week"
label variable health_status "Self-reported health"
label variable married       "Married (1=Yes)"
label variable children      "Number of children (<18) in household"
label variable urban         "Urban residence (1=Yes)"

* Value labels for binary variables
label define yesno 0 "No" 1 "Yes"
foreach v of varlist employed union_member part_time married urban {
    label values `v' yesno
}

* ── 5. Descriptive Statistics ─────────────────────────────────────────────────
* Overall summary (employed observations only)
preserve
keep if employed == 1
summarize wage_hourly ln_wage age education_years experience tenure_years hours_week
restore

* Tabulations
tabulate gender
tabulate race
tabulate region
tabulate industry
tabulate union_member

* Wage by gender
tabulate gender, summarize(wage_hourly)

* Wage by education level
bysort education_years: summarize wage_hourly ln_wage

* Employment rate over time
bysort year: summarize employed

* ── 6. Data Visualisation (via Stata graphs) ──────────────────────────────────
* Wage distribution (log scale)
histogram ln_wage if employed==1, normal xtitle("Log Hourly Wage") ///
    title("Distribution of Log Hourly Wages (2015-2022)") ///
    name(hist_lnwage, replace)
graph export "$out/hist_lnwage.png", replace width(1200)

* Mean wage over time by gender
preserve
keep if employed==1
collapse (mean) mean_wage=wage_hourly, by(year gender)
twoway (line mean_wage year if gender=="Male", lcolor(navy) lwidth(medium)) ///
       (line mean_wage year if gender=="Female", lcolor(cranberry) lwidth(medium) lpattern(dash)), ///
       legend(label(1 "Male") label(2 "Female")) ///
       xtitle("Year") ytitle("Mean Hourly Wage (USD)") ///
       title("Mean Hourly Wage by Gender, 2015-2022") ///
       name(wage_gender_trend, replace)
graph export "$out/wage_gender_trend.png", replace width(1200)
restore

* Wage-experience profile (Mincer curve)
twoway (scatter ln_wage experience if employed==1, mcolor(gs10) msize(small) msymbol(circle_hollow)) ///
       (lfit ln_wage experience if employed==1, lcolor(navy) lwidth(medthick)), ///
       xtitle("Potential Experience (Years)") ytitle("Log Hourly Wage") ///
       legend(off) title("Log Wage–Experience Profile") ///
       name(wage_exp, replace)
graph export "$out/wage_experience.png", replace width(1200)

* ── 7. Cross-Sectional OLS (Mincer Earnings Equation) ────────────────────────
* Baseline: log wage on education and experience
reg ln_wage education_years experience experience_sq ///
    if employed==1, robust

estimates store ols_base

* Full Mincer equation with controls
reg ln_wage education_years experience experience_sq ///
    i.gender i.race i.region i.industry ///
    union_member tenure_years urban married children ///
    if employed==1, robust

estimates store ols_full

* Display comparison table
esttab ols_base ols_full, ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    label title("OLS Mincer Earnings Equations") ///
    mtitles("Baseline" "Full Controls") ///
    stats(N r2_a, labels("N" "Adj. R-squared"))

* ── 8. Panel Data Estimation ──────────────────────────────────────────────────
* Pooled OLS (ignores panel structure)
reg ln_wage education_years experience experience_sq ///
    union_member tenure_years married urban ///
    if employed==1, robust cluster(person_id)
estimates store pooled_ols

* Random Effects (GLS) — assumes individual effects uncorrelated with regressors
xtreg ln_wage education_years experience experience_sq ///
    union_member tenure_years married urban ///
    if employed==1, re

estimates store re_model

* Fixed Effects — within-estimator, controls for time-invariant unobservables
xtreg ln_wage experience experience_sq ///
    union_member tenure_years married urban ///
    if employed==1, fe

* Note: Education and time-invariant characteristics (gender, race) drop from FE
estimates store fe_model

* Hausman test: Fixed vs Random Effects
* H0: RE is consistent (individual effects uncorrelated with regressors)
hausman fe_model re_model

* If p-value < 0.05 → reject H0 → prefer Fixed Effects

* ── 9. First-Difference Estimation ───────────────────────────────────────────
* First-difference removes all time-invariant unobservables
sort person_id year
gen d_ln_wage       = d.ln_wage
gen d_experience    = d.experience
gen d_experience_sq = d.experience_sq
gen d_union         = d.union_member
gen d_tenure        = d.tenure_years
gen d_married       = d.married

reg d_ln_wage d_experience d_experience_sq d_union d_tenure d_married ///
    if employed==1, robust cluster(person_id)
estimates store fd_model

* ── 10. Union Wage Premium ────────────────────────────────────────────────────
* Cross-section union premium
reg ln_wage union_member education_years experience experience_sq ///
    i.gender i.race i.region i.industry ///
    tenure_years urban married children ///
    if employed==1, robust
di "Union wage premium (OLS, %): " (exp(_b[union_member])-1)*100

* Within-person union premium (FE)
xtreg ln_wage union_member experience experience_sq ///
    tenure_years married urban ///
    if employed==1, fe
di "Within-person union premium (FE, %): " (exp(_b[union_member])-1)*100

* ── 11. Gender Wage Gap ───────────────────────────────────────────────────────
* Raw gap
ttest wage_hourly if employed==1, by(gender)

* Adjusted gap (conditional on controls)
reg ln_wage i.gender education_years experience experience_sq ///
    i.race i.region i.industry ///
    union_member tenure_years urban married children ///
    if employed==1, robust
di "Adjusted gender gap (%, Female vs Male): " (exp(_b[2.gender])-1)*100

* ── 12. COVID-19 Shock: Employment Impact (2020) ─────────────────────────────
gen covid = (year == 2020)

* Employment probability — Linear Probability Model with individual FE
xtreg employed covid married children age, fe robust
estimates store emp_lpm_fe

* ── 13. Export Results ────────────────────────────────────────────────────────
* Regression table to Excel
esttab pooled_ols re_model fe_model fd_model ///
    using "$out/panel_results.csv", ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    label title("Panel Earnings Regressions") ///
    mtitles("Pooled OLS" "Random Effects" "Fixed Effects" "First Difference") ///
    stats(N r2_o r2_w, labels("N" "R-sq Overall" "R-sq Within")) ///
    replace

* ── 14. Save Processed Dataset ───────────────────────────────────────────────
save "$data/labor_panel_cleaned.dta", replace

log close
*******************************************************************************/
