# BDHS ANC timing and childhood stunting

Reproducible R analysis for children aged 6-35 months in Bangladesh DHS 2017-18 and 2022. The code follows the supplied manuscript's prevalence-ratio analysis, rather than the logistic/odds-ratio analysis in the original Word code file.

## Files

- `scripts/01_anc_stunting.R`: cleaning, sample flow, descriptive summaries, crude and adjusted Poisson models, age-stratified models, inference, manuscript comparison and forest plot.
- `tests/test_core.R`: runnable tests using synthetic data only.
- `VERIFICATION.md`: results of the local verification and material limitations.
- `.gitignore`: excludes private data, local output and session files.

Raw DHS data, participant-level extracts and manuscript documents are not included.

## Get the data

Apply for access through the [DHS Program](https://dhsprogram.com/data/available-datasets.cfm). Obtain the Birth Recode files `BDBR7RFL.DTA` and `BDBR81FL.DTA`. Keep them outside this repository. Access and redistribution remain subject to DHS conditions.

Required numeric variables: `hw70`, `m13`, `b19`, `b4`, `v012`, `v106`, `v190`, `v025`, `v005`, `v021`, `v022`.

## Run in RStudio

Open the downloaded repository folder as your working directory. Install the Stata reader once if needed:

```r
install.packages('haven')
```

Then run:

```r
source('scripts/01_anc_stunting.R')
main(c('D:/BDHS/BDBR7RFL.DTA', 'D:/BDHS/BDBR81FL.DTA', 'results'))
```

Replace the two input paths with your own. A private numeric CSV with the same variables is also accepted. The script does not install packages automatically or save participant-level data.

Command-line alternative:

```text
Rscript scripts/01_anc_stunting.R D:/BDHS/BDBR7RFL.DTA D:/BDHS/BDBR81FL.DTA results
```

Run synthetic tests from the repository root:

```r
source('tests/test_core.R')
```

or `Rscript tests/test_core.R`.

## Definitions and modelling choices

- Ages 6-35 months, inclusive; valid HW70 within -600 to 600. HW70 is HAZ times 100; stunting is HW70 < -200. Missing/special anthropometry codes are excluded. No observations were at either valid-range endpoint in the checked inputs.
- First ANC in months 1-3 is early; months 4-9 is later (reference). These reproduce the manuscript definition. Other codes, including M13=0, are excluded and audited. Do not interpret M13=0 as absence of ANC without checking the survey-specific label/questionnaire.
- Complete cases for child age, sex, maternal age (15-49), maternal education, wealth, residence and valid sampling information. The comparison does not include all pregnancies without ANC. Prevalence describes the restricted analytical sample, not the national under-five population.
- Weights V005/1,000,000, normalized to mean one within each fitted model. V021 denotes PSU and V022 stratum. Confirm this stratum choice against each survey's documentation before extending the analysis.
- Poisson estimating equations with a log link and stratified PSU sandwich variance. Within-stratum centering and the m/(m-1) correction are implemented explicitly in base R. Poisson model-based standard errors are never used for inference. A finite-population correction is not assumed.
- The default full valid-design Birth Recode frame preserves zero-score PSUs outside the analytic domain. The alternative analytic-record PSU frame is exported as a variance sensitivity check. Their rounded primary confidence intervals agree in these data.
- Primary models adjust for ANC timing, continuous child age, child sex, current maternal age, maternal education, wealth and residence. Current maternal age is not age at delivery.
- Age-stratified models use 6-11, 12-23 and 24-35 months and omit continuous child age, following the manuscript; this does not remove possible residual confounding by age within each band.
- Primary inference uses explicitly specified large-sample normal Wald intervals/p-values, consistent with the manuscript's sandwich approach. Design-df t intervals/p-values are exported separately, not mixed with normal intervals.
- Singleton strata trigger an error; they are not silently adjusted. Rank deficiency and convergence failures also trigger an error. Failed age models are flagged with their reason.
- If `survey` is installed, the script additionally checks each model's coefficients and covariance against `survey::svyglm` using the same domain. This optional check was not executable in the verification environment because `survey` was absent.
- Pooled logistic interactions from the old script are not used to support the revised manuscript's PR findings. No causal effect, IPW correction, E-value or unmeasured-confounding adjustment is claimed.

## Outputs

The results directory contains CSV equivalents of Tables 1-6, all crude/adjusted coefficients, age-model diagnostics, data-cleaning audits, included/excluded descriptive counts, an analytic-frame variance sensitivity table, draft-manuscript comparison, Figure 1 as PDF/PNG and the R session information. Table 2 is a compact sample summary matching the supplied draft, not an exhaustive covariate balance table. Outputs are computed from the input files and are not inserted from manuscript values.

## Method references

- [DHS-7 recode map: HW70 range and special codes](https://www.dhsprogram.com/pubs/pdf/DHSG4/Recode7_Map_31Aug2018_DHSG4.pdf).
- [Survey-weighted GLMs and model-robust variance](https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/svyglm.html).
- [Survey domain/subpopulation analysis](https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/subset.survey.design.html).

## Scope

This repository contains analysis code. Publicly sharing it does not change the unpublished manuscript's publication status or make DHS participant data public. Review aggregate results and coauthor requirements before adding additional files to GitHub.
