# Local verification

Checked on 27 October 2026 using R 4.6.0 (Windows). This records the observed environment; it does not adopt the R 4.6.1 / RStudio version stated in the source Word file.

## Observed reproduction

The two local Birth Recode Stata files were read with bundled Python/pandas, restricted to the required numeric variables, and converted to private temporary CSVs. The delivered R script was then run end to end on those CSVs. Input-row data and temporary CSVs are excluded from the public package. No results were tuned to match the manuscript.

| Measure | BDHS 2017-18 | BDHS 2022 |
|---|---:|---:|
| Raw Birth Recode records | 47,828 | 64,722 |
| Aged 6-35 months | 4,302 | 6,600 |
| Valid HAZ among age-eligible | 3,906 | 2,048 |
| Final analytic sample | 3,423 | 1,783 |
| Stunted, unweighted n | 1,115 | 420 |
| Weighted stunting prevalence | 31.8148% | 22.9439% |
| Early ANC group, n | 1,426 | 726 |
| Later ANC group, n | 1,997 | 1,057 |
| Weighted stunting prevalence, early ANC | 26.8642% | 20.9220% |
| Weighted stunting prevalence, later ANC | 35.1515% | 24.2789% |
| Adjusted PR, early vs later | 0.8956374 | 0.9809018 |
| Normal Wald 95% CI | 0.7992526-1.003646 | 0.8012519-1.200831 |

Primary PRs and confidence intervals match the supplied draft at two decimal places: 0.90 (0.80-1.00) and 0.98 (0.80-1.20). The 2017-18 upper limit is above 1 before rounding; the rounded lower PR should not be described as a statistically clear protective association.

Age-stratified results also match all six draft entries at two decimal places:

| Survey | Age (months) | Model n | Adjusted PR (95% CI) |
|---|---|---:|---|
| 2017-18 | 6-11 | 715 | 0.92 (0.62-1.35) |
| 2017-18 | 12-23 | 1,439 | 0.88 (0.74-1.04) |
| 2017-18 | 24-35 | 1,269 | 0.93 (0.80-1.09) |
| 2022 | 6-11 | 432 | 0.67 (0.39-1.16) |
| 2022 | 12-23 | 699 | 1.03 (0.78-1.36) |
| 2022 | 24-35 | 652 | 1.15 (0.82-1.60) |

## Repairs and audit findings

1. Replaced logistic OR models with log-link Poisson PR models to follow the revised manuscript.
2. Excluded invalid HW70 codes before defining stunting. The original broad `hw70 >= -200` condition would classify nonmissing special codes as not stunted. There were 34 such age-eligible records in 2017-18 and 10 in 2022; their final inclusion also depends on other eligibility criteria.
3. Retained full-frame PSU information for domain variance and exported an analytic-frame variance sensitivity check.
4. Removed pooled interaction analyses that were not carried into the revised manuscript. Separate survey fits avoid cross-round PSU-ID collisions in the old pooled design.
5. Exported exact p-values, unrounded PRs, normal and t inference separately, and recorded convergence, fitted values above one and age-model failures.
6. Found one age-eligible M13=0 record in 2017-18 and none in 2022. It is excluded under the draft's 1-9-month definition. Its meaning must be established from survey-specific documentation before changing that definition.

## Tests actually run

- Full end-to-end R execution on both private numeric CSV extracts, including all primary/age models, aggregate CSVs and PDF/PNG export.
- Synthetic HAZ boundary and invalid-code checks, ANC category checks, zero-weight exclusion, equality of crude model PR to the directly weighted prevalence ratio, invariance to common weight scaling, an independently derived two-group PSU influence-function variance, and a domain fit retaining non-domain PSUs. All passed.
- Comparison of primary sample sizes, weighted prevalence, primary ANC PR/CI and all six age-specific ANC PR/CIs against the supplied manuscript. These selected values match at the reported precision. This is not a full independent audit of every manuscript statement or reference.

## Limits of verification

`haven` and `survey` were not installed in the available R library. Therefore the R script's direct Stata import and optional `survey::svyglm` comparison were not executed. The numeric modelling path was executed in R on private extracts read directly from the Stata files with pandas. The explicit sandwich implementation was checked against an independently derived synthetic two-group variance; no external package equivalence is claimed as an executed test.

The original draft does not specify all finite-sample conventions. This implementation declares normal Wald inference, within-stratum correction and full-frame domains explicitly. Rounded agreement does not establish every computational detail of the original analysis.

The script uses the manuscript's complete-case design and covariate list; reproducing it does not establish causal validity. Included/excluded summaries are descriptive. Selection-bias correction, unmeasured-confounding sensitivity and a clinical/scientific reconsideration of confounder choices remain separate tasks.
