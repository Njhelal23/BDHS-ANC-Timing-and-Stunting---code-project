# BDHS 2017-18 / 2022: first ANC timing and stunting, ages 6-35 months.
# Run from repository root:
# Rscript scripts/01_anc_stunting.R path/to/2017.DTA path/to/2022.DTA results
# Uses base R for modelling; haven is required only to read Stata files.
# Approved private CSV extracts with the same numeric column names are supported.
# This script computes results from data; it never substitutes manuscript values.

required_vars <- c('hw70','m13','b19','b4','v012','v106','v190',
                   'v025','v005','v021','v022')

read_input <- function(path) {
  if (!file.exists(path)) stop('Input not found: ', path)
  if (tolower(tools::file_ext(path)) == 'dta') {
    if (!requireNamespace('haven', quietly=TRUE)) {
      stop('Install haven once: install.packages("haven"). No automatic installation is performed.')
    }
    x <- as.data.frame(haven::read_dta(path))
  } else if (tolower(tools::file_ext(path)) == 'csv') {
    x <- read.csv(path, check.names=FALSE)
  } else stop('Use a .DTA file or an approved private numeric .csv extract.')
  names(x) <- tolower(names(x))
  missing <- setdiff(required_vars, names(x))
  if (length(missing)) stop('Missing columns: ', paste(missing, collapse=', '))
  x <- x[required_vars]
  x[] <- lapply(x, as.numeric)
  x
}

between <- function(x, lo, hi) !is.na(x) & is.finite(x) & x >= lo & x <= hi
prepare <- function(d, year) {
  d$survey <- year
  d$weight <- d$v005 / 1e6
  d$age_eligible <- between(d$b19,6,35)
  # HW70 has two implied decimal places; do not classify special codes as healthy.
  # Endpoints follow the recode range -600:600; endpoint counts are exported.
  d$valid_haz <- between(d$hw70,-600,600)
  d$stunted <- ifelse(d$valid_haz, as.integer(d$hw70 < -200), NA_integer_)
  # Prespecified manuscript comparison: months 1-3 vs 4-9.
  # M13=0 is excluded, not automatically equated with no ANC. Review its label.
  d$anc_timing <- factor(ifelse(d$m13 %in% 1:3,'Early',
                               ifelse(d$m13 %in% 4:9,'Later',NA)),
                         levels=c('Later','Early'))
  d$child_sex <- factor(d$b4,levels=c(1,2),labels=c('Male','Female'))
  d$maternal_education <- factor(d$v106,levels=0:3,
                               labels=c('None','Primary','Secondary','Higher'))
  d$wealth <- factor(d$v190,levels=1:5,
                    labels=c('Poorest','Poorer','Middle','Richer','Richest'))
  d$residence <- factor(d$v025,levels=c(2,1),labels=c('Rural','Urban'))
  d$age_group <- cut(d$b19, breaks=c(6,12,24,36), right=FALSE,
                    labels=c('6-11','12-23','24-35'))
  d$valid_design <- between(d$weight, .Machine$double.xmin, Inf) &
    !is.na(d$v021) & is.finite(d$v021) & !is.na(d$v022) & is.finite(d$v022)
  d$valid_covariates <- between(d$v012,15,49) &
    complete.cases(d[c('child_sex','maternal_education','wealth','residence')])
  d$included <- d$age_eligible & d$valid_haz & !is.na(d$anc_timing) &
    d$valid_covariates & d$valid_design
  d
}

# Stratified ultimate-cluster sandwich, sampling with replacement, no FPC.
# Each PSU score is a sum of row scores; non-domain PSUs contribute zero.
# Within each stratum: m/(m-1) * sum[(U_hj - mean(U_h)) (.)'].
# No independent-observation glm standard errors are used for inference.
psu_meat <- function(scores, rows, frame) {
  skeleton <- unique(frame[c('v022','v021')])
  skeleton <- skeleton[order(skeleton$v022,skeleton$v021),,drop=FALSE]
  key <- function(d) paste(d$v022,d$v021,sep=':')
  idx <- match(key(rows),key(skeleton))
  if (anyNA(idx)) stop('Model PSU absent from survey frame.')
  totals <- matrix(0,nrow(skeleton),ncol(scores))
  sums <- rowsum(scores,idx,reorder=FALSE)
  totals[as.integer(rownames(sums)),] <- sums
  meat <- matrix(0,ncol(scores),ncol(scores))
  dfs <- 0L
  for (h in unique(skeleton$v022)) {
    u <- totals[skeleton$v022 == h,,drop=FALSE]
    m <- nrow(u)
    if (m < 2L) stop('Singleton PSU stratum ',h,
      '. Confirm the sampling design; no automatic lonely-PSU adjustment is applied.')
    centered <- sweep(u,2,colMeans(u),'-')
    meat <- meat + m/(m-1)*crossprod(centered)
    dfs <- dfs + m-1L
  }
  list(meat=meat,design_df=dfs,psus=nrow(skeleton),strata=length(unique(skeleton$v022)))
}

fit_pr <- function(d, formula, frame) {
  if (!nrow(d)) stop('Empty model domain.')
  if (length(unique(d$stunted)) < 2L || length(unique(d$anc_timing)) < 2L)
    stop('Both outcome values and ANC groups are needed.')
  # Check constants before glm; sparse models must not silently change covariates.
  vars <- all.vars(formula)[-1]
  constant <- vars[vapply(d[vars],function(x) length(unique(x)) < 2L,logical(1))]
  if (length(constant)) stop('Constant predictor(s): ',paste(constant,collapse=', '))
  d$w_norm <- d$weight/mean(d$weight)
  fit <- glm(formula,data=d,weights=w_norm,family=quasipoisson(link='log'),
             na.action=na.fail,x=TRUE,y=TRUE,control=glm.control(maxit=100))
  if (!fit$converged || anyNA(coef(fit))) stop('Non-converged or rank-deficient model.')
  X <- fit$x; mu <- fitted(fit); w <- d$w_norm
  bread <- crossprod(X,sweep(X,1,w*mu,'*'))
  if (qr(bread)$rank < ncol(X)) stop('Singular information matrix.')
  scores <- sweep(X,1,w*(d$stunted-mu),'*')
  vc <- psu_meat(scores,d,frame)
  inv <- solve(bread)
  V <- inv %*% vc$meat %*% inv
  se <- sqrt(pmax(diag(V),0)); beta <- coef(fit)
  # Manuscript-style large-sample normal Wald CI/p. Explicitly export t sensitivity.
  z <- qnorm(.975)
  residual_df <- vc$design_df-length(beta)+1L
  critical_t <- if (residual_df > 0) qt(.975,residual_df) else NA_real_
  tab <- data.frame(term=names(beta),log_PR=as.numeric(beta),SE=se,
    PR=exp(beta),CI_lower=exp(beta-z*se),CI_upper=exp(beta+z*se),
    p_value=2*pnorm(-abs(beta/se)),
    CI_lower_t=exp(beta-critical_t*se),CI_upper_t=exp(beta+critical_t*se),
    p_value_t=if(residual_df > 0) 2*pt(-abs(beta/se),residual_df) else NA_real_,
    model_n=nrow(d),design_df=vc$design_df,residual_df=residual_df,
    fitted_above_one=sum(mu>1),row.names=NULL)
  if (any(mu>1)) warning(sum(mu>1),' fitted values >1; inspect model specification.')
  # Optional independent comparison to the established survey package.
  if (requireNamespace('survey',quietly=TRUE)) {
    # Primary-domain cross-check only; valid full-frame data already contain recodes.
    f <- frame
    f$test_domain <- seq_len(nrow(f)) %in% match(rownames(d),rownames(f))
    des <- survey::svydesign(ids=~v021,strata=~v022,weights=~weight,data=f,nest=TRUE)
    dom <- subset(des,test_domain)
    check <- survey::svyglm(formula,design=dom,family=quasipoisson('log'))
    if (!isTRUE(all.equal(unname(coef(check)),unname(beta),tolerance=1e-6)) ||
        !isTRUE(all.equal(unname(vcov(check)),unname(V),tolerance=1e-5)))
      stop('Independent survey-package cross-check failed. Inspect the design.')
  }
  list(table=tab,fit=fit,vcov=V)
}

main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if (length(args)<2L) stop('Usage: Rscript scripts/01_anc_stunting.R 2017.DTA 2022.DTA [results_dir]')
  out <- if(length(args)>=3) args[3] else 'results'
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  emit <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
  primary <- stunted ~ anc_timing + b19 + child_sex + v012 + maternal_education + wealth + residence
  age_formula <- stunted ~ anc_timing + child_sex + v012 + maternal_education + wealth + residence
  years <- c('2017-18','2022')
  all_flow <- all_summary <- all_anc <- all_models <- all_age <- all_audit <- list()
  all_sensitivity <- all_diagnostics <- list()
  for (i in seq_along(years)) {
    year <- years[i]; d <- prepare(read_input(args[i]),year)
    # Stable row identifiers permit domain matching in the optional survey check.
    rownames(d) <- seq_len(nrow(d))
    frame <- d[d$valid_design,,drop=FALSE]
    a <- d[d$included,,drop=FALSE]
    if (!nrow(a)) stop('No analytic records for ',year)
    flow <- c(raw_records=nrow(d),age_6_35=sum(d$age_eligible),
      age_and_valid_HAZ=sum(d$age_eligible & d$valid_haz),
      age_HAZ_and_ANC=sum(d$age_eligible & d$valid_haz & !is.na(d$anc_timing)),
      complete_covariates=sum(d$age_eligible & d$valid_haz & !is.na(d$anc_timing) & d$valid_covariates),
      final_analytical_N=nrow(a))
    all_flow[[year]] <- data.frame(survey=year,step=names(flow),N=as.integer(flow))
    wmean <- function(x,w) sum(x*w)/sum(w)
    all_summary[[year]] <- data.frame(survey=year,N=nrow(a),
      stunted_n=sum(a$stunted),stunting_unweighted_percent=100*mean(a$stunted),
      stunting_weighted_percent=100*wmean(a$stunted,a$weight),
      early_ANC_unweighted_percent=100*mean(a$anc_timing=='Early'),
      early_ANC_weighted_percent=100*wmean(a$anc_timing=='Early',a$weight),
      child_age_weighted_mean=wmean(a$b19,a$weight),maternal_age_weighted_mean=wmean(a$v012,a$weight))
    all_anc[[year]] <- do.call(rbind,lapply(levels(a$anc_timing),function(g) {
      z <- a[a$anc_timing==g,,drop=FALSE]
      data.frame(survey=year,ANC=g,N=nrow(z),stunted_n=sum(z$stunted),
                 weighted_stunting_percent=100*wmean(z$stunted,z$weight))
    }))
    all_audit[[year]] <- data.frame(survey=year,
      invalid_nonmissing_HAZ_in_age=sum(d$age_eligible & !d$valid_haz & !is.na(d$hw70)),
      HAZ_exact_endpoints_in_age=sum(d$age_eligible & d$hw70 %in% c(-600,600)),
      M13_zero_in_age=sum(d$age_eligible & d$m13 %in% 0),
      invalid_design_in_age=sum(d$age_eligible & !d$valid_design))
    emit(as.data.frame(table(M13=d$m13[d$age_eligible],useNA='always')),paste0('audit_M13_',year))
    # Included/excluded comparison is descriptive, not a selection-bias correction.
    eligible <- d[d$age_eligible,,drop=FALSE]
    all_diagnostics[[year]] <- do.call(rbind,lapply(c(FALSE,TRUE),function(flag) {
      z <- eligible[eligible$included==flag,,drop=FALSE]
      data.frame(survey=year,included=flag,N=nrow(z),valid_HAZ_N=sum(z$valid_haz),
        ANC_known_N=sum(!is.na(z$anc_timing)),maternal_age_mean=mean(z$v012,na.rm=TRUE),
        child_age_mean=mean(z$b19,na.rm=TRUE))
    }))
    adj <- fit_pr(a,primary,frame)$table; adj$survey <- year; adj$model <- 'Adjusted'
    crude <- fit_pr(a,stunted~anc_timing,frame)$table; crude$survey <- year; crude$model <- 'Crude'
    all_models[[year]] <- rbind(crude,adj)
    # Explicit variance sensitivity: build the PSU frame from analytic records only.
    # The source manuscript does not fully specify its treatment of non-domain PSUs.
    sens <- fit_pr(a,primary,a)$table; sens$survey <- year
    all_sensitivity[[year]] <- sens
    age_results <- list()
    for (g in levels(a$age_group)) {
      z <- a[a$age_group==g,,drop=FALSE]
      tryCatch({
        result <- fit_pr(z,age_formula,frame)$table
        result$survey <- year;result$age_group <- g;result$status <- 'estimated'
        age_results[[g]] <- result
      },error=function(e) {
        warning('Age group ',g,' / ',year,': ',conditionMessage(e))
        age_results[[g]] <<- data.frame(term='anc_timingEarly',log_PR=NA_real_,SE=NA_real_,
          PR=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,p_value=NA_real_,
          CI_lower_t=NA_real_,CI_upper_t=NA_real_,p_value_t=NA_real_,
          model_n=nrow(z),design_df=NA_integer_,residual_df=NA_integer_,
          fitted_above_one=NA_integer_,survey=year,age_group=g,status=conditionMessage(e))
      })
    }
    all_age[[year]] <- do.call(rbind,age_results)
  }
  flow <- do.call(rbind,all_flow); desc <- do.call(rbind,all_summary)
  anc <- do.call(rbind,all_anc); models <- do.call(rbind,all_models)
  ages <- do.call(rbind,all_age); sens <- do.call(rbind,all_sensitivity)
  emit(flow,'table1_sample_flow');emit(desc,'table2_sample_summary');emit(anc,'table3_ANC_prevalence')
  emit(models[models$model=='Adjusted',],'table4_adjusted_PR')
  emit(ages[ages$survey=='2017-18' & ages$term=='anc_timingEarly',],'table5_age_PR_2017_18')
  emit(ages[ages$survey=='2022' & ages$term=='anc_timingEarly',],'table6_age_PR_2022')
  emit(models,'all_crude_adjusted_models');emit(ages,'all_age_model_coefficients')
  emit(sens,'sensitivity_analytic_frame_variance')
  emit(do.call(rbind,all_audit),'data_cleaning_audit')
  emit(do.call(rbind,all_diagnostics),'included_excluded_descriptive')
  # Published draft values are comparison targets, never used in fitting.
  observed <- models[models$model=='Adjusted' & models$term=='anc_timingEarly',]
  targets <- data.frame(survey=years,manuscript_N=c(3423,1783),
    manuscript_PR=c(.90,.98),manuscript_lower=c(.80,.80),manuscript_upper=c(1.00,1.20))
  comparison <- merge(observed,targets,by='survey')
  comparison$N_matches <- comparison$model_n==comparison$manuscript_N
  comparison$PR_matches_2dp <- round(comparison$PR,2)==comparison$manuscript_PR
  comparison$CI_matches_2dp <- round(comparison$CI_lower,2)==comparison$manuscript_lower &
    round(comparison$CI_upper,2)==comparison$manuscript_upper
  emit(comparison,'manuscript_comparison')
  plot_forest <- function() {
    old_par <- par(mar=c(5,6,4,2)+0.1)
    on.exit(par(old_par))
    z <- observed[match(years,observed$survey),]
    limits <- range(c(z$CI_lower,z$CI_upper,1),finite=TRUE)
    plot(z$PR,2:1,log='x',xlim=limits*c(.9,1.1),ylim=c(.5,2.5),yaxt='n',
      xlab='Adjusted prevalence ratio (log scale)',ylab='',pch=19,
      main='Early versus later ANC and childhood stunting')
    axis(2,at=2:1,labels=years,las=1)
    segments(z$CI_lower,2:1,z$CI_upper,2:1);abline(v=1,lty=2,col='grey50')
    mtext('Stratified PSU sandwich; normal Wald 95% CI',side=3,cex=.75)
  }
  pdf(file.path(out,'figure1_adjusted_PR.pdf'),width=7,height=4);plot_forest();dev.off()
  png(file.path(out,'figure1_adjusted_PR.png'),width=1800,height=1000,res=220);plot_forest();dev.off()
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  message('Analysis complete. Results: ',normalizePath(out))
  print(comparison[c('survey','model_n','PR','CI_lower','CI_upper','PR_matches_2dp','CI_matches_2dp')])
  invisible(list(flow=flow,summary=desc,models=models,ages=ages,comparison=comparison))
}

if (sys.nframe()==0L) main()
