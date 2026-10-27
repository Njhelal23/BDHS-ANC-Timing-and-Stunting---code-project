source('scripts/01_anc_stunting.R')
# Synthetic data only; no DHS data are bundled with these tests.
set.seed(20261027)
n <- 1200L
raw <- data.frame(hw70=sample(c(-350,-201,-200,0,200),n,TRUE),
 m13=sample(1:9,n,TRUE),b19=sample(6:35,n,TRUE),b4=sample(1:2,n,TRUE),
 v012=sample(18:42,n,TRUE),v106=sample(0:3,n,TRUE),v190=sample(1:5,n,TRUE),
 v025=sample(1:2,n,TRUE),v005=runif(n,2e5,2e6),
 v021=rep(1:40,each=30),v022=rep(1:4,each=300))
edge <- raw[1:9,]
edge$hw70 <- c(-601,-600,-201,-200,600,601,9996,9998,NA)
z <- prepare(edge,'synthetic')
stopifnot(identical(z$stunted,c(NA_integer_,1L,1L,0L,0L,NA_integer_,NA_integer_,NA_integer_,NA_integer_)))
edge <- raw[1:6,];edge$m13 <- c(0,1,3,4,9,98)
z <- prepare(edge,'synthetic')
stopifnot(identical(as.character(z$anc_timing),c(NA,'Early','Early','Later','Later',NA)))
x <- prepare(raw,'synthetic');rownames(x) <- seq_len(nrow(x))
f <- stunted~anc_timing
m <- fit_pr(x,f,x)
# For the two-group unadjusted model, PR equals weighted prevalence ratio.
group_mean <- function(g) with(x[x$anc_timing==g,],sum(weight*stunted)/sum(weight))
stopifnot(abs(m$table$PR[m$table$term=='anc_timingEarly']-group_mean('Early')/group_mean('Later'))<1e-7)
# Multiplying every weight by a constant cannot alter coefficients or covariance.
x2 <- x;x2$weight <- x$weight*100
m2 <- fit_pr(x2,f,x2)
stopifnot(isTRUE(all.equal(m$vcov,m2$vcov,tolerance=1e-9)))
# Direct cluster calculation checks the sandwich for the log prevalence ratio.
muE <- group_mean('Early');muL <- group_mean('Later')
wE <- sum(x$weight[x$anc_timing=='Early']);wL <- sum(x$weight[x$anc_timing=='Later'])
inf <- with(x,ifelse(anc_timing=='Early',weight*(stunted-muE)/(wE*muE),
                    -weight*(stunted-muL)/(wL*muL)))
ps <- aggregate(inf,list(h=x$v022,j=x$v021),sum)
expected <- sum(vapply(split(ps$x,ps$h),function(u) length(u)/(length(u)-1)*sum((u-mean(u))^2),numeric(1)))
stopifnot(abs(m$vcov[2,2]-expected)<1e-8)
# Domains retain zero-score PSUs; estimate remains the direct prevalence ratio.
dom <- x[x$v021<=35,]
md <- fit_pr(dom,f,x)
stopifnot(all(is.finite(md$table$SE)))
# Invalid design weight must not enter the analysis.
raw$v005[1] <- 0
stopifnot(!prepare(raw,'synthetic')$included[1])
cat('PASS: HAZ boundaries, special codes, ANC coding, direct prevalence ratio,\n')
cat('weight invariance, independently derived cluster variance, domain fit, invalid weight.\n')
