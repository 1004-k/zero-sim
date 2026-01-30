# ============================================================
# ZERO / Drug Safety simulation (extended grid) — reviewer-proof add-ons
#   A) PS refit sensitivity in confirmed-use cohort (Spec B)
#   B) Target-population shift quantification (initiators vs confirmed)
#   C) ATE weight truncation sensitivity (q = 0.99 by default)
#
# Grid: tau2nd (30/60/90) x p2nd_delta (0.05/0.15/0.30)
# Estimands: ATE (IPTW) vs overlap (ATO)
# Specs:
#   - A_init (time zero at initiation)
#   - B_confirm (time zero at 2nd dispensing; delayed entry; restricted)
# PS handling:
#   - fixedPS: PS fit in A_init, applied to both A_init and B_confirm
#   - refitPS: PS re-fit within B_confirm (sensitivity)
#
# Outputs (written to current working directory):
#   1) perf_ZERO_ext_refitPS_trunc_B{B}_N{N}.csv
#   2) balance_summary_ZERO_ext_refitPS_trunc_B{B}_N{N}.csv
#   3) target_shift_summary_ZERO_ext_B{B}_N{N}.csv
#   4) target_shift_raw_ZERO_ext_B{B}_N{N}.csv
#   5) sessionInfo_ZERO_ext_refitPS_trunc.txt
#   6) sim_progress_ZERO_ext_refitPS_trunc.log
# ============================================================

# ---------------------------
# Packages
# ---------------------------
req <- c("data.table","survival","future.apply","future")
# cobalt is optional (only if you want bal.tab). Default uses fast manual SMD.
req_opt <- c("cobalt")

to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if(length(to_install) > 0) library(data.table)
library(survival)
library(future)
library(future.apply)

USE_COBALT <- FALSE # More faster than "Ture"
if (USE_COBALT) {
  to_install2 <- req_opt[!sapply(req_opt, requireNamespace, quietly = TRUE)]
  if(length(to_install2) > 0) library(cobalt)
}

set.seed(1)

# ---------------------------
# PARALLEL
# ---------------------------
n_cores <- as.integer(Sys.getenv("N_CORES","3"))
future::plan(future::multisession, workers = n_cores)
cat("Using", n_cores, "workers.\n")

# ---------------------------
# LOG FILE (local path)
# ---------------------------
log_dir <- if (.Platform$OS.type == "windows") "C:/temp" else tempdir()
out_dir <- Sys.getenv("OUT_DIR", file.path("output","raw"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, "sim_progress_ZERO_ext_refitPS_trunc.log")
if (file.exists(log_file)) file.remove(log_file)

log_line <- function(msg){
  cat(msg, file = log_file, append = TRUE)
}
log_line(sprintf("[%s] ZERO ext (refitPS + targetShift + trunc) started\n",
                 format(Sys.time(), "%H:%M:%S")))

# ---------------------------
# SIMULATION KNOBS
# ---------------------------
B <- as.integer(Sys.getenv("B","500"))
N <- as.integer(Sys.getenv("N","50000"))

# scenario grid
tau2nd_grid      <- c(30, 60, 90)
p2nd_delta_grid  <- c(0.05, 0.15, 0.30)
p2nd_base        <- 0.75

# DGM parameters
beta_true     <- 0                 # causal null
lambda_event  <- 0.00008
lambda_death  <- 0.00002
followup_days <- 365

# truncation sensitivity for ATE weights
ATE_TRUNC_Q <- 0.99

# numeric stability (avoid div-by-0)
PS_EPS <- 1e-6

# ---------------------------
# Utility
# ---------------------------
ess <- function(w) (sum(w)^2) / sum(w^2)

cap_weights <- function(w, q=0.99){
  cap <- as.numeric(stats::quantile(w, probs=q, na.rm=TRUE, type=7))
  pmin(w, cap)
}

summarize_perf <- function(beta_hat, se_hat, beta_true){
  z  <- qnorm(0.975)
  lo <- beta_hat - z*se_hat
  hi <- beta_hat + z*se_hat
  list(
    bias   = mean(beta_hat - beta_true, na.rm=TRUE),
    rmse   = sqrt(mean((beta_hat - beta_true)^2, na.rm=TRUE)),
    cover  = mean(lo <= beta_true & beta_true <= hi, na.rm=TRUE),
    # under null, P(beta_hat<0)
    sign_rev = mean(beta_hat < 0, na.rm=TRUE),
    # under null, CI exclusion = apparent signal rate (type I error)
    fpr_ci = mean(!(lo <= beta_true & beta_true <= hi), na.rm=TRUE)
  )
}

wmean <- function(x, w) sum(w*x, na.rm=TRUE) / sum(w, na.rm=TRUE)
wvar <- function(x, w){
  mu <- wmean(x, w)
  sum(w*(x-mu)^2, na.rm=TRUE) / sum(w, na.rm=TRUE)
}

weighted_smd_one <- function(x, A, w){
  # A is 0/1; w are analysis weights
  w1 <- w[A==1]; w0 <- w[A==0]
  x1 <- x[A==1]; x0 <- x[A==0]
  mu1 <- wmean(x1, w1); mu0 <- wmean(x0, w0)
  v1  <- wvar(x1, w1);  v0  <- wvar(x0, w0)
  sp  <- sqrt((v1 + v0)/2)
  if (!is.finite(sp) || sp == 0) return(0)
  (mu1 - mu0) / sp
}

balance_fast <- function(d, w, vars=c("age","sex","bmi","egfr","util","gall")){
  smds <- sapply(vars, function(v) weighted_smd_one(d[[v]], d$A, w))
  smd_abs <- abs(smds)
  data.table(
    pct_smd_le_0.1 = mean(smd_abs <= 0.1, na.rm=TRUE),
    max_smd        = max(smd_abs, na.rm=TRUE),
    ess            = ess(w),
    ess_over_n     = ess(w) / nrow(d),
    w_max          = max(w, na.rm=TRUE),
    w_med          = median(w, na.rm=TRUE),
    w_max_over_med = max(w, na.rm=TRUE) / median(w, na.rm=TRUE)
  )
}

balance_cobalt <- function(d, w, estimand=c("ATE","ATO")){
  estimand <- match.arg(estimand)
  bt <- cobalt::bal.tab(
    A ~ age + sex + bmi + egfr + util + gall,
    data     = d,
    weights  = w,
    method   = "weighting",
    estimand = estimand
  )
  smd <- abs(bt$Balance$Diff.Adj)
  smd <- smd[is.finite(smd)]
  data.table(
    pct_smd_le_0.1 = mean(smd <= 0.1),
    max_smd        = if(length(smd)) max(smd) else NA_real_,
    ess            = ess(w),
    ess_over_n     = ess(w) / nrow(d),
    w_max          = max(w, na.rm=TRUE),
    w_med          = median(w, na.rm=TRUE),
    w_max_over_med = max(w, na.rm=TRUE) / median(w, na.rm=TRUE)
  )
}

# Unweighted SMD (target population shift): confirmed(B) vs initiators(A)
unweighted_smd <- function(xB, xA){
  mB <- mean(xB, na.rm=TRUE); mA <- mean(xA, na.rm=TRUE)
  vB <- stats::var(xB, na.rm=TRUE); vA <- stats::var(xA, na.rm=TRUE)
  sp <- sqrt((vB + vA)/2)
  if (!is.finite(sp) || sp == 0) return(0)
  (mB - mA) / sp
}

# ---------------------------
# 1) Data generation
# ---------------------------
gen_covariates <- function(N){
  data.table(
    id  = 1:N,
    age = pmin(pmax(rnorm(N, 62, 10), 18), 90),
    sex = rbinom(N, 1, 0.48),
    bmi = pmin(pmax(rnorm(N, 31, 6), 18), 60),
    egfr= pmin(pmax(rnorm(N, 75, 18), 10), 120),
    util= rgamma(N, shape=2, rate=0.5),
    gall= rbinom(N, 1, plogis(-2.2 + 0.02*(rnorm(N,62,10)-60)))
  )
}

assign_treatment <- function(dt){
  lin <- with(dt,
              -1.2 + 0.02*(age-60) + 0.08*(bmi-30) - 0.015*(egfr-75) + 0.25*util + 0.8*gall)
  ps <- plogis(lin)
  dt[, A := rbinom(.N, 1, ps)]
  dt[]
}

gen_second_dispense <- function(dt, tau2nd, p2nd_base, p2nd_A_effect){
  lin2 <- with(dt,
               qlogis(p2nd_base) +
                 p2nd_A_effect*A +
                 0.15*util - 0.4*gall +
                 0.03*(bmi-30))
  p2 <- plogis(lin2)
  dt[, R2 := rbinom(.N, 1, p2)]
  dt[, t2 := ifelse(R2==1, sample(7:tau2nd, .N, replace=TRUE), NA_integer_)]
  dt[]
}

gen_event_times <- function(dt, beta_true, lambda_event, lambda_death, followup_days){
  xb <- with(dt,
             0.02*(age-60) + 0.10*(bmi-30) - 0.01*(egfr-75) + 0.25*gall + 0.10*util)
  haz_event <- lambda_event * exp(xb + beta_true*dt$A)
  haz_death <- lambda_death * exp(0.03*(dt$age-60) + 0.05*dt$util)
  
  t_event <- rexp(nrow(dt), rate=haz_event)
  t_death <- rexp(nrow(dt), rate=haz_death)
  
  t_obs <- pmin(t_event, t_death, followup_days)
  status <- ifelse(t_obs==t_event & t_obs < pmin(t_death, followup_days), 1,
                   ifelse(t_obs==t_death & t_obs < pmin(t_event, followup_days), 2, 0))
  
  dt[, `:=`(time=t_obs, status=status)]
  dt[]
}

# ---------------------------
# 2) Specs
# ---------------------------
make_spec_A <- function(dt){
  dA <- copy(dt)
  dA[, `:=`(tstart=0, tstop=time, event=as.integer(status==1))]
  dA
}

make_spec_B <- function(dt){
  dB <- dt[R2==1 & !is.na(t2)]
  dB[, `:=`(
    tstart = t2,
    tstop  = time,
    event  = as.integer(status==1 & time > t2)
  )]
  dB <- dB[tstop > tstart]
  dB[]
}

# ---------------------------
# 3) PS models
# ---------------------------
fit_ps <- function(d){
  glm(A ~ age + sex + bmi + egfr + util + gall,
      data=d, family=binomial())
}

predict_ps <- function(fit, newdata){
  ps <- as.numeric(predict(fit, newdata=newdata, type="response"))
  pmin(pmax(ps, PS_EPS), 1-PS_EPS)
}

# ---------------------------
# 4) Weights
# ---------------------------
get_w <- function(d, ps, estimand=c("ate","overlap"), trunc_q=NULL){
  estimand <- match.arg(estimand)
  if (estimand=="ate") {
    w <- ifelse(d$A==1, 1/ps, 1/(1-ps))
    if (!is.null(trunc_q)) w <- cap_weights(w, q=trunc_q)
  } else {
    w <- ifelse(d$A==1, 1-ps, ps)
  }
  w
}

# ---------------------------
# 5) Outcome model
# ---------------------------
fit_cox <- function(d, w){
  fit <- coxph(Surv(tstart, tstop, event) ~ A + cluster(id),
               data=d, weights=w)
  b  <- unname(coef(fit)[["A"]])
  se <- sqrt(vcov(fit)[["A","A"]])
  c(beta=b, se=se)
}

# ---------------------------
# 6) One replicate
# ---------------------------
run_one <- function(tau2nd, p2nd_delta){
  p2nd_A_effect <- -p2nd_delta
  
  dt <- gen_covariates(N)
  dt <- assign_treatment(dt)
  dt <- gen_second_dispense(dt, tau2nd, p2nd_base, p2nd_A_effect)
  dt <- gen_event_times(dt, beta_true, lambda_event, lambda_death, followup_days)
  
  dA <- make_spec_A(dt)
  dB <- make_spec_B(dt)
  
  # ---- B) Target population shift metrics (unweighted) ----
  cr_A1  <- mean(dt[A==1]$R2==1, na.rm=TRUE)
  cr_A0  <- mean(dt[A==0]$R2==1, na.rm=TRUE)
  cr_all <- mean(dt$R2==1, na.rm=TRUE)
  
  n_init <- nrow(dA)
  n_conf <- nrow(dB)
  
  covars <- c("age","sex","bmi","egfr","util","gall")
  smd_pop <- sapply(covars, function(v){
    unweighted_smd(dB[[v]], dA[[v]])
  })
  
  shift_dt <- data.table(
    n_init = n_init,
    n_confirm = n_conf,
    frac_confirm = n_conf / n_init,
    confirm_rate_all = cr_all,
    confirm_rate_A1 = cr_A1,
    confirm_rate_A0 = cr_A0,
    confirm_rate_diff = cr_A1 - cr_A0,
    max_abs_smd_pop = max(abs(smd_pop), na.rm=TRUE)
  )
  for (v in names(smd_pop)) {
    shift_dt[[paste0("abs_smd_pop_", v)]] <- abs(smd_pop[[v]])
  }
  
  # B cohort empty edge-case
  dB_empty <- (n_conf < 10)
  
  # ---- PS fits ----
  ps_fit_A <- fit_ps(dA)
  ps_A <- predict_ps(ps_fit_A, dA)
  ps_B_fixed <- if (!dB_empty) predict_ps(ps_fit_A, dB) else rep(NA_real_, n_conf)
  
  # ---- A) refit PS within B (sensitivity) ----
  ps_B_refit <- rep(NA_real_, n_conf)
  if (!dB_empty) {
    ps_fit_B <- tryCatch(fit_ps(dB), error=function(e) NULL)
    if (!is.null(ps_fit_B)) {
      ps_B_refit <- tryCatch(predict_ps(ps_fit_B, dB),
                             error=function(e) rep(NA_real_, n_conf))
    }
  }
  
  # ---- Define all analysis specs ----
  specs <- list(
    # A cohort
    A_init_ate            = list(d=dA, ps=ps_A,        estimand="ate",     trunc_q=NULL,        estimand_label="ATE"),
    A_init_ate_trunc      = list(d=dA, ps=ps_A,        estimand="ate",     trunc_q=ATE_TRUNC_Q, estimand_label="ATE"),
    A_init_overlap        = list(d=dA, ps=ps_A,        estimand="overlap", trunc_q=NULL,        estimand_label="ATO"),
    
    # B cohort with fixed PS
    B_confirm_ate_fixedPS        = list(d=dB, ps=ps_B_fixed,  estimand="ate",     trunc_q=NULL,        estimand_label="ATE"),
    B_confirm_ate_fixedPS_trunc  = list(d=dB, ps=ps_B_fixed,  estimand="ate",     trunc_q=ATE_TRUNC_Q, estimand_label="ATE"),
    B_confirm_overlap_fixedPS    = list(d=dB, ps=ps_B_fixed,  estimand="overlap", trunc_q=NULL,        estimand_label="ATO"),
    
    # B cohort with refit PS
    B_confirm_ate_refitPS        = list(d=dB, ps=ps_B_refit,  estimand="ate",     trunc_q=NULL,        estimand_label="ATE"),
    B_confirm_ate_refitPS_trunc  = list(d=dB, ps=ps_B_refit,  estimand="ate",     trunc_q=ATE_TRUNC_Q, estimand_label="ATE"),
    B_confirm_overlap_refitPS    = list(d=dB, ps=ps_B_refit,  estimand="overlap", trunc_q=NULL,        estimand_label="ATO")
  )
  
  est_rows <- list()
  bal_rows <- list()
  
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    
    # handle empty B
    if (grepl("^B_confirm", nm) && dB_empty) {
      est_rows[[nm]] <- data.table(spec=nm, beta=NA_real_, se=NA_real_, n=n_conf, ess=NA_real_)
      bal_rows[[nm]] <- data.table(
        spec=nm, n=n_conf, treat_rate=NA_real_, ess=NA_real_, ess_over_n=NA_real_,
        w_max=NA_real_, w_med=NA_real_, w_max_over_med=NA_real_, pct_smd_le_0.1=NA_real_, max_smd=NA_real_
      )
      next
    }
    
    w <- get_w(sp$d, sp$ps, estimand=sp$estimand, trunc_q=sp$trunc_q)
    
    est <- tryCatch(fit_cox(sp$d, w),
                    error=function(e) c(beta=NA_real_, se=NA_real_))
    
    bal <- tryCatch({
      if (USE_COBALT) {
        balance_cobalt(sp$d, w, estimand=sp$estimand_label)
      } else {
        balance_fast(sp$d, w)
      }
    }, error=function(e){
      data.table(
        pct_smd_le_0.1=NA_real_, max_smd=NA_real_, ess=NA_real_, ess_over_n=NA_real_,
        w_max=NA_real_, w_med=NA_real_, w_max_over_med=NA_real_
      )
    })
    
    est_rows[[nm]] <- data.table(
      spec=nm,
      beta=est[["beta"]],
      se=est[["se"]],
      n=nrow(sp$d),
      ess=ess(w)
    )
    
    bal_rows[[nm]] <- data.table(
      spec=nm,
      n=nrow(sp$d),
      treat_rate=mean(sp$d$A),
      ess=bal$ess,
      ess_over_n=bal$ess_over_n,
      w_max=bal$w_max,
      w_med=bal$w_med,
      w_max_over_med=bal$w_max_over_med,
      pct_smd_le_0.1=bal$pct_smd_le_0.1,
      max_smd=bal$max_smd
    )
  }
  
  list(
    est_dt = rbindlist(est_rows, fill=TRUE),
    bal_dt = rbindlist(bal_rows, fill=TRUE),
    shift_dt = shift_dt
  )
}

# ============================================================
# 7) Run simulation: scenario grid
# ============================================================
scen_grid <- data.table(
  tau2nd     = rep(tau2nd_grid, each = length(p2nd_delta_grid)),
  p2nd_delta = rep(p2nd_delta_grid, times = length(tau2nd_grid))
)
scen_grid[, scenario := sprintf("tau%02d_d%.2f", tau2nd, p2nd_delta)]

cat("Scenario grid:\n")
print(scen_grid)
log_line(sprintf("Scenario grid: %s\n", paste(scen_grid$scenario, collapse=", ")))

all_perf  <- list()
all_bal   <- list()
all_shift_raw <- list()
all_shift_sum <- list()

spec_names <- c(
  "A_init_ate","A_init_ate_trunc","A_init_overlap",
  "B_confirm_ate_fixedPS","B_confirm_ate_fixedPS_trunc","B_confirm_overlap_fixedPS",
  "B_confirm_ate_refitPS","B_confirm_ate_refitPS_trunc","B_confirm_overlap_refitPS"
)

for (i in seq_len(nrow(scen_grid))) {
  tau_i  <- scen_grid$tau2nd[i]
  dlt_i  <- scen_grid$p2nd_delta[i]
  scen_i <- scen_grid$scenario[i]
  
  cat("\n--- Scenario", i, "/", nrow(scen_grid), ":", scen_i, "---\n")
  log_line(sprintf("[%s] Start scenario %s (tau2nd=%d, p2nd_delta=%.2f)\n",
                   format(Sys.time(), "%H:%M:%S"), scen_i, tau_i, dlt_i))
  
  t_start <- proc.time()[3]
  
  rr_list <- future_lapply(
    X = 1:B,
    FUN = function(b){
      rr <- run_one(tau2nd=tau_i, p2nd_delta=dlt_i)
      if (b %% 25 == 0) {
        cat(sprintf("[%s] %s rep %d/%d\n",
                    format(Sys.time(), "%H:%M:%S"), scen_i, b, B),
            file=log_file, append=TRUE)
      }
      rr
    },
    future.seed = TRUE
  )
  
  est_rep   <- rbindlist(lapply(rr_list, `[[`, "est_dt"), fill=TRUE)
  bal_rep   <- rbindlist(lapply(rr_list, `[[`, "bal_dt"), fill=TRUE)
  shift_rep <- rbindlist(lapply(rr_list, `[[`, "shift_dt"), fill=TRUE)
  
  est_rep[,   `:=`(scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i)]
  bal_rep[,   `:=`(scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i)]
  shift_rep[, `:=`(scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i)]
  
  # Performance summary per spec
  perf_i <- rbindlist(lapply(spec_names, function(k){
    dt_k <- est_rep[spec==k]
    s <- summarize_perf(dt_k$beta, dt_k$se, beta_true)
    data.table(
      scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i,
      spec=k,
      bias=s$bias,
      rmse=s$rmse,
      cover=s$cover,
      sign_rev=s$sign_rev,
      fpr_ci=s$fpr_ci,
      ess_median=median(dt_k$ess, na.rm=TRUE),
      n_median=median(dt_k$n, na.rm=TRUE)
    )
  }), fill=TRUE)
  
  # Balance/diagnostics summary per spec
  bal_i <- rbindlist(lapply(spec_names, function(k){
    dt_k <- bal_rep[spec==k]
    data.table(
      spec=k,
      n_median=median(dt_k$n, na.rm=TRUE),
      treat_rate_median=median(dt_k$treat_rate, na.rm=TRUE),
      ess_median=median(dt_k$ess, na.rm=TRUE),
      ess_over_n_median=median(dt_k$ess_over_n, na.rm=TRUE),
      w_max_median=median(dt_k$w_max, na.rm=TRUE),
      w_med_median=median(dt_k$w_med, na.rm=TRUE),
      w_max_over_med_median=median(dt_k$w_max_over_med, na.rm=TRUE),
      pct_smd_le_0.1_median=median(dt_k$pct_smd_le_0.1, na.rm=TRUE),
      max_smd_median=median(dt_k$max_smd, na.rm=TRUE),
      scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i
    )
  }), fill=TRUE)
  
  # Target-population shift summary (scenario-level)
  shift_sum_i <- shift_rep[, {
    out <- list(
      n_init_median = median(n_init, na.rm=TRUE),
      n_confirm_median = median(n_confirm, na.rm=TRUE),
      frac_confirm_median = median(frac_confirm, na.rm=TRUE),
      confirm_rate_all_median = median(confirm_rate_all, na.rm=TRUE),
      confirm_rate_A1_median = median(confirm_rate_A1, na.rm=TRUE),
      confirm_rate_A0_median = median(confirm_rate_A0, na.rm=TRUE),
      confirm_rate_diff_median = median(confirm_rate_diff, na.rm=TRUE),
      max_abs_smd_pop_median = median(max_abs_smd_pop, na.rm=TRUE)
    )
    smd_cols <- grep("^abs_smd_pop_", names(.SD), value=TRUE)
    for (cc in smd_cols) {
      out[[paste0(cc, "_median")]] <- median(get(cc), na.rm=TRUE)
    }
    out
  }]
  shift_sum_i[, `:=`(scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i)]
  
  t_end <- proc.time()[3]
  log_line(sprintf("[%s] End scenario %s | %.2f sec\n",
                   format(Sys.time(), "%H:%M:%S"), scen_i, (t_end - t_start)))
  
  all_perf[[i]] <- perf_i
  all_bal[[i]]  <- bal_i
  all_shift_raw[[i]] <- shift_rep
  all_shift_sum[[i]] <- shift_sum_i
}

perf_all <- rbindlist(all_perf, fill=TRUE)
bal_all  <- rbindlist(all_bal,  fill=TRUE)
shift_raw_all <- rbindlist(all_shift_raw, fill=TRUE)
shift_sum_all <- rbindlist(all_shift_sum, fill=TRUE)

# ============================================================
# 8) Save outputs
# ============================================================
perf_file <- file.path(out_dir, sprintf("perf_ZERO_ext_refitPS_trunc_B%d_N%d.csv", B, N))
bal_file <- file.path(out_dir, sprintf("balance_summary_ZERO_ext_refitPS_trunc_B%d_N%d.csv", B, N))
shift_sum_file <- file.path(out_dir, sprintf("target_shift_summary_ZERO_ext_B%d_N%d.csv", B, N))
shift_raw_file <- file.path(out_dir, sprintf("target_shift_raw_ZERO_ext_B%d_N%d.csv", B, N))

fwrite(perf_all, perf_file)
fwrite(bal_all,  bal_file)
fwrite(shift_sum_all, shift_sum_file)
fwrite(shift_raw_all, shift_raw_file)

sink(file.path(out_dir, "sessionInfo_ZERO_ext_refitPS_trunc.txt"))
print(sessionInfo())
sink()

cat("\nSaved files:\n")
cat(" - ", perf_file, "\n", sep="")
cat(" - ", bal_file, "\n", sep="")
cat(" - ", shift_sum_file, "\n", sep="")
cat(" - ", shift_raw_file, "\n", sep="")
cat(" - sessionInfo_ZERO_ext_refitPS_trunc.txt\n")
cat(" - ", log_file, "\n\n", sep="")

print(perf_all)
print(bal_all)
print(shift_sum_all)