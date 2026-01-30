# ============================================================
# ZERO / SER simulation (extended grid + parallel + logging)
#  - Grid: tau2nd (30/60/90) x p2nd_delta (0.05/0.15/0.30)
#  - Fixed PS fit on Spec A, applied to Spec B (same as original)
#  - IPW: ATE vs overlap
#  - Parallel: 3 cores via future.apply
# Outputs:
#   1) perf_ZERO_ext_B{B}_N{N}.csv
#   2) balance_summary_ZERO_ext_B{B}_N{N}.csv
#   3) sessionInfo_ZERO_ext.txt
#   4) sim_progress_ZERO_ext.log
# ============================================================

# ---- Packages ----
req <- c("data.table","survival","cobalt","WeightIt","future.apply","future")
to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if(length(to_install) > 0){
  }
library(data.table)
library(survival)
library(cobalt)
library(WeightIt)      # loaded for workflow consistency; PS computed via glm here
library(future)
library(future.apply)

set.seed(1)

# ---- PARALLEL ----
n_cores <- as.integer(Sys.getenv("N_CORES","3"))
future::plan(future::multisession, workers = n_cores)
cat("Using", n_cores, "workers.\n")

# ---- LOG FILE (force local path, avoid OneDrive) ----
log_dir <- Sys.getenv("LOG_DIR", file.path("output","logs")) 
out_dir <- Sys.getenv("OUT_DIR", file.path("output","raw"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

log_file <- file.path(log_dir, "sim_progress_ZERO_ext.log")

if (file.exists(log_file)) file.remove(log_file)

cat(sprintf("[%s] ZERO extended simulation started\n",
            format(Sys.time(), "%H:%M:%S")),
    file = log_file, append = TRUE)

log_line <- function(msg){
  cat(msg, file = log_file, append = TRUE)
}

# ---- SIMULATION KNOBS (edit here if needed) ----
B <- as.integer(Sys.getenv("B","500"))
N <- as.integer(Sys.getenv("N","50000"))

# scenario grid
tau2nd_grid       <- c(30, 60, 90)
p2nd_delta_grid   <- c(0.05, 0.15, 0.30)  # imbalance magnitude
p2nd_base         <- 0.75

beta_true    <- 0                 # null DGM: "artifact masquerading as signal"
lambda_event <- 0.00008
lambda_death <- 0.00002
followup_days <- 365

# Optional: weight truncation (sensitivity-like)
do_truncate <- FALSE
w_cap_q     <- 0.99

# -----------------------------------------------

# ---- Utility ----
ess <- function(w){
  (sum(w)^2) / sum(w^2)
}

truncate_w <- function(w, q=0.99){
  cap <- as.numeric(stats::quantile(w, probs=q, na.rm=TRUE))
  pmin(w, cap)
}

summarize_perf <- function(beta_hat, se_hat, beta_true){
  z <- qnorm(0.975)
  lo <- beta_hat - z*se_hat
  hi <- beta_hat + z*se_hat
  list(
    bias     = mean(beta_hat - beta_true, na.rm=TRUE),
    rmse     = sqrt(mean((beta_hat - beta_true)^2, na.rm=TRUE)),
    cover    = mean(lo <= beta_true & beta_true <= hi, na.rm=TRUE),
    sign_rev = mean(beta_hat < 0, na.rm=TRUE),
    fpr_ci   = mean(!(lo <= 0 & 0 <= hi), na.rm=TRUE)
  )
}

# ---- 1) Data generation ----
gen_covariates <- function(N){
  dt <- data.table(
    id  = 1:N,
    age = pmin(pmax(rnorm(N, 62, 10), 18), 90),
    sex = rbinom(N, 1, 0.48),
    bmi = pmin(pmax(rnorm(N, 31, 6), 18), 60),
    egfr= pmin(pmax(rnorm(N, 75, 18), 10), 120),
    util= rgamma(N, shape=2, rate=0.5),
    gall= rbinom(N, 1, plogis(-2.2 + 0.02*(rnorm(N,62,10)-60)))
  )
  dt[]
}

assign_treatment <- function(dt){
  lin <- with(dt, -1.2 + 0.02*(age-60) + 0.08*(bmi-30) - 0.015*(egfr-75) + 0.25*util + 0.8*gall)
  ps  <- plogis(lin)
  dt[, A := rbinom(.N, 1, ps)]   # 1=GLP-1RA, 0=SGLT2i
  dt[]
}

gen_second_dispense <- function(dt, tau2nd, p2nd_base, p2nd_A_effect){
  lin2 <- with(dt,
               qlogis(p2nd_base) +
                 p2nd_A_effect*A +
                 0.15*util - 0.4*gall +
                 0.03*(bmi-30)
  )
  p2 <- plogis(lin2)
  dt[, R2 := rbinom(.N, 1, p2)]
  dt[, t2 := ifelse(R2==1, sample(7:tau2nd, .N, replace=TRUE), NA_integer_)]
  dt[]
}

gen_event_times <- function(dt, beta_true, lambda_event, lambda_death, followup_days){
  xb <- with(dt,
             0.02*(age-60) + 0.10*(bmi-30) - 0.01*(egfr-75) + 0.25*gall + 0.10*util
  )
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

# ---- 2) Specs (time-zero) ----
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
  dB
}

# ---- 3) PS model (fixed) and weights (ATE vs overlap) ----
fit_ps_A <- function(dA){
  glm(A ~ age + sex + bmi + egfr + util + gall, data=dA, family=binomial())
}

get_weights_fixedps <- function(d, ps_hat, method=c("ate","overlap")){
  method <- match.arg(method)
  if(method=="ate"){
    w <- ifelse(d$A==1, 1/ps_hat, 1/(1-ps_hat))
  } else {
    # overlap weights
    w <- ifelse(d$A==1, 1-ps_hat, ps_hat)
  }
  if(do_truncate){
    w <- truncate_w(w, w_cap_q)
  }
  w
}

# ---- 4) Outcome model ----
fit_cox <- function(d, w){
  fit <- coxph(
    Surv(tstart, tstop, event) ~ A + cluster(id),
    data = d,
    weights = w
  )
  b  <- unname(coef(fit)[["A"]])
  se <- sqrt(vcov(fit)[["A","A"]])
  c(beta=b, se=se)
}

# ---- 5) Balance/QBS summary (cobalt) ----
balance_summary <- function(d, w, wm){
  estimand_label <- if (wm == "ate") "ATE" else "ATO"
  
  bt <- bal.tab(
    A ~ age + sex + bmi + egfr + util + gall,
    data     = d,
    weights  = w,
    method   = "weighting",
    estimand = estimand_label
  )
  
  smd <- abs(bt$Balance$Diff.Adj)
  smd <- smd[is.finite(smd)]
  
  pct_smd_le_0.1 <- mean(smd <= 0.1)
  
  w_max <- max(w, na.rm=TRUE)
  w_med <- median(w, na.rm=TRUE)
  
  data.table(
    pct_smd_le_0.1 = pct_smd_le_0.1,
    max_smd        = if(length(smd)) max(smd) else NA_real_,
    ess            = ess(w),
    ess_over_n     = ess(w) / nrow(d),
    w_max          = w_max,
    w_med          = w_med,
    w_max_over_med = w_max / w_med
  )
}

# ---- 6) One replicate for a given scenario (tau2nd, p2nd_delta) ----
run_one <- function(tau2nd, p2nd_delta){
  # p2nd_A_effect is negative: A has fewer 2nd dispenses => mismatch engine
  p2nd_A_effect <- -p2nd_delta
  
  dt <- gen_covariates(N)
  dt <- assign_treatment(dt)
  dt <- gen_second_dispense(dt, tau2nd, p2nd_base, p2nd_A_effect)
  dt <- gen_event_times(dt, beta_true, lambda_event, lambda_death, followup_days)
  
  dA <- make_spec_A(dt)
  dB <- make_spec_B(dt)
  
  ps_fit  <- fit_ps_A(dA)
  ps_hat_A <- predict(ps_fit, newdata=dA, type="response")
  ps_hat_B <- predict(ps_fit, newdata=dB, type="response")
  
  out_est <- list()
  out_bal <- list()
  
  for(spec in c("A_init","B_confirm")){
    if(spec=="A_init"){
      d      <- dA
      ps_hat <- ps_hat_A
    } else {
      d      <- dB
      ps_hat <- ps_hat_B
    }
    
    for(wm in c("ate","overlap")){
      w <- get_weights_fixedps(d, ps_hat, wm)
      
      est <- fit_cox(d, w)
      bal <- balance_summary(d, w, wm)
      
      key <- paste0(spec,"_",wm)
      
      out_est[[key]] <- c(
        est,
        ess        = ess(w),
        n          = nrow(d),
        treat_rate = mean(d$A)
      )
      
      out_bal[[key]] <- c(
        n                 = nrow(d),
        treat_rate        = mean(d$A),
        ess               = bal$ess,
        ess_over_n        = bal$ess_over_n,
        w_max             = bal$w_max,
        w_med             = bal$w_med,
        w_max_over_med    = bal$w_max_over_med,
        pct_smd_le_0.1    = bal$pct_smd_le_0.1,
        max_smd           = bal$max_smd
      )
    }
  }
  
  list(
    est = unlist(out_est),
    bal = unlist(out_bal)
  )
}

# ============================================================
# 7) Run simulation: scenario grid (parallel within each scenario)
# ============================================================

scen_grid <- data.table(
  tau2nd     = rep(tau2nd_grid, each = length(p2nd_delta_grid)),
  p2nd_delta = rep(p2nd_delta_grid, times = length(tau2nd_grid))
)
scen_grid[, scenario := sprintf("tau%02d_d%.2f", tau2nd, p2nd_delta)]

cat("Scenario grid:\n")
print(scen_grid)
log_line(sprintf("Scenario grid: %s\n",
                 paste(scen_grid$scenario, collapse=", ")))

all_perf <- list()
all_bal  <- list()

keys <- c("A_init_ate","A_init_overlap","B_confirm_ate","B_confirm_overlap")

extract_spec <- function(res, key){
  data.table(
    beta = res[[paste0(key,".beta")]],
    se   = res[[paste0(key,".se")]],
    ess  = res[[paste0(key,".ess")]],
    n    = res[[paste0(key,".n")]]
  )
}

summ_bal <- function(res_bal, key){
  dt <- data.table(
    n              = res_bal[[paste0(key,".n")]],
    treat_rate     = res_bal[[paste0(key,".treat_rate")]],
    ess            = res_bal[[paste0(key,".ess")]],
    ess_over_n     = res_bal[[paste0(key,".ess_over_n")]],
    w_max          = res_bal[[paste0(key,".w_max")]],
    w_med          = res_bal[[paste0(key,".w_med")]],
    w_max_over_med = res_bal[[paste0(key,".w_max_over_med")]],
    pct_smd_le_0.1 = res_bal[[paste0(key,".pct_smd_le_0.1")]],
    max_smd        = res_bal[[paste0(key,".max_smd")]]
  )
  
  data.table(
    spec                  = key,
    n_median              = median(dt$n, na.rm=TRUE),
    treat_rate_median     = median(dt$treat_rate, na.rm=TRUE),
    ess_median            = median(dt$ess, na.rm=TRUE),
    ess_over_n_median     = median(dt$ess_over_n, na.rm=TRUE),
    w_max_median          = median(dt$w_max, na.rm=TRUE),
    w_max_over_med_median = median(dt$w_max_over_med, na.rm=TRUE),
    pct_smd_le_0.1_median = median(dt$pct_smd_le_0.1, na.rm=TRUE),
    max_smd_median        = median(dt$max_smd, na.rm=TRUE)
  )
}

for (i in seq_len(nrow(scen_grid))) {
  tau_i  <- scen_grid$tau2nd[i]
  dlt_i  <- scen_grid$p2nd_delta[i]
  scen_i <- scen_grid$scenario[i]
  
  cat("\n--- Scenario", i, "/", nrow(scen_grid), ":", scen_i, "---\n")
  log_line(sprintf("[%s] Start scenario %s (tau2nd=%d, p2nd_delta=%.2f)\n",
                   format(Sys.time(), "%H:%M:%S"), scen_i, tau_i, dlt_i))
  
  t_start <- proc.time()[3]
  
  # parallel over replicates
  rr_list <- future_lapply(
    X = 1:B,
    FUN = function(b){
      t0 <- proc.time()[3]
      rr <- run_one(tau2nd = tau_i, p2nd_delta = dlt_i)
      t1 <- proc.time()[3]
      # light log every 25 reps to avoid huge file
      if (b %% 25 == 0) {
        cat(sprintf("[%s] %s rep %d/%d | %.2f sec\n",
                    format(Sys.time(), "%H:%M:%S"),
                    scen_i, b, B, (t1 - t0)),
            file = log_file, append = TRUE)
      }
      rr
    },
    future.seed = TRUE
  )
  
  est_list <- lapply(rr_list, `[[`, "est")
  bal_list <- lapply(rr_list, `[[`, "bal")
  
  res_est <- rbindlist(lapply(est_list, as.list), fill=TRUE)
  res_bal <- rbindlist(lapply(bal_list, as.list), fill=TRUE)
  
  # performance summary for this scenario
  perf_i <- rbindlist(lapply(keys, function(k){
    dt_k <- extract_spec(res_est, k)
    s  <- summarize_perf(dt_k$beta, dt_k$se, beta_true)
    data.table(
      scenario   = scen_i,
      tau2nd     = tau_i,
      p2nd_delta = dlt_i,
      spec       = k,
      bias       = s$bias,
      rmse       = s$rmse,
      cover      = s$cover,
      sign_rev   = s$sign_rev,
      fpr_ci     = s$fpr_ci,
      ess_median = median(dt_k$ess, na.rm=TRUE),
      n_median   = median(dt_k$n,   na.rm=TRUE)
    )
  }))
  
  # balance summary for this scenario
  bal_i <- rbindlist(lapply(keys, function(k){
    tmp <- summ_bal(res_bal, k)
    tmp[, `:=`(scenario=scen_i, tau2nd=tau_i, p2nd_delta=dlt_i)]
    tmp
  }), fill=TRUE)
  
  t_end <- proc.time()[3]
  log_line(sprintf("[%s] End scenario %s | %.2f sec\n",
                   format(Sys.time(), "%H:%M:%S"), scen_i, (t_end - t_start)))
  
  all_perf[[i]] <- perf_i
  all_bal[[i]]  <- bal_i
}

perf_all <- rbindlist(all_perf, fill=TRUE)
bal_all  <- rbindlist(all_bal,  fill=TRUE)

# ============================================================
# 8) Save outputs
# ============================================================

perf_file <- file.path(out_dir, sprintf("perf_ZERO_ext_B%d_N%d.csv", B, N))
bal_file <- file.path(out_dir, sprintf("balance_summary_ZERO_ext_B%d_N%d.csv", B, N))

fwrite(perf_all, perf_file)
fwrite(bal_all,  bal_file)

sink(file.path(out_dir, "sessionInfo_ZERO_ext.txt"))
print(sessionInfo())
sink()

cat("\nSaved files:\n -", perf_file, "\n -", bal_file,
    "\n - sessionInfo_ZERO_ext.txt\n -", log_file, "\n\n")
print(perf_all)
print(bal_all)