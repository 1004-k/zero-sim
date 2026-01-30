# ============================================================
# ZERO simulation — Non-null single-scenario (Supplementary)
#  - Purpose: Option (Time 0/estimand) non-null single-scenario
#  - True effect: HR_true (0.80)
#  - Analytic options: A_init/B_confirm × ATE/overlap (4 Options)
#  - Outputs:
#      1) raw_est_ZERO_nonnull_B{B}_N{N}.csv
#      2) perf_ZERO_nonnull_B{B}_N{N}.csv
# ============================================================

# ---- Packages ----
req <- c("data.table","survival","cobalt","WeightIt","future.apply","future")
to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) library(data.table)
library(survival)
library(cobalt)
library(WeightIt)   # workflow consistency, PS -> glm
library(future)
library(future.apply)

set.seed(2026)

# ---- Parallel ----
n_cores <- as.integer(Sys.getenv("N_CORES","3"))
future::plan(future::multisession, workers = n_cores)

# ---- Log (local only, not OneDrive) ----
log_dir <- Sys.getenv("LOG_DIR", file.path("output","logs"))
out_dir <- Sys.getenv("OUT_DIR", file.path("output","raw"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, "sim_ZERO_nonnull.log")
if (file.exists(log_file)) file.remove(log_file)

log_line <- function(msg){
  cat(msg, file = log_file, append = TRUE)
}

log_line(sprintf("[%s] ZERO non-null single-scenario start\n",
                 format(Sys.time(), "%H:%M:%S")))

# ============================================================
# 1) SETTINGS
# ============================================================
B <- as.integer(Sys.getenv("B","500"))
N <- as.integer(Sys.getenv("N","50000"))
HR_true  <- 0.80     # non-null true effect (HR)
beta_true <- log(HR_true)

lambda_event   <- 0.00008
lambda_death   <- 0.00002
followup_days  <- 365

# single-scenario: tau2nd, p2nd_delta
tau2nd        <- 30          # Confirmation window (days) for the second dispensing
p2nd_delta    <- 0.15        # Differential reduction in second-dispensing probability
p2nd_base     <- 0.75        # Baseline probability of receiving a second dispensing

# ============================================================
# 2) Utility functions
# ============================================================

ess <- function(w){
  (sum(w)^2) / sum(w^2)
}

summarize_perf <- function(beta_hat, se_hat, beta_true){
  z  <- qnorm(0.975)
  lo <- beta_hat - z * se_hat
  hi <- beta_hat + z * se_hat
  
  list(
    bias     = mean(beta_hat - beta_true, na.rm = TRUE),
    rmse     = sqrt(mean((beta_hat - beta_true)^2, na.rm = TRUE)),
    cover    = mean(lo <= beta_true & beta_true <= hi, na.rm = TRUE),
    sign_rev = mean(beta_hat * beta_true < 0, na.rm = TRUE),
    fpr_ci   = mean(!(lo <= beta_true & beta_true <= hi), na.rm = TRUE)
  )
}

# ============================================================
# 3) DGM: covariates, treatment, second dispense, events
# ============================================================

gen_covariates <- function(N){
  dt <- data.table(
    id   = 1:N,
    age  = pmin(pmax(rnorm(N, 62, 10), 18), 90),
    sex  = rbinom(N, 1, 0.48),
    bmi  = pmin(pmax(rnorm(N, 31, 6), 18), 60),
    egfr = pmin(pmax(rnorm(N, 75, 18), 10), 120),
    util = rgamma(N, shape = 2, rate = 0.5),
    gall = rbinom(N, 1, plogis(-2.2 + 0.02*(rnorm(N, 62, 10) - 60)))
  )
  dt[]
}

assign_treatment <- function(dt){
  lin <- with(dt,
              -1.2 +
                0.02*(age - 60) +
                0.08*(bmi - 30) -
                0.015*(egfr - 75) +
                0.25*util +
                0.8*gall)
  ps <- plogis(lin)
  dt[, A := rbinom(.N, 1, ps)]  # 1=GLP-1RA, 0=SGLT2i
  dt[]
}

gen_second_dispense <- function(dt, tau2nd, p2nd_base, p2nd_A_effect){
  lin2 <- with(dt,
               qlogis(p2nd_base) +
                 p2nd_A_effect * A +
                 0.15*util - 0.4*gall +
                 0.03*(bmi - 30))
  p2 <- plogis(lin2)
  dt[, R2 := rbinom(.N, 1, p2)]
  dt[, t2 := ifelse(R2 == 1, sample(7:tau2nd, .N, replace = TRUE),
                    NA_integer_)]
  dt[]
}

gen_event_times <- function(dt, beta_true, lambda_event, lambda_death, followup_days){
  xb <- with(dt,
             0.02*(age - 60) +
               0.10*(bmi - 30) -
               0.01*(egfr - 75) +
               0.25*gall +
               0.10*util)
  
  haz_event <- lambda_event * exp(xb + beta_true * dt$A)
  haz_death <- lambda_death * exp(0.03*(dt$age - 60) + 0.05*dt$util)
  
  t_event <- rexp(nrow(dt), rate = haz_event)
  t_death <- rexp(nrow(dt), rate = haz_death)
  
  t_obs <- pmin(t_event, t_death, followup_days)
  status <- ifelse(t_obs == t_event & t_obs < pmin(t_death, followup_days), 1,
                   ifelse(t_obs == t_death & t_obs < pmin(t_event, followup_days), 2, 0))
  
  dt[, `:=`(time = t_obs, status = status)]
  dt[]
}

# ============================================================
# 4) Time-zero specs & PS/weights & Cox
# ============================================================

make_spec_A <- function(dt){
  dA <- copy(dt)
  dA[, `:=`(
    tstart = 0,
    tstop  = time,
    event  = as.integer(status == 1)
  )]
  dA[]
}

make_spec_B <- function(dt){
  dB <- dt[R2 == 1 & !is.na(t2)]
  dB[, `:=`(
    tstart = t2,
    tstop  = time,
    event  = as.integer(status == 1 & time > t2)
  )]
  dB <- dB[tstop > tstart]
  dB[]
}

fit_ps_A <- function(dA){
  glm(A ~ age + sex + bmi + egfr + util + gall,
      data = dA, family = binomial())
}

get_weights <- function(d, ps_hat, estimand = c("ate","overlap")){
  estimand <- match.arg(estimand)
  if (estimand == "ate") {
    w <- ifelse(d$A == 1, 1/ps_hat, 1/(1 - ps_hat))
  } else {
    w <- ifelse(d$A == 1, 1 - ps_hat, ps_hat)
  }
  w
}

fit_cox <- function(d, w){
  fit <- coxph(
    Surv(tstart, tstop, event) ~ A + cluster(id),
    data    = d,
    weights = w
  )
  b  <- unname(coef(fit)[["A"]])
  se <- sqrt(vcov(fit)[["A","A"]])
  c(beta = b, se = se)
}

# ============================================================
# 5) One replicate for this single scenario
# ============================================================

run_one <- function(){

  p2nd_A_effect <- -p2nd_delta
  
  dt <- gen_covariates(N)
  dt <- assign_treatment(dt)
  dt <- gen_second_dispense(dt, tau2nd, p2nd_base, p2nd_A_effect)
  dt <- gen_event_times(dt, beta_true, lambda_event, lambda_death, followup_days)
  
  dA <- make_spec_A(dt)
  dB <- make_spec_B(dt)
  
  ps_fit   <- fit_ps_A(dA)
  ps_hat_A <- predict(ps_fit, newdata = dA, type = "response")
  ps_hat_B <- predict(ps_fit, newdata = dB, type = "response")
  
  rows <- list()
  
  # 4 options: A_init/B_confirm × ATE/overlap
  for (spec in c("A_init","B_confirm")) {
    d      <- if (spec == "A_init") dA else dB
    ps_hat <- if (spec == "A_init") ps_hat_A else ps_hat_B
    
    for (wm in c("ate","overlap")) {
      w   <- get_weights(d, ps_hat, wm)
      est <- fit_cox(d, w)
      
      rows[[paste0(spec,"_",wm)]] <- data.table(
        spec      = paste0(spec,"_",wm),
        beta      = est[["beta"]],
        se        = est[["se"]],
        ess       = ess(w),
        n         = nrow(d),
        treat_rate= mean(d$A),
        HR_true   = HR_true,
        tau2nd    = tau2nd,
        p2nd_delta= p2nd_delta
      )
    }
  }
  
  rbindlist(rows)
}

# ============================================================
# 6) Run B = 500 replicates
# ============================================================

log_line(sprintf("Scenario: HR_true=%.2f, tau2nd=%d, p2nd_delta=%.2f\n",
                 HR_true, tau2nd, p2nd_delta))

res_list <- future_lapply(
  X = 1:B,
  FUN = function(b){
    dt <- run_one()
    dt[, rep := b]
    dt
  },
  future.seed = TRUE
)

raw_all <- rbindlist(res_list, fill = TRUE)

# ============================================================
# 7) Summarize performance
# ============================================================

perf_list <- lapply(split(raw_all, raw_all$spec), function(dt_k){
  s <- summarize_perf(
    beta_hat  = dt_k$beta,
    se_hat    = dt_k$se,
    beta_true = beta_true
  )
  data.table(
    spec       = unique(dt_k$spec),
    HR_true    = HR_true,
    tau2nd     = tau2nd,
    p2nd_delta = p2nd_delta,
    bias       = s$bias,
    rmse       = s$rmse,
    cover      = s$cover,
    sign_rev   = s$sign_rev,
    fpr_ci     = s$fpr_ci,
    ess_median = median(dt_k$ess, na.rm = TRUE),
    n_median   = median(dt_k$n,   na.rm = TRUE)
  )
})

perf_all <- rbindlist(perf_list, fill = TRUE)

# ============================================================
# 8) Save outputs
# ============================================================

raw_file <- file.path(out_dir, sprintf("raw_est_ZERO_nonnull_B%d_N%d.csv", B, N))
perf_file <- file.path(out_dir, sprintf("perf_ZERO_nonnull_B%d_N%d.csv", B, N))

fwrite(raw_all,  raw_file)
fwrite(perf_all, perf_file)

log_line(sprintf("[%s] Saved: %s, %s\n",
                 format(Sys.time(), "%H:%M:%S"), raw_file, perf_file))

print(perf_all)
cat("\nSaved:\n -", raw_file, "\n -", perf_file, "\n")
