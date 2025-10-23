
###############################################################################
# CEA_Glucorix.R — IGES interview assignment (Markov model in R)
# Author: <Hosein Shabaninejad>
# Purpose: Build a 3‑state cohort Markov CEA (Well → Complication → Dead)
# Horizon: lifetime (stop when alive < 1% or 50y)
# Cycle length: 1 year | Discount: 3% | Half‑cycle correction applied
# Arms: Glucorix (INT) vs Standard of Care (SoC)
#
# INPUT FILES required in project subfolder ./Data/ :
#   1) glucorix_ipd.csv — columns needed:
#        id, sex ("M"/"F"), age (at baseline),
#        event_complication (0/1), event_complication_time (years),
#        death (0/1), death_time (years)
#      (IPD are for Glucorix. If SoC IPD exist, can extend accordingly.)
#   2) lifetable.xlsx — UK life table with sheet "UK" (or first sheet) including:
#        age (integer), sex ("M"/"F"), qx (annual probability of death)
#
# OUTPUTS (written to ./Outputs/):
#   - base_case_results.csv  (per arm + incremental)
#   - psa_results.csv        (one row per PSA draw)
#   - ce_scatter.png         (ΔQALY vs ΔCost)
#   - ceac.png               (Cost-Effectiveness Acceptability Curve)
#   - tornado.png            (Deterministic SA tornado)
#
# Major packages: data.table, flexsurv, heemod (diagram), BCEA, ggplot2, mvtnorm,
#                 readxl, tidyr, patchwork, scales, triangle
###############################################################################





#############################
# Step 1 — Project set-up   #
#############################

## ---- 1A_project_and_folders ----------------------------------------------

proj_dir <- normalizePath(file.path(getwd(), "Glucorix_CEA_IGES"), mustWork = FALSE)
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

# Standard subfolders
dirs <- file.path(proj_dir, c("Data", "Markov", "Outputs", "Scripts"))
invisible(lapply(dirs, function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE)))

# Convenience: path builder
pf <- function(...) file.path(proj_dir, ...)

message("Project directory set at: ", proj_dir)
message("Subfolders present: ", paste(basename(dirs), collapse = ", "))



## ---- 1B_packages ----------------------------------------------------------
needed <- c(
  "readxl", "heemod", "BCEA", "flexsurv", "data.table", "ggplot2",
  "mvtnorm", "tidyr", "patchwork", "scales", "triangle"
)

installed <- rownames(installed.packages())
to_install <- setdiff(needed, installed)
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

# Load (stop on failure for clear errors)
ok <- sapply(needed, require, character.only = TRUE)
if (!all(ok)) stop("Failed to load: ", paste(names(ok)[!ok], collapse = ", "))



## ---- 1C_seed_and_options --------------------------------------------------
# One global seed; keep a named list of seeds if you want to split by component later
SEED <- 20251019L
set.seed(SEED)

options(stringsAsFactors = FALSE, scipen = 999)

# Stash seed & paths so other scripts can source this file and inherit them
globals <- list(
  proj_dir = proj_dir,
  SEED = SEED,
  data_dir = pf("Data"),
  markov_dir = pf("Markov"),
  outputs_dir = pf("Outputs"),
  scripts_dir = pf("Scripts")
)
saveRDS(globals, pf("Outputs", "globals.rds"))
message("Seed set to: ", SEED)




## ---- 1D_import_data_ALT (use this instead of the earlier 1D) -------------
ipd_path <- pf("Data", "glucorix_ipd.csv")
lt_path  <- pf("Data", "lifetable.xlsx")

if (!file.exists(ipd_path)) stop("Missing file: ", ipd_path)
if (!file.exists(lt_path))  stop("Missing file: ", lt_path)

# ---- IPD ----
ipd <- data.table::fread(ipd_path)
names(ipd) <- tolower(names(ipd))
req_ipd <- c("age", "sex", "event_complication", "event_complication_time", "death", "death_time")
missing_ipd <- setdiff(req_ipd, names(ipd))
if (length(missing_ipd)) stop("IPD missing required columns: ", paste(missing_ipd, collapse = ", "))
ipd[, sex := toupper(trimws(sex))]
ipd[!(sex %in% c("M","F")), sex := NA_character_]

# ---- Life table ----
lt_raw <- readxl::read_excel(lt_path)
lt_raw <- data.table::as.data.table(lt_raw)
names(lt_raw) <- tolower(names(lt_raw))

# Accept several common wide formats:
#   age | males | females
#   age | male  | female
wide_candidates <- c("males","male","men","m","females","female","women","f")
has_age <- "age" %in% names(lt_raw)
has_any_wide <- any(wide_candidates %in% names(lt_raw))
has_long <- all(c("age","sex") %in% names(lt_raw))

if (!has_age) stop("Life table must contain an 'age' column (integer years).")

if (has_long) {
  # Already long: expect columns: age, sex, qx or mx
  lt <- data.table::copy(lt_raw)
  if (!("qx" %in% names(lt)) && !("mx" %in% names(lt))) {
    stop("Life table (long) must have either 'qx' (annual death probability) or 'mx' (mortality rate).")
  }
} else if (has_any_wide) {
  # Wide -> Long. Try to find male & female columns.
  male_col <- intersect(names(lt_raw), c("males","male","men","m"))
  fem_col  <- intersect(names(lt_raw), c("females","female","women","f"))
  if (length(male_col) == 0 || length(fem_col) == 0) {
    stop("Could not detect both male and female columns in the life table.")
  }
  male_col <- male_col[1]; fem_col <- fem_col[1]
  
  # Rename to standard then pivot longer
  setnames(lt_raw, c("age", male_col, fem_col), c("age", "males", "females"))
  lt <- tidyr::pivot_longer(
    lt_raw,
    cols = c("males","females"),
    names_to = "sex",
    values_to = "qx"
  )
  lt <- data.table::as.data.table(lt)
  lt[, sex := ifelse(sex == "males", "M", "F")]
} else {
  stop("Life table must be either long (age/sex/ qx|mx) or wide with male/female columns.")
}

# Normalize columns
lt[, age := as.integer(round(age))]
lt[, sex := toupper(trimws(as.character(sex)))]

# If only mx provided, convert to qx; if qx provided, keep it.
if (!("qx" %in% names(lt)) && "mx" %in% names(lt)) {
  lt[, qx := 1 - exp(-mx)]
}
if (!("qx" %in% names(lt))) stop("Life table must contain or allow derivation of 'qx'.")

# Sanity checks
if (any(is.na(lt$qx)) || any(lt$qx < 0) || any(lt$qx > 1)) {
  stop("Life table 'qx' must be in [0,1] with no missing values.")
}
if (!all(lt$sex %in% c("M","F"))) {
  stop("Life table 'sex' must be 'M' or 'F' after processing.")
}

# ---- Quick summaries for later use ----
ipd_age_mean <- mean(ipd$age, na.rm = TRUE)
ipd_age_sd   <- sd(ipd$age,   na.rm = TRUE)
sex_prop <- ipd[, .N, by = sex][, prop := N / sum(N)][order(sex)]

message(sprintf("IPD age: mean=%.2f, sd=%.2f", ipd_age_mean, ipd_age_sd))
print(sex_prop)

# Save cleaned copies
data.table::fwrite(ipd, pf("Outputs", "ipd_clean_preview.csv"))
data.table::fwrite(lt[, .(age, sex, qx)], pf("Outputs", "lifetable_clean_preview.csv"))

saveRDS(list(ipd = ipd, lifetable = lt[, .(age, sex, qx)]), pf("Outputs", "inputs_clean.rds"))
message("Life table reshaped to long format with columns: age, sex, qx")




#############################
# Step 2 — Model Diagram   #
#############################

## ---- 2A_load_env ----------------------------------------------------------
# Reuse globals from Step 1
globals <- readRDS(file.path("Glucorix_CEA_IGES", "Outputs", "globals.rds"))
attach(globals)  # gives proj_dir, outputs_dir, etc.

# We'll use heemod for the transition object; DiagrammeR/htmlwidgets just help save a nice diagram.
for (pkg in c("heemod", "DiagrammeR", "htmlwidgets")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}


## ---- 2B_define_transition --------------------------------------------------
# State names per IGES brief
state_names <- c("Well", "Complication", "Dead")

# Placeholders for probabilities (we'll replace with real, cycle-specific values later).
# For the diagram, we only need structure; rows will still sum to 1 symbolically.
# Well row: stays Well, or moves to Complication, or Death
# Complication row: remains in Complication or Death (no return to Well)
# Dead row: absorbing
tmat <- heemod::define_transition(
  state_names = state_names,
  # From Well
  "1 - p_WC - p_WD", "p_WC",              "p_WD",
  # From Complication
  "0",               "1 - p_CD",          "p_CD",
  # From Dead (absorbing)
  "0",               "0",                 "1"
)

# Quick look
tmat

## ---- 2C_selfloops_plot_and_save -------------------------------------------
# Uses DOT -> SVG -> PNG (robust export path)
for (pkg in c("DiagrammeR", "DiagrammeRsvg", "rsvg")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

#in case of no attached globals in this session:
if (!exists("outputs_dir")) {
  globals <- readRDS(file.path("Glucorix_CEA_IGES", "Outputs", "globals.rds"))
  outputs_dir <- globals$outputs_dir
}

dot <- "
digraph G {
  graph [rankdir = LR, layout = dot]
  node  [shape = circle, fontsize = 12, style = filled, fillcolor = \"#e8f0fe\"]
  edge  [arrowsize = 0.9]

  Well         [label = 'Well']
  Complication [label = 'Complication']
  Dead         [label = 'Dead', fillcolor = '#fdecea']

  # Transitions
  Well         -> Complication [label = 'p_WC']
  Well         -> Dead         [label = 'p_WD']
  Complication -> Dead         [label = 'p_CD']

  # Self-loops (patients can remain in state)
  Well         -> Well         
  Complication -> Complication 

  # Dead is absorbing: no outgoing edges
}
"

g <- DiagrammeR::grViz(dot)   # show in Viewer
g

# Save SVG and PNG
svg_txt <- DiagrammeRsvg::export_svg(g)
out_svg <- file.path(outputs_dir, "state_diagram.svg")
out_png <- file.path(outputs_dir, "state_diagram.png")
writeLines(svg_txt, con = out_svg)
rsvg::rsvg_png(charToRaw(svg_txt), file = out_png)
message("Saved diagram:\n  ", out_svg, "\n  ", out_png)




#############################
# Step 3 — Model Input   #
#############################

## ---- 3A_load --------------------------------------------------------------
globals <- readRDS(file.path("Glucorix_CEA_IGES", "Outputs", "globals.rds"))
attach(globals)  # proj_dir, outputs_dir, etc.

inp <- readRDS(file.path(proj_dir, "Outputs", "inputs_clean.rds"))
ipd <- inp$ipd
lifetable <- inp$lifetable  # columns: age, sex, qx

# Quick sanity message
message(sprintf("Loaded IPD: n=%d; Life table rows=%d", nrow(ipd), nrow(lifetable)))



## ---- 3B_core_settings -----------------------------------------------------
# Treatments
treatments <- c("SoC", "Glucorix")

# Time
cycle_length_yr <- 1          # years
disc_rate_annual <- 0.03      # costs & QALYs
apply_HCC <- TRUE             # half-cycle correction (we'll implement in Step 6)

# Horizon control:
# We'll allow up to 50 cycles, but later stop early when alive < 1%.
n_cycles_cap <- 50L

# Starting age: mean and SD from IPD (PSA will sample from Normal later)
start_age_mean <- as.integer(mean(ipd$age, na.rm = TRUE))
start_age_sd   <- sd(ipd$age,   na.rm = TRUE)

# Sex mix for background mortality weighting
sex_mix <- as.data.frame(ipd[ , .N, by = sex][, prop := N/sum(N)][order(sex)])
names(sex_mix) <- c("sex","n","prop")
sex_mix




## ---- 3C_costs_utils -------------------------------------------------------
# State costs (€/year)
cost_well          <-  500
cost_complication  <- 5000
drug_cost_glucorix <- 2000   # applies only while in Well and only for Glucorix

# Utilities
u_glucorix_well <- 0.1055932
u_soc_well      <- 0.03738858
du_complication <- 0.20      # subtract from the corresponding well utility

# Bounds
u_lower <- -0.594
u_upper <-  1.000

# Helper to bound utilities
bound_u <- function(x, lower = u_lower, upper = u_upper) pmax(lower, pmin(upper, x))

# Arm- & state-specific utilities (bounded)
u_soc_well_bounded      <- bound_u(u_soc_well)
u_glucorix_well_bounded <- bound_u(u_glucorix_well)
u_soc_comp_bounded      <- bound_u(u_soc_well_bounded      - du_complication)
u_glucorix_comp_bounded <- bound_u(u_glucorix_well_bounded - du_complication)

data.frame(
  arm = c("SoC","SoC","Glucorix","Glucorix"),
  state = c("Well","Complication","Well","Complication"),
  utility = c(u_soc_well_bounded, u_soc_comp_bounded,
              u_glucorix_well_bounded, u_glucorix_comp_bounded)
)



## ---- 3D_discount_helpers --------------------------------------------------
# Annual discount factor for a given cycle (0-indexed), mid-cycle if HCC
disc_factor <- function(cycle, rate = disc_rate_annual, hcc = apply_HCC, cycle_len = cycle_length_yr) {
  # cycle 0 spans time [0,1). Mid-cycle occurs at 0.5 years, then 1.5, 2.5, ...
  t <- (cycle + ifelse(hcc, 0.5, 1.0)) * cycle_len
  1 / ((1 + rate) ^ t)
}
# quick test
disc_factor(0:5)


## ---- 3E_param_bundle ------------------------------------------------------
params <- list(
  meta = list(
    created = Sys.time(),
    seed = SEED
  ),
  structure = list(
    states = c("Well","Complication","Dead"),
    treatments = treatments,
    cycle_length_yr = cycle_length_yr,
    horizon_cap_cycles = n_cycles_cap,
    hcc = apply_HCC
  ),
  population = list(
    start_age_mean = start_age_mean,
    start_age_sd   = start_age_sd,
    sex_mix        = sex_mix,
    lifetable      = lifetable   # keep for background mortality lookup
  ),
  costs = list(
    cost_well = cost_well,
    cost_complication = cost_complication,
    drug_cost_glucorix = drug_cost_glucorix
  ),
  utilities = list(
    lower = u_lower, upper = u_upper, du_complication = du_complication,
    u_soc = list(well = u_soc_well_bounded, complication = u_soc_comp_bounded),
    u_glx = list(well = u_glucorix_well_bounded, complication = u_glucorix_comp_bounded)
  ),
  discounting = list(
    annual_rate = disc_rate_annual,
    factor_fun = disc_factor
  ),
  psa_priors = list(
    # Examples we’ll use later:
    start_age = list(dist = "normal", mean = start_age_mean, sd = start_age_sd)
    # (Other param distributions added at Step 8 when TPs are known)
  )
)

# Save to disk for later steps
saveRDS(params, file.path(outputs_dir, "params_step3.rds"))
write.csv(
  data.frame(
    key = c("cost_well","cost_complication","drug_cost_glucorix",
            "u_soc_well","u_soc_comp","u_glx_well","u_glx_comp",
            "disc_rate","hcc"),
    value = c(cost_well, cost_complication, drug_cost_glucorix,
              u_soc_well_bounded, u_soc_comp_bounded,
              u_glucorix_well_bounded, u_glucorix_comp_bounded,
              disc_rate_annual, apply_HCC)
  ),
  file.path(outputs_dir, "basecase_params_preview.csv"),
  row.names = FALSE
)

message("Saved parameter bundle to: ", file.path(outputs_dir, "params_step3.rds"))





#############################
# Step 4 — TP  #
#############################

## ===== RESET BACKGROUND MORTALITY, GLUCORIX p_WC, AND DIAGNOSTICS =====
library(data.table); library(flexsurv)

# Load inputs again (guard against stale objects)
globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)
params <- readRDS(file.path(outputs_dir, "params_step3.rds"))
inp <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
ipd <- as.data.table(inp$ipd)
lifetable <- as.data.table(inp$lifetable)  # age, sex, qx

# ---------- Background mortality (robust) ----------
lt <- copy(lifetable)
lt[, h_bg := -log(pmax(1 - qx, .Machine$double.eps))]
ltM <- lt[sex=="M"][order(age)]
ltF <- lt[sex=="F"][order(age)]
if (nrow(ltM)==0 || nrow(ltF)==0) stop("Life table must contain both sexes (M and F).")

# sex mix (fallback to 50/50 if invalid)
sex_mix <- params$population$sex_mix
w_m <- ifelse(any(sex_mix$sex=="M"), sex_mix$prop[sex_mix$sex=="M"], NA_real_)
w_f <- ifelse(any(sex_mix$sex=="F"), sex_mix$prop[sex_mix$sex=="F"], NA_real_)
if (length(w_m)==0 || length(w_f)==0 || is.na(w_m) || is.na(w_f) || (w_m+w_f)<=0) {
  w_m <- 0.5; w_f <- 0.5
}

h_BG <- function(age_years) {
  a <- as.integer(floor(pmax(0, as.numeric(age_years))))
  aM <- pmin(a, max(ltM$age)); aF <- pmin(a, max(ltF$age))
  hm <- ltM$h_bg[match(aM, ltM$age)]; hm[is.na(hm)] <- ltM$h_bg[.N]
  hf <- ltF$h_bg[match(aF, ltF$age)]; hf[is.na(hf)] <- ltF$h_bg[.N]
  as.numeric(w_m*hm + w_f*hf)
}

# ---------- All-cause hazard fit (if not in env) ----------
## === All-cause hazard h_all(t) – robust extractor ===========================
h_all <- function(tt) {
  # Try vectorized tidy=TRUE first
  df_try <- try(as.data.frame(
    summary(fit_death_best, t = tt, type = "hazard", ci = FALSE, tidy = TRUE)
  ), silent = TRUE)
  
  if (!inherits(df_try, "try-error") && "est" %in% names(df_try)) {
    est <- as.numeric(df_try$est)
    est[is.na(est)] <- 0
    return(pmax(est, 0))
  }
  
  # Fallback: per-time safe pull (handles list/data.frame returns)
  v <- sapply(tt, function(ti) {
    out <- summary(fit_death_best, t = ti, type = "hazard", ci = FALSE)
    if (is.data.frame(out) && "est" %in% names(out)) {
      as.numeric(out$est[1])
    } else if (is.list(out) && length(out) >= 1 &&
               is.data.frame(out[[1]]) && "est" %in% names(out[[1]])) {
      as.numeric(out[[1]]$est[1])
    } else {
      NA_real_
    }
  })
  v[is.na(v)] <- 0
  pmax(v, 0)
}

## === HR lives in params and is read at run-time ============================
params$hazards$hr_comp <- 1   # base-case HR

h_excess <- function(t, age) pmax(0, h_all(t) - h_BG(age))
h_WD     <- function(t, age) h_BG(age) + h_excess(t, age)
h_CD     <- function(t, age) params$hazards$hr_comp * h_WD(t, age)

## === Scalar-safe wrappers ===================================================
get_scalar <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return(default)
  x[1]
}
h_WD_safe <- function(t, age) get_scalar(h_WD(t, age), 0)
h_CD_safe <- function(t, age) get_scalar(h_CD(t, age), 0)





## ===== Refit Glucorix complication model + define p_WC_glx (cumhaz-based) =====
library(data.table); library(flexsurv)

# 1) Reload cleaned inputs to guarantee 'ipd' is present
globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)
inp <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
ipd <- as.data.table(inp$ipd)

# 2) Subset to Glucorix if an arm column exists; otherwise use all rows
if ("arm" %in% names(ipd)) {
  ipd_glx <- ipd[grepl("^glucorix$|^glx$|^g$", arm, ignore.case = TRUE)]
  if (nrow(ipd_glx) == 0) {
    message("No 'Glucorix' rows found in 'arm'; using all IPD for the complication fit.")
    ipd_glx <- copy(ipd)
  }
} else {
  ipd_glx <- copy(ipd)
}

# 3) Fit several parametric forms and pick best by AIC
dists <- c("weibull","gompertz","exp","lnorm","llogis")
fits_comp <- lapply(dists, function(d)
  flexsurvreg(Surv(event_complication_time, event_complication) ~ 1,
              data = ipd_glx, dist = d)
)
names(fits_comp) <- dists
aics <- sapply(fits_comp, AIC)
best_dist_comp <- names(which.min(aics))
fit_comp_best  <- fits_comp[[best_dist_comp]]
message("Refit Glucorix complications: best dist = ", best_dist_comp, 
        " (AIC = ", round(min(aics),1), ")")



## ---- Robust rebuild of p_wc_glx_vec (always yields a vector) --------------
library(flexsurv); library(data.table)

# Ensure the fitted object exists
if (!exists("fit_comp_best")) {
  inp <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
  ipd <- as.data.table(inp$ipd)
  dists <- c("weibull","gompertz","exp","lnorm","llogis")
  fits_comp <- lapply(dists, function(d)
    flexsurvreg(Surv(event_complication_time, event_complication) ~ 1,
                data = ipd, dist = d)
  )
  names(fits_comp) <- dists
  fit_comp_best <- fits_comp[[ which.min(sapply(fits_comp, AIC)) ]]
}

# integer times 0..N+1
n_precompute <- params$structure$horizon_cap_cycles
if (is.null(n_precompute) || !is.finite(n_precompute)) n_precompute <- 50L
t_grid <- 0:(n_precompute + 1L)

# 1) Try tidy=TRUE (usually returns one row per time)
S_grid <- try({
  df <- as.data.frame(summary(fit_comp_best,
                              t = t_grid,
                              type = "survival",
                              ci = FALSE,
                              tidy = TRUE))
  as.numeric(df$est)
}, silent = TRUE)

# 2) If that didn't work or length mismatch, fall back to per-time sapply
if (inherits(S_grid, "try-error") || length(S_grid) != length(t_grid)) {
  S_grid <- sapply(t_grid, function(tt) {
    out <- summary(fit_comp_best, t = tt, type = "survival", ci = FALSE)
    # out can be df or list(df); grab first 'est'
    if (is.data.frame(out) && "est" %in% names(out)) {
      as.numeric(out$est[1])
    } else if (is.list(out) && length(out) >= 1 && is.data.frame(out[[1]]) && "est" %in% names(out[[1]])) {
      as.numeric(out[[1]]$est[1])
    } else {
      NA_real_
    }
  })
}

# Guardrails
S_grid[is.na(S_grid)] <- .Machine$double.eps
S_grid <- pmin(pmax(S_grid, .Machine$double.eps), 1)

# Compute annual p: 1 - S(t+1)/S(t)  for t = 0..N
if (length(S_grid) < 2) stop("Still cannot get survival grid >1 point.")
p_wc_glx_vec <- 1 - (S_grid[-1] / S_grid[-length(S_grid)])
p_wc_glx_vec <- as.numeric(pmin(pmax(p_wc_glx_vec, 0), 1))

cat("Glucorix p_WC (t=0..10): ",
    paste(round(p_wc_glx_vec[1:min(11, length(p_wc_glx_vec))], 6), collapse = ", "),
    "\n")




# SoC annual p (constant) as a scalar function
p_WC_soc <- function(t) {
  # If already computed r_pool earlier, re-use it; otherwise recompute:
  soc_studies <- data.table::data.table(
    study=c("pubA_rate","pubA_ci3y","pubB_ci3y"),
    n    =c(120,120,150),
    type =c("rate","ci3y","ci3y"),
    value=c(0.08,0.25,0.20)
  )
  soc_studies[, r := ifelse(type=="rate", value, -log(1 - value)/3)]
  r_pool <- exp(soc_studies[, sum(n*log(r))/sum(n)])
  as.numeric(1 - exp(-r_pool))
}

# Glucorix precomputed p_WC vector (make sure it exists & numeric)
stopifnot(exists("p_wc_glx_vec"))
stopifnot(is.numeric(p_wc_glx_vec), length(p_wc_glx_vec) >= 1)

# Death hazard helpers should be present from Step 4:
stopifnot(exists("h_WD"), exists("h_CD"))



## ===== Rebuild death hazards (BG + excess) and safe wrappers ===============
library(data.table); library(flexsurv)

# Load inputs (life table + IPD) and params
globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)
inp    <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
params <- readRDS(file.path(outputs_dir, "params_step3.rds"))

ipd <- as.data.table(inp$ipd)
lt  <- as.data.table(inp$lifetable)   # expected columns: age, sex, qx

# --- Background mortality hazard from life table (sex-weighted) -------------
if (!all(c("age","sex","qx") %in% names(lt))) stop("lifetable must have age/sex/qx.")
lt[, h_bg := -log(pmax(1 - qx, .Machine$double.eps))]
ltM <- lt[sex=="M"][order(age)]
ltF <- lt[sex=="F"][order(age)]

# sex mix
sx <- params$population$sex_mix
w_m <- ifelse(any(sx$sex=="M"), sx$prop[sx$sex=="M"], 0.5)
w_f <- ifelse(any(sx$sex=="F"), sx$prop[sx$sex=="F"], 0.5)
if (!is.finite(w_m) || !is.finite(w_f) || w_m + w_f == 0) { w_m <- 0.5; w_f <- 0.5 }

h_BG <- function(age_years) {
  a  <- as.integer(floor(pmax(0, as.numeric(age_years))))
  aM <- pmin(a, max(ltM$age)); aF <- pmin(a, max(ltF$age))
  hm <- ltM$h_bg[match(aM, ltM$age)]; if (anyNA(hm)) hm[is.na(hm)] <- ltM$h_bg[.N]
  hf <- ltF$h_bg[match(aF, ltF$age)]; if (anyNA(hf)) hf[is.na(hf)] <- ltF$h_bg[.N]
  as.numeric(w_m*hm + w_f*hf)
}

# --- All-cause hazard from IPD (parametric) ---------------------------------
if (!exists("fit_death_best")) {
  dists_d <- c("weibull","gompertz","exp","lnorm","llogis")
  fits_d  <- lapply(dists_d, function(d) flexsurvreg(Surv(death_time, death) ~ 1, data = ipd, dist = d))
  names(fits_d) <- dists_d
  fit_death_best <- fits_d[[ which.min(sapply(fits_d, AIC)) ]]
}

h_all <- function(t) {
 est <- as.numeric(summary(fit_death_best, t = t, type = "hazard")[["est"]])
 if (length(est)==0 || is.na(est)) est <- 0
 pmax(est, 0)
}

# --- Excess & state-specific hazards ----------------------------------------


## Base-case HR lives in params and is read dynamically each call
params$hazards$hr_comp <- 1   # set base case here

h_excess <- function(t, age) pmax(0, h_all(t) - h_BG(age))
h_WD     <- function(t, age) h_BG(age) + h_excess(t, age)
h_CD     <- function(t, age) {
 hr <- params$hazards$hr_comp  # <- read current value each call
  hr * h_WD(t, age)
}

# safe wrappers (used in TPM builder)
get_scalar <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return(default)
  x[1]
}
h_WD_safe <- function(t, age) get_scalar(h_WD(t, age), 0)
h_CD_safe <- function(t, age) get_scalar(h_CD(t, age), 0)




# --- Scalar-safe wrappers (always return a single finite number) ------------
get_scalar <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x)==0 || is.na(x) || !is.finite(x)) return(default)
  x[1]
}
h_WD_safe <- function(t, age) get_scalar(h_WD(t, age), 0)
h_CD_safe <- function(t, age) get_scalar(h_CD(t, age), 0)

# --- Quick diagnostics at start age for first few cycles ---------------
start_age <- params$population$start_age_mean
cat("BG hazards ages (start..start+4): ",
    paste(round(sapply(0:4, function(k) h_BG(start_age+k)), 6), collapse=", "), "\n")
cat("h_all(t=0..4): ",
    paste(round(sapply(0:4, h_all), 6), collapse=", "), "\n")
cat("h_WD_safe(t=0..4): ",
    paste(round(sapply(0:4, function(k) h_WD_safe(k, start_age+k)), 6), collapse=", "), "\n")
cat("h_CD_safe(t=0..4): ",
    paste(round(sapply(0:4, function(k) h_CD_safe(k, start_age+k)), 6), collapse=", "), "\n")




############################################################build tpm list#######################

## ===== Safety: scalar getters for hazards (never length-0) ==============
get_scalar <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return(default)
  x[1]
}
h_WD_safe <- function(t, age) get_scalar(h_WD(t, age), 0)  # uses  Step-4 h_WD()
h_CD_safe <- function(t, age) get_scalar(h_CD(t, age), 0)  # uses  Step-4 h_CD()

## ===== Ensure Glucorix p_WC vector exists & numeric =====================
# If already created p_wc_glx_vec earlier, this will just check it.
if (!exists("p_wc_glx_vec") || !is.numeric(p_wc_glx_vec) || length(p_wc_glx_vec) < 2) {
  # Rebuild from the fitted model (works whether summary returns df or list)
  inp <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
  ipd <- data.table::as.data.table(inp$ipd)
  if (!exists("fit_comp_best")) {
    dists <- c("weibull","gompertz","exp","lnorm","llogis")
    fits_comp <- lapply(dists, function(d)
      flexsurv::flexsurvreg(Surv(event_complication_time, event_complication) ~ 1,
                            data = ipd, dist = d))
    names(fits_comp) <- dists
    fit_comp_best <- fits_comp[[ which.min(sapply(fits_comp, AIC)) ]]
  }
  n_pre <- params$structure$horizon_cap_cycles; if (is.null(n_pre)) n_pre <- 50L
  t_grid <- 0:(n_pre + 1L)
  
  sum_obj <- summary(fit_comp_best, t = t_grid, type = "survival", ci = FALSE)
  if (is.data.frame(sum_obj) && "est" %in% names(sum_obj)) {
    S_grid <- as.numeric(sum_obj$est)
  } else if (is.list(sum_obj)) {
    S_grid <- vapply(sum_obj, function(df)
      if (is.data.frame(df) && "est" %in% names(df)) as.numeric(df$est[1]) else NA_real_,
      numeric(1))
  } else {
    # Per-time fallback
    S_grid <- sapply(t_grid, function(tt) {
      out <- summary(fit_comp_best, t = tt, type = "survival", ci = FALSE)
      if (is.data.frame(out) && "est" %in% names(out)) as.numeric(out$est[1])
      else if (is.list(out) && is.data.frame(out[[1]]) && "est" %in% names(out[[1]])) as.numeric(out[[1]]$est[1])
      else NA_real_
    })
  }
  S_grid[is.na(S_grid)] <- .Machine$double.eps
  S_grid <- pmin(pmax(S_grid, .Machine$double.eps), 1)
  stopifnot(length(S_grid) >= 2)
  
  p_wc_glx_vec <- 1 - (S_grid[-1] / S_grid[-length(S_grid)])
  p_wc_glx_vec <- as.numeric(pmin(pmax(p_wc_glx_vec, 0), 1))
}
cat("Glucorix p_WC (t=0..10): ",
    paste(round(p_wc_glx_vec[1:min(11, length(p_wc_glx_vec))], 6), collapse=", "),
    "\n")

## ===== A per-cycle TPM builder (always returns numeric 3×3) =============
tpm_for_cycle <- function(arm = c("SoC","Glucorix"), k, start_age) {
  arm <- match.arg(arm)
  age_k <- as.numeric(start_age + k)
  t_k   <- as.numeric(k)
  
  # Well -> Complication probability (per-year)
  if (arm == "Glucorix") {
    idx <- k + 1L
    if (idx > length(p_wc_glx_vec)) idx <- length(p_wc_glx_vec)
    p_wc <- get_scalar(p_wc_glx_vec[idx], 0)
  } else {
    # SoC constant pooled p
    soc_studies <- data.table::data.table(
      study=c("pubA_rate","pubA_ci3y","pubB_ci3y"),
      n    =c(120,120,150),
      type =c("rate","ci3y","ci3y"),
      value=c(0.08,0.25,0.20)
    )
    soc_studies[, r := ifelse(type=="rate", value, -log(1 - value)/3)]
    r_pool <- exp(soc_studies[, sum(n*log(r))/sum(n)])
    p_wc <- 1 - exp(-r_pool)
  }
  p_wc <- pmin(pmax(p_wc, 0), 1)
  h_wc <- -log(pmax(1 - p_wc, .Machine$double.eps))
  
  # Death hazards (safe scalars)
  h_wd <- h_WD_safe(t_k, age_k)
  h_cd <- h_CD_safe(t_k, age_k)
  
  # Well row via competing risks
  H <- h_wc + h_wd
  if (!is.finite(H) || length(H) == 0) H <- 0
  if (H < .Machine$double.eps) {
    WW <- 1; WC <- 0; WD <- 0
  } else {
    S  <- exp(-H)
    pT <- 1 - S
    WC <- (h_wc / H) * pT
    WD <- (h_wd / H) * pT
    WW <- 1 - WC - WD
  }
  
  # Complication row (to Death only)
  if (!is.finite(h_cd) || h_cd < .Machine$double.eps) {
    CC <- 1; CD <- 0
  } else {
    CC <- exp(-h_cd); CD <- 1 - CC
  }
  
  # Assemble numeric matrix
  tpm <- matrix(c(WW, WC, WD,
                  0,  CC, CD,
                  0,  0,  1),
                nrow = 3, byrow = TRUE)
  storage.mode(tpm) <- "double"
  dimnames(tpm) <- list(from = c("Well","Complication","Dead"),
                        to   = c("Well","Complication","Dead"))
  tpm
}

## ===== Build list of TPMs with lapply (guaranteed matrices) =============
build_tpm_list2 <- function(arm = c("SoC","Glucorix"), start_age, n_cycles) {
  arm <- match.arg(arm)
  tpm_list <- lapply(0:(n_cycles - 1L), function(k) tpm_for_cycle(arm, k, start_age))
  # final sanity coercion
  for (i in seq_along(tpm_list)) {
    if (!is.matrix(tpm_list[[i]])) tpm_list[[i]] <- as.matrix(tpm_list[[i]])
    storage.mode(tpm_list[[i]]) <- "double"
  }
  tpm_list
}

## ===== Preview first 6 TPMs =============================================
start_age <- readRDS(file.path(outputs_dir,"params_step3.rds"))$population$start_age_mean
tpms_glx  <- build_tpm_list2("Glucorix", start_age = start_age, n_cycles = 12L)

for (i in 1:6) {
  cat("\nCycle", i-1, "TPM (Glucorix):\n")
  print(round(tpms_glx[[i]], 7))
}


## === Safe Markov runner that always uses build_tpm_list2 ====================
run_markov_arm <- function(arm = c("SoC","Glucorix"),
                           start_age,
                           n_cycles_cap = params$structure$horizon_cap_cycles) {
  arm <- match.arg(arm)
  
  # Build TPMs with the safe builder
  tpms <- build_tpm_list2(arm, start_age, n_cycles_cap)
  
  # Coerce each TPM to numeric 3x3
  for (k in seq_along(tpms)) {
    if (!is.matrix(tpms[[k]])) tpms[[k]] <- as.matrix(tpms[[k]])
    storage.mode(tpms[[k]]) <- "double"
    if (any(dim(tpms[[k]]) != c(3,3))) stop("TPM at cycle ", k-1, " is not 3x3.")
  }
  
  states <- c("Well","Complication","Dead")
  trace <- matrix(NA_real_, nrow = n_cycles_cap + 1, ncol = length(states),
                  dimnames = list(0:n_cycles_cap, states))
  trace[1, ] <- c(1, 0, 0)
  
  for (k in 1:n_cycles_cap) {
    M <- tpms[[k]]
    storage.mode(M) <- "double"
    trace[k + 1, ] <- as.numeric(trace[k, , drop = TRUE] %*% M)
  }
  
  alive <- rowSums(trace[, c("Well","Complication")])
  stop_k <- which(alive < 0.01)[1]; if (is.na(stop_k)) stop_k <- n_cycles_cap
  eff_n <- stop_k
  
  mid_occ <- (trace[1:eff_n, , drop = FALSE] + trace[2:(eff_n + 1), , drop = FALSE]) / 2
  
  # Utilities / costs
  u_soc <- with(params$utilities$u_soc, c(Well = well, Complication = complication, Dead = 0))
  u_glx <- with(params$utilities$u_glx, c(Well = well, Complication = complication, Dead = 0))
  u_vec <- if (arm == "SoC") u_soc else u_glx
  
  c_state  <- c(Well = params$costs$cost_well,
                Complication = params$costs$cost_complication,
                Dead = 0)
  drug_glx <- params$costs$drug_cost_glucorix
  
  cycle_len   <- params$structure$cycle_length_yr
  apply_HCC   <- params$structure$hcc
  disc_rate   <- params$discounting$annual_rate
  disc_factor <- params$discounting$factor_fun
  
  qaly_cycle <- as.numeric(mid_occ %*% u_vec) * cycle_len
  cost_cycle <- as.numeric(mid_occ %*% c_state) * cycle_len
  if (arm == "Glucorix") cost_cycle <- cost_cycle + (mid_occ[, "Well"] * drug_glx * cycle_len)
  
  df <- sapply(0:(eff_n - 1), function(k) disc_factor(k, rate = disc_rate, hcc = apply_HCC, cycle_len = cycle_len))
  qaly_disc <- qaly_cycle * df
  cost_disc <- cost_cycle * df
  
  list(
    arm = arm,
    tpms = tpms[1:eff_n],
    trace = trace[1:(eff_n + 1), , drop = FALSE],
    per_cycle = data.frame(
      cycle = 0:(eff_n - 1),
      age_start = round(start_age + 0:(eff_n - 1), 2),
      Well_mid = mid_occ[, "Well"],
      Comp_mid = mid_occ[, "Complication"],
      Dead_mid = mid_occ[, "Dead"],
      qaly_cycle = qaly_cycle,
      cost_cycle = cost_cycle,
      disc_factor = df,
      qaly_disc = qaly_disc,
      cost_disc = cost_disc
    ),
    total_qaly = sum(qaly_disc),
    total_cost = sum(cost_disc),
    effective_cycles = eff_n
  )
}


## =====  Re-run per-cycle table for Glucorix ==============================
#  run_markov_arm() uses build_tpm_list(). Point it to build_tpm_list2:
run_markov_arm <- local({
  run_markov_arm_inner <- run_markov_arm  # keep original
  function(arm = c("SoC","Glucorix"), start_age, n_cycles_cap = params$structure$horizon_cap_cycles) {
    arm <- match.arg(arm)
    tpms <- build_tpm_list2(arm, start_age, n_cycles_cap)  # << use the new builder
    # reuse the rest of  original function: copy-paste from Step 5C but replace tpms <- build_tpm_list(...) with the line above
    # --- BEGIN copy of  Step 5C body (shortened) ---
    states <- c("Well","Complication","Dead")
    trace <- matrix(NA_real_, nrow = n_cycles_cap + 1, ncol = length(states),
                    dimnames = list(0:n_cycles_cap, states))
    trace[1, ] <- c(1,0,0)
    for (k in 1:n_cycles_cap) trace[k+1,] <- trace[k,] %*% tpms[[k]]
    
    alive <- rowSums(trace[, c("Well","Complication")])
    stop_k <- which(alive < 0.01)[1]; if (is.na(stop_k)) stop_k <- n_cycles_cap
    eff_n <- stop_k
    
    mid_occ <- (trace[1:eff_n, , drop=FALSE] + trace[2:(eff_n+1), , drop=FALSE]) / 2
    
    u_soc <- with(params$utilities$u_soc, c(Well = well, Complication = complication, Dead = 0))
    u_glx <- with(params$utilities$u_glx, c(Well = well, Complication = complication, Dead = 0))
    u_vec <- if (arm=="SoC") u_soc else u_glx
    
    c_state <- c(Well = params$costs$cost_well, Complication = params$costs$cost_complication, Dead = 0)
    drug_glx <- params$costs$drug_cost_glucorix
    
    cycle_len  <- params$structure$cycle_length_yr
    apply_HCC  <- params$structure$hcc
    disc_rate  <- params$discounting$annual_rate
    disc_factor <- params$discounting$factor_fun
    
    qaly_cycle <- as.numeric(mid_occ %*% u_vec) * cycle_len
    cost_cycle <- as.numeric(mid_occ %*% c_state) * cycle_len
    if (arm == "Glucorix") cost_cycle <- cost_cycle + (mid_occ[, "Well"] * drug_glx * cycle_len)
    
    df <- sapply(0:(eff_n - 1), function(k) disc_factor(k, rate = disc_rate, hcc = apply_HCC, cycle_len = cycle_len))
    qaly_disc <- qaly_cycle * df
    cost_disc <- cost_cycle * df
    
    list(
      arm = arm,
      tpms = tpms[1:eff_n],
      trace = trace[1:(eff_n + 1), , drop = FALSE],
      per_cycle = data.frame(
        cycle = 0:(eff_n - 1),
        age_start = round(start_age + 0:(eff_n - 1), 2),
        Well_mid = mid_occ[, "Well"],
        Comp_mid = mid_occ[, "Complication"],
        Dead_mid = mid_occ[, "Dead"],
        qaly_cycle = qaly_cycle,
        cost_cycle = cost_cycle,
        disc_factor = df,
        qaly_disc = qaly_disc,
        cost_disc = cost_disc
      ),
      total_qaly = sum(qaly_disc),
      total_cost = sum(cost_disc),
      effective_cycles = eff_n
    )
    # --- END copy ---
  }
})

res_glx <- run_markov_arm("Glucorix", start_age = start_age)
print(res_glx$per_cycle[1:10, c("cycle","age_start","Well_mid","Comp_mid","Dead_mid","qaly_cycle","cost_cycle")],
     row.names = FALSE)


#############################
# Step 5 — Apply cost & qaly  #
#############################


#########################step 5#####################

## ---- 5A_load --------------------------------------------------------------
globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)  # proj_dir, outputs_dir, etc.

params <- readRDS(file.path(outputs_dir, "params_step3.rds"))

# make sure 'hazards' exists on the freshly reloaded params
if (is.null(params$hazards)) params$hazards <- list()
params$hazards$hr_comp <- 1  # <-- set  base-case HR here

# make sure 'hazards' exists
if (is.null(params$hazards)) params$hazards <- list()
params$hazards$hr_comp <- 1   # <-- base-case HR

# Use the hardened builder defined earlier in the script
build_tpm_list <- build_tpm_list2


# DO NOT load or use tp_helpers_step4.rds
# tp_help  <- readRDS(file.path(outputs_dir, "tp_helpers_step4.rds"))  # remove
# HR_comp  <- tp_help$HR_comp                                           # remove


# Bring key settings into the environment
states        <- params$structure$states
treatments    <- params$structure$treatments
n_cap         <- params$structure$horizon_cap_cycles
cycle_len     <- params$structure$cycle_length_yr
apply_HCC     <- params$structure$hcc
disc_rate     <- params$discounting$annual_rate
disc_factor   <- params$discounting$factor_fun

# Rebuild any functions that were not saved (RDS can’t store closures reliably)
#source_missing <- function(sym, expr) if (!exists(sym, inherits = FALSE)) eval(parse(text = expr))
#source_missing("build_tpm_list", "build_tpm_list <- tp_help$build_tpm_list")
#HR_comp <- tp_help$HR_comp



## ---- 5B_utils_costs -------------------------------------------------------
# Utilities by arm/state (already bounded in Step 3)
u_soc <- with(params$utilities$u_soc, c(Well = well, Complication = complication, Dead = 0))
u_glx <- with(params$utilities$u_glx, c(Well = well, Complication = complication, Dead = 0))

# Annual state costs (€/yr). Drug cost applies only to Glucorix in Well.
c_state <- c(Well = params$costs$cost_well,
             Complication = params$costs$cost_complication,
             Dead = 0)
drug_glx <- params$costs$drug_cost_glucorix

## ---- 5C_Markov runner with HCC + discounting ------------------------------------------------------------
# Stop rule: end when alive < 1% or hit cap
stop_alive_threshold <- 0.01

run_markov_arm <- function(arm = c("SoC","Glucorix"),
                           start_age = params$population$start_age_mean,
                           n_cycles_cap = n_cap) {
  
  arm <- match.arg(arm)
  # transition matrices (time-inhomogeneous)
  tpms <- build_tpm_list(arm, start_age = start_age, n_cycles = n_cycles_cap)
  
  # Trace: rows = cycles (0..n), cols = states
  trace <- matrix(NA_real_, nrow = n_cycles_cap + 1, ncol = length(states),
                  dimnames = list(0:n_cycles_cap, states))
  trace[1, ] <- c(1, 0, 0)  # start all in Well
  
  # Evolve the cohort
  for (k in 1:n_cycles_cap) trace[k + 1, ] <- trace[k, ] %*% tpms[[k]]
  
  # Determine effective length by alive threshold
  alive <- rowSums(trace[, c("Well","Complication")])
  stop_k <- which(alive < stop_alive_threshold)[1]
  if (is.na(stop_k)) stop_k <- n_cycles_cap
  # we’ll use cycles 0..stop_k
  eff_n <- stop_k
  
  # Mid-cycle (HCC) occupancies: average of start & end-of-cycle
  # For cycle k (0-indexed), mid-state = (trace[k,] + trace[k+1,]) / 2
  mid_occ <- (trace[1:eff_n, , drop = FALSE] + trace[2:(eff_n + 1), , drop = FALSE]) / 2
  
  # Choose utilities vector for the arm
  u_vec <- if (arm == "SoC") u_soc else u_glx
  
  # Per-cycle undiscounted QALYs and state costs (HCC)
  qaly_cycle <- as.numeric(mid_occ %*% u_vec) * cycle_len
  
  cost_cycle <- as.numeric(mid_occ %*% c_state) * cycle_len
  # Add drug cost (Glucorix only, only while in Well)
  if (arm == "Glucorix") {
    cost_cycle <- cost_cycle + (mid_occ[, "Well"] * drug_glx * cycle_len)
  }
  
  # Discount factors (mid-cycle if HCC)
  df <- sapply(0:(eff_n - 1), function(k) disc_factor(k, rate = disc_rate, hcc = apply_HCC, cycle_len = cycle_len))
  
  qaly_disc <- qaly_cycle * df
  cost_disc <- cost_cycle * df
  
  per_cycle <- data.frame(
    cycle = 0:(eff_n - 1),
    age_start = round(start_age + 0:(eff_n - 1), 2),
    Well_mid = mid_occ[, "Well"],
    Comp_mid = mid_occ[, "Complication"],
    Dead_mid = mid_occ[, "Dead"],
    qaly_cycle = qaly_cycle,
    cost_cycle = cost_cycle,
    disc_factor = df,
    qaly_disc = qaly_disc,
    cost_disc = cost_disc
  )
  
  out <- list(
    arm = arm,
    tpms = tpms[1:eff_n],
    trace = trace[1:(eff_n + 1), , drop = FALSE],
    per_cycle = per_cycle,
    total_qaly = sum(qaly_disc),
    total_cost = sum(cost_disc),
    effective_cycles = eff_n
  )
  out
}



## ===== Hardened Markov runner (uses build_tpm_list2 and enforces numeric TPMs) ====
run_markov_arm2 <- function(arm = c("SoC","Glucorix"),
                            start_age,
                            n_cycles_cap = params$structure$horizon_cap_cycles) {
  arm <- match.arg(arm)
  
  # Build TPMs with the safe builder  already created
  tpms <- build_tpm_list2(arm, start_age, n_cycles_cap)
  
  # Coerce every TPM to numeric 3x3 matrix (belt & suspenders)
  for (k in seq_along(tpms)) {
    if (!is.matrix(tpms[[k]])) tpms[[k]] <- as.matrix(tpms[[k]])
    storage.mode(tpms[[k]]) <- "double"
    if (any(dim(tpms[[k]]) != c(3,3))) stop("TPM at cycle ", k-1, " is not 3x3.")
  }
  
  states <- c("Well","Complication","Dead")
  
  trace <- matrix(NA_real_, nrow = n_cycles_cap + 1, ncol = length(states),
                  dimnames = list(0:n_cycles_cap, states))
  trace[1, ] <- c(1, 0, 0)
  
  for (k in 1:n_cycles_cap) {
    # multiply with explicit matrix coercion to avoid type issues
    M <- tpms[[k]]
    storage.mode(M) <- "double"
    trace[k + 1, ] <- as.numeric(trace[k, , drop = TRUE] %*% M)
  }
  
  alive <- rowSums(trace[, c("Well","Complication")])
  stop_k <- which(alive < 0.01)[1]; if (is.na(stop_k)) stop_k <- n_cycles_cap
  eff_n <- stop_k
  
  mid_occ <- (trace[1:eff_n, , drop = FALSE] + trace[2:(eff_n + 1), , drop = FALSE]) / 2
  
  # Utilities and costs (from Step 5B)
  u_soc <- with(params$utilities$u_soc, c(Well = well, Complication = complication, Dead = 0))
  u_glx <- with(params$utilities$u_glx, c(Well = well, Complication = complication, Dead = 0))
  u_vec <- if (arm == "SoC") u_soc else u_glx
  
  c_state  <- c(Well = params$costs$cost_well,
                Complication = params$costs$cost_complication,
                Dead = 0)
  drug_glx <- params$costs$drug_cost_glucorix
  
  cycle_len   <- params$structure$cycle_length_yr
  apply_HCC   <- params$structure$hcc
  disc_rate   <- params$discounting$annual_rate
  disc_factor <- params$discounting$factor_fun
  
  qaly_cycle <- as.numeric(mid_occ %*% u_vec) * cycle_len
  cost_cycle <- as.numeric(mid_occ %*% c_state) * cycle_len
  if (arm == "Glucorix") cost_cycle <- cost_cycle + (mid_occ[, "Well"] * drug_glx * cycle_len)
  
  df <- sapply(0:(eff_n - 1), function(k) disc_factor(k, rate = disc_rate, hcc = apply_HCC, cycle_len = cycle_len))
  qaly_disc <- qaly_cycle * df
  cost_disc <- cost_cycle * df
  
  list(
    arm = arm,
    tpms = tpms[1:eff_n],
    trace = trace[1:(eff_n + 1), , drop = FALSE],
    per_cycle = data.frame(
      cycle = 0:(eff_n - 1),
      age_start = round(start_age + 0:(eff_n - 1), 2),
      Well_mid = mid_occ[, "Well"],
      Comp_mid = mid_occ[, "Complication"],
      Dead_mid = mid_occ[, "Dead"],
      qaly_cycle = qaly_cycle,
      cost_cycle = cost_cycle,
      disc_factor = df,
      qaly_disc = qaly_disc,
      cost_disc = cost_disc
    ),
    total_qaly = sum(qaly_disc),
    total_cost = sum(cost_disc),
    effective_cycles = eff_n
  )
}





## ---- 5D_run_and_totals (fixed print) --------------------------------------
start_age <- params$population$start_age_mean

res_soc <- run_markov_arm2("SoC",      start_age = start_age)
res_glx <- run_markov_arm2("Glucorix", start_age = start_age)

inc_cost <- unname(res_glx$total_cost - res_soc$total_cost)
inc_qaly <- unname(res_glx$total_qaly - res_soc$total_qaly)
icer     <- if (is.finite(inc_qaly) && inc_qaly > 0) inc_cost / inc_qaly else NA_real_

basecase <- data.frame(
  outcome = c("total_cost_SoC",
              "total_qaly_SoC",
              "total_cost_Glucorix",
              "total_qaly_Glucorix",
              "inc_cost",
              "inc_qaly",
              "ICER"),
  value = c(res_soc$total_cost,
            res_soc$total_qaly,
            res_glx$total_cost,
            res_glx$total_qaly,
            inc_cost,
            inc_qaly,
            icer),
  stringsAsFactors = FALSE
)

# Print nicely (round only the numeric column)
basecase_print <- basecase
basecase_print$value <- round(basecase_print$value, 6)
print(basecase_print, row.names = FALSE)



## ---- 5E_save_outputs (unchanged, with rounded preview) --------------------
write.csv(res_soc$per_cycle, file.path(outputs_dir, "per_cycle_SoC.csv"), row.names = FALSE)
write.csv(res_glx$per_cycle, file.path(outputs_dir, "per_cycle_Glucorix.csv"), row.names = FALSE)
write.csv(basecase,         file.path(outputs_dir, "basecase_step5.csv"), row.names = FALSE)

# Optional: a rounded -readable table
basecase_print <- basecase
basecase_print$value <- round(basecase_print$value, 2)
write.csv(basecase_print,  file.path(outputs_dir, "basecase_step5_rounded.csv"), row.names = FALSE)









#############################
# step 6: Diagnostics and validation   #
#############################


##6a) ---- Inspect first  cycles for SoC and Glucorix --------------------------
head_n <- 30

cat("\n===== SoC: First", head_n, "cycles =====\n")
print(
  res_soc$per_cycle[1:head_n, 
                    c("cycle", "age_start", "Well_mid", "Comp_mid", "Dead_mid",
                      "qaly_cycle", "cost_cycle", "qaly_disc", "cost_disc")],
  row.names = FALSE
)

cat("\n===== Glucorix: First", head_n, "cycles =====\n")
print(
  res_glx$per_cycle[1:head_n, 
                    c("cycle", "age_start", "Well_mid", "Comp_mid", "Dead_mid",
                      "qaly_cycle", "cost_cycle", "qaly_disc", "cost_disc")],
  row.names = FALSE
)

## 6B: cumulative over first 5 cycles -----------------------------
cum_soc <- res_soc$per_cycle
cum_glx <- res_glx$per_cycle

cum_soc$cum_cost <- cumsum(cum_soc$cost_disc)
cum_soc$cum_qaly <- cumsum(cum_soc$qaly_disc)

cum_glx$cum_cost <- cumsum(cum_glx$cost_disc)
cum_glx$cum_qaly <- cumsum(cum_glx$qaly_disc)

cat("\nCumulative discounted cost & QALY (first 5 cycles)\n")
df_compare <- data.frame(
  cycle = 0:(head_n - 1),
  cum_cost_SoC = round(cum_soc$cum_cost[1:head_n], 2),
  cum_qaly_SoC = round(cum_soc$cum_qaly[1:head_n], 4),
  cum_cost_Glucorix = round(cum_glx$cum_cost[1:head_n], 2),
  cum_qaly_Glucorix = round(cum_glx$cum_qaly[1:head_n], 4)
)
print(df_compare, row.names = FALSE)




## ---- 6C_state_occupancy_age ------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(dplyr) })

# Base-case runs from Step 5 (uses run_markov_arm2)
start_age <- params$population$start_age_mean
res_soc_bc <- run_markov_arm2("SoC", start_age)
res_glx_bc <- run_markov_arm2("Glucorix", start_age)

# State occupancy over time
soc_occ <- res_soc_bc$per_cycle %>% 
  select(cycle, Well_mid, Comp_mid, Dead_mid) %>% 
  pivot_longer(-cycle, names_to="state", values_to="occupancy") %>%
  mutate(arm="SoC")
glx_occ <- res_glx_bc$per_cycle %>% 
  select(cycle, Well_mid, Comp_mid, Dead_mid) %>% 
  pivot_longer(-cycle, names_to="state", values_to="occupancy") %>%
  mutate(arm="Glucorix")
occ_df <- bind_rows(soc_occ, glx_occ)

p_occ <- ggplot(occ_df, aes(x = cycle, y = occupancy, color = state)) +
  geom_line() +
  facet_wrap(~arm, ncol = 1) +
  labs(title = "State occupancy over time (mid-cycle)", x = "Cycle (years)", y = "Proportion") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "diagnostic_state_occupancy.png"), p_occ, width = 8, height = 6, dpi = 300)

# Age trajectory (just for info)
age_df <- data.frame(cycle = res_soc_bc$per_cycle$cycle,
                     age_start = res_soc_bc$per_cycle$age_start)
p_age <- ggplot(age_df, aes(x=cycle, y=age_start)) +
  geom_line() + labs(title="Age trajectory", x="Cycle", y="Age (years)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "diagnostic_age_trajectory.png"), p_age, width = 6, height = 4, dpi = 300)



## ---- 6D) parametric_vs_KM_complication -------------------------------------
suppressPackageStartupMessages({ library(survival); library(flexsurv); library(ggplot2) })
inp <- readRDS(file.path(outputs_dir, "inputs_clean.rds"))
ipd <- as.data.table(inp$ipd)

# KM for time to first complication (all rows or Glucorix subset if 'arm' exists)
if ("arm" %in% names(ipd)) {
  ipd_glx <- ipd[grepl("^glucorix$|^glx$|^g$", arm, ignore.case = TRUE)]
  if (nrow(ipd_glx)==0) ipd_glx <- copy(ipd)
} else ipd_glx <- copy(ipd)

fit_km <- survfit(Surv(event_complication_time, event_complication) ~ 1, data = ipd_glx)
t_km   <- fit_km$time
S_km   <- fit_km$surv

# Use  chosen parametric fit fit_comp_best
dists <- c("weibull","gompertz","exp","lnorm","llogis")
fits  <- lapply(dists, function(d)
  flexsurvreg(Surv(event_complication_time, event_complication) ~ 1, data = ipd_glx, dist = d))
fit_comp_best <- fits[[ which.min(sapply(fits, AIC)) ]]

t_grid <- seq(0, max(t_km), by = 0.5)
S_par  <- as.data.frame(summary(fit_comp_best, t = t_grid, type = "survival", tidy = TRUE))$est

df_plot <- data.frame(
  time = c(t_km, t_grid),
  S    = c(S_km, S_par),
  type = rep(c("KM","Parametric"), c(length(S_km), length(S_par)))
)
p_km <- ggplot(df_plot, aes(x=time, y=S, color=type)) +
  geom_step(data = subset(df_plot, type=="KM")) +
  geom_line(data = subset(df_plot, type=="Parametric")) +
  labs(title="No-complication survival: KM vs parametric", x="Time (years)", y="S(t)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "diagnostic_km_vs_parametric_complication.png"), p_km, width = 7, height = 5, dpi = 300)


## ---- 6E)_probabilities_in_bounds -------------------------------------------
# Check a handful of TPMs to ensure rows sum to 1 and probs in [0,1]
check_tpm <- function(arm = "Glucorix", n_check = 10) {
  tpms <- build_tpm_list2(arm, start_age = params$population$start_age_mean, n_cycles = n_check)
  errs <- 0
  for (k in seq_len(n_check)) {
    M <- tpms[[k]]
    if (any(M < -1e-10) || any(M > 1 + 1e-10)) errs <- errs + 1
    if (abs(rowSums(M)[1] - 1) > 1e-8 || abs(rowSums(M)[2] - 1) > 1e-8 || abs(rowSums(M)[3] - 1) > 1e-8) errs <- errs + 1
  }
  if (errs == 0) message("TPM checks passed for arm=", arm) else message("TPM checks found ", errs, " issues for arm=", arm)
}
check_tpm("SoC", 12); check_tpm("Glucorix", 12)





#############################
# Step 7- Sensitivity analysis  #
#############################


#####################################################################  SA ######################

## ---- 7A_helpers_overrides (fixed) -----------------------------------------
# Safe clamp
safel <- function(x, lower, upper) pmax(lower, pmin(upper, x))

par0 <- params
on.exit({ params <<- par0 }, add = FALSE)   # always restore

run_ce_with_overrides <- function(
    overrides = list(),
    lambda = 30000,
    hr_arm_specific = TRUE   # <- TRUE: HR applies only in Glucorix; FALSE: both arms
) {
  # keep a copy of params and restore at the end
  par0 <- params
  
  ## ----- apply simple overrides to params -----
  if (!is.null(overrides$cost_well))          params$costs$cost_well <- overrides$cost_well
  if (!is.null(overrides$cost_complication))  params$costs$cost_complication <- overrides$cost_complication
  if (!is.null(overrides$drug_cost_glucorix)) params$costs$drug_cost_glucorix <- overrides$drug_cost_glucorix
  
  if (!is.null(overrides$u_soc_well)) {
    u_soc_well_b <- safel(overrides$u_soc_well, params$utilities$lower, params$utilities$upper)
    params$utilities$u_soc$well <- u_soc_well_b
    params$utilities$u_soc$complication <- safel(u_soc_well_b - params$utilities$du_complication,
                                                 params$utilities$lower, params$utilities$upper)
  }
  if (!is.null(overrides$u_glx_well)) {
    u_glx_well_b <- safel(overrides$u_glx_well, params$utilities$lower, params$utilities$upper)
    params$utilities$u_glx$well <- u_glx_well_b
    params$utilities$u_glx$complication <- safel(u_glx_well_b - params$utilities$du_complication,
                                                 params$utilities$lower, params$utilities$upper)
  }
  if (!is.null(overrides$du_complication)) {
    params$utilities$du_complication <- overrides$du_complication
    params$utilities$u_soc$complication <- safel(params$utilities$u_soc$well - overrides$du_complication,
                                                 params$utilities$lower, params$utilities$upper)
    params$utilities$u_glx$complication <- safel(params$utilities$u_glx$well - overrides$du_complication,
                                                 params$utilities$lower, params$utilities$upper)
  }
  
  if (!is.null(overrides$discount)) params$discounting$annual_rate <- overrides$discount
  
  if (is.null(params$hazards)) params$hazards <- list()
  if (!is.null(overrides$hr_comp)) params$hazards$hr_comp <- overrides$hr_comp
  
  soc_rate_mult <- if (!is.null(overrides$soc_rate_mult)) overrides$soc_rate_mult else 1.0
  excess_mult   <- if (!is.null(overrides$excess_mult))   overrides$excess_mult   else 1.0
  
  # handle possibly NULL input robustly
  stop_drug_after <- overrides$stop_drug_after
  stop_drug_after <- suppressWarnings(as.numeric(stop_drug_after))
  if (length(stop_drug_after) == 0 || !is.finite(stop_drug_after)) stop_drug_after <- NA_real_
  
  # --- local safe helpers (do NOT rely on globals) ---
  safe_num <- function(x, fallback) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0 || is.na(x) || !is.finite(x)) return(fallback)
    x[1]
  }
  safe_h_wd <- function(t_k, age_k) {
    bg <- h_BG(age_k)
    wd <- safe_num(h_WD(t_k, age_k), NA_real_)
    if (is.na(wd) || wd <= 0) bg else wd            # fallback to background
  }
  safe_h_cd <- function(t_k, age_k, arm, hr_arm_specific, hr) {
    bg <- h_BG(age_k)
    wd <- safe_h_wd(t_k, age_k)
    cd <- if (hr_arm_specific) {
      if (arm == "Glucorix") hr * wd else wd
    } else {
      hr * wd
    }
    max(bg, cd)                                     # NEVER below background
  }
  
  
  ## ----- per-cycle TPM builder with multipliers applied --------------------
  tpm_for_cycle_sa <- function(arm = c("SoC","Glucorix"), k, start_age) {
    arm  <- match.arg(arm)
    age_k <- as.numeric(start_age + k)
    t_k   <- as.numeric(k)
    
    # Well -> Complication (as in base model) + optional multipliers
    if (arm == "Glucorix") {
      idx <- k + 1L; if (idx > length(p_wc_glx_vec)) idx <- length(p_wc_glx_vec)
      p_wc <- safel(p_wc_glx_vec[idx], 0, 1)
      # PSA multiplier (if set)
      mult_glx <- getOption("glx_wc_mult_psa", 1)
      h_wc <- -log(pmax(1 - p_wc, .Machine$double.eps)) * mult_glx
      p_wc <- 1 - exp(-h_wc)
    } else {
      # pooled SoC rate then scale by soc_rate_mult
      soc_studies <- data.table::data.table(
        study=c("pubA_rate","pubA_ci3y","pubB_ci3y"),
        n    =c(120,120,150),
        type =c("rate","ci3y","ci3y"),
        value=c(0.08,0.25,0.20)
      )
      soc_studies[, r := ifelse(type=="rate", value, -log(1 - value)/3)]
      r_pool <- exp(soc_studies[, sum(n*log(r))/sum(n)]) * soc_rate_mult
      p_wc <- 1 - exp(-r_pool)
    }
    p_wc <- safel(p_wc, 0, 1)
    h_wc <- -log(pmax(1 - p_wc, .Machine$double.eps))
    
    # Death hazards: background + (scaled) excess; CD via (arm-specific) HR
    # NEW (robust & aligned with base model)
    hr_c <- safe_num(params$hazards$hr_comp, 1)
    # scale excess by excess_mult, but keep ≥ background
    h_bg <- h_BG(age_k)
    h_all_k <- safe_num(h_all(t_k), NA_real_)
    if (is.na(h_all_k) || h_all_k < 0) h_all_k <- h_bg
    h_ex <- pmax(0, h_all_k - h_bg) * excess_mult
    h_wd <- max(h_bg, h_bg + h_ex)                      # floor at background
    h_cd <- safe_h_cd(t_k, age_k, arm, hr_arm_specific, hr_c)
    
    
    
    # Competing risks for Well row
    H <- h_wc + h_wd
    if (!is.finite(H) || H < .Machine$double.eps) {
      WW <- 1; WC <- 0; WD <- 0
    } else {
      S  <- exp(-H)
      pT <- 1 - S
      WC <- (h_wc / H) * pT
      WD <- (h_wd / H) * pT
      WW <- 1 - WC - WD
    }
    
    # Complication row
    if (!is.finite(h_cd) || h_cd < .Machine$double.eps) {
      CC <- 1; CD <- 0
    } else {
      CC <- exp(-h_cd); CD <- 1 - CC
    }
    
    matrix(c(WW, WC, WD,
             0,  CC, CD,
             0,  0,  1),
           nrow = 3, byrow = TRUE,
           dimnames = list(from = c("Well","Complication","Dead"),
                           to   = c("Well","Complication","Dead")))
  }
  
  build_tpm_list_sa <- function(arm = c("SoC","Glucorix"), start_age, n_cycles) {
    arm <- match.arg(arm)
    lapply(0:(n_cycles - 1L), function(k) tpm_for_cycle_sa(arm, k, start_age))
  }
  
  run_arm_sa <- function(arm, start_age = params$population$start_age_mean,
                         n_cycles = params$structure$horizon_cap_cycles) {
    tpms <- build_tpm_list_sa(arm, start_age, n_cycles)
    states <- c("Well","Complication","Dead")
    trace <- matrix(NA_real_, nrow = n_cycles + 1, ncol = 3,
                    dimnames = list(0:n_cycles, states))
    trace[1,] <- c(1,0,0)
    for (k in 1:n_cycles) trace[k+1,] <- as.numeric(trace[k,]) %*% tpms[[k]]
    
    # stop when alive < 1%
    alive <- rowSums(trace[,1:2]); stop_k <- which(alive < 0.01)[1]
    if (is.na(stop_k)) stop_k <- n_cycles
    eff_n <- stop_k
    mid_occ <- (trace[1:eff_n, , drop=FALSE] + trace[2:(eff_n+1), , drop=FALSE]) / 2
    
    cycle_len <- params$structure$cycle_length_yr
    df <- sapply(0:(eff_n - 1),
                 function(k) params$discounting$factor_fun(k,
                                                           rate = params$discounting$annual_rate,
                                                           hcc  = params$structure$hcc,
                                                           cycle_len = cycle_len))
    
    u_soc <- with(params$utilities$u_soc, c(Well=well, Complication=complication, Dead=0))
    u_glx <- with(params$utilities$u_glx, c(Well=well, Complication=complication, Dead=0))
    u_vec <- if (arm=="SoC") u_soc else u_glx
    
    c_state <- c(Well=params$costs$cost_well,
                 Complication=params$costs$cost_complication,
                 Dead=0)
    drug_glx <- params$costs$drug_cost_glucorix
    
    qaly_cycle <- as.numeric(mid_occ %*% u_vec) * cycle_len
    cost_cycle <- as.numeric(mid_occ %*% c_state) * cycle_len
    
    if (arm == "Glucorix") {
      k_idx <- 0:(eff_n - 1)
      if (is.finite(stop_drug_after)) {
        on_tx <- as.numeric(k_idx < stop_drug_after)
        cost_cycle <- cost_cycle + (mid_occ[,"Well"] * drug_glx * cycle_len * on_tx)
      } else {
        cost_cycle <- cost_cycle + (mid_occ[,"Well"] * drug_glx * cycle_len)
      }
    }
    
    list(
      total_qaly = sum(qaly_cycle * df),
      total_cost = sum(cost_cycle * df)
    )
  }
  
  # run both arms
  start_age <- params$population$start_age_mean
  res_soc <- run_arm_sa("SoC", start_age)
  res_glx <- run_arm_sa("Glucorix", start_age)
  
  inc_cost <- res_glx$total_cost - res_soc$total_cost
  inc_qaly <- res_glx$total_qaly - res_soc$total_qaly
  icer <- if (is.finite(inc_qaly) && inc_qaly > 0) inc_cost / inc_qaly else NA_real_
  inb  <- lambda * inc_qaly - inc_cost
  
  # restore params (important)
  params <<- par0
  
  c(total_cost_SoC = res_soc$total_cost,
    total_qaly_SoC = res_soc$total_qaly,
    total_cost_Glucorix = res_glx$total_cost,
    total_qaly_Glucorix = res_glx$total_qaly,
    inc_cost = inc_cost, inc_qaly = inc_qaly, ICER = icer, INB = inb)
}



## ---- 7B_DSA_tornado INB -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

if (!exists("outputs_dir")) {
  globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
  outputs_dir <- globals$outputs_dir
}

lambda_ref <- 30000  # use 30k by default for tornado (edit if needed)

# Base values (read from params)
base_vals <- list(
  cost_well          = params$costs$cost_well,
  cost_complication  = params$costs$cost_complication,
  drug_cost_glucorix = params$costs$drug_cost_glucorix,
  u_soc_well         = params$utilities$u_soc$well,
  u_glx_well         = params$utilities$u_glx$well,
  du_complication    = params$utilities$du_complication,
  hr_comp            = ifelse(is.null(params$hazards$hr_comp), 1.5, params$hazards$hr_comp),
  soc_rate_mult      = 1.0,              # relative change to pooled SoC complication rate
  #excess_mult        = 1.0,              # relative change to excess mortality
  #discount           = params$discounting$annual_rate,
  stop_drug_after    = NA_real_          # NA = never stop; scenario will set 5
)

# Build low/high levels
mk_range <- function(x) c(low = 0.8*x, high = 1.2*x)

ranges <- list(
  cost_well          = mk_range(base_vals$cost_well),
  cost_complication  = mk_range(base_vals$cost_complication),
  drug_cost_glucorix = mk_range(base_vals$drug_cost_glucorix),
  u_soc_well         = mk_range(base_vals$u_soc_well),
  u_glx_well         = mk_range(base_vals$u_glx_well),
  du_complication    = mk_range(base_vals$du_complication),
  hr_comp            = mk_range(base_vals$hr_comp),
  soc_rate_mult      = mk_range(base_vals$soc_rate_mult)
  #excess_mult        = mk_range(base_vals$excess_mult)
)
# Special DSA nodes:
dsa_special <- list(
  #discount = c(low = 0.00, high = 0.05),
  stop_drug_after = c(low = 5, high = NA)  # compare "stop at 5y" vs "never stop"
)

# Run DSA one-way
do_one <- function(name, low, high) {
  o_low  <- base_vals;  o_low[[name]]  <- low
  o_high <- base_vals;  o_high[[name]] <- high
  res_low  <- run_ce_with_overrides(o_low,  lambda = lambda_ref)
  res_high <- run_ce_with_overrides(o_high, lambda = lambda_ref)
  data.frame(
    param = name,
    INB_low  = res_low["INB"],
    INB_high = res_high["INB"],
    ICER_low = res_low["ICER"],
    ICER_high= res_high["ICER"],
    stringsAsFactors = FALSE
  )
}

dsa_results <- do.call(rbind, c(
  # regular ±20%
  lapply(names(ranges), function(k) do_one(k, ranges[[k]]["low"], ranges[[k]]["high"])),
  # specials
  list(
    #do_one("discount", dsa_special$discount["low"], dsa_special$discount["high"]),
    do_one("stop_drug_after", dsa_special$stop_drug_after["low"], dsa_special$stop_drug_after["high"])
  )
))
dsa_results$param <- factor(dsa_results$param)

# Tornado  on  range
dsa_results$INB_min <- pmin(dsa_results$INB_low, dsa_results$INB_high)
dsa_results$INB_max <- pmax(dsa_results$INB_low, dsa_results$INB_high)
dsa_results$width   <- dsa_results$INB_max - dsa_results$INB_min

dsa_plot_df <- dsa_results[order(dsa_results$width, decreasing = TRUE), ]
dsa_plot_df$param <- factor(dsa_plot_df$param, levels = dsa_plot_df$param)

p_tornado <- ggplot(dsa_plot_df, aes(y = param)) +
  geom_segment(aes(x = INB_min, xend = INB_max, yend = param), size = 6, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(labels = label_dollar(prefix = "€")) +
  labs(title = "Tornado diagram (Incremental Net Benefit, λ = €30,000/QALY)",
       x = "INB range (low–high)", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(outputs_dir, "tornado_INB.png"), p_tornado, width = 8, height = 6, dpi = 300)

# Save the DSA table
write.csv(dsa_results, file.path(outputs_dir, "dsa_results.csv"), row.names = FALSE)



# -Tornado on ICER: run base-case once to mark it on the plot
res_base   <- run_ce_with_overrides(base_vals, lambda = lambda_ref)
base_ICER  <- as.numeric(res_base["ICER"])

# --- Build ICER tornado data
dsa_results$ICER_min   <- pmin(dsa_results$ICER_low,  dsa_results$ICER_high)
dsa_results$ICER_max   <- pmax(dsa_results$ICER_low,  dsa_results$ICER_high)
dsa_results$width_ICER <- dsa_results$ICER_max - dsa_results$ICER_min

# Keep only finite ICERs for plotting (avoid Inf/NaN from tiny or negative ΔQALY)
dsa_plot_icer_df <- subset(dsa_results, is.finite(ICER_min) & is.finite(ICER_max))
dsa_plot_icer_df <- dsa_plot_icer_df[order(dsa_plot_icer_df$width_ICER, decreasing = TRUE), ]
dsa_plot_icer_df$param <- factor(dsa_plot_icer_df$param, levels = dsa_plot_icer_df$param)

# --- ICER tornado plot
p_tornado_ICER <- ggplot(dsa_plot_icer_df, aes(y = param)) +
  geom_segment(aes(x = ICER_min, xend = ICER_max, yend = param), size = 6, alpha = 0.8) +
  geom_vline(xintercept = base_ICER, linetype = 2) +
  scale_x_continuous(labels = label_dollar(prefix = "€")) +
  labs(title = "Tornado diagram (ICER)",
       subtitle = sprintf("Base-case ICER = €%s per QALY", scales::comma(round(base_ICER,0))),
       x = "ICER range (low–high; €/QALY)", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(outputs_dir, "tornado_ICER.png"), p_tornado_ICER, width = 8, height = 6, dpi = 300)


#############################
        # step 8 - PSA  #
#############################



## ---- 8A_PSA_setup (helpers) -----------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

# Load paths and params
globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)
params <- readRDS(file.path(outputs_dir, "params_step3.rds"))

# Safeguards
if (!exists("SEED")) SEED <- 20251019L
if (!exists("safel")) safel <- function(x, lower, upper) pmax(lower, pmin(upper, x))

# ---- Distributions ----
# Beta from mean & sd (bounded)
beta_from_mean_sd <- function(mu, sd) {
  mu <- safel(mu, 1e-4, 1 - 1e-4)
  v  <- max(sd, 1e-6)^2
  a  <- ((1 - mu) / v - 1 / mu) * mu^2
  b  <- a * (1/mu - 1)
  a <- pmax(a, 1e-3); b <- pmax(b, 1e-3)
  c(shape1 = a, shape2 = b)
}

# Gamma from mean & sd
gamma_from_mean_sd <- function(mu, sd) {
  sd <- max(sd, 1e-6)
  shape <- (mu/sd)^2
  scale <- (sd^2)/mu
  c(shape=shape, scale=scale)
}

# Lognormal from mean & CV
ln_from_mean_cv <- function(mu, cv = 0.2) {
  cv <- max(cv, 1e-6)
  sigma2 <- log(1 + cv^2)
  mu_log <- log(mu) - 0.5*sigma2
  c(meanlog = mu_log, sdlog = sqrt(sigma2))
}

# ---- SoC complication baseline (pooled) ----
get_soc_rate_pool <- function() {
  soc_studies <- data.table(
    study=c("pubA_rate","pubA_ci3y","pubB_ci3y"),
    n    =c(120,120,150),
    type =c("rate","ci3y","ci3y"),
    value=c(0.08,0.25,0.20)
  )
  soc_studies[, r := ifelse(type=="rate", value, -log(1 - value)/3)]
  as.numeric(exp(soc_studies[, sum(n*log(r))/sum(n)]))
}
base_soc_rate <- get_soc_rate_pool()
base_soc_prob <- 1 - exp(-base_soc_rate)

# ---- Glucorix complication hazards vector from  fitted survival ----
# p_wc_glx_vec must exist from Step 4; convert to hazards to allow multiplicative noise
stopifnot(exists("p_wc_glx_vec"))
base_h_wc_glx <- -log(pmax(1 - p_wc_glx_vec, .Machine$double.eps))



## ---- 8B_PSA_one_draw -------------------------------------------------------
# NOTE: relies on run_ce_with_overrides()  added in Step 8A (uses hazards math)
psa_one_draw <- function(thresholds = c(20000,30000,50000)) {
  # --- Costs (Gamma, assume 20% sd if no SEs) ---
  c_w  <- gamma_from_mean_sd(params$costs$cost_well,         0.2*params$costs$cost_well)
  c_c  <- gamma_from_mean_sd(params$costs$cost_complication, 0.2*params$costs$cost_complication)
  c_dr <- gamma_from_mean_sd(params$costs$drug_cost_glucorix,0.2*params$costs$drug_cost_glucorix)
  cost_well          <- rgamma(1, shape=c_w["shape"],  scale=c_w["scale"])
  cost_complication  <- rgamma(1, shape=c_c["shape"],  scale=c_c["scale"])
  drug_cost_glucorix <- rgamma(1, shape=c_dr["shape"], scale=c_dr["scale"])
  
  # --- Utilities (Beta; use sd = 20%*mean, with floor to avoid degeneracy) ---
  u_sd_floor <- 0.01
  pars_soc <- beta_from_mean_sd(params$utilities$u_soc$well,  max(0.2*abs(params$utilities$u_soc$well), u_sd_floor))
  pars_glx <- beta_from_mean_sd(params$utilities$u_glx$well,  max(0.2*abs(params$utilities$u_glx$well), u_sd_floor))
  pars_duc <- beta_from_mean_sd(params$utilities$du_complication, max(0.2*abs(params$utilities$du_complication), u_sd_floor))
  
  u_soc_well <- rbeta(1, pars_soc["shape1"], pars_soc["shape2"])
  u_glx_well <- rbeta(1, pars_glx["shape1"], pars_glx["shape2"])
  du_comp    <- rbeta(1, pars_duc["shape1"], pars_duc["shape2"])
  
  u_soc_well <- safel(u_soc_well, params$utilities$lower, params$utilities$upper)
  u_glx_well <- safel(u_glx_well, params$utilities$lower, params$utilities$upper)
  # (Complication utilities are computed inside run_ce_with_overrides via du_complication,
  #  but we pass du_comp itself as the override so both arms' comp utilities are updated.)
  
  # --- SoC complication rate ~ lognormal around baseline (CV 20%) ---
  ln_soc   <- ln_from_mean_cv(base_soc_rate, 0.2)
  soc_rate <- rlnorm(1, meanlog = ln_soc["meanlog"], sdlog = ln_soc["sdlog"])
  soc_prob <- 1 - exp(-soc_rate)
  soc_mult <- soc_prob / base_soc_prob   # ratio to scale probabilities (≈ rates) in builder
  
  # --- Glucorix complication hazards multiplier (CV 20%) ---
  ln_mult      <- ln_from_mean_cv(1, 0.2)
  mult_glx_wc  <- rlnorm(1, meanlog = ln_mult["meanlog"], sdlog = ln_mult["sdlog"])
  
  options(glx_wc_mult_psa = mult_glx_wc)
  
  # --- Excess mortality multiplier (CV 20%) ---
  mult_excess <- rlnorm(1, meanlog = ln_mult["meanlog"], sdlog = ln_mult["sdlog"])
  
  # --- HR Complication->Death (lognormal) ---
  base_hr <- ifelse(is.null(params$hazards$hr_comp), 1.5, params$hazards$hr_comp)
  ln_hr   <- ln_from_mean_cv(base_hr, 0.2)
  hr_comp <- rlnorm(1, meanlog = ln_hr["meanlog"], sdlog = ln_hr["sdlog"])
  
  # --- Run CEA with sampled inputs ---
  res <- run_ce_with_overrides(list(
    cost_well          = cost_well,
    cost_complication  = cost_complication,
    drug_cost_glucorix = drug_cost_glucorix,
    u_soc_well         = u_soc_well,
    u_glx_well         = u_glx_well,
    du_complication    = du_comp,
    soc_rate_mult      = soc_mult,      # affects SoC WC
    excess_mult        = mult_excess,   # affects death hazards
    hr_comp            = hr_comp        # CD hazard multiplier
    # discount & stop_drug_after: keep base in PSA (varied in DSA)
  ), lambda = thresholds[2])
  
  c(inc_cost = unname(res["inc_cost"]),
    inc_qaly = unname(res["inc_qaly"]))
}


# default Glucorix WC multiplier = 1 (PSA will overwrite each draw)
options(glx_wc_mult_psa = 1)


## ---- 8C_PSA_run_plots -----------------------------------------------------

set.seed(SEED + 101)
n_sim <- 5000L                # bump to 5000 for submission if time allows
ths   <- c(20000, 30000, 50000)

psa_mat <- replicate(n_sim, psa_one_draw(thresholds = ths))
psa_df  <- data.frame(sim = 1:n_sim,
                      inc_cost = psa_mat["inc_cost",],
                      inc_qaly = psa_mat["inc_qaly",])

# Save PSA results
write.csv(psa_df, file.path(outputs_dir, "psa_results.csv"), row.names = FALSE)





# --- Run base-case to get deterministic ICER point ---
res_base <- run_ce_with_overrides(base_vals, lambda = lambda_ref)

# Extract deterministic incremental cost & QALY for highlight
inc_cost_base <- as.numeric(res_base["inc_cost"])
inc_qaly_base <- as.numeric(res_base["inc_qaly"])
base_ICER     <- as.numeric(res_base["ICER"])

# --- CE plane plot with ICER highlight ---
p_scatter <- ggplot(psa_df, aes(x = inc_qaly, y = inc_cost)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(alpha = 0.4, size = 1) +
  # Add deterministic ICER point in red
  geom_point(aes(x = inc_qaly_base, y = inc_cost_base),
             color = "red", size = 3) +
  #annotate("text",
          # x = inc_qaly_base,
         #  y = inc_cost_base,
          # label = sprintf("ICER = €%s/QALY", scales::comma(round(base_ICER, 0))),
          # color = "red", vjust = -1, hjust = 0.5, size = 3.5) +
  labs(title = "Cost-Effectiveness Plane (PSA)",
       x = "Incremental QALYs (Glucorix – SoC)",
       y = "Incremental Costs (€)") +
  theme_minimal(base_size = 12)

ggsave(file.path(outputs_dir, "ce_scatter.png"), p_scatter, width = 7, height = 6, dpi = 300)



# CEAC
lam_grid  <- seq(0, 60000, by = 1000)
ceac_grid <- sapply(lam_grid, function(lam) mean(lam*psa_df$inc_qaly - psa_df$inc_cost > 0))
ceac_plot_df <- data.frame(lambda = lam_grid, prob_CE = ceac_grid)

p_ceac <- ggplot(ceac_plot_df, aes(x = lambda, y = prob_CE)) +
  geom_line(size = 1) +
  scale_x_continuous(labels = label_dollar(prefix = "€"), breaks = seq(0,60000,10000)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Cost-Effectiveness Acceptability Curve (CEAC)",
       x = "Willingness-to-pay (€/QALY)",
       y = "Pr(Glucorix is cost-effective)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "ceac.png"), p_ceac, width = 7, height = 5, dpi = 300)




#summary PSA
library(scales)

# --- PSA means ---
## ---- 8D_PSA_summaries (with mean ICER) -------------------------------------
stopifnot(exists("psa_df"), nrow(psa_df) > 0)

# 1) Mean incremental cost & QALYs
mean_inc_cost  <- mean(psa_df$inc_cost, na.rm = TRUE)
mean_inc_qaly  <- mean(psa_df$inc_qaly, na.rm = TRUE)

# 2) ICERs across PSA draws
icer_vec <- with(psa_df, inc_cost / inc_qaly)
icer_vec <- icer_vec[is.finite(icer_vec) & inc_qaly > 0]

# 3) Mean ICER (ratio of means)
icer_mean <- mean_inc_cost / mean_inc_qaly

# 4) Median ICER (median of per-simulation ratios)
icer_med <- median(icer_vec, na.rm = TRUE)

# 5) Probability cost-effective at thresholds
ths <- c(20000, 30000, 50000)
prob_ce <- sapply(ths, function(lam) mean(lam*psa_df$inc_qaly - psa_df$inc_cost > 0))

# 6) Print nicely
cat(sprintf("• PSA mean incremental cost: €%s\n", scales::comma(round(mean_inc_cost, 0))))
cat(sprintf("• PSA mean incremental QALYs: %0.3f\n", mean_inc_qaly))
cat(sprintf("• Mean ICER from PSA (ratio of means): €%s per QALY\n", 
            scales::comma(round(icer_mean, 0))))
cat(sprintf("• Median ICER from PSA (median of ratios): €%s per QALY\n", 
            scales::comma(round(icer_med, 0))))
cat(sprintf("• Probability of cost-effectiveness at thresholds:\n"))
cat(sprintf("   o €20,000/QALY: %s\n", scales::percent(prob_ce[1], accuracy = 0.1)))
cat(sprintf("   o €30,000/QALY: %s\n", scales::percent(prob_ce[2], accuracy = 0.1)))
cat(sprintf("   o €50,000/QALY: %s\n", scales::percent(prob_ce[3], accuracy = 0.1)))

# 7) Save a compact summary table
psa_summary <- data.frame(
  metric   = c("PSA mean incremental cost",
               "PSA mean incremental QALYs",
               "Mean ICER from PSA (ratio of means)",
               "Median ICER from PSA (median of ratios)",
               "Pr(CE) at €20,000/QALY",
               "Pr(CE) at €30,000/QALY",
               "Pr(CE) at €50,000/QALY"),
  value    = c(mean_inc_cost,
               mean_inc_qaly,
               icer_mean,
               icer_med,
               prob_ce[1],
               prob_ce[2],
               prob_ce[3]),
  formatted = c(
    paste0("€", scales::comma(round(mean_inc_cost, 0))),
    sprintf("%0.3f", mean_inc_qaly),
    paste0("€", scales::comma(round(icer_mean, 0)), " per QALY"),
    paste0("€", scales::comma(round(icer_med, 0)), " per QALY"),
    scales::percent(prob_ce[1], accuracy = 0.1),
    scales::percent(prob_ce[2], accuracy = 0.1),
    scales::percent(prob_ce[3], accuracy = 0.1)
  )
)
write.csv(psa_summary, file.path(outputs_dir, "psa_summary.csv"), row.names = FALSE)





## ---- 8D_sanity_tests -------------------------------------------------------
# 1) Drug cost to 0 => should lower incremental cost (maybe dominance)
res_zero_drug <- run_ce_with_overrides(list(drug_cost_glucorix = 0))
print(res_zero_drug[c("inc_cost","inc_qaly","ICER")])
# 3) Hazard -> 0: set SoC & Glucorix complication hazards extremely low (approx everyone stays Well)
#    Do this by scaling SoC prob tiny and Glucorix hazards tiny
res_haz0 <- run_ce_with_overrides(list(soc_rate_mult = 1e-6, excess_mult = 1e-6))
print(res_haz0[c("inc_cost","inc_qaly","ICER")])

# 3) HR_comp = 1.5 => Complication->Death equal to Well->Death hazard
res_hr1 <- run_ce_with_overrides(list(hr_comp = 1.5))
print(res_hr1[c("inc_cost","inc_qaly","ICER")])







## ---- 10A_export_tables -----------------------------------------------------
# Base-case table was saved earlier as basecase_step5.csv
# Also export a tidy base-case table (rounded)
basecase <- read.csv(file.path(outputs_dir, "basecase_step5.csv"))
basecase$value_round <- round(basecase$value, 4)
write.csv(basecase, file.path(outputs_dir, "basecase_final.csv"), row.names = FALSE)

# Ensure DSA, PSA, and figures exist (created above)
stopifnot(file.exists(file.path(outputs_dir, "psa_results.csv")))
stopifnot(file.exists(file.path(outputs_dir, "ce_scatter.png")))
stopifnot(file.exists(file.path(outputs_dir, "ceac.png")))
stopifnot(file.exists(file.path(outputs_dir, "tornado_INB.png")))




## ---- 10B_save_master_script -----------------------------------------------
# Write an executable script that sources nothing else (lightweight wrapper).
master_script <- '
# CEA_Glucorix.R — master runner
# Usage: source("CEA_Glucorix.R") from project root or RStudio

globals <- readRDS(file.path("Glucorix_CEA_IGES","Outputs","globals.rds"))
attach(globals)
message("Project: ", proj_dir)

# Re-run PSA & CE plots quickly
psa_df <- read.csv(file.path(outputs_dir, "psa_results.csv"))
library(ggplot2); library(scales)

p_scatter <- ggplot(psa_df, aes(x = inc_qaly, y = inc_cost)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(alpha = 0.4, size = 1) +
  labs(title = "Cost-Effectiveness Plane (PSA)",
       x = "Incremental QALYs",
       y = "Incremental Costs (€)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "ce_scatter.png"), p_scatter, width = 7, height = 6, dpi = 300)

lam_grid <- seq(0, 60000, by = 1000)
ceac_grid <- sapply(lam_grid, function(lam) mean(lam*psa_df$inc_qaly - psa_df$inc_cost > 0))
ceac_plot_df <- data.frame(lambda = lam_grid, prob_CE = ceac_grid)
p_ceac <- ggplot(ceac_plot_df, aes(x = lambda, y = prob_CE)) +
  geom_line(size = 1) +
  scale_x_continuous(labels = label_dollar(prefix = "€"), breaks = seq(0,60000,10000)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "CEAC", x = "WTP (€/QALY)", y = "Pr(CE)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outputs_dir, "ceac.png"), p_ceac, width = 7, height = 5, dpi = 300)

message("CEA_Glucorix.R finished. Outputs in: ", outputs_dir)
'
writeLines(master_script, con = file.path(outputs_dir, "CEA_Glucorix.R"))





## ---- 10C_readme ------------------------------------------------------------
readme <- c(
  "# Glucorix_CEA_IGES – HOW TO RUN",
  "",
  "## Prerequisites",
  "- R 4.2+; packages: readxl, heemod, BCEA, flexsurv, data.table, ggplot2, mvtnorm, tidyr, patchwork, scales, triangle, survival.",
  "",
  "## Data",
  "- Place `Data/glucorix_ipd.csv` and `Data/lifetable.xlsx` in the project.",
  "",
  "## Running the model",
  "1. Execute Steps 1–7 to build inputs, hazards, and the Markov runner.",
  "2. For sensitivity analyses and figures:",
  "   - Run the Step 8 code (DSA + PSA).",
  "   - Run the Step 9 code (diagnostics).",
  "3. Optional: run `Outputs/CEA_Glucorix.R` to regenerate CE plots from saved PSA.",
  "",
  "## Outputs (in `Outputs/`)",
  "- `basecase_final.csv` – base-case totals and ICER.",
  "- `per_cycle_SoC.csv`, `per_cycle_Glucorix.csv` – cycle-level results.",
  "- `tornado_INB.png` – tornado diagram (INB at λ = €30,000).",
  "- `psa_results.csv` – PSA Δcost and ΔQALY per simulation.",
  "- `ce_scatter.png` – CE plane.",
  "- `ceac.png` – CE acceptability curve.",
  "- Diagnostics: `diagnostic_state_occupancy.png`, `diagnostic_age_trajectory.png`,",
  "  `diagnostic_km_vs_parametric_complication.png`.",
  "",
  "## Notes",
  "- DSA uses ±20% ranges; discount (0%/5%), and a treatment discontinuation scenario (5 years).",
  "- PSA uses: Gamma (costs), Beta (utilities), Lognormal (rates/HRs), with 20% CV where SEs are unknown.",
  "- Competing risks handled in hazards; half-cycle correction applied to costs & QALYs; 3% discount base case.",
  ""
)
writeLines(readme, con = file.path(outputs_dir, "README.md"))

