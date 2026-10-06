# ============================================================
# Model 4 - measured affinities, and the dose question
#
# Questions
#   1. What erlotinib concentration restores mutant occupancy
#      to the level reached in the drug-sensitive mutant?
#   2. Is that concentration plausible?
#   3. Does the competitive model reproduce the measured
#      potency shift?
#
# This is the first script that uses measured values. They are
# read from data/parameters.csv, never written here, so the
# numbers in the output can always be traced to a row with a
# citation.
#
# A note on which variants are compared
#   The source paper measures L858R and L858R/T790M, not
#   wild-type against T790M. T790M arises clinically on an
#   already-mutant receptor, so L858R against L858R/T790M is
#   the comparison the data support. The wild-type ATP value
#   is reported too and is carried along for reference.
#
# Run with:  Rscript R/04_mutant_dose.R
# Requires:  the Ki, Km and IC50 rows of data/parameters.csv
#
# Output
#   figures/04_mutant_dose.png
# ============================================================


# ---- 1. Read the parameter table -------------------------
# Values live in the CSV with their citations. Any row still
# marked unverified is reported loudly, because an unchecked
# number that reaches a figure is worse than a missing one.

params <- read.csv("data/parameters.csv", stringsAsFactors = FALSE)

get_value <- function(symbol) {
  row <- params[params$symbol == symbol, ]
  if (nrow(row) != 1) {
    stop(sprintf("expected exactly one row for '%s', found %d", symbol, nrow(row)))
  }
  if (is.na(row$value) || row$value == "") {
    stop(sprintf(paste("no value for '%s' in data/parameters.csv.",
                       "Fill it from a published measurement and cite it."),
                 symbol))
  }
  as.numeric(row$value)
}

needed <- c("Ki_erlotinib_L858R", "Ki_erlotinib_L858R_T790M",
            "Km_ATP_L858R", "Km_ATP_L858R_T790M", "Km_ATP_WT",
            "A_assay", "IC50_erlotinib_L858R", "IC50_erlotinib_L858R_T790M")

unverified <- params$symbol[params$symbol %in% needed & params$verified != "yes"]
if (length(unverified) > 0) {
  cat("------------------------------------------------------------\n")
  cat("WARNING: these values have not been checked against the paper\n")
  for (s in unverified) cat(sprintf("  %s\n", s))
  cat("Every number below inherits that caveat. Set verified = yes\n")
  cat("in data/parameters.csv once each has been read off the source.\n")
  cat("------------------------------------------------------------\n\n")
}

Ki_sens <- get_value("Ki_erlotinib_L858R")
Ki_res  <- get_value("Ki_erlotinib_L858R_T790M")
Km_sens <- get_value("Km_ATP_L858R")
Km_res  <- get_value("Km_ATP_L858R_T790M")
A_assay <- get_value("A_assay") * 1000   # mM to uM, to match Km


# ---- 2. Occupancy and dose -------------------------------
# From model 1, with Ki in place of Kd and the ATP Michaelis
# constant in place of Ka. Treating Km,ATP as an apparent
# dissociation constant for ATP is an approximation: Km equals
# Kd only when catalysis is slow next to dissociation. It is
# the standard approximation and it is also a limitation to
# state in the write-up.

occupancy <- function(D, Ki, Km, A) {
  (D / Ki) / (1 + D / Ki + A / Km)
}

dose_for <- function(f, Ki, Km, A) {
  Ki * f * (1 + A / Km) / (1 - f)
}


# ---- 3. The dose answer ----------------------------------

cat("ATP competition term at the assay ATP concentration\n")
cat(sprintf("  L858R        : A/Km = %5.1f  (Km = %.0f uM)\n", A_assay / Km_sens, Km_sens))
cat(sprintf("  L858R/T790M  : A/Km = %5.1f  (Km = %.0f uM)\n", A_assay / Km_res, Km_res))
cat(sprintf("  wild-type    : A/Km = %5.1f  (Km = %.0f uM, reference only)\n\n",
            A_assay / get_value("Km_ATP_WT"), get_value("Km_ATP_WT")))

cat("Erlotinib concentration for 90% occupancy\n")
d_sens <- dose_for(0.9, Ki_sens, Km_sens, A_assay)
d_res  <- dose_for(0.9, Ki_res,  Km_res,  A_assay)
cat(sprintf("  L858R        : %8.0f nM\n", d_sens))
cat(sprintf("  L858R/T790M  : %8.0f nM\n", d_res))
cat(sprintf("  ratio        : %8.1f fold\n\n", d_res / d_sens))

cat("Where that factor comes from\n")
cat(sprintf("  affinity (Ki)      : %5.1f fold\n", Ki_res / Ki_sens))
cat(sprintf("  ATP competition    : %5.2f fold\n",
            (1 + A_assay / Km_res) / (1 + A_assay / Km_sens)))
cat("  In this dataset the shift is carried almost entirely by the\n")
cat("  inhibition constant. The ATP term contributes little, because\n")
cat("  the double mutant's ATP affinity is only modestly tighter.\n")
cat("  Yun et al. 2008 report a much larger ATP shift; the sensitivity\n")
cat("  analysis below is where that disagreement is handled, and it\n")
cat("  needs the Yun values in the CSV first.\n\n")


# ---- 4. Does the model match the measurement? ------------
# The model's prediction of an IC50 is the concentration
# giving 50% occupancy, Ki*(1 + A/Km). The paper measured
# IC50 at the same ATP concentration, so the two are directly
# comparable. This is the only real test in the project of
# whether the model describes anything.

cat("Model against measurement, IC50 at the assay ATP concentration\n")
pred_sens <- Ki_sens * (1 + A_assay / Km_sens)
pred_res  <- Ki_res  * (1 + A_assay / Km_res)
meas_sens <- get_value("IC50_erlotinib_L858R")
meas_res  <- get_value("IC50_erlotinib_L858R_T790M")
cat(sprintf("  L858R        : predicted %7.0f nM   measured %7.0f nM\n", pred_sens, meas_sens))
cat(sprintf("  L858R/T790M  : predicted %7.0f nM   measured %7.0f nM\n", pred_res, meas_res))
cat(sprintf("  ratio        : predicted %7.1f      measured %7.1f\n\n",
            pred_res / pred_sens, meas_res / meas_sens))

cat("Reading this honestly\n")
cat("  For L858R the prediction lands close to the measured value.\n")
cat("  For the double mutant it is far too low: the model accounts\n")
cat("  for only part of the measured loss of potency. Simple\n")
cat("  competitive binding with these two constants is therefore not\n")
cat("  sufficient, and the gap is the interesting result rather than a\n")
cat("  defect to hide. Candidate explanations, none tested here:\n")
cat("  Km,ATP is not the ATP Kd; the Ki and IC50 measurements come\n")
cat("  from different assay formats; or the mutant's resistance\n")
cat("  involves something this model omits.\n\n")


# ---- 5. Is the dose plausible? ---------------------------
# Deliberately not answered with a number yet. The comparison
# needs the free plasma concentration of erlotinib at the
# standard dose, which means both C_ss_erlotinib and
# fu_erlotinib, and both are still empty rows in the CSV.
# Printing a dose in nM next to nothing invites the reader to
# supply their own comparison from memory, so the script says
# what is missing instead.

missing <- c("C_ss_erlotinib", "fu_erlotinib")[
  vapply(c("C_ss_erlotinib", "fu_erlotinib"), function(s) {
    v <- params$value[params$symbol == s]
    length(v) != 1 || is.na(v) || v == ""
  }, logical(1))]

if (length(missing) > 0) {
  cat("Plausibility of the dose: not yet assessable.\n")
  cat("  Missing from data/parameters.csv:", paste(missing, collapse = ", "), "\n")
  cat("  A concentration in nM cannot be called achievable or not\n")
  cat("  until the free plasma concentration at the licensed dose is\n")
  cat("  in the table with its citation.\n\n")
}


# ---- 6. Figure -------------------------------------------

dose_grid <- 10^seq(0, 5, length.out = 200)
dir.create("figures", showWarnings = FALSE)
png("figures/04_mutant_dose.png", width = 1600, height = 1200, res = 200)
plot(NA, xlim = range(dose_grid), ylim = c(0, 1), log = "x",
     xlab = "Free erlotinib concentration (nM)",
     ylab = "Fractional receptor occupancy",
     main = "Erlotinib occupancy, L858R against L858R/T790M")
lines(dose_grid, occupancy(dose_grid, Ki_sens, Km_sens, A_assay),
      col = "#1b6ca8", lwd = 2)
lines(dose_grid, occupancy(dose_grid, Ki_res, Km_res, A_assay),
      col = "#c05621", lwd = 2)
abline(h = 0.9, lty = 3, col = "grey40")
points(c(d_sens, d_res), c(0.9, 0.9), pch = 19,
       col = c("#1b6ca8", "#c05621"))
legend("topleft",
       legend = c(sprintf("L858R (Ki = %.1f nM)", Ki_sens),
                  sprintf("L858R/T790M (Ki = %.0f nM)", Ki_res),
                  "90% occupancy"),
       col = c("#1b6ca8", "#c05621", "grey40"),
       lwd = 2, lty = c(1, 1, 3), bty = "n")
mtext(sprintf("ATP %.0f uM. Values from data/parameters.csv%s",
              A_assay,
              if (length(unverified) > 0) " (unverified)" else ""),
      side = 1, line = 4, cex = 0.7)
invisible(dev.off())

cat("Wrote figures/04_mutant_dose.png\n")


# ============================================================
# STILL OPEN
# - Confirm every value against the source paper and set
#   verified = yes in data/parameters.csv.
# - Add the Yun 2008 ATP values and repeat the dose
#   calculation under both labs' ATP affinities, as the
#   sensitivity analysis the plan calls for.
# - Add the free plasma concentration before making any claim
#   about whether a dose is tolerable.
# - The predicted and measured IC50 ratios disagree by about
#   tenfold for the double mutant. That gap is a result, and
#   the write-up should present it as one.
# ============================================================
