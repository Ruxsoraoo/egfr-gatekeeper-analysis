# ============================================================
# Model 1 - reversible inhibitor binding with ATP competition
#
# Purpose
#   Relate free inhibitor concentration to fractional receptor
#   occupancy when ATP competes for the same site.
#
#   Erlotinib is ATP-competitive, and T790M acts mainly by
#   raising the receptor's affinity for ATP (Yun et al. 2008),
#   not only by steric hindrance. A binding model without ATP
#   therefore cannot represent the resistance mechanism, so the
#   competition is present from the first model rather than
#   added later.
#
# Scope
#   Binding only. No receptor synthesis, degradation or
#   downstream signal: those are models 2 and 3. Free drug and
#   free ATP are treated as fixed inputs, not state variables.
#
# Run with:  Rscript R/01_binding_atp_competition.R
# Requires:  deSolve
#
# Output
#   figures/01_occupancy_vs_dose.png
#   Printed check of the numerical solution against the
#   closed-form steady state.
#
# Axes are dimensionless (D/Kd and A/Ka) so that this script
# runs before the parameter table is filled in. Model 4
# substitutes measured wild-type and T790M values from
# data/parameters.csv.
# ============================================================

library(deSolve)


# ---- 1. Equations ----------------------------------------
# State variables, all concentrations:
#   R  free receptor
#   DR inhibitor-receptor complex
#   AR ATP-receptor complex
#
# Mass action in both directions. D (free inhibitor) and
# A (free ATP) are parameters, not states: the drug is assumed
# buffered by the circulation and ATP by cellular metabolism,
# so binding does not deplete either. This is the standard
# assumption and it is also the model's first limitation.
#
#   dR/dt  = -konD*D*R + koffD*DR - konA*A*R + koffA*AR
#   dDR/dt =  konD*D*R - koffD*DR
#   dAR/dt =  konA*A*R - koffA*AR
#
# Total receptor is conserved: R + DR + AR is constant.

binding_atp <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR  <- -konD * D * R + koffD * DR - konA * A * R + koffA * AR
    dDR <-  konD * D * R - koffD * DR
    dAR <-  konA * A * R - koffA * AR
    list(c(dR, dDR, dAR))
  })
}


# ---- 2. Closed-form steady state -------------------------
# At equilibrium Kd = koffD/konD and Ka = koffA/konA, and
# fractional inhibitor occupancy is
#
#   DR / (R + DR + AR) = (D/Kd) / (1 + D/Kd + A/Ka)
#
# The ATP term A/Ka sits in the denominator, so raising ATP
# affinity (lowering Ka) lowers occupancy at the same dose.
# That is the mechanism in one line, and it is why the
# numerical solution below is checked against this expression
# rather than trusted on its own.

occupancy_analytic <- function(D_over_Kd, A_over_Ka) {
  D_over_Kd / (1 + D_over_Kd + A_over_Ka)
}


# ---- 3. Numerical solution -------------------------------
# Solved with lsoda. Rate constants are set to 1 in
# dimensionless form so that Kd = Ka = 1 and the x axis reads
# directly as a multiple of Kd. Only the ratios matter for
# the steady state; the absolute rates set how fast it is
# reached, which models 2 and 3 use.
#
# The system is integrated to t = 1e4, far beyond equilibrium,
# and the final value is taken as the steady state.

occupancy_numeric <- function(D_over_Kd, A_over_Ka, R_total = 1, t_end = 1e4) {
  parms <- c(konD = 1, koffD = 1, konA = 1, koffA = 1,
             D = D_over_Kd, A = A_over_Ka)
  state <- c(R = R_total, DR = 0, AR = 0)
  out <- ode(y = state, times = c(0, t_end), func = binding_atp,
             parms = parms, method = "lsoda")
  final <- out[nrow(out), ]
  final[["DR"]] / (final[["R"]] + final[["DR"]] + final[["AR"]])
}


# ---- 4. Check the solver against the algebra -------------
# A numerical result that is not checked against something
# independent is not evidence. If these disagree, the
# equations and the closed form have diverged and nothing
# below should be believed.

dose_grid <- 10^seq(-2, 3, length.out = 60)
atp_levels <- c(0, 1, 10, 100)

cat("Checking numerical solution against closed form\n")
max_err <- 0
for (a in atp_levels) {
  for (d in dose_grid) {
    err <- abs(occupancy_numeric(d, a) - occupancy_analytic(d, a))
    max_err <- max(max_err, err)
  }
}
cat(sprintf("  maximum absolute difference: %.2e\n", max_err))
stopifnot(max_err < 1e-6)
cat("  agreement within tolerance\n\n")


# ---- 5. What dose restores occupancy? --------------------
# Rearranging the closed form, the dose giving occupancy f is
#
#   D/Kd = f * (1 + A/Ka) / (1 - f)
#
# so the dose needed scales linearly with 1 + A/Ka. This is
# the quantity model 4 evaluates with measured values, and
# the reason the ATP term cannot be left out of the dose
# argument.

dose_for_occupancy <- function(f, A_over_Ka) {
  f * (1 + A_over_Ka) / (1 - f)
}

cat("Dose (as a multiple of Kd) for 50% occupancy\n")
for (a in atp_levels) {
  cat(sprintf("  A/Ka = %6.1f : D/Kd = %8.2f\n", a, dose_for_occupancy(0.5, a)))
}
cat("\n")


# ---- 6. Figure -------------------------------------------
# Occupancy against dose at several ATP levels, log x axis.
# Raising A/Ka shifts the curve right without changing its
# shape: the same occupancy is still reachable in principle,
# but only at a higher concentration.

dir.create("figures", showWarnings = FALSE)
png("figures/01_occupancy_vs_dose.png", width = 1600, height = 1200, res = 200)
plot(NA, xlim = range(dose_grid), ylim = c(0, 1), log = "x",
     xlab = "Free inhibitor concentration (multiples of Kd)",
     ylab = "Fractional receptor occupancy",
     main = "Inhibitor occupancy under ATP competition")
palette_atp <- c("black", "#1b6ca8", "#c05621", "#8b2252")
for (i in seq_along(atp_levels)) {
  lines(dose_grid, occupancy_analytic(dose_grid, atp_levels[i]),
        col = palette_atp[i], lwd = 2)
}
abline(h = 0.9, lty = 3, col = "grey40")
legend("topleft", legend = sprintf("A/Ka = %g", atp_levels),
       col = palette_atp, lwd = 2, bty = "n")
invisible(dev.off())

cat("Wrote figures/01_occupancy_vs_dose.png\n")


# ============================================================
# NEXT
# Model 2 adds receptor synthesis and degradation.
# Model 4 substitutes measured Kd and ATP Km values for
# wild-type and T790M from data/parameters.csv, which must be
# filled from published measurements first. No affinity value
# is assumed or computed anywhere in this script.
# ============================================================
