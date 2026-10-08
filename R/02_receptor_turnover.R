# ============================================================
# Model 2 - binding with ATP competition and receptor turnover
#
# Question
#   How does receptor turnover change steady-state occupancy?
#
# Model 1 held the receptor pool fixed. A real receptor is
# synthesised and degraded, and the drug-bound complex is
# internalised. This model adds those three processes and asks
# what difference they make.
#
# Scope
#   No downstream signal: that is model 3. Free drug and free
#   ATP remain fixed inputs.
#
# Run with:  Rscript R/02_receptor_turnover.R
# Requires:  deSolve
#
# Output
#   figures/02_turnover_shifts_apparent_kd.png
#   Printed check of the numerical solution against the
#   closed-form steady state.
#
# Rates are dimensionless, as in model 1. Model 4 substitutes
# measured values from data/parameters.csv.
# ============================================================

library(deSolve)


# ---- 1. Equations ----------------------------------------
# States: R free receptor, DR inhibitor complex,
#         AR ATP complex.
#
#   dR/dt  = ksyn - kdeg*R - konD*D*R + koffD*DR
#                          - konA*A*R + koffA*AR
#   dDR/dt = konD*D*R - koffD*DR - kint*DR
#   dAR/dt = konA*A*R - koffA*AR - kdeg*AR
#
# Three additions to model 1:
#   ksyn  zero-order synthesis of new receptor
#   kdeg  first-order degradation of free receptor, applied
#         to the ATP complex as well, since ATP binding does
#         not protect the receptor
#   kint  loss of the inhibitor-bound complex, by
#         internalisation or degradation
#
# Total receptor is no longer conserved. It reaches a steady
# state set by synthesis against the three loss routes, so the
# conservation check used in model 1 does not apply here and
# the closed form below takes its place.

turnover <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR  <- ksyn - kdeg * R - konD * D * R + koffD * DR -
           konA * A * R + koffA * AR
    dDR <- konD * D * R - koffD * DR - kint * DR
    dAR <- konA * A * R - koffA * AR - kdeg * AR
    list(c(dR, dDR, dAR))
  })
}


# ---- 2. Closed-form steady state -------------------------
# Setting the derivatives to zero and solving:
#
#   Kd_app = (koffD + kint) / konD
#   Ka_app = (koffA + kdeg) / konA
#   R  = ksyn / (kdeg + kint*D/Kd_app + kdeg*A/Ka_app)
#   DR = R * D/Kd_app
#   AR = R * A/Ka_app
#
# and therefore
#
#   occupancy = (D/Kd_app) / (1 + D/Kd_app + A/Ka_app)
#
# This is the result worth pausing on. The expression has the
# same form as model 1, but the dissociation constants are
# replaced by apparent ones that include the loss rates. Loss
# of the complex acts exactly like faster dissociation.
#
# So turnover does not change the shape of the occupancy curve.
# It shifts it, by however much kint adds to koffD. If kint is
# small next to koffD the shift is negligible; if the complex
# is cleared as fast as it dissociates, the apparent Kd
# doubles.

steady_state <- function(parms) {
  with(as.list(parms), {
    Kd_app <- (koffD + kint) / konD
    Ka_app <- (koffA + kdeg) / konA
    R  <- ksyn / (kdeg + kint * D / Kd_app + kdeg * A / Ka_app)
    c(R = R, DR = R * D / Kd_app, AR = R * A / Ka_app,
      Kd_app = Kd_app, Ka_app = Ka_app)
  })
}

occupancy_analytic <- function(parms) {
  ss <- steady_state(parms)
  ss[["DR"]] / (ss[["R"]] + ss[["DR"]] + ss[["AR"]])
}


# ---- 3. Numerical solution -------------------------------
# Integrated from an empty cell to t = 1e6, long enough for
# synthesis and loss to balance. Starting from zero receptor
# rather than a guessed steady state means the approach to
# equilibrium is also visible if it is ever plotted.

occupancy_numeric <- function(parms, t_end = 1e6) {
  out <- ode(y = c(R = 0, DR = 0, AR = 0), times = c(0, t_end),
             func = turnover, parms = parms, method = "lsoda")
  final <- out[nrow(out), ]
  final[["DR"]] / (final[["R"]] + final[["DR"]] + final[["AR"]])
}


# ---- 4. Check the solver against the algebra -------------
# Over a grid of rate combinations, not a single case: a check
# at one parameter set can pass by coincidence.

base <- c(konD = 1, koffD = 1, konA = 1, koffA = 1,
          ksyn = 1, kdeg = 0.1, kint = 0.1, D = 1, A = 10)

cat("Checking numerical solution against closed form\n")
max_err <- 0
for (d in 10^seq(-2, 2, length.out = 12)) {
  for (ki in c(0.01, 0.1, 1)) {
    for (a in c(0, 10, 100)) {
      parms <- base
      parms[["D"]] <- d; parms[["kint"]] <- ki; parms[["A"]] <- a
      max_err <- max(max_err,
                     abs(occupancy_numeric(parms) - occupancy_analytic(parms)))
    }
  }
}
cat(sprintf("  maximum absolute difference: %.2e\n", max_err))
stopifnot(max_err < 1e-6)
cat("  agreement within tolerance\n\n")


# ---- 5. How much does turnover matter? -------------------
# The answer is the ratio Kd_app / Kd, which is
# (koffD + kint) / koffD. Printed for three rates of complex
# loss so the size of the effect is explicit rather than
# asserted.

cat("Effect of complex loss on the apparent dissociation constant\n")
for (ki in c(0.01, 0.1, 1)) {
  parms <- base
  parms[["kint"]] <- ki
  ss <- steady_state(parms)
  cat(sprintf("  kint = %5.2f (koffD = %.2f) : Kd_app/Kd = %.2f\n",
              ki, parms[["koffD"]], ss[["Kd_app"]] / parms[["koffD"]]))
}
cat("\n")


# ---- 6. Figure -------------------------------------------
# Occupancy against dose with and without turnover, at fixed
# ATP. The curves are parallel on a log axis, which is the
# point: turnover costs dose, not achievable occupancy.

dose_grid <- 10^seq(-2, 3, length.out = 60)
dir.create("figures", showWarnings = FALSE)
png("figures/02_turnover_shifts_apparent_kd.png",
    width = 1600, height = 1200, res = 200)
plot(NA, xlim = range(dose_grid), ylim = c(0, 1), log = "x",
     xlab = "Free inhibitor concentration (multiples of Kd)",
     ylab = "Fractional receptor occupancy",
     main = "Receptor turnover shifts the curve, not its shape")
kint_levels <- c(0, 0.1, 1)
palette_kint <- c("black", "#1b6ca8", "#c05621")
for (i in seq_along(kint_levels)) {
  parms <- base
  parms[["kint"]] <- kint_levels[i]
  occ <- vapply(dose_grid, function(d) {
    p <- parms; p[["D"]] <- d; occupancy_analytic(p)
  }, numeric(1))
  lines(dose_grid, occ, col = palette_kint[i], lwd = 2)
}
legend("topleft", legend = sprintf("kint = %g", kint_levels),
       col = palette_kint, lwd = 2, bty = "n",
       title = sprintf("A/Ka = %g, koffD = %g", base[["A"]], base[["koffD"]]))
invisible(dev.off())

cat("Wrote figures/02_turnover_shifts_apparent_kd.png\n")


# ============================================================
# ANSWER TO THE QUESTION
# Turnover changes steady-state occupancy only through the
# apparent dissociation constant, Kd_app = (koffD + kint)/konD.
# The functional form is unchanged, so the curve shifts along
# the dose axis without changing shape. Whether the shift
# matters is an empirical question about kint against koffD,
# and both are rows in data/parameters.csv.
#
# NEXT
# Model 3 adds a downstream signal and asks whether 90%
# occupancy produces a 90% loss of signal.
# ============================================================
