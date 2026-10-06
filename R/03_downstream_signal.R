# ============================================================
# Model 3 - turnover with a downstream signal
#
# Question
#   Does 90% receptor occupancy produce a 90% loss of signal?
#
# Occupancy is what a binding model predicts, but it is not
# what a patient experiences. What matters is the signal the
# receptor transmits. This model adds that signal and shows
# that the two numbers are not interchangeable.
#
# Scope
#   Free drug and free ATP remain fixed inputs. The signal is
#   a single lumped variable, not a mapped pathway.
#
# Run with:  Rscript R/03_downstream_signal.R
# Requires:  deSolve
#
# Output
#   figures/03_occupancy_vs_signal.png
#   Printed check of the numerical solution against the
#   closed-form steady state.
# ============================================================

library(deSolve)


# ---- 1. Equations ----------------------------------------
# Model 2 plus one state:
#
#   dS/dt = ktrans*R - kout*S
#
# The signal is driven by FREE receptor, because an
# inhibitor-bound receptor is inactive. ATP-bound receptor is
# counted as free for this purpose only in the sense that it
# is not inhibited; treating it as signalling would require a
# catalytic step this model does not have, which is a
# simplification worth stating in the write-up.

signal_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR  <- ksyn - kdeg * R - konD * D * R + koffD * DR -
           konA * A * R + koffA * AR
    dDR <- konD * D * R - koffD * DR - kint * DR
    dAR <- konA * A * R - koffA * AR - kdeg * AR
    dS  <- ktrans * R - kout * S
    list(c(dR, dDR, dAR, dS))
  })
}


# ---- 2. Closed-form steady state -------------------------
# As model 2, with S = ktrans*R/kout.

steady_state <- function(parms) {
  with(as.list(parms), {
    Kd_app <- (koffD + kint) / konD
    Ka_app <- (koffA + kdeg) / konA
    R  <- ksyn / (kdeg + kint * D / Kd_app + kdeg * A / Ka_app)
    c(R = R, DR = R * D / Kd_app, AR = R * A / Ka_app,
      S = ktrans * R / kout)
  })
}

occupancy <- function(parms) {
  ss <- steady_state(parms)
  ss[["DR"]] / (ss[["R"]] + ss[["DR"]] + ss[["AR"]])
}

# Signal relative to the untreated steady state, which is the
# same model with D = 0.
signal_fraction <- function(parms) {
  untreated <- parms
  untreated[["D"]] <- 0
  steady_state(parms)[["S"]] / steady_state(untreated)[["S"]]
}


# ---- 3. The answer, algebraically ------------------------
# Substituting the dose that gives occupancy f, which is
# D/Kd_app = f*(1 + A/Ka_app)/(1 - f), into the expression for
# S and dividing by the untreated signal, the ATP terms cancel
# and what is left is
#
#   S/S0 = 1 / (1 + (kint/kdeg) * f/(1 - f))
#
# Three things follow, and they are the content of this model.
#
# 1. Signal loss is NOT equal to occupancy. At f = 0.9 the
#    retained signal is 52.6% when kint/kdeg = 0.1 and 1.1%
#    when kint/kdeg = 10.
# 2. The two coincide only when kint = kdeg, where
#    S/S0 = 1 - f exactly. That is an arithmetic coincidence
#    at one parameter ratio, not a general rule, and reading
#    occupancy as signal loss silently assumes it.
# 3. The ratio does not depend on ATP. ATP competition sets
#    the dose needed to reach a given occupancy, which is
#    model 1's result, but not what that occupancy does to the
#    signal. The two effects are separable, which is useful:
#    the mutation changes the first and not the second.

signal_from_occupancy <- function(f, kint_over_kdeg) {
  1 / (1 + kint_over_kdeg * f / (1 - f))
}


# ---- 4. Check the solver against the algebra -------------

base <- c(konD = 1, koffD = 1, konA = 1, koffA = 1, ksyn = 1,
          kdeg = 0.1, kint = 0.1, ktrans = 1, kout = 1,
          D = 1, A = 10)

cat("Checking numerical solution against closed form\n")
max_err <- 0
for (d in 10^seq(-2, 2, length.out = 10)) {
  for (ki in c(0.01, 0.1, 1)) {
    for (a in c(0, 10, 100)) {
      parms <- base
      parms[["D"]] <- d; parms[["kint"]] <- ki; parms[["A"]] <- a
      out <- ode(y = c(R = 0, DR = 0, AR = 0, S = 0),
                 times = c(0, 1e6), func = signal_model,
                 parms = parms, method = "lsoda")
      num <- out[nrow(out), ]
      ana <- steady_state(parms)
      max_err <- max(max_err,
                     max(abs(c(num[["R"]], num[["DR"]], num[["AR"]], num[["S"]]) -
                             c(ana[["R"]], ana[["DR"]], ana[["AR"]], ana[["S"]]))))
    }
  }
}
cat(sprintf("  maximum absolute difference: %.2e\n", max_err))
stopifnot(max_err < 1e-6)

# And check the two routes to the answer agree: the full model
# solved at a dose, against the one-line expression in part 3.
cat("Checking the closed form against the occupancy expression\n")
max_err2 <- 0
for (ki in c(0.01, 0.1, 1, 10)) {
  for (a in c(0, 10, 100)) {
    parms <- base
    parms[["kint"]] <- ki; parms[["A"]] <- a
    Kd_app <- (parms[["koffD"]] + ki) / parms[["konD"]]
    Ka_app <- (parms[["koffA"]] + parms[["kdeg"]]) / parms[["konA"]]
    for (f in c(0.5, 0.9, 0.99)) {
      parms[["D"]] <- f * (1 + a / Ka_app) / (1 - f) * Kd_app
      max_err2 <- max(max_err2,
                      abs(signal_fraction(parms) -
                          signal_from_occupancy(f, ki / parms[["kdeg"]])))
    }
  }
}
cat(sprintf("  maximum absolute difference: %.2e\n", max_err2))
stopifnot(max_err2 < 1e-9)
cat("  agreement within tolerance\n\n")


# ---- 5. The table that answers the question --------------

cat("Signal retained at 90% occupancy\n")
for (r in c(0.1, 0.5, 1, 2, 10)) {
  cat(sprintf("  kint/kdeg = %5.1f : S/S0 = %5.1f%%  (occupancy-matched loss would be 10.0%%)\n",
              r, 100 * signal_from_occupancy(0.9, r)))
}
cat("\n")


# ---- 6. Figure -------------------------------------------
# Retained signal against occupancy for several values of
# kint/kdeg, with the line S/S0 = 1 - f drawn for reference.
# Only the kint/kdeg = 1 curve lies on it.

f_grid <- seq(0.001, 0.995, length.out = 200)
ratios <- c(0.1, 1, 10)
palette_ratio <- c("#1b6ca8", "black", "#c05621")

dir.create("figures", showWarnings = FALSE)
png("figures/03_occupancy_vs_signal.png",
    width = 1600, height = 1200, res = 200)
plot(NA, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Fractional receptor occupancy",
     ylab = "Signal retained, relative to untreated",
     main = "Occupancy is not signal loss")
lines(f_grid, 1 - f_grid, lty = 2, col = "grey50", lwd = 2)
for (i in seq_along(ratios)) {
  lines(f_grid, signal_from_occupancy(f_grid, ratios[i]),
        col = palette_ratio[i], lwd = 2)
}
legend("topright",
       legend = c(sprintf("kint/kdeg = %g", ratios), "S/S0 = 1 - occupancy"),
       col = c(palette_ratio, "grey50"), lwd = 2,
       lty = c(1, 1, 1, 2), bty = "n")
invisible(dev.off())

cat("Wrote figures/03_occupancy_vs_signal.png\n")


# ============================================================
# ANSWER TO THE QUESTION
# No. At 90% occupancy the retained signal is
# 1/(1 + (kint/kdeg)*9), which is 52.6% when the complex is
# cleared ten times more slowly than free receptor is
# degraded, and 1.1% when it is cleared ten times faster.
# The two agree only when kint = kdeg.
#
# Reporting occupancy as though it were pathway inhibition
# therefore assumes a parameter ratio that has to be measured.
# kint and kdeg are rows in data/parameters.csv.
#
# NEXT
# Model 4 substitutes measured wild-type and T790M affinities
# and asks what dose restores occupancy in the mutant, and
# whether that dose is tolerable. It needs the parameter table
# filled first.
# ============================================================
