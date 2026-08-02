"""
prion_model.py
==============
A quantitative model of prion replication and PrP-lowering therapy, built to
answer one specific question:

    How deeply, and how early, must PrP-C be lowered to *halt* or *reverse*
    prion disease -- as opposed to merely slowing it down?

STRUCTURE
---------
1. Replication layer: the Masel-Jansen-Nowak (1999) nucleated polymerization
   model (NPM). Polymers of PrP-Sc elongate by recruiting PrP-C, fragment into
   new seeds, and are cleared. This layer has an exact, analytic threshold:
   below a critical substrate concentration x_crit, the prion cannot sustain
   itself and existing load DECLINES. That threshold is the mathematical
   definition of "cure by substrate reduction."

2. Toxicity layer: deliberately SEPARATED from prion load, following
   Mallucci 2003/2007 (Science 302:871; Neuron 53:325) and Sandberg 2011
   (Nature 470:540). Damage is driven by the *neuronal conversion flux*, not by
   the standing PrP-Sc burden. Dysfunction (D) is reversible and self-repairing;
   only sustained dysfunction above a threshold converts into irreversible
   neuron loss (N). This is what makes "reversal" a well-posed question.

REFERENCES FOR THE CALIBRATION TARGETS
--------------------------------------
  Buelet 1993/1994    Prnp-/- mice never develop disease; Prnp+/- roughly 2x incubation
  Fischer 1996        Tga20 (~8-10x PrP) incubation shortened to ~60-70 d
  Mallucci 2003/2007  Neuronal PrP depletion mid-infection reverses spongiosis
                      and behavioural/electrophysiological deficits
  Minikel 2020        PrP-lowering ASOs: ~50% lowering -> up to ~3x survival
  An/Davis 2025       Base editing, ~50% PrP reduction -> +52% lifespan
                      (humanized PRNP mice, sCJD MM1 and E200K isolates)
  Gentile 2026        Divalent siRNA: 49% residual PrP -> 2.7x survival
                      (chronic presymptomatic) or +64% (single dose at onset);
                      17% residual PrP achievable with a single 348 ug dose

IMPORTANT CAVEAT
----------------
This is a hypothesis-generating toy model with lumped parameters, not a
validated simulator. Its value is that it makes the threshold question
quantitative and falsifiable, and it says exactly which mouse experiment
would pin down the answer. Treat every number as a prediction to be tested.

Author's note: parameters are in units of days.
"""

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------

@dataclass
class PrionParams:
    """Lumped parameters for the replication + toxicity model.

    Replication (NPM):
        x0    normalised PrP-C steady state in an untreated animal (== 1.0)
        d     PrP-C turnover rate (1/day). PrP-C half-life ~3-5 h in neurons.
        beta  elongation rate constant (per polymer, per unit PrP-C, per day)
        b     fragmentation rate per peptide bond per day
        a     polymer clearance rate per day
        n     minimum stable polymer size (subcritical fragments dissolve)

    Toxicity:
        kappa  dysfunction generated per unit conversion flux
        rho    repair rate of reversible dysfunction (1/day)
        D_tox  dysfunction level above which irreversible neuron loss begins
        mu     rate of neuron loss above threshold
        D_sx   dysfunction level at which clinical signs appear
        N_death fraction of neurons lost that is terminal
    """
    x0: float = 1.0
    d: float = 4.0
    beta: float = 6.0
    b: float = 2.0e-4
    a: float = 0.02
    n: int = 6

    kappa: float = 1.0
    rho: float = 0.06
    D_tox: float = 0.35
    mu: float = 0.02
    D_sx: float = 0.30
    N_death: float = 0.30

    # inoculum: number of seeds delivered by intracerebral injection
    y_inoculum: float = 1e-6


# ----------------------------------------------------------------------------
# 1. Analytic replication layer
# ----------------------------------------------------------------------------

def growth_rate(x, p: PrionParams):
    """Exponential growth rate r of prion load at fixed PrP-C level x.

    Linearising the NPM around zero prion load gives the 2x2 system

        d/dt [y, z] = [[-(a + b(2n-1)),  b        ],   [y]
                       [ beta*x - b*n(n-1), -a    ]] * [z]

    r is its dominant eigenvalue. r > 0 means the prion replicates;
    r < 0 means existing prion load is cleared faster than it is made.
    """
    A = p.a + p.b * (2 * p.n - 1)
    B = p.b
    C = p.beta * x - p.b * p.n * (p.n - 1)
    Dd = p.a
    disc = (A - Dd) ** 2 + 4.0 * B * C
    if disc < 0:
        return -(A + Dd) / 2.0
    return (-(A + Dd) + np.sqrt(disc)) / 2.0


def x_crit(p: PrionParams):
    """Critical PrP-C level: the substrate concentration below which prion
    replication cannot be sustained (r(x) = 0).

        x_crit = [a*(a + b*(2n-1)) + b^2 * n * (n-1)] / (b * beta)

    Expressed as a fraction of the untreated level x0, 1 - x_crit/x0 is the
    minimum fractional PrP knockdown required for true suppression rather
    than mere slowing. This is the single most decision-relevant number in
    the whole PrP-lowering programme.
    """
    num = p.a * (p.a + p.b * (2 * p.n - 1)) + p.b ** 2 * p.n * (p.n - 1)
    return num / (p.b * p.beta)


def required_knockdown(p: PrionParams):
    """Fractional PrP lowering needed to push r below zero. >1 means
    'impossible even at complete knockout' (should not happen for sane params)."""
    return 1.0 - x_crit(p) / p.x0


# ----------------------------------------------------------------------------
# 2. Full nonlinear system (replication + toxicity)
# ----------------------------------------------------------------------------

def _rhs(t, s, p: PrionParams, knockdown_fn):
    """State vector s = [x, y, z, D, N].

    x  PrP-C concentration
    y  number concentration of PrP-Sc polymers (seeds)
    z  total PrP-Sc monomer mass in polymers
    D  reversible neuronal dysfunction (spongiosis / synaptic deficit)
    N  cumulative irreversible neuron loss
    """
    x, y, z, D, N = s
    x = max(x, 0.0)
    y = max(y, 0.0)
    z = max(z, 0.0)

    # therapy: fractional suppression of PrP-C synthesis at time t
    kd = knockdown_fn(t)
    lam = p.d * p.x0 * (1.0 - kd)

    conv_flux = p.beta * x * y          # PrP-C -> PrP-Sc conversion
    dissolve = p.b * p.n * (p.n - 1) * y  # subcritical fragments returning monomer

    dx = lam - p.d * x - conv_flux + dissolve
    dy = p.b * (z - (2 * p.n - 1) * y) - p.a * y
    dz = conv_flux - p.a * z - dissolve

    # toxicity is driven by ongoing conversion in neurons, NOT by standing load
    # (Mallucci 2003: PrP-Sc kept accumulating extraneuronally while the animal
    #  recovered, once neuronal conversion was switched off)
    dD = p.kappa * conv_flux - p.rho * D
    dN = p.mu * max(0.0, D - p.D_tox)

    return [dx, dy, dz, dD, dN]


def simulate(p: PrionParams, t_end=1200.0, knockdown=0.0, t_treat=0.0,
             ramp_days=7.0, max_step=2.0):
    """Integrate an infected animal, optionally treated.

    knockdown : fractional PrP-C lowering achieved by therapy (0-1)
    t_treat   : day post-inoculation when therapy starts
    ramp_days : time for the drug to reach full effect (ASO/siRNA onset)
    """
    def kd_fn(t):
        if t < t_treat:
            return 0.0
        return knockdown * min(1.0, (t - t_treat) / max(ramp_days, 1e-9))

    s0 = [p.x0, p.y_inoculum, p.y_inoculum * p.n, 0.0, 0.0]

    def terminal(t, s, *args):
        return s[4] - p.N_death
    terminal.terminal = True
    terminal.direction = 1

    sol = solve_ivp(_rhs, (0.0, t_end), s0, args=(p, kd_fn),
                    method="LSODA", rtol=1e-8, atol=1e-12,
                    max_step=max_step, dense_output=True, events=terminal)
    return sol


def survival_time(p: PrionParams, t_end=3000.0, **kw):
    """Day of terminal disease, or np.inf if the animal survives the window."""
    sol = simulate(p, t_end=t_end, **kw)
    if sol.t_events and len(sol.t_events[0]) > 0:
        return float(sol.t_events[0][0])
    return np.inf


def onset_time(sol, p: PrionParams):
    """First day on which clinical signs appear (D crosses D_sx)."""
    D = sol.y[3]
    idx = np.where(D >= p.D_sx)[0]
    return float(sol.t[idx[0]]) if len(idx) else np.inf


# ----------------------------------------------------------------------------
# 3. Analytic survival model used for calibration
# ----------------------------------------------------------------------------

def relative_survival_analytic(x_rel, p: PrionParams, t_tox_frac=0.25):
    """Survival at PrP level x_rel, relative to untreated, from the analytic
    growth rate.

    Total time to death = replication phase + toxic phase (Sandberg 2011).
    The replication phase scales as 1/r(x); the toxic phase is treated as a
    roughly fixed fraction of untreated survival that does not shorten
    indefinitely with slower replication.
    """
    r1 = growth_rate(p.x0, p)
    rx = growth_rate(x_rel * p.x0, p)
    if r1 <= 0:
        return np.nan
    if rx <= 0:
        return np.inf
    return (1.0 - t_tox_frac) * (r1 / rx) + t_tox_frac


def solve_beta_for_target(p: PrionParams, x_rel, target_ratio, t_tox_frac=0.25,
                          bounds=(0.05, 5000.0)):
    """Find the elongation rate beta such that lowering PrP to `x_rel`
    multiplies survival by `target_ratio`.

    This is the key inference step: the observed dose-response at 50% lowering
    constrains where the threshold sits.
    """
    def f(beta):
        q = replace(p, beta=beta)
        if growth_rate(q.x0, q) <= 0:      # must be able to cause disease at all
            return 1e6
        val = relative_survival_analytic(x_rel, q, t_tox_frac)
        if not np.isfinite(val):
            return -1e6
        return val - target_ratio
    lo, hi = bounds
    flo, fhi = f(lo), f(hi)
    if np.sign(flo) == np.sign(fhi):
        return None
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-12)


# ----------------------------------------------------------------------------
# 4. Published calibration anchors
# ----------------------------------------------------------------------------

# (label, residual PrP fraction, survival ratio vs untreated, source)
ANCHORS = [
    ("Prnp-/- knockout",              0.00, np.inf, "Bueler 1993"),
    ("di-siRNA, chronic pre-sx",      0.49, 2.70,   "Gentile 2026 NAR"),
    ("ASO, chronic early",            0.50, 3.00,   "Minikel 2020 (upper)"),
    ("Prnp+/- heterozygote",          0.50, 2.00,   "Bueler 1994"),
    ("base editing, prophylactic",    0.50, 1.52,   "An/Davis 2025 Nat Med"),
    ("di-siRNA, single dose at onset",0.49, 1.64,   "Gentile 2026 NAR"),
    ("wild type",                     1.00, 1.00,   "reference"),
    ("Tga20 overexpressor",           8.00, 0.42,   "Fischer 1996 (approx)"),
]


if __name__ == "__main__":
    p = PrionParams()
    print(f"growth rate at x=1.0 : {growth_rate(1.0, p):.4f} /day")
    print(f"x_crit               : {x_crit(p):.4f}")
    print(f"required knockdown   : {100*required_knockdown(p):.1f}%")
    print(f"untreated survival   : {survival_time(p):.0f} d")
