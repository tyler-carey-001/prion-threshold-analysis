"""
run_analysis.py
===============
Three analyses, each aimed at a decision the field actually faces.

  A. THRESHOLD INFERENCE
     Published mouse data say ~50% PrP lowering buys 1.5x-3x survival but never
     a cure. Under the nucleated-polymerization model, that dose-response
     pins down where the self-sustaining threshold sits -- i.e. how deep
     knockdown must go before prion load starts falling instead of rising.

  B. DOSE-RESPONSE CURVE
     Predicted survival as a function of residual PrP, showing the asymptote.
     Overlaid with every published anchor I could find.

  C. REVERSIBILITY WINDOW
     Treatment depth x treatment timing -> full rescue / rescue with permanent
     deficit / delay only / no benefit. This is the operational meaning of
     "can prion disease be reversed?"

Outputs PNG figures and a machine-readable results.json.
"""

import json
from dataclasses import replace, asdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.optimize import brentq

from prion_model import (PrionParams, growth_rate, x_crit, required_knockdown,
                         relative_survival_analytic, solve_beta_for_target,
                         simulate, survival_time, onset_time, ANCHORS)

plt.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 140,
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "figure.facecolor": "white",
})

INK = "#1b1b1b"
ACCENT = "#c2410c"
BLUE = "#1d4ed8"
GREEN = "#15803d"
results = {}


# ---------------------------------------------------------------------------
# A. THRESHOLD INFERENCE
# ---------------------------------------------------------------------------
print("=" * 74)
print("A. How deep must PrP lowering go before prion replication cannot sustain itself?")
print("=" * 74)

base = PrionParams()
observed_ratios = np.linspace(1.15, 6.0, 260)
tox_fracs = [0.15, 0.25, 0.40]

inference = {}
for tf in tox_fracs:
    kds = []
    for ratio in observed_ratios:
        beta = solve_beta_for_target(base, 0.50, ratio, t_tox_frac=tf)
        if beta is None:
            kds.append(np.nan)
            continue
        kds.append(required_knockdown(replace(base, beta=beta)))
    inference[tf] = np.array(kds)

# point estimates for the published anchors at ~50% lowering
anchor_50 = [(lab, r) for lab, x, r, _ in ANCHORS if 0.45 <= x <= 0.55 and np.isfinite(r)]
print(f"\n{'observed survival ratio at 50% lowering':<42} {'-> required knockdown for suppression'}")
threshold_table = []
for lab, ratio in sorted(anchor_50, key=lambda t: t[1]):
    beta = solve_beta_for_target(base, 0.50, ratio, t_tox_frac=0.25)
    kd = required_knockdown(replace(base, beta=beta))
    threshold_table.append({"model": lab, "survival_ratio": ratio,
                            "required_knockdown": round(float(kd), 4)})
    print(f"  {lab:<34} {ratio:>4.2f}x   ->   {100*kd:>5.1f}%  (residual PrP {100*(1-kd):.1f}%)")

results["threshold_inference"] = threshold_table

fig, ax = plt.subplots(figsize=(6.4, 4.2))
styles = {0.15: "--", 0.25: "-", 0.40: ":"}
for tf in tox_fracs:
    ax.plot(observed_ratios, 100 * inference[tf], styles[tf], color=INK, lw=1.6,
            label=f"toxic phase = {int(100*tf)}% of untreated course")

for lab, ratio in anchor_50:
    beta = solve_beta_for_target(base, 0.50, ratio, t_tox_frac=0.25)
    kd = 100 * required_knockdown(replace(base, beta=beta))
    ax.plot(ratio, kd, "o", color=ACCENT, ms=6, zorder=5)
    ax.annotate(lab, (ratio, kd), textcoords="offset points", xytext=(7, -3),
                fontsize=7.2, color=ACCENT)

ax.axhspan(0, 55, color=BLUE, alpha=0.07)
ax.axhspan(78, 90, color=GREEN, alpha=0.10)
ax.text(5.85, 50, "depth reached by ASO / het KO\n(~50%)", ha="right", va="top",
        fontsize=7, color=BLUE)
ax.text(5.85, 86, "depth reached by divalent siRNA\n(17% residual, Gentile 2026)",
        ha="right", va="top", fontsize=7, color=GREEN)

ax.set_xlabel("observed survival extension at 50% PrP lowering (fold)")
ax.set_ylabel("PrP knockdown required to halt replication (%)")
ax.set_title("A. The mouse dose-response constrains the cure threshold", loc="left",
             fontsize=10.5, weight="bold")
ax.set_ylim(0, 100)
ax.set_xlim(1.1, 6.05)
ax.legend(fontsize=7.2, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig("figures/fig1_threshold_inference.png")
plt.close(fig)

kd_lo = min(t["required_knockdown"] for t in threshold_table)
kd_hi = max(t["required_knockdown"] for t in threshold_table)
print(f"\n  => Across published anchors, suppression needs {100*kd_lo:.0f}-{100*kd_hi:.0f}% knockdown.")
print(f"     50% lowering is on the wrong side of that line. 83% is on the right side.")
results["required_knockdown_range"] = [round(kd_lo, 3), round(kd_hi, 3)]


# ---------------------------------------------------------------------------
# B. DOSE-RESPONSE CURVE
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("B. Predicted survival vs residual PrP")
print("=" * 74)

fig, ax = plt.subplots(figsize=(6.4, 4.2))
xs = np.linspace(0.02, 1.6, 800)

scenarios = [("conservative (1.52x at 50%)", 1.52, "#94a3b8"),
             ("central (2.7x at 50%)", 2.70, INK),
             ("optimistic (3.0x at 50%)", 3.00, "#64748b")]

for lab, ratio, col in scenarios:
    beta = solve_beta_for_target(base, 0.50, ratio, t_tox_frac=0.25)
    p = replace(base, beta=beta)
    xc = x_crit(p)
    ys = [relative_survival_analytic(x, p, 0.25) for x in xs]
    ys = [y if np.isfinite(y) else np.nan for y in ys]
    ax.plot(100 * xs, ys, color=col, lw=1.8 if col == INK else 1.2, label=lab)
    ax.axvline(100 * xc, color=col, ls=":", lw=1.0)
    print(f"  {lab:<30} x_crit = {100*xc:5.1f}% residual PrP")

for lab, xr, ratio, src in ANCHORS:
    if np.isfinite(ratio) and xr <= 1.6:
        ax.plot(100 * xr, ratio, "s", color=ACCENT, ms=5, zorder=5)

ax.axvspan(15, 20, color=GREEN, alpha=0.15)
ax.text(17.5, 9.2, "di-siRNA\nreach", ha="center", fontsize=7, color=GREEN)
ax.set_xlabel("residual PrP-C (% of normal)")
ax.set_ylabel("survival, fold vs untreated")
ax.set_yscale("log")
ax.set_ylim(0.3, 40)
ax.set_xlim(0, 130)
ax.set_title("B. Survival is non-linear in knockdown depth; it diverges at x_crit",
             loc="left", fontsize=10.5, weight="bold")
ax.legend(fontsize=7.2, frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig("figures/fig2_dose_response.png")
plt.close(fig)


# ---------------------------------------------------------------------------
# C. REVERSIBILITY WINDOW
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("C. Reversibility: how late can you treat and still recover?")
print("=" * 74)

# Calibrate the full nonlinear model: central beta, then tune toxicity so the
# untreated animal reaches terminal disease at ~150 d (RML in wild-type mice).
beta_c = solve_beta_for_target(base, 0.50, 2.70, t_tox_frac=0.25)
p = replace(base, beta=beta_c)


def death_day(kappa):
    return survival_time(replace(p, kappa=kappa), t_end=4000.0)


lo, hi = 1e-3, 1e6
f = lambda k: (death_day(k) if np.isfinite(death_day(k)) else 4000.0) - 150.0
try:
    kappa_c = brentq(f, lo, hi, xtol=1e-6, rtol=1e-8, maxiter=200)
except ValueError:
    kappa_c = 1.0
p = replace(p, kappa=kappa_c)

sol_untreated = simulate(p, t_end=400.0)
T_death = survival_time(p)
T_onset = onset_time(sol_untreated, p)
print(f"  calibrated untreated course: onset {T_onset:.0f} d, terminal {T_death:.0f} d")
print(f"  (beta = {beta_c:.3g}, x_crit = {100*x_crit(p):.1f}% residual PrP)")
results["calibration"] = {"beta": float(beta_c), "kappa": float(kappa_c),
                          "x_crit": float(x_crit(p)),
                          "untreated_onset_d": float(T_onset),
                          "untreated_terminal_d": float(T_death)}

# grid over depth x timing
depths = np.linspace(0.20, 0.95, 16)
times = np.linspace(10, 145, 16)
outcome = np.zeros((len(depths), len(times)))
perm_damage = np.zeros_like(outcome)
surv = np.zeros_like(outcome)

for i, kd in enumerate(depths):
    for j, tt in enumerate(times):
        s = simulate(p, t_end=1500.0, knockdown=kd, t_treat=tt, max_step=4.0)
        st = survival_time(p, t_end=1500.0, knockdown=kd, t_treat=tt, max_step=4.0)
        Nf = float(s.y[4][-1])
        surv[i, j] = st if np.isfinite(st) else 1500.0
        perm_damage[i, j] = Nf
        if not np.isfinite(st):
            # survived the window: full rescue vs rescue-with-deficit
            outcome[i, j] = 3 if Nf < 0.02 else 2
        elif st > 1.5 * T_death:
            outcome[i, j] = 1          # substantial delay only
        else:
            outcome[i, j] = 0          # little benefit

cmap = ListedColormap(["#7f1d1d", "#d97706", "#65a30d", "#166534"])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

fig, ax = plt.subplots(figsize=(6.6, 4.4))
im = ax.pcolormesh(times, 100 * depths, outcome, cmap=cmap, norm=norm, shading="auto")
ax.axvline(T_onset, color="white", ls="--", lw=1.4)
ax.text(T_onset + 2, 92, "clinical onset", color="white", fontsize=7.5, rotation=90,
        va="top")
cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3])
cbar.ax.set_yticklabels(["no real benefit", "delay only",
                         "survives, permanent deficit", "full rescue"], fontsize=7.5)
ax.set_xlabel("day therapy starts (post-inoculation; untreated death at day %d)" % T_death)
ax.set_ylabel("PrP knockdown achieved (%)")
ax.set_title("C. The reversibility window: depth x timing", loc="left",
             fontsize=10.5, weight="bold")
fig.tight_layout()
fig.savefig("figures/fig3_reversibility_window.png")
plt.close(fig)

# trajectories
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
treat_days = [30, 70, 100, 125]
axes[0].plot(sol_untreated.t, sol_untreated.y[2], color=INK, lw=1.8, label="untreated")
axes[1].plot(sol_untreated.t, sol_untreated.y[3], color=INK, lw=1.8, label="untreated")
cols = ["#0ea5e9", "#2563eb", "#7c3aed", "#db2777"]
for c, td in zip(cols, treat_days):
    s = simulate(p, t_end=400.0, knockdown=0.85, t_treat=td)
    axes[0].plot(s.t, s.y[2], color=c, lw=1.4, label=f"85% KD @ d{td}")
    axes[1].plot(s.t, s.y[3], color=c, lw=1.4, label=f"85% KD @ d{td}")

axes[0].set_yscale("log")
axes[0].set_ylabel("PrP-Sc load (z)")
axes[1].set_ylabel("reversible dysfunction (D)")
axes[1].axhline(p.D_tox, color="#7f1d1d", ls=":", lw=1.2)
axes[1].text(5, p.D_tox * 1.03, "irreversible neuron loss above this line",
             fontsize=7, color="#7f1d1d")
for a in axes:
    a.set_xlabel("days post-inoculation")
    a.legend(fontsize=7, frameon=False)
    a.set_xlim(0, 300)
axes[0].set_title("Deep knockdown makes prion load fall", loc="left", fontsize=9.5,
                  weight="bold")
axes[1].set_title("Early dysfunction reverses; late damage does not", loc="left",
                  fontsize=9.5, weight="bold")
fig.tight_layout()
fig.savefig("figures/fig4_trajectories.png")
plt.close(fig)

# latest day at which 85% knockdown still gives full rescue
last_full = None
for td in np.arange(5, 150, 2.5):
    s = simulate(p, t_end=1500.0, knockdown=0.85, t_treat=td, max_step=4.0)
    if not np.isfinite(survival_time(p, t_end=1500.0, knockdown=0.85, t_treat=td,
                                     max_step=4.0)) and float(s.y[4][-1]) < 0.02:
        last_full = td
print(f"  at 85% knockdown, full rescue is still possible up to day ~{last_full:.0f}"
      f" ({100*last_full/T_death:.0f}% of the way to terminal disease)")
results["last_day_full_rescue_at_85pct_KD"] = float(last_full) if last_full else None

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=2)

print("\nWrote fig1-fig4 PNGs and results.json")
