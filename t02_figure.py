"""fig_t02_penetrance_shift.png — 2016 vs gnomAD v4.1.1, per variant."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t02_power_check import penetrance_confint
from t02_power_check_pervariant import N_CASE, VARIANTS

N_EXAC = 60706
V4 = {"P102L": (2, 1461704), "A117V": (0, 1461738),
      "D178N": (1, 1461876), "E200K": (13, 1461878)}

fig, ax = plt.subplots(figsize=(9.5, 5.0))
ys = range(len(VARIANTS))[::-1]
for y, (name, ac_case) in zip(ys, VARIANTS):
    lo16, b16, hi16 = penetrance_confint(ac_case, N_CASE, 0, N_EXAC)
    ac4, an4 = V4[name]
    lo4, b4, hi4 = penetrance_confint(ac_case, N_CASE, ac4, an4 // 2)
    ax.plot([lo16 * 100, hi16 * 100], [y + .16] * 2, lw=5, color="#9aa5b1",
            solid_capstyle="butt")
    col = "#B4436C" if name == "A117V" else "#1A6966"
    ax.plot([lo4 * 100, hi4 * 100], [y - .16] * 2, lw=5, color=col,
            solid_capstyle="butt")
    # A clamped point estimate is a censored bound, not a measurement. Draw it
    # as a right-pointing open marker so an external reader cannot misread it
    # as an estimated value sitting at 100%.
    if b4 >= 1.0:
        ax.plot(b4 * 100, y - .16, ">", ms=9, mfc="white", mec=col, mew=2, zorder=3)
    else:
        ax.plot(b4 * 100, y - .16, "o", ms=8, color=col, zorder=3)
    if b16 >= 1.0:
        ax.plot(b16 * 100, y + .16, ">", ms=9, mfc="white", mec="#5b6670",
                mew=2, zorder=4)
    tag = f"AC={ac4}" + ("  (zero persists)" if ac4 == 0 else "")
    ax.text(101.5, y - .16, tag, va="center", fontsize=8.5, color=col)

ax.axvspan(60, 90, color="#f0c419", alpha=.16, zorder=0)
ax.set_ylim(-0.95, len(VARIANTS) - 0.45)
ax.text(75, -0.88, "published E200K survival estimates (~60–90%)",
        ha="center", va="bottom", fontsize=8, color="#8a6d00")
ax.set_yticks(list(ys))
ax.set_yticklabels([f"{n}\n{c} case alleles" for n, c in VARIANTS], fontsize=9)
ax.set_xlabel("Lifetime risk in heterozygotes (%)  —  95% CI, corner method")
ax.set_xlim(0, 118)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_title("PRNP penetrance: 2016 ExAC (grey, n=60,706) vs gnomAD v4.1.1 "
             "(teal, n≈730,939)\ncrude estimates; case counts frozen at 2016\n▷ = censored bound (clamped at 100%), not a point estimate",
             fontsize=10.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(axis="x", alpha=.25, lw=.6)
fig.tight_layout()
fig.savefig("fig_t02_penetrance_shift.png", dpi=170)
print("wrote fig_t02_penetrance_shift.png")
