"""
plots_relatorio.py
==================
Gera as figuras do Trabalho Final de Aeroacústica (PME3558).
Metodologia de Pegg (1979) -- NASA TM-80200.

As figuras são salvas em figures/ (PNG, 150 dpi).
Resultados numéricos extraídos das saídas verificadas de
example_I.py, example_II.py e example_III.py.
"""

import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Garantir que o diretório de figuras existe
os.makedirs("figures", exist_ok=True)

# Configuração global de estilo
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.35,
    "grid.linestyle": "--",
})

# -------------------------------------------------------------
# EXEMPLO I -- Ruído rotacional: harmônicos calculados vs. experimentais
# -------------------------------------------------------------

# Rotor principal (Me=0.54, beta=80deg)
mr_mB   = [5,   10,   15,   20  ]
mr_freq = [17.5, 35.0, 52.5, 70.0]
mr_calc = [85.7, 80.8, 77.5, 75.7]
mr_exp  = [85.0, 79.5, 76.0, 73.0]

# Rotor de cauda (Me=0.54, beta=80deg)
tr_mB   = [5,     10,    15,    20   ]
tr_freq = [105.5, 211.0, 316.5, 422.0]
tr_calc = [85.1,  79.9,  76.6,  74.5 ]
tr_exp  = [82.0,  76.0,  75.0,  73.0 ]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

for ax, freqs, calcs, exps, mbs, title in [
    (axes[0], mr_freq, mr_calc, mr_exp, mr_mB, "Rotor Principal"),
    (axes[1], tr_freq, tr_calc, tr_exp, tr_mB, "Rotor de Cauda"),
]:
    ax.plot(freqs, calcs, "o-", color="#1f77b4", linewidth=1.6,
            markersize=6, label="Calculado")
    ax.plot(freqs, exps, "s--", color="#d62728", linewidth=1.4,
            markersize=6, label="Experimental")
    for f, c, mb in zip(freqs, calcs, mbs):
        ax.annotate(f"$mB={mb}$", xy=(f, c), xytext=(4, 4),
                    textcoords="offset points", fontsize=8)
    ax.set_title(title)
    ax.set_xlabel("Frequência [Hz]")
    ax.set_ylabel("NPS [dB re 20 uPa]")
    ax.set_ylim(68, 92)
    ax.legend(loc="upper right")

fig.suptitle("Exemplo I -- Ruído Rotacional (Harmônicos)", fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig("figures/fig_ex1_rotacional.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex1_rotacional.png")

# -------------------------------------------------------------
# EXEMPLO I -- Espectro combinado em bandas de 1/3 de oitava
# -------------------------------------------------------------

# Contribuições por banda (NPS em dB; None = ausente)
bands_hz = [16, 20, 31.5, 40, 50, 63, 80, 100, 160, 200, 315, 400, 630, 1250, 2500]

mr_rot_vals = {16: 85.7, 31.5: 80.8, 50: 77.5, 63: 75.7}
tr_rot_vals = {100: 85.1, 200: 79.9, 315: 76.6, 400: 74.5}
mr_bb_vals  = {40: 64.8, 80: 69.0, 160: 72.6, 315: 76.8, 630: 72.8, 1250: 72.2, 2500: 67.8}
tr_bb_vals  = {20: 56.1, 40: 58.6, 80: 63.6, 160: 67.8, 315: 71.4, 630: 75.6, 1250: 71.6, 2500: 71.0}

exp_total   = {16: 85.0, 31.5: 79.5, 50: 76.0, 80: 73.0,
               100: 82.0, 200: 76.0, 315: 75.0, 400: 73.0}

def db_sum(*vals):
    v = [x for x in vals if x is not None]
    if not v:
        return None
    return 10.0 * math.log10(sum(10.0**(x/10.0) for x in v))

totals = {}
for fc in bands_hz:
    parts = [mr_rot_vals.get(fc), tr_rot_vals.get(fc),
             mr_bb_vals.get(fc),  tr_bb_vals.get(fc)]
    totals[fc] = db_sum(*parts)

# Separar bandas com contribuição
def get_series(d):
    xs = sorted(d.keys())
    return xs, [d[x] for x in xs]

fig, ax = plt.subplots(figsize=(11, 5))

def plot_contrib(d, marker, color, label, linestyle="-"):
    xs, ys = get_series(d)
    ax.plot(xs, ys, marker=marker, linestyle=linestyle, color=color,
            markersize=5, linewidth=1.4, label=label)

plot_contrib(mr_rot_vals, "o",  "#1f77b4", "RP -- Rotacional")
plot_contrib(tr_rot_vals, "s",  "#ff7f0e", "RC -- Rotacional")
plot_contrib(mr_bb_vals,  "^",  "#2ca02c", "RP -- Banda Larga")
plot_contrib(tr_bb_vals,  "v",  "#9467bd", "RC -- Banda Larga")

# Total calculado
tot_xs = [fc for fc in bands_hz if totals[fc] is not None]
tot_ys = [totals[fc] for fc in tot_xs]
ax.plot(tot_xs, tot_ys, "k-", linewidth=2.0, marker="D", markersize=5, label="Total Calculado")

# Experimental
exp_xs, exp_ys = get_series(exp_total)
ax.plot(exp_xs, exp_ys, "r--", linewidth=1.8, marker="*", markersize=8, label="Experimental (Ref. 13)")

ax.set_xscale("log")
ax.set_xlabel("Frequência central da banda [Hz]")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title("Exemplo I -- Espectro Combinado em Bandas de 1/3 de Oitava", fontweight="bold")
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
ax.set_xticks([16, 20, 31.5, 40, 50, 63, 80, 100, 160, 200, 315, 400, 630, 1250, 2500])
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_ylim(50, 95)
ax.legend(loc="upper right", ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig("figures/fig_ex1_combinado.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex1_combinado.png")

# -------------------------------------------------------------
# EXEMPLO II -- Espectro de ruído de interação pá-vórtice (BVI)
# -------------------------------------------------------------

bvi_mB   = [4,    10,   16,   20,   24,    30,    36,    40,    44  ]
bvi_freq = [228,  455,  683,  910,  1138,  1365,  1593,  1820,  2048]
bvi_spl  = [75.8, 87.3, 92.8, 95.8, 97.1,  96.5,  94.3,  89.3,  76.5]

fig, ax = plt.subplots(figsize=(9, 4.5))
markerline, stemlines, baseline = ax.stem(bvi_freq, bvi_spl, linefmt="C0-",
                                           markerfmt="C0o", basefmt="k-")
plt.setp(stemlines, linewidth=1.4)
plt.setp(markerline, markersize=7)

# Envoltória
ax.plot(bvi_freq, bvi_spl, "r--", linewidth=1.2, alpha=0.7, label="Envoltória")

for f, s, mb in zip(bvi_freq, bvi_spl, bvi_mB):
    ax.annotate(f"$mB={mb}$", xy=(f, s), xytext=(0, 6),
                textcoords="offset points", ha="center", fontsize=8)

ax.set_xlabel("Frequência [Hz]")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title("Exemplo II -- Espectro Harmônico de Interação Pá-Vórtice (BVI)", fontweight="bold")
ax.legend(["Harmônicos calculados", "Envoltória"])
ax.set_ylim(65, 105)
ax.set_xlim(0, 2300)
fig.tight_layout()
fig.savefig("figures/fig_ex2_bvi.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex2_bvi.png")

# -------------------------------------------------------------
# EXEMPLO III -- Contribuição de cada fonte de ruído
#   Me=0.9 (esquerda) e Me=0.8 (direita)
# -------------------------------------------------------------

mB_vals  = [2,  4,   8,   12,  16,  20,  30 ]

# Me=0.9 (T=485 N, N=47.2 rps, f0=94.4 Hz)
freq_09  = [94,  188, 376, 564, 752, 940, 1410]
spl_C_09 = [99.1, 100.5, 101.2, 100.9, 100.1, 99.1,  95.4]
spl_T_09 = [98.4, 104.1, 108.4, 110.1, 110.8, 111.0, 110.0]
spl_R_09 = [94.7,  95.9,  96.5,  95.8,  94.7,  93.7,  91.6]
spl_tot09 = [102.6, 106.1, 109.4, 110.7, 111.3, 111.4, 110.2]

# Me=0.8 (T=676 N, N=41.95 rps, f0=83.9 Hz; compressibilidade inativa)
freq_08   = [84,  168, 336, 503, 671, 839, 1258]
spl_T_08  = [94.7, 99.6, 100.6, 99.7, 98.1, 95.9, 89.2]
spl_R_08  = [93.5, 93.4,  92.5,  89.8, 86.9, 85.1, 81.5]
spl_tot08 = [97.2, 100.5, 101.2, 100.1, 98.4, 96.3, 89.9]

x = np.arange(len(mB_vals))
width = 0.22

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

# Painel esquerdo: Me=0.9
ax = axes[0]
ax.bar(x - width, spl_C_09, width, label="Compressibilidade", color="#e377c2", alpha=0.85)
ax.bar(x,         spl_T_09, width, label="Espessura",          color="#1f77b4", alpha=0.85)
ax.bar(x + width, spl_R_09, width, label="Rotacional",         color="#ff7f0e", alpha=0.85)
ax.plot(x, spl_tot09, "k-o", linewidth=1.8, markersize=5, label="Total (energético)", zorder=5)
ax.set_xticks(x)
ax.set_xticklabels([f"{mb}\n({f} Hz)" for mb, f in zip(mB_vals, freq_09)], fontsize=8)
ax.set_xlabel("Harmônico mB (Frequência [Hz])")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title("$M_e = 0{,}9$  --  $T = 485$ N", fontweight="bold")
ax.set_ylim(80, 120)
ax.legend(loc="lower left", fontsize=8)

# Painel direito: Me=0.8
ax = axes[1]
ax.bar(x - width/2, spl_T_08, width, label="Espessura",  color="#1f77b4", alpha=0.85)
ax.bar(x + width/2, spl_R_08, width, label="Rotacional", color="#ff7f0e", alpha=0.85)
ax.plot(x, spl_tot08, "k-o", linewidth=1.8, markersize=5, label="Total (energético)", zorder=5)
ax.set_xticks(x)
ax.set_xticklabels([f"{mb}\n({f} Hz)" for mb, f in zip(mB_vals, freq_08)], fontsize=8)
ax.set_xlabel("Harmônico mB (Frequência [Hz])")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title("$M_e = 0{,}8 = M_{dd}$  --  $T = 676$ N  (compressibilidade inativa)", fontweight="bold")
ax.set_ylim(75, 115)
ax.legend(loc="lower left", fontsize=8)
ax.text(0.02, 0.97, "Compressibilidade: $M_e = M_{dd}$ -> inativa",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

fig.suptitle("Exemplo III -- Contribuição Individual das Fontes de Ruído\nModelo 1/7 do Rotor UH-1 em Pairado",
             fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig("figures/fig_ex3_componentes.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex3_componentes.png")

print("\nTodas as figuras geradas com sucesso.")
