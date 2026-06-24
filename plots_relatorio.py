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
mr_calc = [85.8, 80.9, 77.6, 75.8]
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
    for i, (f, c, mb) in enumerate(zip(freqs, calcs, mbs)):
        yoff = 7 if i % 2 == 0 else -14
        ax.annotate(f"$mB={mb}$", xy=(f, c), xytext=(5, yoff),
                    textcoords="offset points", fontsize=8)
    ax.set_title(title)
    ax.set_xlabel("Frequência [Hz]")
    ax.set_ylabel("NPS [dB re 20 uPa]")
    ax.set_ylim(66, 93)
    ax.legend(loc="upper right", framealpha=0.95)

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

mr_rot_vals = {16: 85.8, 31.5: 80.9, 50: 77.6, 63: 75.8}
tr_rot_vals = {100: 85.1, 200: 79.9, 315: 76.6, 400: 74.5}
mr_bb_vals  = {40: 65.0, 80: 69.2, 160: 72.8, 315: 77.0, 630: 73.0, 1250: 72.4, 2500: 68.0}
tr_bb_vals  = {20: 53.4, 40: 57.9, 80: 62.9, 160: 67.1, 315: 70.7, 630: 74.9, 1250: 70.9, 2500: 70.3}

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
ax.legend(loc="lower right", ncol=2, fontsize=8, framealpha=0.95)
fig.tight_layout()
fig.savefig("figures/fig_ex1_combinado.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex1_combinado.png")

# -------------------------------------------------------------
# EXEMPLO I -- Contribuição de cada fonte por banda de 1/3 de oitava
# -------------------------------------------------------------

bands_src = [16, 20, 31.5, 40, 50, 63, 80, 100, 160, 200, 315, 400, 630, 1250, 2500]
band_lbls = ["16", "20", "31,5", "40", "50", "63", "80",
             "100", "160", "200", "315", "400", "630", "1250", "2500"]

src_data_ex1 = [
    ("RP -- Rotacional",  mr_rot_vals, "#1f77b4"),
    ("RC -- Rotacional",  tr_rot_vals, "#ff7f0e"),
    ("RP -- Banda Larga", mr_bb_vals,  "#2ca02c"),
    ("RC -- Banda Larga", tr_bb_vals,  "#9467bd"),
]

n_src1 = len(src_data_ex1)
bw1    = 0.18
x_b1   = np.arange(len(bands_src))
offs1  = np.linspace(-(n_src1 - 1) / 2, (n_src1 - 1) / 2, n_src1) * bw1

fig, ax = plt.subplots(figsize=(14, 5))
labeled1 = set()

for i, (label, data, color) in enumerate(src_data_ex1):
    for j, fc in enumerate(bands_src):
        val = data.get(fc)
        if val is not None:
            lbl = label if label not in labeled1 else "_nolegend_"
            ax.bar(x_b1[j] + offs1[i], val, bw1, color=color, alpha=0.85, label=lbl)
            labeled1.add(label)

# Linha de total energético
tot_x1, tot_y1 = [], []
for j, fc in enumerate(bands_src):
    parts = [d.get(fc) for _, d, _ in src_data_ex1 if d.get(fc) is not None]
    if parts:
        tot_x1.append(x_b1[j])
        tot_y1.append(db_sum(*parts))
ax.plot(tot_x1, tot_y1, "k-o", linewidth=2.0, markersize=5,
        label="Total (energético)", zorder=5)

# Pontos experimentais
exp_items1 = sorted(exp_total.items())
exp_xb1 = [x_b1[bands_src.index(fc)] for fc, _ in exp_items1 if fc in bands_src]
exp_yb1 = [val for fc, val in exp_items1 if fc in bands_src]
ax.plot(exp_xb1, exp_yb1, "r--*", linewidth=1.6, markersize=8,
        label="Experimental (Ref. 13)")

ax.set_xticks(x_b1)
ax.set_xticklabels(band_lbls, rotation=45, ha="right", fontsize=8)
ax.set_xlabel("Frequência central da banda [Hz]")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title(
    "Exemplo I -- Contribuição Individual das Fontes de Ruído\n"
    "por Banda de 1/3 de Oitava  (RP = Rotor Principal, RC = Rotor de Cauda)",
    fontweight="bold")
ax.set_ylim(50, 95)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8, ncol=3, framealpha=0.9)
fig.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig("figures/fig_ex1_fontes.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex1_fontes.png")

# -------------------------------------------------------------
# EXEMPLO II -- Espectro de ruído de interação pá-vórtice (BVI)
# -------------------------------------------------------------

bvi_mB   = [4,    10,   16,   20,   24,    30,    36,    40,    44  ]
bvi_freq = [228,  455,  683,  910,  1138,  1365,  1593,  1820,  2048]
bvi_spl  = [75.8, 87.3, 92.8, 95.8, 97.1,  96.5,  94.3,  89.3,  76.5]

# Alternância de offset das anotações: índices pares acima, ímpares abaixo
_bvi_yoff = [8, 8, 8, -16, 8, -16, 8, -16, 8]

fig, ax = plt.subplots(figsize=(9, 4.5))
markerline, stemlines, baseline = ax.stem(bvi_freq, bvi_spl, linefmt="C0-",
                                           markerfmt="C0o", basefmt="k-")
plt.setp(stemlines, linewidth=1.4)
plt.setp(markerline, markersize=7)
markerline.set_label("Harmônicos calculados")

# Envoltória
envelope, = ax.plot(bvi_freq, bvi_spl, "r--", linewidth=1.2, alpha=0.7, label="Envoltória")

for i, (f, s, mb) in enumerate(zip(bvi_freq, bvi_spl, bvi_mB)):
    ax.annotate(f"$mB={mb}$", xy=(f, s), xytext=(0, _bvi_yoff[i]),
                textcoords="offset points", ha="center", fontsize=8)

ax.set_xlabel("Frequência [Hz]")
ax.set_ylabel("NPS [dB re 20 uPa]")
ax.set_title("Exemplo II -- Espectro Harmônico de Interação Pá-Vórtice (BVI)", fontweight="bold")
ax.legend(handles=[markerline, envelope], labels=["Harmônicos calculados", "Envoltória"],
          loc="upper right", framealpha=0.95)
ax.set_ylim(65, 108)
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
spl_T_09 = [98.5, 104.1, 108.5, 110.2, 110.9, 111.1, 110.1]
spl_R_09 = [94.7,  95.9,  96.5,  95.8,  94.7,  93.7,  91.6]
spl_tot09 = [102.6, 106.1, 109.4, 110.8, 111.3, 111.4, 110.3]

# Me=0.8 (T=676 N, N=41.95 rps, f0=83.9 Hz; compressibilidade inativa)
freq_08   = [84,  168, 336, 503, 671, 839, 1258]
spl_T_08  = [94.7, 98.9, 100.6, 99.7, 98.1, 95.9, 89.2]
spl_R_08  = [96.4, 96.3,  95.4,  92.7, 89.8, 88.0, 84.4]
spl_tot08 = [98.6, 100.8, 101.8, 100.5, 98.7, 96.6, 90.4]

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
ax.legend(loc="upper left", fontsize=8)

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
ax.legend(loc="upper right", fontsize=8)
ax.text(0.02, 0.97, "$M_e = M_{dd}$ → Compressibilidade inativa",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

fig.suptitle("Exemplo III -- Contribuição Individual das Fontes de Ruído\nModelo 1/7 do Rotor UH-1 em Pairado",
             fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig("figures/fig_ex3_componentes.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex3_componentes.png")

# -------------------------------------------------------------
# EXEMPLO III -- Espectro calculado vs. referência Pegg (1979)
# -------------------------------------------------------------

ref_09 = [102.6, 106.1, 109.4, 110.7, 111.3, 111.4, 110.2]
ref_08 = [97.2,  100.5, 101.2, 100.1,  98.4,  96.3,  89.9]

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

# Painel esquerdo: M_e = 0.9 — captura handles para legenda compartilhada
ax = axes[0]
h_C,   = ax.plot(mB_vals, spl_C_09,  "^--", color="#e377c2", linewidth=1.4, markersize=6)
h_T,   = ax.plot(mB_vals, spl_T_09,  "s--", color="#1f77b4", linewidth=1.4, markersize=6)
h_R,   = ax.plot(mB_vals, spl_R_09,  "v--", color="#ff7f0e", linewidth=1.4, markersize=6)
h_tot, = ax.plot(mB_vals, spl_tot09, "k-o", linewidth=2.0,   markersize=7)
h_ref, = ax.plot(mB_vals, ref_09,    "r*",  markersize=11,   zorder=5)
ax.set_xlabel("Harmônico $mB$")
ax.set_ylabel("NPS [dB re 20 $\\mu$Pa]")
ax.set_title("$M_e = 0{,}9$  —  $T = 485$ N", fontweight="bold")
ax.set_ylim(85, 120)
ax.set_xticks(mB_vals)

# Painel direito: M_e = 0.8
ax = axes[1]
ax.plot(mB_vals, spl_T_08,  "s--", color="#1f77b4", linewidth=1.4, markersize=6)
ax.plot(mB_vals, spl_R_08,  "v--", color="#ff7f0e", linewidth=1.4, markersize=6)
ax.plot(mB_vals, spl_tot08, "k-o", linewidth=2.0,   markersize=7)
ax.plot(mB_vals, ref_08,    "r*",  markersize=11,   zorder=5)
ax.text(0.03, 0.97,
        "Compressibilidade inativa ($M_e = M_{dd}$)\n"
        "$T_{\\mathrm{cód.}} = 676$ N  ·  Ref.: $T = 485$ N",
        transform=ax.transAxes, fontsize=7.5, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.85))
ax.set_xlabel("Harmônico $mB$")
ax.set_ylabel("NPS [dB re 20 $\\mu$Pa]")
ax.set_title("$M_e = 0{,}8 = M_{dd}$  —  $T_{\\mathrm{cód.}} = 676$ N", fontweight="bold")
ax.set_ylim(80, 115)
ax.set_xticks(mB_vals)

# Legenda compartilhada abaixo dos dois painéis
_leg_labels = ["Compressibilidade", "Espessura", "Rotacional",
               "Total calculado", "Referência Pegg (1979)"]
fig.legend([h_C, h_T, h_R, h_tot, h_ref], _leg_labels,
           loc="upper center", bbox_to_anchor=(0.5, 0.01), ncol=5,
           fontsize=8, framealpha=0.95)

fig.suptitle(
    "Exemplo III -- NPS Total e Fontes Individuais vs. Referência Pegg (1979)\n"
    "Modelo 1/7 do Rotor UH-1 em Pairado",
    fontsize=11, fontweight="bold",
)
fig.tight_layout(rect=[0, 0.10, 1, 1])
fig.savefig("figures/fig_ex3_espectro.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_ex3_espectro.png")

# =============================================================
# PROPOSTAS DE ENGENHARIA — Gráficos paramétricos (§5.1–§5.3)
# =============================================================

# --- §5.1 Sensibilidade ao Mach de ponta (espessura + banda larga) ---

M_t_sweep = np.linspace(0.60, 0.95, 300)
M_t_ref_  = 0.9
delta_T_  = 40.0 * np.log10(M_t_sweep / M_t_ref_)
delta_BB_ = 60.0 * np.log10(M_t_sweep / M_t_ref_)

M_t_81_  = 0.81
dT_81_   = 40.0 * math.log10(M_t_81_ / M_t_ref_)
dBB_81_  = 60.0 * math.log10(M_t_81_ / M_t_ref_)

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(M_t_sweep, delta_T_,  "-",  color="#1f77b4", linewidth=2.0,
        label="Espessura  ($\\propto 40\\log M_t$)  [Eq. (5)]")
ax.plot(M_t_sweep, delta_BB_, "--", color="#2ca02c", linewidth=2.0,
        label="Banda Larga  ($\\propto 60\\log M_t$)  [Eq. (8)]")
ax.axvline(M_t_ref_, color="k",       linewidth=0.9, linestyle=":",  alpha=0.7)
ax.axvline(M_t_81_,  color="gray",    linewidth=1.0, linestyle="--", alpha=0.6)
ax.axhline(0.0,      color="k",       linewidth=0.8, alpha=0.4)
ax.axvline(0.80,     color="#d62728", linewidth=1.2, linestyle="-.", alpha=0.85)

# Marcadores nos pontos de interesse (M_t = 0.81)
ax.scatter([M_t_81_], [dT_81_],  s=55, color="#1f77b4", zorder=6)
ax.scatter([M_t_81_], [dBB_81_], s=55, color="#2ca02c", marker="^", zorder=6)

# Rótulos das linhas verticais (coordenadas de dados, sem transform misto)
ax.text(0.903, -13.5,  "Ref.\n$M_t{=}0{,}9$",
        fontsize=7.5, color="k", va="center", ha="left")
ax.text(0.77, -13.5, "Limiar\ncompressib.\n($M_{dd}{=}0{,}8$)",
        fontsize=7.5, color="#d62728", va="center", ha="left")
ax.text(0.813, -13.5, "$M_t{=}0{,}81$\n($-10\\%$)",
        fontsize=7.5, color="gray", va="bottom", ha="left")

# Caixa de valores-resumo no canto inferior esquerdo (abaixo das curvas)
ax.text(0.615, -14.2,
        f"$\\Delta$NPS em $M_t = 0{{,}}81$ ($-10\\%$):\n"
        f"  • Espessura: ${dT_81_:.1f}$ dB\n"
        f"  • Banda larga: ${dBB_81_:.1f}$ dB",
        fontsize=8.5, va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.92,
                  edgecolor="gray", linewidth=0.8))

ax.set_xlabel("Mach de ponta $M_t$")
ax.set_ylabel("$\\Delta$NPS relativo a $M_t = 0{,}9$ [dB]")
ax.set_title("Proposta §5.1 — Sensibilidade do NPS à Redução de $M_t$\n"
             "(componente dominante de cada fonte)", fontweight="bold")
ax.set_xlim(0.60, 0.95)
ax.set_ylim(-15, 3)
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig("figures/fig_eng_velocidade.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_eng_velocidade.png")

# --- §5.2 NPS total de compressibilidade vs M_dd ---

_ME09   = 0.9
_RE09   = 0.945
_r09    = 3.14
_c09    = 0.0762
_R09    = 1.05
_dC     = {2: 173.9, 4: 175.3, 8: 176.0, 12: 175.7, 16: 174.9, 20: 173.9, 30: 170.2}
_mBC    = sorted(_dC.keys())
_baseC  = 20.0 * math.log10(_RE09 / _r09) - 21.6


def _total_compress(mdd):
    exc = _ME09 - mdd
    if exc <= 0.0:
        return float("nan")
    fac = 20.0 * math.log10(exc * _c09 / _R09)
    return 10.0 * math.log10(sum(10.0 ** ((_baseC + fac + _dC[mb]) / 10.0) for mb in _mBC))


M_dd_sweep = np.linspace(0.72, 0.895, 400)
C_curve    = [_total_compress(v) for v in M_dd_sweep]

val_C_080 = _total_compress(0.80)
val_C_085 = _total_compress(0.85)

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(M_dd_sweep, C_curve, "-", color="#e377c2", linewidth=2.2,
        label="NPS total de compressibilidade\n(soma energética, 7 harmônicos — Eq. (4))")

ax.axvline(0.80, color="#d62728", linestyle="--", linewidth=1.3, alpha=0.8)
ax.plot(0.80, val_C_080, "ro", markersize=9, zorder=5)
ax.annotate(f"$M_{{dd}}=0{{,}}80$\n(Exemplo III)\n{val_C_080:.1f} dB",
            xy=(0.80, val_C_080), xytext=(0.815, val_C_080 - 6),
            fontsize=8, color="#d62728", va="top",
            arrowprops=dict(arrowstyle="->", color="#d62728",
                            connectionstyle="arc3,rad=-0.15"))

ax.axvline(0.85, color="#2ca02c", linestyle="--", linewidth=1.3, alpha=0.8)
ax.plot(0.85, val_C_085, "g^", markersize=9, zorder=5)
ax.annotate(f"$M_{{dd}}=0{{,}}85$\n(perfil avançado)\n{val_C_085:.1f} dB"
            f"  ({val_C_085 - val_C_080:+.1f} dB)",
            xy=(0.85, val_C_085), xytext=(0.886, val_C_085+8),
            fontsize=8, color="#2ca02c", va="top", ha="right",
            arrowprops=dict(arrowstyle="->", color="#2ca02c",
                            connectionstyle="arc3,rad=0.2"))

ax.set_xlabel("Mach de divergência de arrasto $M_{dd}$")
ax.set_ylabel("NPS total de compressibilidade [dB re 20 $\\mu$Pa]")
ax.set_title("Proposta §5.2 — Redução do Ruído de Compressibilidade por Aumento de $M_{dd}$\n"
             "(Exemplo III: $M_e = 0{,}9$, $r = 3{,}14$ m, 7 harmônicos)",
             fontweight="bold")
ax.legend(fontsize=9, loc="lower left")
ax.set_xlim(0.72, 0.90)
fig.tight_layout()
fig.savefig("figures/fig_eng_mdd.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_eng_mdd.png")

# --- §5.3 NPS de BVI vs ΔL/L₀ ---

_T_BVI   = 458.35
_N_BVI   = 22.75
_RHO_BVI = 1.228
_C0_BVI  = 344.0
_r_BVI   = 1.02
_BETA_BVI = 49.0
_dBV     = {4: -23.0, 10: -11.5, 16: -6.0, 20: -3.0, 24: -1.7,
            30: -2.3,  36: -4.5,  40: -9.5,  44: -22.3}
_mBBVI   = sorted(_dBV.keys())

_baseTN  = 20.0 * math.log10(_T_BVI * _N_BVI / (_RHO_BVI * _C0_BVI**3 * _r_BVI))
_cosBeta = math.cos(math.radians(_BETA_BVI))


def _calc_bvi(dl, mb):
    return 20.0 * math.log10(dl * _cosBeta) + _baseTN + _dBV[mb] + 190.2


def _total_bvi(dl):
    return 10.0 * math.log10(sum(10.0 ** (_calc_bvi(dl, mb) / 10.0) for mb in _mBBVI))


DL_sweep     = np.linspace(0.04, 0.35, 300)
spl_pk24     = [_calc_bvi(dl, 24) for dl in DL_sweep]
spl_tot_bvi_ = [_total_bvi(dl) for dl in DL_sweep]

val_pk_020 = _calc_bvi(0.20, 24)
val_pk_010 = _calc_bvi(0.10, 24)

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(DL_sweep, spl_pk24,     "-",  color="#9467bd", linewidth=2.0,
        label="NPS pico ($mB=24$, $f\\approx1138$ Hz)  [Eq. (6)]")
ax.plot(DL_sweep, spl_tot_bvi_, "--", color="#8c564b", linewidth=2.0,
        label="NPS total energético (9 harmônicos)")

ax.axvline(0.20, color="#d62728", linestyle="--", linewidth=1.3, alpha=0.8)
ax.plot(0.20, val_pk_020, "ro", markersize=9, zorder=5)
ax.annotate(f"$\\Delta L/L_0=0{{,}}20$\n(Exemplo II)\n{val_pk_020:.1f} dB",
            xy=(0.20, val_pk_020), xytext=(0.25, 95),
            fontsize=8, color="#d62728", va="center",
            arrowprops=dict(arrowstyle="->", color="#d62728",
                            connectionstyle="arc3,rad=0.2"))

ax.axvline(0.10, color="#2ca02c", linestyle="--", linewidth=1.3, alpha=0.8)
ax.plot(0.10, val_pk_010, "g^", markersize=9, zorder=5)
ax.annotate(f"$\\Delta L/L_0=0{{,}}10$\n(alvo IBC/HHC)\n{val_pk_010:.1f} dB"
            f"\n({val_pk_010 - val_pk_020:+.1f} dB)",
            xy=(0.10, val_pk_010), xytext=(0.15, 90.0),
            fontsize=8, color="#2ca02c", va="top",
            arrowprops=dict(arrowstyle="->", color="#2ca02c",
                            connectionstyle="arc3,rad=0.3"))

ax.set_xlabel("Razão de carga impulsiva $\\Delta L / L_0$")
ax.set_ylabel("NPS de BVI [dB re 20 $\\mu$Pa]")
ax.set_title("Proposta §5.3 — Redução do Ruído de BVI pelo Controle de $\\Delta L / L_0$\n"
             "(Parâmetros do Exemplo II: $B=2$, $M_t=0{,}433$, $\\beta=49°$)",
             fontweight="bold")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figures/fig_eng_bvi.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvo: figures/fig_eng_bvi.png")

print("\nTodas as figuras geradas com sucesso.")
