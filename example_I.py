"""
example_I.py - Exemplo I de Pegg (1979) NASA TM-80200.
Helicóptero de grande porte em pairado: ruído rotacional e de banda larga
do rotor principal e do rotor de cauda. Equações (2), (3), (7), (8).

ΔSPL_mr e ΔSPL_tr: leituras das Figs. 2 e 3 em M_e = 0.54, β = 80°.
S_1/3: interpolação linear em log(f/f_p), ancorada em leituras da Fig. 9
       (f_p,RP = 289 Hz; f_p,RC = 544 Hz).

Nota: ambos os rotores usam V_T = 2π·N·R na Eq. (8). O rotor principal
resulta em base ~0.2 dB acima da tabela do relatório, que parece ter usado
V_T = M_tip·c_0 = 206.4 m/s. Para o rotor de cauda persiste uma discrepância
de ~0.7 dB sem origem identificável nos dados do relatório.
"""

import math
import os
import sys
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pegg1979 import (
    RotorCondition,
    rotational_noise,
    broadband_peak_frequency,
    broadband_spl,
)

# Valores de ΔSPL — ruído rotacional, M_e = 0.54, β = 80°
_DSPL_MR = {5: 148.4, 10: 143.5, 15: 140.2, 20: 138.4}  # Fig. 2 — rotor principal
_DSPL_TR = {5: 154.4, 10: 149.2, 15: 145.9, 20: 143.8}  # Fig. 3 — rotor de cauda


def _dspl_mr(mB: int):
    return lambda M_e, beta_deg: _DSPL_MR[mB]


def _dspl_tr(mB: int):
    return lambda M_e, beta_deg: _DSPL_TR[mB]


# Espectro S_1/3 (Fig. 9 de Schlegel) — 19 pontos (f/f_p, S)
_FP_MR_NOM = 289
_FP_TR_NOM = 544

_FIG9_PTS = sorted([
    (20   / _FP_TR_NOM, -29.0),
    (40   / _FP_TR_NOM, -24.5),
    (40   / _FP_MR_NOM, -19.5),
    (80   / _FP_TR_NOM, -19.5),
    (80   / _FP_MR_NOM, -15.3),
    (160  / _FP_TR_NOM, -15.3),
    (160  / _FP_MR_NOM, -11.7),
    (315  / _FP_TR_NOM, -11.7),
    (315  / _FP_MR_NOM,  -7.5),
    (630  / _FP_TR_NOM,  -7.5),
    (630  / _FP_MR_NOM, -11.5),
    (1250 / _FP_TR_NOM, -11.5),
    (1250 / _FP_MR_NOM, -12.1),
    (2500 / _FP_TR_NOM, -12.1),
    (2500 / _FP_MR_NOM, -16.5),
    (5000 / _FP_TR_NOM, -16.5),
    (5000 / _FP_MR_NOM, -17.0),
    (10000/ _FP_TR_NOM, -17.0),
    (10000/ _FP_MR_NOM, -21.8),
])

_LOG_PTS = [(math.log10(r), s) for r, s in _FIG9_PTS]


def _s13(f_over_fp: float) -> float:
    """Interpolação linear por partes de S_1/3 no espaço log(f/f_p)."""
    lx = math.log10(f_over_fp)
    if lx <= _LOG_PTS[0][0]:
        return _LOG_PTS[0][1]
    if lx >= _LOG_PTS[-1][0]:
        return _LOG_PTS[-1][1]
    for i in range(len(_LOG_PTS) - 1):
        x0, y0 = _LOG_PTS[i]
        x1, y1 = _LOG_PTS[i + 1]
        if x0 <= lx <= x1:
            t = (lx - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return _LOG_PTS[-1][1]


# Parâmetros do rotor e do observador
C_0 = 344.0
RHO = 1.228
M_TIP = 0.6

# Rotor principal
B_MR  = 5
R_MR  = 9.45
N_MR  = 3.5
T_MR  = 69420.0
A_B_MR = 18.58
CL_MR  = 0.438
r_MR   = 61.6
BETA_MR   = 5.0
THETA1_MR = 90.0 - BETA_MR
R_E_MR = 0.9 * R_MR

# Rotor de cauda
B_TR  = 5
R_TR  = 1.52
N_TR  = 21.1
T_TR  = 5206.0
A_B_TR = 3.44
CL_TR  = 0.182
r_TR   = 62.2
BETA_TR   = 80.0
THETA1_TR = 90.0 - BETA_TR
R_E_TR = 0.9 * R_TR

# M_e em 0.9R (pairado: M_e = M_tip × 0.9)
M_E = M_TIP * 0.9   # 0.54

# Velocidades de ponta para banda larga (Eqs. 7 e 8)
V_T_MR = 2.0 * math.pi * N_MR * R_MR
V_T_TR = 2.0 * math.pi * N_TR * R_TR


# Para ruído rotacional, M_t = M_E = 0.54 para que effective_mach() retorne 0.54
# em pairado (M_F = 0), satisfazendo a faixa [0.5, 0.9] e a entrada das Figs. 2 e 3.
# Para banda larga, M_t = V_T / c_0 para que cond.V_T seja fisicamente correto na Eq. (8).
def _make_cond(B, R, R_e, N, T, r, beta, theta1, A_b, CL, M_t=M_E):
    return RotorCondition(
        B=B, R=R, R_e=R_e,
        c=A_b / (B * R),
        h=0.0,
        M_t=M_t, M_F=0.0, M_dd=0.82, N=N,
        T=T, r=r, beta=beta, psi=90.0,
        rho=RHO, c_0=C_0,
        A_b=A_b, CL_bar=CL,
        L_0=0.0, DeltaL=0.0,
        theta_1=theta1,
    )


COND_MR_ROT = _make_cond(B_MR, R_MR, R_E_MR, N_MR, T_MR, r_MR, BETA_MR, THETA1_MR, A_B_MR, CL_MR)
COND_TR_ROT = _make_cond(B_TR, R_TR, R_E_TR, N_TR, T_TR, r_TR, BETA_TR, THETA1_TR, A_B_TR, CL_TR)

COND_MR_BB = replace(COND_MR_ROT, M_t=V_T_MR / C_0)
COND_TR_BB = replace(COND_TR_ROT, M_t=V_T_TR / C_0)


def _db_sum(*levels: float) -> float:
    return 10.0 * math.log10(sum(10.0 ** (L / 10.0) for L in levels))


_THIRD_OCT_BANDS = [
    16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160,
    200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500,
]
_BAND_FACTOR = 10 ** (1.0 / 20.0)


def _assign_band(f_hz: float) -> float | None:
    for fc in _THIRD_OCT_BANDS:
        if fc / _BAND_FACTOR <= f_hz <= fc * _BAND_FACTOR:
            return fc
    return None


def _fc_label(fc: float) -> str:
    if fc >= 1000:
        return f"{fc / 1000:.1f} K"
    return f"{fc:>5.1f}" if fc % 1 else f"{int(fc):>5d}"


# Dados experimentais de referência (Ref. 13, Fig. 17)
EXP_ROT_MR = {5: 85.0, 10: 79.5, 15: 76.0, 20: 73.0}
EXP_ROT_TR = {5: 82,   10: 76,   15: 75,   20: 73  }
EXP_TOTAL  = {16: 85,  31.5: 79.5, 50: 76, 80: 73,
              100: 82,  200: 76,  315: 75, 400: 73}


def run():
    sys.stdout.reconfigure(encoding="utf-8")

    W = 70
    SEP = "=" * W

    print(SEP)
    print("EXEMPLO I — Pegg (1979) NASA TM-80200")
    print("Helicóptero em pairado: ruído rotacional + banda larga")
    print(SEP)

    print(f"\n{'Parâmetro':<26}{'Rotor Principal':>15}{'Rotor de Cauda':>15}")
    print("-" * 56)
    rows = [
        ("Mach na ponta",              f"{M_TIP:.3f}",       f"{M_TIP:.3f}"),
        ("Veloc. de rotação (rps)",    f"{N_MR:.1f}",        f"{N_TR:.1f}"),
        ("Raio da pá (m)",             f"{R_MR:.2f}",        f"{R_TR:.2f}"),
        ("Ângulo do observador β (°)", f"{BETA_MR:.0f}°",    f"{BETA_TR:.0f}°"),
        ("Número de pás",              f"{B_MR}",            f"{B_TR}"),
        ("Empuxo (N)",                 f"{T_MR:.0f}",        f"{T_TR:.0f}"),
        ("Área das pás (m²)",          f"{A_B_MR:.2f}",      f"{A_B_TR:.2f}"),
        ("Densidade do ar (kg/m³)",    f"{RHO:.3f}",         f"{RHO:.3f}"),
        ("Distância ao obs. (m)",      f"{r_MR:.1f}",        f"{r_TR:.1f}"),
        ("Coef. méd. sustent. C̄_L",   f"{CL_MR:.3f}",      f"{CL_TR:.3f}"),
    ]
    for label, mr, tr in rows:
        print(f"  {label:<24}{mr:>15}{tr:>15}")

    # --- Ruído rotacional ---
    print(f"\n{SEP}")
    print("RUÍDO ROTACIONAL  [Eq. (3)]")
    print(SEP)

    base_mr = (20.0 * math.log10(R_E_MR / r_MR)
               + 20.0 * math.log10(T_MR / (RHO * C_0**2 * R_MR**2)))
    base_tr = (20.0 * math.log10(R_E_TR / r_TR)
               + 20.0 * math.log10(T_TR / (RHO * C_0**2 * R_TR**2)))

    print(f"\n  M_e = {M_E:.2f}  |  Base (RP) = {base_mr:.1f} dB,  f₀ = {N_MR * B_MR:.1f} Hz\n")

    rot_mr = {}
    hdr = f"  {'mB':>4}  {'f (Hz)':>8}  {'ΔSPL_mr*':>9}  {'Calc. (dB)':>10}  {'Exp. (dB)':>9}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for m in [1, 2, 3, 4]:
        mB  = m * B_MR
        res = rotational_noise(COND_MR_ROT, m, lookup_fn=_dspl_mr(mB))
        f, spl = res["f_Hz"], res["SPL_dB"]
        rot_mr[mB] = (f, spl)
        print(f"  {mB:>4}  {f:>8.1f}  {_DSPL_MR[mB]:>9.1f}  {spl:>10.1f}  {EXP_ROT_MR[mB]:>9.1f}")
    print("  * ΔSPL em M_e = 0.54, β = 80°  (Ref. 13, Fig. 17)")

    print(f"\n  Base (RC) = {base_tr:.1f} dB,  f₀ = {N_TR * B_TR:.1f} Hz\n")

    rot_tr = {}
    print(hdr.replace("ΔSPL_mr*", "ΔSPL_tr*"))
    print("  " + "-" * (len(hdr) - 2))
    for m in [1, 2, 3, 4]:
        mB  = m * B_TR
        res = rotational_noise(COND_TR_ROT, m, lookup_fn=_dspl_tr(mB))
        f, spl = res["f_Hz"], res["SPL_dB"]
        rot_tr[mB] = (f, spl)
        print(f"  {mB:>4}  {f:>8.1f}  {_DSPL_TR[mB]:>9.1f}  {spl:>10.1f}  {EXP_ROT_TR[mB]:>9.1f}")
    print("  * ΔSPL em M_e = 0.54, β = 80°  (Ref. 13, Fig. 17)")

    # --- Ruído de banda larga ---
    print(f"\n{SEP}")
    print("RUÍDO DE BANDA LARGA  [Eq. (8)]")
    print(SEP)

    fp_mr = broadband_peak_frequency(T_MR, V_T_MR, si=True)
    fp_tr = broadband_peak_frequency(T_TR, V_T_TR, si=True)

    print(f"\n  f_p (RP) = {fp_mr:.0f} Hz  |  f_p (RC) = {fp_tr:.0f} Hz\n")

    print(f"  Rotor principal  [V_T = {V_T_MR:.1f} m/s, θ₁ = {THETA1_MR:.0f}°, C̄_L = {CL_MR}]")
    bb_mr = {}
    hdr_bb = f"  {'f (Hz)':>8}  {'f/f_p':>7}  {'S_1/3':>6}  {'SPL_1/3 (dB)':>13}"
    print(hdr_bb)
    print("  " + "-" * (len(hdr_bb) - 2))
    for f in [40, 80, 160, 315, 630, 1250, 2500, 5000, 10000]:
        spl = broadband_spl(COND_MR_BB, f=float(f), f_p=fp_mr, lookup_S13=_s13)
        s13 = _s13(f / fp_mr)
        bb_mr[f] = spl
        print(f"  {f:>8d}  {f/fp_mr:>7.4f}  {s13:>6.1f}  {spl:>13.1f}")

    print(f"\n  Rotor de cauda  [V_T = {V_T_TR:.1f} m/s, θ₁ = {THETA1_TR:.0f}°, C̄_L = {CL_TR}]")
    print(f"  (discrepância de ~0.7 dB na base em relação ao relatório; causa não identificada)\n")
    bb_tr = {}
    print(hdr_bb)
    print("  " + "-" * (len(hdr_bb) - 2))
    for f in [20, 40, 80, 160, 315, 630, 1250, 2500, 5000, 10000]:
        spl = broadband_spl(COND_TR_BB, f=float(f), f_p=fp_tr, lookup_S13=_s13)
        s13 = _s13(f / fp_tr)
        bb_tr[f] = spl
        print(f"  {f:>8d}  {f/fp_tr:>7.4f}  {s13:>6.1f}  {spl:>13.1f}")

    # --- Resumo combinado ---
    print(f"\n{SEP}")
    print("RESUMO EM BANDAS DE TERÇO DE OITAVA")
    print(SEP)

    mr_rot_band: dict[float, float] = {}
    for mB, (f, spl) in rot_mr.items():
        fc = _assign_band(f)
        if fc is not None:
            mr_rot_band[fc] = spl

    tr_rot_band: dict[float, float] = {}
    for mB, (f, spl) in rot_tr.items():
        fc = _assign_band(f)
        if fc is not None:
            tr_rot_band[fc] = spl

    def _col(val: float | None) -> str:
        return f"{val:>7.1f}" if val is not None else f"{'':>7}"

    hdr_sum = (f"\n  {'Freq.':>6}  {'RC Rot':>7}  {'RC BL':>7}  "
               f"{'RP Rot':>7}  {'RP BL':>7}  {'Total':>7}  {'Exp.':>6}")
    print(hdr_sum)
    print("  " + "-" * (len(hdr_sum) - 3))

    for fc in _THIRD_OCT_BANDS:
        tr_r = tr_rot_band.get(fc)
        tr_b = bb_tr.get(fc)
        mr_r = mr_rot_band.get(fc)
        mr_b = bb_mr.get(fc)
        parts = [v for v in [tr_r, tr_b, mr_r, mr_b] if v is not None]
        total = _db_sum(*parts) if parts else None
        exp   = EXP_TOTAL.get(fc)
        exp_s = f"{exp:>6.1f}" if exp is not None else f"{'':>6}"
        print(f"  {_fc_label(fc):>6}  {_col(tr_r)}  {_col(tr_b)}  "
              f"{_col(mr_r)}  {_col(mr_b)}  {_col(total)}  {exp_s}")


if __name__ == "__main__":
    run()
