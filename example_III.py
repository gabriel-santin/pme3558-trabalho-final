"""
example_III.py - Exemplo III de Pegg (1979) NASA TM-80200.
Ruído a alta velocidade de ponta: modelo 1/7 do rotor UH-1 em pairado.
Compressibilidade (Eq. 4), espessura (Eq. 5) e rotacional (Eq. 3)
para dois casos: M_e = 0.9 e M_e = 0.8.

Erros tipográficos aparentes no relatório original:
  1. Base de compressibilidade = −10.4 − 42.8 − 21.6 = −74.8 dB; o relatório
     imprime "+74.8 + ΔSPL_C" (sinal negativo omitido na impressão).
  2. ΔSPL_T(mB=4, M_e=0.8) = 122.1 (lido na Fig. 5) resulta em SPL = 98.9 dB,
     0.7 dB abaixo do valor da tabela (99.6 dB) — provável erro tipográfico.

Nota de arredondamento: a espessura em M_e=0.9 apresenta +0.1 dB sistemático
em relação ao relatório. O relatório arredonda a base intermediária para −24.6
antes de somar ΔSPL_T + ΔSPL_K; o código usa a base exata (≈ −24.641).

Observações:
  - N e M_e são dados primários; V_T = 2π × N × R.
  - M_e é passado como M_t no RotorCondition (pairado: M_F=0 → effective_mach() retorna M_t).
  - Ruído rotacional: forças em 0.9R → M_e_rot = 0.9 × M_e (RotorCondition dedicado).
  - Cada caso usa seu empuxo físico: M_e=0.9 → T=485 N; M_e=0.8 → T=676 N.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pegg1979 import RotorCondition, compressibility_noise, thickness_noise, rotational_noise
from dataclasses import replace

# Dados de ΔSPL — leituras das Figs. 4, 5, 6 e 2 de Pegg (1979)

# ΔSPL_C (compressibilidade) em (mB, M_e=0.9), β=0° — Fig. 4
_DSPL_C = {
     2: 173.9,
     4: 175.3,
     8: 176.0,
    12: 175.7,
    16: 174.9,
    20: 173.9,
    30: 170.2,
}

# ΔSPL_T (espessura) em (mB, M_e), β=0° — Fig. 5
_DSPL_T = {
    ( 2, 0.9): 119.6,
    ( 4, 0.9): 125.2,
    ( 8, 0.9): 129.7,
    (12, 0.9): 131.6,
    (16, 0.9): 132.4,
    (20, 0.9): 132.7,
    (30, 0.9): 132.2,
    ( 2, 0.8): 117.9,
    ( 4, 0.8): 122.1,   # Fig. 5; tabela do relatório indica SPL=99.6 (discrepância de 0.7 dB)
    ( 8, 0.8): 123.9,
    (12, 0.8): 123.2,
    (16, 0.8): 121.7,
    (20, 0.8): 119.6,
    (30, 0.8): 113.4,
}

# ΔSPL_K (correção corda-diâmetro) em mB — Fig. 6; c/D = 0.0762/(2×1.05) ≈ 0.036
_DSPL_K = {
     2: 3.5,
     4: 3.5,
     8: 3.4,
    12: 3.2,
    16: 3.1,
    20: 3.0,
    30: 2.5,
}

# ΔSPL_mr (rotacional) em (mB, M_e_rot), β=0° — Fig. 2
# M_e_rot = 0.9 × M_e: 0.81 para M_e=0.9, 0.72 para M_e=0.8
_DSPL_MR = {
    ( 2, 0.81): 155.5,
    ( 4, 0.81): 156.7,
    ( 8, 0.81): 157.3,
    (12, 0.81): 156.6,
    (16, 0.81): 155.5,
    (20, 0.81): 154.5,
    (30, 0.81): 152.4,
    ( 2, 0.72): 154.3,
    ( 4, 0.72): 154.2,
    ( 8, 0.72): 153.3,
    (12, 0.72): 150.6,
    (16, 0.72): 147.7,
    (20, 0.72): 145.9,
    (30, 0.72): 142.3,
}

# Callables de consulta
def _lookup_C(mB: int, M_e: float, beta_deg: float) -> float:
    return _DSPL_C[mB]

def _lookup_T(mB: int, M_e: float, beta_deg: float) -> float:
    return _DSPL_T[(mB, round(M_e, 1))]

def _lookup_K(mB: int, c_over_D: float) -> float:
    return _DSPL_K[mB]

# Parâmetros físicos
B    = 2
R    = 1.05
R_E  = 0.945         # 0.9R
c    = 0.0762
h    = 0.12 * c      # espessura = 12% da corda
r    = 3.14
BETA = 0.0           # no plano
C_0  = 344.0
RHO  = 1.228
M_DD = 0.8

D    = 2.0 * R

ME_09  = 0.9
N_09   = 47.2
T_09   = 485.0

ME_08  = 0.8
N_08   = 41.95
T_08   = 676.0

ME_ROT_09 = 0.9 * ME_09   # 0.81
ME_ROT_08 = 0.9 * ME_08   # 0.72

# mB, m_index, freq_09_Hz, freq_08_Hz (frequências da tabela do relatório, sem Doppler)
_ROWS = [
    ( 2,  1,   94,   84),
    ( 4,  2,  188,  168),
    ( 8,  4,  376,  336),
    (12,  6,  564,  503),
    (16,  8,  752,  671),
    (20, 10,  940,  839),
    (30, 15, 1410, 1258),
]

# RotorConditions — compressibilidade e espessura usam M_t = M_e;
# ruído rotacional usa M_t = M_e_rot = 0.9 × M_e (forças em 0.9R)
def _make_cond(M_t, N, T):
    return RotorCondition(
        B=B, R=R, R_e=R_E,
        c=c, h=h,
        M_t=M_t, M_F=0.0, M_dd=M_DD, N=N,
        T=T, r=r, beta=BETA, psi=90.0,
        rho=RHO, c_0=C_0,
        A_b=B * c * R, CL_bar=0.30,
        L_0=T / B, DeltaL=0.0,
        Delta_psi=0.0,
        theta_1=90.0 - BETA,
    )

COND_09 = _make_cond(M_t=ME_09, N=N_09, T=T_09)
COND_08 = _make_cond(M_t=ME_08, N=N_08, T=T_08)

COND_09_ROT = _make_cond(M_t=ME_ROT_09, N=N_09, T=T_09)
COND_08_ROT = _make_cond(M_t=ME_ROT_08, N=N_08, T=T_08)

def _fig2(mB: int, M_e_key: float):
    key = round(M_e_key, 2)
    def _fn(M_e, beta_deg):
        return _DSPL_MR[(mB, key)]
    return _fn


def _energetic_sum(spls):
    return 10.0 * math.log10(sum(10.0 ** (s / 10.0) for s in spls))


def run():
    sys.stdout.reconfigure(encoding="utf-8")

    W   = 72
    SEP = "=" * W

    print(SEP)
    print("EXEMPLO III — Pegg (1979) NASA TM-80200")
    print("Ruído a alta velocidade de ponta: modelo 1/7 do rotor UH-1 em pairado")
    print("Compressibilidade (Eq. 4), Espessura (Eq. 5), Rotacional (Eq. 3)")
    print(SEP)

    print("\nParâmetros de entrada:")
    rows_in = [
        ("Número de pás B",                    f"{B}"),
        ("Raio da pá R (m)",                   f"{R}"),
        ("Raio efetivo R_e (m)",               f"{R_E}  (= 0.9R)"),
        ("Corda c (m)",                        f"{c}"),
        ("Espessura h (m)",                    f"{h:.6f}  (= 0.12c)"),
        ("Mach de divergência de arrasto M_dd",f"{M_DD}"),
        ("Distância ao observador r (m)",      f"{r}"),
        ("Elevação do observador β (°)",       f"{BETA:.0f}°  (no plano)"),
        ("Velocidade do som c₀ (m/s)",        f"{C_0:.0f}"),
        ("Densidade do ar ρ (kg/m³)",         f"{RHO}"),
        ("M_e=0.9: N (rps), T (N)",           f"{N_09}  |  {T_09:.0f}"),
        ("M_e=0.8: N (rps), T (N)",            f"{N_08}  |  {T_08:.0f}"),
    ]
    print(f"  {'Parâmetro':<38}{'Valor'}")
    print("  " + "-" * 54)
    for label, val in rows_in:
        print(f"  {label:<38}{val}")

    # --- Compressibilidade (M_e = 0.9 somente) ---
    print(f"\n{SEP}")
    print("RUÍDO DE COMPRESSIBILIDADE  [Eq. (4)]")
    print(SEP)

    t1_C = 20.0 * math.log10(R_E / r)
    t2_C = 20.0 * math.log10((ME_09 - M_DD) * c / R)
    base_C = t1_C + t2_C - 21.6

    f0_09 = B * N_09

    print(f"\nM_e = 0.9  (ativo, M_e > M_dd = {M_DD})  |  Base = {base_C:.1f} dB  |  f₀ = {f0_09:.1f} Hz")
    print(f"  (relatório imprime \"+74.8\" — sinal negativo omitido na impressão)")

    print(f"\n  {'mB':>4}  {'Freq (Hz)':>10}  {'ΔSPL_C':>8}  {'SPL_C (dB)':>12}")
    print("  " + "-" * 43)

    spls_C = []
    for mB, m, f09, _ in _ROWS:
        res = compressibility_noise(COND_09, m, _lookup_C)
        if res["active"]:
            spl = res["SPL_dB"]
            dspl = _DSPL_C[mB]
            spls_C.append(spl)
            print(f"  {mB:>4}  {f09:>10d}  {dspl:>8.1f}  {spl:>12.1f}")

    total_C = _energetic_sum(spls_C)
    print(f"  {'Total (energético)':>4}  {'':>10}  {'':>8}  {total_C:>12.1f}")
    print(f"\n  M_e = 0.8  (M_e = M_dd → compressibilidade inativa)")

    # --- Espessura ---
    print(f"\n{SEP}")
    print("RUÍDO DE ESPESSURA  [Eq. (5)]")
    print(SEP)

    t1_T09 = 40.0 * math.log10(ME_09)
    t1_T08 = 40.0 * math.log10(ME_08)
    t2_T   = 20.0 * math.log10(R / r)
    t3_T   = 20.0 * math.log10(h / c)
    t4_T   = 20.0 * math.log10(B)
    base_T09 = t1_T09 + t2_T + t3_T + t4_T - 0.9
    base_T08 = t1_T08 + t2_T + t3_T + t4_T - 0.9

    f0_08 = B * N_08

    print(f"\n  Base (M_e=0.9) = {base_T09:.1f} dB  |  Base (M_e=0.8) = {base_T08:.1f} dB")
    print(f"  f₀ (M_e=0.9) = {f0_09:.1f} Hz  |  f₀ (M_e=0.8) = {f0_08:.1f} Hz")
    print(f"  Nota: ΔSPL_T(mB=4, M_e=0.8) = 122.1 → SPL=98.9 dB; relatório: 99.6 dB (erro tipográfico).")
    print(f"  Nota: espessura M_e=0.9 +0.1 dB em relação ao relatório (arredondamento intermediário).")

    print(f"\n  {'mB':>4}  {'f_09':>6}  {'ΔSPL_T09':>9}  {'ΔK':>5}  {'SPL09':>8}"
          f"  {'f_08':>6}  {'ΔSPL_T08':>9}  {'SPL08':>8}")
    print("  " + "-" * 72)

    spls_T09, spls_T08 = [], []
    for mB, m, f09, f08 in _ROWS:
        res09 = thickness_noise(COND_09, m, _lookup_T, _lookup_K)
        spl09 = res09["SPL_dB"]
        res08 = thickness_noise(COND_08, m, _lookup_T, _lookup_K)
        spl08 = res08["SPL_dB"]
        spls_T09.append(spl09)
        spls_T08.append(spl08)
        dspl_t09 = _DSPL_T[(mB, 0.9)]
        dspl_t08 = _DSPL_T[(mB, 0.8)]
        dspl_k   = _DSPL_K[mB]
        print(f"  {mB:>4}  {f09:>6d}  {dspl_t09:>9.1f}  {dspl_k:>5.1f}  {spl09:>8.1f}"
              f"  {f08:>6d}  {dspl_t08:>9.1f}  {spl08:>8.1f}")

    total_T09 = _energetic_sum(spls_T09)
    total_T08 = _energetic_sum(spls_T08)
    print(f"  {'Total':>4}  {'':>6}  {'':>9}  {'':>5}  {total_T09:>8.1f}"
          f"  {'':>6}  {'':>9}  {total_T08:>8.1f}")

    # --- Ruído rotacional ---
    print(f"\n{SEP}")
    print("RUÍDO ROTACIONAL  [Eq. (3)]")
    print(SEP)

    t1_R   = 20.0 * math.log10(R_E / r)
    t2_R09 = 20.0 * math.log10(T_09 / (RHO * C_0**2 * R**2))
    t2_R08 = 20.0 * math.log10(T_08 / (RHO * C_0**2 * R**2))
    base_R09 = t1_R + t2_R09
    base_R08 = t1_R + t2_R08

    print(f"\n  M_e_rot: {ME_ROT_09:.2f} (M_e=0.9)  |  {ME_ROT_08:.2f} (M_e=0.8)")
    print(f"  Base (M_e=0.9) = {base_R09:.1f} dB  |  Base (M_e=0.8) = {base_R08:.1f} dB")

    print(f"\n  {'mB':>4}  {'f_09':>6}  {'ΔSPL(0.81)':>11}  {'SPL09':>8}"
          f"  {'f_08':>6}  {'ΔSPL(0.72)':>11}  {'SPL08':>8}")
    print("  " + "-" * 68)

    spls_R09, spls_R08 = [], []
    for mB, m, f09, f08 in _ROWS:
        res09 = rotational_noise(COND_09_ROT, m, _fig2(mB, ME_ROT_09))
        spl09 = res09["SPL_dB"]

        res08 = rotational_noise(COND_08_ROT, m, _fig2(mB, ME_ROT_08))
        spl08 = res08["SPL_dB"]

        spls_R09.append(spl09)
        spls_R08.append(spl08)

        dspl09 = _DSPL_MR[(mB, round(ME_ROT_09, 2))]
        dspl08 = _DSPL_MR[(mB, round(ME_ROT_08, 2))]

        print(f"  {mB:>4}  {f09:>6d}  {dspl09:>11.1f}  {spl09:>8.1f}"
              f"  {f08:>6d}  {dspl08:>11.1f}  {spl08:>8.1f}")

    total_R09 = _energetic_sum(spls_R09)
    total_R08 = _energetic_sum(spls_R08)
    print(f"  {'Total':>4}  {'':>6}  {'':>11}  {total_R09:>8.1f}"
          f"  {'':>6}  {'':>11}  {total_R08:>8.1f}")

    # --- Soma dos componentes ---
    print(f"\n{SEP}")
    print("SOMA DOS COMPONENTES DE RUÍDO")
    print(SEP)

    _REF_08 = {2:97.2, 4:100.5, 8:101.2, 12:100.1, 16:98.4, 20:96.3, 30:89.9}
    _REF_09 = {2:102.6, 4:106.1, 8:109.4, 12:110.7, 16:111.3, 20:111.4, 30:110.2}

    print(f"\n  Soma energética (potência) de compressibilidade + espessura + rotacional:")
    print(f"  (Compressibilidade inativa para M_e=0.8)\n")
    print(f"  {'mB':>4}  {'f_09':>6}  {'M_T=0.9':>9}  {'Ref(0.9)':>10}"
          f"  {'f_08':>6}  {'M_T=0.8':>9}  {'Ref(0.8)':>10}")
    print("  " + "-" * 70)

    for i, (mB, m, f09, f08) in enumerate(_ROWS):
        spl09 = _energetic_sum([spls_C[i], spls_T09[i], spls_R09[i]])
        spl08 = _energetic_sum([spls_T08[i], spls_R08[i]])
        ref09 = _REF_09[mB]
        ref08 = _REF_08[mB]
        print(f"  {mB:>4}  {f09:>6d}  {spl09:>9.1f}  {ref09:>10.1f}"
              f"  {f08:>6d}  {spl08:>9.1f}  {ref08:>10.1f}")

    print(f"\n  Totais por fonte (soma energética sobre todos os mB):")
    print(f"  {'':>36} M_e=0.9    M_e=0.8")
    print(f"  {'Compressibilidade (dB)':<36} {total_C:>6.1f}   inativa")
    print(f"  {'Espessura (dB)':<36} {total_T09:>6.1f}   {total_T08:>6.1f}")
    print(f"  {'Rotacional (dB)':<36} {total_R09:>6.1f}   {total_R08:>6.1f}")

    combined_09 = _energetic_sum(spls_C + spls_T09 + spls_R09)
    combined_08 = _energetic_sum(spls_T08 + spls_R08)
    print(f"  {'Combinado (dB)':<36} {combined_09:>6.1f}   {combined_08:>6.1f}")



if __name__ == "__main__":
    run()
