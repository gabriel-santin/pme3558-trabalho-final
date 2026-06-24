"""
example_II.py - Exemplo II de Pegg (1979) NASA TM-80200.
Ruído de IPV de um modelo de rotor de duas pás em túnel de vento. Equação (6).

Um vórtice de ponta desprendido por uma pá anterior é encontrado pela pá seguinte
ao longo de um arco Δψ = 18° (5% de uma revolução). A variação de carga ΔL = 20%
de L₀ gera o espectro harmônico de IPV avaliado aqui.

Os harmônicos rotacional e de espessura NÃO são calculados neste exemplo — a IPV
domina o espectro nesta condição, conforme indicado por Pegg (1979).

Nota sobre frequências: as frequências da tabela (228, 455, 683, … Hz) são os
5.º, 10.º, … 45.º múltiplos da FPP = B·N = 45.5 Hz, sem correção Doppler.
Os índices mB (4, 10, 16, …) passados para a Fig. 8 correspondem a m × B para
m = 2, 5, 8, …, conforme a numeração de Wright (Ref. 8) em Pegg (1979).
As frequências são portanto codificadas diretamente da tabela do relatório.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pegg1979 import RotorCondition, bvi_noise

# ΔSPL_bv em Δψ/ψ₀ = 0.05 — leituras da Fig. 8 de Pegg (1979)
_DSPL_BV = {
     4: -23.0,
    10: -11.5,
    16:  -6.0,
    20:  -3.0,
    24:  -1.7,
    30:  -2.3,
    36:  -4.5,
    40:  -9.5,
    44: -22.3,
}


def _lookup_bv(mB: int, psi_ratio: float) -> float:
    return _DSPL_BV[mB]


# Parâmetros de rotor / observador / atmosfera
B     = 2
T     = 458.35       # empuxo total, N (103 lb)
M_T   = 0.433
N_RPM = 22.75        # rps
BETA  = 49.0         # graus
r     = 1.02         # m (3.35 ft)
V_F   = 21.3         # m/s (70 ft/s)
C_0   = 344.0        # m/s
RHO   = 1.228        # kg/m³
DELTA_PSI = 18.0     # graus (5% de 360°)

V_T   = M_T * C_0
R     = V_T / (2.0 * math.pi * N_RPM)
M_F   = V_F / C_0

L_0     = T / B
DELTA_L = 0.20 * L_0

COND = RotorCondition(
    B=B, R=R, R_e=0.9 * R,
    c=0.10,
    h=0.012,
    M_t=M_T, M_F=M_F, M_dd=0.82, N=N_RPM,
    T=T, r=r, beta=BETA, psi=90.0,
    rho=RHO, c_0=C_0,
    A_b=B * 0.10 * R, CL_bar=0.30,
    L_0=L_0, DeltaL=DELTA_L,
    Delta_psi=DELTA_PSI,
    theta_1=90.0 - BETA,
)

# Cada linha: (m_call, mB, freq_Hz)
# m_call é o índice m passado para bvi_noise; mB = m × B é o índice para a Fig. 8.
# freq_Hz vem diretamente da tabela do relatório (múltiplos da FPP sem Doppler).
_TABLE_ROWS = [
    (  2,  4,   228),
    (  5, 10,   455),
    (  8, 16,   683),
    ( 10, 20,   910),
    ( 12, 24,  1138),
    ( 15, 30,  1365),
    ( 18, 36,  1593),
    ( 20, 40,  1820),
    ( 22, 44,  2048),
]


def run():
    sys.stdout.reconfigure(encoding="utf-8")

    W   = 68
    SEP = "=" * W

    print(SEP)
    print("EXEMPLO II — Pegg (1979) NASA TM-80200")
    print("Ruído de IPV: modelo de rotor de duas pás em túnel de vento  [Eq. (6)]")
    print(SEP)

    print("\nParâmetros de entrada:")
    rows_in = [
        ("Empuxo total (N)",              f"{T:.2f}  (103 lb)"),
        ("Mach na ponta",                 f"{M_T:.3f}"),
        ("Velocidade de rotação (rps)",   f"{N_RPM:.2f}"),
        ("Número de pás",                 f"{B}"),
        ("Ângulo do observador β (°)",    f"{BETA:.0f}°"),
        ("Distância ao observador (m)",   f"{r:.2f}  (3.35 ft)"),
        ("Velocidade de avanço (m/s)",    f"{V_F:.1f}  (70 ft/s)"),
        ("Velocidade do som (m/s)",       f"{C_0:.0f}"),
        ("Arco de azimute IPV Δψ (°)",   f"{DELTA_PSI:.0f}°  (Δψ/ψ₀ = 0.05)"),
        ("Razão ΔL / L₀",                "0.20  (20% da sustentação de uma pá)"),
    ]
    print(f"  {'Parâmetro':<32}{'Valor'}")
    print("  " + "-" * 46)
    for label, val in rows_in:
        print(f"  {label:<32}{val}")

    cos_beta   = math.cos(math.radians(BETA))
    lift_ratio = DELTA_L / L_0

    term1 = 20.0 * math.log10(lift_ratio * cos_beta)
    term2 = 20.0 * math.log10(T * N_RPM / (RHO * C_0**3 * r))
    base  = term1 + term2 + 190.2

    print(f"\nBase [Eq. (6)]:  {term1:.1f} + {term2:.1f} + 190.2 = {base:.1f} dB  "
          f"(Δψ/ψ₀ = {DELTA_PSI/360:.2f})")

    print(f"\nEspectro harmônico de IPV:")
    print(f"\n  {'mB':>4}  {'Freq (Hz)':>10}  {'ΔSPL_bv':>9}  {'SPL (dB)':>9}")
    print("  " + "-" * 39)

    for m_call, mB, freq_hz in _TABLE_ROWS:
        res = bvi_noise(COND, m_call, _lookup_bv)
        spl = res["SPL_dB"]
        dspl = _DSPL_BV[mB]
        print(f"  {mB:>4}  {freq_hz:>10d}  {dspl:>9.1f}  {spl:>9.1f}")

    print(f"\n  * ΔSPL_bv em Δψ/ψ₀ = {DELTA_PSI/360:.2f}  |  pico em mB = 24 (f ≈ 1138 Hz)")
    print(f"    Nota: mB=44 → 76.5 dB calc. vs 76.3 dB no relatório (provável erro tipográfico).")


if __name__ == "__main__":
    run()
