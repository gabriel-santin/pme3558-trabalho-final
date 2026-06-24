"""
pegg1979.py - Previsão semi-empírica de ruído de rotor de helicóptero.
Implementação das Equações (1)–(8) de Pegg (1979), NASA TM-80200.

Mecanismos implementados:
    §1  Ruído rotacional (Eq. 3)           — Lowson & Ollerhead
    §2  Ruído de compressibilidade (Eq. 4) — Arndt & Borgman
    §3  Ruído de espessura (Eq. 5)         — Hawkings & Lowson
    §4  Ruído de IPV (Eq. 6)               — Wright
    §5  Ruído de instalação                — sem procedimento disponível
    §6  Ruído de banda larga (Eq. 7–8)     — Lowson / Schlegel / Munch
"""

import math
from dataclasses import dataclass
from typing import Callable


@dataclass
class RotorCondition:
    """
    Parâmetros de entrada do rotor, voo e observador para os procedimentos de Pegg (1979).
    Todas as grandezas em SI, salvo indicação.
    """
    B: int          # número de pás
    R: float        # raio do rotor, m
    R_e: float      # raio efetivo da pá (= x·R), m
    c: float        # corda da pá, m
    h: float        # espessura máxima da pá, m
    M_t: float      # Mach na ponta (= V_T / c_0)
    M_F: float      # Mach de avanço (= V / c_0)
    M_dd: float     # Mach de divergência de arrasto
    N: float        # velocidade de rotação, rev/s
    T: float        # empuxo total, N
    r: float        # distância cubo-observador, m
    beta: float     # ângulo de elevação do observador, graus (0–80°)
    psi: float      # ângulo de azimute da pá, graus
    rho: float      # densidade do ar, kg/m³
    c_0: float      # velocidade do som, m/s
    A_b: float      # área total das pás, m²
    CL_bar: float   # coeficiente médio de sustentação C̄_L
    L_0: float      # sustentação média da pá, N
    DeltaL: float   # variação de sustentação por interação com vórtice, N
    Delta_psi: float = 18.0   # arco de azimute de IPV, graus
    psi_0: float    = 360.0   # azimute completo do rotor, graus
    theta_1: float  = 90.0    # ângulo eixo de empuxo negativo-observador, graus

    @property
    def V_T(self) -> float:
        return self.M_t * self.c_0

    @property
    def D(self) -> float:
        return 2.0 * self.R


def effective_mach(cond: RotorCondition) -> float:
    """Mach efetivo na ponta da pá — Eq. (1): M_e = (M_t + M_F·sen ψ) / (1 − M_F·cos β)."""
    psi  = math.radians(cond.psi)
    beta = math.radians(cond.beta)
    return (cond.M_t + cond.M_F * math.sin(psi)) / (1.0 - cond.M_F * math.cos(beta))


def blade_passage_frequency(cond: RotorCondition, m: int) -> float:
    """Frequência de passagem de pá no harmônico m com correção Doppler — Eq. (2)."""
    beta = math.radians(cond.beta)
    return m * cond.B * cond.N / (1.0 - cond.M_F * math.cos(beta))


# §1 — Ruído rotacional (Lowson & Ollerhead)

def rotational_noise(
    cond: RotorCondition,
    m: int,
    lookup_fn: Callable[[float, float], float],
) -> dict:
    """
    NPS de ruído rotacional no harmônico m — Eq. (3).
    lookup_fn(M_e, beta_deg) deve retornar ΔSPL_mr ou ΔSPL_tr.
    """
    M_e = effective_mach(cond)

    if not (0.5 <= M_e <= 0.9):
        raise ValueError(
            f"M_e = {M_e:.4f} fora da faixa [0.5, 0.9] (Pegg §1)."
        )
    if not (0.0 <= cond.beta <= 80.0):
        raise ValueError(
            f"β = {cond.beta}° fora da faixa [0°, 80°] (Pegg §1)."
        )

    f    = blade_passage_frequency(cond, m)
    dSPL = lookup_fn(M_e, cond.beta)

    spl = (
        20.0 * math.log10(cond.R_e / cond.r)
        + 20.0 * math.log10(cond.T / (cond.rho * cond.c_0**2 * cond.R**2))
        + dSPL
    )
    return {"M_e": M_e, "f_Hz": f, "SPL_dB": spl}


# §2 — Ruído de compressibilidade por arrasto de perfil (Arndt & Borgman)

def compressibility_noise(
    cond: RotorCondition,
    m: int,
    lookup_C: Callable[[int, float, float], float],
) -> dict:
    """
    NPS de ruído de compressibilidade no harmônico m — Eq. (4).
    Retorna active=False quando M_e ≤ M_dd.
    lookup_C(mB, M_e, beta_deg) deve retornar ΔSPL_C.
    """
    M_e = effective_mach(cond)

    if M_e <= cond.M_dd:
        return {"M_e": M_e, "active": False, "SPL_dB": None}

    if M_e > 0.95:
        raise ValueError(
            f"M_e = {M_e:.4f} excede o limite 0.95 (Pegg §2)."
        )

    mB   = m * cond.B
    dSPL = lookup_C(mB, M_e, cond.beta)

    spl = (
        20.0 * math.log10(cond.R_e / cond.r)
        + 20.0 * math.log10((M_e - cond.M_dd) * cond.c / cond.R)
        + dSPL
        - 21.6
    )
    return {"M_e": M_e, "active": True, "SPL_dB": spl}


# §3 — Ruído de espessura (Hawkings & Lowson)

def thickness_noise(
    cond: RotorCondition,
    m: int,
    lookup_T: Callable[[int, float, float], float],
    lookup_K: Callable[[int, float], float],
) -> dict:
    """
    NPS de ruído de espessura no harmônico m — Eq. (5).
    lookup_T(mB, M_e, beta_deg) → ΔSPL_T; lookup_K(mB, c_over_D) → ΔSPL_K.
    """
    M_e = effective_mach(cond)

    if not (0.8 <= M_e <= 0.95):
        raise ValueError(
            f"M_e = {M_e:.4f} fora da faixa [0.8, 0.95] (Pegg §3)."
        )

    mB       = m * cond.B
    c_over_D = cond.c / cond.D
    dSPL_T   = lookup_T(mB, M_e, cond.beta)
    dSPL_K   = lookup_K(mB, c_over_D)

    spl = (
        40.0 * math.log10(M_e)
        + 20.0 * math.log10(cond.R / cond.r)
        + 20.0 * math.log10(cond.h / cond.c)
        + 20.0 * math.log10(cond.B)
        + dSPL_T
        + dSPL_K
        - 0.9
    )
    return {"M_e": M_e, "SPL_dB": spl}


# §4 — Ruído de interação pá-vórtice / IPV (Wright)

def bvi_noise(
    cond: RotorCondition,
    m: int,
    lookup_bv: Callable[[int, float], float],
) -> dict:
    """
    NPS de ruído de IPV no harmônico m — Eq. (6).
    lookup_bv(mB, delta_psi_ratio) → ΔSPL_bv.
    """
    psi_ratio = cond.Delta_psi / cond.psi_0

    if not (0.02 <= psi_ratio <= 0.05):
        raise ValueError(
            f"Δψ/ψ₀ = {psi_ratio:.4f} fora da faixa [0.02, 0.05] (Pegg §4)."
        )

    mB         = m * cond.B
    beta_rad   = math.radians(cond.beta)
    lift_ratio = cond.DeltaL / cond.L_0
    dSPL       = lookup_bv(mB, psi_ratio)

    spl = (
        20.0 * math.log10(lift_ratio * math.cos(beta_rad))
        + 20.0 * math.log10(cond.T * cond.N / (cond.rho * cond.c_0**3 * cond.r))
        + dSPL
        + 190.2
    )
    return {"psi_ratio": psi_ratio, "SPL_dB": spl}


# §5 — Ruído de instalação

def installation_noise(*args, **kwargs) -> None:
    """Pegg (1979) §5 — sem procedimento verificado disponível."""
    raise NotImplementedError(
        "Pegg (1979) §5: nenhum procedimento verificado disponível para ruído de instalação."
    )


# §6 — Ruído de banda larga (Lowson / Hubbard / Schlegel / Munch)

def broadband_peak_frequency(T: float, V_T: float, si: bool = True) -> float:
    """
    Frequência de pico do ruído de banda larga — Eq. (7).
    Equação dimensionalmente inconsistente: si=True para entradas em SI (T em N, V_T em m/s),
    si=False para unidades inglesas (T em lbf, V_T em ft/s).
    """
    if si:
        return -240.0 * math.log10(T) + 2.448 * V_T + 942.0
    else:
        return -240.0 * math.log10(T) + 0.746 * V_T + 786.0


def _lift_coeff_correction(CL_bar: float) -> float:
    """Correção f(C̄_L) = 10 log(C̄_L / 0.4), válida para C̄_L ≤ 0.48."""
    if CL_bar <= 0.48:
        return 10.0 * math.log10(CL_bar / 0.4)
    raise NotImplementedError(
        "Pegg (1979): sem expressão analítica para C̄_L > 0.48."
    )


def broadband_spl(
    cond: RotorCondition,
    f: float,
    f_p: float,
    lookup_S13: Callable[[float], float],
) -> float:
    """
    NPS em banda de terço de oitava para ruído de banda larga — Eq. (8).
    lookup_S13(f_over_fp) → S_1/3 (espectro normalizado de Schlegel).
    """
    theta1_rad = math.radians(cond.theta_1)

    vel_term = 20.0 * math.log10((cond.V_T / cond.c_0) ** 3)
    dir_term = 10.0 * math.log10(
        (cond.A_b / cond.r**2) * (math.cos(theta1_rad) ** 2 + 0.1)
    )
    s13  = lookup_S13(f / f_p)
    f_CL = _lift_coeff_correction(cond.CL_bar)

    return vel_term + dir_term + s13 + f_CL + 130.0
