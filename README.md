# Pegg (1979) — Previsão de Ruído de Rotor de Helicóptero

Implementação em Python dos procedimentos semi-empíricos de previsão de ruído de
helicópteros de Pegg (1979), NASA TM-80200.

---

## Mecanismos implementados

| # | Mecanismo | Fonte | Equação(ões) |
|---|-----------|--------|--------------|
| §1 | Ruído rotacional | Lowson & Ollerhead (Ref. 3) | (1), (2), (3) |
| §2 | Ruído de compressibilidade | Arndt & Borgman (Ref. 6) | (1), (4) |
| §3 | Ruído de espessura | Hawkings & Lowson (Ref. 7) | (1), (5) |
| §4 | Interação pá-vórtice (IPV) | Wright (Ref. 8) | (6) |
| §5 | Ruído de instalação | — | sem procedimento |
| §6 | Ruído de banda larga | Lowson / Hubbard / Schlegel / Munch | (7), (8) |

---

## Requisitos

- Python 3.10 ou superior
- Apenas biblioteca padrão (`math`, `dataclasses`)

---

## Estrutura dos arquivos

- `pegg1979.py` — biblioteca principal com `RotorCondition` e as funções de ruído
- `example_I.py` — helicóptero em pairado: ruído rotacional + banda larga (Eqs. 2, 3, 7, 8)
- `example_II.py` — rotor de duas pás em túnel de vento: IPV (Eq. 6)
- `example_III.py` — modelo 1/7 do UH-1 em pairado: compressibilidade, espessura, rotacional (Eqs. 3–5)

Para executar qualquer exemplo:

```bash
python example_I.py
python example_II.py
python example_III.py
```

---

## Uso da biblioteca

Cada função de ruído exige uma callable de consulta fornecida pelo chamador
que retorna o valor de ΔSPL correspondente (lido dos gráficos de Pegg).

```python
from pegg1979 import RotorCondition, rotational_noise, broadband_peak_frequency

cond = RotorCondition(
    B=5, R=9.45, R_e=8.505, c=0.74, h=0.0,
    M_t=0.54, M_F=0.0, M_dd=0.82, N=3.5,
    T=69420.0, r=61.6, beta=5.0, psi=90.0,
    rho=1.228, c_0=344.0,
    A_b=18.58, CL_bar=0.438,
    L_0=0.0, DeltaL=0.0,
)

# Fornecer os ΔSPL lidos dos gráficos do relatório
def dspl_mr(M_e, beta_deg):
    return 148.4  # Fig. 2, harmônico mB=5, M_e=0.54, β=80°

result = rotational_noise(cond, m=1, lookup_fn=dspl_mr)
# result = {"M_e": 0.54, "f_Hz": 17.5, "SPL_dB": 85.7}
```

---

## Referência de parâmetros

| Parâmetro | Símbolo | Unidade | Descrição |
|-----------|---------|---------|-----------|
| `B` | $B$ | — | Número de pás |
| `R` | $R$ | m | Raio do rotor |
| `R_e` | $R_e$ | m | Raio efetivo da pá ($= x R$) |
| `c` | $c$ | m | Corda da pá |
| `h` | $h$ | m | Espessura máxima da pá |
| `M_t` | $M_t$ | — | Mach na ponta ($= V_T / c_0$) |
| `M_F` | $M_F$ | — | Mach de avanço ($= V / c_0$) |
| `M_dd` | $M_{dd}$ | — | Mach de divergência de arrasto |
| `N` | $N$ | rev/s | Velocidade de rotação |
| `T` | $T$ | N | Empuxo total |
| `r` | $r$ | m | Distância rotor–observador |
| `beta` | $\beta$ | graus | Elevação do observador (0–80°) |
| `psi` | $\psi$ | graus | Azimute do rotor |
| `theta_1` | $\theta_1$ | graus | Ângulo eixo de empuxo–observador |
| `rho` | $\rho$ | kg/m³ | Densidade do ar |
| `c_0` | $c_0$ | m/s | Velocidade do som |
| `A_b` | $A_b$ | m² | Área total das pás |
| `CL_bar` | $\bar{C}_L$ | — | Coeficiente médio de sustentação (≤ 0.48) |
| `L_0` | $L_0$ | N | Sustentação média da pá |
| `DeltaL` | $\Delta L$ | N | Variação de sustentação por IPV |
| `Delta_psi` | $\Delta\psi$ | graus | Arco de azimute de IPV (padrão: 18°) |

---

## Faixas de validade

| Procedimento | Restrição |
|---|---|
| `rotational_noise` | $0.5 \le M_e \le 0.9$ ; $0° \le \beta \le 80°$ |
| `compressibility_noise` | $M_{dd} < M_e \le 0.95$ (inativo se $M_e \le M_{dd}$) |
| `thickness_noise` | $0.8 \le M_e \le 0.95$ |
| `bvi_noise` | $0.02 \le \Delta\psi/\psi_0 \le 0.05$ |

---

## Limitações

- Procedimentos semi-empíricos calibrados com dados de 1979; destinam-se a
  projeto preliminar e estudos paramétricos, não a análises de certificação.
- A equação de frequência de pico (Eq. 7) é dimensionalmente inconsistente —
  use `si=True` para entradas em SI e `si=False` para unidades inglesas.
- Sem procedimento disponível para ruído de instalação (§5).
