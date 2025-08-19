
# README: Modified SWAT Model – MCN and MGW Enhancements

## Overview

This project introduces two key modifications to the SWAT (Soil and Water Assessment Tool) model (version 664) to improve hydrological simulation accuracy:

- **MCN (Modified Curve Number)**: Enhances the SCS-CN runoff estimation method by allowing flexible initial abstraction ratios (λ).
- **MGW (Modified Groundwater Recharge)**: Introduces seasonality into the aquifer recharge delay time (δ) using a sinusoidal function.

These modifications aim to better represent local hydrological conditions and seasonal groundwater dynamics, especially in regions with distinct wet and dry seasons.

---

## 1. Modified SCS-CN Theory (MCN)

### Traditional SCS-CN Equation

Runoff depth Q is calculated as:

- If P > Ia:
  Q = (P - Ia)^2 / (P - Ia + S)
- Else:
  Q = 0

Where:
- Q: runoff depth (mm)
- P: rainfall depth (mm)
- S: maximum potential retention (mm)
- Ia = 0.2S: initial abstraction

### CN-S Relationship:
CN = 25,400 / (254 + S)

### Modified CN for λ = 0.05:
Using empirical relationships from Woodward et al. (2003):

- Modified retention:
  S_0.05 = 0.8187 * S_0.2^1.15

- Modified curve number:
  CN_0.05 = 100 / ((1.1879 * (100 / CN_0.2 - 1))^1.15 + 1)

### Generalized Form for Any λ:

- Retention:
  S_λ = α * S_0.2^β

- Curve number:
  CN_λ = 100 / ((α / 25.4 * (25,400 / CN_0.2 - 254)^β) + 1)

This formulation ensures that 0 ≤ CN_λ ≤ 100, preserving theoretical bounds.

---

## 2. Modified Groundwater Recharge (MGW)

SWAT calculates aquifer recharge using:

w_rchrg,i = (1 - exp[-1/δ]) * w_seep + exp[-1/δ] * w_rchrg,i-1

Where:
- w_rchrg,i: recharge on day i
- w_seep: seepage from soil profile
- δ: delay time (days)

### Proposed Seasonal Delay Time Model

To reflect seasonal variation in groundwater levels, a sinusoidal model is proposed:

δ(i) = ((δ_max - δ_min)/2) * sin((2π(i - i_max)/365) + π/2) + ((δ_max + δ_min)/2)

Where:
- i: current day of the year
- i_max: day of lowest groundwater level
- δ_max, δ_min: max/min delay times

---

## 3. Implementation Details

### Modified Files in SWAT Source Code

- Parameter Definitions: modparm.f, allocate_parms.f
- Subroutines Modified: curno.f, surq_daycn.f, dailycn.f

### Calibration Parameters

These new parameters can be calibrated manually or via automated tools like SWAT-CUP:

| Parameter | Distribution | Range | Mean (μ) |
|-----------|--------------|-------|----------|
| α         | Uniform      | [0.20, 1.50] | 0.85 |
| β         | Uniform      | [0.50, 2.00] | 1.25 |
| λ         | Uniform      | [0.00, 0.40] | 0.20 |

For MGW, the sinusoidal model introduces three new variables: δ_min, δ_max, i_max. Only one additional parameter needs calibration, as i_max can be estimated from streamflow records.

---

## 4. Model Naming Convention

- **MCN**: Modified Curve Number method
- **MGW**: Modified Groundwater recharge method
- **MCN/GW**: Combined implementation of both MCN and MGW
