# The Transmission of U.S. Monetary Policy Shocks across Volatility Regimes

Project for **Monetary Economics** (MSc Quantitative Finance and Insurance, University of Turin, February 2026). The code is in MATLAB.

Does a monetary tightening work differently when markets are stressed? This project estimates a monthly U.S. VAR separately for **high- and low-volatility regimes** and compares how a Fed Funds shock is transmitted to activity, prices, credit spreads and equities.

## Data and regimes

- Monthly FRED data from 1962: industrial production growth, CPI inflation, the **Federal Funds Rate**, the **BAA–AAA credit spread** and **S&P 500 returns**.
- The VIX only starts in 1990, so a **VIX proxy** is built from the realised volatility of the S&P 500. It is rescaled by the average VIX / realised-volatility ratio over 1990–2020.
- A month is in the high-volatility regime when the proxy is above its 75th percentile.

## Identification

1. **Recursive (Cholesky)**, following Christiano–Eichenbaum–Evans (1998). The slow macro variables come first, then the policy rate, then the fast financial variables.
2. **Zero and sign restrictions**, following the agnostic approach of Uhlig (2005).
   - IP and inflation do not respond on impact. This is imposed with a block-diagonal rotation `H = diag(I₂, Q)`.
   - A tightening raises the FFR, widens the spread and lowers equity returns.

The lag order is chosen by AIC on the full sample (p = 10) and kept the same in both regimes. Confidence bands (16–84%) come from a residual bootstrap with 1,000 draws. IRFs run to 36 months, and the FEVD is reported at 12, 24 and 36 months.

## Findings

- A monetary tightening works **more strongly through financial conditions when volatility is high**. The credit spread response is larger and more persistent, and equities fall more.
- The responses of IP and inflation are smaller, slower and less precisely estimated, especially in the high-volatility regime.
- **The FEVD changes with the regime.** At 12 months, the share of the spread's forecast-error variance due to its own shock falls from 0.79 in low volatility to 0.44 in high volatility, while equity shocks explain about 0.42 of it.

## Files

| File | Content |
|---|---|
| `Project3.mlx` | Main live script: data, AIC, Cholesky and sign–zero IRFs for each regime, FEVD |
| `SignZeroRestrictBoot.mlx` | Bootstrap with zero and sign restrictions (block rotation, acceptance sampling) |
| `CholeskyBoot.m`, `IRF_Chol_VAR.mlx` | Cholesky IRFs with bootstrap bands |
| `VARTopicsOLS.m`, `VarStr.mlx`, `companion.mlx` | OLS VAR, regression matrices, companion form |
| `TopicsVD2.m`, `RotationMatrix.m` | FEVD, random orthonormal rotations |
| `*_export.m` | Plain-text copies of the `.mlx` live scripts, readable on GitHub |
| `excelprogetto.xlsx` | Data: FRED series, VIX proxy and regime dummy |
| `Report.pdf` | Project report (11 pages) |

## Running it

Open `Project3.mlx` in MATLAB from this folder and run it. The data file is included.
