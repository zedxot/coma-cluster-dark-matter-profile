# Dynamically Mapping the Coma Cluster: Testing $\Lambda\text{CDM}$ Halo Models & Virial Mass with DESI Spectroscopy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![DESI Iron](https://img.shields.io/badge/Data-DESI%20Iron-orange.svg)](https://data.desi.lbl.gov/)
[![Paper PDF](https://img.shields.io/badge/Paper-AAS%20TeX%20PDF-red.svg)](paper/aas_template.pdf)

## Overview

This repository presents a complete dynamical mass analysis of the **Coma Cluster (Abell 1656)** using spectroscopic galaxy kinematics from the **Dark Energy Spectroscopic Instrument (DESI)** Iron Data Release. By isolating a high-purity sample of $N = 2,098$ confirmed member galaxies within a projected angular radius of $\theta \leq 2.0^\circ$ ($\approx 3.45\text{ Mpc}$), we dynamically measure the cluster's line-of-sight velocity dispersion profile $\sigma_v(\theta)$, test competing dark matter halo density profiles ($\Lambda\text{CDM}$ cuspy NFW vs. Singular Isothermal Sphere vs. Soft-Core), and determine Coma's enclosed mass distribution and total virial mass.

---

## Key Results & Figures

### 1. Dark Matter Halo Model Selection ($\chi^2$ & BIC Analysis)

By partitioning the kinematic sample into 8 radial annuli ($\Delta \theta = 0.25^\circ$) and subtracting spectrographic measurement uncertainties in quadrature, we extract an intrinsic velocity dispersion profile declining from $\sigma_{\text{int}} = 1075.94 \pm 40.72\text{ km s}^{-1}$ in the cluster core down to $\sigma_{\text{int}} = 807.01 \pm 85.07\text{ km s}^{-1}$ at the virial boundary.

Spherical Jeans modeling and Bayesian model selection conclusively demonstrate that the **Navarro-Frenk-White (NFW)** cuspy density profile provides an optimal fit to DESI galaxy kinematics, outperforming flat and cored profiles:

| Halo Density Model | Free Parameters ($k$) | Best-Fit Parameters | $\chi^2$ | d.o.f. | $\chi_\nu^2$ | BIC | $\Delta \text{BIC}$ | Model Verdict |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Singular Isothermal Sphere (SIS)** | 1 | $\sigma_0 = 946.99 \pm 14.69\text{ km s}^{-1}$ | 22.61 | 7 | 3.23 | 24.69 | $+16.42$ | Decisively Ruled Out |
| **Soft-Core (Plummer-like)** | 2 | $\sigma_{0,c} = 1027.35\text{ km s}^{-1}, \theta_c = 1.414^\circ$ | 7.47 | 6 | 1.24 | 11.63 | $+3.36$ | Disfavored |
| **NFW (Cuspy CDM)** | 2 | $\sigma_{\text{max}} = 1093.30 \pm 42.66\text{ km s}^{-1}, \theta_s = 1.047^\circ \pm 0.353^\circ$ | **4.11** | **6** | **0.69** | **8.27** | **0.00** | **Optimal Fit ($\Lambda\text{CDM}$)** |

![Coma Cluster Model Selection](model_comparison.png)

*Figure 1: Observed line-of-sight velocity dispersion profile of Coma member galaxies from DESI Iron (black points with red error bars) compared against best-fit SIS baseline (blue dashed), Soft-Core profile (orange dash-dotted), and cuspy NFW profile (green solid).*

---

### 2. Dynamical Mass & Virial Radius Estimation

Integrating the best-fit NFW gravitational potential assuming a cluster distance of $D \approx 98.9\text{ Mpc}$ ($z_{\text{coma}} = 0.0231, H_0 = 70\text{ km s}^{-1}\text{Mpc}^{-1}$) yields a physical scale radius of $r_s \approx 1.81\text{ Mpc}$ ($1810\text{ kpc}$).

- **Virial Mass ($r_{\text{vir}} \approx 2.0\text{ Mpc}$, $R_{200}$)**:
  $$\mathbf{M_{\text{vir}} = (5.11 \pm 0.45) \times 10^{14} \, \mathbf{M}_\odot}$$
- **Total Aperture Enclosed Mass ($\theta \le 2.0^\circ, r \approx 3.45\text{ Mpc}$)**:
  $$\mathbf{M(<2.0^\circ) = 9.57 \times 10^{14} \, \mathbf{M}_\odot}$$

![Coma Cluster Enclosed Mass Profile](coma_enclosed_mass_profile.png)

*Figure 2: Enclosed dynamical mass profile $M(<\theta)$ of the Coma Cluster derived from NFW Jeans modeling of DESI galaxy kinematics. The vertical dotted line indicates the scale radius $r_s \approx 1.05^\circ$ ($1.81\text{ Mpc}$).*

---

## Data Reduction & Methodology Pipeline

1. **Target Selection (`data_prep.py`)**:
   - Extract galaxy targets from DESI Iron primary catalogs (`zall-pix-iron.fits`).
   - Quality cuts: `SPECTYPE == 'GALAXY'`, `ZWARN == 0`, `ZCAT_PRIMARY == True`.
2. **Relativistic Kinematics & Interloper Clipping (`data_filter_calc.py`)**:
   - Convert redshifts to exact relativistic line-of-sight Doppler velocities:
     $$v_{\text{LOS}} = c \cdot \frac{(1+z)^2 - 1}{(1+z)^2 + 1}$$
   - Spatial cone centered on $(\alpha_0, \delta_0) = (194.953^\circ, +27.981^\circ)$ within $\theta \le 2.0^\circ$ and $z \in [0.013, 0.033]$.
   - Apply iterative $3\sigma$ velocity clipping ($\|v_{\text{LOS}} - v_{\text{median}}\| \le 3\sigma_v$) to reject unbound background/foreground interlopers, yielding $N = 2,098$ confirmed member galaxies.
3. **Jeans Modeling & Mass Profiling (`calc_n_model_plot.py`)**:
   - Measure intrinsic velocity dispersion $\sigma_{\text{int}}$ across 8 radial bins by subtracting mean spectrograph measurement errors in quadrature.
   - Perform non-linear least-squares fitting (`scipy.optimize.curve_fit`) to evaluate SIS, Soft-Core, and NFW halo profiles.
   - Calculate Bayesian Information Criterion ($\text{BIC} = \chi^2 + k \ln N$) and integrate the NFW potential to compute $M(<r)$.

---

## Repository Structure

```text
├── data_prep.py                # Extracts DESI Iron FITS catalog & applies primary quality flags
├── data_filter_calc.py         # Relativistic velocity calculations, spatial selection & 3σ clipping
├── calc_n_model_plot.py        # Radial profile extraction, halo model fitting & plot generation
├── coma_galaxy_data.csv        # Cleaned catalog of N=2,098 confirmed Coma member galaxies
├── model_comparison.png        # Figure 1: Observed dispersion profile vs. SIS/Core/NFW fits
├── coma_enclosed_mass_profile.png # Figure 2: NFW enclosed mass profile M(<r)
├── paper/
│   ├── aas_template.tex        # AAS TeX 6.31 LaTeX source file for the manuscript
│   └── aas_template.pdf        # Compiled publication-ready PDF paper
├── .gitignore                  # Git tracking rules (only paper/aas_template.tex tracked in paper/)
├── LICENSE                     # MIT License
└── README.md                   # Project overview & documentation
```

---

## Quickstart & Usage

### Prerequisites
- Python 3.10+
- `numpy`, `pandas`, `scipy`, `matplotlib`, `astropy`

### Environment Setup
```bash
# Clone repository
git clone https://github.com/zedxot/coma-cluster-dark-matter-profile.git
cd coma-cluster-dark-matter-profile

# Create & activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install requirements
pip install numpy pandas scipy matplotlib astropy
```

### Running the Analysis Pipeline
```bash
# 1. Perform velocity calculations and interloper rejection
python data_filter_calc.py

# 2. Run Jeans halo fitting, BIC model selection, and generate figures
python calc_n_model_plot.py
```

### Compiling the Manuscript
```bash
cd paper
pdflatex -shell-escape -interaction=nonstopmode aas_template.tex
```

---

## Citation & Paper

If you use this codebase or data products in your research, please cite the manuscript in [`paper/aas_template.tex`](paper/aas_template.tex):

```bibtex
@article{Hossain2026ComaDESI,
  title={Dynamical Mapping of the Coma Cluster: Testing $\Lambda$CDM Halo Models and Virial Mass with DESI Spectroscopy},
  author={Hossain, Jubair},
  journal={AAS Journal Submission Draft},
  year={2026}
}
```

---

## License
This project is licensed under the [MIT License](LICENSE).
