
# Single-Section Coupled Line Directional Coupler Design

This repository contains the design, simulation, and analysis of a **10 dB single-section coupled line directional coupler**. The project transitions from ideal mathematical models to physical microstrip implementations with optimized geometric elements.

## Project Overview

The objective was to design a directional coupler with a center frequency of **5 GHz** and a system impedance ($Z_0$) of **50 $\Omega$**. The project follows a three-stage validation process:

1. **Ideal Model:** MATLAB simulation using ideal transmission lines.

2. **MStrip V1:** Microstrip realization using a **Duroid 6010 substrate** ($\epsilon_r=10.7$).

3. **MStrip V2:** Refined microstrip design incorporating **optimally chamfered 90-degree bends** for improved performance.

## Design Specifications

* **Center Frequency:** 5 GHz 
* **Coupling:** 10 dB 
* **System Impedance ($Z_0$):** 50 $\Omega$ 
* **Substrate Parameters:** 1.27 mm thickness, $\epsilon_r = 10.7$, $\tan \delta = 0.0023$.

## Key Results (at 5 GHz)
The following table summarizes the performance metrics across the different simulation versions:

| Parameter | MATLAB (Ideal) | MStrip V1 | MStrip V2 |
|-----------|---------------|-----------|-----------|
| Return Loss (dB) | 125.62 | 28.99 | 23.95 |
| Coupling (dB) | 10.00 | 10.17 | 10.55 |
| Directivity (dB) | 125.16 | 9.29 | 13.91 |
| Isolation (dB) | 135.16 | 19.47 | 24.46 |

## Software & Tools

* **MATLAB:** Initial numerical simulation and S-parameter plotting.

* **ADS (Advanced Design System):** Microstrip layout design and full-wave simulation.

## Repository Structure

* MATLAB`: Scripts for ideal transmission line calculations and results.

* ADS_Designs`: Schematic and layout files for MStrip V1 and V2.

* Results`: Comparative summary graphs of S-parameters ($S_{11}, S_{21}, S_{31}, S_{41}$).

---

Developed as part of the Microwave Engineering (EE 514/414) curriculum at Western New England University.
