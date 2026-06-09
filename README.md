# THMM: Penalized Threshold Hidden Markov Model for Narwhal (Monodon monoceros) Disturbance Response
This repository contains the code associated with the paper:
**Estimating the distance at which narwhal respond to disturbance: a penalized threshold hidden Markov model**
by Fanny Dupont, Marianne Marcoux, Nigel E. Hussey, Jackie Dawson, Marie Auger-Méthé
---
## Repository Structure
| File Name                                  | Description                                                                                                      |
|--------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| `Fit_THMM_casestudy.R`                     | Script to reproduce the results of the narwhal case study                                                        |
| `Fit_THMM_casestudy_sensitivity.R`         | Sensitivity analysis version of the case study script, identical to `Fit_THMM_casestudy.R` but with a modified threshold |
| `Fit_THMM_Simulation.R`                    | Script to reproduce the simulation study                                                                         |
| `SourceFunctions.R`                        | Core functions for fitting the THMM (gamma state-dependent or von Mises)                                         |
| `Tutorial.html`                            | Step-by-step tutorial of the case study with detailed explanations and more heavily commented code               |
| `data_casestudy.csv`                       | Data used in the case study                                                                                      |
---
## Key Features
- **Penalized Threshold Hidden Markov Model (THMM):** Designed to estimate disturbance response thresholds in narwhal behavior.
- **Flexible Distributions:** Currently supports gamma (state-dependent) and von Mises distributions, but is easily adaptable to others.
- **Reproducibility:** All results from the paper can be reproduced using the provided scripts.
