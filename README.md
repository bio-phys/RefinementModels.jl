RefinementModels.jl
===================

# Overview

The module provides simple models for ensemble refinement. 

The models are mainly intended for method devolpment. For these models, we generate synthetic data for calculated observables and experimental data. 

Currently, the following models are included:

* Discrete double-well [1]
* Continuous double-well
* Random Gaussian data [2]

Note that a Python implementation of the von Mises polymer model [3] is availabe at https://github.com/bio-phys/BioFF

# Usage 

The sub-folder /example/ contains scripts to generate models and save the output for subsequent reweighting with BioEn [1]. 

See https://github.com/bio-phys/BioEn for a Python/C implementation and https://github.com/bio-phys/BioEn.jl for a Julia implementation of methods to solve the BioEn optimization problem [2]. 

# References

- [1] Köfinger and Hummer, J. Chem. Phys. 143 (2015) https://aip.scitation.org/doi/10.1063/1.4937786
- [2] Köfinger et al. J. Chem. Theory and Comput. 15 (2019) https://doi.org/10.1021/acs.jctc.8b01231 
- [3] Köfinger and Hummer, Eur. Phys. J. B. 94 (2021) https://doi.org/10.1140/epjb/s10051-021-00234-4
