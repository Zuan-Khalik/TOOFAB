# TOOFAB(ulous) (Toolbox for FAst Battery simulation)
Version 1.5.1 (issued 11-10-2021)

A fast implementation of the Doyle-Fuller-Newman (DFN) battery model usable for analysis and control. 

## Features

### Current features
- Accurate and computationally efficient simulation of the DFN model describing the internal electrochemical processes of a battery. 
- Ability to easily use and change concentration-dependent parameters.
- Ability to easily configure the toolbox for a desired trade-off between accuracy and computational efficiency through selection of model simplifications and adjusting the coarseness of the discretization grid. 
- Option to use a battery ageing model describing the side reactions that lead to Li-ion loss due to the build-up SEI (solid-electrolyte interphase) layer.
- Option to use a lumped thermal model describing the thermal dynamics.  
- Ability to efficiently use the model in a closed-loop setting. 
- NEW: Estimate DFN model parameters based on experimental current/voltage data (work in progress)

### Coming features
- Documentation for the parameter_determination function

Any other desired features can still be requested. In that case, please contact the author (z.khalik@tue.nl).

## Getting Started
These instructions will set you up to use TOOFAB.

### Prerequisites 
This toolbox only requires a working version of MATLAB. 
The toolbox has been tested with MATLAB R2020b, but should work with any MATLAB version equal to or newer than MATLAB R2016b. This compatibility requirement comes from the feature that allows local functions, added to MATLAB since version R2016b. A legacy version compatible with older MATLAB versions is planned to be added in the future, or upon request. 

### Using the toolbox
Simulation of the DFN model using TOOFAB can be done with the DFN function defined as

out = DFN(input_current,tf,init_cond,param)

where

- out: contains all the output variables, such as the output voltage, the concentrations and the potentials.
- input_current: contains information about the current profile. This field can be provided either as a scalar representing the desired applied current from time 0 to final_time, an array which contains the current levels at each specified sample time, or as a function which takes the output voltage, current, concentration and potentials, and the parameters as input and mainly provides the current as output. The latter form is especially useful when the battery is desired to be controlled in closed-loop. Example functions for input_current are provided with the toolbox.
- final_time: specifies the simulation time
- init_cond: specifies the initial condition, which can be either an initial state-of-charge, as a value between 0 and 1, an initial voltage, or a MATLAB struct where the initial condition for a non-steady-state c_s, c_e, and T can be specified. Further details on how init_cond can be specified can be found in the documentation of the toolbox. 
- param: can be used to change user-configurable parameters, such as all the model parameters, and simulation parameters, e.g., the temporal and spatial grid discretization variables. Note that this field is optional, and a default set of parameters is already contained in the DFN function. 

Parameter estimation of the DFN model using TOOFAB can be done with the parameter_determination function defined as

[p, ph, results] = parameter_determination(equil_data,est_data,options_in) 

where

- p: contains the estimated DFN model parameters.
- ph: contains the estimated reformulated DFN model parameters (see [2] for more details on the reformulated DFN model).
- results: contains additional information on the estimation results.
- equil_data is a struct that contains the information on the equilbrium potential of the cell (or the EMF). When the EMF is used, there should be 2 fields in equil: EMF and Cbat. EMF should be function of SOC and Cbat should be given in C. 
Several example files are included that show how to use the toolbox. 
- est_data is a cell that can contain multiple structs that in itself contain the measurement data. In each of the structs, there should be at least 3 fields: t, I, V. t is the time vector, I contains the applied current measurements and V the terminal voltage measurements. For simplicity, make sure that you interpolate the data such that the sample time is 1 s (see the example data). 
- options: specify additional options for the parameter estimation routine. 

## Issues
It could happen that you encounter errors or unexpected behavior under certain circumstances. If you find any problems with the toolbox, please contact me (z.khalik@tue.nl). I am open to assist you with the use of the toolbox for your usecase. 

## Authors
Zuan Khalik (https://www.tue.nl/en/research/researchers/zuan-khalik/)

## References
[1] Z. Khalik, M.C.F. Donkers, H.J. Bergveld, "Model Simplifications and Their Impact on Computational Complexity for an Electrochemistry-Based Battery Modeling Toolbox", in Journal of Power Sources, 2021
[2] Z. Khalik, M.C.F. Donkers, H.J. Bergveld, "Parameter estimation of the Doyle–Fuller–Newman model for Lithium-ion batteries by parameter normalization, grouping, and sensitivity analysis", in Journal of Power Sources, 2021

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE.md](LICENSE.md) file for details


