# BattEry-Simulation-Toolbox (BEST)
Version 0.1.0 (issued 12-12-2019)

A fast implementation of the Doyle-Fuller-Newman (DFN) battery model usable for analysis and control. 

## Features
### Current features
The implementation has the following currently (tested) features:
- Simulation of the Simplified DFN (SDFN) model as presented in [1]. 

### Coming features
BEST will be continuously upgraded, and the current planned features are as follows: 
- Fast simulation of the full DFN model (without simplifications), including concentration-dependent parameters
- The ability to apply several simplifications to allow for a trade-off between accuracy and computation time. 
- Thermal dynamics through a lumped thermal model. 
- Battery ageing model 

## Getting Started
These instructions will set you up to use the BattEry Simulation Toolbox (BEST).

### Prerequisites 
This toolbox only requires a working version of MATLAB. 
The toolbox has been tested with MATLAB R2019b, but should work with any MATLAB version equal to or newer than MATLAB R2016b. This compatibility requirement comes from the feature that allows local functions, added to MATLAB since version R2016b. A legacy version compatible with older MATLAB versions is planned to be added in the future, or upon request. 

### Using the toolbox
The battery simulation toolbox can be interfaced with the DFN function defined as

out = DFN(input_current,tf,init_cond,param)

where

- out is a MATLAB struct containing all the output variables, such as the output voltage, the concentrations and the potentials.
- input_current contains information about the current profile. This field can be provided either as an array which contains the current levels at each specified time sample, or as a function which takes the output voltage, current, concentration and potentials, and the parameters as input and mainly provides the current as output. The latter form is especially useful when the battery is desired to be controlled in closed-loop.
- tf specifies the simulation time
- init_cond specifies the initial condition, which can be either an initial State-of-Charge (SoC), as a value between 0 and 1, or an initial voltage. If init_cond is chosen to be larger than 2, the toolbox assumes that init_cond specifies an initial voltage. 
- param is a struct containing all the model parameters, and simulation parameters, such as the temporal and spatial grid discretization variables.

An example file is included that shows how to use the toolbox. 

## Authors
Zuan Khalik (https://www.tue.nl/en/research/researchers/zuan-khalik/)

## References
[1] Z. Khalik, H.J. Bergveld, M.C.F. Donkers, "On Trade-offs Between Computational Complexity and Accuracy of Electrochemistry-based Battery Models", in Conference on Decision & Control, 2019

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE.md](LICENSE.md) file for details


