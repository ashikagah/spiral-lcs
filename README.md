# Hidden structures of information transport underlying spiral wave dynamics

The original paper is [here](http://aip.scitation.org/doi/full/10.1063/1.4973542).

## Citation

Please cite the following paper when you use the code in this repo:

Ashikaga H and James RG. Hidden structures of information transport underlying spiral wave dynamics. _Chaos_ 27: 013106, 2017

## Installation

Clone the github repository.
```
$ git clone https://github.com/ashikagah/spiral-lcs
```
## Dependencies

1. Function `rm_spirals.m`, a MATLAB implementation of the Rogers-McCulloch model in two dimensions (2-D) from [Rogers-McCulloch repo](https://github.com/ashikagah/Rogers-McCulloch).

2. Java Information Dynamics Toolkit [JIDT](https://github.com/jlizier/jidt).

## Usage

1. Generate time series of spiral waves

In MATLAB command window, 

```
>> generate_data
```
It uses the function `rm_spirals.m` to create random stimulations to induce four spiral waves in a 2-D lattice. It uses a stimulation file `s4_stim.mat`. It saves the time series of the excitation variable (_ts_) in a file `orig60.mat` and its binarized time series in a file `bi60.mat`. 

2. Convert to phase map

The function `phase_map.m` will map the excitation variable (_ts_) to phase [-pi, pi] using Hilbert transform. It will save the time series of the phase (_p_) in a file `phase.mat` and a movie file `phase_movie.avi`.

3. Identify phase singularities

The function `phase_singularity.m` will detect the phase singularities (= rotors) using topological charge. It will save the time series of the phase singularities (_ps_) in a file `singularity.mat` and a movie file `singularity_movie.avi`.

## Variables
The Rogers-McCulloch model involves two normalized state variables:
- **_v_**: excitation variable (~ transmembrane potential)
- _**r**_: relaxation variable 

## Spatial domain
- Matrix size: 120 x 120
- Grid spacing: 0.99 mm
- Grid size: 11.9 x 11.9 cm

## Numerical solution
Model equations are solved using a finite difference method for spatial derivatives and explicit Euler integration for time derivatives. Newman boundary conditions are assumed. 

## Licence
MIT

## References
1. Rogers JM and McCulloch AD. A collocation-Galerkin finite element model of cardiac action potential propagation. _IEEE Trans Biomed Eng_ 41: 743-57, 1994
2. Bray MA and Wikswo JP. Use of topological charge to determine filament location and dynamics in a numerical model of scroll wave activity. _IEEE Trans Biomed Eng_ 49: 1086-93, 2002
3. Hammer P. Spiral waves in monodomain reaction-diffusion model. [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/22492-spiral-waves-in-monodomain-reaction-diffusion-model)
4. Spottiswoode BS. 2D phase unwrapping algorithms. [MATLAB File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/22504-2d-phase-unwrapping-algorithms?focused=5111677&tab=function)
5. Atienza FA et al., A probabilistic model of cardiac electrical activity based on a cellular automata system. _Rev Esp Cardiol_ 58:41-7, 2005


