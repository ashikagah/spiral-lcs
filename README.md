# Hidden structures of information transport underlying spiral wave dynamics

The original paper is [here](http://aip.scitation.org/doi/full/10.1063/1.4973542).

## How to cite

Please cite the following paper when you use the code in this repo:

Ashikaga H and James RG. Hidden structures of information transport underlying spiral wave dynamics. _Chaos_ 27: 013106, 2017

## Installation

Clone the github repository.
```
$ git clone https://github.com/ashikagah/spiral-lcs
```
## Dependencies

1. Function `rm_spirals.m`, a MATLAB implementation of the Rogers-McCulloch model in two dimensions (2-D) from [Rogers-McCulloch repo](https://github.com/ashikagah/Rogers-McCulloch).

2. Java Information Dynamics Toolkit [JIDT](https://github.com/jlizier/jidt/wiki/Installation).

## Usage

1. Generate time series of spiral waves and information flow

In MATLAB command window, 

```
>> generate_data
```
It uses the function `rm_spirals.m` to create sequential stimulations according to the stimulation file `s4_stim.mat` to induce four spiral waves in a 2-D lattice. It saves the time series of the excitation variable (_ts_) in a file `orig60.mat` and its binarized time series in a file `bi60.mat`. It also creates and saves the time series of information flow [uo, vo] in a file `uvo60.mat`. The whole process will take several hours, depending on the system used.

2. Eulerian analysis

In MATLAB command window, 

```
>> eulerian_analysis
```
It shows Shannon entropy, instantaneous information flow and total information flow over time, all in an Eulerian perspective.  

3. Lagrangian analysis

In MATLAB command window, 

```
>> lagrangian_analysis
```
It shows the forward trajectory of information flow and the Lagrangian coherent structures in an Lagrangian perspective.

## Spatial domain
- Matrix size: 120 x 120
- Grid spacing: 0.99 mm
- Grid size: 11.9 x 11.9 cm

## Licence
MIT

## References
1. Lizier JT. JIDT: An information-theoretic toolkit for studying the dynamics of complex systems. _Frontiers in Robotics and AI_ 1:11, 2014 [html](https://journal.frontiersin.org/article/10.3389/frobt.2014.00011/full)
