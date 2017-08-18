# Vegetation-Banding

## Overview

The following are code excerpts as part of a Mathematics REU. Numerical simulations written in MATLAB include arclength continuation methods, direct simulation of a system of partial differential equations using finite differences, and programs for generating figures.

## Background on Project

There are many instances in nature where vegetation growth in dryland environments displays a band-like pattern formation. Consider a closed two species diffusion-reaction system inspired by the work of Turing and Klausmeier. The goal is to expand on the work done for a simple model between *w* the water concentration and *b* the biomass concentration. 

<img src="https://latex.codecogs.com/gif.latex?b_t=b_{xx}+wb^2-b"/>
<img src="https://latex.codecogs.com/gif.latex?w_t=cw_x-wb^2+b"/>

At the end of the project a LaTeX document will accompany the files for full exposition.

## Layman's Terms

I'm looking at a model of water and plants on a hill with the water going down the hill due to the slope and the plants spreading out. There's a more complicated nonlinear relationship between how plants and water interact with each other, and there are two main ways to simulate this whole system: 

1.   Directly using a mesh for time and space.

2.   Transforming to a simpler system and following the solution carefully along a (potentially twisted and sharp) path. 

As it turns out, both cases show that there's an interesting pattern where the vegetation and water concentrations move up the hill as a block with a speed *s* that depends on how much vegetation/water we started with.

## Main Components
```
- simulation
   - continuation
   - direct_sim
- continuation
- sim_code
```
### simulation (most updated folder)

**continuation:** In the continuation folder there are implmementations of numerical arclength continuation used to find heteroclinic orbits for a scaled version of the system above.  

**direct\_sim:** In the direct\_sim folder there are direct simulations of the model above for both heteroclinic orbits we consider.




