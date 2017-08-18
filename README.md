# Vegetation-Banding

## Overview

The following are code excerpts as part of a Mathematics REU. Numerical simulations written in MATLAB include arclength continuation methods, direct simulation of a system of partial differential equations using finite differences, and programs for generating figures.

## Background on Project

There are many instances in nature where vegetation growth in dryland environments displays a band-like pattern formation. Our goal is to expand on the work done for a simple model between *w* the water concentration and *b* the biomass concentration. 

<img src="https://latex.codecogs.com/gif.latex?b_t=b_{xx}+wb^2-b"/>
<img src="https://latex.codecogs.com/gif.latex?w_t=cw_x-wb^2+b"/>

At the end of the project a LaTeX document will accompany the files for full exposition.

## Main Components

-simulation
   -continuation
   -direct\_sim

### simulation (most updated folder)

**continuation:** In the continuation folder there are implmementations of numerical arclength continuation used to find heteroclinic orbits for a scaled version of the system above.  

**direct\_sim:** In the direct\_sim folder there are direct simulations of the model above for both heteroclinic orbits we consider.




