# Various MATLAB Files to be used for Vegetation Bands Numerics. 

How-to-use description should be within individual files, however this is a broad overview of what each file is for and what it is capable of producing. 

## continuation_h1.m
	(Note : Since this was done before continuation_h2.m and earlier in the project, is not as complete as future files)
	- relies on cont_df.m, cont_dz.m, dz_solver.m, and use of fsolve function from optimization toolbox. 
	- contains many modes but overall purpose is to use arclength continuation to find **upper edge** heteroclinic. 
	- analyzes scaled system to ultimately produce theta (tilde) vs. s graph (mode 4). Mode 5 transforms back to original coordinates and plots discrete points from direct simulation produced from veg_dir_sim_h1.m
	- FIGURES GENERATED: 
		- Mode 1: error of jacobian components
		- Mode 2: single computation
		- Mode 3: baby continuation updates
		- Mode 4: arc length continuation
		- Mode 5: transformed data for het1

	### cont_df.m 
		- function for fsolve, produces jacobian as well for u = [b; v; s]. Does not yet include parameter theta. can debug (look for irregular errors) in continuation.m mode 1.

	### cont_dz.m 
		- function for fsolve, produces jacobian for continuation in z = [b; v; s; theta]. 

	### dz_solver.m
		- (small) function for fsolve, finds new direction dz to search in (is normalized after).


## continuation_h2.m
	-Similar to continuation_h1.m but for the **lower edge** heteroclinic. Same modes, but contains mode 6 to replot heteroclinics (potentially both 1 & 2) on a single plot and with legend.
	- Requires cont_df_h2.m, cont_dz_h2.m, and dz_solver.m.
	- Has similar functions for fsolve, they do the same things as in continuation.m but have different boundary conditions.
	- Pretty long and complex so may take a while to understand (same holds for continuation.m)

	### cont_df_h2.m 
		- same ideas as cont_df.m
	
	### cont_dz_h2.m 
		- same ideas as cont_dz.m


## veg_dir_sim_h1.m 
	- Direct simulation of original model for **upper edge** heteroclinic. 
	- Relies on the functions veg_model_f and veg_model_j for ode15s functions.
	- Plots space vs. interfaces of biomass and water. 
	- Calculates front speed s directly.
	- FIGURES GENERATED: 
		- can save initial state or final state. nothing special besides what you see when you run it.
	- data can be saved for other use (i.e. continuation.m).

	### veg_model_f 
		- uses upwind finite differences.

	### veg_model_j 
		- should be straight forward, but could always have bugs potentially.


## veg_dir_sim_h2.m 
	- Direct simulation of original model for **lower edge** heteroclinic.
	- Relies on the functions veg_model_f2 and veg_model_j for ode15s functions.
	- Same idea as veg_dir_sim_h1.m

	### veg_model_f2 
		- similar.




