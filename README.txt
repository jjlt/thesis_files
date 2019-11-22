# thesis_files
All files relevant to my Honours Thesis.
This README.md is the documentation to GitHub repository
#######################################################################
#      Name         #                  Description                    #
#######################################################################

01_simulate.R		# For rudamentory testing of complex normal 
					# generation

02_simulate.R		# Testing for cases where there is no 
					# psuedo-covariance

03_simulate.R		# Testing covariate matrices size p > 1 work

04_simulate.R		# Inspecting the results of the iterative estimation
					# method

05_simulate.R		# Initial comparison with testing taking the modulus
					# of complex data

06_simulate.R		# First important file. Simulating complex normal
					# and checking for normality of the beta estimates.
					# This file produced the graphics for the density
					# histograms. Warning that large values of n, will
					# take ages.

07_simulate.R		# Another important file where graphics were made,
					# this one created the numerical testing of the 
					# components of the regression parameters having
					# the consistency, eg they converge to the true 
					# values.

08_simulate.R		# Testing file for statistical tests. 
					# Non-substaintial file.

09_simulate.R		# Looking at the covariance of the beta estimate,
					# using true values of gamma and estimated values,
					# seeing if they line up with the theoretical amount.

10_simulate.R		# Used to produce the graphics for showing the
					# the consistency of the real regression parameters
					# numerically.

11_simulate.R		# Non-substaintial file inspecting the convergence 
					# of gamma. Turns out it does converge, just not as
					# fast.

12_simulate.R		# Attempting to do statistical tests and transforming
					# the beta estimates. Seeing if they results looked 
					# correct. Tried making insignificant data, does work.
					# p-values can be higher when there isn't a linear
					# relationship.

13_simulate.R		# Trying out the wind data for the first time. 30 
					# minute beforehand wind on current wind.

14_simulate.R		# We investigate if we can actually non-linearly solve
					# the MLE, but it looks way too disgusting, because it 
					# has to be split into real and imaginary parts as the
					# function cannot deal with them.

15_simulate.R		# Interesting file. Comparisons with other linear 
					# regression models on their predictive power on
					# wind data. lm on the modulus works pretty well,
					# actually better than clm, mlm and gls. I think it is
					# because the wind data is very simple and wind from
					# thirty minutes ago is always going to be pretty close.
					# gls worked okay, it was implemented in two ways, but 
					# it surpisingly wasn't close to clm or mlm. This would
					# be the last moment gls was relevant.

16_simulate.R		# We introduce the new time covariates and to tests
					# with clm.

17_simulate.R		# Important file. We do comparisons with the other 
					# linear regression methods, similar to 15_simulate,
					# but this time on simulated data and we test the 
					# closeness of regression co-efficient beta. clm wins
					# and gls is no where near, no worth putting into the
					# graphics to be honest. Be warned when running this
					# use smaller sample sizes for speed and compare the 
					# gls for yourself.

18_simulate.R		# Changing the previous winds to 24 hours beforehand
					# rather than 30 minutes beforehand. Makes it a more
					# interesting story.

1912794.csv			# Raw data from NOAA, had to do a lot of work finding
					# access to data, and finding access to data in a csv
					# form so I can easy clean it.

new_data.csv		# Cleaned data without time variables

new_data_1.csv		# Cleaned data with time. 

simulation_headers.R# Important file. This file is used as a header 
					# containing all the useful functions including the
					# whole clm function. Inside the code is very well
					# documented, except for the last few functions,
					# which are obvious in function. Careful the clm and
					# fgls functions are sensitive to the correct imput.

wind_data.csv		# Cleaned and rearranged data.

wind_data2.csv		# See new_data_1.csv

wind_script.pl		# The perl script which cleans the raw data. 

/06_simulate/		# Contains graphics and R workspaces, so the very long
					# processes don't have to be repeated. Works spaces with
					# error in their name contained data generated with		
					# a bug, so don't use those.

/07_simulate/		# See above, has a a section of older graphics.
					# again use working R workspace

/08_simulate/		# Just a workspace

/09_simulate/		# Just a workspace

/10_simulate/		# Graphics

/17_simulate/		# Workspace

/graphics/			# The final set of graphics used to compile the
					# thesis and presentation
