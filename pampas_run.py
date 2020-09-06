# script to fit the mass and radius of an individual kepler planet
# by building a grid
import sys
import planet_prop as pp
import mesa_build as mb
import numpy as np
import planet_prop as pp

I=int(sys.argv[1])
J=int(sys.argv[2])
NI=int(sys.argv[3])
NJ=int(sys.argv[4])

planet_dir="/rds/general/user/jeowen/ephemeral/Young_planets_2/Model_0/planets"

Mcore_grid=np.linspace(2.5,9.625,NI)
Rock_fraction=2./3.
Rcore=(0.0592*Rock_fraction+0.0975)*(np.log10(Mcore_grid))**2.+(0.2337*Rock_fraction+0.4938)*np.log10(Mcore_grid)+(0.3102*Rock_fraction+0.7932)

period_grid=np.logspace(np.log10(11.567),np.log10(11.567),NJ) # days
day_to_sec=24.*60.*60.
G=6.67408e-8
separation=((period_grid*day_to_sec/(2.*np.pi))**2.*G*pp.msun)**(1./3.)

#planet properties
Mcore=Mcore_grid[I]
Mcore=(np.float64(Mcore).item())
rho_core=Mcore*pp.mearth/(4.0/3.0*np.pi*(Rcore[I]*pp.rearth)**3.0)
rho_core=(np.float64(rho_core).item())
sep=separation[J]/pp.au_to_cm
Tstar=4280.4 # at 3e6 years from Baraffe model for Sun-like star
Rstar=1.6776*pp.rsun # at 3e6 years from Baraffe model for Sun-like star
min_mass=(1+1e-5)*Mcore*pp.mearth/pp.msun

#calculate Teq
Teq=Tstar*(Rstar/(2.0*sep*pp.au_to_cm))**0.5

L_M_whittle=6e-8
L_M_min=6e-8
L_M_max=1e2
Xmax=10.
Xmin=0.001
NL=96
Nenv=128
Xbondi=0.333333
Tmax=6.e9
min_tcool=1.e4
max_tcool=1.e11
ES=True # evolve the star
do_not_skip_L0=True

#build starting models
mb.build_models(Mcore,rho_core,L_M_whittle,L_M_min,L_M_max,NL,Teq,Xmax,Xmin,Nenv,planet_dir,Xbondi,min_tcool,do_not_skip_L0)

#run models
mb.run_models(Teq,sep,NL,Nenv,planet_dir,min_mass,Xbondi,Tmax,ES,min_tcool,max_tcool,L_M_min,L_M_max)
