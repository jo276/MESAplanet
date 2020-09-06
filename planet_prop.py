# python script containing info and functions for planetary calculations

import math as mt

# constants
msun=1.9891e33 #grams
rsun=6.955e10  # cm
lsun=3.846e33  # erg/s
mjup=1.8986e30 # grams
rjup=6.9911e9  # cm
mearth=5.97219e27 # grams
rearth=6.371e8 # cm
mneptune=1.0243e29 #grams
rneptune=2.4622e9 # cm
au_to_cm=1.49597871e13 #cm
a=[0.0912, 0.0592]


def get_core_prop(Mass,Xiron,Xice):

# uses fortney et al. 2007ab fitting formula

# Mass is core mass in earth masses
# Radius is in earth masses
# density is in core density in cgs units

#fitting constants
	a=[0.0912, 0.0592]
	b=[0.1603, 0.0975]
	c=[0.3330, 0.2337]
	d=[0.7387, 0.4938]
	e=[0.4639, 0.3102]
	f=[1.1193, 0.7932]

	prop=[0]*2

	Xrock=1.0-Xiron

	if Xice > 0.0:
		if Xiron > 0.0:
			print("Error both ice and iron frac cannot be > 0")
			return prop;
		else:
			prop[0]=(a[0]*Xice+b[0])*(mt.log10(Mass))**2.0+(c[0]*Xice+d[0])*(mt.log10(Mass))+(e[0]*Xice+f[0])
	else:
		prop[0]=(a[1]*Xrock+b[1])*(mt.log10(Mass))**2.0+(c[1]*Xrock+d[1])*(mt.log10(Mass))+(e[1]*Xrock+f[1])

	# now calculate radius
	prop[1]=Mass*mearth/(4.0/3.0*mt.pi*(prop[0]*rearth)**3)

	return prop;
