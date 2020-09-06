# collection of mesa functions for building and running the code

def build_starting_models(Mcore,rho_core,L_M_start,L_M_target,Teq,Xstart,Xend,Nenv,planet_dir):

	# builds Nenv models with Mcore and rho_core between Xstart and Xend
	# with inital luminosity of L_M_start

	import os
	import shutil
	import f90nml #fortran namelist tools
	import mesa_namelist as mn # for generating mesa namelist
	import planet_prop as pp   # tools for planetary properties
	import numpy as np
	import math

	start_file="3.00_Mj.mod"
	file_in="model_in.mod"
	file_out="model_out.mod"
	core_dir="put_core"
	make_dir="../make_planets"
	curr_dir=os.getcwd()

	#usful constants
	sb=5.6704e-5
	kb=1.38e-16
	curr_mass=5.695799690399999e+30 # exact starting mass of 3Mjup model
	Flux=4.0*sb*Teq**4.0
	ES=False

	Sigma=250 # heating depth in F-Sigma routine
	Thold=1e14 # holding time
	Twhittle=1e14 # mass removal time-scale

	shutil.copyfile(start_file,"%s/%s" %(core_dir,file_in))

	os.chdir(core_dir)

	os.system("./mk")

	nml=mn.build_namelist_pc(rho_core,Mcore*pp.mearth/(pp.msun))
	nml.write('inlist', force=True )

	#now run MESA
	os.system("./create_planet")

	#now copy output to make_planet
	shutil.copyfile(file_out,"%s/%s" %(make_dir,file_in))
	#move to make make_planets
	os.chdir(make_dir)
	#build mesa
	os.system("./mk")

	nml=mn.build_namelist_mp_Fsigma_hold(Thold,L_M_start,Flux,Sigma)
	nml.write('inlist', force=True )
	os.system("./create_planet")

	Xfrac=np.logspace(math.log10(Xstart),math.log10(Xend),Nenv)
	final_mass=curr_mass

	for i in range(0, Nenv):

		curr_mass=final_mass
		final_mass=(1.0+Xfrac[i])*Mcore*pp.mearth
		mdot=(final_mass-curr_mass)/Twhittle/pp.msun
		mdot=(np.float64(mdot).item())
		final_mass=(np.float64(final_mass).item())

		shutil.copyfile(file_out,file_in)
		nml=mn.build_namelist_mp_whittle(Twhittle,L_M_start,final_mass/pp.msun,mdot,Flux,Sigma)
		nml.write('inlist', force=True )
		os.system("./create_planet")

		shutil.copyfile(file_out,file_in)
		nml=mn.build_namelist_mp_hold(Thold,L_M_start,Flux,Sigma)
		nml.write('inlist', force=True )
		os.system("./create_planet")

		shutil.copyfile(file_out,file_in)
		nml=mn.build_namelist_ev_Fsigma(1e4,Flux,Sigma,ES)
		nml['star_job']['years_for_initial_dt']=5e3
		nml.write('inlist', force=True )
		os.system("./create_planet")

		shutil.copyfile(file_out,file_in)
		nml=mn.build_namelist_mp_hold(Thold,L_M_target,Flux,Sigma)
		nml.write('inlist', force=True )
		os.system("./create_planet")

		shutil.copyfile(file_out,file_in)
		nml=mn.build_namelist_ev_Fsigma(1e4,Flux,Sigma,ES)
		nml['star_job']['years_for_initial_dt']=5e3
		nml.write('inlist', force=True )
		os.system("./create_planet")

		file_store="model_%d_L0.mod" %(i)
		shutil.copyfile(file_out,"%s/%s" %(planet_dir,file_store))
		his_out="his_%d_L0.data" %(i)
		shutil.copyfile("LOGS/history.data","%s/%s" %(planet_dir,his_out))

	os.chdir(curr_dir)

def heat_up_model(planet_dir,L_M_want,Teq,i,j):

	# function to heat up a mesa model to a desired L_M

	import os
	import shutil
	import f90nml #fortran namelist tools
	import mesa_namelist as mn # for generating mesa namelist

	file_in="model_in.mod"
	file_out="model_out.mod"
	run_dir="make_planets"
	curr_dir=os.getcwd()

	start_file ="model_%d_L%d.mod" %(i,j-1)
	finish_file="model_%d_L%d.mod" %(i,j)

	#usful constants
	sb=5.6704e-5
	Flux=4.0*sb*Teq**4.0
	Sigma=250 # heating depth in F-Sigma routine
	Thold=1e14 # holding time
	ES=False

	os.chdir(run_dir)

	shutil.copyfile("%s/%s" %(planet_dir,start_file),file_in)

	nml=mn.build_namelist_mp_hold(Thold,L_M_want,Flux,Sigma)
	nml.write('inlist', force=True )
	os.system("./create_planet")

	#evolve forward for short 1e4 year time-scale to get stable model
	shutil.copyfile(file_out,file_in)
	nml=mn.build_namelist_ev_Fsigma(1e4,Flux,Sigma,ES)
	nml['star_job']['years_for_initial_dt']=5e3
	nml.write('inlist', force=True )
	os.system("./create_planet")

	shutil.copyfile(file_out,"%s/%s" %(planet_dir,finish_file))
	his_out="his_%d_L%d.data" %(i,j)
	shutil.copyfile("LOGS/history.data","%s/%s" %(planet_dir,his_out))
	os.chdir(curr_dir)

def build_models(Mcore,rho_core,L_M_whittle,L_M_start,L_M_end,NL,Teq,Xstart,Xend,Nenv,planet_dir,Xbondi,tcool_min,do_not_skip_L0):
	# this function builds many models with envolope mass fractions and luminosities
	#create initial population

	import os
	import numpy as np
	import math
	from nugridpy import mesa as ms
	import shutil
	import planet_prop as pp

	curr_dir=os.getcwd()


	cs=(1.38e-16*Teq/(2.35*1.67e-24))**0.5

	if (do_not_skip_L0):
		build_starting_models(Mcore,rho_core,L_M_whittle,L_M_start,Teq,Xstart,Xend,Nenv,planet_dir)

	L_M_grid=np.logspace(math.log10(L_M_start),math.log10(L_M_end),NL)

	for i in range(0,Nenv):
		for j in range(1,NL):
			print("***123", i, j)
			#read in history file and check planet radius and KH timescale
			his_file="his_%d_L%d.data" %(i,j-1)
			shutil.copyfile("%s/%s" %(planet_dir,his_file),"history.data")
			his=ms.history_data('.',clean_starlog=True)
			rad=10**(his.get('log_R'))*pp.rsun
			mass=his.get('star_mass')*pp.msun
			lum=his.get('luminosity')*pp.lsun
			khtime=his.get('kh_timescale') # in years
			# now adjust khtime
			cooling_time=khtime[-1]*lum[-1]/(L_M_grid[j-1]*mass[-1])
			Rbondi=6.67e-8*mass[-1]/(2.0*cs**2.0)
			XB=rad[-1]/Rbondi
			print(cooling_time/1e6,"Myr",XB)
			if XB > Xbondi:
				break
			if cooling_time < tcool_min:
				break
			# else heat-up to next value
			L_M_curr=(np.float64(L_M_grid[j]).item())
			heat_up_model(planet_dir,L_M_curr,Teq,i,j)
	os.chdir(curr_dir)


def run_models(Teq,sep,NL,Nenv,planet_dir,min_mass,Xbondi,Tevolve,Evolve_star,tcool_min,tcool_max,L_M_start,L_M_end):

	import os
	import numpy as np
	import math
	from nugridpy import mesa as ms
	import shutil
	import planet_prop as pp
	import mesa_namelist as mn

	curr_dir=os.getcwd()
	file_in="model_in.mod"
	file_out="model_out.mod"
	run_dir="make_planets"

	sb=5.6704e-5
	Flux=4.0*sb*Teq**4.0
	Sigma=250 # heating depth in F-Sigma routine

	L_M_grid=np.logspace(math.log10(L_M_start),math.log10(L_M_end),NL)

	cs=(1.38e-16*Teq/(2.35*1.67e-24))**0.5
	for i in range(0,Nenv):
		for j in range(0,NL):
			#read in history file and check planet radius
			his_file="his_%d_L%d.data" %(i,j)
			shutil.copyfile("%s/%s" %(planet_dir,his_file),"history.data")
			his=ms.history_data('.',clean_starlog=True)
			rad=10**(his.get('log_R'))*pp.rsun
			mass=his.get('star_mass')*pp.msun
			lum=his.get('luminosity')*pp.lsun
			khtime=his.get('kh_timescale') # in years
			# now adjust khtime
			cooling_time=khtime[-1]*lum[-1]/(L_M_grid[j]*mass[-1])
			Rbondi=6.67e-8*mass[-1]/(2.0*cs**2.0)
			XB=rad[-1]/Rbondi
			print(i,j, XB, cooling_time/1.e6, "Myr")
			if XB > Xbondi:
				break
			if cooling_time < tcool_min:
				break
			if cooling_time < tcool_max:
				# can run model
				os.chdir(run_dir)
				mod_file="model_%d_L%d.mod" %(i,j)
				shutil.copyfile("%s/%s" %(planet_dir,mod_file),file_in)
				nml=mn.build_namelist_ev_Fsigma_we(Tevolve,Flux,Sigma,Teq,sep,min_mass,Evolve_star)
				nml.write('inlist', force=True )
				os.system("./create_planet")
				his_out="e_his_%d_L%d.data" %(i,j)
				shutil.copyfile("LOGS/history.data","%s/%s" %(planet_dir,his_out))
			os.chdir(curr_dir)
