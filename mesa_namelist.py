def build_namelist_pc(rho_core,Mcore):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'


	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = False
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['relax_core'] = True
	sj['new_core_mass'] = Mcore
	sj['core_avg_rho'] = rho_core
	sj['dlg_core_mass_per_step'] = 0.05
	sj['relax_core_years_for_dt'] = 0.1
	sj['core_avg_eps'] = 0

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 1
	cn['terminal_cnt'] = 10
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 0.01
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = 1e3
	cn['max_model_number']=2000

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;

def build_namelist_mp_hold(Thold,L_M_ratio,Flux,Sigma):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=100 # sep AU
	x_ctrl[12]=1 # stellar mass Msun
	x_ctrl[21]=L_M_ratio # luminisoity to core mass ratio

	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[5]=True # core/envelope heating

	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'

	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = False
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 0.01
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Thold
	cn['use_other_energy']=True
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl
	cn['max_model_number']=2000

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;

def build_namelist_mp_whittle(Thold,L_M_ratio,new_M,mdot,Flux,Sigma):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=100 # sep AU
	x_ctrl[12]=1 # stellar mass Msun
	x_ctrl[21]=L_M_ratio # luminisoity to core mass ratio

	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[5]=True # core/envelope heating

	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'

	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = False
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 0.01
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Thold
	cn['RGB_wind_scheme']=''
	cn['AGB_wind_scheme']=''
	cn['mass_change']=mdot
	cn['star_mass_min_limit']=new_M
	cn['use_other_energy']=True
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl
	cn['max_model_number']=2000

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;

def build_namelist_mp_Fsigma_hold(Thold,L_M_ratio,Flux,Sigma):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=100 # sep AU
	x_ctrl[12]=1 # stellar mass Msun
	x_ctrl[21]=L_M_ratio # luminisoity to core mass ratio

	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[5]=True # core/envelope heating

	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'

	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = False
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_initial_dt']=True
	sj['years_for_initial_dt']=1e8
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 0.01
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Thold
	cn['use_other_energy']=True
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl
	cn['max_model_number']=2000

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;

def build_namelist_ev_Fsigma(Tmax,Flux,Sigma,ES):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=100 # sep AU
	x_ctrl[12]=1 # stellar mass Msun

	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[5]=False # core/envelope heating
	x_logical_ctrl[7]=ES # evolve star

	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'

	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['set_initial_dt']=True
	sj['years_for_initial_dt']=1e5
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = True
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 1e-4
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Tmax
	cn['max_years_for_timestep']=1e8
	cn['use_other_energy']=True
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;

def build_namelist_ev_Fsigma_we(Tmax,Flux,Sigma,Temp,sep,stop_mass,ES):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=sep # sep AU
	x_ctrl[12]=1 # stellar mass Msun
	x_ctrl[18]=Temp
	x_ctrl[23]=6.9183e7 #saturation for X-rays
	x_ctrl[24]=10**-3.6
	x_ctrl[25]=1.19


	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[5]=False # core/envelope heating
	x_logical_ctrl[6]=True #evaporation
	x_logical_ctrl[7]=ES # include evolution of star or not


	pg['read_extra_pgstar_inlist1']=True
	pg['extra_pgstar_inlist1_name']='inlist_pgstar'

	sj['read_extra_star_job_inlist1']=True
	sj['extra_star_job_inlist1_name']='inlist_cache'
	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 3.e6
	sj['set_initial_dt']=True
	sj['years_for_initial_dt']=1e5
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = True
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 3e-4
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Tmax
	cn['max_years_for_timestep']=5e7
	cn['use_other_energy']=True
	cn['use_other_wind']=True
	cn['star_mass_min_limit']=stop_mass
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl
	cn['min_q_for_k_below_const_q'] = 0.9
	cn['min_q_for_k_const_mass'] = 0.5

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;
