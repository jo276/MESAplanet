! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module evolve

      use star_private_def
      use const_def
      use star_utils

      implicit none

      contains
      

      integer function do_evolve_step(id, first_try, just_did_backup)
         use num_def
         use do_one_utils, only: write_terminal_header
         use report, only: do_report
         use winds, only: set_mdot
         use adjust_mass, only: do_adjust_mass         
         use utils_lib, only: is_bad_num, has_bad_num
         use element_diffusion, only: do_element_diffusion
         use alloc, only: check_sizes
         use evolve_support, only: set_current_to_old
         use solve_hydro, only: set_L_burn_by_category
         use struct_burn_mix, only: do_struct_burn_mix
         use hydro_vars, only: set_vars, set_final_vars
         use hydro_mtx, only: dump_struct
         use star_utils, only: start_time, update_time
         use solve_omega_mix, only: do_solve_omega_mix
         use mix_info, only: set_mixing_info
         use hydro_rotation, only: set_rotation_info
         use profile
         
         logical, intent(in) :: first_try, just_did_backup
         integer, intent(in) :: id

         type (star_info), pointer :: s            
         integer :: ierr, time0, clock_rate, &
            j, k, j_cnt, mdot_redo_cnt, max_mdot_redo_cnt
         logical :: okay, trace, skip_global_corr_coeff_limit, &
            have_too_large_wind_mdot, have_too_small_wind_mdot, &
            ignored_first_step, was_in_implicit_wind_limit
         real(dp) :: J_tot1, J_tot2, &
            w_div_w_crit, w_div_w_crit_prev, mstar_dot, mstar_dot_prev, abs_mstar_delta, &
            explicit_mdot, max_wind_mdot, wind_mdot, r_phot, kh_timescale, dmskhf, dmsfac, &
            too_large_wind_mdot, too_small_wind_mdot, boost, mstar_dot_nxt, &
            surf_w_div_w_crit_limit, dt, total
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         do_evolve_step = terminate
         
         ierr = 0         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return         
         
         !call dump_struct(s)
         !stop 'evolve'

         call system_clock(s% system_clock_at_start_of_step, clock_rate)
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'evolve starting xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
            do k=1,s% nz
               write(*,3) 'evolve starting j_rot(k)', &
                  s% model_number, k, s% j_rot(k)
            end do
         end if
         
         mstar_dot = 0
         w_div_w_crit = -1
         surf_w_div_w_crit_limit = s% surf_w_div_w_crit_limit
         mdot_redo_cnt = 0
         max_mdot_redo_cnt = s% max_mdot_redo_cnt
         
         max_wind_mdot = 10*Msun/secyer
         have_too_large_wind_mdot = .false.
         have_too_small_wind_mdot = .false.
         too_large_wind_mdot = 0
         too_small_wind_mdot = 0

         trace = s% trace_evolve
         s% just_did_backup = just_did_backup
         s% doing_newton_iterations = .false.
         s% have_mixing_info = .false.
         ignored_first_step = .false.
         
         if (dbg) then
            call check_sizes(s, ierr)
            if (ierr /= 0) then
               write(*,*) 'do_evolve_step: check_sizes returned ierr', ierr
               return
            end if
         end if

         if (s% doing_first_model_of_run) then
            if (s% do_history_file) then
               if (first_try) then
                  call write_terminal_header(s)
               else
                  write(*,1) '1st model retry log10(dt/yr)', log10_cr(s% dt/secyer)
               end if
            end if
            call system_clock(time0,clock_rate)
            s% starting_system_clock_time = time0
            s% system_clock_rate = clock_rate
            s% initial_timestep = s% dt_next
            s% last_backup = -s% backup_hold
            s% timestep_hold = -111
            if (first_try) s% model_number_old = s% model_number
         end if

         call debug('before prepare_for_new_step')
         
         if (first_try) then
            if (trace) write(*,'(/,a,i8)') 'call prepare_for_new_step', s% model_number
            do_evolve_step = prepare_for_new_step(s)
            if (do_evolve_step /= keep_going) return
         end if

         call debug('before prepare_for_new_try')

         if (trace) write(*,'(/,a,i8)') 'call prepare_for_new_try', s% model_number
         do_evolve_step = prepare_for_new_try(s)
         if (do_evolve_step /= keep_going) return

         call debug('before set_mdot')

         if (trace) write(*,'(/,a,i8)') 'call set_mdot', s% model_number
         
         ! set mdot for the step
         call set_mdot(s, s% L_phot*Lsun, s% mstar, s% Teff, ierr)
         if (ierr /= 0) then
            do_evolve_step = retry
            s% result_reason = nonzero_ierr
            if (s% report_ierr) write(*, *) 'do_evolve_step: set_mdot ierr'
            return
         end if
         
         if (s% use_other_adjust_mdot) then
            call s% other_adjust_mdot(s% id, ierr)
            if (ierr /= 0) then
               do_evolve_step = retry
               s% result_reason = nonzero_ierr
               if (s% report_ierr) write(*, *) 'do_evolve_step: other_adjust_mdot ierr'
               return
            end if
         end if

         explicit_mdot = s% mstar_dot
         
         was_in_implicit_wind_limit = s% was_in_implicit_wind_limit
         if (was_in_implicit_wind_limit .and. &
             s% generations >= 2 .and. &
             abs(s% mstar_dot_old) > 0 .and. &
             abs((s% mstar_dot-s% mstar_dot_old)/s% mstar_dot_old)+1 > s% mdot_revise_factor) then
             write(*,*) "Skipping first step in implicit mdot"
             s% mstar_dot = s% mstar_dot_old
             mdot_redo_cnt = 1
             ignored_first_step = .true.
         end if
         
         if (s% dt <= 0) then
            do_evolve_step = terminate
            s% termination_code = t_dt_is_zero
            s% result_reason = dt_is_zero
            return
         end if
      
         abs_mstar_delta = 0 
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'evolve after remesh xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
         end if
         
      implicit_mdot_loop: do
         
            s% time = s% time_old + s% dt

            call debug('before do_adjust_mass')

            if (trace) write(*,'(/,a)') 'call do_adjust_mass'         
            if (s% doing_timing) call start_time(s, time0, total)
            call do_adjust_mass(s, s% species, ierr)
            if (s% doing_timing) call update_time(s, time0, total, s% time_adjust_mass)         
            if (ierr /= 0) then
               do_evolve_step = retry
               s% result_reason = adjust_mass_failed
               if (s% report_ierr) write(*, *) 'do_evolve_step: do_adjust_mass ierr'
               return
            end if

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'evolve after do_adjust_mass xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            dt = s% dt

            call debug('before do_set_vars')
            
            if (.not. s% do_element_diffusion) then

               call do_set_vars(0, dt, ierr)
               if (ierr /= 0) return

            else
            
               call do_set_vars(1, dt, ierr)
               if (ierr /= 0) return
               
               call debug('before do_element_diffusion')
               if (trace) write(*,'(/,a)') 'call do_element_diffusion'
               if (s% doing_timing) call start_time(s, time0, total)
               okay = do_element_diffusion(s, s% dt)
               if (s% doing_timing) call update_time(s, time0, total, s% time_element_diffusion)         
               if (.not. okay) then
                  if (s% report_ierr) then
                     write(*, *) 'element diffusion failed: retry', s% model_number
                  end if
                  do_evolve_step = retry
                  s% result_reason = diffusion_failed
                  return
               end if
                  
               call do_set_vars(2, dt, ierr)
               if (ierr /= 0) return
               
            end if

            s% need_to_adjust_J_lost = .true.
            if (s% rotation_flag .and. s% premix_omega) then
               do_evolve_step = do_solve_omega_mix(s, 0.5d0*dt)
               if (do_evolve_step /= keep_going) return
               call set_rotation_info(s, ierr)
               if (ierr /= 0) return
               call set_mixing_info(s, ierr)
               if (ierr /= 0) return            
            end if

            call save_start_values(s, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'save_start_values failed'
               return
            end if
            
            call set_Eulerian_Lagrangian_for_eps_grav

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'evolve before do_struct_burn_mix xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            call debug('before do_struct_burn_mix')
            if (trace) write(*,'(/,a,i8)') 'call do_struct_burn_mix', s% model_number            
            skip_global_corr_coeff_limit = &
               (first_try .or. (just_did_backup .and. s% number_of_backups_in_a_row == 1))
            if (s% doing_timing) call start_time(s, time0, total)
            do_evolve_step = do_struct_burn_mix( &
               s, skip_global_corr_coeff_limit, dt)
            if (s% doing_timing) call update_time(s, time0, total, s% time_struct_burn_mix)         
            if (do_evolve_step /= keep_going) return
            call debug('after do_struct_burn_mix')
        
            if (.not. s% rotation_flag) exit implicit_mdot_loop
            if (s% mstar_dot == 0) exit implicit_mdot_loop
            if (max_mdot_redo_cnt <= 0) exit implicit_mdot_loop
                        
            mstar_dot_prev = mstar_dot
            mstar_dot = s% mstar_dot
            wind_mdot = -s% mstar_dot 
               ! for a wind, wind_mdot > 0, but we don't require this in the following.
               ! i.e., the following does not assume that mstar_dot < 0
            
            if (mdot_redo_cnt == 1 .or. ignored_first_step) then ! this is the 1st correction to mdot
               r_phot = sqrt(s% L(1)/(pi*crad*clight*pow4(s% Teff)))
               if (s% M_center >= 0.0) then
                  call new_eval_kh_timescale(s,kh_timescale)
               else
                  kh_timescale = eval_kh_timescale(s% cgrav(1), s% mstar, r_phot, s% L(1))
               endif
               dmskhf = s% rotational_mdot_kh_fac
               dmsfac = s% rotational_mdot_boost_fac
               max_wind_mdot = dmskhf*s% mstar/kh_timescale
               if (wind_mdot > 0) max_wind_mdot = min(max_wind_mdot, wind_mdot*dmsfac)
            end if
            
            w_div_w_crit_prev = w_div_w_crit
            ! check the new w_div_w_crit to make sure not too large
            call set_surf_avg_rotation_info(s)
            w_div_w_crit = s% w_div_w_crit_avg_surf
            
            !write(*,2) 'log wind_mdot', mdot_redo_cnt, log10(abs(wind_mdot)/(Msun/secyer))
            
            if (wind_mdot >= max_wind_mdot) then
               if (mdot_redo_cnt == 0) then
                  write(*,*) 'cannot fix omega >= omega_crit -- wind mass loss already at max'
               else
                  write(*,2) 'retry: at max wind mass loss', s% model_number, &
                     log10_cr(max_wind_mdot/(Msun/secyer))
                  do_evolve_step = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               write(*,*)
               if (w_div_w_crit > surf_w_div_w_crit_limit) then
                  write(*,1) 'retry: w_div_w_crit > surf_w_div_w_crit_limit', &
                     w_div_w_crit, surf_w_div_w_crit_limit
                  do_evolve_step = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               exit implicit_mdot_loop
            end if
               
            ! NOTE: we assume that if surface omega/omega_crit (w_div_w_crit) is too large,
            ! then mass loss needs to be made larger to fix the problem.
            ! if that assumption is wrong, i.e. if bigger mass loss makes w_div_w_crit worse,
            ! then we are in an unstable situation and will remove mass until regain stability.
            
            if (w_div_w_crit <= surf_w_div_w_crit_limit &
                  .and. mdot_redo_cnt == 0) then
               s% was_in_implicit_wind_limit = .false.
               exit implicit_mdot_loop 
            end if
               ! normal case; no problem; no redo required.
            
            if (w_div_w_crit <= surf_w_div_w_crit_limit &
                  .and. s% mstar_dot == explicit_mdot) exit implicit_mdot_loop 
               ! implicit scheme reached the limit setted by the explicit_mdot; no problem; no redo required.

            s% was_in_implicit_wind_limit = .true.
               
            if (s% dt/secyer < s% min_years_dt_for_redo_mdot) then
               if (.true.) write(*,1) &
                  'dt too small for fix to fix w > w_crit; min_years_dt_for_redo_mdot', &
                  s% dt/secyer, s% min_years_dt_for_redo_mdot
               exit implicit_mdot_loop
            end if
            
            ! if get here, need to revise mdot to fix w_div_w_crit
            
            mdot_redo_cnt = mdot_redo_cnt + 1
            
            if (mdot_redo_cnt == 1) then ! this is the 1st correction to mdot

               call set_current_to_old(s)
               do_evolve_step = prepare_for_new_try(s)
               if (do_evolve_step /= keep_going) return
               
               have_too_small_wind_mdot = .true.
               too_small_wind_mdot = wind_mdot
               if (s% mstar_dot < 0) then
                  s% mstar_dot = mstar_dot*s% mdot_revise_factor
               else
                  s% mstar_dot = mstar_dot/s% mdot_revise_factor
               end if
               
               if (-s% mstar_dot > max_wind_mdot) s% mstar_dot = -max_wind_mdot
                  
               write(*,3) 'w > w_crit: revise mdot and redo', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))

               !abs_mstar_delta = max(abs(s% mstar_dot), 1d-6*Msun/secyer)
               abs_mstar_delta = abs(s% mstar_dot)
               
               cycle implicit_mdot_loop
               
            else if (mdot_redo_cnt == 2 .and. ignored_first_step) then
               abs_mstar_delta = abs(s% mstar_dot_old)
            end if
            
            ! have already done at least one correction -- check if okay now
            if (w_div_w_crit <= surf_w_div_w_crit_limit .and. &
                  have_too_small_wind_mdot .and. &
                  abs((wind_mdot-too_small_wind_mdot)/wind_mdot) < s% surf_w_div_w_crit_tol) then
               write(*,3) 'OKAY', s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
               write(*,*)
               exit implicit_mdot_loop ! in bounds so accept it
            end if

            if (mdot_redo_cnt >= max_mdot_redo_cnt) then
               if (max_mdot_redo_cnt > 0) then
                  write(*,3) 'failed to fix w > w_crit: too many tries', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer))
                  do_evolve_step = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               exit implicit_mdot_loop
            end if
            
            if (w_div_w_crit > surf_w_div_w_crit_limit &
                  .and. w_div_w_crit_prev >= surf_w_div_w_crit_limit &
                  .and. -mstar_dot >= max_wind_mdot) then
               write(*,3) 'failed to fix w > w_crit', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
               write(*,*)
               if (.true.) then
                  do_evolve_step = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               exit implicit_mdot_loop ! give up
            end if
            
            if (w_div_w_crit >= surf_w_div_w_crit_limit) then ! wind too small
               !write(*,*) "entering too small wind mdot"
               if (.not. have_too_small_wind_mdot) then
                  !write(*,*) "setting too small wind mdot"
                  too_small_wind_mdot = wind_mdot
                  have_too_small_wind_mdot = .true.
               else if (wind_mdot > too_small_wind_mdot) then
                  !write(*,*) "changing too small wind mdot"
                  too_small_wind_mdot = wind_mdot
               else if (.false.) then ! oops. not converging.
                  write(*,3) 'failed to fix w > w_crit: not converging', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer))
                  write(*,*)
                  if (.true.) then
                     do_evolve_step = retry
                     s% result_reason = nonzero_ierr
                     return
                  end if
                  exit implicit_mdot_loop ! give up
               end if
            else if (w_div_w_crit < surf_w_div_w_crit_limit) then ! wind too large
               !write(*,*) "entering too large wind mdot"
               if (.not. have_too_large_wind_mdot) then
                  !write(*,*) "setting too large wind mdot"
                  too_large_wind_mdot = wind_mdot
                  have_too_large_wind_mdot = .true.
               else if (wind_mdot < too_large_wind_mdot) then
                  !write(*,*) "changing too large wind mdot"
                  too_large_wind_mdot = wind_mdot
               else if (.false.) then ! oops. not converging.
                  write(*,3) 'failed to fix w > w_crit: not converging', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer))
                  write(*,*)
                  if (.true.) then
                     do_evolve_step = retry
                     s% result_reason = nonzero_ierr
                     return
                  end if
                  exit implicit_mdot_loop ! give up
               end if
            end if

            call set_current_to_old(s)
            do_evolve_step = prepare_for_new_try(s)
            if (do_evolve_step /= keep_going) return
            
            if (have_too_large_wind_mdot .and. have_too_small_wind_mdot) then
               if (abs((too_large_wind_mdot-too_small_wind_mdot)/too_large_wind_mdot) &
                   < s% surf_w_div_w_crit_tol) then
                  write(*,*) "too_large_wind_mdot good enough, using it"
                  s% mstar_dot = -too_large_wind_mdot
               else
                  ! have bracketing mdots; bisect for next one.
                  s% mstar_dot = -0.5d0*(too_large_wind_mdot + too_small_wind_mdot)
                  write(*,3) 'fix w > w_crit: bisect mdots and redo', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer)), &
                     log10_cr(abs(too_large_wind_mdot)/(Msun/secyer)), &
                     log10_cr(abs(too_small_wind_mdot)/(Msun/secyer))
               end if
                  
            else ! still have wind too small so boost it again
               if (have_too_small_wind_mdot) then
                  if (mod(mdot_redo_cnt,2) == 1) then
                     boost = s% implicit_mdot_boost
                     ! increase mass loss
                     mstar_dot_nxt = mstar_dot - boost*abs_mstar_delta
                  else
                     if (mstar_dot < 0) then ! increase mass loss
                        mstar_dot_nxt = mstar_dot*s% mdot_revise_factor
                     else ! decrease mass gain
                        mstar_dot_nxt = mstar_dot/s% mdot_revise_factor
                     end if
                  end if
               else
                  if (mod(mdot_redo_cnt,2) == 1) then
                     boost = s% implicit_mdot_boost
                     ! decrease mass loss
                     mstar_dot_nxt = mstar_dot + boost*abs_mstar_delta
                  else
                     if (mstar_dot < 0) then ! decrease mass loss
                        mstar_dot_nxt = mstar_dot/s% mdot_revise_factor
                     else ! increase mass gain
                        mstar_dot_nxt = mstar_dot*s% mdot_revise_factor
                     end if
                  end if
               end if
               if (mstar_dot_prev /= explicit_mdot) &
                  mstar_dot_nxt = min(mstar_dot_nxt, explicit_mdot)
               if (mstar_dot_nxt == explicit_mdot) &
                  write(*,*) "implicit mdot: reached explicit_mdot"
               if (.false.) then
                  write(*,*)
                  write(*,1) 'log mstar_dot_prev', log10_cr(abs(mstar_dot_prev)/(Msun/secyer))
                  write(*,1) 'log mstar_dot', log10_cr(abs(mstar_dot)/(Msun/secyer))
                  write(*,1) 'log mstar_dot_nxt', log10_cr(abs(mstar_dot_nxt)/(Msun/secyer))
                  write(*,1) 'boost*abs_mstar_delta', boost*abs_mstar_delta
                  write(*,1) 'boost', boost
                  write(*,1) 'abs_mstar_delta', abs_mstar_delta
                  write(*,1) 'abs_mstar_delta/mstar_dot', abs_mstar_delta/mstar_dot
                  write(*,1) 'w_div_w_crit', w_div_w_crit
                  write(*,*)
               end if
               s% mstar_dot = mstar_dot_nxt
               if (-s% mstar_dot > max_wind_mdot) s% mstar_dot = -max_wind_mdot
               !abs_mstar_delta = max(abs_mstar_delta, abs(s% mstar_dot))
               write(*,3) 'fix w > w_crit: change mdot and redo', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
            end if

         end do implicit_mdot_loop

         call debug('before set_final_vars')
         if (trace) write(*,'(/,a)') 'call set_final_vars'
         call set_final_vars(s, s% dt, ierr) ! will use these for next step
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_evolve_step: set_final_vars ierr'
            do_evolve_step = retry
            s% result_reason = nonzero_ierr
            return
         end if
         
         if (s% max_timestep_hi_T_limit > 0 .and. &
               s% max_years_for_timestep /= s% hi_T_max_years_for_timestep) then
            if (maxval(s% T(1:s% nz)) >= s% max_timestep_hi_T_limit) then
               write(*,1) 'switch to high T max timesteps'
               if (.false.) then
                  write(*,1) 's% max_timestep_hi_T_limit', s% max_timestep_hi_T_limit
                  write(*,1) 's% max_timestep', s% max_timestep
                  write(*,1) 's% max_years_for_timestep', s% max_years_for_timestep
                  write(*,1) 's% hi_T_max_years_for_timestep', s% hi_T_max_years_for_timestep
                  write(*,1) 'maxval(s% T(1:s% nz))', maxval(s% T(1:s% nz))
                  write(*,*) 's% max_years_for_timestep /= s% hi_T_max_years_for_timestep', &
                     s% max_years_for_timestep /= s% hi_T_max_years_for_timestep
                  !stop
               end if
               s% max_years_for_timestep = s% hi_T_max_years_for_timestep
               s% max_timestep = secyer*s% max_years_for_timestep
            end if
         end if
         
         call debug('before do_report')

         if (trace) write(*,'(/,a)') 'call do_report'
         call system_clock(time0,clock_rate)
         s% current_system_clock_time = time0
         s% total_elapsed_time = dble(time0 - s% starting_system_clock_time)/dble(clock_rate)
         call do_report(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_evolve_step: do_report ierr'
            do_evolve_step = retry
            s% result_reason = nonzero_ierr
            return
         end if
         
         call set_L_burn_by_category(s) ! final values for use in selecting timestep
         
         s% total_angular_momentum = total_angular_momentum(s)
         !write(*,2) 'total_angular_momentum after do_evolve_step', &
         !   s% model_number, s% total_angular_momentum
         
         if (trace) write(*,'(/,a)') 'done do_evolve_step'
         
         call debug('done do_evolve_step')
         
         
         contains
         
         
         subroutine debug(str) 
            use chem_def
            character (len=*), intent(in) :: str
            integer :: k, j
            include 'formats'
            
            !write(*,2) trim(str) // ' s% xa(1,2)', s% model_number, s% xa(1,2)
            return
            
            if (.not. s% rotation_flag) return
            k = 1
            write(*,3) trim(str) // ' s% omega(k)', k, s% model_number, s% omega(k)
            return
            j = 2
            !do j=1,1 !s% species
               if (.true. .or. s% xa(j,k) > 1d-9) &
                  write(*,1) trim(str) // ' xin(net_iso(i' // &
                     trim(chem_isos% name(s% chem_id(j))) // '))= ', &
                     s% xa(j,k), s% abar(k)
            !end do
         end subroutine debug
         
         
         subroutine do_set_vars(which_case, dt, ierr)
            use hydro_vars, only: set_vars, set_vars_before_diffusion, reset_vars_after_diffusion
            integer, intent(in) :: which_case
            real(dp), intent(in) :: dt
            integer, intent(out) :: ierr
            integer :: nz
            nz = s% nz
            ierr = 0
            if (trace) write(*,'(/,a)') 'call set_vars'
            select case(which_case)
            case (0)
               call set_vars(s, dt, ierr)
            case (1)
               call set_vars_before_diffusion(s, dt, ierr)
            case (2)
               call reset_vars_after_diffusion(s, dt, ierr)
            end select
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*, *) 'do_evolve_step: set_vars ierr: retry', s% model_number
               end if
               do_evolve_step = retry
               s% result_reason = nonzero_ierr
               return
            end if
         end subroutine do_set_vars
         
         
         subroutine set_Eulerian_Lagrangian_for_eps_grav
            use star_utils, only: set_k_CpTMdot_lt_L
            real(dp) :: dxm_CpTMdot_lt_L, dxm_kA, dxm_target
            integer :: k, kA, kB, nz, transition
            
            include 'formats'
            
            if (s% mstar_dot == 0d0) then
               s% k_below_Eulerian_eps_grav = 1
               s% k_Lagrangian_eps_grav = 1
               return
            end if
            
            nz = s% nz
            kA = 1
            dxm_kA = 0
            
            if (s% min_dxm_Eulerian_div_dxm_CpTMdot_lt_L > 0) then
               call set_k_CpTMdot_lt_L(s)
               dxm_CpTMdot_lt_L = sum(s% dm(1:s% k_CpTMdot_lt_L))
               dxm_target = s% min_dxm_Eulerian_div_dxm_CpTMdot_lt_L*dxm_CpTMdot_lt_L
               do k = kA, nz
                  if (dxm_kA >= dxm_target) exit
                  kA = k
                  dxm_kA = dxm_kA + s% dm(k)
               end do
            end if
            
            if (s% mstar_dot > 0 .and. s% min_dxm_Eulerian_div_dxm_added > 0) then
               dxm_target = (s% xmstar - s% xmstar_old)*s% min_dxm_Eulerian_div_dxm_added
               do k = kA, nz
                  if (dxm_kA >= dxm_target) exit
                  kA = k
                  dxm_kA = dxm_kA + s% dm(k)
               end do
            end if
            
            if (s% mstar_dot < 0 .and. s% min_dxm_Eulerian_div_dxm_removed > 0) then
               dxm_target = (s% xmstar_old - s% xmstar)*s% min_dxm_Eulerian_div_dxm_removed
               do k = kA, nz
                  if (dxm_kA >= dxm_target) exit
                  kA = k
                  dxm_kA = dxm_kA + s% dm(k)
               end do
            end if
            
            if (kA == 1) then ! pure Lagrangian
               kB = 1
            else if (kA >= nz) then ! pure Eulerian
               kA = nz+1
               kB = kA
            else ! transition zone between
               transition = max(0, s% min_cells_for_Eulerian_to_Lagrangian_transition)
               kB = kA
               kA = kB-transition
               if (kA < 1) then
                  kA = 1
                  kB = kA+transition
               end if
            end if

            s% k_below_Eulerian_eps_grav = kA ! pure Eulerian for k < this
            s% k_Lagrangian_eps_grav = kB ! pure Lagrangian for k >= this
      
         end subroutine set_Eulerian_Lagrangian_for_eps_grav  


      end function do_evolve_step
                     

      subroutine save_start_values(s, ierr)
         use solve_hydro, only: set_L_burn_by_category
         use chem_def, only: num_categories
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, j
         include 'formats'    
         ierr = 0         
         call set_L_burn_by_category(s)
         do k=1,s% nz
            s% lnd_start(k) = s% lnd(k)
            s% lnP_start(k) = s% lnP(k)
            s% lnPgas_start(k) = s% lnPgas(k)
            s% lnT_start(k) = s% lnT(k)
            s% lnR_start(k) = s% lnR(k)
            s% L_start(k) = s% L(k)
            s% omega_start(k) = s% omega(k)
            s% ye_start(k) = s% ye(k)
            s% Z_start(k) = min(1d0, max(0d0, 1d0 - (s% X(k) + s% Y(k))))
            s% i_rot_start(k) = s% i_rot(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% non_nuc_neu_start(k) = s% non_nuc_neu(k)
            s% mass_correction_start(k) = s% mass_correction(k)
            s% P_div_rho_start(k) = s% P(k)/s% rho(k)
            do j=1,s% species
               s% dxdt_nuc_start(j,k) = s% dxdt_nuc(j,k)
            end do
            do j=1,num_categories
               s% luminosity_by_category_start(j,k) = &
                  s% luminosity_by_category(j,k)
            end do
            s% grada_start(k) = s% grada(k)
            s% gradr_start(k) = s% gradr(k)
            s% grada_at_face_start(k) = s% grada_at_face(k)
            s% chiT_start(k) = s% chiT(k)
            s% chiRho_start(k) = s% chiRho(k)
            s% cp_start(k) = s% cp(k)
            s% cv_start(k) = s% cv(k)
            s% gam_start(k) = s% gam(k)
            s% eta_start(k) = s% eta(k)
            s% T_start(k) = s% T(k)
            s% abar_start(k) = s% abar(k)
            s% zbar_start(k) = s% zbar(k)
            s% mu_start(k) = s% mu(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% csound_start(k) = s% csound(k)
            s% burn_num_iters(k) = 0
            s% burn_num_substeps(k) = 0
         end do
         s% have_start_values = .true.
      end subroutine save_start_values

      
      integer function prepare_for_new_step(s)
         use utils_lib, only: has_bad_num, is_bad_num
         use evolve_support, only: new_generation
         use chem_def

         type (star_info), pointer :: s
         
         integer :: ierr, nz
         logical :: trace
         
         include 'formats'

         ierr = 0
         trace = s% trace_evolve
         nz = s% nz

         prepare_for_new_step = keep_going

         if (s% dt_next <= 0) then
            prepare_for_new_step = terminate
            if (s% time >= s% max_age*secyer) then
               s% result_reason = result_reason_normal
               s% termination_code = t_max_age
            else
               s% result_reason = dt_is_zero
               s% termination_code = t_dt_is_zero
            end if
            return
         end if

         if (s% dt_next < s% min_timestep_limit) then
            write(*, *) 's% dt_next', s% dt_next
            write(*, *) 's% min_timestep_limit', s% min_timestep_limit
            prepare_for_new_step = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if
            
         if (trace) write(*,*) 'call set_start_of_step_info'
         call set_start_of_step_info(s, ierr)
         if (failed('set_start_of_step_info ierr')) return
                  
         if (.not. s% doing_first_model_of_run) then
            if (trace) write(*,*) 'call do_mesh'
            prepare_for_new_step = do_mesh(s)
            if (prepare_for_new_step /= keep_going) then
               write(*,*) 'failed in do_mesh'
               prepare_for_new_step = terminate
               return
            end if
         end if

         if (trace) write(*,*) 'call new_generation'
         call new_generation(s, ierr)         
         if (failed('new_generation ierr')) return

         s% dt = s% dt_next 
         s% dt_start = s% dt         
         s% retry_cnt = 0
         s% generations = min(max_generations, s% generations+1)
         
         s% need_to_save_profiles_now = .false.
         s% need_to_update_history_now = s% doing_first_model_of_run
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            prepare_for_new_step = terminate
            if (s% report_ierr) write(*, *) 'prepare_for_new_step: ' // trim(str)
            s% result_reason = nonzero_ierr
         end function failed


      end function prepare_for_new_step


      integer function do_mesh(s)
         use adjust_mesh, only: remesh
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s      
         logical, parameter :: okay_to_merge = .true.
         integer :: time0, clock_rate
         real(dp) :: total
         include 'formats'
         if (s% doing_timing) call start_time(s, time0, total)
         do_mesh = remesh(s, okay_to_merge)
         if (s% doing_timing) call update_time(s, time0, total, s% time_remesh)         
         if (do_mesh /= keep_going) then
            s% result_reason = adjust_mesh_failed
            if (s% report_ierr) write(*, *) 'do_mesh: remesh failed'
            return
         end if
         s% have_start_values = .false.
      end function do_mesh
      

      integer function prepare_for_new_try(s) 
         ! return keep_going, terminate, retry, or backup
         ! if don't return keep_going, then set result_reason to say why.
         use utils_lib, only: is_bad_num, has_bad_num
         use net_lib, only: clean_up_fractions
         use net, only: get_screening_mode
         use star_utils, only: save_for_d_dt

         type (star_info), pointer :: s
         
         integer :: ierr, i, j, k, nz, nvar, nvar_hydro
         real(dp), parameter :: max_sum_abs = 10d0, xsum_tol = 1d-2
         real(dp) :: r00, r003, rp13, rm13, r_in, r_out, screening 
         logical :: okay

         include 'formats'
         
         ierr = 0
         nvar = s% nvar
         nvar_hydro = s% nvar_hydro
         nz = s% nz

         s% result_reason = result_reason_normal
         s% model_number = s% model_number_old + 1
         s% termination_code = 0
         s% newton_iter = 0
         s% newton_adjust_iter = 0
         
         screening = get_screening_mode(s,ierr)
         if (ierr /= 0) then
            write(*,*) 'bad value for screening_mode ' // trim(s% screening_mode)
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         
         ! check dimensions
         if (size(s% xh_old,dim=1) /= nvar_hydro .or. size(s% xh_old,dim=2) < nz) then
            write(*,*) 'bad dimensions for xh_old', size(s% xh_old,dim=1), nvar_hydro, &
               size(s% xh_old,dim=2), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% xa_old,dim=1) /= s% species .or. size(s% xa_old,dim=2) < nz) then
            write(*,*) 'bad dimensions for xa_old', size(s% xa_old,dim=1), s% species, &
               size(s% xa_old,dim=2), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% q_old,dim=1) < nz) then
            write(*,*) 'bad dimensions for q_old', size(s% q_old,dim=1), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% dq_old,dim=1) < nz) then
            write(*,*) 'bad dimensions for dq_old', size(s% dq_old,dim=1), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         
         do k = 1, nz
            do j=1,nvar_hydro
               s% xh(j,k) = s% xh_old(j,k) ! start from copy of old structure
            end do
            do j=1,s% species
               s% xa(j,k) = s% xa_old(j,k) ! start from copy of old composition
            end do
            s% q(k) = s% q_old(k) ! start with same q's
            s% dq(k) = s% dq_old(k) ! start with same dq's
         end do

         call set_m_and_dm(s)
         call set_dm_bar(nz, s% dm, s% dm_bar)
         
         if (s% rotation_flag) then
            okay = .true.
            do k=1,nz
               s% omega(k) = s% omega_old(k)
               if (is_bad_num(s% omega(k)) .or. &
                     s% omega(k) > 1d50 .or. &
                     (s% omega_old(nz) /= 0 .and. s% omega(k) < 1d-50)) then
                  write(*,2) 's% omega(k)', k, s% omega(k)
                  okay = .false.
               end if
            end do
            if (.not. okay) then
               write(*,2) 'model_number', s% model_number
               stop 'prepare_for_new_try'
            end if
            call use_xh_to_update_i_rot_and_j_rot(s)
            s% total_angular_momentum = total_angular_momentum(s)
            !write(*,2) 'total_angular_momentum after use_xh_to_update_i_rot_and_j_rot', &
            !   s% model_number, s% total_angular_momentum
            if (s% total_angular_momentum < 0) then
               write(*,*) 'ERROR: doing rotation, but total_angular_momentum < 0'
               prepare_for_new_try = terminate
               s% termination_code = t_negative_total_angular_momentum
               return
               stop 'prepare_for_new_try'
            end if
         end if
         
         if (s% just_did_backup) then
            !write(*,*) 'just_did_backup call save_for_d_dt'
            call save_for_d_dt(s)
            !write(*,*) 'just_did_backup call set_start_of_step_info'
            call set_start_of_step_info(s, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*, *) 'prepare_for_new_try: set_start_of_step_info ierr'
               s% result_reason = nonzero_ierr
               prepare_for_new_try = terminate; return
            end if
         end if
         
         prepare_for_new_try = keep_going
         
         
         !write(*,2) 'done prepare_for_new_try omega(1)', s% model_number, s% omega(1)
         
         
      end function prepare_for_new_try
      
      
      integer function pick_next_timestep(id)
         ! determine what we want for the next timestep
         ! if don't return keep_going, then set result_reason to say why.
         use timestep, only: timestep_controller
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: i, j, n
         real(dp) :: max_timestep, remaining_years, min_max, prev_max_years
         include 'formats'
         
         pick_next_timestep = terminate         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% max_years_for_timestep > 0) then
            max_timestep = secyer*s% max_years_for_timestep
         else
            max_timestep = s% max_timestep
         end if
         pick_next_timestep = timestep_controller(s, max_timestep)
         if (pick_next_timestep /= keep_going) then
            if (s% trace_evolve) &
               write(*,*) 'pick_next_timestep: timestep_controller /= keep_going'
            return
         end if
         s% dt_next_unclipped = s% dt_next 
               ! write out the unclipped timestep in saved models
         
         if ((s% time + s% dt_next) > s% max_age*secyer) then
            s% dt_next = max(0d0, s% max_age*secyer - s% time)
         else if (s% num_adjusted_dt_steps_before_max_age > 0 .and. &
                  s% max_years_for_timestep > 0) then
            remaining_years = s% max_age - s% star_age
            if (s% using_revised_max_yr_dt) &
               s% max_years_for_timestep = s% revised_max_yr_dt
            n = floor(remaining_years/s% max_years_for_timestep + 1d-6)
            j = s% num_adjusted_dt_steps_before_max_age
            if (remaining_years <= s% max_years_for_timestep) then
               s% max_years_for_timestep = remaining_years
               s% using_revised_max_yr_dt = .true.
               s% revised_max_yr_dt = s% max_years_for_timestep
               s% dt_next = s% max_years_for_timestep*secyer
               write(*,3) 'remaining steps and years until max age', &
                  s% model_number, 1, remaining_years
            else if (n <= j) then
               prev_max_years = s% max_years_for_timestep
               i = floor(remaining_years/s% dt_years_for_steps_before_max_age + 1d-6)
               if ((i+1d-9)*s% dt_years_for_steps_before_max_age < remaining_years) then
                  s% max_years_for_timestep = remaining_years/(i+1)
               else
                  s% max_years_for_timestep = remaining_years/i
               end if
               min_max = prev_max_years*s% reduction_factor_for_max_timestep
               if (s% max_years_for_timestep < min_max) &
                  s% max_years_for_timestep = min_max
               if (.not. s% using_revised_max_yr_dt) then
                  s% using_revised_max_yr_dt = .true.
                  write(*,2) 'begin reducing max timestep prior to max age', &
                     s% model_number, remaining_years
               else if (s% revised_max_yr_dt > s% max_years_for_timestep) then
                  write(*,2) 'reducing max timestep prior to max age', &
                     s% model_number, remaining_years
               else if (s% max_years_for_timestep <= s% dt_years_for_steps_before_max_age) then
                  i = floor(remaining_years/s% max_years_for_timestep + 1d-6)
                  write(*,3) 'remaining steps and years until max age', &
                     s% model_number, i, remaining_years
               else 
                  write(*,2) 'remaining_years until max age', &
                     s% model_number, remaining_years
               end if
               s% revised_max_yr_dt = s% max_years_for_timestep
               if (s% dt_next/secyer > s% max_years_for_timestep) &
                  s% dt_next = s% max_years_for_timestep*secyer
            end if
            
         end if
         
      end function pick_next_timestep
      
      
      integer function prepare_to_redo(id)
         use evolve_support, only: set_current_to_old
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr
         include 'formats'        
         ierr = 0 
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            prepare_to_redo = terminate
            return
         end if         
         prepare_to_redo = keep_going         
         if (s% trace_evolve) write(*,'(/,a)') 'prepare_to_redo'        
         call set_current_to_old(s)             
      end function prepare_to_redo
      
      
      integer function prepare_to_retry(id)
         use evolve_support, only: set_current_to_old
         integer, intent(in) :: id
         
         real(dp) :: retry_factor
         type (star_info), pointer :: s
         integer :: ierr, k
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            prepare_to_retry = terminate
            return
         end if
         
         if (s% trace_evolve) write(*,'(/,a)') 'prepare_to_retry'
                  
         s% retry_cnt = s% retry_cnt + 1
         if (s% retry_limit > 0 .and. s% retry_cnt > s% retry_limit) then
            s% dt_start = sqrt(s% dt*s% dt_start)
            prepare_to_retry = backup
            write(*,*) 'have reached retry limit so now backup'
            !stop
            return
         end if

         prepare_to_retry = keep_going
         
         retry_factor = s% timestep_factor_for_retries
         s% dt = s% dt*retry_factor
         write(*,'(a50,2i6,2f16.6)') 'retry cnt, step, log10(dt/yr), retry_factor', &
            s% retry_cnt, s% model_number, &
            log10_cr(s% dt*retry_factor/secyer), retry_factor
         if (s% dt <= max(s% min_timestep_limit,0d0)) then
            write(*,1) 'dt', s% dt
            write(*,1) 'min_timestep_limit', s% min_timestep_limit
            call report_convergence_problems(s, 'dt < min_timestep_limit')
            prepare_to_retry = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if

         if (s% max_years_for_timestep > 0) &
            s% max_timestep = secyer*s% max_years_for_timestep
         if (s% max_timestep > 0) s% dt = min(s% dt, s% max_timestep)
         
         call set_current_to_old(s)

         s% num_retries = s% num_retries+1
         !write(*,2) 'prepare_to_retry s% num_retries', s% num_retries
         !write(*,2) 's% max_number_retries', s% max_number_retries
         if (s% num_retries > s% max_number_retries .and. s% max_number_retries >= 0) then
            write(*,2) 'num_retries', s% num_retries
            write(*,2) 'max_number_retries', s% max_number_retries
            call report_convergence_problems(s, '-- too many retries')
            s% termination_code = t_max_number_retries
            prepare_to_retry = terminate; return
         end if

         s% model_number_for_last_retry = s% model_number
            
      end function prepare_to_retry
      
      
      subroutine report_convergence_problems(s,str)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: str
         write(*,*)
         write(*,*) 'stopping because of convergence problems ' // trim(str)
         write(*,*)
      end subroutine report_convergence_problems
      
      
      integer function do1_backup(id)
         ! returns either keep_going or terminate
         ! at end of this routine, must have same vars set as at end of prepare_for_new_step.
         
         use evolve_support, only: restore_older, set_current_to_old
         use alloc, only: free_star_info_arrays, allocate_star_info_arrays, &
            set_var_info, set_chem_names
         
         integer, intent(in) :: id
         
         integer :: ierr
         type (star_info), pointer :: s
         
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if
         
         if (s% trace_evolve) write(*,'(/,a)') 'do1_backup'
         
         do1_backup = keep_going
         
         s% dt = s% dt_start * &
            (s% timestep_factor_for_backups**(1+2*s% number_of_backups_in_a_row))         
         if (s% doing_first_model_of_run) s% dt = 0.01d0*s% dt

         if (s% dt < s% min_timestep_limit) then
            write(*,1) 'dt', s% dt
            write(*,1) 'min_timestep_limit', s% min_timestep_limit
            call report_convergence_problems(s, 'dt < min_timestep_limit')
            do1_backup = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if
         
         if (s% generations > 2) then        
            !write(*,2) 'do1_backup omega_older(1)', s% model_number, s% omega_older(1)
            !write(*,2) 'do1_backup omega_old(1)', s% model_number, s% omega_old(1)
            call restore_older(s) ! set old = older
            s% generations = s% generations-1        
            call set_var_info(s, ierr)
            if (ierr /= 0) then
               write(*, *) 'do1_backup: set_var_info ierr', ierr
               s% result_reason = nonzero_ierr
               do1_backup = terminate
               return
            end if            
            call set_chem_names(s)            
            call set_current_to_old(s)
            call free_star_info_arrays(s)
            call allocate_star_info_arrays(s, ierr)
            if (ierr /= 0) then
               write(*, *) 'do1_backup: allocate_star_info_arrays ierr'
               s% result_reason = nonzero_ierr
               do1_backup = terminate
               return
            end if            
         else
            call set_current_to_old(s)         
         end if
         
         !write(*,2) 'do1_backup omega(1)', s% model_number, s% omega(1)
      
         s% dt_limit_ratio_old = 0 ! don't use predictive timestep control after backup
         s% dt_start = s% dt                   
         s% last_backup = s% model_number               
         s% num_backups = s% num_backups + 1
         if (s% num_backups > s% max_number_backups .and. s% max_number_backups >= 0) then
            write(*,2) 'num_backups', s% num_backups
            write(*,2) 'max_number_backups', s% max_number_backups
            call report_convergence_problems(s, 'num_backups > max_number_backups')
            s% termination_code = t_max_number_backups
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if

         s% number_of_backups_in_a_row = s% number_of_backups_in_a_row + 1
         if (s% number_of_backups_in_a_row > s% max_backups_in_a_row &
               .and. s% max_backups_in_a_row > 0) then
            write(*,2) 'number_of_backups_in_a_row', s% number_of_backups_in_a_row
            write(*,2) 'max_backups_in_a_row', s% max_backups_in_a_row
            call report_convergence_problems(s, 'too many backups in a row')
            s% termination_code = t_max_backups_in_a_row
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if
         
         s% timestep_hold = s% model_number + s% backup_hold         
         s% retry_cnt = 0
         s% model_number_for_last_retry = s% model_number
         s% have_start_values = .false.

         if (s% report_ierr) write(*, *) 'backup model_number', &
            s% model_number, s% num_backups
         
      end function do1_backup
      
      
      subroutine set_start_of_step_info(s, ierr)
         use report, only: do_report
         use hydro_vars, only: set_vars
         use mlt_info, only: set_gradT_excess_alpha
         use alloc, only: non_crit_get_work_array, non_crit_return_work_array

         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         logical :: trace
         integer :: nz
         
         include 'formats'
         
         ierr = 0
         trace = s% trace_evolve
         nz = s% nz
                  
         if (trace) write(*,*) 'call set_vars'
         call set_vars(s, s% dt, ierr) ! this does set_mixing_info too
         if (failed('set_vars')) return
         
         if (trace) write(*,*) 'call do_report'
         call do_report(s, ierr) ! set values in case used during step
         if (failed('do_report ierr')) return
         
         ! save a few things from start of step that will need later
         s% prev_Lmax = maxval(abs(s% L(1:nz)))
         if (s% rotation_flag) then
            s% surf_r_equatorial = s% r_equatorial(1)
         else
            s% surf_r_equatorial = s% r(1)
         end if
         s% starting_T_center = s% T(nz)
         s% surf_opacity = s% opacity(1)
         s% surf_csound = s% csound(1)
         s% surf_rho = s% rho(1)
         s% prev_Ledd = eval_Ledd(s)
         
         if (s% generations == 1) then
            s% surf_accel_grav_ratio = 0
         else
            s% surf_accel_grav_ratio = &
               (s% v_surf - s% v_surf_old)/(s% dt*s% grav(1))
         end if
         
         if (trace) write(*,*) 'call set_gradT_excess_alpha'
         call set_gradT_excess_alpha(s, ierr)
         if (failed('set_gradT_excess_alpha ierr')) return
         
         if (trace) write(*,*) 'call save_prev_mesh_info'
         call save_prev_mesh_info(ierr)
         if (failed('save_prev_mesh_info ierr')) return
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            if (s% report_ierr) write(*, *) 'set_start_of_step_info: ' // trim(str)
            s% result_reason = nonzero_ierr
         end function failed
         
         
         subroutine save_prev_mesh_info(ierr)
            integer, intent(out) :: ierr
            real(dp) :: xm
            integer :: k           
            ierr = 0         
            call do1_alloc(s% prev_mesh_xm, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            call do1_alloc(s% prev_mesh_lnS, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            call do1_alloc(s% prev_mesh_mu, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            xm = 0
            do k=1,nz
               s% prev_mesh_xm(k) = xm
               xm = xm + s% dm(k)
               s% prev_mesh_lnS(k) = s% lnS(k)
               s% prev_mesh_mu(k) = s% mu(k)
            end do
            s% prev_mesh_nz = nz
            s% have_prev_lnS = .false.
            s% have_prev_mu = .false.         
         end subroutine save_prev_mesh_info
         
         
         subroutine do1_alloc(p, sz, ierr)
            use alloc, only: non_crit_do1_alloc_if_necessary
            real(dp), pointer :: p(:)
            integer, intent(in) :: sz
            integer, intent(out) :: ierr
            call non_crit_do1_alloc_if_necessary( &
               s, p, sz, 'prepare_for_new_step', ierr)
         end subroutine do1_alloc

      
      end subroutine set_start_of_step_info


      integer function finish_step( &
            id, id_extra, do_photo, &
            how_many_extra_profile_columns, data_for_extra_profile_columns, &
            how_many_extra_history_columns, data_for_extra_history_columns, ierr)
         ! returns keep_going or terminate
         ! if don't return keep_going, then set result_reason to say why.
         use evolve_support, only: output
         use do_one_utils, only: do_save_profiles
         use history, only: write_history_info
         use utils_lib, only: free_iounit, number_iounits_allocated
         use alloc, only: size_work_arrays

         integer, intent(in) :: id, id_extra
         logical, intent(in) :: do_photo ! if true, then save "photo" for restart
         interface
            include 'extra_profile_cols.inc'
            include 'extra_history_cols.inc'
         end interface
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer, parameter :: nvals = 1, n_ivals = 0
         integer :: j, k, nz, &
            current_num_iounits_in_use, prev_num_iounits_in_use
         integer :: ivals(n_ivals)
         real(dp) :: vals(nvals)
         
         include 'formats'
         
         finish_step = terminate
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         nz = s% nz
         
         s% h1_czb_mass_prev = s% h1_czb_mass
         prev_num_iounits_in_use = number_iounits_allocated()
         
         finish_step = keep_going
         s% result_reason = result_reason_normal
                  
         if (s% need_to_save_profiles_now .and. s% write_profiles_flag) then
            call do_save_profiles( &
               s, id_extra, how_many_extra_profile_columns, data_for_extra_profile_columns, ierr)
            s% need_to_save_profiles_now = .false.
         end if
         
         call check(1)
         
         if (s% need_to_update_history_now .and. s% do_history_file) then
            call write_history_info( &
               s, id_extra, how_many_extra_history_columns, data_for_extra_history_columns, ierr)
            if (ierr /= 0) then
               finish_step = terminate
               if (s% report_ierr) write(*, *) 'finish_step: write_history_info ierr', ierr
               s% result_reason = nonzero_ierr
               return
            end if
            s% need_to_update_history_now = .false.
         end if
         
         call check(2)
         
         if (do_photo .or. &
               (s% photostep > 0 .and. mod(s% model_number, s% photostep) == 0)) then
               
            call output(id, ierr)

            if (ierr /= 0) then
               finish_step = terminate
               if (s% report_ierr) write(*, *) 'finish_step: output ierr', ierr
               s% result_reason = nonzero_ierr
               return
            end if
         
            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'finish_step after save photo xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do k=1,s% nz
                  write(*,3) 'lnr', s% model_number, k, s% xh(s% i_lnR, k)
                  write(*,3) 'i_rot', s% model_number, k, s% i_rot(k)
                  write(*,3) 'j_rot', s% model_number, k, s% j_rot(k)
                  write(*,3) 'omega', s% model_number, k, s% omega(k)
               end do
            end if
            
            
         end if
         
         call check(3)

         s% screening_mode_value = -1 ! force a new lookup for next step         
         s% doing_first_model_of_run = .false.
         s% number_of_backups_in_a_row = 0         
         
                  
         contains
         
         
         subroutine check(i)
            integer, intent(in) :: i
            include 'formats'
            !return
            
            current_num_iounits_in_use = number_iounits_allocated()
            if (current_num_iounits_in_use > 3 .and. &
                  current_num_iounits_in_use > prev_num_iounits_in_use) then
               write(*,2) 's% model_number', s% model_number
               write(*,2) 'prev_num_iounits_in_use', prev_num_iounits_in_use
               write(*,2) 'current_num_iounits_in_use', current_num_iounits_in_use
               write(*,2) 'i', i
               stop 'finish_step' 
            end if
            prev_num_iounits_in_use = current_num_iounits_in_use
         end subroutine check
         
         
      end function finish_step
      
      
      subroutine set_age(id, age, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: age
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% time = age*secyer
         s% star_age = age
         s% profile_age = age
         s% post_he_age = age
      end subroutine set_age



      end module evolve


