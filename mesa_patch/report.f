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

      module report

      use star_private_def
      use chem_def
      use utils_lib
      use star_utils
      use num_lib, only: find0
      
      use const_def, only: avo, kerg, pi, clight, crad, Rsun, Lsun, Msun, &
         secyer, ln10, mev_amu, ev2erg, two_thirds

      implicit none
      

      contains


      subroutine do_report(s, ierr)
         use rates_def, only: &
            i_rate, i_rate_dRho, i_rate_dT, std_reaction_Qs, std_reaction_neuQs
         use chem_def, only: ipp, icno, i3alf
         use star_utils, only: get_r_phot
         type (star_info), pointer :: s
         integer, intent(out) :: ierr 
      
         integer :: k, i, j, ic, nz, kcore
         real(dp) :: w1, radius, dr, dm, hpc, cur_m, cur_r, prev_r, &
            twoGmrc2, cur_h, prev_h, cur_he, non_fe_core_mass, &
            prev_he, cur_c, prev_c, luminosity, mstar, pdg, pdg_prev, &
            prev_m, cell_mass, wf, conv_time, mv, bminv, uminb, kh_time
         logical, parameter :: new_only = .false.
         integer, pointer :: net_iso(:)
         integer :: h1, h2, he3, he4, c12, n14, o16, ne20, si28
         
         include 'formats'
         
         ierr = 0

         nz = s% nz
         net_iso => s% net_iso

         h1 = net_iso(ih1)
         h2 = net_iso(ih2)
         he3 = net_iso(ihe3)
         he4 = net_iso(ihe4)
         c12 = net_iso(ic12)
         n14 = net_iso(in14)
         o16 = net_iso(io16)
         ne20 = net_iso(ine20)
         si28 = net_iso(isi28)

         luminosity = s% L(1) ! in ergs per second
         s% log_surface_luminosity = log10_cr(luminosity/Lsun) 
            ! log10(stellar luminosity in solar units)         
         s% L_phot = luminosity/Lsun
         s% photosphere_L = s% L_phot
         s% photosphere_r = get_r_phot(s)/Rsun
         radius = s% r(1)  !  radius in cm
         s% log_surface_radius = log10_cr(radius/Rsun) 
            ! log10(stellar radius in solar units)
         s% log_center_density = &
            log10_cr((s% xmstar*s% center_avg_value_dq)/ &
                     (volume_at_q(s% center_avg_value_dq) - volume_at_q(0d0)))
         s% log_max_temperature = maxval(s% lnT(1:nz))/ln10
         s% log_center_temperature = center_value(s, s% lnT)/ln10
         s% log_center_pressure = center_value(s, s% lnP)/ln10
         s% center_degeneracy = center_value(s, s% eta)    
         
         s% center_eps_grav = center_value(s, s% eps_grav)    
         s% center_eps_nuc = center_value(s, s% eps_nuc)
         
         s% d_center_eps_nuc_dlnT = center_value(s, s% d_epsnuc_dlnT)    
         s% d_center_eps_nuc_dlnd = center_value(s, s% d_epsnuc_dlnd)    
         s% center_non_nuc_neu = center_value(s, s% non_nuc_neu)    
         s% center_dL_dm = center_value(s, s% dL_dm)    
         
         s% center_gamma = center_value(s, s% gam)
         s% center_abar = center_value(s, s% abar)
         s% center_zbar = center_value(s, s% zbar)
         s% center_mu = center_value(s, s% mu)
         s% center_ye = center_value(s, s% ye)
         s% center_entropy = exp_cr(center_value(s, s% lnS))*amu/kerg
         s% max_entropy = exp_cr(maxval(s% lnS(1:nz)))*amu/kerg
         
         if (.not. s% rotation_flag) then
            s% omega(1:nz) = 0
            s% center_omega = 0
            s% center_omega_div_omega_crit = 0
         else
            s% center_omega = center_value(s, s% omega)
            s% center_omega_div_omega_crit = center_omega_div_omega_crit()
         end if
         
         s% log_surface_temperature = s% lnT(1)/ln10 ! log10(temperature at surface)
         s% log_surface_pressure = s% lnP(1)/ln10 ! log10(pressure at surface)
         s% log_surface_density = s% lnd(1)/ln10 ! log10(density at surface)
         s% log_surface_gravity = safe_log10_cr(s% grav(1)) ! log10(gravity at surface)
         
         do j=1,num_categories 
            s% L_by_category(j) = &
               dot_product(s% dm(1:nz), s% eps_nuc_categories(i_rate,j,1:nz))/Lsun
            s% center_eps_burn(j) = center_value_eps_burn(j,i_rate)
            s% center_eps_burn_dT(j) = center_value_eps_burn(j,i_rate_dT)
            s% center_eps_burn_dRho(j) = center_value_eps_burn(j,i_rate_dRho)
         end do

         s% power_nuc_burn = dot_product(s% dm(1:nz), s% eps_nuc(1:nz))/Lsun
         s% power_h_burn = s% L_by_category(ipp) + s% L_by_category(icno)
         
         s% power_he_burn = s% L_by_category(i3alf)
         s% power_nuc_neutrinos = &
            dot_product(s% dm(1:nz),s% eps_nuc_neu_total(1:nz))/Lsun
         s% power_nonnuc_neutrinos = &
            dot_product(s% dm(1:nz),s% non_nuc_neu(1:nz))/Lsun
         s% power_neutrinos = s% power_nuc_neutrinos + s% power_nonnuc_neutrinos

         s% L_nuc_burn_total = s% power_nuc_burn
         
         if (s% v_flag) then
            s% v_surf = s% v(1)
         else
            s% v_surf = s% r(1)*s% dlnR_dt(1)
         end if
         
         call set_surf_avg_rotation_info(s)
         
         s% time_step = s% dt/secyer         ! timestep in years
         s% star_age = s% time/secyer
         if ( s% model_number <= 0 ) s% star_age = 0d0
         mstar = s% mstar
         s% star_mass = mstar/Msun             ! stellar mass in solar units
         s% star_mdot = s% mstar_dot/(Msun/secyer)      ! dm/dt in msolar per year
         
         if (s% M_center >= 0.0) then
            call new_eval_kh_timescale(s,kh_time)
            s% kh_timescale = kh_time / secyer
         else
            s% kh_timescale = eval_kh_timescale(s% cgrav(1), mstar, radius, luminosity)/secyer
         endif 
         ! kelvin-helmholtz timescale in years (about 1.6x10^7 for the sun)         
         s% nuc_timescale = 1d10*s% star_mass/(luminosity/Lsun) 
         ! nuclear timescale in years (e.g., about 10^10 years for sun)         
         s% dynamic_timescale = 2*pi*sqrt(radius*radius*radius/(s% cgrav(1)*mstar))

         if (h1 /= 0) then
            s% center_h1 = center_avg_x(s,h1)
            s% surface_h1 = surface_avg_x(s,h1)
         end if
         if (he3 /= 0) then
            s% center_he3 = center_avg_x(s,he3)
            s% surface_he3 = surface_avg_x(s,he3)
         end if
         if (he4 /= 0) then
            s% center_he4 = center_avg_x(s,he4)
            s% surface_he4 = surface_avg_x(s,he4)
         end if
         if (c12 /= 0) then
            s% center_c12 = center_avg_x(s,c12)
            s% surface_c12 = surface_avg_x(s,c12)
         end if
         if (n14 /= 0) then
            s% center_n14 = center_avg_x(s,n14)
            s% surface_n14 = surface_avg_x(s,n14)
         end if
         if (o16 /= 0) then
            s% center_o16 = center_avg_x(s,o16)
            s% surface_o16 = surface_avg_x(s,o16)
         end if
         if (ne20 /= 0) then
            s% center_ne20 = center_avg_x(s,ne20)
            s% surface_ne20 = surface_avg_x(s,ne20)
         end if
         if (si28 /= 0) then
            s% center_si28 = center_avg_x(s,si28)
         end if
         
         ! FYI profile stuff
         do k=1,nz
         
            s% entropy(k) = exp_cr(s% lnS(k))/(avo*kerg)
            if (is_bad_num(s% entropy(k))) then
               ierr = -1
               write(*,2) 'report: s% entropy(k)', k, s% entropy(k)
               return
            end if
            
            if (k == nz) then
               dr = s% r(k) - s% R_center
            else
               dr = s% r(k) - s% r(k+1)
            end if
            s% dr_div_csound(k) = dr/s% csound(k)
            if (is_bad_num(s% dr_div_csound(k))) then
               ierr = -1
               write(*,2) 'report: s% dr_div_csound(k)', &
                  k, s% dr_div_csound(k), dr, s% csound(k)
               return
            end if

            if (.not. s% v_flag) then
               s% velocity(k) = s% r(k) * s% dlnR_dt(k)
               if (is_bad_num(s% velocity(k))) then
                  ierr = -1
                  write(*,2) 'report: s% velocity(k)', &
                     k, s% velocity(k), s% r(k), s% dlnR_dt(k), s% dt
                  return
               end if
            else
               s% velocity(k) = s% v(k)
               if (is_bad_num(s% velocity(k))) then
                  ierr = -1
                  write(*,2) 'report: s% velocity(k)', k, s% velocity(k)
                  return
               end if
            end if
         
            s% v_div_csound(k) = s% velocity(k)/s% csound_at_face(k)
            if (is_bad_num(s% v_div_csound(k))) then
               ierr = -1
               write(*,2) 'report: s% v_div_csound(k)', k, s% v_div_csound(k), &
                  s% velocity(k), s% csound_at_face(k)
               return
            end if

         end do

         if (s% photosphere_r*Rsun >= s% r(1)) then
            s% photosphere_acoustic_r = sum(s% dr_div_csound(1:nz)) + &
               (s% photosphere_r*Rsun - s% r(1))/s% csound(1)
         else
            do k=2,nz
               if (s% photosphere_r*Rsun > s% r(k)) then
                  s% photosphere_acoustic_r = sum(s% dr_div_csound(k:nz)) + &
                     (s% photosphere_r*Rsun - s% r(k))/s% csound(k-1)
                  exit
               end if
            end do
         end if
         
         if (.not. s% get_delta_nu_from_scaled_solar) then
            s% delta_nu = 1d6/(2*s% photosphere_acoustic_r) ! microHz
         else
            s% delta_nu = &
         		s% delta_nu_sun*sqrt(s% star_mass)*pow3(s% Teff/s% Teff_sun)/s% L_phot
         end if
         
         s% nu_max = &
            s% nu_max_sun*s% star_mass/(pow2(s% photosphere_r)*sqrt(s% Teff/s% Teff_sun))
         s% acoustic_cutoff = &
            0.5d6*s% grav(1)*sqrt(s% gamma1(1)*s% rho(1)/s% P(1))
         call get_delta_Pg(s, s% nu_max, s% delta_Pg)

         call get_tau(s, s% tau)
         
         call get_mass_info(s, s% dm, ierr)
         if (failed('get_mass_info')) return
         
         call get_power_info(s, s% dm, ierr)
         if (failed('get_power_info')) return
         
         call get_mixing_regions(s, ierr)
         if (failed('get_mixing_regions')) return
         
         call get_neutrino_fluxes(s, ierr)
         if (failed('get_neutrino_fluxes')) return

         call check_tau(10d0, s% tau10_radius, s% tau10_mass, &
               s% tau10_lgP, s% tau10_lgT, s% tau10_lgRho, s% tau10_L)
         call check_tau(100d0, s% tau100_radius, s% tau100_mass, &
               s% tau100_lgP, s% tau100_lgT, s% tau100_lgRho, s% tau100_L)
         
         call set_mass_conv_core
         call find_conv_mx_regions
         call find_mx_regions
         call get_burn_zone_info(s, ierr)
         if (failed('get_burn_zone_info')) return

         s% fe_core_infall = 0
         s% non_fe_core_infall = 0
         if (s% v_flag) then
            if (s% fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*s% fe_core_mass) cycle
                  if (-s% v(k) > s% fe_core_infall) &
                     s% fe_core_infall = -s% v(k)
               end do
            end if
            non_fe_core_mass = max( &
               s% he_core_mass, s% c_core_mass, s% o_core_mass, s% si_core_mass)
            if (non_fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*non_fe_core_mass) cycle
                  if (s% m(k) < Msun*s% fe_core_mass) exit
                  if (-s% v(k) > s% non_fe_core_infall) &
                     s% non_fe_core_infall = -s% v(k)
               end do
            end if
         end if

         
         contains

         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed
         
         
         subroutine set_mass_conv_core
            integer :: j, nz
            s% mass_conv_core = 0
            nz = s% nz
            do j = 1, s% n_conv_regions
               if (s% cz_bot_mass(j) <= s% m(nz)) then
                  s% mass_conv_core = s% cz_top_mass(j)/Msun
                  exit
               end if
            end do
         end subroutine set_mass_conv_core

         
         subroutine check_tau(check, tau_radius, tau_mass, tau_lgP, tau_lgT, tau_lgd, tau_L)
            use num_lib, only: binary_search
            real(dp), intent(in) :: check
            real(dp), intent(out) :: tau_radius, tau_mass, tau_lgP, tau_lgT, tau_lgd, tau_L
            integer :: k, hint, nz
            real(dp) :: frac, v00, vp1, taup1, tau00
            
            include 'formats'
            
            hint = 0
            nz = s% nz
            k = binary_search(nz, s% tau, hint, check)
            
            if (k < 1 .or. k > nz-2) then
               tau_radius=0; tau_mass=0; tau_lgP=0; tau_lgT=0; tau_lgd=0; tau_L=0
               return
            end if
            
            tau00 = s% tau(k)
            taup1 = s% tau(k+1)

            ! tau00 <= check < taup1
            frac = (check - tau00)/(taup1 - tau00)
            
            if (frac > 1 .or. frac < 0) then
               tau_radius=0; tau_mass=0; tau_lgP=0; tau_lgT=0; tau_lgd=0; tau_L=0
               return
            end if

            tau_lgP = (s% lnP(k) + (s% lnP(k+1)-s% lnP(k))*frac)/ln10
            tau_lgT = (s% lnT(k) + (s% lnT(k+1)-s% lnT(k))*frac)/ln10
            tau_lgd = (s% lnd(k) + (s% lnd(k+1)-s% lnd(k))*frac)/ln10
            
            v00 = (s% r(k+1) + s% r(k))/2
            vp1 = (s% r(k+2) + s% r(k+1))/2
            tau_radius = (v00 + (vp1 - v00)*frac)/Rsun
            
            v00 = (s% q(k+1) + s% q(k))/2
            vp1 = (s% q(k+2) + s% q(k+1))/2
            tau_mass = (s% M_center + s% xmstar*(v00 + (vp1 - v00)*frac))/Msun
            
            v00 = (s% L(k+1) + s% L(k))/2
            vp1 = (s% L(k+2) + s% L(k+1))/2
            tau_L = (v00 + (vp1 - v00)*frac)/Lsun

         end subroutine check_tau


         real(dp) function volume_at_q(q)
            use interp_1d_def
            use interp_1d_lib
            real(dp), intent(in) :: q
            real(dp) :: vp2, vp1, v00, vm1
            integer, parameter :: n_old = 4, n_new = 1, nwork = pm_work_size
            real(dp) :: qlo, x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
            integer :: k, nz, k00, ierr
            real(dp), target :: work_ary(n_old*nwork)
            real(dp), pointer :: work(:)
            work => work_ary
            
            include 'formats'
            nz = s% nz

            !write(*,1) 'q', q
            !write(*,1) 's% q(nz)', s% q(nz)
            
            if (q == 0d0) then
               volume_at_q = 0
               return
            end if
            
            if (q <= s% q(nz)) then
               volume_at_q = (q/s% q(nz))*(4*pi/3)*s% r(nz)*s% r(nz)*s% r(nz)
               return
            end if
            k00 = 1
            do k=nz-1,2,-1
               if (s% q(k) >= q) then
                  k00 = k; exit
               end if
            end do
            if (k00 == 1) then
               volume_at_q = (4*pi/3)*q*s% r(1)*s% r(1)*s% r(1)
               write(*,1) 'volume_at_q', volume_at_q
               return
            end if
            
            x_old(1) = 0
            if (k00+1 == nz) then
               v_old(1) = s% R_center*s% R_center*s% R_center
               qlo = 0
            else
               v_old(1) = s% r(k00+2)*s% r(k00+2)*s% r(k00+2)
               qlo = s% q(k00+2)
            end if
            
            x_old(2) = s% dq(k00+1)
            v_old(2) = s% r(k00+1)*s% r(k00+1)*s% r(k00+1)

            x_old(3) = x_old(2) + s% dq(k00)
            v_old(3) = s% r(k00)*s% r(k00)*s% r(k00)

            x_old(4) = x_old(3) + s% dq(k00-1)
            v_old(4) = s% r(k00-1)*s% r(k00-1)*s% r(k00-1)
            
            ierr = 0
            x_new(1) = q-qlo
            call interpolate_vector( &
               n_old, x_old, n_new, x_new, v_old, v_new, interp_pm, nwork, work, &
               'report volume_at_q', ierr)
            if (ierr /= 0) then
               write(*,*) 'volume_at_q: failed in interpolate_vector'
               volume_at_q = (v_old(2) + v_old(3))/2
               return
            end if
            
            volume_at_q = (4*pi/3)*v_new(1)
            
            return

            write(*,1) 'x_old', x_old
            write(*,1) 'v_old', v_old
            write(*,1) 'qlo', qlo
            write(*,1) 'x_new', x_new
            write(*,1) 'v_new', v_new
            write(*,1) 'volume_at_q', volume_at_q
            write(*,1) 'mass_at_q', q*s% mstar
            write(*,1) 'mass_at_q/volume_at_q', q*s% mstar/volume_at_q
            write(*,1) 'rho(k00)', s% rho(k00)
            write(*,1) 'rho(k00+1)', s% rho(k00+1)

         end function volume_at_q


         real(dp) function center_value_eps_burn(j,ir)
            integer, intent(in) :: j, ir
            real(dp) :: sum_x, sum_dq, dx, dq
            integer :: k
            sum_x = 0
            sum_dq = 0
            do k = s% nz, 1, -1
               dq = s% dq(k)
               dx = s% eps_nuc_categories(ir,j,k)*dq
               if (sum_dq+dq >= s% center_avg_value_dq) then
                  sum_x = sum_x + dx*(s% center_avg_value_dq - sum_dq)/dq
                  sum_dq = s% center_avg_value_dq
                  exit
               end if
               sum_x = sum_x + dx
               sum_dq = sum_dq + dq
            end do
            center_value_eps_burn = sum_x/sum_dq
         end function center_value_eps_burn         


         real(dp) function center_omega_div_omega_crit()
            real(dp) :: sum_x, sum_dq, dx, dq
            integer :: k
            center_omega_div_omega_crit = 0
            if (.not. s% rotation_flag) return
            sum_x = 0
            sum_dq = 0
            do k = s% nz, 1, -1
               dq = s% dq(k)
               dx = dq*s% omega(k)/omega_crit(s,k)
               if (sum_dq+dq >= s% center_avg_value_dq) then
                  sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
                  sum_dq = s% center_avg_value_dq
                  exit
               end if
               sum_x = sum_x + dx
               sum_dq = sum_dq + dq
            end do
            center_omega_div_omega_crit = min(1d0,sum_x/sum_dq)
         end function center_omega_div_omega_crit         
         
         
         subroutine find_conv_mx_regions
            use mlt_def, only: convective_mixing
            real(dp) :: conv_mx1_dq, conv_mx2_dq, mx_dq
            integer :: i, ktop, kbot, conv_mx1, conv_mx2
            
            include 'formats'
            
            s% largest_conv_mixing_region = 0
            
            s% conv_mx1_top = 0
            s% conv_mx1_bot = 0
            s% conv_mx2_top = 0
            s% conv_mx2_bot = 0
            
            s% conv_mx1_top_r = 0
            s% conv_mx1_bot_r = 0
            s% conv_mx2_top_r = 0
            s% conv_mx2_bot_r = 0
            
            if (s% num_mixing_regions == 0) return
            
            conv_mx1 = 0
            conv_mx2 = 0
            conv_mx1_dq = 0
            conv_mx2_dq = 0
            
            do i = 1, s% num_mixing_regions
               if (s% mixing_region_type(i) /= convective_mixing) cycle
               mx_dq = s% q(s% mixing_region_top(i)) - s% q(s% mixing_region_bottom(i))
               if (mx_dq > conv_mx1_dq) then
                  conv_mx2_dq = conv_mx1_dq; conv_mx1_dq = mx_dq
                  conv_mx2 = conv_mx1; conv_mx1 = i
               else if (mx_dq > conv_mx2_dq) then
                  conv_mx2_dq = mx_dq
                  conv_mx2 = i
               end if
            end do

            if (conv_mx1 > 0) then
               s% largest_conv_mixing_region = conv_mx1
               ktop = s% mixing_region_top(conv_mx1)
               kbot = s% mixing_region_bottom(conv_mx1)
               s% conv_mx1_top = s% q(ktop) - s% mixing_type_change_dq(ktop)
               s% conv_mx1_bot = s% q(kbot) - s% mixing_type_change_dq(kbot)
               s% conv_mx1_top_r = s% r(ktop)/Rsun
               s% conv_mx1_bot_r = s% r(kbot)/Rsun
            end if
            
            if (conv_mx2 > 0) then
               ktop = s% mixing_region_top(conv_mx2)
               kbot = s% mixing_region_bottom(conv_mx2)
               s% conv_mx2_top = s% q(ktop) - s% mixing_type_change_dq(ktop)
               s% conv_mx2_bot = s% q(kbot) - s% mixing_type_change_dq(kbot)
               s% conv_mx2_top_r = s% r(ktop)/Rsun
               s% conv_mx2_bot_r = s% r(kbot)/Rsun
            end if

         end subroutine find_conv_mx_regions
         
         
         subroutine find_mx_regions
            real(dp) :: mx1_dq, mx2_dq, mx_dq
            integer :: i, ktop, kbot, mx1_top_region, mx1_bottom_region, &
               mx2_top_region, mx2_bottom_region, &
               current_top_region, current_bottom_region, &
               current_top_point, current_bottom_point
            
            include 'formats'
            
            s% mx1_top = 0
            s% mx1_bot = 0
            s% mx2_top = 0
            s% mx2_bot = 0
            
            s% mx1_top_r = 0
            s% mx1_bot_r = 0
            s% mx2_top_r = 0
            s% mx2_bot_r = 0
            
            if (s% num_mixing_regions == 0) return
            
            mx1_top_region = 0
            mx1_bottom_region = 0
            mx1_dq = 0
            
            mx2_top_region = 0
            mx2_bottom_region = 0
            mx2_dq = 0
            
            i = 1
            do
               if (i > s% num_mixing_regions) exit
               current_top_region = i
               current_top_point = s% mixing_region_top(current_top_region)
               do
                  current_bottom_region = i
                  current_bottom_point = s% mixing_region_bottom(current_bottom_region)
                  i = i+1
                  if (i > s% num_mixing_regions) exit
                  if (s% mixing_region_top(i) /= current_bottom_point+1) exit
               end do
               mx_dq = s% q(current_top_point) - s% q(current_bottom_point)
               if (mx_dq > mx1_dq) then
                  mx2_dq = mx1_dq; mx1_dq = mx_dq
                  mx2_top_region = mx1_top_region
                  mx1_top_region = current_top_region
                  mx2_bottom_region = mx1_bottom_region
                  mx1_bottom_region = current_bottom_region
               else if (mx_dq > mx2_dq) then
                  mx2_dq = mx_dq
                  mx2_top_region = current_top_region
                  mx2_bottom_region = current_bottom_region
               end if
            end do
            
            if (mx1_top_region > 0) then
               ktop = s% mixing_region_top(mx1_top_region)
               kbot = s% mixing_region_bottom(mx1_bottom_region)
               s% mx1_top = s% q(ktop) - s% mixing_type_change_dq(ktop)
               s% mx1_bot = s% q(kbot) - s% mixing_type_change_dq(kbot)
               s% mx1_top_r = s% r(ktop)/Rsun
               s% mx1_bot_r = s% r(kbot)/Rsun
            end if
            
            if (mx2_top_region > 0) then
               ktop = s% mixing_region_top(mx2_top_region)
               kbot = s% mixing_region_bottom(mx2_bottom_region)
               s% mx2_top = s% q(ktop) - s% mixing_type_change_dq(ktop)
               s% mx2_bot = s% q(kbot) - s% mixing_type_change_dq(kbot)
               s% mx2_top_r = s% r(ktop)/Rsun
               s% mx2_bot_r = s% r(kbot)/Rsun
            end if

         end subroutine find_mx_regions


      end subroutine do_report


      real(dp) function surface_avg_x(s,j)
         use utils_lib, only: is_bad_num
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq
         integer :: k
         include 'formats'
         sum_x = 0
         sum_dq = 0
         do k = 1, s% nz
            sum_x = sum_x + s% xa(j,k)*s% dq(k)
            sum_dq = sum_dq + s% dq(k)
            if (sum_dq >= s% surface_avg_abundance_dq) exit
         end do
         surface_avg_x = sum_x/sum_dq
         
         return
         
         if (is_bad_num(surface_avg_x)) then
            write(*,1) 'surface_avg_x ' // trim(chem_isos% name(s% chem_id(j))), surface_avg_x
            write(*,1) 'sum_x ' // trim(chem_isos% name(s% chem_id(j))), sum_x
            write(*,1) 'sum_dq', sum_dq
            stop 'debug: surface_avg_x'
         end if
      end function surface_avg_x         


      real(dp) function center_avg_x(s,j)
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = s% xa(j,k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_avg_x = sum_x/sum_dq
      end function center_avg_x


      subroutine get_burn_zone_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: burn_min1, burn_min2
         integer :: i, i_start
         include 'formats'
         burn_min1 = s% burn_min1; burn_min2 = s% burn_min2
         ! up to 3 zones where eps_nuc > burn_min1 erg/g/s
         ! for each zone have 4 numbers: start1, start2, end2, end1
         ! start1 is mass of inner edge where first goes > burn_min1 (or -20 if none such)
         ! start2 is mass of inner edge where first zone reaches burn_min2 erg/g/sec (or -20 if none such)
         ! end2 is mass of outer edge where first zone drops back below burn_min2 erg/g/s
         ! end1 is mass of outer edge where first zone ends (i.e. eps_nuc < burn_min1)
         ! similar for second and third zones
         i_start = s% nz
         do i=1,3
            call find_epsnuc_zone(s, i_start, &
               s% burn_zone_mass(1,i), s% burn_zone_mass(2,i), &
               s% burn_zone_mass(3,i), s% burn_zone_mass(4,i), &
               s% burn_min1, s% burn_min2, ierr)
         end do
      end subroutine get_burn_zone_info


      subroutine find_epsnuc_zone( &
            s, i_start, bzm_1, bzm_2, bzm_3, bzm_4, burn_min1, burn_min2, ierr)
         use const_def, only:Msun
         type (star_info), pointer :: s
         integer, intent(inout) :: i_start
         real(dp), intent(out) :: bzm_1, bzm_2, bzm_3, bzm_4
         real(dp), intent(in) :: burn_min1, burn_min2
         integer, intent(out) :: ierr

         real(dp), parameter :: null_zone = -20
         integer :: i, burn_zone
         real(dp) :: prev_m, prev_x, cur_m, cur_x
         ierr = 0
         bzm_1 = null_zone; bzm_2 = null_zone; bzm_3 = null_zone; bzm_4 = null_zone
         burn_zone = 0 ! haven't entered the zone yet
         if (i_start .ne. s% nz) then
            i = i_start+1
            prev_m = s% m(i)
            prev_x = s% eps_nuc(i)
         else ! keep the compiler happy
            prev_m = 0
            prev_x = 0
         end if
         do i = i_start, 1, -1
            cur_m = s% m(i)
            cur_x = s% eps_nuc(i)
            select case (burn_zone)
               case (0)
                  if ( cur_x .gt. burn_min2 ) then
                     if ( i .eq. s% nz ) then ! use star center as start of zone
                        bzm_2 = 0d0
                     else ! interpolate to estimate where rate reached burn_min1 
                        bzm_2 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     end if
                     bzm_1 = bzm_2
                     burn_zone = 2
                  elseif ( cur_x .gt. burn_min1 ) then
                     if ( i .eq. s% nz ) then ! use star center as start of zone
                        bzm_1 = 0d0
                     else ! interpolate to estimate where rate reached burn_min1 
                        bzm_1 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     end if
                     burn_zone = 1
                  end if
               case (1) ! in the initial eps > burn_min1 region
                  if ( cur_x .gt. burn_min2 ) then
                     bzm_2 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     burn_zone = 2
                  else if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     i_start = i
                     return
                  end if
               case (2) ! in the initial eps > burn_min2 region
                  if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     bzm_3 = bzm_4
                     i_start = i
                     return
                  end if
                  if ( cur_x .lt. burn_min2 ) then
                     bzm_3 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     burn_zone = 3
                  end if
               case (3) ! in the final eps > burn_min1 region
                  if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     i_start = i
                     return
                  end if
               case default
                  ierr = -1
                  write(*,*) 'error in find_eps_nuc_zone'
                  return
            end select
            prev_m = cur_m; prev_x = cur_x
         end do
         i_start = 0
         select case (burn_zone)
            case (0)
               return
            case (1)
               bzm_4 = cur_m
            case (2)
               bzm_3 = cur_m
               bzm_4 = cur_m
            case (3)
               bzm_4 = cur_m
            case default
               ierr = -1
               write(*,*) 'error in find_eps_nuc_zone'
               return
         end select
      end subroutine find_epsnuc_zone

      
      
      subroutine get_mass_info(s, cell_masses, ierr)
         use utils_lib, only: is_bad_num
         type (star_info), pointer :: s
         real(dp), pointer, intent(in) :: cell_masses(:) 
         integer, intent(out) :: ierr 
         
         integer :: k, nz, j, nzlo, nzhi, kbdy, nzlo_prev
         real(dp) :: cell_mass
         integer, pointer :: net_iso(:)    
         
         include 'formats'
         
         ierr = 0
         nz = s% nz
         net_iso => s% net_iso
         
         s% star_mass_h1=0d0; s% star_mass_he3=0d0; s% star_mass_he4=0d0; s% star_mass_c12 = 0d0
         s% star_mass_n14 = 0d0; s% star_mass_o16 = 0d0; s% star_mass_ne20 = 0d0
         
         do k = 1, nz
            cell_mass = cell_masses(k)
            do j=1, s% species
               if (s% chem_id(j) == ih1) then
                  s% star_mass_h1 = s% star_mass_h1 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ihe3) then
                  s% star_mass_he3 = s% star_mass_he3 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ihe4) then
                  s% star_mass_he4 = s% star_mass_he4 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ic12) then
                  s% star_mass_c12 = s% star_mass_c12 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == in14) then
                  s% star_mass_n14 = s% star_mass_n14 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == io16) then
                  s% star_mass_o16 = s% star_mass_o16 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ine20) then
                  s% star_mass_ne20 = s% star_mass_ne20 + cell_mass*s% xa(j, k)
               end if
            end do
         end do

         s% star_mass_h1 = s% star_mass_h1 / Msun
         s% star_mass_he3 = s% star_mass_he3 / Msun
         s% star_mass_he4 = s% star_mass_he4 / Msun
         s% star_mass_c12 = s% star_mass_c12 / Msun
         s% star_mass_n14 = s% star_mass_n14 / Msun
         s% star_mass_o16 = s% star_mass_o16 / Msun
         s% star_mass_ne20 = s% star_mass_ne20 / Msun
         
         call get_core_info(s)
         call get_trace_mass_location_info(s)
         call get_max_T_location_info(s)
      
      end subroutine get_mass_info
      
      
      subroutine get_max_T_location_info(s)
         type (star_info), pointer :: s
         integer :: k, kk
         k = maxloc(s% lnT(1:s% nz),dim=1)
         s% max_T_lgT = s% lnT(k)/ln10
         s% max_T_lgRho = s% lnd(k)/ln10
         s% max_T_lgP = s% lnP(k)/ln10
         if (k == s% nz) then
            s% max_T_mass = s% M_center/Msun ! (Msun)
            s% max_T_radius = s% R_center/Rsun ! (Rsun)
            s% max_T_L = s% L_center/Lsun! (Lsun)
         else
            s% max_T_mass = s% m(k)/Msun ! (Msun)
            s% max_T_radius = s% r(k)/Rsun ! (Rsun)
            s% max_T_L = s% L(k)/Lsun! (Lsun)
         end if
         s% max_T_entropy = s% entropy(k)
         s% max_T_eps_nuc = s% eps_nuc(k) ! (erg/g/s)
         s% max_T_lgP_thin_shell = &
            log10_cr(s% cgrav(k)*s% m(k)*(s% m(1)-s% m(k))/(4*pi*pow4(s% r(k))))
         s% max_T_shell_binding_energy = 0d0
         do kk = 1, k
            s% max_T_shell_binding_energy = s% max_T_shell_binding_energy + &
               s% dm(kk)*(exp_cr(s% lnE(kk)) + &
                     s% P(kk)/s% rho(kk) - s% cgrav(k)*s% m_grav(kk)/s% r(kk) + &
                     0.5d0*s% velocity(kk)*s% velocity(kk))
         end do
      end subroutine get_max_T_location_info
      
      
      subroutine get_trace_mass_location_info(s)
         type (star_info), pointer :: s
         real(dp) :: q, m
         integer :: kbdy
         include 'formats'
         q = max(0d0, min(1d0, s% trace_mass_location / s% star_mass))
         call get_info_at_q(s, q, &
            kbdy, m, s% trace_mass_radius, s% trace_mass_lgT, &
            s% trace_mass_lgRho, s% trace_mass_L, s% trace_mass_v, &
            s% trace_mass_lgP, s% trace_mass_g, s% trace_mass_X, s% trace_mass_Y, &
            s% trace_mass_edv_H, s% trace_mass_edv_He, s% trace_mass_scale_height, &
            s% trace_mass_dlnX_dr, s% trace_mass_dlnY_dr, s% trace_mass_dlnRho_dr, &
            s% trace_mass_omega, s% trace_mass_omega_div_omega_crit)    
      end subroutine get_trace_mass_location_info


      subroutine get_info_at_surface( &
            s, bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, bdy_time, &
            bdy_omega, bdy_omega_div_omega_crit)
         type (star_info), pointer :: s
         real(dp), intent(out) :: &
            bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_omega, bdy_omega_div_omega_crit, bdy_time
         
         real(dp) :: bdy_omega_crit
         
         bdy_time = 0
         bdy_m = s% star_mass
         bdy_r = s% r(1)/Rsun
         bdy_lgT = s% lnT(1)/ln10
         bdy_lgRho = s% lnd(1)/ln10
         bdy_L = s% L(1)/Lsun
         if (s% v_flag) then
            bdy_v = s% v(1)
         else
            bdy_v = s% r(1)*s% dlnR_dt(1)
         end if
         if (s% rotation_flag) then
            bdy_omega = s% omega_avg_surf
            bdy_omega_crit = s% omega_crit_avg_surf
            bdy_omega_div_omega_crit = s% w_div_w_crit_avg_surf
         else
            bdy_omega = 0
            bdy_omega_div_omega_crit = 0
         end if

      end subroutine get_info_at_surface
      

      subroutine get_core_info(s)
         use num_lib, only: find0
         type (star_info), pointer :: s
         
         integer :: k, j, A_max, h1, he4, c12, o16, si28, species, nz
         integer, pointer :: net_iso(:)
         logical :: have_he, have_c, have_o, have_si, have_fe
         real(dp) :: sumx, min_x
         integer, parameter :: &
            A_max_fe = 47, &
            A_max_si = 28, &
            A_max_o = 16, &
            A_max_c = 12, &
            A_max_he = 4
         
         include 'formats'

         net_iso => s% net_iso
         species = s% species
         nz = s% nz
         h1 = net_iso(ih1)
         he4 = net_iso(ihe4)
         c12 = net_iso(ic12)
         o16 = net_iso(io16)
         si28 = net_iso(isi28)
         min_x = s% min_boundary_fraction
         
         call clear_core_info( &
            s% neutron_rich_core_mass, s% neutron_rich_core_radius, s% neutron_rich_core_lgT, &
            s% neutron_rich_core_lgRho, s% neutron_rich_core_L, s% neutron_rich_core_v, &
            s% neutron_rich_core_omega, s% neutron_rich_core_omega_div_omega_crit)
         
         do k=1,nz
            if (s% Ye(k) <= s% neutron_rich_core_boundary_Ye_max) then
               call set_core_info(s, k, &
                  s% neutron_rich_core_mass, s% neutron_rich_core_radius, s% neutron_rich_core_lgT, &
                  s% neutron_rich_core_lgRho, s% neutron_rich_core_L, s% neutron_rich_core_v, &
                  s% neutron_rich_core_omega, s% neutron_rich_core_omega_div_omega_crit)
               exit
            end if
         end do

         call clear_core_info( &
            s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
            s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
            s% fe_core_omega, s% fe_core_omega_div_omega_crit)
         call clear_core_info( &
            s% si_core_mass, s% si_core_radius, s% si_core_lgT, &
            s% si_core_lgRho, s% si_core_L, s% si_core_v, &
            s% si_core_omega, s% si_core_omega_div_omega_crit)
         call clear_core_info( &
            s% o_core_mass, s% o_core_radius, s% o_core_lgT, &
            s% o_core_lgRho, s% o_core_L, s% o_core_v, &
            s% o_core_omega, s% o_core_omega_div_omega_crit)
         call clear_core_info( &
            s% c_core_mass, s% c_core_radius, s% c_core_lgT, &
            s% c_core_lgRho, s% c_core_L, s% c_core_v, &
            s% c_core_omega, s% c_core_omega_div_omega_crit)
         call clear_core_info( &
            s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
            s% he_core_lgRho, s% he_core_L, s% he_core_v, &
            s% he_core_omega, s% he_core_omega_div_omega_crit)

         have_he = .false.
         have_c = .false.
         have_o = .false.
         have_si = .false.
         have_fe = .false.
         
         do k=1,nz
         
            j = maxloc(s% xa(1:species,k),dim=1)
            A_max = chem_isos% Z_plus_N(s% chem_id(j))
            
            if (.not. have_fe) then
               if (s% fe_core_boundary_si28_fraction < 0) then
                  if (A_max >= A_max_fe) then
                     call set_core_info(s, k, &
                        s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
                        s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
                        s% fe_core_omega, s% fe_core_omega_div_omega_crit)
                     exit
                  end if
               else if (si28 /= 0) then
                  if (s% xa(si28,k) <= s% fe_core_boundary_si28_fraction) then
                     sumx = 0
                     do j=1,species
                        if (chem_isos% Z_plus_N(s% chem_id(j)) >= A_max_fe) &
                           sumx = sumx + s% xa(j,k)
                     end do
                     if (sumx >= min_x) then
                        call set_core_info(s, k, &
                           s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
                           s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
                           s% fe_core_omega, s% fe_core_omega_div_omega_crit)
                        exit
                     end if
                  end if
               end if
            end if
            
            if (.not. have_si) then
               if (s% si_core_boundary_o16_fraction < 0) then
                  if (A_max >= A_max_si) then
                     call set_core_info(s, k, &
                        s% si_core_mass, s% si_core_radius, s% si_core_lgT, &
                        s% si_core_lgRho, s% si_core_L, s% si_core_v, &
                        s% si_core_omega, s% si_core_omega_div_omega_crit)
                     have_si = .true.
                  end if
               else if (o16 /= 0 .and.si28 /= 0) then
                  if (s% xa(o16,k) <= s% si_core_boundary_o16_fraction .and. &
                      s% xa(si28,k) >= min_x) then
                     call set_core_info(s, k, &
                        s% si_core_mass, s% si_core_radius, s% si_core_lgT, &
                        s% si_core_lgRho, s% si_core_L, s% si_core_v, &
                        s% si_core_omega, s% si_core_omega_div_omega_crit)
                     have_si = .true.
                  end if
               end if
            end if
            
            if (.not. have_o) then
               if (s% o_core_boundary_c12_fraction < 0) then
                  if (A_max >= A_max_o) then
                     call set_core_info(s, k, &
                        s% o_core_mass, s% o_core_radius, s% o_core_lgT, &
                        s% o_core_lgRho, s% o_core_L, s% o_core_v, &
                        s% o_core_omega, s% o_core_omega_div_omega_crit)
                     have_o = .true.
                  end if
               else if (c12 /= 0 .and. o16 /= 0) then
                  if (s% xa(c12,k) <= s% o_core_boundary_c12_fraction .and. &
                      s% xa(o16,k) >= min_x) then
                     call set_core_info(s, k, &
                        s% o_core_mass, s% o_core_radius, s% o_core_lgT, &
                        s% o_core_lgRho, s% o_core_L, s% o_core_v, &
                        s% o_core_omega, s% o_core_omega_div_omega_crit)
                     have_o = .true.
                  end if
               end if
            end if            
            
            if (.not. have_c) then
               if (s% c_core_boundary_he4_fraction < 0) then
                  if (A_max >= A_max_c) then
                     call set_core_info(s, k, &
                        s% c_core_mass, s% c_core_radius, s% c_core_lgT, &
                        s% c_core_lgRho, s% c_core_L, s% c_core_v, &
                        s% c_core_omega, s% c_core_omega_div_omega_crit)
                     have_c = .true.
                  end if
               else if (he4 /= 0 .and. c12 /= 0) then
                  if (s% xa(he4,k) <= s% c_core_boundary_he4_fraction .and. &
                      s% xa(c12,k) >= min_x) then
                     call set_core_info(s, k, &
                        s% c_core_mass, s% c_core_radius, s% c_core_lgT, &
                        s% c_core_lgRho, s% c_core_L, s% c_core_v, &
                        s% c_core_omega, s% c_core_omega_div_omega_crit)
                     have_c = .true.
                  end if
               end if
            end if
            
            if (.not. have_he) then
               if (s% he_core_boundary_h1_fraction < 0) then
                  if (A_max >= A_max_he) then
                     call set_core_info(s, k, &
                        s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
                        s% he_core_lgRho, s% he_core_L, s% he_core_v, &
                        s% he_core_omega, s% he_core_omega_div_omega_crit)
                     have_he = .true.
                  end if
               else if (h1 /= 0 .and. he4 /= 0) then
                  if (s% xa(h1,k) <= s% he_core_boundary_h1_fraction .and. &
                      s% xa(he4,k) >= min_x) then
                     call set_core_info(s, k, &
                        s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
                        s% he_core_lgRho, s% he_core_L, s% he_core_v, &
                        s% he_core_omega, s% he_core_omega_div_omega_crit)
                     have_he = .true.
                  end if
               end if
            end if
 
         end do
         
      end subroutine get_core_info


      subroutine clear_core_info( &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit)
         real(dp), intent(out) :: &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit
         core_m = 0
         core_r = 0
         core_lgT = 0
         core_lgRho = 0
         core_L = 0
         core_v = 0
         core_omega = 0
         core_omega_div_omega_crit = 0
      end subroutine clear_core_info
      

      subroutine set_core_info(s, k, &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit
         
         integer :: kcore, j, jm1, j00
         real(dp) :: dm1, d00, qm1, q00, core_q, &
            core_lgP, core_g, core_X, core_Y, core_edv_H, core_edv_He, &
            core_scale_height, core_dlnX_dr, core_dlnY_dr, core_dlnRho_dr
            
         include 'formats'
         
         if (k == 1) then
            core_q = 1d0
         else
            jm1 = maxloc(s% xa(:,k-1), dim=1)
            j00 = maxloc(s% xa(:,k), dim=1)
            qm1 = s% q(k-1) - 0.5d0*s% dq(k-1) ! center of k-1
            q00 = s% q(k) - 0.5d0*s% dq(k) ! center of k
            dm1 = s% xa(j00,k-1) - s% xa(jm1,k-1)
            d00 = s% xa(j00,k) - s% xa(jm1,k)
            if (dm1*d00 > 0d0) then
               write(*,2) 'bad args for set_core_info', k, dm1, d00
               stop 1
               core_q = 05d0*(qm1 + q00)
            else
               core_q = find0(qm1, dm1, q00, d00)
            end if
         end if
            
         call get_info_at_q(s, core_q, &
            kcore, core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_lgP, core_g, core_X, core_Y, core_edv_H, core_edv_He, &
            core_scale_height, core_dlnX_dr, core_dlnY_dr, core_dlnRho_dr, &
            core_omega, core_omega_div_omega_crit)

      end subroutine set_core_info
               


      subroutine get_info_at_q(s, bdy_q, &
            kbdy, bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_lgP, bdy_g, bdy_X, bdy_Y, bdy_edv_H, bdy_edv_He, &
            bdy_scale_height, bdy_dlnX_dr, bdy_dlnY_dr, bdy_dlnRho_dr, &
            bdy_omega, bdy_omega_div_omega_crit)

         type (star_info), pointer :: s
         real(dp), intent(in) :: bdy_q
         integer, intent(out) :: kbdy
         real(dp), intent(out) :: &
            bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_lgP, bdy_g, bdy_X, bdy_Y, bdy_edv_H, bdy_edv_He, &
            bdy_scale_height, bdy_dlnX_dr, bdy_dlnY_dr, bdy_dlnRho_dr, &
            bdy_omega, bdy_omega_div_omega_crit
         
         real(dp) :: x, x0, x1, x2, alfa, beta, bdy_omega_crit
         integer :: k, ii, klo, khi
         
         include 'formats'
         
         bdy_m=0; bdy_r=0; bdy_lgT=0; bdy_lgRho=0; bdy_L=0; bdy_v=0
         bdy_lgP=0; bdy_g=0; bdy_X=0; bdy_Y=0; bdy_edv_H=0; bdy_edv_He=0
         bdy_scale_height=0; bdy_dlnX_dr=0; bdy_dlnY_dr=0; bdy_dlnRho_dr=0
         bdy_omega=0; bdy_omega_div_omega_crit=0
         kbdy = 0

         if (bdy_q <= 0) return
         k = k_for_q(s,bdy_q)
         if (k >= s% nz) return
         if (k <= 1) then
            bdy_m = s% star_mass
            bdy_r = s% r(1)/Rsun
            bdy_lgT = s% lnT(1)/ln10
            bdy_lgRho = s% lnd(1)/ln10
            bdy_L = s% L(1)/Lsun
            if (s% v_flag) bdy_v = s% v(1)
            bdy_lgP = s% lnP(1)/ln10
            bdy_g = s% grav(1)
            bdy_X = s% X(1)
            bdy_Y = s% Y(1)
            if (s% do_element_diffusion) then
            ii = s% net_iso(ih1)
            if (ii > 0) bdy_edv_H = s% edv(ii,1)
            ii = s% net_iso(ihe4)
            if (ii > 0) bdy_edv_He = s% edv(ii,1)
            end if
            bdy_scale_height = s% scale_height(1)
            bdy_omega = s% omega(k)
            bdy_omega_crit = omega_crit(s,1)
            bdy_omega_div_omega_crit = min(1d0,bdy_omega/bdy_omega_crit)
            kbdy = 1
            return
         end if
         
         kbdy = k+1
                  
         bdy_m = (s% M_center + s% xmstar*bdy_q)/Msun
         
         x = s% q(k-1) - bdy_q
         x0 = s% dq(k-1)/2
         x1 = s% dq(k)/2 + s% dq(k-1)
         x2 = s% dq(k+1)/2 + s% dq(k) + s% dq(k-1)
         
         alfa = max(0d0, min(1d0, (bdy_q - s% q(k+1))/s% dq(k)))
         
         bdy_lgT = interp3(s% lnT(k-1), s% lnT(k), s% lnT(k+1))/ln10
         bdy_lgRho = interp3(s% lnd(k-1), s% lnd(k), s% lnd(k+1))/ln10
         bdy_lgP = interp3(s% lnP(k-1), s% lnP(k), s% lnP(k+1))/ln10
         bdy_X = interp3(s% X(k-1), s% X(k), s% X(k+1))
         bdy_Y = interp3(s% Y(k-1), s% Y(k), s% Y(k+1))
                  
         bdy_r = pow_cr( &
            interp2(s% r(k)*s% r(k)*s% r(k), s% r(k+1)*s% r(k+1)*s% r(k+1)),1d0/3d0)/Rsun
         bdy_L = interp2(s% L(k), s% L(k+1))/Lsun         
         bdy_g = interp2(s% grav(k), s% grav(k+1))
         bdy_scale_height = interp2(s% scale_height(k), s% scale_height(k+1))

         if (s% do_element_diffusion) then
            ii = s% net_iso(ih1)
            if (ii > 0) bdy_edv_H = interp2(s% edv(ii,k), s% edv(ii,k+1))
            ii = s% net_iso(ihe4)
            if (ii > 0) bdy_edv_He = interp2(s% edv(ii,k), s% edv(ii,k+1))
         end if

         klo = k-1
         khi = k+1
         bdy_dlnX_dr = log_cr(max(1d-99,max(1d-99,s% X(klo))/max(1d-99,s% X(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
         bdy_dlnY_dr = log_cr(max(1d-99,max(1d-99,s% Y(klo))/max(1d-99,s% Y(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
         bdy_dlnRho_dr = (s% lnd(klo) - s% lnd(khi))/(s% rmid(klo) - s% rmid(khi))
         
         if (s% v_flag) bdy_v = interp2(s% v(k), s% v(k+1))
         if (s% rotation_flag) then
            bdy_omega = interp2(s% omega(k), s% omega(k+1))
            bdy_omega_crit = interp2(omega_crit(s,k), omega_crit(s,k+1))
            bdy_omega_div_omega_crit = min(1d0,bdy_omega/bdy_omega_crit)
         end if
            
         contains
         
         real(dp) function interp3(f0, f1, f2)
            real(dp), intent(in) :: f0, f1, f2
            real(dp) :: fmin, fmax
            fmin = min(f0,f1,f2)
            fmax = max(f0,f1,f2)
            interp3 = (f1*(x-x0)*(x-x2)*(x0-x2)-&
                  (x-x1)*(f0*(x-x2)*(x1-x2) + (x-x0)*(x0-x1)*x2))/ &
                     ((x0-x1)*(x0-x2)*(-x1+x2))
            interp3 = min(fmax, max(fmin, interp3))
         end function interp3
         
         real(dp) function interp2(f0, f1)
            real(dp), intent(in) :: f0, f1
            interp2 = alfa*f0 + (1-alfa)*f1
         end function interp2
         
      end subroutine get_info_at_q

      
      subroutine get_mixing_regions(s, ierr)
         use mlt_def, only: no_mixing
         type (star_info), pointer :: s
         integer, intent(out) :: ierr 
         integer :: prev_type, cur_type, cur_top, n, k
         include 'formats'
         ierr = 0         
         cur_type = s% mixing_type(1)
         cur_top = 1
         n = 0
         do k = 2, s% nz
            prev_type = cur_type
            cur_type = s% mixing_type(k)
            if (cur_type == prev_type .and. k < s% nz) cycle
            ! change of type from k-1 to k
            if (prev_type /= no_mixing) then
               n = n + 1
               s% mixing_region_type(n) = prev_type
               s% mixing_region_top(n) = cur_top
               if (k == s% nz) then
                  s% mixing_region_bottom(n) = k
               else
                  s% mixing_region_bottom(n) = k-1
               end if
               if (n == max_num_mixing_regions) exit
            end if
            cur_top = k
         end do
         
         
         s% num_mixing_regions = n
         
         return
         write(*,2) 's% num_mixing_regions', s% num_mixing_regions
         do n = 1, s% num_mixing_regions
            write(*,'(a20,99i6)') 'mixing region', n, &
               s% mixing_region_type(n), s% mixing_region_top(n), s% mixing_region_bottom(n)
         end do
         
      end subroutine get_mixing_regions
      
      
      subroutine get_power_info(s, cell_masses, ierr)
         use rates_def, only: i_rate, std_reaction_Qs, std_reaction_neuQs
         use chem_def, only: ipp, icno, i3alf
         type (star_info), pointer :: s
         real(dp), pointer, intent(in) :: cell_masses(:) 
         integer, intent(out) :: ierr 
         
         integer :: i, j, k, nz
         real(dp), pointer :: eps_h(:), eps_he(:), eps_z(:)
         
         include 'formats'
         
         ierr = 0
         nz = s% nz
         
         call do_alloc(ierr)
         if (ierr /= 0) return
         
         do k=1,nz
            eps_h(k) = s% eps_nuc_categories(i_rate,ipp,k) + s% eps_nuc_categories(i_rate,icno,k)
            eps_he(k) = s% eps_nuc_categories(i_rate,i3alf,k)
            eps_z(k) = s% eps_nuc(k) - (eps_h(k) + eps_he(k))
         end do
         
         call set_max_burn_info(s, nz, eps_h, &
            s% max_eps_h, s% max_eps_h_lgT, s% max_eps_h_lgRho, s% max_eps_h_m, &
            s% max_eps_h_lgR, s% max_eps_h_lgP, s% max_eps_h_opacity, s% max_eps_h_cp,  &
            s% max_eps_h_k, .false., ierr)
         if (failed('set_max_burn_info h')) return
         
         call set_max_burn_info(s, nz, eps_he, &
            s% max_eps_he, s% max_eps_he_lgT, s% max_eps_he_lgRho, s% max_eps_he_m, &
            s% max_eps_he_lgR, s% max_eps_he_lgP, s% max_eps_he_opacity, s% max_eps_he_cp, &
            s% max_eps_he_k, .false., ierr)
         if (failed('set_max_burn_info he')) return
         
         call set_max_burn_info(s, nz, eps_z, &
            s% max_eps_z, s% max_eps_z_lgT, s% max_eps_z_lgRho, s% max_eps_z_m, &
            s% max_eps_z_lgR, s% max_eps_z_lgP, s% max_eps_z_opacity, s% max_eps_z_cp, &
            s% max_eps_z_k, .false., ierr)
         if (failed('set_max_burn_info z')) return
         
         call check_for_bad_num(s% max_eps_z, 's% max_eps_z')
         call check_for_bad_num(s% max_eps_z_lgT, 's% max_eps_z_lgT')
         call check_for_bad_num(s% max_eps_z_lgRho, 's% max_eps_z_lgRho')
         call check_for_bad_num(s% max_eps_z_m, 's% max_eps_z_m')
         call check_for_bad_num(s% max_eps_z_lgP, 's% max_eps_z_lgP')
         call check_for_bad_num(s% max_eps_z_lgR, 's% max_eps_z_lgR')
         call check_for_bad_num(s% max_eps_z_opacity, 's% max_eps_z_opacity')
         if (ierr /= 0) return
         
         call set_max_burn_info(s, nz, s% eps_nuc, &
            s% max_eps_nuc, s% max_eps_nuc_lgT, s% max_eps_nuc_lgRho, s% max_eps_nuc_m, &
            s% max_eps_nuc_lgR, s% max_eps_nuc_lgP, s% max_eps_nuc_opacity, s% max_eps_nuc_cp, &
            s% max_eps_nuc_k, .false., ierr)
         if (failed('set_max_burn_info eps_nuc')) return         

         call dealloc
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
               call dealloc
            end if
         end function failed
         
         
         subroutine do_alloc(ierr)
            use alloc
            integer, intent(out) :: ierr
            call non_crit_get_work_array(s, eps_h, nz, nz_alloc_extra, 'report', ierr)
            if (ierr /= 0) return            
            call non_crit_get_work_array(s, eps_he, nz, nz_alloc_extra, 'report', ierr)
            if (ierr /= 0) return            
            call non_crit_get_work_array(s, eps_z, nz, nz_alloc_extra, 'report', ierr)
            if (ierr /= 0) return            
         end subroutine do_alloc

         subroutine dealloc
            use alloc
            call non_crit_return_work_array(s, eps_h, 'report')            
            call non_crit_return_work_array(s, eps_he, 'report')            
            call non_crit_return_work_array(s, eps_z, 'report')            
         end subroutine dealloc
         
         subroutine check_for_bad_num(v, str)
            real(dp), intent(in) :: v
            character (len=*), intent(in) :: str
            if (is_bad_num(v)) then
               write(*,*) 'get_power_info ' // trim(str) , v
               ierr = -1
            end if
         end subroutine check_for_bad_num
         
      
      end subroutine get_power_info
      
      
      subroutine set_max_burn_info( &
            s, nz, eps, eps_max, max_lgT, max_lgRho, max_m, max_lgR, max_lgP, &
            max_opacity, max_cp, max_k, show, ierr)
         use interp_1d_def
         use interp_1d_lib
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), intent(in) :: eps(:) ! (nz)
         real(dp), intent(out) :: &
            eps_max, max_lgT, max_lgRho, max_m, max_lgR, max_lgP, max_opacity, max_cp
         logical, intent(in) :: show
         integer, intent(out) :: max_k, ierr 

         integer, parameter :: n_old = 5, n_new = 1, nwork = pm_work_size
         real(dp) :: x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new), &
               delta, delta_2, mx_2, delta_3, mx_3
         integer :: k, i
         logical :: dbg
         real(dp), target :: work_ary(n_old*nwork), f1_ary(4*n_old)
         real(dp), pointer :: work(:), f1(:), f(:,:)
         work => work_ary
         
         include 'formats'
         
         dbg = show
         
         f1 => f1_ary
         f(1:4,1:n_old) => f1(1:4*n_old)
         
         ierr = 0
         k = maxloc(eps(1:nz), dim=1)
         max_k = k
         if (k < 3 .or. k > nz-2) then
            k = nz ! move to center
            call skip_interpolation
            return
         end if
         
         ! place x's at cell midpoints
         x_old(1) = s% dq(k-2)/2
         x_old(2) = s% dq(k-2) + s% dq(k-1)/2
         x_old(3) = s% dq(k-2) + s% dq(k-1) + s% dq(k)/2
         x_old(4) = s% dq(k-2) + s% dq(k-1) + s% dq(k) + s% dq(k+1)/2
         x_old(5) = s% dq(k-2) + s% dq(k-1) + s% dq(k) + s% dq(k+1) + s% dq(k+2)/2
         f(1,1:5) = eps(k-2:k+2)
         call interp_pm(x_old, n_old, f1, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interp_pm for eps'
            return
         end if
         
         if (dbg) write(*,*) 'for 2'
         call find_delta_and_max(2,delta_2,mx_2,ierr)
         if (ierr /= 0) then
            if (dbg) then
               write(*,2) 'mx_2', k, mx_2
               do i=k-2,k+2
                  write(*,2) 'eps', i, eps(i), s% dq(i)
               end do
               write(*,*)
               do i=1,5
                  write(*,2) 'x_old', i, x_old(i)
               end do
               write(*,*)
               stop 'debug: set_max_burn_info'
            end if
            write(*,*) 'set_max_burn_info: failed in find_delta_and_max 2'
            return
         end if
         if (dbg) then
            write(*,1) 'fraction', delta_2/(x_old(3)-x_old(2))
            write(*,1) 'mx_2', mx_2
            write(*,*)
         end if

         if (dbg) write(*,*) 'for 3'
         call find_delta_and_max(3,delta_3,mx_3,ierr)
         if (ierr /= 0) then
            if (dbg) then
               write(*,2) 'mx_3', k, mx_3
               do i=k-2,k+2
                  write(*,2) 'eps', i, eps(i), s% dq(i)
               end do
               do i=1,5
                  write(*,2) 'x_old', i, x_old(i)
               end do
               stop 'debug: set_max_burn_info'
            end if
            write(*,*) 'set_max_burn_info: failed in find_delta_and_max 3'
            return
         end if
         if (dbg) then
            write(*,1) 'fraction', delta_3/(x_old(4)-x_old(3))
            write(*,1) 'mx_3', mx_3
            write(*,*)
         end if

         if (mx_2 > mx_3 .and. delta_2 < x_old(3)-x_old(2)) then
            delta = delta_2; eps_max = mx_2; x_new(1) = x_old(2) + delta
            max_m = s% M_center + s% xmstar*(s% q(k-2) - (x_old(2)+delta))
            if (dbg) write(*,*) 'use 2', delta/(x_old(3)-x_old(2))
         else
            delta = delta_3; eps_max = mx_3; x_new(1) = x_old(3) + delta
            max_m = s% M_center + s% xmstar*(s% q(k-2) - (x_old(3)+delta))
            if (dbg) write(*,*) 'use 3', delta/(x_old(4)-x_old(3))
         end if
                  
         ! interpolate max_lgT, max_lgRho, max_lgP, max_opacity, max_cp
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% lnT(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for lnT'
            return
         end if
         max_lgT = v_new(1)/ln10
         
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% cp(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for cp'
            return
         end if
         max_cp = v_new(1)
         
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% lnd(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for lnd'
            return
         end if
         max_lgRho = v_new(1)/ln10
         
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% lnP(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for lnP'
            return
         end if
         max_lgP = v_new(1)/ln10
         
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% opacity(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for opacity'
            return
         end if
         max_opacity = v_new(1)

         ! change to interpolating values at cell boundaries
         x_new(1) = (max_m - s% M_center)/s% xmstar ! q at max
         x_old(1:5) = s% q(k-2:k+2)
         ! interpolate max_lgR
         call interpolate_vector( &
            n_old, x_old, n_new, x_new, s% lnr(k-2:k+2), v_new, interp_pm, nwork, work, &
            'report set_max_burn_info', ierr)
         if (ierr /= 0) then
            write(*,*) 'set_max_burn_info: failed in interpolate_vector for lgR'
            return
         end if
         max_lgR = v_new(1)/ln10
         
         if (dbg) then
            write(*,*) 'k', k
            write(*,1) 'x_old(1:5)', x_old(1:5)
            write(*,1) 'mass(k-2:k+2)', s% M_center + s% xmstar*(s% q(k-2:k+2) + s% q(k-1:k+3))/2
            write(*,1) 'max_m', max_m
            write(*,1) 's% M_center', s% M_center
            write(*,1) 's% xmstar', s% xmstar
            write(*,1) 'q(k-2:k+2)', s% q(k-2:k+2)
            write(*,1) 'dq(k-2:k+2)', s% dq(k-2:k+2)
            write(*,1) 'm(k-2:k+2)', s% m(k-2:k+2)
            write(*,*)
            write(*,1) 'eps(k-2:k+2)', eps(k-2:k+2)
            write(*,1) 'eps_max', eps_max
            write(*,*)
            write(*,1) 'cp(k-2:k+2)', s% cp(k-2:k+2)/ln10
            write(*,1) 'max_cp', max_cp
            write(*,*)
            write(*,1) 'lgT(k-2:k+2)', s% lnT(k-2:k+2)/ln10
            write(*,1) 'max_lgT', max_lgT
            write(*,*)
            write(*,1) 'lgd(k-2:k+2)', s% lnd(k-2:k+2)/ln10
            write(*,1) 'max_lgRho', max_lgRho
            write(*,*)
            write(*,1) 'lgR(k-2:k+2)', s% lnR(k-2:k+2)/ln10
            write(*,1) 'max_lgR', max_lgR
            write(*,*)
            if (is_bad_num(max_lgR)) stop
         end if 
         
         contains         
         
         subroutine find_delta_and_max(i,delta,mx,ierr)
            integer, intent(in) :: i
            real(dp), intent(out) :: delta, mx
            integer, intent(out) :: ierr
            real(dp) :: s2, s, delta1, delta2, v1, v2
            include 'formats'
            ierr = 0
            if (abs(f(4,i)) < 1d-99) then ! treat as quadratic
               call quad(i,delta,mx)
               if (is_bad_num(mx)) ierr = -1
               return
            end if
            s2 = f(3,i)*f(3,i) - 3*f(2,i)*f(4,i)
            if (s2 < 0 .or. f(4,i) == 0 .or. f(3,i) == 0) then
               call quad(i,delta,mx)
               if (is_bad_num(mx)) ierr = -1
               return
            end if
            s = sqrt(s2)
            delta1 = (-f(3,i) + s)/(3*f(4,i))
            delta2 = -(f(3,i) + s)/(3*f(4,i))
            v1 = f(1,i) + delta1*(f(2,i) + delta1*(f(3,i) + delta1*f(4,i)))
            v2 = f(1,i) + delta2*(f(2,i) + delta2*(f(3,i) + delta2*f(4,i)))
            if (v1 > v2) then
               delta = delta1; mx = v1
            else
               delta = delta2; mx = v2
            end if
            if (is_bad_num(mx)) then
               call quad(i,delta,mx)
               if (is_bad_num(mx)) then
                  ierr = -1
                  if (dbg) write(*,1) 'find_delta_and_max mx', mx
               end if
            end if
         end subroutine find_delta_and_max
         
         
         subroutine quad(i,delta,mx)
            integer, intent(in) :: i
            real(dp), intent(out) :: delta, mx
            if (f(3,i) == 0) then
               delta = 0; mx = f(1,i); return
            end if
            delta = -f(2,i)/(2*f(3,i))
            if (delta < 0 .or. delta + x_old(i) > x_old(i+1)) then ! no extremum
               if (f(1,i) > f(1,i+1)) then
                  delta = 0; mx = f(1,i)
               else
                  delta = x_old(i+1)-x_old(i); mx = f(1,i+1)
               end if
            else
               mx = f(1,i) + delta*(f(2,i) + delta*f(3,i))
            end if
         end subroutine quad
         
         
         subroutine skip_interpolation
            eps_max = eps(k)
            max_lgT = s% lnT(k)/ln10
            max_lgRho = s% lnd(k)/ln10
            max_lgP = s% lnp(k)/ln10
            max_opacity = s% opacity(k)
            max_lgR = s% lnR(k)/ln10
            max_cp = s% cp(k)
            if (k < nz) then
               max_m = s% M_center + s% xmstar*(s% q(k) + s% q(k+1))/2/Msun ! midpoint of cell k
            else
               max_m = s% M_center
            end if
         end subroutine skip_interpolation


      end subroutine set_max_burn_info
      
      
      subroutine get_neutrino_fluxes(s, ierr)
         use rates_def
         use net_lib, only: get_net_reaction_table_ptr
         type (star_info), pointer :: s
         integer, intent(out) :: ierr 
      
         real(dp) :: Lpp, Lpep, Lbe7ec, Lb8ep
         real(dp) :: Lpp_neu, Lpep_neu, Lbe7ec_neu, Lb8ep_neu
         real(dp) :: npp, npep, nbe7ec, nb8ep
         integer, pointer :: net_reaction_ptr(:) ! maps reaction id to net reaction number
         integer :: r_h1_h1_wk_h2, r_h1_h1_ec_h2, r_be7_ec_li7, r_b8_wk_he4_he4, nz
         
         include 'formats'
         
         ierr = 0
         if (.not. s% trace_solar_neutrinos) return
         
         nz = s% nz
         
         call get_net_reaction_table_ptr(s% net_handle, net_reaction_ptr, ierr)
         if (ierr /= 0) then
            write(*,*) 'report get_neutrino_fluxes failed in get_net_reaction_table_ptr'
            return
         end if
         
         r_h1_h1_wk_h2 = net_reaction_ptr(ir_h1_h1_wk_h2)
         r_h1_h1_ec_h2 = net_reaction_ptr(ir_h1_h1_ec_h2)
         r_be7_ec_li7 = net_reaction_ptr(ir_be7_ec_li7)
         r_b8_wk_he4_he4 = net_reaction_ptr(ir_b8_wk_he4_he4)
         
         if (r_h1_h1_wk_h2 == 0) then
            r_h1_h1_wk_h2 = net_reaction_ptr(irpp_to_he3)
            write(*,*) 'use irpp_to_he3'
         else
            write(*,*) 'use irpp'
         end if
         if (r_h1_h1_ec_h2 == 0) then
            r_h1_h1_ec_h2 = net_reaction_ptr(irpep_to_he3)
            write(*,*) 'use irpep_to_he3'
         else
            write(*,*) 'use irpep'
         end if
         if (r_be7_ec_li7 == 0) then
            r_be7_ec_li7 = net_reaction_ptr(ir34_pp2)
            write(*,*) 'use ir34_pp2'
         else
            write(*,*) 'use ir_be7_ec_li7'
         end if
         if (r_b8_wk_he4_he4 == 0) then
            r_b8_wk_he4_he4 = net_reaction_ptr(ir34_pp3)
            write(*,*) 'use ir34_pp3'
         else
            write(*,*) 'use irb8ep'
         end if
         
         call do1(r_h1_h1_wk_h2, ir_h1_h1_wk_h2, Lpp, npp, Lpp_neu, s% flux_pp)
         call do1(r_h1_h1_wk_h2, ir_h1_h1_wk_h2, Lpep, npep, Lpep_neu, s% flux_pep)
         call do1(r_be7_ec_li7, ir_be7_ec_li7, Lbe7ec, nbe7ec, Lbe7ec_neu, s% flux_be7ec)
         call do1(r_b8_wk_he4_he4, ir_b8_wk_he4_he4, Lb8ep, nb8ep, Lb8ep_neu, s% flux_b8ep)
         
         write(*,1) 'r_h1_h1_wk_h2', Lpp/Lsun, npp, Lpp_neu, s% flux_pp
         write(*,1) 'r_h1_h1_ec_h2', Lpep/Lsun, npep, Lpep_neu, s% flux_pep
         write(*,1) 'r_be7_ec_li7', Lbe7ec/Lsun, nbe7ec, Lbe7ec_neu, s% flux_be7ec
         write(*,1) 'r_b8_wk_he4_he4', Lb8ep/Lsun, nb8ep, Lb8ep_neu, s% flux_b8ep
         write(*,*)
         
         contains
         
         subroutine do1(r, ir, L, n, L_neu, flux)
            use rates_def, only: i_rate, std_reaction_Qs, std_reaction_neuQs
            integer, intent(in) :: r, ir
            real(dp), intent(out) :: L, n, L_neu, flux
            real(dp) :: eps
            if (r == 0) then
               L = 0
            else
               L = dot_product(s% dm(1:nz), s% reaction_eps_nuc(i_rate, r, 1:nz))
            end if
            eps = ev2erg*1d6*(std_reaction_Qs(ir)-std_reaction_neuQs(ir)) ! ergs per reaction
            n = L/eps ! number per second produced
            flux = n/(4*pi*au*au) ! number per cm^2 per second at earth
            L_neu = n*ev2erg*1d6*std_reaction_neuQs(ir)
         end subroutine do1
         
         
      end subroutine get_neutrino_fluxes
      

      end module report
