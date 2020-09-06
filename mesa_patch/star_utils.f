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

      module star_utils

      use star_private_def
      use const_def
      use num_lib

      implicit none

      
      contains
      
      
      subroutine foreach_cell(s,nzlo,nzhi,use_omp,do1,ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: use_omp
         interface
            subroutine do1(s,k,ierr)
               use star_private_def
               type (star_info), pointer :: s         
               integer, intent(in) :: k
               integer, intent(out) :: ierr
            end subroutine do1
         end interface
         integer, intent(out) :: ierr
         
         integer :: k, op_err
         logical :: okay
         ierr = 0
         
         if (nzlo == nzhi) then
            call do1(s,nzlo,ierr)
            return
         end if
         
         if (use_omp) then
            okay = .true.
!xe$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(STATIC,10)
!$OMP PARALLEL DO PRIVATE(k,op_err)
            do k = nzlo, nzhi
               if (.not. okay) cycle
               op_err = 0
               call do1(s,k,op_err)
               if (op_err /= 0) okay = .false. ! cannot just exit from a parallel loop
            end do
!$OMP END PARALLEL DO
            if (.not. okay) ierr = -1
         else
            do k = nzlo, nzhi
               call do1(s,k,ierr)
               if (ierr /= 0) exit
            end do
         end if
      
      end subroutine foreach_cell
      
      
      subroutine foreach_cell_dynamic(s,nzlo,nzhi,use_omp,do1,ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: use_omp
         interface
            subroutine do1(s,k,ierr)
               use star_private_def
               type (star_info), pointer :: s         
               integer, intent(in) :: k
               integer, intent(out) :: ierr
            end subroutine do1
         end interface
         integer, intent(out) :: ierr
         
         integer :: k, op_err
         ierr = 0
         
         if (nzlo == nzhi) then
            call do1(s,nzlo,ierr)
            return
         end if
         
         if (use_omp) then
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(DYNAMIC,10)
            do k = nzlo, nzhi
               if (ierr /= 0) cycle
               op_err = 0
               call do1(s,k,op_err)
               if (op_err /= 0) ierr = op_err
            end do
!$OMP END PARALLEL DO
         else
            do k = nzlo, nzhi
               call do1(s,k,ierr)
               if (ierr /= 0) exit
            end do
         end if
      
      end subroutine foreach_cell_dynamic


      real(dp) function sum_Egrav(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: dq_prev, dq_cur
         dq_cur = 0
         sum_Egrav = 0
         do k = 1, s% nz-1
            dq_prev = dq_cur
            dq_cur = s% dq(k)
            sum_Egrav = sum_Egrav + 0.5d0*(dq_prev + dq_cur)*s% q(k)/s% r(k)
         end do
         k = s% nz
         sum_Egrav = sum_Egrav + (0.5d0*s% dq(k-1) + s% dq(k))*s% q(k)/s% r(k)
         sum_Egrav = -sum_Egrav*s% cgrav(k)*s% mstar*s% mstar
      end function sum_Egrav


      real(dp) function sum_Etherm(s)
         type (star_info), pointer :: s
         integer :: nz, k
         nz = s% nz
         sum_Etherm = 0d0
         do k=1,nz
            sum_Etherm = sum_Etherm + s% dm(k)*exp_cr(s% lnE(k))
         enddo
      end function sum_Etherm


      real(dp) function sum_Ebinding(s)
         use chem_def
         type (star_info), pointer :: s
         integer :: cid, j
         real(dp) :: m1, n1, E1, Etotal
         integer :: nz
         nz = s% nz
         Etotal = 0
         do j=1, s% species
            cid = s% chem_id(j)
            m1 = dot_product(s% xa(j,1:nz),s% dm(1:nz)) ! grams of species j
            n1 = m1*avo/chem_isos% W (cid) ! number of species j nuclei
            E1 = ev2erg*1d6*chem_isos% binding_energy(cid) ! ergs binding energy per nuclei of species j
            Etotal = Etotal + E1*n1
         end do
         sum_Ebinding = Etotal
      end function sum_Ebinding
      
      
      real(dp) function sum_L_nuc(s)
         type (star_info), pointer :: s
         integer :: nz
         nz = s% nz
         sum_L_nuc = dot_product(s% dm(1:nz), s% eps_nuc(1:nz))
      end function sum_L_nuc
      
      
      real(dp) function sum_L_grav(s)
         type (star_info), pointer :: s
         integer :: nz
         nz = s% nz
         sum_L_grav = s% eps_grav_factor*dot_product(s% dm(1:nz), s% eps_grav(1:nz))
      end function sum_L_grav
      
      
      subroutine get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         use chem_def
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(out) :: y_avg, z_avg
         integer, intent(out) :: ierr
         
         integer :: k, nz,  h1, h2, he3, he4
         real(dp) :: total_mass_h, total_mass_he, total_mass_z, &
            cell_mass, total_mass
         
         ierr = 0
         
         nz = s% nz
         h1 = s% net_iso(ih1)
         h2 = s% net_iso(ih2)
         he3 = s% net_iso(ihe3)
         he4 = s% net_iso(ihe4)
         total_mass=0; total_mass_h=0; total_mass_he=0; total_mass_z=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            total_mass_h = total_mass_h + cell_mass*s% xa(h1, k)
            if (h2 /= 0) total_mass_h = total_mass_h + cell_mass*s% xa(h2, k)
            total_mass_he = total_mass_he + cell_mass*s% xa(he4, k)
            if (he3 /= 0) total_mass_he = total_mass_he + cell_mass*s% xa(he3, k)
         end do
         total_mass_z = total_mass - (total_mass_h + total_mass_he)
         z_avg = total_mass_z / total_mass
         y_avg = total_mass_he / total_mass
         
      end subroutine get_average_Y_and_Z


      real(dp) function eval_current_y(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_y = y_avg
      end function eval_current_y


      real(dp) function eval_current_z(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_z = z_avg
      end function eval_current_z


      real(dp) function eval_current_abundance(s, j, nzlo, nzhi, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: j, nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: cell_mass, jmass, total_mass
         
         ierr = 0
         
         if (j == 0) then
            eval_current_abundance = 0
            return
         end if
         
         nz = s% nz
         total_mass=0; jmass=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            jmass = jmass + cell_mass*s% xa(j, k)
         end do
         eval_current_abundance = jmass / total_mass
         
      end function eval_current_abundance

      
      subroutine smooth_abundances(s, cnt, nzlo, nzhi, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: cnt ! make this many passes
         integer, intent(in) :: nzlo, nzhi ! only smooth zones nzlo to nzhi inclusive
         integer, intent(out) :: ierr
         integer :: k, j, nz
         ierr = 0
         nz = s% nz
         do j = 1, cnt
            do k = max(nzlo,2), min(nzhi, nz)
               s% xa(:,k) = (s% xa(:,k-1) + s% xa(:,k) + s% xa(:,k+1))/3
            end do
            if (nzhi == nz) s% xa(:,nz) = (s% xa(:,nz-1) + s% xa(:,nz) + s% xa(:,nz))/3
            if (nzlo == 1) s% xa(:,1) = (s% xa(:,2) + s% xa(:,1) + s% xa(:,1))/3
         end do
      end subroutine smooth_abundances
      
      
      integer function k_for_q(s, q)
         ! return k s.t. q(k) >= q > q(k)-dq(k)
         type (star_info), pointer :: s         
         real(dp), intent(in) :: q
         integer :: k, nz
         nz = s% nz
         if (q >= 1) then
            k_for_q = 1; return
         else if (q <= s% q(nz)) then
            k_for_q = nz; return
         end if
         do k = 1, nz-1
            if (q > s% q(k+1)) then
               k_for_q = k; return
            end if
         end do
         k_for_q = nz
      end function k_for_q
  
  
      subroutine get_name_for_restart_file(n, num_digits, num)
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         call get_string_for_model_number('x', n, num_digits, num)
      end subroutine get_name_for_restart_file
  
  
      subroutine get_string_for_model_number(prefix, n, num_digits, num)
         character (len=*), intent(in) :: prefix
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         integer :: val
         character (len=32) :: fstring
         include 'formats'
         val = mod(n, 10**num_digits) ! wrap around
         if (val == 0) then
            write(num,*) n
            num = adjustl(num)
            return
         end if
        write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits
        write(num,fstring) trim(prefix), val
      end subroutine get_string_for_model_number
         
         
      subroutine report_xa_bad_nums(s,ierr)         
         use utils_lib, only: is_bad_num
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, j
         ierr = 0
         do k=1,s% nz
            do j=1,s% species
               if (is_bad_num(s% xa(j,k))) then
                  ierr = -1
                  write(*,*) j, k, s% xa(j,k)
               end if
            end do
         end do
      end subroutine report_xa_bad_nums
      
      
      real(dp) function eval_csound(s,k,ierr) result(cs)
         use utils_lib, only: is_bad_num
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (s% use_sr_sound_speed) then
            cs = clight*sqrt( &
               s% gamma1(k)/(1d0 + (exp_cr(s% lnE(k)) + clight*clight)*s% rho(k)/s% P(k)))
         else
            cs = sqrt(s% gamma1(k)*s% P(k)/s% rho(k))
         end if
         if (is_bad_num(cs)) then
            if (s% report_ierr) &
               write(*,2) 'bad csound', k, cs, s% gamma1(k), s% P(k), s% rho(k)
            ierr = -1
         end if
      end function eval_csound
      
      
      subroutine set_m_grav_and_grav(s) ! using mass_corrections
         type (star_info), pointer :: s
         integer :: k, nz
         real(dp) :: twoGmrc2
         include 'formats'
         nz = s% nz
         if (.not. s% use_mass_corrections) then
            do k=1,nz
               s% m_grav(k) = s% m(k)
            end do
         else
            s% m_grav(nz) = &
               s% M_center + s% dm(nz)*s% mass_correction(nz)
            do k=nz-1,1,-1
               s% m_grav(k) = &
                  s% m_grav(k+1) + s% dm(k)*s% mass_correction(k)
            end do
         end if
         do k=1,nz
            s% grav(k) = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
            if (s% use_gr_factors) then ! GR gravity factor = 1/sqrt(1-2Gm/(rc^2))
               twoGmrc2 = 2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight)
               s% grav(k) = s% grav(k)/sqrt(1d0 - twoGmrc2)
            end if
         end do
      end subroutine set_m_grav_and_grav
      
      
      subroutine set_m_and_dm(s)
         type (star_info), pointer :: s
         integer :: k
         do k = 1, s% nz
            s% m(k) = s% M_center + s% q(k)*s% xmstar
            s% dm(k) = s% dq(k)*s% xmstar
         end do
      end subroutine set_m_and_dm
      
      
      subroutine set_dm_bar(nz, dm, dm_bar)
         integer, intent(in) :: nz
         real(dp), intent(in) :: dm(:) ! (nz)
         real(dp), intent(out) :: dm_bar(:) ! (nz)
         integer :: k
         do k=2,nz-1
            dm_bar(k) = 0.5d0*(dm(k-1) + dm(k))
         end do
         dm_bar(1) = 0.5d0*dm(1)
         dm_bar(nz) = 0.5d0*dm(nz-1) + dm(nz)
      end subroutine set_dm_bar
      
      
      subroutine normalize_dqs(nz, dq, ierr) 
         ! rescale dq's so that add to 1.000
         ! work in from boundaries to meet at largest dq
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k, midq
         real(dp) :: dqsum1, dqsum2
         include 'formats'
         midq = maxloc(dq(1:nz),dim=1)
         ! surface inward
         dqsum1 = 0
         do k=1, midq
            dqsum1 = dqsum1 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               !write(*,2) 'normalize_dqs: bad dq(k)', k, dq(k)
               return
            end if
         end do
         ! center outward
         dqsum2 = 0
         do k=nz, midq+1, -1
            dqsum2 = dqsum2 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               !write(*,2) 'normalize_dqs: bad dq(k)', k, dq(k)
               return
            end if
         end do
         !write(*,1) 'normalize_dqs: dqsum1+dqsum2', dqsum1+dqsum2
         dq(1:nz) = dq(1:nz)/(dqsum1 + dqsum2)
      end subroutine normalize_dqs
      
      
      subroutine set_qs(nz, q, dq, ierr) ! set q's using dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(out) :: q(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k   
         include 'formats'      
         ierr = 0
         q(1) = 1
         do k=2,nz-1
            q(k) = q(k-1) - dq(k-1)
            if (q(k) < 0d0 .or. q(k) > 1d0) then
               write(*,2) 'set_qs q', k-1, q(k-1)
               write(*,2) 'set_qs dq', k-1, dq(k-1)
               write(*,2) 'set_qs q-dq', k-1, q(k-1) - dq(k-1)
               write(*,2) 'set_qs q', k, q(k)
               !stop
               ierr = -1
               return
            end if
         end do
         q(nz) = dq(nz)
         if (q(nz) >= q(nz-1)) then
            q(nz) = q(nz-1) - dq(nz-1)
            dq(nz) = q(nz)
            if (dq(nz) <= 0) then
               write(*,2) 'set_qs q', nz-1, q(nz-1)
               write(*,2) 'set_qs dq', nz-1, dq(nz-1)
               write(*,2) 'set_qs q-dq', nz-1, q(nz-1) - dq(nz-1)
               write(*,2) 'set_qs q', nz, q(nz)
               !stop
               ierr = -1
               return
            end if
         end if
      end subroutine set_qs
      
      
      subroutine set_xqs(nz, xq, dq, ierr) ! set xq's using dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(out) :: xq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k
         include 'formats'      
         ierr = 0
         xq(1) = 0
         do k=2,nz-1
            xq(k) = xq(k-1) + dq(k-1)
         end do
         xq(nz) = 1 - dq(nz)
         if (xq(nz) < xq(nz-1)) then
            xq(nz) = xq(nz-1) + dq(nz-1)
            dq(nz) = 1 - xq(nz)
            if (dq(nz) <= 0) then
               ierr = -1
               return
            end if
         end if
      end subroutine set_xqs

         
      real(dp) function interp_val_to_pt(v,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         integer, intent(in) :: k, sz
         real(dp), pointer :: v(:), dq(:)
         character (len=*), intent(in) :: str
         integer :: ierr
         include 'formats'
         if (k == 1) then
            interp_val_to_pt = v(k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               v(k-2), v(k-1), v(k), v(k+1), &
               interp_val_to_pt, str, ierr)
            if (ierr == 0) return
            write(*,1) '0.5d0*(dq(k-2)+dq(k-1))', 0.5d0*(dq(k-2)+dq(k-1))
            write(*,1) '0.5d0*(dq(k-1)+dq(k))', 0.5d0*(dq(k-1)+dq(k))
            write(*,1) '0.5d0*(dq(k)+dq(k+1))', 0.5d0*(dq(k)+dq(k+1))
            write(*,2) 'dq(k-2)', k-2, dq(k-2)
            write(*,2) 'dq(k-1)', k-1, dq(k-1)
            write(*,2) 'dq(k)', k, dq(k)
            write(*,2) 'dq(k+1)', k+1, dq(k+1)
            
            stop 'interp_val_to_pt'
         endif
         interp_val_to_pt = (v(k)*dq(k-1) + v(k-1)*dq(k))/(dq(k-1) + dq(k))
      end function interp_val_to_pt
      
      
      real(dp) function interp_xa_to_pt(xa,j,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         real(dp), pointer :: xa(:,:), dq(:)
         character (len=*), intent(in) :: str
         integer, intent(in) :: j, k, sz
         integer :: ierr
         include 'formats'
         if (j == 0) then
            interp_xa_to_pt = 0
            return
         end if
         if (k == 1) then
            interp_xa_to_pt = xa(j,k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               xa(j,k-2), xa(j,k-1), xa(j,k), xa(j,k+1), &
               interp_xa_to_pt, str, ierr)
            if (ierr == 0) return
         endif
         interp_xa_to_pt = (xa(j,k)*dq(k-1) + xa(j,k-1)*dq(k))/(dq(k-1) + dq(k))
      end function interp_xa_to_pt

      
      real(dp) function get_dtau1(s)
         type (star_info), pointer :: s  
         get_dtau1 = s% dm(1)*s% opacity(1)/(4*pi*s% rmid(1)*s% rmid(1))
      end function get_dtau1
      
      
      subroutine get_tau(s, tau)
         type (star_info), pointer :: s  
         real(dp), pointer :: tau(:) 
         ! tau(k) is optical depth at outer boundary of cell k
         real(dp) :: dtau, dr
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         dtau = get_dtau1(s)
         tau(1) = s% tau_factor*s% tau_base 
         do k = 2, s% nz
            tau(k) = tau(k-1) + dtau
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
         end do
      end subroutine get_tau
         
         
      integer function find_cell_for_mass(s, m)
         type (star_info), pointer :: s
         real(dp), intent(in) :: m
         integer :: k
         find_cell_for_mass = s% nz
         do k = 1, s% nz-1
            if (s% m(k) >= m .and. m > s% m(k+1)) then
               find_cell_for_mass = k
               return
            end if
         end do
      end function find_cell_for_mass
      
      
      subroutine get_delta_Pg(s, nu_max, delta_Pg)
         use utils_lib, only: is_bad_num
         type (star_info), pointer :: s  
         real(dp), intent(in) :: nu_max ! microHz
         real(dp), intent(out) :: delta_Pg ! seconds
         ! g-mode period spacing for l=1
         real(dp) :: integral, N2, omega2, kr2, L2, el, &
            dr, r, r2, cs2, sl2, I_integral, I_integral_limit
         integer :: k, k_sl2
         logical, parameter :: dbg = .false.
         include 'formats'
         if (dbg) then
            write(*,2) 'nu_max', s% model_number, nu_max
            write(*,2) 's% star_mass', s% model_number, s% star_mass
            write(*,2) 's% photosphere_r', s% model_number, s% photosphere_r
            write(*,2) 's% Teff', s% model_number, s% Teff
         end if
         delta_Pg = 0
         integral = 0
         I_integral = 0
         I_integral_limit = 0.5
         omega2 = pow2(2*pi*nu_max/1d6)
         if (dbg) write(*,1) 'log omega2', log10_cr(omega2)
         el = 1
         L2 = el*(el+1)
         k_sl2 = 0
         do k = 2, s% nz
            N2 = s% brunt_N2(k)
            r = s% r(k)
            r2 = r*r
            cs2 = s% csound_at_face(k)*s% csound_at_face(k)
            sl2 = L2*cs2/r2
            dr = s% rmid(k-1) - s% rmid(k)
            !if (N2 > 0) then
            !   write(*,*) 'omega2 >= sl2 omega2 < N2', omega2 >= sl2, omega2 < N2
            !   write(*,2) '      omega2, sl2, N2', k, omega2, sl2, N2
            !end if
            if (omega2 >= sl2) then
               !write(*,3) 'omega2 >= sl2', k, s% nz, omega2, sl2, cs2, r2
               cycle
            end if
            if (k_sl2 == 0) then
               k_sl2 = k
               if (dbg) write(*,2) 'k_sl2', k
            end if
            !if (dbg) write(*,2) 'N2 - omega2', k, N2 - omega2, N2, omega2
            if (N2 > omega2) then ! in g-cavity
               if (dbg .and. integral == 0) write(*,2) 'enter g-cavity', k
               integral = integral + sqrt(N2)*dr/r
               !write(*,3) 'integral', k, s% nz, integral, sqrt(N2)*dr/r, N2, dr, r
            else ! in decay region
               if (integral == 0) cycle ! ! haven't been in g-cavity yet
               if (dbg .and. I_integral == 0) write(*,2) 'enter decay', k
               !write(*,3) 'omega2 < N2', k, s% nz, omega2 - N2, omega2, N2
               ! in decay region below g-cavity; I_integral estimates decay
               kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)*omega2/cs2
               I_integral = I_integral + sqrt(-kr2)*dr
               if (I_integral > I_integral_limit) exit
            end if
         end do
         
         if (dbg) write(*,2) 'omega2 nu_max integral I_integral', &
            s% model_number, omega2, nu_max, integral, I_integral

         if (integral == 0) return
         delta_Pg = sqrt(2d0)*pi*pi/integral
         if (is_bad_num(delta_Pg)) delta_Pg = 0
         
         if (dbg) write(*,2) 'delta_Pg', s% model_number, delta_Pg
         !stop 'get_delta_Pg'
         
      end subroutine get_delta_Pg
      
      
      real(dp) function get_fraction_NSE_burn(s, k) result(alfa)
         use chem_def, only: isi28
         type (star_info), pointer :: s  
         integer, intent(in) :: k
         integer :: si28
         if (s% T(k) <= s% T_NSE_full_off) then
            alfa = 0d0
            return
         end if
         if (s% T(k) >= s% T_NSE_full_on) then               
            alfa = 1d0
            return
         end if
!         si28 = s% net_iso(isi28)
!         if (si28 > 0 .and. s% si28_NSE_full_off > 0 .and. &
!               s% si28_NSE_full_off > s% si28_NSE_full_off) then
!            alfa = (s% xa(si28,k) - s% si28_NSE_full_off)/ &
!               (s% si28_NSE_full_on - s% si28_NSE_full_off)
!            alfa = max(0d0, min(1d0, alfa))
!         else
!            alfa = 1d0
!         end if
         alfa = (s% T(k) - s% T_NSE_full_off)/ &
            (s% T_NSE_full_on - s% T_NSE_full_off)
      end function get_fraction_NSE_burn
      
      
      real(dp) function get_tau_at_r(s, r)
         type (star_info), pointer :: s  
         real(dp), intent(in) :: r
         real(dp) :: dtau, dr, tau_m1, tau_00
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         dtau = get_dtau1(s)
         tau_00 = s% tau_factor*s% tau_base 
         get_tau_at_r = tau_00
         if (r >= s% r(1)) return
         do k = 2, s% nz
            tau_m1 = tau_00
            tau_00 = tau_m1 + dtau
            if (r < s% r(k-1) .and. r >= s% r(k)) then
               get_tau_at_r = &
                  (tau_00*(s% r(k-1)-r) + tau_m1*(r-s% r(k)))/(s% r(k-1)-s% r(k))
               return
            end if
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
         end do
      end function get_tau_at_r
      
      
      integer function find_tau_phot(s, tau00, taup1, ierr) 
         ! return k for the cell containing optical depth = 2/3
         type (star_info), pointer :: s  
         real(dp), intent(out) :: tau00, taup1
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dtau
         
         real(dp), parameter :: tau_phot = 2d0/3d0
         
         include 'formats'
         ierr = 0
         tau00 = 0
         taup1 = 0
         find_tau_phot = 1
         if (s% tau_factor >= 1) return
         tau00 = s% tau_factor*s% tau_base
         do k = 1, s% nz
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot) then
               find_tau_phot = k
               return
            end if
            tau00 = taup1
         end do
         ierr = -1
      end function find_tau_phot
      
      
      real(dp) function get_r_phot(s) ! return r where optical depth = 2/3
         type (star_info), pointer :: s  

         integer :: k
         real(dp) :: tau00, taup1, dtau, r003, rp13, r3
         
         real(dp), parameter :: tau_phot = 2d0/3d0
         
         include 'formats'

         tau00 = 0
         taup1 = 0
         get_r_phot = s% r(1)
         if (s% tau_factor >= 1) return
         tau00 = s% tau_factor*s% tau_base
         if (tau00 >= tau_phot) return
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               r003 = s% r(k)*s% r(k)*s% r(k)
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
               r3 = r003 + (rp13 - r003)*(tau_phot - tau00)/dtau
               get_r_phot = pow_cr(r3,1d0/3d0)
               return
            end if
            tau00 = taup1
         end do

      end function get_r_phot


      real(dp) function center_value(s, p)
         type (star_info), pointer :: s  
         real(dp), intent(in) :: p(:)
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = p(k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x + dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_value = sum_x/sum_dq
      end function center_value         

            
      subroutine interp_q( &
            nz2, nvar_hydro, species, qval, xh, xa, q, dq, struct, comp, ierr)
         use num_lib, only: binary_search
         integer, intent(in) :: nz2, nvar_hydro, species
         !real(dp), intent(in) :: qval, xh(nvar_hydro,nz2), xa(species,nz2), q(nz2), dq(nz2)
         real(dp), intent(in) :: qval
         real(dp), intent(in), pointer :: xh(:,:), xa(:,:), q(:), dq(:)
         !real(dp), intent(out) :: struct(nvar_hydro), comp(species)
         real(dp), intent(out) :: struct(:), comp(:)
         integer, intent(out) :: ierr         
         integer :: k
         real(dp) :: alfa         
         ierr = 0         
         if (qval <= q(nz2)) then
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = xh(1:nvar_hydro,nz2)
            if (species > 0) &
               comp(1:species) = xa(1:species,nz2)
            return
         end if
         k = binary_search(nz2, q, 0, qval)
         if (k < 1 .or. k >= nz2) then
            ierr = -1
            return
         end if
         if (qval <= q(k) .and. qval > q(k+1)) then
            alfa = (qval - q(k+1)) / dq(k)
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = &
                  alfa*xh(1:nvar_hydro,k) + (1-alfa)*xh(1:nvar_hydro,k+1)
            if (species > 0) &
               comp(1:species) = alfa*xa(1:species,k) + (1-alfa)*xa(1:species,k+1)
            return
         end if
         ierr = -1         
      end subroutine interp_q

      
      subroutine std_write_internals_to_file(id, num)
         integer, intent(in) :: num, id
         character (len=256) :: fname
         integer :: ierr
         ierr = 0
         write(fname, '(a, i1, a)') 'plot_data/internals', mod(abs(num), 10), '.data'
         write(*,*) 'call write_internals_to_file ' // trim(fname)
         call write_internals_to_file(id, fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in write_internals_to_file ' // trim(fname)
         end if
      end subroutine std_write_internals_to_file
      
      
      subroutine write_internals_to_file(id, filename, ierr)
         use utils_lib
         character (len=*), intent(in) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: iounit
         ierr = 0
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr == 0) then
            call write_internals(id, iounit, ierr)
            close(iounit)
         else
            write(*, *) 'failed to open internals file ' // trim(filename)
         end if
         call free_iounit(iounit)
      end subroutine write_internals_to_file
      
      
      subroutine write_internals(id, iounit, ierr)
         use chem_def
         integer, intent(in) :: iounit, id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s 
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'write_internals: get_star_ptr ierr', ierr
            return
         end if
         call write_model_info(s, iounit, ierr)
      end subroutine write_internals
      
      
      subroutine write_model_info(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s 
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, i, nz, species
         
         integer :: he4

         ierr = 0
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         he4 = s% net_iso(ihe4)
         
         write(iounit,'(a)') '            mass         initial_z       n_shells'
         write(iounit,'(2x,2f15.4,i15)') s% star_mass, s% initial_z, nz
         write(iounit,fmt='(i5)',advance='no') 1
         do i=2,88
            write(iounit,fmt='(i12,15x)',advance='no') i
         end do
         write(iounit,*)
         write(iounit,fmt='(a5,1x,99(a26,1x))',advance='no') &
               'grid', 'r', 'm', 'log_dq', &
               'log10d', 'log10T', 'log10m', 'log10r', 'L', 'r_div_rstar', &
               'log10P', 'log10Pgas', 'chiT', 'chiRho', &
               'dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas', &
               'd_dlnd_const_lnT', 'd_dlnT_const_lnd', 'xscale_lnd', 'xscale_lnT', &
               'dequP_dlnPgas', 'term2', 'term1', 'Prad_div_P'
         ! composition info
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') trim(chem_isos% name(chem_id(i)))
         end do
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') 'lg_' // trim(chem_isos% name(chem_id(i)))
         end do
         write(iounit,fmt=*)
               
         do k=1, nz
            write(iounit,'(i5,1x,99(1pd26.16,1x))',advance='no') k,  &
               s% r(k)/Rsun, s% m(k)/Msun, safe_log10_cr(s% dq(k)), &
               s% lnd(k)/ln10, s% lnT(k)/ln10, log10_cr(s% m(k)),  &
               s% lnR(k)/ln10, s% L(k)/Lsun, s% r(k)/s% r(1), &
               s% lnP(k)/ln10, s% lnPgas(k)/ln10, s% chiT(k), s% chiRho(k), &
               s% dlnRho_dlnPgas_const_T(k), s% dlnRho_dlnT_const_Pgas(k), &
               s% profile_extra(k,1), s% profile_extra(k,2), &
               s% profile_extra(k,3), s% profile_extra(k,4), &
               s% profile_extra(k,5), s% profile_extra(k,6), &
               s% profile_extra(k,7), s% profile_extra(k,8)
            do i=1, species
               write(iounit, fmt='(1pd26.16, 1x)', advance='no') s% xa(i, k)
            end do               
            do i=1, species
               write(iounit, fmt='(1pd26.16, 1x)', advance='no') safe_log10_cr(s% xa(i, k))
            end do               
            write(iounit,*)
         end do
      
      end subroutine write_model_info

      
      subroutine std_dump_model_info_for_ndiff(s, num, ierr)
         use utils_lib
         type (star_info), pointer :: s 
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         character (len=256) :: fname
         integer :: iounit
         ierr = 0
         write(fname, '(a, i1)') 'n', mod(abs(num), 10)
         write(*,*) 'dump_model_info_for_ndiff ' // trim(fname)
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(iounit, file=trim(fname), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            call free_iounit(iounit)
            return
         end if
         call dump_model_info_for_ndiff(s, iounit, ierr)
         call free_iounit(iounit)
      end subroutine std_dump_model_info_for_ndiff

      
      subroutine dump_model_info_for_ndiff(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s 
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, j, nz, species
         include 'formats'
         ierr = 0
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         write(iounit,*) 'nz', nz
         write(iounit,1) 'star_mass', s% star_mass
         write(iounit,1) 'initial_z', s% initial_z
         do k=1, nz
            do j=1, s% nvar_hydro
               write(iounit,2) trim(s% nameofvar(j)), k, s% xh(j,k)
            end do
            do j=1,species
               write(iounit,2) trim(chem_isos% name(chem_id(j))), k, clip(s% xa(j,k))
            end do
         end do
         
         contains
         
         real(dp) function clip(x)
            real(dp), intent(in) :: x
            if (.true. .or. x > 1d-30) then
               clip = x
            else
               clip = 0d0
            end if
         end function clip
         
      end subroutine dump_model_info_for_ndiff


      subroutine set_tau_base(id, tau_base, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: tau_base
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% tau_base = tau_base
      end subroutine set_tau_base
      

      real(dp) function dt_Courant(s)
         type (star_info), pointer :: s 
         integer :: k
         real(dp) :: dt
         dt_Courant = s% r(s% nz)/s% csound(s% nz)
         do k=1, s% nz-1
            dt = (s% r(k) - s% r(k+1))/s% csound(k)
            if (dt < dt_Courant) dt_Courant = dt
         end do
      end function dt_Courant


      ! largest k s.t. for all k' < k, cell k' has Cp(k')*T(k')*mstar_dot < L(k).
      subroutine set_k_CpTMdot_lt_L(s)
         type (star_info), pointer :: s 
         integer :: k, nz    
         if (s% mstar_dot <= 0d0) then
            s% k_CpTMdot_lt_L = 1
            return      
         end if    
         nz = s% nz
         do k = 2, nz
            if (s% Cp(k)*s% T(k)*s% mstar_dot >= max(1d-99,s% L(k))) then
               s% k_CpTMdot_lt_L = k-1
               return
            end if
         end do
         s% k_CpTMdot_lt_L = nz      
      end subroutine set_k_CpTMdot_lt_L
      

      subroutine set_scale_height(s)
         type (star_info), pointer :: s 
         real(dp) :: Hp, alt_Hp, alfa, beta, rho_face, P_face
         integer :: k
         include 'formats'
         do k=1,s% nz
            if (k == 1) then
               alfa = 1
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            end if
            beta = 1 - alfa
            if (alfa == 1) then
               rho_face = s% rho(k)
               P_face = s% P(k)
            else
               rho_face = alfa*s% rho(k) + beta*s% rho(k-1)
               P_face = alfa*s% P(k) + beta*s% P(k-1)
            end if
            Hp = P_face/(rho_face*s% grav(k))
            alt_Hp = sqrt(P_face / s% cgrav(k)) / rho_face
            s% scale_height(k) = min(Hp, alt_Hp)
         end do
      end subroutine set_scale_height


      real(dp) function eval_Ledd(s)
         type (star_info), pointer :: s
         real(dp) :: dtau1, dtau, dr, tau, dqsum, Ledd_sum
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         dtau1 = get_dtau1(s)
         dtau = dtau1
         tau = s% tau_factor*s% tau_base
         dqsum = s% dq(1)
         Ledd_sum = s% dq(1)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(1)
         do k = 2, s% nz
            tau = tau + dtau
            if (tau > s% surf_avg_tau) exit
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            dqsum = dqsum + s% dq(k)
            Ledd_sum = Ledd_sum + &
               s% dq(k)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(k)
         end do
         eval_Ledd = Ledd_sum/dqsum
      end function eval_Ledd
         
         
      real(dp) function total_angular_momentum(s) result(J)
         type (star_info), pointer :: s 
         include 'formats'
         if (.not. s% rotation_flag) then
            J = 0
         else
            J = dot_product(s% dm_bar(1:s% nz), s% j_rot(1:s% nz))
         end if
      end function total_angular_momentum


      real(dp) function eval_irradiation_heat(s,k)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp) :: irradiation_dq, xq, eps
         eval_irradiation_heat = 0
         if (s% irradiation_flux /= 0) then
            irradiation_dq = s% area(1)*s% column_depth_for_irradiation/s% xmstar
            xq = 1 - s% q(k)
            if (irradiation_dq > xq) then ! add irradiation heat for cell k
               eps = 0.25d0 * s% irradiation_flux / s% column_depth_for_irradiation
               if (irradiation_dq < xq + s% dq(k)) then ! only part of cell gets heated
                  eval_irradiation_heat = eps*(irradiation_dq - xq)/s% dq(k)
               else ! all of cell gets heated
                  eval_irradiation_heat = eps
               end if
            end if
         end if
      end function eval_irradiation_heat

      
      subroutine start_time(s, time0, total_all_before)
         type (star_info), pointer :: s 
         integer, intent(out) :: time0
         real(dp), intent(out) :: total_all_before
         integer :: clock_rate
         if (.not. s% doing_timing) return
         total_all_before = total_times(s)
         call system_clock(time0,clock_rate)
      end subroutine start_time

      
      subroutine update_time(s, time0, total_all_before, total)
         type (star_info), pointer :: s 
         integer, intent(in) :: time0
         real(dp), intent(in) :: total_all_before
         real(dp), intent(inout) :: total
         real(dp) :: total_all_after, other_stuff
         integer :: time1, clock_rate
         if (.not. s% doing_timing) return
         call system_clock(time1,clock_rate)
         total_all_after = total_times(s)
         other_stuff = total_all_after - total_all_before 
            ! don't double count any other stuff
         total = total + (dble(time1-time0)/clock_rate - other_stuff)
      end subroutine update_time
      
      
      real(dp) function total_times(s)
         type (star_info), pointer :: s 
         total_times = &
            s% time_evolve_step + &
            s% time_remesh + &
            s% time_adjust_mass + &
            s% time_element_diffusion + &
            s% time_struct_burn_mix + &
            s% time_solve_mix + &
            s% time_solve_burn + &
            s% time_solve_struct + &
            s% time_solve_omega_mix + &
            s% time_eos + &
            s% time_neu_kap + &
            s% time_nonburn_net + &
            s% time_mlt + &
            s% time_set_hydro_vars + &
            s% time_set_vars
      end function total_times
      
      
      subroutine dump_model_for_diff(s,io)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: io         
         integer :: k, j, i
         include 'formats'      
         write(io,1) 'dump_model_for_diff'   
         do k = 1, s% nz
            do j = 1, s% nvar_hydro
               write(io,3) s% nameofvar(j), k, j, s% xh(j,k)
            end do
            do j = s% nvar_hydro+1, s% nvar
               i = j-s% nvar_hydro
               write(io,3) 'var_' // trim(chem_isos% name(s% chem_id(i))), k, i, s% xa(i,k)
            end do
         end do      
      end subroutine dump_model_for_diff
      
      
      subroutine smooth(dc, sz)
         real(dp), dimension(:), pointer :: dc
         integer, intent(in) :: sz
         integer :: k
         k = 1
         dc(k) = (3*dc(k) + dc(k+1))/4
         do k=2,sz-1
            dc(k) = (dc(k-1) + 2*dc(k) + dc(k+1))/4
         end do
         k = sz
         dc(k) = (dc(k-1) + 3*dc(k))/4
      end subroutine smooth

      
      subroutine dq_smooth_nonconv(s, dqsm, dc, work)
         use mlt_def, only: convective_mixing
         type (star_info), pointer :: s
         real(dp), intent(in) :: dqsm
         real(dp), dimension(:) :: dc, work
      
         real(dp) :: q0, qcntr, qsurf, dq, dqsum
         integer :: nz, k, kcntr, ksurf
         
         if (dqsm <= 0) return
         nz = s% nz
         work(1:nz) = dc(1:nz)
         dq = dqsm
         do k = 1, nz
            if (s% mixing_type(k) == convective_mixing) cycle
            q0 = s% q(k)
            qcntr = q0 - dq*0.5
            qsurf = q0 + dq*0.5
            kcntr = k
            do while (kcntr < nz .and. s% q(kcntr) > qcntr)
               if (s% mixing_type(kcntr+1) == convective_mixing) exit
               kcntr = kcntr+1
            end do
            ksurf = k
            do while (ksurf > 1 .and. s% q(ksurf) < qsurf)
               if (s% mixing_type(ksurf-1) == convective_mixing) exit
               ksurf = ksurf-1
            end do
            dqsum = sum(s% dq(ksurf:kcntr))
            dc(k) = dot_product(s% dq(ksurf:kcntr),work(ksurf:kcntr))/dqsum
         end do
         
      end subroutine dq_smooth_nonconv

      
      subroutine dr_div_R_smooth_nonconv(s, dr_div_R_width, cell_dr, v, work)
         use mlt_def, only: convective_mixing
         use utils_lib, only: is_bad_num
         
         type (star_info), pointer :: s
         real(dp), intent(in) :: dr_div_R_width
         real(dp), dimension(:) :: cell_dr, v, work
      
         real(dp) :: r0, rcntr, rsurf, dr, drsum, rstar
         integer :: nz, k, j, kcntr, ksurf
         
         include 'formats'
         
         if (dr_div_R_width <= 0) return
         nz = s% nz
         rstar = s% r(1)
         work(1:nz) = v(1:nz)
         dr = dr_div_R_width*rstar
         do k = 1, nz
            if (s% mixing_type(k) == convective_mixing) cycle
            r0 = s% rmid(k)
            rcntr = r0 - dr*0.5
            rsurf = r0 + dr*0.5
            kcntr = k
            do while (kcntr < nz .and. s% rmid(kcntr) > rcntr)
               if (s% mixing_type(kcntr+1) == convective_mixing) exit
               kcntr = kcntr+1
            end do
            ksurf = k
            do while (ksurf > 1 .and. s% rmid(ksurf) < rsurf)
               if (s% mixing_type(ksurf-1) == convective_mixing) exit
               ksurf = ksurf-1
            end do
            drsum = sum(cell_dr(ksurf:kcntr))
            v(k) = dot_product(cell_dr(ksurf:kcntr),work(ksurf:kcntr))/drsum
            if (is_bad_num(v(k))) then
               v(k) = work(k)
               !write(*,2) 'v(k)', k, v(k)
               !write(*,2) 'drsum', k, drsum
               !do j=ksurf,kcntr
               !   write(*,2) 'work(j)', j, work(j)
               !   write(*,2) 'cell_dr(j)', j, cell_dr(j)
               !end do
               !stop 'debug dr_div_R_smooth_nonconv'
            end if
         end do
         
      end subroutine dr_div_R_smooth_nonconv

      
      subroutine dr_div_R_smooth(s, preserve_sign, dr_div_R_width, cell_dr, v, work)
         use utils_lib, only: is_bad_num
         
         type (star_info), pointer :: s
         logical, intent(in) :: preserve_sign
         real(dp), intent(in) :: dr_div_R_width
         real(dp), dimension(:) :: cell_dr, v, work
      
         real(dp) :: r0, rcntr, rsurf, dr, drsum, rstar, v0
         integer :: nz, k, j, kcntr, ksurf
         
         include 'formats'
         
         if (dr_div_R_width <= 0) return
         nz = s% nz
         rstar = s% r(1)
         work(1:nz) = v(1:nz)
         dr = dr_div_R_width*rstar
         do k = 1, nz
            r0 = s% rmid(k)
            v0 = work(k)
            rcntr = r0 - dr*0.5
            rsurf = r0 + dr*0.5
            kcntr = k
            do while (kcntr < nz .and. s% rmid(kcntr) > rcntr)
               if (preserve_sign .and. v(kcntr)*v0 <= 0) exit
               kcntr = kcntr+1
            end do
            ksurf = k
            do while (ksurf > 1 .and. s% rmid(ksurf) < rsurf)
               if (preserve_sign .and. v(ksurf)*v0 <= 0) exit
               ksurf = ksurf-1
            end do            
            drsum = sum(cell_dr(ksurf:kcntr))
            v(k) = dot_product(cell_dr(ksurf:kcntr),work(ksurf:kcntr))/drsum
            if (is_bad_num(v(k))) then
               v(k) = work(k)
            else if (preserve_sign .and. v(k)*work(k) <= 0) then
               v(k) = 0d0
            end if
         end do
         
      end subroutine dr_div_R_smooth

      
      subroutine fraction_scale_height_smooth( &
            s, preserve_sign, num_times, fraction, cell_dr, v, work)
         use utils_lib, only: is_bad_num
         
         type (star_info), pointer :: s
         logical, intent(in) :: preserve_sign
         integer, intent(in) :: num_times
         real(dp), intent(in) :: fraction
         real(dp), dimension(:) :: cell_dr, v, work
      
         real(dp) :: r00, rp1, r0, rcntr, rsurf, dr, drsum, &
            vsum, rstar, v0, scale_height
         integer :: nz, i, k, kk, j, kcntr, ksurf
         
         include 'formats'
         
         if (fraction <= 0 .or. num_times <= 0) return
         nz = s% nz
         rstar = s% r(1)
         do i=1,num_times
            work(1:nz) = v(1:nz)
            do k = 1, nz
               if (k < nz) then
                  scale_height = 0.5d0*(s% scale_height(k) + s% scale_height(k+1))
               else
                  scale_height = s% scale_height(k)
               end if
               dr = fraction*scale_height
               r0 = s% rmid(k)
               rcntr = r0 - dr*0.5
               rsurf = r0 + dr*0.5
               r00 = s% r(k)
               if (k < nz) then
                  rp1 = s% r(k+1)
               else
                  rp1 = s% R_center
               end if
               dr = r00 - rp1
               v0 = work(k)
               vsum = v0*dr
               drsum = dr
               do kk = k+1, nz
                  if (preserve_sign .and. v(kk)*v0 <= 0) exit
                  r00 = s% r(kk)
                  if (kk < nz) then
                     rp1 = s% r(kk+1)
                  else
                     rp1 = s% R_center
                  end if
                  dr = r00 - max(rp1,rcntr)
                  if (dr > 0) then
                     vsum = vsum + work(kk)*dr
                     drsum = drsum + dr
                  end if
                  if (rp1 <= rcntr) exit
               end do
               do kk = k-1, 1, -1
                  if (preserve_sign .and. v(kk)*v0 <= 0) exit
                  r00 = s% r(kk)
                  rp1 = s% r(kk+1)
                  dr = min(r00,rsurf) - rp1
                  if (dr > 0) then
                     vsum = vsum + work(kk)*dr
                     drsum = drsum + dr
                  end if
                  if (r00 >= rsurf) exit
               end do
               v(k) = vsum/drsum
               if (is_bad_num(v(k))) then
                  v(k) = work(k)
               else if (preserve_sign .and. v(k)*work(k) <= 0) then
                  v(k) = 0d0
               end if
            end do
         end do
         
      end subroutine fraction_scale_height_smooth
      
      
      subroutine get_XYZ(s, xa, X, Y, Z)
         use chem_def, only: ih1, ih2, ihe3, ihe4
         type (star_info), pointer :: s
         real(dp), intent(in) :: xa(:)
         real(dp), intent(out) :: X, Y, Z
         X = 0d0
         if (s% net_iso(ih1) /= 0) X = X + xa(s% net_iso(ih1))
         if (s% net_iso(ih2) /= 0) X = X + xa(s% net_iso(ih2))
         X = min(1d0, max(0d0, X))
         Y = 0d0
         if (s% net_iso(ihe3) /= 0) Y = Y + xa(s% net_iso(ihe3))
         if (s% net_iso(ihe4) /= 0) Y = Y + xa(s% net_iso(ihe4))
         Y = min(1d0, max(0d0, Y))
         Z = min(1d0, max(0d0, 1d0 - (X + Y)))
      end subroutine get_XYZ
      
      
      subroutine get_face_values(s, v_mid, v_face, ierr)
         ! simple interpolation by mass
         type (star_info), pointer :: s         
         real(dp), pointer, intent(in) :: v_mid(:)
         real(dp), pointer, intent(out) :: v_face(:)
         integer, intent(out) :: ierr         
         integer :: k
         real(dp) :: dq_sum
         ierr = 0
         v_face(1) = v_mid(1)
         do k=2, s% nz
            dq_sum = s% dq(k-1) + s% dq(k)
            v_face(k) = (v_mid(k)*s% dq(k-1) + v_mid(k-1)*s% dq(k))/dq_sum
         end do
      end subroutine get_face_values
         
         
      real(dp) function get_Ledd(s,k) result(Ledd)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp) :: kap_face
         integer :: j
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         kap_face = interp_val_to_pt(s% opacity,j,s% nz,s% dq,'get_Ledd')
         Ledd = pi4*clight*s% cgrav(j)*s% m_grav(j)/kap_face
      end function get_Ledd

      
      real(dp) function get_L_rad(s,k)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer :: j
         real(dp) :: kap_face, del_m, del_T4
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         kap_face = interp_val_to_pt(s% opacity,j,s% nz,s% dq,'get_L_rad')
         del_m = 0.5d0*(s% dm(j-1) + s% dm(j))
         del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))
         get_L_rad = -s% area(j)*s% area(j)*crad*clight/(3*kap_face)*(del_T4/del_m)
      end function get_L_rad

      
      real(dp) function get_L_rad_div_Ledd(s,k) result(L_rad_div_Ledd)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer :: j
         real(dp) :: del_m, del_T4
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         del_m = 0.5d0*(s% dm(j-1) + s% dm(j))
         del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))
         L_rad_div_Ledd = &
            -(s% area(j)*s% area(j)*crad*(del_T4/del_m)/3)/(pi4*s% cgrav(j)*s% m_grav(j))
      end function get_L_rad_div_Ledd

      
      real(dp) function get_log_concentration(s,j,k) result(log_c)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s         
         integer, intent(in) :: j, k
         ! concentration = number density / number density of electrons
         !  Ci = (Xi/Ai) / sum(Zi*Xi/Ai)   [see Thoul et al, ApJ 421:828-842, 1994]         
         integer :: i, cid, species
         real(dp) :: tmp, c         
         log_c = -1d99
         if (s% chem_id(j) == 0) return
         species = s% species
         tmp = 0d0
         do i=1,species
            cid = s% chem_id(i)
            tmp = tmp + chem_isos% Z(cid)*s% xa(i,k)/chem_isos% Z_plus_N(cid)
         end do
         cid = s% chem_id(j)
         c = (s% xa(j,k)/chem_isos% Z_plus_N(cid))/tmp
         log_c = safe_log10_cr(c)
      end function get_log_concentration

      
      real(dp) function get_phi_Joss(s,k) result(phi)
         use eos_def, only: i_lnPgas
         ! Joss, Salpeter, Ostriker, 1973. density inversion when Lrad/Ledd > phi.
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         phi = 1d0/(1d0 + (s% Pgas(k)/(4* s% Prad(k)))*s% d_eos_dlnT(i_lnPgas,k))
      end function get_phi_Joss
      
      
      logical function after_He_burn(s, he4_limit)
         use chem_def
         type (star_info), pointer :: s   
         real(dp), intent(in) :: he4_limit      
         integer :: nz, h1, he4
         real(dp) :: small = 1d-4
         after_He_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         if (h1 == 0 .or. he4 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > he4_limit) return
         after_He_burn = .true.
      end function after_He_burn
      
      
      logical function after_C_burn(s, c12_limit)
         use chem_def
         type (star_info), pointer :: s   
         real(dp), intent(in) :: c12_limit      
         integer :: nz, h1, he4, c12
         real(dp) :: small = 1d-4
         after_C_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > small .or. &
             s% xa(c12,nz) > c12_limit) return
         after_C_burn = .true.
      end function after_C_burn


      real(dp) function omega_crit(s, k)  ! result always is > 0
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp) :: Ledd, gamma_factor
         include 'formats'
         if (.not. s% rotation_flag) then
            omega_crit = 0
            return
         end if
         Ledd = 4*pi*clight*s% cgrav(k)*s% m_grav(k)/s% opacity(k)
         gamma_factor = 1d0 - min(s% L(k)/Ledd, 0.9999d0)
         omega_crit = sqrt(gamma_factor*s% cgrav(k)*s% m_grav(k)/pow3(s% r(k)))
      end function omega_crit
      
      
      subroutine set_surf_avg_rotation_info(s)
         type (star_info), pointer :: s
         real(dp) :: &
            dm, dmsum, omega_sum, omega_crit_sum, omega_div_omega_crit_sum, &
            v_rot_sum, v_crit_sum, v_div_v_crit_sum, L_div_Ledd_sum, &
            kap_face, Ledd, gamma_factor, omega_crit, omega, kap_sum, &
            j_rot_sum, j_rot, v_rot, v_crit, L_div_Ledd, dtau, tau, &
            cgrav, kap, mmid, Lmid, rmid, logT_sum, logRho_sum
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         
         if (.not. s% rotation_flag) then
            s% omega_avg_surf = 0
            s% omega_crit_avg_surf = 0
            s% w_div_w_crit_avg_surf = 0
            s% j_rot_avg_surf = 0
            s% v_rot_avg_surf = 0
            s% v_crit_avg_surf = 0
            s% v_div_v_crit_avg_surf = 0
            s% L_div_Ledd_avg_surf = 0
            s% opacity_avg_surf = 0
            s% logT_avg_surf = 0
            s% logRho_avg_surf = 0
            return
         end if

         tau = s% tau_factor*s% tau_base
         dmsum = 0d0
         L_div_Ledd_sum = 0d0
         rmid = 0d0
         
         do k = 1, s% nz - 1
            kap = s% opacity(k)
            rmid = s% rmid(k)
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            Lmid = 0.5d0*(s% L(k) + s% L(k+1))
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            dm = s% dm(k)
            dtau = dm*kap/(4*pi*rmid*rmid)
            if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
               dm = dm*(s% surf_avg_tau - tau)/dtau
               !write(*,2) 'tau limit', k, (s% surf_avg_tau - tau)/dtau
            end if
            dmsum = dmsum + dm
            Ledd = 4*pi*clight*cgrav*mmid/kap
            L_div_Ledd_sum = L_div_Ledd_sum + dm*Lmid/Ledd
            tau = tau + dtau
            if (tau >= s% surf_avg_tau) exit
         end do
         
         !write(*,2) 'L_div_Ledd_sum', s% model_number, L_div_Ledd_sum
         
         s% L_div_Ledd_avg_surf = L_div_Ledd_sum/dmsum
         if (s% generations > 2) & ! time average
            s% L_div_Ledd_avg_surf = &
               0.5d0*(s% L_div_Ledd_avg_surf + s% L_div_Ledd_avg_surf_old)
         L_div_Ledd = s% L_div_Ledd_avg_surf

         gamma_factor = 1d0 - min(L_div_Ledd, 0.9999d0)

         tau = s% tau_factor*s% tau_base
         dmsum = 0
         j_rot_sum = 0
         omega_sum = 0
         omega_crit_sum = 0
         omega_div_omega_crit_sum = 0
         v_rot_sum = 0
         v_crit_sum = 0
         v_div_v_crit_sum = 0
         kap_sum = 0
         logT_sum = 0
         logRho_sum = 0
         
         do k = 1, s% nz - 1
         
            kap = s% opacity(k)
            rmid = s% rmid(k)            
            dm = s% dm(k)
            dtau = dm*kap/(4*pi*rmid*rmid)
            
            if (tau + dtau <= s% surf_avg_tau_min) then 
               tau = tau + dtau
               cycle
            end if
            
            if (tau < s% surf_avg_tau_min) then ! only use part of this cell
               dm = dm*(tau + dtau - s% surf_avg_tau_min)/dtau
            else if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
               dm = dm*(s% surf_avg_tau - tau)/dtau
            end if
            
            dmsum = dmsum + dm
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            omega = 0.5d0*(s% omega(k) + s% omega(k+1))
            j_rot = 0.5d0*(s% j_rot(k) + s% j_rot(k+1))
            
            omega_crit = sqrt(gamma_factor*cgrav*mmid/pow3(rmid))
            v_rot = omega*rmid
            v_crit = omega_crit*rmid
            
            kap_sum = kap_sum + dm*kap
            j_rot_sum = j_rot_sum + dm*j_rot
            omega_sum = omega_sum + dm*omega
            omega_crit_sum = omega_crit_sum + dm*omega_crit
            omega_div_omega_crit_sum = omega_div_omega_crit_sum + dm*omega/omega_crit
            v_rot_sum = v_rot_sum + dm*v_rot
            v_crit_sum = v_crit_sum + dm*v_crit
            v_div_v_crit_sum = v_div_v_crit_sum + dm*v_rot/v_crit
            logT_sum = logT_sum + dm*s% lnT(k)/ln10
            logRho_sum = logRho_sum + dm*s% lnd(k)/ln10
            kap_sum = kap_sum + dm*kap
            tau = tau + dtau
            if (tau >= s% surf_avg_tau) exit

         end do

         s% logT_avg_surf = logT_sum/dmsum
         s% logRho_avg_surf = logRho_sum/dmsum
         s% opacity_avg_surf = kap_sum/dmsum
         s% j_rot_avg_surf = j_rot_sum/dmsum
         s% omega_avg_surf = omega_sum/dmsum
         s% omega_crit_avg_surf = omega_crit_sum/dmsum
         s% w_div_w_crit_avg_surf = omega_div_omega_crit_sum/dmsum
         s% v_rot_avg_surf = v_rot_sum/dmsum
         s% v_crit_avg_surf = v_crit_sum/dmsum
         s% v_div_v_crit_avg_surf = v_div_v_crit_sum/dmsum
   
      end subroutine set_surf_avg_rotation_info   
      

      subroutine median_smoothing(dd, n, ns, dmed)
         use num_lib, only: qsort
         real(dp), pointer, intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         real(dp), pointer :: dmed(:) ! (n) work array

         real(dp) :: x(2*ns+1)
         integer :: i, j, k, nmed, index(2*ns+1)

         nmed = 2*ns+1

         do i=1,n
            if ((i > 1+ns) .and. (i < n-ns)) then
               k = 1
               do j = i-ns, i+ns
                  x(k) = dd(j)
                  k = k+1
               end do
               call qsort(index,nmed,x)
               dmed(i) = x(index(ns+1))
            else
               dmed(i) = dd(i)
            end if
         end do

         do i=1,n
            if (dmed(i) /= 0) dd(i) = dmed(i)
         end do

      end subroutine median_smoothing
      
      
      subroutine weighed_smoothing(dd, n, ns, preserve_sign, ddold)
      !     based on routine written by S.-C. Yoon, 18 Sept. 2002
      !     for smoothing  any variable (dd) with size n over 2*ns+1 cells. 
         real(dp), pointer, intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         logical, intent(in) :: preserve_sign
         real(dp), pointer :: ddold(:) ! (n) work array
         
         integer :: nweight, mweight, i, j, k
         real(dp) :: weight(2*ns+1), sweight, v0
         
         include 'formats'

         do i = 1,n
           ddold(i) = dd(i)
         end do

         !--preparation for smoothing --------
         nweight = ns
         mweight = 2*nweight+1
         do i = 1,mweight
            weight(i) = 0d0
         end do
         weight(1) = 1d0
         do i = 1,mweight-1
            do j = i+1,2,-1
               weight(j) = weight(j) + weight(j-1)
            end do
         end do

         !--smoothing ------------------------
         do i=2,n-1
            sweight=0d0
            dd(i)=0d0
            v0 = ddold(i)
            do j = i, max(1,i-nweight), -1
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            do j = i+1, min(n,i+nweight)
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            if (sweight > 0) then
               sweight = 1d0/sweight
               dd(i) = dd(i)*sweight
            end if
         end do
      
      end subroutine weighed_smoothing
        
      
      ! inner radius of shell ri
      ! outer radius of shell ra
      real(dp) function eval_i_rot(s,ri,r00,ra) result(i_rot)
         type (star_info), pointer :: s         
         real(dp), intent(in) :: ri,r00,ra
         real(dp) :: rai,ra2,ri2,rm2
         
         if (s% simple_i_rot_flag) then
            i_rot = (2d0/3d0)*r00*r00
            return
         end if
         
      ! expression for evaluation without subtraction from Langer code
         rai=ra*ri
         ra2=ra*ra
         ri2=ri*ri
         rm2=ri2+rai+ra2
         i_rot=0.4D0*(ri2*ri2+rai*rm2+ra2*ra2)/rm2
         
      end function eval_i_rot
      
      
      subroutine set_i_rot(s)
         type (star_info), pointer :: s         
         integer :: k, nz
         include 'formats'
         s% i_rot(1) = eval_i_rot(s, s% rmid(1), s% r(1), s% r(1))
         do k=2,s% nz
            s% i_rot(k) = eval_i_rot(s, s% rmid(k), s% r(k), s% rmid(k-1))
         end do
         !write(*,2) 's% i_rot(1)', s% model_number, s% i_rot(1)
      end subroutine set_i_rot
      
      
      subroutine set_j_rot(s)
         type (star_info), pointer :: s         
         integer :: k
         include 'formats'
         do k=1,s% nz
   			s% j_rot(k) = s% i_rot(k)*s% omega(k)
         end do
      end subroutine set_j_rot
      
      
      subroutine set_omega(s, str)
         type (star_info), pointer :: s 
         character (len=*) :: str       
         integer :: k
         include 'formats'
         do k=1,s% nz
   			s% omega(k) = s% j_rot(k)/s% i_rot(k)
         end do
         !write(*,2) trim(str) // ' s% omega(1)', s% model_number, s% omega(1)
      end subroutine set_omega
      
      
      subroutine check_omega(s, str)
         type (star_info), pointer :: s  
         character (len=*) :: str       
         integer :: k
         logical :: okay
         include 'formats'
         okay = .true.
         do k=1,s% nz
   			if (abs(s% omega(k) - s% j_rot(k)/s% i_rot(k)) > 1d-14) then
   			   write(*,2) 'omega error', k, s% omega(k) - s% j_rot(k)/s% i_rot(k)
   			   okay = .false.
   			   exit
   			end if
         end do
         if (okay) return
         write(*,*) trim(str)
         stop 'check_omega'
      end subroutine check_omega
         
         
      subroutine use_xh_to_update_i_rot(s)
         type (star_info), pointer :: s         
         real(dp) :: r00, r003, rp1, rp13, rm13, r_in, r_out
         integer :: k, nz
         include 'formats'
         
         nz = s% nz

         if (s% simple_i_rot_flag) then
            do k=1,nz
               r00 = exp_cr(s% xh(s% i_lnR,k))
      			s% i_rot(k) = (2d0/3d0)*r00*r00
            end do
            return
         end if
         
         r00 = exp_cr(s% xh(s% i_lnR,1))
			r003 = r00*r00*r00
         rm13 = r003
         r_out = r00
         do k=1,nz
   			if (k == nz) then
   			   rp1 = s% R_center
   			else
   			   rp1 = exp_cr(s% xh(s% i_lnR,k+1))
   			end if
            rp13 = rp1*rp1*rp1
            r_in = pow_cr((r003 + rp13)/2,1d0/3d0)
   			s% i_rot(k) = eval_i_rot(s,r_in,r00,r_out)
   			rm13 = r003
   			r003 = rp13
   			r00 = rp1
            r_out = r_in
         end do
         
      end subroutine use_xh_to_update_i_rot
      
         
      subroutine use_xh_to_update_i_rot_and_j_rot(s)
         type (star_info), pointer :: s 
         call use_xh_to_update_i_rot(s)
         call set_j_rot(s)
      end subroutine use_xh_to_update_i_rot_and_j_rot


      real(dp) function eval_kh_timescale(G,M,R,L) result(kh)
         real(dp), intent(in) :: G,M,R,L
         kh = 0.75d0*G*M*M/(R*L) ! 0.75 is based on sun.  Hansen & Kawaler eqn 1.30
      end function eval_kh_timescale
      
      subroutine new_eval_kh_timescale(s,kh)
        !added by JO to calculate actual KH time-scale      
         use const_def
         type (star_info), pointer :: s
         integer                   :: k
         real(dp)                  :: Menc,U
         real(dp)                  :: kh
         real(dp)                  :: maxL
         ! if core mass greater than zero evaluate properly
         if (s% M_center .gt. 0.0) then
            kh=0.0
            Menc=(s% M_center)
            U=0.0
            do k= (s% nz),1,-1
               Menc=Menc+(s% dq(k))*((s% mstar)-(s% M_center))
               U=U+standard_cgrav*Menc/(s% r(k))*(s% dq(k))*((s% mstar)-(s% M_center))
            end do
            kh=abs(U)/(s% L(1))
         endif
         return
      end subroutine new_eval_kh_timescale
      
      subroutine set_log_abs_shear_using_xh(s, log_abs_shear)
         type (star_info), pointer :: s         
         real(dp) :: log_abs_shear(:)
         
         integer :: i_omega, i_lnR, k, nz
         real(dp), pointer :: r(:)
         real(dp) :: domega_dr, dlnomega_dlnr
         
         i_lnR = s% i_lnR
         nz = s% nz
         allocate(r(nz))
         
         do k=1,nz
            r(k) = exp_cr(s% xh(i_lnR,k))
         end do
         
         do k=2,nz-1
            domega_dr = 2*(s% omega(k-1) - s% omega(k))/(r(k-1) - r(k+1))
            dlnomega_dlnr = domega_dr*(r(k-1) + r(k+1))/(s% omega(k-1) + s% omega(k))
            log_abs_shear(k) = log10_cr(max(1d-30,min(1d30,abs(dlnomega_dlnr))))
         end do
         log_abs_shear(1) = log_abs_shear(2)
         log_abs_shear(nz) = log_abs_shear(nz-1)
         
         deallocate(r)
      
      end subroutine set_log_abs_shear_using_xh
      
      
      real(dp) function yrs_for_init_timestep(s)
         type (star_info), pointer :: s
         if (s% initial_mass <= 1) then
            yrs_for_init_timestep = 1d5
         else
            yrs_for_init_timestep = 1d5 / pow_cr(s% initial_mass,2.5d0)
         end if
      end function yrs_for_init_timestep

      
      
      subroutine set_phase_of_evolution(s) ! at start of run
         use rates_def, only: i_rate
         use chem_def, only: i_burn_c
         type (star_info), pointer :: s
         real(dp) :: power_he_burn, power_c_burn, power_neutrinos
         integer :: nz
         include 'formats'
         nz = s% nz
         if (.not. arrived_main_seq(s) .or. s% phase_of_evolution == phase_carbon_burning) return
         power_he_burn = s% power_he_burn
         power_c_burn = dot_product(s% dm(1:nz), s% eps_nuc_categories(i_rate,i_burn_c,1:nz))/Lsun
         power_neutrinos = s% power_neutrinos  
         if (s% phase_of_evolution == phase_helium_burning .and. power_c_burn > power_neutrinos) then
            !write(*, *) 'set_phase_of_evolution: phase_carbon_burning'
            s% phase_of_evolution = phase_carbon_burning
         else if (power_c_burn + power_he_burn > power_neutrinos) then
            !write(*, *) 'set_phase_of_evolution: phase_helium_burning'
            s% phase_of_evolution = phase_helium_burning
         else if (s% center_he4 < center_he_going) then
            !write(*, *) 'set_phase_of_evolution: phase_helium_burning'
            s% phase_of_evolution = phase_helium_burning
         else if (s% center_h1 < center_h_gone) then
            !write(*, *) 'set_phase_of_evolution: phase_wait_for_he'
            s% phase_of_evolution = phase_wait_for_he
         else if (s% center_h1 < center_h_going) then
            !write(*, *) 'set_phase_of_evolution: phase_mid_main_seq'
            s% phase_of_evolution = phase_mid_main_seq
         else
            !write(*, *) 'set_phase_of_evolution: phase_early_main_seq'
            s% phase_of_evolution = phase_early_main_seq
         end if
      end subroutine set_phase_of_evolution
      
      
      subroutine show_phase_of_evolution(s)
         type (star_info), pointer :: s
         include 'formats'         
         select case ( s% phase_of_evolution )
         case ( phase_starting )
            write(*, *) 'phase_starting'
         case ( phase_early_main_seq )
            write(*, *) 'phase_early_main_seq'
         case ( phase_mid_main_seq )
            write(*, *) 'phase_mid_main_seq'
         case ( phase_wait_for_he )
            write(*, *) 'phase_wait_for_he'
         case ( phase_he_igniting )
            write(*, *) 'phase_he_igniting'
         case ( phase_he_ignition_over )
            write(*, *) 'phase_he_ignition_over'
         case ( phase_carbon_burning )
            write(*, *) 'phase_carbon_burning'
         case ( phase_helium_burning )
            write(*, *) 'phase_helium_burning'
         end select
      end subroutine show_phase_of_evolution


      logical function arrived_main_seq(s)
         type (star_info), pointer :: s
         include 'formats'   
         arrived_main_seq = &
            (s% L_nuc_burn_total >= s% L_phot) .and. &
            (s% power_h_burn >= s% L_nuc_burn_total/2)
         return
         write(*,1) 's% L_nuc_burn_total', s% L_nuc_burn_total
         write(*,1) 's% L_phot', s% L_phot 
         write(*,1) 's% power_h_burn', s% L_phot 
         write(*,*) 'arrived_main_seq',  arrived_main_seq
         write(*,*)
      end function arrived_main_seq
                 
      
      subroutine save_for_d_dt(s)
         ! these values will be modified as necessary by adjust mass
         type (star_info), pointer :: s
         integer :: k, nz, i_lnR, i_lnT, i_xlnd, i_lnPgas, i_vel
         nz = s% nz
         i_lnR = s% i_lnR
         i_lnT = s% i_lnT
         i_xlnd = s% i_xlnd
         i_lnPgas = s% i_lnPgas
         i_vel = s% i_vel
         do k=1, nz
            s% lnT_for_d_dt_const_m(k) = s% xh(i_lnT, k)
            s% lnT_for_d_dt_const_q(k) = s% xh(i_lnT, k)
            s% lnR_for_d_dt_const_m(k) = s% xh(i_lnR, k)
            s% del_t_for_just_added(k) = 0d0
         end do
         if (i_xlnd /= 0) then
            do k=1, nz
               s% lnd_for_d_dt_const_m(k) = s% xh(i_xlnd, k)
               s% lnd_for_d_dt_const_q(k) = s% xh(i_xlnd, k)
            end do
         end if
         if (i_lnPgas /= 0) then
            do k=1, nz
               s% lnPgas_for_d_dt_const_m(k) = s% xh(i_lnPgas, k)
               s% lnPgas_for_d_dt_const_q(k) = s% xh(i_lnPgas, k)
            end do
         end if
         if (i_vel /= 0) then
            do k=1, nz
               s% v_for_d_dt_const_m(k) = s% xh(i_vel, k)
            end do
         end if
      end subroutine save_for_d_dt
      
      
      ! e00(i,j,k) is partial of equ(i,k) wrt var(j,k)
      subroutine e00(s,xscale,i,j,k,nvar,v)
         use utils_lib, only: is_bad_num
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, v00
         logical, parameter :: dbg = .false.
         include 'formats'
         
         if (abs(v) < 1d-250) return
         
         if (dbg) then
            if (is_bad_num(v)) then
               write(*,4) 'e00(i,j,k) ' // &
                  trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
               stop 'debug: e00'
            end if
         end if
         
         if (i == -3 .and. k == 1921) then
            write(*,4) 'e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if
         
         if (.false. .and. k == s% trace_k) then
            write(*,4) 'e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if
         
         if (j == 0) then
            write(*,*) 'called e00 with j=0 for ' // s% nameofequ(i)
            write(*,*) 's% lnPgas_flag', s% lnPgas_flag
            write(*,*) 's% i_xlnd', s% i_xlnd
            write(*,*) 's% i_lnPgas', s% i_lnPgas
            stop 'e00'
         end if
         
         if (s% hydro_matrix_type == block_tridiag_dble_matrix_type .or. &
             s% hydro_matrix_type == block_tridiag_quad_matrix_type) then
            s% dblk(i,j,k) = s% dblk(i,j,k) + v*xscale(j,k)
            return
         end if
         b = nvar*(k-1)
         q = s% idiag + b + i
         v00 = b + j
         s% jacobian(q-v00,v00) = s% jacobian(q-v00,v00) + v*xscale(j,k)
         
      end subroutine e00

      
      ! em1(i,j,k) is partial of equ(i,k) wrt var(j,k-1)
      subroutine em1(s,xscale,i,j,k,nvar,v)
         use utils_lib, only: is_bad_num
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, vm1
         logical, parameter :: dbg = .false.
         if (k == 1) return
         include 'formats'
         
         if (abs(v) < 1d-250) return

         if (dbg) then
            if (is_bad_num(v)) then
               write(*,4) 'em1(i,j,k) ' // &
                  trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
               stop 'debug: em1'
            end if
         end if
         
         if (.false. .and. s% newton_iter <= 1 .and. i == 3 .and. v /= 0) then
            write(*,4) 'em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if
         
         if (.false. .and. k == s% trace_k) then
            write(*,4) 'em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if


         if (s% hydro_matrix_type == block_tridiag_dble_matrix_type .or. &
             s% hydro_matrix_type == block_tridiag_quad_matrix_type) then
            s% lblk(i,j,k) = s% lblk(i,j,k) + v*xscale(j,k-1)
            return
         end if
         b = nvar*(k-1)
         q = s% idiag + b + i
         vm1 = b + j - nvar
         s% jacobian(q-vm1,vm1) = s% jacobian(q-vm1,vm1) + v*xscale(j,k-1)
         
      end subroutine em1
      
      
      ! ep1(i,j,k) is partial of equ(i,k) wrt var(j,k+1)
      subroutine ep1(s,xscale,i,j,k,nvar,v)
         use utils_lib, only: is_bad_num
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, vp1
         logical, parameter :: dbg = .false.
         include 'formats'
         if (k == s% nz) return
         
         if (abs(v) < 1d-250) return
         
         if (dbg) then
            if (is_bad_num(v)) then
               write(*,4) 'ep1(i,j,k) ' // &
                  trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
               stop 'debug: ep1'
            end if
         end if
         
         if (.false. .and. s% newton_iter <= 1 .and. i == 3 .and. v /= 0) then
            write(*,4) 'ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if
         
         if (.false. .and. k == s% trace_k) then
            write(*,4) 'ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         end if

         if (s% hydro_matrix_type == block_tridiag_dble_matrix_type .or. &
             s% hydro_matrix_type == block_tridiag_quad_matrix_type) then
            s% ublk(i,j,k) = s% ublk(i,j,k) + v*xscale(j,k+1)
            return
         end if
         b = nvar*(k-1)
         q = s% idiag + b + i
         vp1 = b + j + nvar
         s% jacobian(q-vp1,vp1) = s% jacobian(q-vp1,vp1) + v*xscale(j,k+1)
         
      end subroutine ep1

      
      real(dp) function current_min_xa_hard_limit(s) result(min_xa_hard_limit)
         type (star_info), pointer :: s      
         real(dp) :: logTc, alfa
         logTc = s% lnT(s% nz)/ln10
         if (logTc <= s% logT_max_for_xa_hard_limit) then
            min_xa_hard_limit = s% min_xa_hard_limit
         else if (logTc >= s% logT_min_for_xa_hard_limit) then
            min_xa_hard_limit = s% min_xa_hard_limit_for_highT
         else
            alfa = (logTc - s% logT_max_for_xa_hard_limit) / &
                   (s% logT_min_for_xa_hard_limit - s% logT_max_for_xa_hard_limit)
            min_xa_hard_limit = &
               alfa*s% min_xa_hard_limit_for_highT + (1d0 - alfa)*s% min_xa_hard_limit
         end if
      end function current_min_xa_hard_limit



      end module star_utils
