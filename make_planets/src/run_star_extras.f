! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c EDITED by JO Dec '12 to include planetary evaaporation from the 
!c Owen & Jackson 2012 models
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def

      implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     Extra run_star_extra module variables
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! dust opacity variables
      double precision, dimension(:), allocatable, save ::  kn,an,bn
      ! AU to cm
      double precision, parameter                       ::  au_to_cm = 1.49597871d13

      logical                                       ::  read_tables
      logical                                       ::  read_star
      real(dp), allocatable, dimension(:), save     ::  Rp_evap,Mp_evap
      real(dp), allocatable, dimension(:,:,:), save ::  Mdot_evap
      real(dp), allocatable, dimension(:,:), save   ::  Mdot_table
      real(dp), allocatable, dimension(:), save     ::  si_age, si_rad, si_T
      integer      :: NR,NM
      real(dp), dimension(6)     :: tableflux=(/ 3.18831e3, 3.5368e3, &
            8.8419e3, 3.5368e4, 7.9577e4, 3.5368e5 /)
      real(dp)                   :: Tflux=3.53677651e6 ! flux in evaporation table


      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
        ierr = 0
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         s% other_energy => core_envelope_heating
         s% other_atm    => planetary_atm
         s% other_kap_get_Type1 => planetary_kap1
         s% other_wind   => planet_evap         
      end subroutine extras_controls

      subroutine core_envelope_heating(id, ierr)        

        !written by JO to include extra_heating terms necessary for planetary
        !evolution

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

	 ! MY VARIABLES
         integer :: k
         double precision :: S_set, t_cool, S_env, T_env
         
         !core heating
         real(dp)              ::  dTc_dt !time derivative of core temperature
         real(dp)              ::  radio_fact=1 ! adjust the radioactive heating from nominal
         real(dp)              ::  heatc_fact=1 ! adjust the heat capacity factor from nominal
         real(dp), parameter   ::  edot_rad=3.3481e-8 !(erg/s/g) from valencia et al. (2010)
         real(dp), parameter   ::  cv_core=1d7 !(erg/g/K) from valencia et al. (2010)
         real(dp)              ::  Rstar_int
         real(dp)              ::  Teff_int
         real(dp)              ::  Teq_new

         integer               ::  star_length=79
         integer               ::  i,i_gt

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

	 ! INITALIZE
         s% extra_heat = 0.d0
         s% d_extra_heat_dlnd(:) = 0d0
         s% d_extra_heat_dlnT(:) = 0d0

         S_set = s% x_ctrl(6) * kerg / amu
         t_cool = s% x_ctrl(7) * secyer

         if ((s% x_logical_ctrl(3))) then
            do k = 1, s%nz
               S_env = exp( s% lnS(k) )
               T_env = s% T(k) 
               s% extra_heat(k) = s% extra_heat(k) + T_env * ( S_set - &
                    S_env )/ t_cool
            end do
         endif

        !core's luminosity
        if ((s% x_logical_ctrl(1)).or.(s% x_logical_ctrl(2))) then
           !pre-asign
           s% L_center = 0d0
        endif
        
        if ((s% x_logical_ctrl(1))) then
           ! add core heating to centre
           ! from radio-active decay
           s% L_center = s% L_center + radio_fact*edot_rad*(s% M_center)           
        endif
        if ((s% x_logical_ctrl(2))) then
           ! add core heating to centre
           ! from heat capacity
           dTc_dt=(s% dlnT_dt(s% nz))*(s% T(s% nz))/secyer
           if (dTc_dt .lt. 0) then
              s% L_center = s% L_center - heatc_fact*cv_core*(s% M_center)*dTc_dt
           endif
        endif
        if ((s% x_logical_ctrl(6))) then
           !add extra core heating to puff up envelope
           !in terms L/M in cgs units
           s% L_center = s% L_center +(s% x_ctrl(22))*((s% star_mass)*msol)
        endif

        if ((s% x_logical_ctrl(8))) then
          ! need to read in stellar evolution file
          if (read_star .eqv. .false.) then

            allocate(si_age(star_length))
            allocate(si_rad(star_length))
            allocate(si_T(star_length))

            open(78,file='star_info1Msun.dat')
            do i=1,star_length
              read(78,*) si_age(i), si_rad(i), si_T(i)
            end do
            close(78)
            read_star=.true.
          endif

          !now need to compute new irradiation flux
          if (s% star_age .gt. si_age(1)) then 
            !need to adjust irradiation flux
            do i=1,star_length
              if (si_age(i) .gt. s% star_age) then
                i_gt=i
                exit
              endif
            end do 
            Rstar_int=si_rad(i_gt-1)+(si_rad(i_gt)-si_rad(i_gt-1))* & 
            & (s% star_age-si_age(i_gt-1))/(si_age(i_gt)-si_age(i_gt-1))

            Teff_int=si_T(i_gt-1)+(si_T(i_gt)-si_T(i_gt-1))* & 
            & (s% star_age-si_age(i_gt-1))/(si_age(i_gt)-si_age(i_gt-1))

            ! now calculate new flux
            Teq_new=Teff_int*(Rstar_int*rsol/(2.0*(s% x_ctrl(12))*1.5e13))**0.5
            s% irradiation_flux=4.0*5.6704e-5*Teq_new**4.0

          endif
        endif


      end subroutine core_envelope_heating

      subroutine planet_evap(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         
         implicit none

         integer, intent(in)  :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         real(dp)             :: L10R,L10M,L10Mdot_u, L10Mdot_b,L10Mdot
         real(dp)             :: Tstat, PLfoff, Lx_Lbol, Lbol
         real(dp)             :: adjust,smooth,mdot,mdot_EUV
         real(dp)             :: HEflux,l10HEflux,l10F1,l10F2
         real(dp)             :: dif1,dif2
         integer              :: n_up, n_below,n
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         type (star_info), pointer     :: s

         integer               :: which_evap

         which_evap=2
         
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! USED to calculate photoevaporation rates of planets


          if (read_tables .eqv. .false.) then
            !need to read evpaoration tables
            call read_evap
            read_tables=.true.
          endif

          !use seperation and Teq to calculate Lbol
          !Lbol=16.0*pi*((s% x_ctrl(12))*(au))**2.0*boltz_sigma*(s% x_ctrl(19))**(4.0)

          Tstat=(s% x_ctrl(24))
          Lx_Lbol=(s% x_ctrl(25))
          PLfoff=(s% x_ctrl(26))

          !compute X-ray flux
          HEflux=Lx_Lbol*(s% irradiation_flux)

         if ((s% x_logical_ctrl(7))) then 

          
          if (s% star_age .gt.  Tstat) HEflux=HEflux*((s% star_age)/Tstat)**(-PLfoff)

          if (which_evap == 1) then

           n_up=6
           n_below=6
           do n=1,6
             if (HEflux .lt. tableflux(n)) then
               n_up=n
               n_below=n-1
               exit
             endif
           end do
           if (n_below .eq. 0) n_below=1

           L10R=log10(10.**(s% log_surface_radius)*Rsol)
           L10M=log10((s% star_mass)*Msol)

           if (n_below .ne. n_up) then
            !now do interpolation to find mdot
            !call interpolate_evap(L10R,L10M,L10Mdot_u,n_up)
            !call interpolate_evap(L10R,L10M,L10Mdot_b,n_below)
            !now interpolate in flux
            if (L10Mdot_b*L10Mdot_u .gt. 0) then !both mdots present
              !linear interpolation
              L10Mdot=(L10Mdot_u-L10Mdot_b)/(tableflux(n_up)-tableflux(n_below))*(HEflux - tableflux(n_below)) &
              +L10Mdot_b
              mdot=10.0**L10Mdot
            else
              if (L10Mdot_u .eq. 0) then 
                mdot=0.0
              else
                !smoothly interpolate to zero using ~ Flux * exp cut_off
                smooth=(Heflux/tableflux(n_up)) !linear scaling due to flux
                smooth=smooth*exp(-0.1*(tableflux(n_up)-Heflux)/(Heflux - tableflux(n_below)))
                mdot=(10.0**L10Mdot_u)*smooth
              endif
            endif
           else
            !call interpolate_evap(L10R,L10M,L10Mdot,n_up)
            mdot=10.0**L10Mdot*(HEflux/tableflux(n_up))
           endif

           ! check against EUV rate
           call get_EUV_rate(HEflux,mdot_EUV,(s% star_mass)*Msol,10.**(s% log_surface_radius)*Rsol)
           if (mdot_EUV .gt. mdot) mdot=mdot_EUV

           w=mdot

           !correct units to Msun/yr
         
           w=w/Msol*secyer !Msun/year

          elseif (which_evap==2) then
          	L10R=log10(10.**(s% log_surface_radius)*Rsol)
            L10M=log10((s% star_mass)*Msol)
            call interpolate_evap(L10R,L10M,L10Mdot)

            mdot=10.0**L10Mdot*(HEflux/Tflux)

            !print*, L10R,L10M,mdot, HEflux

            w=mdot/Msol*secyer !in Msun/year
			call get_EUV_rate(HEflux,mdot_EUV,(s% star_mass)*Msol,10.**(s% log_surface_radius)*Rsol)
			if (mdot_EUV .gt. mdot) mdot=mdot_EUV


          endif
         else
             w=0d0
         endif

       end subroutine planet_evap

       subroutine get_EUV_rate(HEflux,mdot,Mp,Rp)

        implicit none

        real(dp), intent(in) :: HEflux,Mp,Rp !cgs units
        real(dp), intent(out):: mdot !cgs units
        real(dp), parameter  :: eff=0.15 ! efficiency

        ! calculates a lower bound to the mdot due to pure EUV mass-loss

        mdot=eff*pi*HEflux*Rp**3.0/(4.0*standard_cgrav*Mp) ! low EUV rate Equation 1 of Owen & Jackson (2012)

        return
       end subroutine get_EUV_rate

       subroutine read_evap

       !reads in photoevaporation tables

       character*32          ::      evap_table
       integer               ::      table_length
       integer               ::      i,j,k,n
       real(dp), allocatable, dimension(:) :: evap_M, evap_R, evap_Mdot

       ! Hard code this for now
       table_length=5600
       NR=140
       NM=40

       allocate(Rp_evap(NR))
       allocate(Mp_evap(NM))
       allocate(Mdot_evap(NM,NR,6))
       allocate(Mdot_table(NM,NR))

       allocate(evap_M(table_length))
       allocate(evap_R(table_length))
       allocate(evap_Mdot(table_length))

       evap_M=0.0
       evap_R=0.0
       evap_Mdot=0.0

       do n=1,6
        write(evap_table,'(a,I1.1,a)') "Owen_Rates/Mdot",n,".out"      
        !read in the file
        open(unit=73,file=evap_table)

        do i=1,table_length
           read(73,*) evap_R(i), evap_M(i), evap_Mdot(i)
        end do

        close(73)
        !put in store tables
        k=1
        do i=1,NR
          do j=1,NM
             Rp_evap(i)=evap_R(k)
             Mp_evap(j)=evap_M(k)
             Mdot_evap(j,i,n)=evap_Mdot(k) !reversal to be consistent with MATLAB style indexing
             k=k+1
          end do
        end do
  
       end do

       ! read in new table
       evap_table='Owen_Rates/Mdot_out_Mar15.out'
       open(unit=74, file=evap_table)
       do i=1,table_length
       	read(74,*) evap_R(i), evap_M(i), evap_Mdot(i)
       end do
       close(74)

       !put in store tables
       k=1
       do i=1,NR
          do j=1,NM
             Rp_evap(i)=evap_R(k)
             Mp_evap(j)=evap_M(k)
             Mdot_table(j,i)=evap_Mdot(k) !reversal to be consistent with MATLAB style indexing
             k=k+1
          end do
       end do

     
       end subroutine read_evap

      subroutine interpolate_evap(LR,LM,LMD)
      
      !performs bi-linear interpolation of the evaporation tables
      
      real(dp), intent(in)     :: LR,LM !logarithmic values of Radius and Mass [cgs]
      real(dp), intent(out)    :: LMD   !lograithmic values of the mass-loss rate [cgs]

      integer                  :: NRad,NMass !lower index of evaporation table
      integer                  :: i,j,n
      
      !variables for linear interpolation
      real(dp)                 :: int1, int2, int3

      !variables for nearest neighbour interpolation
      real(dp)                 :: store_dist
      real(dp), dimension(4,2) :: distance !distance to point, mdot at that point

      logical                  :: use_nearest

      use_nearest=.false.

      Nrad=0
      NMass=0
      do i=1,NR
         if (LR .lt. Rp_evap(i)) then
            NRad=i-1
            exit
         endif
      end do

      do j=1,NM
         if (LM .lt. Mp_evap(j)) then
            NMass=j-1
            exit
         endif
      end do

!       if (Nrad==0) then
!          print*, "Planet Radius not on evaporation table!"
!          print*, 10**LR
!          return
!       endif
!       if (Nmass==0) then
!          print*, "Planet Mass not on evaporation table!"
!          print*, 10**LM
!          return
!       endif
      
      !now check not near boundary
      
      if (Mdot_table(NMass  ,Nrad  ) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass+1,Nrad  ) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass  ,Nrad+1) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass+1,Nrad+1) .lt. 1) use_nearest=.true.

      !Now perform interpolation
      if (use_nearest) then
         !use closetest non-zero value
         !Also throw warning to terminal
         distance(1,1)=((LR-Rp_evap(Nrad  ))**2+(LM-Mp_evap(NMass  ))**2)
         distance(2,1)=((LR-Rp_evap(Nrad+1))**2+(LM-Mp_evap(NMass  ))**2)
         distance(3,1)=((LR-Rp_evap(Nrad  ))**2+(LM-Mp_evap(NMass+1))**2)
         distance(4,1)=((LR-Rp_evap(Nrad+1))**2+(LM-Mp_evap(NMass+1))**2)

         distance(1,2)=Mdot_table(NMass  ,NRad  )
         distance(2,2)=Mdot_table(NMass+1,NRad  )
         distance(3,2)=Mdot_table(NMass  ,NRad+1)
         distance(4,2)=Mdot_table(NMass+1,NRad+1)
         
         store_dist=1d100
         LMD=0.
         do i=1,4
            if (distance(i,1).lt.store_dist) then
               if (distance(i,2) .gt. 1) then
                  LMD=distance(i,2)
                  store_dist=distance(i,1)
               endif
            endif
         end do   
!          if (LMD .gt. 0) then
!             print*, "Warning Using Nearest Neighbour Evaporation Interpolation"
!          endif
      else
         !linear interpolataion

         !Mass first then Radius
         int1=Mdot_table(Nmass,NRad)+(Mdot_table(Nmass+1,NRad)-Mdot_table(Nmass,Nrad)) &
         &  *((LM-Mp_evap(Nmass))/(Mp_evap(Nmass+1)-Mp_evap(Nmass)))
      
         int2=Mdot_table(Nmass,NRad+1)+(Mdot_table(Nmass+1,NRad+1)-Mdot_table(Nmass,Nrad+1)) &
         &  *((LM-Mp_evap(Nmass))/(Mp_evap(Nmass+1)-Mp_evap(Nmass)))
      
         int3=int1+(int2-int1)/(Rp_evap(NRad+1)-Rp_evap(Nrad))*(LR-Rp_evap(Nrad))
         LMD=int3
      endif

      end subroutine interpolate_evap




      subroutine planetary_atm( &
          id, M, R, L, X, Z, kap, Teff, &
          lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
          lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
          which_atm_option, switch_to_grey_as_backup, ierr)

          !calculates various planetary atmosphere's 
          !implemented by JO August 2013
          use star_lib, only: star_ptr
      
          implicit none

          type (star_info), pointer :: s

          integer, intent(in) :: id ! star id if available; 0 otherwise
          !     cgs units
          double precision, intent(in) :: M ! enclosed mass at atmosphere/interior boundary
          double precision, intent(in) :: R ! radius at atmosphere/interior boundary
          double precision, intent(in) :: L ! luminosity at atmosphere/interior boundary
          double precision, intent(in) :: X ! hydrogen mass fraction
          double precision, intent(in) :: Z ! metallicity
          double precision, intent(in) :: kap
          ! opacity above photosphere (average by mass)
          ! used to estimate gamma for tau boundary

          double precision, intent(out) :: Teff ! temperature at photosphere (essentially meaningless here)
          double precision, intent(out) :: lnT ! natural log of temperature at base of atmosphere
          double precision, intent(out) :: lnP ! natural log of pressure at base of atmosphere (Pgas + Prad)

          ! partial derivatives of lnT and lnP
          double precision, intent(out) :: dlnT_dL, dlnT_dlnR, dlnT_dlnM,dlnT_dlnkap
          double precision, intent(out) :: dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap

          integer, intent(in) :: which_atm_option
          ! atm_simple_photosphere or one of the options for tables
          ! T(tau) integration is not supported in this routine -- call atm_int_T_tau instead
          logical, intent(in) :: switch_to_grey_as_backup
          ! if you select a table option, but the args are out of the range of the tables,
          ! then this flag determines whether you get an error or the code automatically
          ! switches to option = atm_simple_photosphere as a backup.

          integer, intent(out) :: ierr ! == 0 means AOK



          !local variables
          !opacity variables, see Rogers & Seager 2010 
          double precision, parameter   :: C_atm= 4.786300923226380e-08
          double precision, parameter   :: a_atm=0.68
          double precision, parameter   :: b_atm=0.45
          double precision, parameter   :: b_4=0.1125 !b_kap/4

          double precision              :: gamma, Kappa_v, tau_IR, bound_fact
          double precision              :: Tint4, Teq4

          double precision              :: Patm, Tatm ! cgs valuse of P & T at base 

          double precision              :: Tmin4,Tmax4,Tmiddle4 !max and min temperature to fourth power
          double precision              :: gmin, gmax, gmiddle ! gamma at min and max values
          double precision              :: fmin,fmax,fmiddle !function at max or min

          double precision              :: gatm ! g=GM/R^2
          
          double precision              :: En(0:2) !exponetial integral 0,1,2
          double precision              :: big_bracket

          integer                       :: Max_its,its
          double precision              :: conv_tol ! relative convergance tolerance
          logical                       :: converged

          ! access x_ctrl(1) to see which type of atmosphere is required
          !option  1 - fixed P & T boundary (aka disc atmosphere)
          !option  2 - Guillot (2010) T(tau) relation Eqn (49) with RS10 Opacity fit

          ierr = 0
          call get_star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (s% x_integer_ctrl(1) .eq. 1) then
          !fixed P & T boundary
             !Pressure in Xctrl(2), Temperature in Xctrl(3)
             !read in values from pointers
             Patm= s% x_ctrl(2)
             Tatm= s% x_ctrl(3)

             s% tau_factor = (s% x_ctrl(1)) / (s% tau_base)

             ! assign outputs
             Teff=(L/(4d0*pi*boltz_sigma*R**2))**0.25 ! effective temperature meaninless here
             Teff=maxval((/ Tatm, Teff/))
             lnP=log(Patm)
             lnT=log(Tatm)
             ! all derivatives zero in this implementation
             dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
             dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
             


          elseif (s% x_integer_ctrl(1) .eq. 2) then
             !Guillot (2010) T(tau) relation
             !solves for a given Temperature and Pressure at the atmosphere
             !at fixed optical depth to the out going irradiation
             
             !Equation 49 of Guillot (2010) A&A 520,A27 specifies T(tau)

             !solves this using bi-section method where max and min are controlled using x_ctrl inputs

             !set exit condition to false initially
             converged=.false.

             !calculate gravity
             gatm=standard_cgrav*M/R**2d0 


             !import opacity to optical
             kappa_v=(s% x_ctrl(15))
             !find optical depth at boundary (estimate gamma from kap)
             bound_fact=(s% x_ctrl(16))
             tau_IR=bound_fact*(2d0/3d0)
             if (tau_IR .lt. 2d0/3d0) tau_IR=2d0/3d0
             ! calculate new tau factor 
             s% tau_factor= tau_IR / (s% tau_base)

             !import max iterations
             max_its=(s% x_integer_ctrl(2))
             its=0 ! number of iterations initially zero
             !import convergance tolerance
             conv_tol=(s% x_ctrl(20))

             !evaluate Tint4 & Teq4
             Teq4=(s% x_ctrl(19))**4d0
             Tint4=(L/(4*pi*boltz_sigma*R**2))
             

             !set bi-section max and min
             Tmin4=(s% x_ctrl(17))**4d0
             Tmax4=(s% x_ctrl(18))**4d0


             do while (converged .eqv. .false.)
                !calculate gamma for these values
                gmin=kappa_v/(C_atm*(gatm*tau_IR/(C_atm*Tmin4**b_4))**(a_atm/(1+a_atm))*Tmin4**b_4)
                gmax=kappa_v/(C_atm*(gatm*tau_IR/(C_atm*Tmin4**b_4))**(a_atm/(1+a_atm))*Tmin4**b_4)
                !now evaluate guillot function at max and min
                call guillot_eval(tau_IR,gmin,Teq4,Tint4,Tmin4,fmin)
                call guillot_eval(tau_IR,gmax,Teq4,Tint4,Tmax4,fmax)

                if ((fmin*fmax) > 0.) then
                   !failed in bi-section method
                   ierr=-1
                   write(*,*) "failed in guillot atm at its",its
                   write(*,*) "failed because Tmin and Tmax do not straddle roots"
                   write(*,*) "fmin=", fmin, "fmax=", fmax, "fmin*fmax", fmin*fmax
                   write(*,*) "Tmin=", Tmin4**0.25, "Tmax=", Tmax4**0.25
                   write(*,*) "gamma_min=", gmin, "gamma_max=", gmax
                   write(*,*) "FAILED in other_atm"                   
                   write(*,*) "Star Info Below"
                   write(*,*) "Mass", M
                   write(*,*) "Radius", R
                   write(*,*) "Luminosity", L
                   write(*,*) "tau_IR", tau_IR
                   
                   Teff = 0; lnT = 0; lnP = 0
                   dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
                   dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
                   return
                endif
                !now find middle Tmean in linear T space
                Tmiddle4=(Tmin4**0.25+(Tmax4**0.25-Tmin4**0.25)/2)**4d0
                gmiddle=kappa_v/(C_atm*(gatm*tau_IR/(C_atm*Tmiddle4**b_4))**(a_atm/(1+a_atm))*Tmiddle4**b_4)
                call guillot_eval(tau_IR,gmiddle,Teq4,Tint4,Tmiddle4,fmiddle)

                if (fmax*fmiddle < 0. ) then
                   ! root sits between gmax and gmiddle - relabel middle to min
                   Tmin4=Tmiddle4                   
                elseif (fmin*fmiddle < 0. ) then
                   ! root sits between gmin and gmiddle - relabel middle to max
                   Tmax4=Tmiddle4
                else
                   !failed somewhere
                   ierr=-1
                   write(*,*) "Failed to find Tmiddle in guillot atm at its",its
                   write(*,*) "fmin,fmiddle,fmax below"
                   write(*,*) fmin,fmiddle,fmax
                   write(*,*) "Tmin,Tmiddle,Tmax below"
                   write(*,*) Tmin4**0.25,Tmiddle4**0.25,Tmax4**0.25
                   write(*,*) "gamma_min,gamma_middle,gamma_max below"
                   write(*,*) gmin,gmiddle,gmax
                   write(*,*) "FAILED in other_atm"
                   write(*,*) "Star Info Below"
                   write(*,*) "Mass", M
                   write(*,*) "Radius", R
                   write(*,*) "Luminosity", L
                   write(*,*) "tau_IR", tau_IR
                   return                
                endif

                !convergence testing
                if ((Tmax4**0.25-Tmin4**0.25)/(Tmax4**0.25) .lt. conv_tol) then
                   ! converged
                   converged=.true.
                   exit
                endif

                its=its+1
                !max_its check
                if (its .gt. max_its) then
                   converged=.true.
                   write(*,*) "Atmosphere not converged: failed"
                   write(*,*) "failed to converge in",its,"iterations"
                   write(*,*) fmin,fmiddle,fmax
                   write(*,*) "Tmin,Tmiddle,Tmax below"
                   write(*,*) Tmin4**0.25,Tmiddle4**0.25,Tmax4**0.25
                   write(*,*) "gamma_min,gamma_middle,gamma_max below"
                   write(*,*) gmin,gmiddle,gmax
                   write(*,*) "FAILED in other_at"
                   write(*,*) "Star Info Below"
                   write(*,*) "Mass", M
                   write(*,*) "Radius", R
                   write(*,*) "Luminosity", L
                   write(*,*) "tau_IR", tau_IR
                   ierr=-1
                   return
                endif
             end do

             !now have converged temperature
             !use Tmin as actual value 
             
             Tatm=Tmin4**0.25
             Patm=(gatm*tau_IR/(C_atm*Tatm**b_atm))**(1/(1+a_atm))

             !just output Teff as Tatm
             Teff=Tatm
             lnT=log(Tatm)
             lnP=log(Patm)
          
             !partials of T
             dlnT_dL=((2d0/3d0+tau_IR)*0.75*Tint4/(4d0*Tatm**4d0))/L !dlnT_dlnL / L
             dlnT_dlnR=0. !zero for fixed optical depth at constant pressure
             dlnT_dlnM=0. !zero for fixed optical depth at constant pressure
             !evaluate exonetial integral
             call enxa ( 2, gmin*tau_IR, En )
             big_bracket=-2d0/(3d0*gmin**2)-gmin*tau_IR/3d0*exp(-gmin*tau_IR)+2d0/(3d0*gmin**2d0)*exp(-gmin*tau_IR) &
                  +2d0/3d0*exp(-gmin*tau_IR)+2d0/3d0*(1-(tau_IR**2d0/2d0))*En(2)  & 
                  -2d0*gmin*tau_IR/3d0*(1-(tau_IR**2d0)/2d0)*En(1)

             dlnT_dlnkap=-big_bracket*gmin*(3d0*Teq4/(16d0*Tatm**4d0))

             !partials of P
             ! just come from P=(gatm*tau/(kappa))
             dlnP_dL=0. !change in pressure with luminosity at fixed T is zero
             dlnP_dlnR=-2d0
             dlnP_dlnM=1d0
             dlnP_dlnkap=-1d0

          elseif ((s% x_integer_ctrl(1)) .eq. 3) then
             !Guillot (2010) T(tau) relation
             !solves for a given Temperature and Pressure at the atmosphere
             !at fixed optical depth to the out going radiation
             !unlike option 2 this assumes that the gamma factor is specified in the inlists
             !calculate gravity
             gatm=standard_cgrav*M/R**2d0 

             bound_fact=(s% x_ctrl(16))
             tau_IR=bound_fact*(2d0/3d0)
             tau_IR=2d0
             if (tau_IR .lt. 2d0/3d0) tau_IR=2d0/3d0
             ! calculate new tau factor 
             s% tau_factor= tau_IR / (s% tau_base)

             !evaluate Tint4 & Teq4
             Teq4=(s% x_ctrl(19))**4d0
             Tint4=(L/(4*pi*boltz_sigma*R**2))
             gmiddle=(s% x_ctrl(15)) ! get gamma from x_ctrl
             
             !get the temperature
             call guillot_eval(tau_IR,gmiddle,Teq4,Tint4,0d0,Tmiddle4)
             !assign atmosphere pressure and temperature
             Tatm=(-Tmiddle4)**0.25 ! Tmiddle4 is given minus sign by guillot_eval
             Patm=(gatm*tau_IR/(C_atm*Tatm**b_atm))**(1/(1+a_atm))
             !outputs
             !just output Teff as Tatm
             Teff=Tatm
             lnT=log(Tatm)
             lnP=log(Patm)

             !assing derivatives
             !partials of T
             dlnT_dL=((2d0/3d0+tau_IR)*0.75*Tint4/(4d0*Tatm**4d0))/L !dlnT_dlnL / L
             dlnT_dlnR=0. !zero for fixed optical depth at constant pressure
             dlnT_dlnM=0. !zero for fixed optical depth at constant pressure

             !evaluate exonetial integral
             call enxa ( 2, gmiddle*tau_IR, En )
             big_bracket=-2d0/(3d0*gmiddle**2)-gmiddle*tau_IR/3d0*exp(-gmiddle*tau_IR)+2d0/(3d0*gmiddle**2d0)*exp(-gmiddle*tau_IR) &
                  +2d0/3d0*exp(-gmiddle*tau_IR)+2d0/3d0*(1-(tau_IR**2d0/2d0))*En(2)  & 
                  -2d0*gmiddle*tau_IR/3d0*(1-(tau_IR**2d0)/2d0)*En(1)

             dlnT_dlnkap=-big_bracket*gmiddle*(3d0*Teq4/(16d0*Tatm**4d0))

             !partials of P
             ! just come from P=(gatm*tau/(kappa))
             dlnP_dL=0. !change in pressure with luminosity at fixed T is zero
             dlnP_dlnR=-2d0
             dlnP_dlnM=1d0
             dlnP_dlnkap=-1d0

          else


            

             write(*,*) "no atmosphere option selected"
             write(*,*) "option set by x_integer_ctrl(1), which is currently",s% x_integer_ctrl(1)
             write(*,*) "FAILED IN other_atm"
             ierr=-1
             Teff = 0; lnT = 0; lnP = 0
             dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
             dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
         
          endif

          return
      
        end subroutine planetary_atm

        subroutine guillot_eval(tau,gamma,Teq4,Tint4,Tmean4,func)

          implicit none

          !evaluates Equation 49 of Guillot et al. (2010) for input values
          
          double precision, intent(in)    ::  tau !optical depth in IR
          double precision, intent(in)    ::  gamma ! ratio of opacities Kappa_v/Kappa_IR
          double precision, intent(in)    ::  Teq4 !equilibrium temperature to the fourth power
          double precision, intent(in)    ::  Tint4 !Internal temperature to the fourth power
          double precision, intent(in)    ::  Tmean4 !mean atmosphere temperature to the fourth power
          double precision, intent(out)   ::  func  ! value of function

          !internals
          double precision                ::  En(0:2) !exponetial integral 0,1,2
          double precision                ::  big_bracket ! curly brackets of EQN 49
          double precision                ::  gt !gamma*tau
          
          !evaluate gamma*tau
          gt=gamma*tau
          !evaluate exonetial integral
          call enxa ( 2, gt, En )

          big_bracket=(2d0/3d0)+(2d0/(3d0*gamma))*(1d0+(gt/2d0-1d0)*exp(-gt))+2d0*gamma/3d0*(1-(tau**2d0)/3d0)*En(2)

          func=Tmean4-0.75*Tint4*(2d0/3d0+tau)-0.75*Teq4*big_bracket
          
        end subroutine guillot_eval


        subroutine planetary_kap1( &
             id, k, handle, zbar, X, Zbase, log10_rho, log10_T,  &
             species, chem_id, net_iso, xa, &
             lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
             kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         
          use chem_def, only: num_chem_isos
 
          ! INPUT
          integer, intent(in) :: id ! star id if available; 0 otherwise
          integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
          integer, intent(in) :: handle ! from alloc_kap_handle
          real(dp), intent(in) :: zbar ! average ion charge
          real(dp), intent(in) :: X ! the hydrogen mass fraction
          real(dp), intent(in) :: Zbase ! the metallicity
          real(dp), intent(in) :: log10_rho ! the density
          real(dp), intent(in) :: log10_T ! the temperature
          integer, intent(in) :: species
          integer, pointer :: chem_id(:) ! maps species to chem id
          ! index from 1 to species
          ! value is between 1 and num_chem_isos         
          integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
          real(dp), intent(in) :: xa(:) ! mass fractions
          double precision, intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
          real(dp), intent(out) :: kap ! opacity
          real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
          real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
          integer, intent(out) :: ierr ! 0 means AOK.
         

          !internals
          double precision, parameter   :: C_atm=4.786300923226380e-08
          double precision, parameter   :: a_atm=0.68
          double precision, parameter   :: b_atm=0.45
          double precision, parameter   :: dV_s=1d-4 ! small percentage difference in variable (T/rho)
          double precision              :: den, temp, mu_atm ! density, temperature, mean particle weigth
          double precision              :: kap_cons ! prefactor in front (C*(kb/mu)^a)
          double precision              :: d_kap,g_kap
          double precision              :: dtgr=0.01
          double precision              :: T1,T2,rho1,rho2,kap1,kap2
          type (star_info), pointer     :: s

          !Use the Rogers & Seager (2010) fit to the opacity table
          !kappa=C*P^b*T^b

          ierr=0
          call get_star_ptr(id, s, ierr)
          if (ierr /= 0) return

          ! find mu
          mu_atm=(s% mu(1))*amu ! assume atmosphere is atomic 
          ! find kap_cons
          kap_cons=(C_atm*(kerg/mu_atm)**a_atm)
          
          !calculate opacity and derivatives
          den=10d0**log10_rho
          temp=10d0**log10_T
          
          kap=kap_cons*(den**a_atm)*(temp**(a_atm+b_atm))

          !derivatives
          dln_kap_dlnRho=a_atm
          dln_kap_dlnT  =a_atm+b_atm

          if ((s% x_logical_ctrl(4))) then
             ! include PPD dust opacity
             !read in dust to gas ratio
             dtgr=(s% x_ctrl(4))
             call get_BL_opacity(den,temp,d_kap,dtgr)
             g_kap=kap
             kap=kap+d_kap
             !now calculate update derivatives
             ! use simple finite difference method
             !temperature first
             T1=temp*(1d0-dV_s/2d0)
             T2=temp*(1d0+dV_s/2d0)
             kap1=kap_cons*(den**a_atm)*(T1**(a_atm+b_atm))
             kap2=kap_cons*(den**a_atm)*(T2**(a_atm+b_atm))
             call get_BL_opacity(den,T1,d_kap,dtgr)
             kap1=kap1+d_kap
             call get_BL_opacity(den,T2,d_kap,dtgr)
             kap2=kap2+d_kap
             dln_kap_dlnT=(temp/kap)*((kap2-kap1)/(temp*dV_s))
             !now density
             rho1=den*(1d0-dV_s/2d0)
             rho2=den*(1d0+dV_s/2d0)
             kap1=kap_cons*(rho1**(a_atm))*(temp**(a_atm+b_atm))
             kap2=kap_cons*(rho2**(a_atm))*(temp**(a_atm+b_atm))
             call get_BL_opacity(rho1,temp,d_kap,dtgr)
             kap1=kap1+d_kap
             call get_BL_opacity(rho2,temp,d_kap,dtgr)
             kap2=kap2+d_kap
             dln_kap_dlnRho=(den/kap)*((kap2-kap1)/(den*dV_s))
          endif


       end subroutine planetary_kap1

       subroutine get_BL_opacity(rho,T,kappa,dtgr)

        implicit none

        integer                        ::      i,j

        double precision               ::      rho,T,kappa
        double precision, dimension(8) ::      Ttran ! transition temperature between different power-laws
        double precision, parameter    ::      dtgr_BL=0.01 ! dust to gas ratio of Bell & Lin opacities
        double precision               ::      dtgr !actual dust to gas ratio

        ! calculates the opacity for a given density and Temperature

        !find transition temperatures for this density
        ! moving up the Bell & Lin table to postion 8 is dust only

        !check to see if opacities allocated
        if (.not. allocated(kn)) then
           allocate(kn(12))
           allocate(an(12))
           allocate(bn(12))
           call setup_BL_opacity
        endif

        Ttran=0.
        j=8
        do i=1,7
           Ttran(i)=(kn(i)/kn(i+1))**(1d0/(bn(i+1)-bn(i)))*rho**((an(i)-an(i+1))/(bn(i+1)-bn(i)))
           if (T .lt. Ttran(i)) then
              j=i
              exit
           endif
        end do

        !calculate opacity
        kappa=kn(j)*rho**an(j)*T**bn(j)
        !adjust due to dust-to-gas ratio
        kappa=kappa*(dtgr/dtgr_BL)

      end subroutine get_BL_opacity



      integer function extras_startup(s, id, restart, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         ierr = 0
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         double precision    :: safety_fact1,safety_fact2
         double precision    :: Mp,Rp,Tsurf, Rs
         double precision    :: Rhill,Mstellar,sep


         extras_check_model = keep_going         
  
         !pre-assign constants for Hill Sphere
         sep=1.5d15 !100 AU
         Mstellar=msol 

         if ((s% x_ctrl(12)) .gt. 0d0) then
            sep = (s% x_ctrl(12))*au_to_cm
         endif
         if ((s% x_ctrl(13)) .gt. 0d0) then
            Mstellar= (s% x_ctrl(13)) * msol
         endif
        
         extras_check_model = keep_going         
         !need to check if surface temperature close to escape temperature
         ! do this by comparing planet radius to parker/bondi radius
         Mp=(s% mstar)
         Rp=(10**(s% log_surface_radius))*rsol
         Tsurf=10**(s% log_surface_temperature)
        
         Rs = standard_cgrav*Mp*amu*(s% mu(1))/(2*kerg*Tsurf)

         safety_fact1=1d0
         safety_fact2=1d0
        
         if ((s% x_ctrl(9)) .gt. 0d0) then
            safety_fact1=s% x_ctrl(9)
         endif
         if ((s% x_ctrl(10)) .gt. 0d0) then
            safety_fact2=s% x_ctrl(10)
         endif
         

         if (Rp .gt. safety_fact1*Rs) then
            !stop model
            extras_check_model = terminate
            write(*,*) "Planet not hydrostatic, escape temperature too low, STOP"
            return
         endif
         
         if ((s% x_logical_ctrl(5))) then
            ! then also check against Rp > Rhill
            Rhill=sep*(Mp/(3*Mstellar))**(1d0/3d0)
            if (Rp .gt. safety_fact2*Rhill) then
               !stop model
               extras_check_model = terminate
               write(*,*) "Planet bigger than its Hill Sphere, STOP"
               return
            endif
         endif
         
      end function extras_check_model


      integer function how_many_extra_history_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         
         !note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         
         ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         
         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: ierr
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_logfile_now = .true.
            
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(s, id, id_extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info


!cccccccccccccc Subroutines associated with special functions copied from 
! http://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
! using two ENXA for evalution of En(X) and e1xb for evaluation of E1(x)

      subroutine enxa ( n, x, en )

!*****************************************************************************80
!
!! ENXA computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EN(0:N), the function values.
!
        implicit none
        
        integer ( kind = 4 ) n

        real ( kind = 8 ) e1
        real ( kind = 8 ) ek
        real ( kind = 8 ) en(0:n)
        integer ( kind = 4 ) k
        real ( kind = 8 ) x
        
        en(0) = exp ( - x ) / x 
        call e1xb ( x, e1 )
        
        en(1) = e1
        do k = 2, n
           ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
           en(k) = ek
           e1 = ek
        end do

        return
      end subroutine enxa

      subroutine e1xb ( x, e1 )

!*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
        implicit none

        real ( kind = 8 ) e1
        real ( kind = 8 ) ga
        integer ( kind = 4 ) k
        integer ( kind = 4 ) m
        real ( kind = 8 ) r
        real ( kind = 8 ) t
        real ( kind = 8 ) t0
        real ( kind = 8 ) x
        
        if ( x == 0.0D+00 ) then
           
           e1 = 1.0D+300

        else if ( x <= 1.0D+00 ) then

           e1 = 1.0D+00
           r = 1.0D+00

           do k = 1, 25
              r = -r * k * x / ( k + 1.0D+00 )**2
              e1 = e1 + r
              if ( abs ( r ) <= abs ( e1 ) * 1.0D-15 ) then
                 exit
              end if
           end do
    
           ga = 0.5772156649015328D+00
           e1 = - ga - log ( x ) + x * e1

        else

           m = 20 + int ( 80.0D+00 / x )
           t0 = 0.0D+00
           do k = m, 1, -1
              t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
           end do
           t = 1.0D+00 / ( x + t0 )
           e1 = exp ( -x ) * t
    
        end if

        return
      end subroutine e1xb

      subroutine setup_BL_opacity

        implicit none

        !reads the opacity varibles into memory from the Bell & Lin (1997) disc opacities

        ! of the form kappa=kn*rho^an*T^bn
        ! upto 8 is mainly dust dominated

        kn(1)=1e-4
        an(1)=0.
        bn(1)=2.1
        kn(2)=3d0
        an(2)=0.
        bn(2)=-0.01
        kn(3)=1e-2
        an(3)=0d0
        bn(3)=1.1
        kn(4)=5e4
        an(4)=0d0
        bn(4)=-1.5
        kn(5)=1e-1
        an(5)=0d0
        bn(5)=0.7
        kn(6)=2e15
        an(6)=0d0
        bn(6)=-5.2
        kn(7)=2e-2
        an(7)=0d0
        bn(7)=0.8
        kn(8)=2d81
        an(8)=1d0
        bn(8)=-24d0
        kn(9)=1e-8
        an(9)=2./3.
        bn(9)=3d0
        kn(10)=1e-36
        an(10)=1./3.
        bn(10)=10d0
        kn(11)=1.5e20
        an(11)=1d0
        bn(11)=-2.5
        kn(12)=0.348
        an(12)=0.0
        bn(12)=0.0

      end subroutine setup_BL_opacity
      

      end module run_star_extras
      
! list of X_ctrls

! logicals
! 1- core heating from radio-active decay
! 2- core heating from heat-capacity
! 3- envelope heating towards adiabat
! 4- include PPD dust opacity
! 5- check for exit on Hill sphere criterion
! 6- use core luminosity as heating mechnism
! 7- use evaporation
! 8- follow stellar evolution in Tstar and Rstar

! reals
! 1- Atmosphereic optical depth for PPD boundary 
! 2- Atmosphreic Pressure for PPD boundary [cgs]
! 3- Atmosphereic Temperature for PPD boundary [K]
! 4- dust-to-gas mass ratio for using PPD dust opacties
! 5- LEAVE BLANK
! 6- desired entropy for envelope heating [kerg per Baryon]
! 7- cooling time for envelope heating [years]
! 8- LEAVE BLANK
! 9- safety_factor for planet not being bound
! 10- safety_factor for planet not fitting in Hill sphere
! 11- LEAVE BLANK
! 12- sep [AU]
! 13- Mstar [Msun]
! 14- LEAVE BLANK
! 15- kappa_v ! optical depth to incoming stellar radiation (option2) or gamma (option3) 
! 16- atmos/inter tau boundary factor set boundary at tau=fact*(2/(3*gamma))
! 17- Tmin for bi-section method to solve for boundary
! 18- Tmax for bi-section method to solve for boundary
! 19- Equilibruim temperature of planet
! 20- convergence criterion for atmos/interior boundary
! 21- Leave blank
! 22- L/M in cgs units for core luminosity
! 23- Leave blank
! 24- Saturation Time for X-rays [years]
! 25- Lx/Lbol for X-rays during saturated period
! 26- Fall off power-law after Tsat 

! integers
! 1 - type of atmosphere to use 1- disc, 2- guillot two stream
! 2 - max number of iterations for bi-section atmosphere solver
