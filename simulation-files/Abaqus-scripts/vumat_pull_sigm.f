      ! author: Karan Taneja
      ! unit system: mm-N-MPa
      ! dimension: plane strain
      ! documentation: 'Astrocytes in white matter ...'

      !****************************************************************

      module GlobalStorage
      ! this number must >= the number of elements in the mesh
      ! otherwise, the simulation will crash   
      real*8 inicoord(46000,3)
      end module  

      !****************************************************************

      subroutine vumat (
            ! read only
     +      jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, 
     +      stepTime, totalTime, dt, cmname, coordMp, charLength, 
     +      props, density, strainInc, relSpinInc,
     +      tempOld, stretchOld, defgradOld, fieldOld,
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      tempNew, stretchNew, defgradNew, fieldNew,
            ! write only
     +      stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

      dimension jblock(*), coordMp(*), charLength(*), 
     +      props(nprops), density(*), strainInc(*),
     +      relSpinInc(*), tempOld(*), stretchOld(*),
     +      defgradOld(*), fieldOld(*), stressOld(*),
     +      stateOld(*), enerInternOld(*), enerInelasOld(*), 
     +      tempNew(*), stretchNew(*), defgradNew(*),
     +      fieldNew(*), stressNew(*), stateNew(*),
     +      enerInternNew(*), enerInelasNew(*)

      character*80 cmname

      parameter (     
     +      i_umt_nblock = 1,
     +      i_umt_npt    = 2,
     +      i_umt_layer  = 3,
     +      i_umt_kspt   = 4,
     +      i_umt_noel   = 5 )

      ! call particular user material to perform the analysis 
      if (cmname(1:4) .eq. 'WHIT') then
         call  vumat_white (jblock(i_umt_nblock),
     +      ndir, nshr, nstatev, nprops, 
     +      stepTime, totalTime, dt, coordMp,
     +      props, density, strainInc, 
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      stretchNew, defgradNew, 
     +      stressNew, stateNew, enerInternNew, enerInelasNew,
     +      jblock(i_umt_noel))
      else if(cmname(1:4) .eq. 'GRAY') then
         call  vumat_gray (jblock(i_umt_nblock),
     +      ndir, nshr, nstatev, nprops, 
     +      stepTime, totalTime, dt, coordMp,
     +      props, density, strainInc, 
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      stretchNew, defgradNew, 
     +      stressNew, stateNew, enerInternNew, enerInelasNew,
     +      jblock(i_umt_noel))
      endif

      end subroutine vumat

      !****************************************************************

      subroutine vumat_white (
            ! read only
     +      nblock, ndir, nshr, nstatev, nprops, 
     +      stepTime, totalTime, dt, coordMp,
     +      props, density, strainInc, 
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      stretchNew, defgradNew, 
            ! write only
     +      stressNew, stateNew, enerInternNew, enerInelasNew,
            ! extra argument - array of internal element numbers
     +      nElement )

      use GlobalStorage
      include 'vaba_param.inc'
      ! implicit none 

      dimension coordMp(nblock,*), 
     +      props(nprops), density(nblock), strainInc(nblock,ndir+nshr),
     +      stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev), 
     +      enerInternOld(nblock), enerInelasOld(nblock), 
     +      stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr),
     +      stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +      enerInternNew(nblock), enerInelasNew(nblock), 
     +      nElement(nblock)

      ! local quantities
      integer i, km
      real*8 F_tau(3,3), detF, stress_power
      real*8 R_tau(3,3), U_tau(3,3), U_inv(3,3)
      real*8 N_R(2,1), zeta, rot_matrix(3,3)
      real*8 sigma_tau(3,3), sigma_rot(3,3), sigma_rad, sigma_tan
      real*8 theta_dot_1, f_2, thetag_t, thetag_tau
      real*8 coordx, coordy, coordz, maj_axis, min_axis
      real*8 zero, one, two, three, half, third, Pi
      real*8 W_e, Je

      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0,
     +   half=0.5d0, third=1.d0/3.d0, Pi=3.1415926d0)

      ! WM ellipse parameters
      maj_axis = props(8) ! WM major axis (mm)
      min_axis = props(9) ! WM minor axis (mm)



      ! pour initial coordinates into the global variable matrix 
      if (totalTime.lt.1e-4) then
         do km=1,nblock
            inicoord(nElement(km),1) = coordMp(km,1)
            inicoord(nElement(km),2) = coordMp(km,2)
            inicoord(nElement(km),3) = coordMp(km,3)
         enddo
      end if 

      ! start loop over material points:
      do km=1,nblock

         ! copy new deformation gradient (at end of current increment)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         
         if(nshr .lt. 2) then ! 2D case
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else ! 3D case
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if

         ! initialization (dummy) step: set growth parameter at t=0
         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            stateOld(km,1) = one 
         endif

         ! read old growth parameter and coordinates
         thetag_t = stateOld(km,1) ! growth parameter at beginning of increment
         coordx = inicoord(nElement(km),1)
         coordy = inicoord(nElement(km),2)
         coordz = inicoord(nElement(km),3)



         !-------------------------------------------------------------
         ! perform the time integration and compute the
         ! constitutive response based on the material model

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! initialization (dummy) step: pass zero timestep
            call integ_white(props, nprops, F_tau, -1.0d0, totalTime,
     +                       coordx, coordy, coordz, thetag_t,
     +                       theta_dot_1, f_2, thetag_tau, sigma_tau, W_e, Je)
         else
            ! perform explicit time integration procedure
            call integ_white(props, nprops, F_tau, dt, totalTime,
     +                       coordx, coordy, coordz, thetag_t,
     +                       theta_dot_1, f_2, thetag_tau, sigma_tau, W_e, Je)
         endif
         !-------------------------------------------------------------

         ! define stress at end of current increment
         ! (ABAQUS/Explicit uses Cauchy stress in corotational frame, R^T . sigma . R
         ! see ABAQUS User Subroutines Manual, VUMAT, Special Considerations for Hyperelasticity

         call m3inv(U_tau, U_inv)
         R_tau = matmul(F_tau, U_inv)
         sigma_tau = matmul(transpose(R_tau), matmul(sigma_tau, R_tau))

         do i=1,ndir
            stressNew(km,i) = sigma_tau(i,i)
         end do
         
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = sigma_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = sigma_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = sigma_tau(1,3)
               endif
            endif
         endif

         ! rotate stress tensor and get radial and tangential outputs
         ! S' = R.S.R^T
         N_R(1,1) = 2.0*coordx/maj_axis**2.0
         N_R(2,1) = 2.0*coordy/min_axis**2.0
         zeta = atan(N_R(2,1)/N_R(1,1))
         ! write(6, *) 'zeta: ',zeta

         rot_matrix = zero
         rot_matrix(1,1) =  cos(zeta)
         rot_matrix(1,2) =  sin(zeta)
         rot_matrix(2,1) = -sin(zeta)
         rot_matrix(2,2) =  cos(zeta)
         rot_matrix(3,3) =  one

         sigma_rot = matmul(matmul(rot_matrix, sigma_tau), transpose(rot_matrix))
         sigma_rad = sigma_rot(1,1)
         sigma_tan = sigma_rot(2,2)

         ! update state variables
         call mdet(F_tau,detF)

         stateNew(km,1) = thetag_tau ! growth parameter at end of increment
         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz 
         stateNew(km,5) = detF   
         stateNew(km,6) = (sigma_tau(1,1)+sigma_tau(2,2)+sigma_tau(3,3))/three
         stateNew(km,7) = theta_dot_1
         stateNew(km,8) = f_2
         stateNew(km,9) = sigma_rad
         stateNew(km,10) = sigma_tan   

         ! update the specific internal energy
         stress_power = zero
         do i = 1,ndir
            stress_power = stress_power + half*(
     +         (stressOld(km,i) + stressNew(km,i))*strainInc(km,i) )
         enddo
         
         select case (nshr)
         case(1) ! 2D analysis
            stress_power = stress_power + half*(
     +         (stressOld(km,ndir+1) + stressNew(km,ndir+1))*strainInc(km,ndir+1) )
         case(3) ! 3D analysis
            stress_power = stress_power + half*(
     +         (stressOld(km,ndir+1) + stressNew(km,ndir+1))*strainInc(km,ndir+1) +
     +         (stressOld(km,ndir+2) + stressNew(km,ndir+2))*strainInc(km,ndir+2) +
     +         (stressOld(km,ndir+3) + stressNew(km,ndir+3))*strainInc(km,ndir+3) )
         end select
           
         enerInternNew(km) = W_e/(density(km)*Je)
         enerInelasNew(km) = enerInelasOld(km)

      enddo ! end loop over material points
      end subroutine vumat_white

      !****************************************************************

      subroutine vumat_gray (
            ! read only
     +      nblock, ndir, nshr, nstatev, nprops, 
     +      stepTime, totalTime, dt, coordMp, 
     +      props, density, strainInc, 
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      stretchNew, defgradNew,
            ! write only
     +      stressNew, stateNew, enerInternNew, enerInelasNew,
            ! read only extra arguments
     +      nElement )

      use GlobalStorage
      include 'vaba_param.inc'
      ! implicit none 

      dimension coordMp(nblock,*),
     +      props(nprops), density(nblock), strainInc(nblock,ndir+nshr),
     +      stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev), 
     +      enerInternOld(nblock), enerInelasOld(nblock), 
     +      stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr),
     +      stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +      enerInternNew(nblock), enerInelasNew(nblock),
     +      nElement(nblock)

      ! local quantities
      integer i, km
      real*8 F_tau(3,3), detF, stress_power
      real*8 R_tau(3,3), U_tau(3,3), U_inv(3,3)
      real*8 N_R(2,1), zeta, rot_matrix(3,3)
      real*8 sigma_tau(3,3), sigma_rot(3,3), sigma_rad, sigma_tan
      real*8 thetag_t, thetag_tau
      real*8 coordx, coordy, coordz, maj_axis, min_axis
      real*8 zero, one, two, three, half, third
      real*8 W_e, Je

      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, 
     +     half=0.5d0, third=1.d0/3.d0)

      ! pour initial coordinates into the global variable
      if (totalTime.lt.1e-4) then
         do km=1,nblock
            inicoord(nElement(km),1) = coordMp(km,1)
            inicoord(nElement(km),2) = coordMp(km,2)
            inicoord(nElement(km),3) = coordMp(km,3)
         enddo
      end if

      maj_axis = props(5)   ! pial major axis
      min_axis = props(6)   ! pial minor axis

      ! start loop over material points:
      do km=1,nblock

         ! copy new deformation gradient (at end of current increment)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         
         if(nshr .lt. 2) then ! 2D case
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else ! 3D case
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if

         ! initialization (dummy) step: set growth parameter at t=0
         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            stateOld(km,1) = one 
         endif

         ! read old state variables and coordinates
         thetag_t = stateOld(km,1) ! growth parameter at beginning of increment
         coordx = inicoord(nElement(km),1)
         coordy = inicoord(nElement(km),2)
         coordz = inicoord(nElement(km),3)

         !---------------------------------------------------------------
         ! perform the time integration and compute the
         ! constitutive response based on the material model

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! initialization (dummy) step: pass zero timestep
            call integ_gray(props, nprops, F_tau, -1.0d0, totalTime,
     +                      coordx, coordy, coordz, thetag_t,
     +                      thetag_tau, sigma_tau, W_e, Je)
         else
            ! perform explicit time integration procedure
            call integ_gray(props, nprops, F_tau, dt, totalTime,
     +                      coordx, coordy, coordz, thetag_t,
     +                      thetag_tau, sigma_tau, W_e, Je)
         endif
         !---------------------------------------------------------------
         
         ! define stress at end of current increment
         ! (ABAQUS/Explicit uses Cauchy stress in corotational frame, R^T . sigma . R
         ! see ABAQUS User Subroutines Manual, VUMAT, Special Considerations for Hyperelasticity

         call m3inv(U_tau, U_inv)
         R_tau = matmul(F_tau, U_inv)
         sigma_tau = matmul(transpose(R_tau), matmul(sigma_tau, R_tau))

         do i=1,ndir
            stressNew(km,i) = sigma_tau(i,i)
         end do
         
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = sigma_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = sigma_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = sigma_tau(1,3)
               endif
            endif
         endif

         ! rotate stress tensor and get radial and tangential outputs
         ! S' = R.S.R^T
         N_R(1,1) = 2.0*coordx/maj_axis**2.0
         N_R(2,1) = 2.0*coordy/min_axis**2.0
         zeta = atan(N_R(2,1)/N_R(1,1))
         ! write(6, *) 'zeta: ',zeta

         rot_matrix = zero
         rot_matrix(1,1) =  cos(zeta)
         rot_matrix(1,2) =  sin(zeta)
         rot_matrix(2,1) = -sin(zeta)
         rot_matrix(2,2) =  cos(zeta)
         rot_matrix(3,3) =  one

         sigma_rot = matmul(matmul(rot_matrix, sigma_tau), transpose(rot_matrix))
         sigma_rad = sigma_rot(1,1)
         sigma_tan = sigma_rot(2,2)

         ! update state variables
         call mdet(F_tau,detF)

         stateNew(km,1) = thetag_tau ! growth parameter at end of increment
         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz
         stateNew(km,5) = detF
         stateNew(km,6) = (sigma_tau(1,1)+sigma_tau(2,2)+sigma_tau(3,3))/three
         if (nstatev.ge.10) then
            stateNew(km,7)  = zero
            stateNew(km,8)  = zero
            stateNew(km,9)  = sigma_rad
            stateNew(km,10) = sigma_tan
         endif
         ! update the specific internal energy
         stress_power = zero
         do i = 1,ndir
            stress_power = stress_power +
     +         half*((stressOld(km,i)+stressNew(km,i))*
     +         strainInc(km,i))
         enddo
         
         select case (nshr)
         case(1) ! 2D analysis
            stress_power = stress_power + half*(
     +         (stressOld(km,ndir+1) + stressNew(km,ndir+1))*strainInc(km,ndir+1) )
         case(3) ! 3D analysis
            stress_power = stress_power + half*(
     +         (stressOld(km,ndir+1) + stressNew(km,ndir+1))*strainInc(km,ndir+1) +
     +         (stressOld(km,ndir+2) + stressNew(km,ndir+2))*strainInc(km,ndir+2) +
     +         (stressOld(km,ndir+3) + stressNew(km,ndir+3))*strainInc(km,ndir+3) )
         end select
           
         enerInternNew(km) = W_e/(density(km)*Je)
         enerInelasNew(km) = enerInelasOld(km)

      enddo ! end loop over material points
      end subroutine vumat_gray

      !****************************************************************

      subroutine integ_white(props, nprops, F_tau, dtime, totalTime,
     +                       coordx, coordy, coordz, thetag_t,
     +                       theta_dot_1, f_2, thetag_tau, sigma_tau, W_e, Je
     +                       )

      implicit none

      integer nitl, nprops

      ! arguments passed in to read
      real*8 props(nprops), F_tau(3,3), dtime, totalTime
      real*8 coordx, coordy, coordz, thetag_t
      ! arguments passed in to write
      real*8 theta_dot_1, f_2, thetag_tau, sigma_tau(3,3), W_e, Je
      ! properties
      real*8 lambda, mu, G_GM, gamma_1, gamma_hat, T_1, T_2
      real*8 maj_axis, min_axis, b_tilde, N_gyri, alpha, delta
      ! local quantities
      real*8 Fg_tau(3,3), Jg, Fginv(3,3)
      real*8 Fe_tau(3,3), Be_tau(3,3)
      real*8 detF, lnJe, trMe, f_phi
      real*8 a_tilde, rad, r_tilde, psi
      real*8 growth_crit, res, dres, phig, dphig, xtol
      real*8 thetag_dum, theta_dot_3_fac
      real*8 Iden(3,3), zero, one, two, three, four, nine, half, third
      real*8 rho_sig, R_g, theta_max, g_log
      real*8 day_offset, days_sim


      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0, nine=9.d0,
     +   half=0.5d0, third=1.d0/3.d0)
      parameter(days_sim = 1.0d0)    ! 18 days of step time
      parameter(day_offset = -2.0d0)  ! t = 0 is E40 (P-2)
      xtol = 1.d-10 ! tolerance for local Newton iterations in growth parameter calculations
      growth_crit = zero ! critical growth criterion (purely tensile growth in this study)
      call onem(Iden)

      ! obtain material properties
      mu        = props(1)
      lambda    = props(2)
      G_GM      = props(3)  ! gray matter growth rate
      gamma_1   = props(4)  ! G_wm/G_gm used in Phase 1
      gamma_hat = props(5)  ! G_wm/G_gm in Phase 3 (pulling)
      T_1       = props(6)  ! end of Phase 1
      T_2       = props(7)  ! end of Phase 2
      maj_axis  = props(8)  ! WM major axis (a)
      min_axis  = props(9)  ! WM minor axis (b)
      b_tilde   = props(10) ! scaling factor used in calculating r_tilde
      N_gyri    = props(11) ! number of proliferation zones
      alpha     = props(12) ! standard deviation of Gauss function
      delta     = props(13) ! scaled threshold for Gauss function

      ! calculate growth rate in Phase 3 (pulling)
      call mdet(F_tau,detF)
!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (dtime.lt.zero) then
         thetag_tau = thetag_t
         Fg_tau  = (thetag_tau**third)*Iden
         ! Fg_tau  = dsqrt(thetag_tau)*Iden
         ! Fg_tau(3,3) = one
         call m3inv(Fg_tau, Fginv)
         Fe_tau = matmul(F_tau, Fginv)
         Be_tau = matmul(Fe_tau, transpose(Fe_tau))
         call mdet(Fg_tau, Jg)
         Je = detF/Jg
         sigma_tau = ((lambda*dlog(Je) - mu)*Iden + mu*Be_tau)/Je
         lnJe = dlog(Je)
         W_e  = half*mu*(Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)-three)
     +          - mu*lnJe + half*lambda*lnJe*lnJe
         theta_dot_1 = zero
         f_2 = zero
         return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculate other dimensions
      a_tilde = b_tilde*(maj_axis/min_axis)
      ! a_tilde=zero !maj_min_ratio = maj_axis/min_axis is not assigned in the original
      psi = atan(coordy/coordx)
      ! write(6, *) 'psi: ',psi
      rad = sqrt(coordx**2.0 + coordy**2.0)
      r_tilde = rad/sqrt(
     +     ((maj_axis - a_tilde)*cos(psi))**2.0 
     +   + ((min_axis - b_tilde)*sin(psi))**2.0 )

      ! ---- Garcia-calibrated logistic envelope (identical to push) --
      ! Garcia et al. 2024 whole brain: R_g = 0.134/d, t_mid = 20.5 d
      ! 1.0 s of step time <-> 18 postnatal days
      rho_sig   = 0.47d0
      R_g       = 0.134d0*days_sim  ! 1/days
      theta_max = (one + exp(0.134d0*(20.5d0 - day_offset)))
     +            **(two/three)     ! 7.7054

      g_log = (two/three)*R_g*thetag_t
     +        *(one - (thetag_t/theta_max)**1.5d0)
      if (g_log .lt. zero) g_log = zero

      call gauss(r_tilde,delta,alpha,f_phi)
      f_2 = sin(four*psi*(N_gyri - half)) + one

      ! theta_dot_1 = (G_GM*gamma_1)*half*f_phi*f_2+G_const 

      ! phase 3 coefficient
      ! theta_dot_3_fac = gamma_hat*G_GM/mu 

      ! Phases 1+2 (E40 -> T_2 = P6): progenitor (GZ) growth,
      ! patterned peaks f_phi*f_2 plus a constant zone rho_sig
      if (totalTime.le.T_2) then ! Phase 1 and extend to Phase 2
         theta_dot_1 = g_log*gamma_1*(half*f_phi*f_2 + rho_sig)
         thetag_tau  = thetag_t + theta_dot_1*dtime

      else ! Phase 3 (pulling)
         ! thetag_dum = thetag_t
         ! envelope replaces the constant G_GM in the pull gain
         theta_dot_3_fac = gamma_hat*g_log/mu
         thetag_dum = thetag_t

         do nitl = 1, 50

            ! kinematics with current iterate
            Fg_tau  = (thetag_dum**third)*Iden
            ! Fg_tau  = dsqrt(thetag_dum)*Iden
            ! Fg_tau(3,3) = one
            call m3inv(Fg_tau, Fginv)
            Fe_tau = matmul(F_tau, Fginv)
            Be_tau = matmul(Fe_tau, transpose(Fe_tau))
            call mdet(Fg_tau, Jg)
            Je = detF/Jg
            lnJe = log(Je)

            ! growth criterion based on trace of Kirchoff/Mandel stress
            trMe = two*(lambda*lnJe - mu) + mu*(Be_tau(1,1) + Be_tau(2,2) )
            phig = trMe - growth_crit

            if (phig.le.zero) then ! criterion not met: no tensile growth
               thetag_tau = thetag_dum
               exit
            endif

            ! nonlinear residual and derivative
            dphig = -third/thetag_dum*(two*mu*(Be_tau(1,1) + Be_tau(2,2))
     +              + two*three*lambda)
            res  = thetag_dum - thetag_t - (theta_dot_3_fac)*phig*dtime
            dres = one - (theta_dot_3_fac)*dphig*dtime

            ! update iterate
            thetag_dum = thetag_dum - res/dres
            thetag_tau = thetag_dum

            if (dabs(res).le.xtol) exit

            if (nitl.eq.50) print *, '>no convergence!!!! |r|=', dabs(res)

         enddo
         ! SDV7 now reports the ACTUAL pull growth rate
         theta_dot_1 = (thetag_tau - thetag_t)/dtime
      endif

      ! update kinematics
      Fg_tau  = (thetag_tau**third)*Iden
      ! Fg_tau  = dsqrt(thetag_tau)*Iden
      ! Fg_tau(3,3) = one 
      call m3inv(Fg_tau, Fginv)
      Fe_tau = matmul(F_tau, Fginv)
      Be_tau = matmul(Fe_tau, transpose(Fe_tau))
      call mdet(Fg_tau, Jg)
      Je = detF/Jg

      sigma_tau = ((lambda*dlog(Je) - mu)*Iden + mu*Be_tau)/Je
      lnJe = dlog(Je)
      W_e  = half*mu*(Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)-three)
     +          - mu*lnJe + half*lambda*lnJe*lnJe

      end subroutine integ_white

      !****************************************************************

      subroutine integ_gray(props, nprops, F_tau, dtime, totalTime,
     +                      coordx, coordy, coordz, thetag_t,
     +                      thetag_tau, sigma_tau, W_e, Je)

      ! Cortex growth, Eq. (GMGrowth) of the draft:
      !   F^g = theta I + (1 - theta) (N_R x N_R + E_Z x E_Z)
      !   stretch theta along the in-plane tangent, 1 along N_R and E_Z
      !   theta = J^g = sqrt(S(t)/S(P0)), S = cortical surface area
      !   (Garcia et al. 2024)

      implicit none

      ! arguments passed in to read
      integer nprops
      real*8 props(nprops), F_tau(3,3), dtime, totalTime
      real*8 coordx, coordy, coordz, thetag_t
      ! arguments passed in to write
      real*8 thetag_tau, sigma_tau(3,3), W_e, Je
      ! properties
      real*8 mu, lambda, G_GM, T_1, maj_axis, min_axis
      ! local quantities
      real*8 detF, N_R(3,1), E_Z(3,1), NE(3,3), tmp
      real*8 Fg_tau(3,3), Jg, Fginv(3,3)
      real*8 Fe_tau(3,3), Be_tau(3,3), lnJe
      real*8 Iden(3,3), zero, one, two, three, half, third
      real*8 R_sa, tm_sa, t_d, t_on, theta_S
      real*8 day_offset, days_sim


      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, half=0.5d0,
     +          third=1.d0/3.d0)
      ! growth calendar (E40 clock): day = days_sim*t + day_offset
      parameter(days_sim = 1.0d0)     ! step time is in days
      parameter(day_offset = -2.0d0)  ! t = 0 is E40 (P-2)
      call onem(Iden)

      ! obtain material properties
      mu       = props(1)
      lambda   = props(2)
      G_GM     = props(3) ! not used
      T_1      = props(4) ! end of Phase 1 (P0), start of cortical growth
      maj_axis = props(5) ! pial major axis
      min_axis = props(6) ! pial minor axis

      ! N_R: outward unit normal of the ellipse through the point that
      ! has the shape of the pial surface.  E_Z: out-of-plane unit vector.
      N_R(1,1) = two*coordx/maj_axis**2.0
      N_R(2,1) = two*coordy/min_axis**2.0
      N_R(3,1) = zero
      tmp = sqrt(N_R(1,1)**2.0 + N_R(2,1)**2.0 + N_R(3,1)**2.0)
      if (tmp.gt.1.d-12) then
         N_R = N_R/tmp
      else                    ! no position (packaging call at 0,0)
         N_R = zero
         N_R(1,1) = one
      endif
      E_Z = zero
      E_Z(3,1) = one
      NE = matmul(N_R,transpose(N_R)) + matmul(E_Z,transpose(E_Z))

      call mdet(F_tau,detF)
!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (dtime.lt.zero) then
         thetag_tau = thetag_t
         Fg_tau = Iden        ! theta = 1 at the start, no growth yet
         call m3inv(Fg_tau, Fginv)
         Fe_tau = matmul(F_tau, Fginv)
         Be_tau = matmul(Fe_tau, transpose(Fe_tau))
         call mdet(Fg_tau, Jg)
         Je = detF/Jg
         sigma_tau = ((lambda*dlog(Je) - mu)*Iden + mu*Be_tau)/Je
         lnJe = dlog(Je)
         W_e  = half*mu*(Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)-three)
     +          - mu*lnJe + half*lambda*lnJe*lnJe
         return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ---- cortical surface area S(t)/S(P0), Garcia et al. (2024) -----
      ! PLACEHOLDERS read off the figure (they reproduce SA(P8)=351 and
      ! SA(P38)=2260 mm^2); replace with the fitted supplement values.
      R_sa  = 0.161d0                        ! 1/day
      tm_sa = 19.5d0                         ! day
      t_d   = days_sim*totalTime + day_offset   ! postnatal day now
      t_on  = days_sim*T_1 + day_offset         ! = 0 (P0) for T_1 = 2
      if (totalTime.le.T_1) then             ! Phase 1: no GM growth
         theta_S = one
      else
         theta_S = (one + exp(-R_sa*(t_on - tm_sa)))
     +           / (one + exp(-R_sa*(t_d  - tm_sa)))
      endif

      ! growth parameter = J^g = cortical length ratio in the section
      thetag_tau = dsqrt(theta_S)
      Fg_tau = thetag_tau*Iden + (one - thetag_tau)*NE

      ! update kinematics
      call m3inv(Fg_tau, Fginv)
      Fe_tau = matmul(F_tau, Fginv)
      Be_tau = matmul(Fe_tau, transpose(Fe_tau))
      call mdet(Fg_tau, Jg)
      Je = detF/Jg

      sigma_tau = ((lambda*dlog(Je) - mu)*Iden + mu*Be_tau)/Je
      lnJe = dlog(Je)
      W_e  = half*mu*(Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)-three)
     +          - mu*lnJe + half*lambda*lnJe*lnJe
      end subroutine integ_gray


      !****************************************************************

      subroutine gauss(x, mean, sigma, y)
      ! Gaussian curve evaluated at x with given mean and sigma

      implicit none
      real*8 x, mean, sigma, y
      real*8 half, two, Pi
      parameter(half=0.5d0, two=2.d0, Pi=3.1415926d0)

      y = exp(-half * ((x - mean) / sigma)**two) / (sigma * dsqrt(two*Pi))

      end subroutine gauss      

      !****************************************************************

      subroutine onem(A)
      ! set A to the 3x3 identity matrix

      implicit none
      real*8 A(3,3)
        
      A=0d0
      A(1,1)=1d0
      A(2,2)=1d0
      A(3,3)=1d0

      end

      !****************************************************************

      subroutine mtrans(A, Atrans)
      ! compute the transpose of 3x3 matrix A into Atrans.

      implicit none
      real*8 A(3,3), Atrans(3,3)
      integer i, j

      do i = 1, 3
         do j = 1, 3
            Atrans(j,i) = A(i,j)
         end do
      end do

      end

      !****************************************************************

      subroutine mdet(A, det)
      ! compute the determinant of 3x3 matrix A

      implicit none
      real*8 A(3,3), det

      det =   A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     +      + A(1,3)*A(2,1)*A(3,2)
     +      - A(3,1)*A(2,2)*A(1,3) - A(3,2)*A(2,3)*A(1,1)
     +      - A(3,3)*A(2,1)*A(1,2)

      end

      !****************************************************************

      subroutine m3inv(A, Ainv)
      ! compute the inverse of 3x3 matrix A into Ainv; stops if singular

      implicit none
      real*8 A(3,3), Ainv(3,3), det, Acofac(3,3), Aadj(3,3)
      integer i, j

      call mdet(A, det)
      if (det .eq. 0.d0) then
         write(*,*) 'm3inv: matrix is singular'
         stop
      end if

      call mcofac(A, Acofac)
      call mtrans(Acofac, Aadj)

      do i = 1, 3
         do j = 1, 3
            Ainv(i,j) = Aadj(i,j) / det
         end do
      end do

      end

      !****************************************************************

      subroutine mcofac(A, Acofac)
      ! compute the cofactor matrix of 3x3 matrix A for inversion

      implicit none
      real*8 A(3,3), Acofac(3,3)

      Acofac(1,1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
      Acofac(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      Acofac(1,3) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)
      Acofac(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      Acofac(2,2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
      Acofac(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      Acofac(3,1) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)
      Acofac(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      Acofac(3,3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

      end