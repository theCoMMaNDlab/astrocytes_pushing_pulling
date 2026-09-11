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
         ! write(6, *) 'Material: ', cmname(1:8), totalTime, dt
         call  vumat_white (jblock(i_umt_nblock),
     +      ndir, nshr, nstatev, nprops, 
     +      stepTime, totalTime, dt, coordMp,
     +      props, density, strainInc, 
     +      stressOld, stateOld, enerInternOld, enerInelasOld,
     +      stretchNew, defgradNew, 
     +      stressNew, stateNew, enerInternNew, enerInelasNew,
     +      jblock(i_umt_noel))
      else if(cmname(1:4) .eq. 'GRAY') then
         ! write(6, *) 'Material: ', cmname(1:8), totalTime, dt
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


      !  if (totalTime.lt.0.1) then
      !     do km=1,nblock
      !        inicoord(nElement(km),1) = coordMp(km,1)
      !        inicoord(nElement(km),2) = coordMp(km,2)
      !        inicoord(nElement(km),3) = coordMp(km,3)
      !        if (km.eq.1 .or. km.eq.nblock) then
      !           write(6,*) 'POPULATE: nElement=',nElement(km),
      ! +              ' coordMp=',coordMp(km,1),coordMp(km,2),coordMp(km,3)
      !           flush(6)
      !        endif
      !     enddo
      !  end if


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
         ! write(6, *) 'zeta: ',coordx
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
      real*8 thetag_t, thetag_tau, sigma_tau(3,3)
      real*8 coordx, coordy, coordz
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

         ! update state variables
         call mdet(F_tau,detF)
         stateNew(km,1) = thetag_tau ! growth parameter at end of increment
         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz
         stateNew(km,5) = detF
         stateNew(km,6) = (sigma_tau(1,1)+sigma_tau(2,2)+sigma_tau(3,3))/three

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
           
         ! enerInternNew(km) = enerInternOld(km) + stress_power/density(km)
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

      integer nprops

      ! arguments passed in to read
      real*8 props(nprops), F_tau(3,3), dtime, totalTime
      real*8 coordx, coordy, coordz, thetag_t
      ! arguments passed in to write
      real*8 theta_dot_1, f_2, thetag_tau, sigma_tau(3,3), W_e, Je
      ! properties
      real*8 lambda, mu, G_GM, gamma_1, gamma, T_1, T_2
      real*8 maj_axis, min_axis, b_tilde, N_gyri
      real*8 alpha, delta, alpha_bar, delta_bar
      ! local quantities
      real*8 a_tilde, rad, coordiff, r_tilde, psi
      real*8 Fg_tau(3,3), Jg, Fginv(3,3)
      real*8 Fe_tau(3,3), Be_tau(3,3), lnJe
      real*8 detF, f_phi, f_H
      real*8 Iden(3,3), zero, one, two, three, four, half, third

      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three = 3.d0, four=4.d0,
     +     half=0.5d0, third=1.d0/3.d0)
      call onem(Iden)

      ! obtain material properties
      mu        = props(1)
      lambda    = props(2)
      G_GM      = props(3)  ! gray matter growth rate
      gamma_1   = props(4)  ! G_wm/G_gm used in Phase 1
      gamma     = props(5)  ! G_wm/G_gm in Phase 3 (pushing)
      T_1       = props(6)  ! end of Phase 1
      T_2       = props(7)  ! end of Phase 2
      maj_axis  = props(8)  ! WM major axis (a)
      min_axis  = props(9)  ! WM minor axis (b)
      b_tilde   = props(10) ! scaling factor used in calculating r_tilde
      N_gyri    = props(11) ! number of proliferation zones
      alpha     = props(12) ! standard deviation of Gauss function
      delta     = props(13) ! scaled threshold for Gauss function
      alpha_bar = props(14) ! smoothing for Heaviside function used in Phase 3
      delta_bar = props(15) ! scaled threshold for Heaviside function


      ! write(6,*) 'PROPS CHECK: mu=',mu,' lambda=',lambda,' G_GM=',G_GM
      ! write(6,*) '  gamma_1=',gamma_1,' T_1=',T_1,' T_2=',T_2
      ! write(6,*) '  maj_axis=',maj_axis,' min_axis=',min_axis,' b_tilde=',b_tilde
      ! write(6,*) '  N_gyri=',N_gyri,' alpha=',alpha,' delta=',delta
      ! flush(6)
      if (totalTime .gt. 0.20d0 .and. totalTime .lt. 0.20002d0) then
         write(6,*) 'IW nprops =', nprops
         write(6,*) 'IW props  =', props(1), props(2), props(3),
     +                             props(10), props(11), props(12),
     +                             props(13), props(14), props(15)
         flush(6)
      endif
      call mdet(F_tau,detF)

!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (dtime.lt.zero) then
         thetag_tau = thetag_t
         Fg_tau  = (thetag_tau**third)*Iden
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

      psi = atan(coordy/coordx)
      rad = sqrt(coordx**2.0 + coordy**2.0)
      r_tilde = rad/sqrt(
     +     ((maj_axis - a_tilde)*cos(psi))**2.0 
     +   + ((min_axis - b_tilde)*sin(psi))**2.0 )

      ! calculate growth rate in Phase 1
      call gauss(r_tilde,delta,alpha,f_phi)
      f_2 = sin(four*psi*(N_gyri - half)) + one
      theta_dot_1 = (G_GM*gamma_1)*half*f_phi*f_2 

      if (totalTime.le.T_1) then ! Phase 1
         thetag_tau = thetag_t + (theta_dot_1)*dtime

      else if (totalTime.le.T_2) then ! Phase 2 (WM doesn't grow)
         thetag_tau = thetag_t

      else ! Phase 3 (pushing)

         ! calculate growth rate in Phase 3
         coordiff = (r_tilde - delta_bar)*one
         call Hhat(coordiff,alpha_bar,f_H)
         theta_dot_1 = (G_GM*gamma)*half*f_H*f_2
         thetag_tau = thetag_t + (theta_dot_1)*dtime

      endif
      if (thetag_tau .le. 1.0d-6) thetag_tau = 1.0d-6 
      ! update kinematics
      Fg_tau  = (thetag_tau**third)*Iden
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
     +                      thetag_tau, sigma_tau,W_e, Je)

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
      real*8 detF, N_R(3,1), tmp
      real*8 Fg_tau(3,3), Jg, Fginv(3,3)
      real*8 Fe_tau(3,3), Be_tau(3,3), lnJe
      real*8 Iden(3,3), zero, one, two, three, half, third

      ! constants
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, half=0.5d0, third=1.d0/3.d0)
      call onem(Iden)

      ! obtain material properties
      mu       = props(1)
      lambda   = props(2)
      G_GM     = props(3) ! gray matter growth rate
      T_1      = props(4) ! end of Phase 1, beginning of Phase 2
      maj_axis = props(5) ! WM major axis (a)
      min_axis = props(6) ! WM minor axis (b)


      call mdet(F_tau,detF)

!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (dtime.lt.zero) then
         thetag_tau = thetag_t
         Fg_tau = thetag_tau**third*Iden
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
         return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! obtain referential surface outnormal of an elliptical surface
      N_R(1,1) = 2.0*coordx/maj_axis**2.0
      N_R(2,1) = 2.0*coordy/min_axis**2.0
      N_R(3,1) = zero
      tmp = sqrt(N_R(1,1)**2.0 + N_R(2,1)**2.0 + N_R(3,1)**2.0)
      N_R = N_R/tmp

      if (totalTime.le.T_1) then ! Phase 1 (no GM growth)
         thetag_tau = one
      else ! Phase 2 and 3
         thetag_tau = thetag_t + (G_GM)*dtime
      endif
      if (thetag_tau .le. 1.0d-6) thetag_tau = 1.0d-6 
      ! update kinematics
      Fg_tau = dsqrt(thetag_tau)*Iden
     +         +(one - dsqrt(thetag_tau))*matmul(N_R,transpose(N_R))
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
      end subroutine integ_gray

      !****************************************************************

      subroutine Hhat(x, alpha, y)
      ! smooth Heaviside approximation: y = exp(ax)/(1+exp(ax))

      implicit none
      real*8 x,alpha,y

      y = (exp(alpha*x))/(1.0 + exp(alpha*x))

      end subroutine Hhat

      !****************************************************************

      subroutine gauss(x, mean, sigma, y)
      ! Gaussian curve evaluated at x with given mean and sigma

      implicit none
      real*8 x, mean, sigma, y
      real*8 half, two, Pi
      parameter(half=0.5d0, two=2.d0, Pi=3.1415926d0)

      y = exp(-half * ((x - mean) / sigma)**2) / (sigma * dsqrt(two*Pi))

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
