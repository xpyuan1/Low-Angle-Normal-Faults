c----- Mechanics of LANF
c-----------------------------------------------------------------------
      program LANF
      implicit double precision (a-h,o-z)
      logical found
      data tol/1.0d-3/
c
c---- input material parameters
      read(5,*) cb,phib,beta
      read(5,*) cfmin,cfmax,ncf
      read(5,*) phifmin,phifmax,nphif
      read(5,*) rhof,dlambB,dlambF
      read(5,*) nz,ntheta
c
c---- convert angles in rad
      phib = phib/180.d0*3.14159d0
      beta = beta/180.d0*3.14159d0
      phifmin = phifmin/180.d0*3.14159d0
      phifmax = phifmax/180.d0*3.14159d0

c---- increment in phif
      dphif = (phifmax - phifmin)/real(nphif-1)
c---- loop over all phif
      do 1000 iphif = 1,nphif
c---- current value of phif
      phif = phifmin + dphif*real(iphif-1)
c
c---- the force for locked depth D = 1
      call compext_35(quD1,phib,phif,cb,rhof,dlambB,dlambF)
c
c---- the values of theta1 and gamma1 
      theta1 = phib/2.d0+3.14159d0/4.d0
      gamma1 = phib/2.d0+3.14159d0/4.d0
c
c---- increment in cf
      dcf = (cfmax - cfmin)/real(ncf-1)
c---- loop over all cf
      do 2000 icf = 1,ncf
c---- current value of cf
      cf = cfmin + dcf*real(icf-1)
c
c---- initialize the mini of the function
      qulower = 1.0d+12
      thetalower = 0.d0
c---- set the flag to false
      found = .false.
c
c---- increment of depth z, root of hanging wall
      dz = 1.d0/real(nz-1)
c---- loop over dz from 0.d0 to 1.d0
      do 3000 iz = 1,nz
c---- current value of z
      xz = dz*real(iz-1)
c
c---- define the effective region for theta
      thetamin = (phib-beta+phif)*(1.d0+tol)
      thetamax = (3.14159d0-phib-beta+phif)*(1.d0-tol)
c---- increment in theta
      dtheta = (thetamax - thetamin)/real(ntheta-1)
c---- loop over all theta
      do 4000 itheta = 1,ntheta
c---- current value of theta
      theta = thetamin + dtheta*real(itheta-1)
c---- test the angle constraint
      term1 = (1.d0-xz)/dtan(beta) + 1.d0/dtan(theta)
      term2 = xz/dtan(theta1)
      if(term2.gt.term1) goto 4000
c
c---- velocity of the detachment layer
      if(beta.ge.phif) then
        xUd = dsin(theta-phib)/dsin(theta-phib+beta-phif)
        xJ = dsin(beta-phif)/dsin(theta-phib+beta-phif)
      else
        xUd = dsin(theta+phib)/dsin(theta+phib-phif+beta)
        xJ = dsin(phif-beta)/dsin(theta+phib-phif+beta)
      endif
      if((xUd.lt.0.d0).or.(xJ.lt.0.d0)) goto 4000
c
c---- size of each segment
      xLg1e1 = xz/dsin(gamma1)
      xLg1f1 = xz/dsin(theta1)
      xLgg1 = (1.d0-xz)/dsin(beta)
      xLgf = 1.d0/dsin(theta)
c---- the size of Shw, Sd and Sbs
      Shw = xz*xLg1e1*dcos(gamma1)/2.d0
     1    + xz*xLg1f1*dcos(theta1)/2.d0
      Sd = (xz+1.d0)*xLgg1*dcos(beta)/2.d0
     1   + 1.d0*xLgf*dcos(theta)/2.d0
     1   - xz*xLg1f1*dcos(theta1)/2.d0
      Sbs = -1.d0*xLgf*dcos(theta)/2.d0
c---- the size for fluid overpressure power
      Sg1e1 = xz*xLg1e1*dcos(gamma1)/2.d0
      Sg1f1 = xz*xLg1f1*dcos(theta1)/2.d0
      Sgg1 = (xz+1.d0)*xLgg1*dcos(beta)/2.d0
      Sgf = 1.d0*xLgf*dcos(theta)/2.d0
c---- velocity of the hanging-wall layer
      xUhw = xUd*dsin(theta1-phib+beta-phif)
     1     /dsin(theta1+gamma1-2.d0*phib)
      xJ1 = xUd*dsin(gamma1-phib-beta+phif)
     1    /dsin(theta1+gamma1-2.d0*phib)
c
c---- the function Qu(qu,gamma1,theta1,theta,beta)
      call compext_34(qu,gamma1,theta1,theta,
     1     beta,phib,phif,cb,cf,rhof,dlambB,dlambF,
     1     Sg1e1,Sg1f1,Sgg1,Sgf,
     1     Shw,Sd,Sbs,xLg1e1,xLg1f1,xLgg1,xLgf,
     1     xz,gravHW,gravD,fopHW,fopGG1,fopGF,
     1     xUd,xJ,xUhw,xJ1,pmr,grav,fop)
c
c---- to get the MIN[Qu(qu,gamma1,theta1,theta,beta)]
          if(qu.lt.qulower) then
	    qulower = qu
        theta1lower = theta1
        gamma1lower = gamma1
         thetalower = theta
            xzlower = xz
              found = .true.
          end if
c
c---- end of loop for theta
 4000 continue
c---- end of loop for xz
 3000 continue
c
c---- determine which mechanism is better
      if(quD1.lt.qulower) then
c---- LANF is fully locked
        write(60,*) dtan(phif),cf,quD1
        write(61,*) dtan(phif),cf,1.d0
      else
c---- LANF is slipping or particially slipping
        write(60,*) dtan(phif),cf,qulower
        write(61,*) dtan(phif),cf,xzlower
c---- end of determination 
      endif
c
c---- end of loop for cf
 2000 continue
c---- end of loop for phif
 1000 continue
c
      end
c
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
      subroutine compext_34(qu,gamma1,theta1,theta,
     1        beta,phib,phif,cb,cf,rhof,dlambB,dlambF,
     1        Sg1e1,Sg1f1,Sgg1,Sgf,
     1        Shw,Sd,Sbs,xLg1e1,xLg1f1,xLgg1,xLgf,
     1        xz,gravHW,gravD,fopHW,fopGG1,fopGF,
     1        xUd,xJ,xUhw,xJ1,pmr,grav,fop)
c--------E-----F---F--------
c         \   /   /        |
c          \/   /          |
c            \ /           |
c--------------G------------
      implicit double precision (a-h,o-z)
c
c---- maximum resisting power
      pmr = cb*xUhw*dcos(phib)*xLg1e1 
     1    + cb*xJ1*dcos(phib)*xLg1f1
     1    + cb*xJ*dcos(phib)*xLgf
     1    + cf*xUd*dcos(phif)*xLgg1
c---- gravity power
      gravHW = Shw*xUhw*dsin(gamma1-phib)
      gravD = Sd*xUd*dsin(beta-phif)
      grav = (1.d0-rhof)*(gravHW+gravD)
c---- fluid overpressure power 
      fopHW = dlambB*Sg1e1*xUhw*dsin(phib)/dcos(gamma1)
     1      + dlambB*Sg1f1*xJ1*dsin(phib)/dcos(theta1)
c---- theta tends to be pi/2.d0 for infinite 1/dcos(theta)
      fopGF = dlambB*xLgf*xJ*dsin(phib)/2.d0
      fopGG1 = dlambF*Sgg1*xUd*dsin(phif)/dcos(beta)
      fop = fopHW + fopGG1 + fopGF
c---- external power
      pext = grav + fop
c---- upper bound force
      qu = pmr-pext
c
      return
      end 
c
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
      subroutine compext_35(quD1,phib,phif,
     1           cb,rhof,dlambB,dlambF)
c-------E-------F---------
c        \     /         |
c         \   /          |
c          \ /           |
c-----------G-------------
      implicit double precision (a-h,o-z)
c
c---- the values of theta and gamma
      theta = phib/2.d0+3.14159d0/4.d0
      gamma = phib/2.d0+3.14159d0/4.d0
c---- size of each segment
      xLge = 1.d0/dsin(gamma)
      xLgf = 1.d0/dsin(theta)
c---- the size of Shw
      Shw = 1.d0*xLge*dcos(gamma)/2.d0
     1      + 1.d0*xLgf*dcos(theta)/2.d0
c---- the size for fluid overpressure power
      Sge = 1.d0*xLge*dcos(gamma)/2.d0
      Sgf = 1.d0*xLgf*dcos(theta)/2.d0
c---- velocity of the hanging-wall layer
      xUhw = dsin(theta-phib)/dsin(theta+gamma-2.d0*phib)
      xJ = dsin(gamma-phib)/dsin(theta+gamma-2.d0*phib)
c
c---- maximum resisting power
      pmr = cb*xUhw*dcos(phib)*xLge + cb*xJ*dcos(phib)*xLgf
c---- gravity power
      grav = (1.d0-rhof)*Shw*xUhw*dsin(gamma-phib)
c---- fluid overpressure power 
      fop = dlambB*Sge*xUhw*dsin(phib)/dcos(gamma)
     1      + dlambB*Sgf*xJ*dsin(phib)/dcos(theta)
c---- external power
      pext = grav + fop
c---- upper bound force
      quD1 = pmr-pext
c
      return
      end
