!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  this module does setup for galactic discs with live or analytic stars
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: datafiles, dim, dust, extern_spiral, externalforces, io,
!    options, part, physcon, prompting, random, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!--------------------------------------------------------------------------
!
! This subroutine is a utility for setting up isolated galactic discs
!
! The code has two main options. The code can setup a gas disc alone,
! which should then be exposed to some galactic potentials.
! Alternatively it will assume a live stellar component (with optional
! live bulge and dark matter) where stars are simply N-body particles.
!
! If setting a live disc then the code reads asciifiles of the
! galsetup python routines included with phantom.
! Seperate files are needded for each component (gas, stars, bulge, halo).
! This routine can read in ANY IC files providing they are in the format
!     x y z m vx vy vz
! so long as the different components are stored in the different asciifiles.
! Important note: the setup needs to know the number of each type of
! galactic component which requires the small galsetic.txt file.
!
! Created by Alex Pettitt from a similar routine by Clare Dobbs modified by
! Daniel Price for use in phantom.
!
!--------------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,     only:maxp,maxvxyzu,use_dust
 use setup_params, only:rhozero,npart_total
 use physcon, only:Rg,pi,solarm,pc,kpc,au,years,km,kboltz,mass_proton_cgs
 use units,   only:umass,udist,utime,set_units
 use ptmass,  only:h_acc,r_crit,rho_crit_cgs,icreate_sinks
 use mpiutils,       only:bcast_mpi
 use random,  only:ran2
 use eos,     only:gmw
 use extern_spiral, only:idisk,ibulg,ihalo,iarms,ibar,iread
 use part,    only:h2chemistry,abundance,iHI,dustfrac,istar,igas,ibulge,&
                   idarkmatter,iunknown,set_particle_type,ndusttypes
 use options, only:iexternalforce,icooling,nfulldump,use_dustfrac
 use externalforces, only:externalforce,initialise_externalforces
 use io,      only:master,fatal
 use prompting,      only:prompt
 use set_dust,       only:set_dustfrac
 use unifdis,        only:set_unifdis
 use stretchmap,     only:set_density_profile
 use timestep,       only:dtmax,tmax

 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real(kind=8),      intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: iseed,i,itot,it,ierr
 real :: rcyl,rcylin,faclod,xmin,xmax,ymin,ymax,zmin,zmax
 real :: rmax,rcyl2,rcylin2,totmass,totvol
 real :: xmax5,ymax5,zmax5,xi,yi,zi,r2
 real :: radius2,vx1,vy1,vz1,prob,thermal,disp,h2ratio,cs,Tinit,psep
 real :: fextxi,fextyi,fextzi,poti,dr,r,rand
 real(kind=8) :: uergg,angvel,vcirc,phi
 real :: dust_to_gas_ratio
 character(100)        :: partdist,vset,galsetupic
 integer, parameter, dimension(5) :: lenAr=(/1000,3000,2670,1830,4500/)
 real,    dimension(5) :: mratios=(/0.0008,0.0788,0.2293,0.1750,0.5161/)
 integer :: mratiosi(5)
 integer :: ierrf(5)
 real :: cdf1(lenAr(1)),rp1(lenAr(1))
 real :: cdf2(lenAr(2)),rp2(lenAr(2))
 real :: cdf3(lenAr(3)),rp3(lenAr(3))
 real :: cdf4(lenAr(4)),rp4(lenAr(4))
 real :: cdf5(lenAr(5)),rp5(lenAr(5))

 real:: totmassD,totmassG,totmassB,totmassH,totvolB,totvolH
 real:: rhozero1,rhozero2,rhozero3,rhozero4,h2,h3,h4
 real:: xis,yis,zis,mis,vxis,vyis,vzis,phaseis
 character(20) :: yn_gas,yn_star,yn_bulge,yn_halo
 character(30) :: sometext
 logical       :: use_live_stars


 time = 0.0
!
!--initialising units and flags
!
! set code units
 call set_units(dist=1.*kpc,mass=1.d10*solarm,G=1.)
 uergg = udist**2/utime**2
!
!--set input file options
!--maxvxyzu(3-4) and therefore ieos(1-2) are set in dim_galdisc
 iexternalforce = 8   !8=galactic disk potentials
 icooling       = 1   !1=cooling on, 0=off
 nfulldump      = 1
 idisk          = 0
 ibar           = 9 
 iarms          = 0 
 iread          = 0 
 tmax          = 200.
 dtmax         = 0.2

 hfact = 1.2 !-- cubic 
 hfact = 1.0 !-- quintic

!
!-------------------------Setting-energies-------------------------
!
 if (maxvxyzu >= 4) then
    !Non-Isothermal PV^gamma=K
    gamma = 1.6666667
 else
    !Isothermal PV = K as T=const.
    gamma = 1.0
 endif

 if (maxvxyzu >= 4) then

    !--H2 chemistry simulation
    !
    icooling       = 1   !1=cooling on, 0=off
    h2ratio = 0.
    gmw = (2.*h2ratio+(1.-2.*h2ratio)+0.4)/ &
          (0.1+h2ratio+(1.-2.*h2ratio))     ! should be 1.27 when h2ratio = 0
    Tinit = 1.0e4
    call prompt('Enter initial temperature in Kelvin',Tinit,0.0)
    cs = sqrt(kboltz*Tinit/(gmw*mass_proton_cgs))/(udist/utime) ! In code units
    polyk = cs*cs

    print*,' assuming mean molecular weight is ',gmw
    print*,'isothermal sound speed [km/s]  = ',cs*udist/utime/100000.

 else

    !-- isothermal simulation
    !thermal default at 10000K, 1K min and 1000000K max:
    icooling = 0   !1=cooling on, 0=off
    thermal = 15156.  
    call prompt('Enter initial temperature in Kelvin',thermal,0.1,1000000.)
    Tinit = thermal
    gmw = 1.26
    thermal = 3.0/2.0*thermal*Rg/gmw/uergg
    polyk = 2./3.*thermal
    print*,' assuming mean molecular weight is ',gmw
    print*,'isothermal sound speed [km/s]  = ', sqrt(polyk*(udist/utime)**2)/100000.


 endif
!
!--------------Setting-positions/velocities-for-gas-only-disc--------------
!
  rcylin= 0.
  rcyl  = 10.
  zmax  = 0.1

  call prompt('Enter inner radius for particle setup [kpc]',rcylin,0.)
  call prompt('Enter outer radius for particle setup [kpc]',rcyl,0.)
  call prompt('Enter disk height for particle setup  [kpc]',zmax,0.)

  rcylin = rcylin*kpc/udist
  rcyl   = rcyl*kpc/udist
  !--Height dist. comp to width of disk.

  zmax = zmax*kpc/udist
  zmin = -zmax
  xmax = rcyl
  ymax = rcyl
  xmin = -xmax
  ymin = -ymax
  rmax = sqrt(rcyl*rcyl + zmax*zmax)
  rcyl2 = rcyl*rcyl
  rcylin2 = rcylin*rcylin

  !--initialise random number generator
  iseed = -6485
  print "(a,i10)",' random seed = ',iseed
  print "(3(a,f10.3),a)",' galactic disc setup... rmin = ',rcylin,' rmax = ',rcyl,' in units of ',udist/kpc,' kpc'
  xi = ran2(iseed)
  npart = 1000000
  npart_total = 0
  call prompt('Enter number of particles ',npart,1,maxp)
  if (npart > maxp) call fatal('setup','npart > maxp')
  npartoftype(1) = npart

  print "(a,es10.3,a,1pg10.3,a)",'Mass is in units of ',umass,' g (',umass/solarm,' solar masses)'
  print "(a,es10.3,a,1pg10.3,a)",'Distance is in units of ',udist,' cm (',udist/kpc,' kpc)'
  print "(a,es10.3,a,1pg10.3,a)",'Time is in units of ',utime,' s (',utime/(1.0e9*years),' Gyr)'

  totvol = pi*(rcyl2 - rcylin2)
  rhozero = 31.830988618379067           ! 1.0e10 Msun / (pi* 10kpc**2)
  psep = 1.0
  call prompt('Enter initial surface density [Msun/pc^2] ',rhozero, 0.)
  call prompt('Enter particle separetion coeff',psep, 0.,10.)
  rhozero = rhozero*(solarm/(pc*pc))/(umass/(udist*udist))
  print "(a,es10.3)",'initial surface density = ',rhozero

  totmass = rhozero*totvol
  massoftype(igas) = totmass/real(npart) 

  rhozero = totmass/totvol/(zmax-zmin)
  psep = psep*hfact*(massoftype(1)/rhozero)**(1./3.)
  print "(a,es10.3)",'initial mean density  = ',rhozero
  print "(a,es10.3)",'particle separation  = ',psep
  print "(a,es10.3)",'particle mass = ',totmass/real(npart)*umass/solarm

  !
  !------------------------Setting-positions------------------------
  !


  npart = 0
  call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                hfact,npart,xyzh,.False.,rcylmin=rcylin,rcylmax=rcyl,nptot=npart_total)

  !shift z to center
  xyzh(3,:) = xyzh(3,:) - (maxval(xyzh(3,:)) + minval(xyzh(3,:)))/2.

!  call set_density_profile(npart,xyzh,min=0.0,max=rcyl,rhofunc=rhofuncr,geom=2,coord=1)
!  call set_density_profile(npart,xyzh,min=-0.4,max=0.4,rhofunc=rhofuncz,geom=1,coord=3)

  npartoftype(:)    = 0
  npartoftype(igas) = npart
  do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
  enddo

  do i=1,npart

     if (xyzh(1,i)>=1.9  .and. xyzh(1,i)<=2.1 .and. &
         xyzh(2,i)>=-0.1 .and. xyzh(2,i)<=0.1 .and. &
         xyzh(3,i)>=-0.05 .and. xyzh(3,i)<=0.05) then

         print*,'select particle', i
         print*,xyzh(1,i),xyzh(2,i),xyzh(3,i)

     endif

  enddo

  print*,'Final number of particles set = ',npart
  !
  !------------------------Setting-velocities------------------------
  !

  call initialise_externalforces(iexternalforce,ierr)
  do i=1,npart
     fextxi=0.
     fextyi=0.
     fextzi=0.
     poti  =0.
     radius2 = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
     !--Pull an initial velocity from the actual rotation curve defined by .in file
     call externalforce(iexternalforce,xyzh(1,i),xyzh(2,i),0.0,xyzh(4,i), &
                        0.0,fextxi,fextyi,fextzi,poti)
     
     angvel= sqrt( sqrt(fextxi**2+fextyi**2) *sqrt(radius2) )

     vxyzu(1,i) = -angvel*xyzh(2,i)/sqrt(radius2)
     vxyzu(2,i) = +angvel*xyzh(1,i)/sqrt(radius2)
     vxyzu(3,i) = 0.

     !--Add Gaussian random velocity dispersion
     !--For 5 km/s dispersion, use Gaussian with sigma=disp=5
     !--Set velocity dispersion parameter:
     disp=5.

!     rand = huge(rand)
!     prob = 0.
!     do while (prob < rand)
!        vx1=40.*(ran2(iseed)-0.5)
!        prob=exp(-(vx1/disp)**2./2.)
!        rand = ran2(iseed)
!     enddo
!     vxyzu(1,i)=vxyzu(1,i)+vx1*100000.*utime/udist
!
!     rand = huge(rand)
!     prob = 0.
!     do while (prob < rand)
!        vy1=40.*(ran2(iseed)-0.5)
!        prob=exp(-(vy1/disp)**2./2.)
!        rand = ran2(iseed)
!     enddo
!     vxyzu(2,i)=vxyzu(2,i)+vy1*100000.*utime/udist
!
!     rand = huge(rand)
!     prob = 0.
!     do while (prob < rand)
!        vz1=40.*(ran2(iseed)-0.5)
!        prob=exp(-(vz1/disp)**2./2.)
!        rand = ran2(iseed)
!    enddo
!     vxyzu(3,i)=vxyzu(3,i)+vz1*100000.*utime/udist
!

    if (maxvxyzu  >= 4) then

        vxyzu(4,i) = (1./(gamma-1.)*kboltz*Tinit/(gmw*mass_proton_cgs))&
                    /(udist**2/utime**2)     ! code unit
    endif

 !
 !------------------------Setting-abundances------------------------
 !


    if (h2chemistry) then

       abundance(:,i)   = 0.
       abundance(iHI,i) = 1.  ! assume all atomic hydrogen initially

    endif

 !
 !------------------------Dust-fractions----------------------------
 !
     if (use_dustfrac) then
     
         call set_dustfrac(dust_to_gas_ratio,dustfrac(:,i))

     endif


  enddo

! for sink particles 
 h_acc = 0.1*pc/udist
 r_crit= 2.*h_acc
 icreate_sinks = 1
 rho_crit_cgs  = 100.


contains

real function rhofuncr(r)
 real, intent(in) :: r

 rhofuncr = 71.76*exp(-r/4.8)

end function rhofuncr

real function rhofuncz(z)
 real, intent(in) :: z

 rhofuncz = 64000.3/(cosh(z/0.05))**2

end function rhofuncz

end subroutine setpart

end module setup
