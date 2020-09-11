!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  This will reset the particles sound speed for isothermal runs
!
!  REFERENCES: None
!
!  OWNER: Zhi Li
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, centreofmass, dim, part
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use eos,          only: polyk,gmw
 use units,        only: umass,udist,utime
 use physcon,      only: Rg,pi,solarm,pc,kpc,au,years,km,kboltz,mass_proton_cgs
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: thermal
 integer                :: i
 !


 print*,'previous iso. sound speed [km/s]  = ', sqrt(polyk*(udist/utime)**2)/100000.
 thermal = 60623.
 thermal = 3.0/2.0*thermal*Rg/gmw/((udist/utime)**2)
 polyk = 2./3.*thermal
 print*,'updated iso. sound speed [km/s]  = ', sqrt(polyk*(udist/utime)**2)/100000.

!  do i = 1,npart
!       if (xyzh(1,i) < xmin)         xyzh(1,i) = xyzh(1,i) + dxbound
!       if (xyzh(1,i) > xmin+dxbound) xyzh(1,i) = xyzh(1,i) - dxbound
!       if (xyzh(2,i) < ymin)         xyzh(2,i) = xyzh(2,i) + dybound
!       if (xyzh(2,i) > ymin+dybound) xyzh(2,i) = xyzh(2,i) - dybound
!       if (xyzh(3,i) < zmin)         xyzh(3,i) = xyzh(3,i) + dzbound
!       if (xyzh(3,i) > zmin+dzbound) xyzh(3,i) = xyzh(3,i) - dzbound
!  enddo
!

 return
end subroutine modify_dump

end module moddump

