!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine to calculate A.M. loss and torques from different part
!
!  REFERENCES: None
!
!  OWNER: Zhi Li
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, io, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'A.M. loss'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use eos,     only:get_spsound
 use part,    only:fxyzu,fext,rhoh
 use externalforces, only:externalforce
 use forces,  only:fxyzu_press,fxyzu_visco

 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=20) :: filename
 integer :: i,num
 integer :: partid(9)
 integer, parameter :: iparams = 10
 integer, parameter :: iecc    = 23
 real :: fg(3),vel(3),pos(3),rad,zad,Li(3),Qg(3),Qp(3),Qv(3),pot

 !Print the analysis being done

 if (time == 100.0) then

   open(unit=24,file='TorquesPV',form='unformatted',access='stream',status='replace')

   do i=1,npart

    if (xyzh(1,i) >= -1. .and. xyzh(1,i) <= 1. .and. &
        xyzh(2,i) >= -1. .and. xyzh(2,i) <= 1. ) then

    call externalforce(8,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i), &
                          time,fg(1),fg(2),fg(3),pot)    

    write(24),xyzh(1,i),xyzh(2,i),xyzh(3,i), &
              pmass*(xyzh(1,i)*         fg(2)  -xyzh(2,i)*           fg(1)), &
              pmass*(xyzh(1,i)*fxyzu_press(2,i)-xyzh(2,i)*fxyzu_press(1,i)), &
              pmass*(xyzh(1,i)*fxyzu_visco(2,i)-xyzh(2,i)*fxyzu_visco(1,i)) 

    endif


   enddo

   close(24)

 endif

 !------------------------------------------------------------------------!

 !CHOOSE YOUR PARTICLE

 partid = (/1070274,1070275,1070276,1070542,1070543,1070544,1070809,1070810,1070811/)


 do num = 1, size(partid)

   i = partid(num)

   write(filename,"(a,i7.7)")"tracing_",i

   if (i > npart) print*,'Particle chosen does not exist.'

   pos = xyzh(1:3,i)
   vel = vxyz(1:3,i)
   rad = sqrt(pos(1)**2+pos(2)**2)
   zad = pos(3)

   call externalforce(8,xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i), &
                          time,fg(1),fg(2),fg(3),pot)

   Li(1) = pmass*(xyzh(2,i)*vxyz(3,i)-xyzh(3,i)*vxyz(2,i))
   Li(2) = pmass*(xyzh(3,i)*vxyz(1,i)-xyzh(1,i)*vxyz(3,i))
   Li(3) = pmass*(xyzh(1,i)*vxyz(2,i)-xyzh(2,i)*vxyz(1,i))

  ! Qg(3) = pmass*(xyzh(1,i)* fext(2,i)-xyzh(2,i)* fext(1,i))
  ! Qp(3) = pmass*(xyzh(1,i)*fxyzu(2,i)-xyzh(2,i)*fxyzu(1,i))
  ! Qv(3) = 0.0
   Qg(3) = pmass*(xyzh(1,i)*         fg(2)  -xyzh(2,i)*           fg(1))
   Qp(3) = pmass*(xyzh(1,i)*fxyzu_press(2,i)-xyzh(2,i)*fxyzu_press(1,i))
   Qv(3) = pmass*(xyzh(1,i)*fxyzu_visco(2,i)-xyzh(2,i)*fxyzu_visco(1,i))

   if (time==0.0) then
      open(unit=iecc,file=filename,status="unknown")
      write(iecc,'("# A.M. loss and torques for particle ",i10.10)') i
      write(iecc,"('#',9(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time', &
           2,'x', &
           3,'y', &
           4,'z', &
           5,'rho', &
           6,'Lz', &
           7,'Torque-g', &
           8,'Torque-p', &
           9,'Torque-v'
   else
      open(unit=iecc,file=filename,status="old",position="append")
   endif

   write(iecc,'(9(es18.10,1X))') & 
         time,pos(1),pos(2),pos(3),rhoh(xyzh(4,i),pmass),Li(3),Qg(3),Qp(3),Qv(3)

   close(unit=iecc)

 enddo

end subroutine do_analysis

end module analysis

