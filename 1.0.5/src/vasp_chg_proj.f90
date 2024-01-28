program vasp_chg_proj

use constants,  only :  dp, LINELEN

use util,       only :  cross_product,  &
                        time_formatting

use iotool,     only :  io_printerr

use vasp,       only :  chgcar,         &
                        poscar,         &
                        chg_get,        &
                        poscar_read

implicit none

type(chgcar) :: DatChg
type(poscar) :: DatPos

character(LEN=LINELEN) :: comm, f_name_in, f_name_out
character(LEN=20)      :: tform

character(LEN=1), dimension(3) :: Axis

integer :: u = 10
integer :: ii, ix, iy, iz, iGrid1, iGrid2, iAxis, iAxis1, iAxis2

logical :: is_1d
logical :: is_2d

real(kind=dp) :: Loc_1

real(kind=dp), dimension(3)                    :: mult, LenVec
real(kind=dp), dimension(:,:,:),   allocatable :: dat_1d
real(kind=dp), dimension(:,:,:,:), allocatable :: dat_2d

call init()

open (unit=u, file=f_name_in, status='old')

call poscar_read(u, DatPos)
call chg_get(u, DatChg)

close(u)

call time_formatting(DatChg%eTime, tform)

write (*,'(1x," read vasp LOCPOT successfully - ",A)') trim(adjustl(tform))

if (is_1d) then
   allocate (dat_1d(3,maxval(DatChg%nGrid(:)),3))
   dat_1d(1,:,:) = 0d0
   dat_1d(2,:,:) = 9d4
   dat_1d(3,:,:) = -9d4
endif
if (is_2d) then
   allocate (dat_2d(3,maxval(DatChg%nGrid(:)),maxval(DatChg%nGrid(:)),3)) 
   dat_2d(1,:,:,:) = 0d0
   dat_2d(2,:,:,:) = 9d4
   dat_2d(3,:,:,:) = -9d4
endif

do iz = 1, DatChg%nGrid(3)
do iy = 1, DatChg%nGrid(2)
do ix = 1, DatChg%nGrid(1)
   if (is_1d) then
      dat_1d(1,ix,1) = dat_1d(1,ix,1) + DatChg%chg(ix,iy,iz)
      dat_1d(1,iy,2) = dat_1d(1,iy,2) + DatChg%chg(ix,iy,iz)
      dat_1d(1,iz,3) = dat_1d(1,iz,3) + DatChg%chg(ix,iy,iz)
      if (dat_1d(2,ix,1) .gt. DatChg%chg(ix,iy,iz)) dat_1d(2,ix,1) = DatChg%chg(ix,iy,iz)
      if (dat_1d(2,iy,2) .gt. DatChg%chg(ix,iy,iz)) dat_1d(2,iy,2) = DatChg%chg(ix,iy,iz)
      if (dat_1d(2,iz,3) .gt. DatChg%chg(ix,iy,iz)) dat_1d(2,iz,3) = DatChg%chg(ix,iy,iz)
      if (dat_1d(3,ix,1) .lt. DatChg%chg(ix,iy,iz)) dat_1d(3,ix,1) = DatChg%chg(ix,iy,iz)
      if (dat_1d(3,iy,2) .lt. DatChg%chg(ix,iy,iz)) dat_1d(3,iy,2) = DatChg%chg(ix,iy,iz)
      if (dat_1d(3,iz,3) .lt. DatChg%chg(ix,iy,iz)) dat_1d(3,iz,3) = DatChg%chg(ix,iy,iz)
   endif
   if (is_2d) then
      dat_2d(1,ix,iy,3) = dat_2d(1,ix,iy,3) + DatChg%chg(ix,iy,iz)
      dat_2d(1,iy,iz,1) = dat_2d(1,iy,iz,1) + DatChg%chg(ix,iy,iz)
      dat_2d(1,iz,ix,2) = dat_2d(1,iz,ix,2) + DatChg%chg(ix,iy,iz)
      if (dat_2d(2,ix,iy,3) .gt. DatChg%chg(ix,iy,iz)) dat_2d(2,ix,iy,3) = DatChg%chg(ix,iy,iz)
      if (dat_2d(2,iy,iz,1) .gt. DatChg%chg(ix,iy,iz)) dat_2d(2,iy,iz,1) = DatChg%chg(ix,iy,iz)
      if (dat_2d(2,iz,ix,2) .gt. DatChg%chg(ix,iy,iz)) dat_2d(2,iz,ix,2) = DatChg%chg(ix,iy,iz)
      if (dat_2d(3,ix,iy,3) .lt. DatChg%chg(ix,iy,iz)) dat_2d(3,ix,iy,3) = DatChg%chg(ix,iy,iz)
      if (dat_2d(3,iy,iz,1) .lt. DatChg%chg(ix,iy,iz)) dat_2d(3,iy,iz,1) = DatChg%chg(ix,iy,iz)
      if (dat_2d(3,iz,ix,2) .lt. DatChg%chg(ix,iy,iz)) dat_2d(3,iz,ix,2) = DatChg%chg(ix,iy,iz)
   endif
enddo
enddo
enddo

do iAxis = 1, 3
   LenVec(iAxis) = sqrt(dot_product(DatPos%LattVec(:,iAxis),DatPos%LattVec(:,iAxis)))
enddo

if (is_1d) then
   do iAxis = 1, 3
      mult(iAxis) = dble(DatChg%nGrid(iAxis))/dble(DatChg%nGrid(1)*DatChg%nGrid(2)*DatChg%nGrid(3))
      f_name_out = "proj_1d_"//Axis(iAxis)//".dat"
      open (unit=u, file=f_name_out, action='write')
      write (u,'("#",A9,4A16)') 'Loc.','Avgerage','Min','Max'
      do iGrid1 = 1, DatChg%nGrid(iAxis)
         write (u,'(F10.5,3E16.8E2)') &
                LenVec(iAxis)*float(iGrid1-1)/float(DatChg%nGrid(iAxis)), &
                dat_1d(1,iGrid1,iAxis)*mult(iAxis), &
                dat_1d(2,iGrid1,iAxis), &
                dat_1d(3,iGrid1,iAxis)
      enddo
      close (u)
   enddo
endif

if (is_2d) then
   do iAxis = 1, 3
      iAxis1 = 1 + mod(iAxis,3)
      iAxis2 = 1 + mod(iAxis+1,3)
      f_name_out = "proj_2d_"//Axis(iAxis1)//Axis(iAxis2)//".dat"
      open (unit=u, file=f_name_out, action='write')
      write (u,'("#",A9,A10,4A16)') 'Loc.1','Loc.2','Sum','Avgerage','Min','Max'
      do iGrid1 = 1, DatChg%nGrid(iAxis1)
         Loc_1 = LenVec(iAxis1)*float(iGrid1-1)/float(DatChg%nGrid(iAxis1))
         do iGrid2 = 1, DatChg%nGrid(iAxis2)
            write (u,'(2F10.5,4E16.8E2)') &
               Loc_1, &
               LenVec(iAxis2)*float(iGrid2-1)/float(DatChg%nGrid(iAxis2)), &
               dat_2d(1,iGrid1,iGrid2,iAxis), & 
               dat_2d(1,iGrid1,iGrid2,iAxis)/dble(DatChg%nGrid(iAxis)), & 
               dat_2d(2,iGrid1,iGrid2,iAxis), & 
               dat_2d(3,iGrid1,iGrid2,iAxis)
         enddo
      enddo
      close (u)
   enddo
endif

contains

!-----------------------------------------------------------------------------!

 subroutine init()
   implicit none

   integer :: ik, k
   logical :: getFile, getOpt, isFile
   character(LEN=LINELEN) :: line

   k         = iargc()
   f_name_in = 'None'

   call getarg(0, comm)

   is_1d  = .true.
   is_2d  = .false.
   getOpt = .false.

   if (k .ge. 1) then
      LK : do ik = 1, k
         call getarg(ik, line)
         select case (trim(line))
         case ('-h')
            call help()
         case ('-1d')
            if (getOpt) cycle LK
            is_1d  = .true.
            is_2d  = .false.
            getOpt = .true.
            cycle LK
         case ('-2d')
            if (getOpt) cycle LK
            is_1d  = .false.
            is_2d  = .true.
            getOpt = .true.
            cycle LK
         case ('-a')
            if (getOpt) cycle LK
            is_1d  = .true.
            is_2d  = .true.
            getOpt = .true.
            cycle LK
         end select

         if (getFile) cycle LK
         inquire (file=trim(line), exist=isFile)
         if (isFile) then
            getFile   = .true.
            f_name_in = trim(line)
         endif
      enddo LK
   endif

   if (trim(f_name_in) .eq. 'None') call help()
   Axis(1) = 'x'
   Axis(2) = 'y'
   Axis(3) = 'z'

 end subroutine init

!-----------------------------------------------------------------------------!

 subroutine help()
   implicit none

   write (*,'(A)') " usage    : "//trim(comm)//" [CHGCAR]"
   write (*,'(A)') " default  : [OPTION] = None"
   write (*,'(A)') "            [CHGCAR] = None" 
   write (*,*)
   write (*,'(A)') " [OPTION] : -h  : call this help"
   write (*,'(A)') "          : -1d : 1-dim projection only, default"
   write (*,'(A)') "                  (OUTPUTs : x, y, z)"
   write (*,'(A)') "          : -2d : 2-dim projection only"
   write (*,'(A)') "                  (OUTPUTs : xy, yz, zx)"
   write (*,'(A)') "          : -a  : both 1-dim and 2-dim projctions"
   write (*,'(A)') " [CHGCAR] : vasp CHGCAR file"
   stop

 end subroutine help

end program vasp_chg_proj
