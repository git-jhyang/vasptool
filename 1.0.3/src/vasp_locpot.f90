program vasp_locpot

use constants,  only :  dp, LINELEN

use util,       only :  cross_product,  &
                        call_err,       &
                        time_formatting

use vasp,       only :  chgcar,         &
                        poscar,         &
                        chg_get,        &
                        poscar_read

implicit none

type(chgcar) :: DatPot
type(poscar) :: DatPos

character(LEN=LINELEN) :: comm, f_name, line, lineGrid
character(LEN=20)      :: tform
integer                :: k, u = 10
integer                :: ii, ix, iy, iz, iGrid, iAxis
logical                :: isFile

real(kind=dp), dimension(3)                  :: mult, LenVec
real(kind=dp), dimension(:,:,:), allocatable :: pot

k = iargc()
call getarg(0, comm)
if (k .eq. 0) then
   f_name = 'LOCPOT'
else
   call getarg(1, f_name)
endif

inquire (file=f_name, exist=isFile)
if (.not. isFile) then
   call call_err("Input error","file is missing")
   call help(comm)
endif

open (unit=u, file=f_name, status='old')
call poscar_read(u, DatPos)

read (u, *)
read (u, '(A)') lineGrid
read (lineGrid, *) DatPot%nGrid

allocate (DatPot%chg(DatPot%nGrid(1), DatPot%nGrid(2), DatPot%nGrid(3)))
allocate (pot(3,3,maxval(DatPot%nGrid(:))))

pot(1,:,:) = 0d0
pot(2,:,:) = 9d2
pot(3,:,:) = -9d2

call chg_get(u, DatPot)
call time_formatting(DatPot%eTime, tform)

write (*,'(2A)') " read vasp LOCPOT successfully - ",trim(adjustl(tform))

close(u)

do iz = 1, DatPot%nGrid(3)
do iy = 1, DatPot%nGrid(2)
do ix = 1, DatPot%nGrid(1)
   pot(1,1,ix) = pot(1,1,ix) + DatPot%chg(ix,iy,iz)
   pot(1,2,iy) = pot(1,2,iy) + DatPot%chg(ix,iy,iz)
   pot(1,3,iz) = pot(1,3,iz) + DatPot%chg(ix,iy,iz)
   if (pot(2,1,ix) .gt. DatPot%chg(ix,iy,iz)) pot(2,1,ix) = DatPot%chg(ix,iy,iz)
   if (pot(2,2,iy) .gt. DatPot%chg(ix,iy,iz)) pot(2,2,iy) = DatPot%chg(ix,iy,iz)
   if (pot(2,3,iz) .gt. DatPot%chg(ix,iy,iz)) pot(2,3,iz) = DatPot%chg(ix,iy,iz)
   if (pot(3,1,ix) .lt. DatPot%chg(ix,iy,iz)) pot(3,1,ix) = DatPot%chg(ix,iy,iz)
   if (pot(3,2,iy) .lt. DatPot%chg(ix,iy,iz)) pot(3,2,iy) = DatPot%chg(ix,iy,iz)
   if (pot(3,3,iz) .lt. DatPot%chg(ix,iy,iz)) pot(3,3,iz) = DatPot%chg(ix,iy,iz)
enddo
enddo
enddo

do iAxis = 1, 3
   LenVec(iAxis) = sqrt(dot_product(DatPos%LattVec(:,iAxis),DatPos%LattVec(:,iAxis)))
   mult(iAxis)   = dble(DatPot%nGrid(iAxis))/dble(DatPot%nGrid(1)*DatPot%nGrid(2)*DatPot%nGrid(3))
enddo

write (f_name,'(A)') "LOCPOT_proj.dat"
open (unit=u, file=f_name)
write (u,'("#",A9,2(48x,A10))') 'Axis 1','Axis 2','Axis 3'
write (u,'("#",A9,3A16,2(A10,3A16))') 'Loc','Avg','Min','Max','Loc','Avg','Min','Max','Loc','Avg','Min','Max'
do iGrid = 1, maxval(DatPot%nGrid)
   do iAxis = 1, 3
      if (iGrid .gt. DatPot%nGrid(iAxis)) then
         write (u,'(58x)',advance='no')
      else
         write (u,'(F10.5)',advance='no') LenVec(iAxis)*float(iGrid-1)/float(DatPot%nGrid(iAxis))
         write (u,'(3ES16.8E2)',advance='no') pot(1,iAxis,iGrid)*mult(iAxis), pot(2:3,iAxis,iGrid)
      endif
   enddo
   write (u,*)
enddo
close (u)

contains

!-----------------------------------------------------------------------------!

 subroutine help(comm_in)
   implicit none
   character(LEN=*), intent(in) :: comm_in

   write (*,'(3A)') " usage    : ", trim(adjustl(comm_in)), " [LOCPOT]"
   write (*,'(A)')  " default  : ./LOCPOT " 
   write (*,'(A)')  " [LOCPOT] : vasp LOCPOT file"
   stop

 end subroutine help

end program vasp_locpot
