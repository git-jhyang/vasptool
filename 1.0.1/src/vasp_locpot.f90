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
integer                :: ix, iy, iz
logical                :: isFile
real(kind=dp)          :: mult, LenVec

real(kind=dp), dimension(:,:), allocatable :: pot

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
allocate (pot(maxval(DatPot%nGrid(:)),3))

pot(:,:) = 0d0

call chg_get(u, DatPot)
call time_formatting(DatPot%eTime, tform)

write (*,'(2A)') " read vasp LOCPOT successfully - ",trim(adjustl(tform))

close(u)

do iz = 1, DatPot%nGrid(3)
do iy = 1, DatPot%nGrid(2)
do ix = 1, DatPot%nGrid(1)
   pot(ix,1) = pot(ix,1) + DatPot%chg(ix,iy,iz)
   pot(iy,2) = pot(iy,2) + DatPot%chg(ix,iy,iz)
   pot(iz,3) = pot(iz,3) + DatPot%chg(ix,iy,iz)
enddo
enddo
enddo

do ix = 1, 3
   LenVec = sqrt(dot_product(DatPos%LattVec(:,ix),DatPos%LattVec(:,ix)))/float(DatPot%nGrid(ix))
   mult   = dble(DatPot%nGrid(ix))/dble(DatPot%nGrid(1)*DatPot%nGrid(2)*DatPot%nGrid(3))
   write (f_name,'("prj_LOCPOT_axis_",I1,".dat")') ix
   open (unit=u, file=f_name)
   write (u, '("#",A9,A20)') 'Location','Potential'
   do iy = 1, DatPot%nGrid(ix)
      write (u, '(F10.5, ES20.12E2)') LenVec*float(iy-1), pot(iy,ix)*mult
   enddo
   close(u)
enddo

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
