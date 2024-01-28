program vasp_locpot

use constants,  only :  dp, LINELEN

use util,       only :  call_err,       &
                        assert_data,    &
                        time_formatting

use vasp,       only :  chgcar,         &
                        poscar,         &
                        chg_get,        &
                        chg_sum,        &
                        chg_write,      &
                        poscar_read,    &
                        poscar_write

implicit none

type(chgcar), dimension(2) :: DatChgIn
type(chgcar)               :: DatChgOut
type(poscar), dimension(2) :: DatPos

character(LEN=LINELEN) :: comm, line, lineGrid
integer                :: k, u, iFile
logical                :: isFile, isErr, isErrG

real(kind=dp)          :: eTime = 0d0

character(LEN=LINELEN), dimension(2) :: f_name
character(LEN=20),      dimension(2) :: tform
real(kind=dp),          dimension(2) :: M

k = iargc()
call getarg(0, comm)
if (k .ne. 4) then
   call call_err("Input error","need 4 args")
   call help(comm)
endif

call getarg(1, f_name(1))
inquire (file=f_name(1), exist=isFile)
if (.not. isFile) then
   call call_err("Input error","file 1 is missing")
   call help(comm)
endif

call getarg(2, f_name(2))
inquire (file=f_name(2), exist=isFile)
if (.not. isFile) then
   call call_err("Input error","file 2 is missing")
   call help(comm)
endif

call getarg(3, line)
read (line, *) M(1)
call getarg(4, line)
read (line, *) M(2)

do iFile = 1, 2
   open (unit=u, file=f_name(iFile), status='old')
   call poscar_read(u, DatPos(iFile))

   read (u, *)
   read (u, '(A)') lineGrid
   read (lineGrid, *) DatChgIn(iFile)%nGrid
   allocate (DatChgIn(iFile)%chg(DatChgIn(iFile)%nGrid(1), &
             DatChgIn(iFile)%nGrid(2), DatChgIn(iFile)%nGrid(3)))
   call chg_get(u, DatChgIn(iFile))
   eTime = eTime + DatChgIn(iFile)%eTime
   call time_formatting(DatChgIn(iFile)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"read vasp CHGCAR ",I1," - ",A,"( ",A," )")') &
        iFile, trim(adjustl(tform(1))), trim(adjustl(tform(2)))
   close(u)
enddo

call assert_data(DatPos(1)%LattVec, DatPos(2)%LattVec, 3, 3, isErr)
isErrG = isErr
call assert_data(DatPos(1)%nSpec, DatPos(2)%nSpec, isErr)
isErrG = isErrG .and. isErr
call assert_data(DatPos(1)%NameAtoms, DatPos(2)%NameAtoms, DatPos(1)%nSpec, isErr)
isErrG = isErrG .and. isErr
call assert_data(DatPos(1)%nAtomSpec, DatPos(2)%nAtomSpec, DatPos(1)%nSpec, isErr)
isErrG = isErrG .and. isErr
call assert_data(DatPos(1)%nAtom, DatPos(2)%nAtom, isErr)
isErrG = isErrG .and. isErr
call assert_data(DatPos(1)%coo, DatPos(2)%coo, 3, DatPos(1)%nAtom, isErr)
isErrG = isErrG .and. isErr
if (isErrG) call call_err('input error','two systems are different dim and coord')
if (isErrG) stop

call assert_data(DatChgIn(1)%nGrid, DatChgIn(2)%nGrid, 3, isErr)
if (isErr) call call_err('input error','two CHGCARs have different grid')
if (isErr) stop

DatChgOut%nGrid = DatChgIn(1)%nGrid
allocate (DatChgOut%chg(DatChgOut%nGrid(1), DatChgOut%nGrid(2), DatChgOut%nGrid(3)))

write (*,*) M
call chg_sum(DatChgIn(1), DatChgIn(2), DatChgOut, M)

eTime = eTime + DatChgOut%eTime
call time_formatting(DatChgOut%eTime, tform(1))
call time_formatting(eTime, tform(2))
write (*,'(2x,"two CHGCARs are merged - ",A," ( ",A," )")') &
        trim(adjustl(tform(1))), trim(adjustl(tform(2)))

f_name = "CHGCAR_sum"
open (unit=u, file=f_name)
call poscar_write(u, DatPos(1))
write (u,*)
call chg_write(u, DatChgOut)
eTime = eTime + DatChgOut%eTime
call time_formatting(DatChgOut%eTime, tform(1))
call time_formatting(eTime, tform(2))
write (*,'(2x,"charge density is written - ",A,"( ",A," )")') &
        trim(adjustl(tform(1))), trim(adjustl(tform(2))) 
close (u)

contains

!-----------------------------------------------------------------------------!

 subroutine help(comm_in)
   implicit none
   character(LEN=*), intent(in) :: comm_in

   write (*,'(3A)') " usage    : ",trim(adjustl(comm_in)), " [CHGCAR 1] [CHGCAR 2] [RATE 1] [RATE 2]"
   write (*,'(A)')  " default  : None, need 4 args"
   write (*,'(A)')  " [CHGCAR] : vasp CHGCAR file, no default"
   write (*,'(A)')  " [RATE]   : mixing rate of CHGCAR 1 and 2, no default"
   stop

 end subroutine help

end program vasp_locpot
