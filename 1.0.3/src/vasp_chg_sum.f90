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

character(LEN=LINELEN) :: comm, line, lineGrid, f_name_out
integer                :: k, u, iFile
logical                :: isFile, isErr, isErrG

real(kind=dp)          :: eTime = 0d0, deTime

character(LEN=LINELEN), dimension(2) :: f_name_in
character(LEN=20),      dimension(2) :: tform
real(kind=dp),          dimension(2) :: M

k = iargc()
M(:) = 1.0d0

call getarg(0, comm)
if (k .lt. 2) then
   call call_err("Input error","need more than 2 args")
   call help(comm)
else
   do iFile = 1, 2
      call getarg(iFIle, f_name_in(iFile))
      inquire (file=f_name_in(iFile), exist=isFile)
      if (.not. isFile) then
         call call_err("Input error - file is missing",trim(f_name_in(iFile)))
         call help(comm)
      endif
   enddo
endif

if (k .ge. 3) then
   call getarg(3, line)
   read (line, *) M(1)
endif
if (k .ge. 4) then
   call getarg(4, line)
   read (line, *) M(2)
endif

do iFile = 1, 2
   open (unit=u, file=f_name_in(iFile), status='old')
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
isErrG = isErrG .or. isErr
call assert_data(DatPos(1)%NameAtoms, DatPos(2)%NameAtoms, DatPos(1)%nSpec, isErr)
isErrG = isErrG .or. isErr
call assert_data(DatPos(1)%nAtomSpec, DatPos(2)%nAtomSpec, DatPos(1)%nSpec, isErr)
isErrG = isErrG .or. isErr
call assert_data(DatPos(1)%nAtom, DatPos(2)%nAtom, isErr)
isErrG = isErrG .or. isErr
call assert_data(DatPos(1)%coo, DatPos(2)%coo, 3, DatPos(1)%nAtom, isErr)
isErrG = isErrG .or. isErr
if (isErrG) write (*,'(2x," !!! WARNING : two CHGCARs have different crystal information !!!")')
call assert_data(DatChgIn(1)%isLong, DatChgIn(2)%isLong, isErr)
if (isErr) write (*,'(2x," !!! WARNING : two CHGCARs have different format !!!")')
call assert_data(DatChgIn(1)%nGrid, DatChgIn(2)%nGrid, 3, isErr)
if (isErr) call call_err('input error','two CHGCARs have different grid')
if (isErr) stop

DatChgOut%nGrid = DatChgIn(1)%nGrid
allocate (DatChgOut%chg(DatChgOut%nGrid(1), DatChgOut%nGrid(2), DatChgOut%nGrid(3)))

call chg_sum(DatChgIn(1), DatChgIn(2), DatChgOut, M)

eTime = eTime + DatChgOut%eTime
call time_formatting(DatChgOut%eTime, tform(1))
call time_formatting(eTime, tform(2))
write (*,'(2x,"two CHGCARs are merged ( ratio : ",F5.2," and ",F5.2," ) - ",A," ( ",A," )")') &
        M, trim(adjustl(tform(1))), trim(adjustl(tform(2)))

call assert_data(DatChgIn(1)%isLong, DatChgIn(2)%isLong, isErr)
if (      isErr) DatChgOut%isLong = .false.
if (.not. isErr) DatChgOut%isLong = DatChgIn(1)%isLong

if (      DatChgOut%isLong) f_name_out = 'CHGCAR_sum'
if (.not. DatChgOut%isLong) f_name_out = 'CHG_sum'

open (unit=u, file=f_name_out)
call poscar_write(u, DatPos(1), .false.)
write (u,*)
call chg_write(u, DatChgOut)
close (u)

deTime = DatChgOut%eTime
eTime  = eTime + deTime
call time_formatting(deTime, tform(1))
call time_formatting(eTime, tform(2))
write (*,'(2x,"charge density (",A,") is written - ",A,"( ",A," )")') &
        trim(f_name_out), trim(adjustl(tform(1))), trim(adjustl(tform(2))) 

contains

!-----------------------------------------------------------------------------!

 subroutine help(comm_in)
   implicit none
   character(LEN=*), intent(in) :: comm_in

   write (*,'(3A)') " usage    : ",trim(adjustl(comm_in)), " [CHGCAR 1] [CHGCAR 2] [RATE 1] [RATE 2]"
   write (*,'(A)')  " default  : [CHGCAR] = None"
   write (*,'(A)')  "            [RATE]   = 1.0"
   write (*,*)
   write (*,'(A)')  " [CHGCAR] : vasp CHGCAR file, no default"
   write (*,'(A)')  " [RATE]   : mixing rate of CHGCAR 1 and 2, default = 1.0"
   write (*,*)
   stop

 end subroutine help

end program vasp_locpot
