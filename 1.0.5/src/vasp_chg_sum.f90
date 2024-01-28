program vasp_chg_sum

use constants,  only :  dp, LINELEN,    &
                        CHKOK, CHKERR

use util,       only :  assert_data,    &
                        time_formatting

use iotool,     only :  io_printerr,    &
                        io_linesplit

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

character(LEN=LINELEN) :: comm, f_name_out
integer                :: u, iFile
logical                :: isErr, isErrG

real(kind=dp)          :: eTime = 0d0, deTime

character(LEN=LINELEN), dimension(2) :: f_name_in
character(LEN=20),      dimension(2) :: tform
real(kind=dp),          dimension(2) :: M

call init()

do iFile = 1, 2
   open (unit=u, file=f_name_in(iFile), status='old')
   call poscar_read(u, DatPos(iFile))
   call chg_get(u, DatChgIn(iFile))
   eTime = eTime + DatChgIn(iFile)%eTime
   call time_formatting(DatChgIn(iFile)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"read vasp CHGCAR ",I1," - ",A,"( ",A," )")') &
        iFile, trim(adjustl(tform(1))), trim(adjustl(tform(2)))
   close(u)
enddo

isErrG = assert_data(DatPos(1)%LattVec, DatPos(2)%LattVec, 3, 3)
isErrG = isErrG .or. assert_data(DatPos(1)%nSpec, DatPos(2)%nSpec)
isErrG = isErrG .or. assert_data(DatPos(1)%NameAtoms, DatPos(2)%NameAtoms, DatPos(1)%nSpec)
isErrG = isErrG .or. assert_data(DatPos(1)%nAtomSpec, DatPos(2)%nAtomSpec, DatPos(1)%nSpec)
isErrG = isErrG .or. assert_data(DatPos(1)%nAtom, DatPos(2)%nAtom)
isErrG = isErrG .or. assert_data(DatPos(1)%coo, DatPos(2)%coo, 3, DatPos(1)%nAtom)
if (isErrG) call io_printerr(trim(comm)//' - crystal information','two CHGCARs have different crystal information')
isErr = assert_data(DatChgIn(1)%isLong, DatChgIn(2)%isLong)
if (isErr) call io_printerr(trim(comm)//' - chgcar format','two CHGCARs have different format')
isErr = assert_data(DatChgIn(1)%nGrid, DatChgIn(2)%nGrid, 3)
if (isErr) call io_printerr(trim(comm)//' - chgcar grid','two CHGCARs have different grid')
if (isErr) stop

call chg_sum(DatChgIn(1), DatChgIn(2), DatChgOut, M)

eTime = eTime + DatChgOut%eTime
call time_formatting(DatChgOut%eTime, tform(1))
call time_formatting(eTime, tform(2))
write (*,'(2x,"two CHGCARs are merged ( ratio : ",F5.2," and ",F5.2," ) - ",A," ( ",A," )")') &
        M, trim(adjustl(tform(1))), trim(adjustl(tform(2)))

if (      DatChgOut%isLong) f_name_out = 'CHGCAR_sum'
if (.not. DatChgOut%isLong) f_name_out = 'CHG_sum'

open (unit=u, file=f_name_out)

DatPos(1)%isLong = .false.
call poscar_write(u, DatPos(1))
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

 subroutine init()
   implicit none

   integer :: ik, k, CHK
   logical :: getFile1, getFile2, getRate1, getRate2, isFile

   character(LEN=LINELEN) :: line

   k        = iargc()
   M(:)     = 1.0d0
   getFile1 = .false.
   getFile2 = .false.
   getRate1 = .false.
   getRate2 = .false.

   call getarg(0, comm)

   if (k .lt. 2) then
      call io_printerr(trim(comm)//" - Input error","need more than 2 args")
      call help()
   else
      LK : do ik = 1, k
         call getarg(ik, line)
         select case (trim(line))
         case ('-h')
            call help()
         end select
         if (.not. getFile1) then
            inquire (file=trim(line), exist=isFile)
            if (.not. isFile) then
               call io_printerr(trim(comm)//" - Input error","file 1 is missing")
               call help()
            endif
            getFile1     = .true.
            f_name_in(1) = trim(line)
         elseif (.not. getFile2) then
            inquire (file=trim(line), exist=isFile)
            if (.not. isFile) then
               call io_printerr(trim(comm)//" - Input error","file 2 is missing")
               call help()
            endif
            getFile2     = .true.
            f_name_in(2) = trim(line)
         elseif (.not. getRate1) then
            call io_linesplit(line, M(1), CHK, PrintErr_in = .true.)
            if (CHK .ne. CHKOK) then
               call io_printerr(trim(comm)//" - Input error","[RATE 1] - "//trim(line))
               call help()
            endif
            getRate1 = .true.
         elseif (.not. getRate2) then
            call io_linesplit(line, M(2), CHK, PrintErr_in = .true.)
            if (CHK .ne. CHKOK) then
               call io_printerr(trim(comm)//" - Input error","[RATE 2] - "//trim(line))
               call help()
            endif
            getRate2 = .true.
         else
            cycle LK
         endif
      enddo LK
   endif

 end subroutine init

!-----------------------------------------------------------------------------!

 subroutine help()
   implicit none

   write (*,'(A)') " usage    : "//trim(comm)//" [OPTION] [CHGCAR 1] [CHGCAR 2] [RATE 1] [RATE 2]"
   write (*,'(A)') " default  : [OPTION] = None"
   write (*,'(A)') "            [CHGCAR] = None"
   write (*,'(A)') "            [RATE]   = 1.0"
   write (*,*)
   write (*,'(A)') " [OPTION] : -h : call this help"
   write (*,'(A)') " [CHGCAR] : vasp CHGCAR file, no default"
   write (*,'(A)') " [RATE]   : mixing rate of CHGCAR 1 and 2, default = 1.0"
   write (*,*)
   stop

 end subroutine help

end program vasp_chg_sum
