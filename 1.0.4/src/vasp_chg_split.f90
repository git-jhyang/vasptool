program vasp_chg_split

use constants,  only :  dp,             &
                        CHKOK,          &
                        CHKERR,         &
                        LINELEN

use util,       only :  time_formatting

use iotool,     only :  io_printerr

use vasp,       only :  chgcar,         &
                        poscar,         &
                        chg_get,        &
                        chg_sum,        &
                        chg_write,      &
                        poscar_read,    &
                        poscar_write

implicit none

type(chgcar), dimension(4) :: DatChg
type(poscar)               :: DatPos

character(LEN=LINELEN) :: comm, f_name_in, f_name_out
integer                :: u, cnt, nBlock, iBlock
logical                :: isCollin

real(kind=dp)          :: eTime = 0d0, deTime

character(LEN=20), dimension(2) :: tform
real(kind=dp),     dimension(2) :: M

call init()

open (unit=u, file=f_name_in, status='old')
call poscar_read(u, DatPos)

do iBlock = 1, nBlock
   call chg_get(u, DatChg(iBlock))
   if (DatChg(iBlock)%CHK .ne. CHKOK) then
      call io_printerr(trim(comm)//" - data is missing","wrong input file or option")
      write (*,'(5x,"Required block : ",I1," (file : ",I1," blocks)")') nBlock, iBlock-1
      stop
   endif
   eTime = eTime + DatChg(iBlock)%eTime
   call time_formatting(DatChg(iBlock)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"read vasp CHGCAR blocks - ",I1," - ",A,"( ",A," )")') &
        iBlock, trim(adjustl(tform(1))), trim(adjustl(tform(2)))
enddo
close(u)

if (isCollin) then
   M = [0.5d0, 0.5d0]
   call chg_sum(DatChg(1), DatChg(2), DatChg(3), M_in=M)
   M(2) = -0.5d0
   call chg_sum(DatChg(1), DatChg(2), DatChg(4), M_in=M)
   eTime = eTime + DatChg(3)%eTime + DatChg(4)%eTime
   call time_formatting(DatChg(3)%eTime+DatChg(4)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"spin up and down are seperated - ",A,"( ",A," )")') &
        trim(adjustl(tform(1))), trim(adjustl(tform(2)))
endif

do iBlock = 1, 4
   write (f_name_out,'(A,"_",I1)') trim(f_name_in), iBlock

   open (unit=u, file=f_name_out)
   DatPos%isLong = .false.
   call poscar_write(u, DatPos)
   call chg_write(u, DatChg(iBlock))
   close (u)

   deTime = DatChg(iBlock)%eTime
   eTime  = eTime + deTime
   call time_formatting(deTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"charge density ",I1," (",A,") is written - ",A,"( ",A," )")') &
        iBlock, trim(f_name_out), trim(adjustl(tform(1))), trim(adjustl(tform(2))) 
enddo

write (*,*)
if (isCollin) then
   write (*,'(2x,A,"_1 : UP+DOWN")') trim(f_name_in)
   write (*,'(2x,A,"_2 : UP-DOWN")') trim(f_name_in)
   write (*,'(2x,A,"_3 : UP")') trim(f_name_in)
   write (*,'(2x,A,"_4 : Down")') trim(f_name_in)
else
   write (*,'(2x,A,"_1 : SUM")') trim(f_name_in)
   write (*,'(2x,A,"_2 : X")') trim(f_name_in)
   write (*,'(2x,A,"_3 : Y")') trim(f_name_in)
   write (*,'(2x,A,"_4 : Z")') trim(f_name_in)
endif
write (*,*)

contains

!-----------------------------------------------------------------------------!

 subroutine init()
   implicit none

   integer :: ik, k
   logical :: getFile, getOpt, isFile

   character(LEN=LINELEN) :: line

   k        = iargc()
   isCollin = .true.
   nBlock   = 2
   getFile  = .false.
   getopt   = .false.

   call getarg(0, comm)

   if (k .eq. 0) then
      call io_printerr(trim(comm)//" - Input error","need more than 1 arg")
      call help()
   else
      LK : do ik = 1, k
         call getarg(ik, line)
         select case (trim(line))
         case ('-h')
            call help()
         case ('-ncl')
            if (getOpt) cycle LK
            isCollin = .false.
            nBlock   = 4
            getOpt   = .true.
            cycle LK
         case ('-std')
            if (getOpt) cycle LK
            isCollin = .true.
            nBlock   = 2
            getOpt   = .true.
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

   if (.not. getFile) then
      call io_printerr(trim(comm)//" - Input error","file is missing")
      call help()
   endif

 end subroutine init

!-----------------------------------------------------------------------------!

 subroutine help()

   implicit none

   write (*,'(A)') " usage    : "//trim(comm)//" [OPTION] [CHGCAR]"
   write (*,'(A)') " default  : [OPTION] = -std"
   write (*,'(A)') "            [CHGCAR] = None"
   write (*,*)
   write (*,'(A)') " [OPTION] : -h   : call this help"
   write (*,'(A)') "          : -std : collinear spin polarized calculation, default"
   write (*,'(A)') "                   (OUTPUTs : up + down, up - down, up, down)"
   write (*,'(A)') "          : -ncl : non-collinear calculation"
   write (*,'(A)') "                   (OUTPUTs : sum, x, y, z)"
   write (*,'(A)') " [CHGCAR] : vasp CHGCAR file, no default"
   write (*,*)
   stop

 end subroutine help

end program vasp_chg_split
