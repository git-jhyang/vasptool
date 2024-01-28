program vasp_locpot

use constants,  only :  dp, LINELEN

use util,       only :  call_err,       &
                        time_formatting

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

character(LEN=LINELEN) :: comm, f_name_in, f_name_out, line, lineGrid
integer                :: k, u, cnt, nBlock, iBlock
logical                :: isFile, isCollin

real(kind=dp)          :: eTime = 0d0, deTime

character(LEN=20), dimension(2) :: tform
real(kind=dp),     dimension(2) :: M

k = iargc()
call getarg(0, comm)
if (k .eq. 0) then
   call call_err("Input error","need more than 1 arg")
   call help(comm)
elseif (k .eq. 1) then
   call getarg(1, f_name_in)
   u = 1
else
   call getarg(1, f_name_in)
   call getarg(2, line)
   read (line, *) u
endif

if (u .eq. 1) then
   isCollin = .True.
   nBlock   = 2
elseif (u .eq. 2) then
   isCollin = .False.
   nBlock   = 4
else
   call call_err("Input error","[CASE] should be 1 or 2")
   call help(comm)
endif

inquire (file=f_name_in, exist=isFile)
if (.not. isFile) then
   call call_err("Input error","file is missing")
   call help(comm)
endif

open (unit=u, file=f_name_in, status='old')
call poscar_read(u, DatPos)

do iBlock = 1, nBlock
   if (iBlock .eq. 1) then
      read (u, *)
      read (u, '(A)') lineGrid
   else
      cnt = 0
      do while (.True.)
         read (u, '(A)') line
         if (line .eq. lineGrid) exit
         cnt = cnt + 1
         if (cnt .gt. 1000) then
            call call_err("Data is missing","wrong Input file or option")
            write (*,'(5x,"Required block : ",I1," (file : ",I1," blocks)")') nBlock, iBlock-1
            stop
         endif
      enddo
   endif
   read (lineGrid, *) DatChg(iBlock)%nGrid

   allocate (DatChg(iBlock)%chg(DatChg(iBlock)%nGrid(1), &
             DatChg(iBlock)%nGrid(2), DatChg(iBlocK)%nGrid(3)))
   call chg_get(u, DatChg(iBlock))
   eTime = eTime + DatChg(iBlock)%eTime
   call time_formatting(DatChg(iBlock)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"read vasp CHGCAR blocks - ",I1," - ",A,"( ",A," )")') &
        iBlock, trim(adjustl(tform(1))), trim(adjustl(tform(2)))
enddo
close(u)

if (isCollin) then
   do iBlock = 3, 4
      DatChg(iBlock)%nGrid = DatChg(1)%nGrid
      allocate (DatChg(iBlock)%chg(DatChg(iBlock)%nGrid(1), &
                DatChg(iBlock)%nGrid(2), DatChg(iBlocK)%nGrid(3)))
   enddo
   M = [0.5d0, 0.5d0]
   call chg_sum(DatChg(1), DatChg(2), DatChg(3), M)
   M(2) = -0.5d0
   call chg_sum(DatChg(1), DatChg(2), DatChg(4), M)
   eTime = eTime + DatChg(3)%eTime + DatChg(4)%eTime
   call time_formatting(DatChg(3)%eTime+DatChg(4)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"spin up and down are seperated - ",A,"( ",A," )")') &
        trim(adjustl(tform(1))), trim(adjustl(tform(2)))
endif

do iBlock = 1, 4
   write (f_name_out,'(A,"_",I1)') trim(f_name_in), iBlock

   open (unit=u, file=f_name_out)
   call poscar_write(u, DatPos, .false.)
   write (u,*)
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

contains

!-----------------------------------------------------------------------------!

 subroutine help(comm_in)
   implicit none
   character(LEN=*), intent(in) :: comm_in

   write (*,'(3A)') " usage    : ", trim(adjustl(comm_in)), " [CHGCAR] [CASE]"
   write (*,'(A)')  " default  : [CHGCAR] = None"
   write (*,'(A)')  "            [CASE]   = 1"
   write (*,*)
   write (*,'(A)')  " [CHGCAR] : vasp CHGCAR file, no default"
   write (*,'(A)')  " [CASE]   : CHGCAR calculation method, default = 1 "
   write (*,'(A)')  "          : 1 - collinear spin polarized calculation"
   write (*,'(A)')  "                (OUTPUTs : up + down, up - down, up, down)"
   write (*,'(A)')  "          : 2 - non-collinear calculation"
   write (*,'(A)')  "                (OUTPUTs : sum, x, y, z)"
   write (*,*)
   stop

 end subroutine help

end program vasp_locpot
