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

character(LEN=LINELEN) :: comm, f_name, line, lineGrid
integer                :: k, u, cnt, nBlock, iBlock
logical                :: isFile, isCollin

real(kind=dp)          :: eTime = 0d0

character(LEN=20), dimension(2) :: tform
real(kind=dp),     dimension(2) :: M

k = iargc()
call getarg(0, comm)
if (k .ne. 2) then
   call call_err("Input error","need 2 args")
   call help(comm)
endif

call getarg(1, line)
call getarg(2, f_name)

read (line, *) u

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

inquire (file=f_name, exist=isFile)
if (.not. isFile) then
   call call_err("Input error","file is missing")
   call help(comm)
endif

open (unit=u, file=f_name, status='old')
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
         if (cnt .gt. 100) then
            call call_err("Data is missing","wrong Input file or option")
            write (*,'(5x,"Required block : ",I1," (file : ",I1," blocks")') nBlock, iBlock
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
   M = [0.5d0, -0.5d0]
   call chg_sum(DatChg(1), DatChg(2), DatChg(4), M)
   eTime = eTime + DatChg(3)%eTime + DatChg(4)%eTime
   call time_formatting(DatChg(3)%eTime+DatChg(4)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"spin up and down are seperated - ",A,"( ",A," )")') &
        trim(adjustl(tform(1))), trim(adjustl(tform(2)))
endif

do iBlock = 1, 4
   write (f_name,'("CHGCAR_",I1)') iBlock
   open (unit=u, file=f_name)
   call poscar_write(u, DatPos)
   write (u,*)
   call chg_write(u, DatChg(iBlock))
   eTime = eTime + DatChg(1)%eTime
   call time_formatting(DatChg(iBlock)%eTime, tform(1))
   call time_formatting(eTime, tform(2))
   write (*,'(2x,"charge density ",I1," is written - ",A,"( ",A," )")') &
        iBlock, trim(adjustl(tform(1))), trim(adjustl(tform(2))) 
   close (u)
enddo

contains

!-----------------------------------------------------------------------------!

 subroutine help(comm_in)
   implicit none
   character(LEN=*), intent(in) :: comm_in

   write (*,'(3A)') " usage    : ", trim(adjustl(comm_in)), " [option] [VASP_CHGCAR]"
   write (*,'(A)')  " default  : None, need 2 args"
   write (*,'(A)')  " [option] : CHGCAR calculation method, no default "
   write (*,'(A)')  "          : 1 - collinear spin polarized calculation"
   write (*,'(A)')  "                (OUTPUTs : up + down, up - down, up, down)"
   write (*,'(A)')  "          : 2 - non-collinear calculation"
   write (*,'(A)')  "                (OUTPUTs : sum, x, y, z)"
   write (*,'(A)')  " [CHGCAR] : vasp CHGCAR file, no default"
   stop

 end subroutine help

end program vasp_locpot
