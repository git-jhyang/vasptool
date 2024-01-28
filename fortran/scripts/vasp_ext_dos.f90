program vasp_ext_dos
implicit none

integer, parameter :: fd = 10, fo = 12, fp = 13

logical :: isfile

integer :: i, j, k, cnt, chk, ispin, lorbit, lsorbit, natom, ngrid, ncol, bi4(4), sp_num
integer, allocatable :: at_num(:), at_bnum(:)

character :: bc300*300, bc200*200, bc50*50, bc40*40, bc4(4)*50, aname*2
character, allocatable :: at_name(:)*2

double precision :: br, br3(3), br4(4), br7(7), br10(10), br19(19), &
&                   ef
double precision, allocatable :: dos(:,:,:), energy(:)


inquire (file = 'DOSCAR', exist = isfile)
if (.not. isfile) then
   write (*,*) " !!! ERROR : 'DOSCAR' does not exist !!! "
   stop
endif

bc50 = 'ISPIN'
call ext_in(bc50,bc40)
read (bc40, *) ispin
if (ispin .eq. 999) then
   ispin = 1
endif


bc50 = 'LSORBIT'
call ext_in(bc50,bc40)
if ((bc40(1:3) .eq. '.T.') .or. (bc40(1:6) .eq. '.TRUE.')) then
   lsorbit = 1
   write (*,*) " !!! ERROR : 'LSORBIT' is not supported yet, SORRY !!!"
   stop
else
   lsorbit = 0
endif

bc50 = 'LORBIT'
call ext_in(bc50,bc40)
read (bc40,*) lorbit
if (lorbit .eq. 999) then
   lorbit = 0
endif

open (fd, file = 'DOSCAR')
read (fd, '(A)') bc200
read (bc200, *) bi4
natom = bi4(1)

read (fd, *)
read (fd, *)
read (fd, *)
read (fd, *)

read (fd, '(A)') bc200
read (bc200, *) br4

ngrid = br4(3)
ef = br4(4)

open (fo, file = 'dos_total.dat', action = 'write')
if (ispin .eq. 1) then
   write (fo, *) '#   E-EF      TDOS'
else
   write (fo, *) '#   E-EF     TDOSu     TDOSd'
endif
write(fo,*) '####################################################################'

read(fd,*)
allocate (energy(ngrid-1))
do i = 1, ngrid-1
   read(fd, '(A)') bc200
   read(bc200,*) br3
   energy(i) = br3(1) - ef
   write(fo, '(2F10.4)', advance = 'no') energy(i), br3(2)
   if (ispin .eq. 2) then
      write(fo, '(F10.4)', advance = 'no') br3(3)
   endif
   write(fo,*)
enddo

close(fo)
if (lorbit .eq. 0) then
   close(fd)
   stop
endif

! inquire( directory = 'dos', exist = isfile)          ! work for ifort
inquire( file = trim('dos')//'/.', exist = isfile)   ! work for gfortran


if (.not. isfile) then
   call system('mkdir dos')
   write (*,*) " 'dos' directory is created "
endif

inquire(file = "POSCAR", exist = isfile)
if (.not. isfile) then
   write (*,*) " !!! ERROR : 'POSCAR' does not exist !!!"
   stop
endif

open (fp, file = 'POSCAR')
read (fp,*)
read (fp,*)
read (fp,*)
read (fp,*)
read (fp,*)
read (fp, '(A)') bc200

sp_num = 0
do i = 1, 199
   if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
      sp_num = sp_num + 1
   endif
enddo

ncol = 3*ispin
if (lorbit .ne. 10) then
   ncol = ncol * 3
endif
ncol = ncol + 1

allocate (at_num(sp_num), at_bnum(sp_num), at_name(sp_num))
allocate (dos(ncol,ngrid-1,0:sp_num))

dos(:,:,:) = 0.000000000000000000

read (bc200, *) at_name
read (fp, '(A)') bc200
read (bc200, *) at_bnum

j = 0
do i = 1, sp_num
   at_num(i) = at_bnum(i) + j
   j = j + at_bnum(i)
enddo

if (natom .ne. j) then
   write (*,*) " !!! ERROR : check your 'DOSCAR' and 'POSCAR' files !!!"
   write (*,*) "             # of total atoms are different in two files "
   write (*,*)
   stop
endif

cnt = 1
chk = 2

write (*,*)
write (*,*) " Writing lm-decomposed DOSCARs in 'dos' directory ..."
write (*,*) "  0   10   20   30   40   50   60   70   80   90  100 % "
write (*,*) "  |    |    |    |    |    |    |    |    |    |    |"
write (*,'(A)', advance = 'no') "   *"


do i = 1, natom
   if (i .gt. at_num(cnt)) then
      cnt = cnt + 1
   endif
   bc40 = ' '
   bc40(1:8) = 'dos/dos_'
   write(bc40(9:11),'(I3.3)') i
   bc40(12:12) = '_'
   aname = at_name(cnt)
   bc40(13:14) = aname(1:2)
   if (bc40(14:14) .eq. ' ') then
      bc40(14:17) = '.dat'
   else
      bc40(15:18) = '.dat'
   endif
   read (fd,*)
   read (fd,*)
   open (fo, file = bc40, action = 'write')
   if ((lorbit .eq. 10) .and. (ispin .eq. 1))then
      write (fo, '(A)')  '#     E-EF     s-DOS     p-DOS     d-DOS'
   elseif ((lorbit .eq. 10) .and. (ispin .eq. 2))then
      write (fo, '(A)')  '#     E-EF    s-DOSu    s-DOSd    p-DOSu    p-DOSd    d-DOSu    d-DOSd'
   elseif ((lorbit .ne. 10) .and. (ispin .eq. 1))then
      write (fo, '(A)', advance = 'no') '#     E-EF     s-DOS    px-DOS    py-DOS    pz-DOS   dxy-DOS   dyz-DOS'
      write (fo, '(A)')  '   dz2-DOS   dxz-DOS   dx2-DOS'
   elseif ((lorbit .ne. 10) .and. (ispin .eq. 2))then
      write (fo, '(A)', advance = 'no') '#     E-EF    s-DOSu    s-DOSd   px-DOSu   px-DOSd   py-DOSu   py-DOSd'
      write (fo, '(A)', advance = 'no') '   pz-DOSu   pz-DOSd  dxy-DOSu  dxy-DOSd  dyz-DOSu  dyz-DOSd  dz2-DOSu'
      write (fo, '(A)')  '  dz2-DOSd  dxz-DOSu  dxz-DOSd  dx2-DOSu  dx2-DOSd'
   endif
   write (fo, '(A)', advance = 'no') '#########################################################################'
   write (fo, '(A)') '#########################################################################'
   do j = 1, ngrid-1
      write (fo, '(F10.4)', advance = 'no') energy(j)
      read (fd,'(A)') bc300
      read (bc300, *) dos(:, j, 0)
      do k = 2, ncol
         write (fo, '(F10.4)', advance = 'no') dos(k,j,0)
      enddo
      dos(:,j,cnt) = dos(:,j,cnt) + dos(:,j,0)
      write (fo,*)
   enddo
   if (int(i*100/natom) .eq. chk) then
      write (*, '(A1)', advance = 'no') "*"
      chk = chk + 2
   elseif (int(i*100/natom) .gt. chk) then
      do while (int(i*100/natom) .ge. chk)
         write (*, '(A1)', advance = 'no') "*"
         chk = chk + 2
      enddo
   endif
   close (fo)
enddo
write (*,*)
write (*,*)

do i = 1, sp_num
   bc40 = ' '
   bc40(1:10) = 'dos/dos_#_'
   aname = at_name(i)
   bc40(11:12) = aname(1:2)
   if (bc40(12:12) .eq. ' ') then
      bc40(12:15) = '.dat'
   else
      bc40(13:16) = '.dat'
   endif
   open (fo, file = bc40, action = 'write')
   if ((lorbit .eq. 10) .and. (ispin .eq. 1))then
      write (fo, '(A)')  '#     E-EF     s-DOS     p-DOS     z-DOS'
   elseif ((lorbit .eq. 10) .and. (ispin .eq. 2))then
      write (fo, '(A)')  '#     E-EF    s-DOSu    s-DOSd    p-DOSu    p-DOSd    z-DOSu    z-DOSd'
   elseif ((lorbit .ne. 10) .and. (ispin .eq. 1))then
      write (fo, '(A)', advance = 'no') '#     E-EF     s-DOS    px-DOS    py-DOS    pz-DOS   dxy-DOS   dyz-DOS'
      write (fo, '(A)')  '   dz2-DOS   dxz-DOS   dx2-DOS'
   elseif ((lorbit .ne. 10) .and. (ispin .eq. 2))then
      write (fo, '(A)', advance = 'no') '#     E-EF    s-DOSu    s-DOSd   px-DOSu   px-DOSd   py-DOSu   py-DOSd'
      write (fo, '(A)', advance = 'no') '   pz-DOSu   pz-DOSd  dxy-DOSu  dxy-DOSd  dyz-DOSu  dyz-DOSd  dz2-DOSu'
      write (fo, '(A)')  '  dz2-DOSd  dxz-DOSu  dxz-DOSd  dx2-DOSu  dx2-DOSd'
   endif
   write (fo, '(A)', advance = 'no') '#########################################################################'
   write (fo, '(A)')  '#########################################################################'
   do j = 1, ngrid-1
      write (fo, '(F10.4)', advance = 'no') energy(j)
      do k = 2, ncol
         write (fo, '(F10.4)', advance = 'no') dos(k,j,i)
      enddo
      write (fo,*)
   enddo
   close (fo)
enddo

deallocate (at_num, at_bnum, at_name, dos, energy)
close(fd)

contains

subroutine ext_in(command, val)
implicit none

character, intent(in)  :: command*50
character, intent(out) :: val*40

integer, parameter :: fp = selected_real_kind(15), fin = 11
character :: bc200*200, bc100*100, bc50(3)*50
integer :: i
logical :: isin

inquire (file = 'INCAR', exist = isin)
if (.not. isin) then
   write (*,*) " !!! ERROR : 'INCAR' does not exist !!!"
   stop
endif

do i = 1, 49
   if ((command(i:i) .ne. ' ') .and. (command(i+1:i+1) .eq. ' ')) then
      exit
   endif
enddo

val = '999'
open (fin, file = 'INCAR')
do 
100 continue
   read(fin,'(A)', end = 101) bc100
   read(bc100,*, err = 100, end = 100) bc50
!   write (*,*) '2', bc50(1), command(1:i)
   if (bc50(1) .eq. command(1:i)) then
      read(bc50(3),*) val
      exit
   endif
enddo
101 continue

close (fin)

end subroutine
end program
