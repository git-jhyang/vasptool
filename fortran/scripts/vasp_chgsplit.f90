program chgsplit

implicit none

integer, parameter :: fp = selected_real_kind(15), &
&                     fc = 10, fo1 = 11, fo2 = 12, fo3 = 13, fo4 = 14

character :: fname*50, bc100*100, dis*100, at_name*100, co_typ*100, bc200*200, fout*50
logical :: isfile
double precision :: lc, lvec(3,3), tl, tt !, vol
integer :: i, j, k, sp_num, at_tnum, tpoint, grid(3), gridref(3), &
&          chk, t1, t2, trate, col, low, rem, nblock 

integer, allocatable :: at_num(:)
double precision, allocatable :: chgval(:,:,:), coord(:,:), bufchgval(:)

call getarg(1, fname)
call system_clock(t1, trate)

inquire (file = fname, exist = isfile)
if (.not. isfile) then
   if (fname(1:2) .eq. '  ') then
      write (*,*) " !!! ERROR : Invalid input format !!! "
      write (*,*) "             chgsplit [filename]"
   else
      write (*,*) " !!! ERROR : File does not exist !!!"
   endif
   stop
endif

open (fc, file = fname)

read(fc,'(A)') dis            ! file description
read(fc,'(A)') bc100            ! lattice constant
read(bc100,*) lc

do i = 1,3
   read(fc,'(A)') bc100         ! lattice vectors
   read(bc100,*) lvec(i,:)
enddo
lvec(:,:) = lvec(:,:)*lc
call lvec_refin(lvec)
!vol = lvec(1,1)*lvec(2,2)*lvec(3,3)

read(fc,'(A)') at_name            ! atom names
sp_num = 0
do i = 1, 99                    ! # of atom spices
   if ((at_name(i:i) .ne. ' ') .and. (at_name(i+1:i+1) .eq. ' ')) then
      sp_num = sp_num + 1
   endif
enddo
allocate (at_num(sp_num))

read (fc,'(A)') bc100           ! # of atoms
read (bc100,*) at_num(:)

at_tnum = 0
do i = 1, sp_num
   at_tnum = at_tnum + at_num(i)
enddo
allocate (coord(3,at_tnum))

read(fc,'(A)') co_typ            ! cartesian or direct

do i = 1, at_tnum
   read(fc,'(A)') bc100         ! coordination of atom
   read(bc100,*) coord(:,i)
enddo

read(fc,'(A)') bc100            ! empty line
read(fc,'(A)') bc100            ! grid
read(bc100, *) gridref(:)
nblock = 1

tpoint = gridref(1)*gridref(2)*gridref(3)        ! total data number

write (*,*)
write (*,*) " Total grid points : ",tpoint

read(fc,'(A)') bc200            ! read first line
col = 0

do i = 1, 199                   ! count column
   if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
      col = col + 1
   endif
enddo

low = ceiling(float(tpoint)/float(col))         ! calculate low
rem = mod(tpoint,col)                           ! remainder

!write (*,*) low, col, rem
!write (*,*) tpoint, at_tnum, gridref(:)

allocate (chgval(col, 4, low))

if (rem .eq. 0) then
   rem = col
endif

allocate (bufchgval(rem))

read(bc200,*) chgval(:,1,1)             ! put first data

write (*,*)
write (*,*) " Reading data block 1 ..."
write (*,*) "  0   |   20    |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2
do i = 2, low-1
   read(fc,'(A)') bc200
   read(bc200,*) chgval(:,1,i)          ! data read,  THE FIRST BLOCK
!   write (*,'(I12,A,I12,4ES15.6E3)') i,'/',low, chgval(i,1:4,nblock)
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol
write (*,*)
read(fc,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval(1:rem,1,low) = bufchgval(:)!/vol

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tl
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

check1: do
101 continue
    read(fc,'(A)', end = 100) bc100
    j = 0
    do i = 1,99
       if ((bc100(i:i) .ne. ' ') .and. (bc100(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then
       read(bc100, *, err = 101) grid(:)
       if ((gridref(1) .eq. grid(1)) .and. (gridref(2) .eq. grid(2)) &
&         .and. (gridref(3) .eq. grid(3))) then
          nblock = 2
          exit check1
       endif
    endif
end do check1
call system_clock(t1, trate)

write (*,*)
write (*,*) " Reading data block 2 ..."
write (*,*) "  0   |   20    |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc,'(A)') bc200
   read(bc200,*) chgval(:,2,i)          ! data read, THE SECOND BLOCK
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval(1:rem,2,low) = bufchgval(:)!/vol

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tt + tl
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

!write (*,*) 'done?'
check2: do
102 continue
    read(fc,'(A)', end = 100) bc100
    j = 0
!    write (*,*) 'here?'
    do i = 1,99
       if ((bc100(i:i) .ne. ' ') .and. (bc100(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then
       read(bc100, *, err = 102) grid(:)
       if ((gridref(1) .eq. grid(1)) .and. (gridref(2) .eq. grid(2)) & 
&         .and. (gridref(3) .eq. grid(3))) then
          nblock = 3
          exit check2
       endif
    endif
end do check2
call system_clock(t1,trate)

write (*,*)
write (*,*) " Reading data block 3 ..."
write (*,*) "  0   |   20    |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc,'(A)') bc200
   read(bc200,*) chgval(:,3,i)          ! data read, THE THIRD BLOCK
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval(1:rem,3,low) = bufchgval(:)!/vol

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tt + tl
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'


check3: do
103 continue
    read(fc,'(A)', end = 100) bc100
    j = 0
    do i = 1,99
       if ((bc100(i:i) .ne. ' ') .and. (bc100(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then
       read(bc100, *, err = 103) grid(:)
       if ((gridref(1) .eq. grid(1)) .and. (gridref(2) .eq. grid(2)) &
&         .and. (gridref(3) .eq. grid(3))) then
          nblock = 4
          exit check3
       endif
    endif
end do check3

call system_clock(t1, trate)
write (*,*)
write (*,*) " Reading data block 4 ..."
write (*,*) "  0   |    20   |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc,'(A)') bc200
   read(bc200,*) chgval(:,4,i)          ! data read, THE FOURTH BLOCK
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)

!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval(1:rem,4,low) = bufchgval(:)!/vol

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tt + tl
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

100 continue
close(fc)

call system_clock(t1, trate)
if (nblock .eq. 1) then
   write (*,*)           " !!! WARNING : There is 1 data block !!!"
   write (*,'(A,A20,A)') "               Check your '",fname(1:20),"' file"
   stop
endif

if (nblock .eq. 3) then
   write (*,*)           " !!! WARNING : There are 3 data block !!!"
   write (*,'(A,A20,A)') "               Check your '",fname(1:20),"' file"
   stop
endif

do i = 1, 49
   if (fname(i:i) .eq. '.') then
      fname(i:i) = '_'
   endif
   if (fname(i:i) .eq. ' ') then
      exit
   endif
enddo

fout(1:i) = fname(1:i)
if (nblock .eq. 2) then
   fout(i:i+3) = '_tot'
   open (fo1, file = fout(1:i+3), action = 'write')
   fout(i:i+3) = '_mag'
   open (fo2, file = fout(1:i+3), action = 'write')
   fout(i:i+2) = '_up'
   open (fo3, file = fout(1:i+2), action = 'write')
   fout(i:i+4) = '_down'
   open (fo4, file = fout(1:i+4), action = 'write')
elseif (nblock .eq. 4) then
   fout(i:i+3) = '_tot'
   open (fo1, file = fout(1:i+3), action = 'write')
   fout(i:i+4) = '_magX'
   open (fo2, file = fout(1:i+4), action = 'write')
   fout(i:i+4) = '_magY'
   open (fo3, file = fout(1:i+4), action = 'write')
   fout(i:i+4) = '_magZ'
   open (fo4, file = fout(1:i+4), action = 'write')
endif

write(fo1,'(A)') dis
write(fo2,'(A)') dis
write(fo3,'(A)') dis
write(fo4,'(A)') dis

write(fo1,'(A)') ' 1.00000000000000'
write(fo2,'(A)') ' 1.00000000000000'
write(fo3,'(A)') ' 1.00000000000000'
write(fo4,'(A)') ' 1.00000000000000'

do i = 1,3
   write(fo1, '(3F12.6)') lvec(i,:)
   write(fo2, '(3F12.6)') lvec(i,:)
   write(fo3, '(3F12.6)') lvec(i,:)
   write(fo4, '(3F12.6)') lvec(i,:)
enddo

write (fo1,'(A)') at_name
write (fo2,'(A)') at_name
write (fo3,'(A)') at_name
write (fo4,'(A)') at_name

do i = 1, sp_num
   write (fo1,'(I6)', advance = 'no') at_num(i)
   write (fo2,'(I6)', advance = 'no') at_num(i)
   write (fo3,'(I6)', advance = 'no') at_num(i)
   write (fo4,'(I6)', advance = 'no') at_num(i)
enddo

write (fo1,*)
write (fo2,*)
write (fo3,*)
write (fo4,*)

write (fo1,'(A)') co_typ
write (fo2,'(A)') co_typ
write (fo3,'(A)') co_typ
write (fo4,'(A)') co_typ

do i = 1, at_tnum
   write(fo1, '(3F12.6)') coord(:,i)
   write(fo2, '(3F12.6)') coord(:,i)
   write(fo3, '(3F12.6)') coord(:,i)
   write(fo4, '(3F12.6)') coord(:,i)
enddo

write (fo1,*)
write (fo2,*)
write (fo3,*)
write (fo4,*)

write (fo1, '(3I12)') gridref(:)
write (fo2, '(3I12)') gridref(:)
write (fo3, '(3I12)') gridref(:)
write (fo4, '(3I12)') gridref(:)

write (*,*)
write (*,*)
write (*,*) " Writing datas ..."
write (*,*) "  0   |    |    |    |    |    |    |    |    |  100 % "
write (*,'(A)', advance = 'no') "   "

k = 2
do i = 1, low-1
   if (int(i*100/(low-1)) .eq. k) then
      write (*, '(A1)', advance = 'no') "*"
      k = k + 2
   endif
   if (nblock .eq. 2) then
      do j = 1, col
         write(fo1,'(ES12.4E2)', advance = 'no') chgval(j,1,i) 
         write(fo2,'(ES12.4E2)', advance = 'no') chgval(j,2,i)
         write(fo3,'(ES12.4E2)', advance = 'no') (chgval(j,1,i) + chgval(j,2,i))*0.5_fp
         write(fo4,'(ES12.4E2)', advance = 'no') (chgval(j,1,i) - chgval(j,2,i))*0.5_fp
      enddo
   elseif (nblock .eq. 4) then
      do j = 1, col
         write(fo1,'(ES12.4E2)', advance = 'no') chgval(j,1,i)
         write(fo2,'(ES12.4E2)', advance = 'no') chgval(j,2,i)
         write(fo3,'(ES12.4E2)', advance = 'no') chgval(j,3,i)
         write(fo4,'(ES12.4E2)', advance = 'no') chgval(j,4,i)
      enddo
   endif
   write (fo1,*)
   write (fo2,*)
   write (fo3,*)
   write (fo4,*)
enddo

if (nblock .eq. 2) then
   do j = 1, rem
      write(fo1,'(ES12.4E2)', advance = 'no') chgval(j,1,low)
      write(fo2,'(ES12.4E2)', advance = 'no') chgval(j,2,low)
      write(fo3,'(ES12.4E2)', advance = 'no') (chgval(j,1,low) + chgval(j,2,low))*0.5_fp
      write(fo4,'(ES12.4E2)', advance = 'no') (chgval(j,1,low) - chgval(j,2,low))*0.5_fp
   enddo
elseif (nblock .eq. 4) then
   do j = 1, rem
      write(fo1,'(ES12.4E2)', advance = 'no') chgval(j,1,low)
      write(fo2,'(ES12.4E2)', advance = 'no') chgval(j,2,low)
      write(fo3,'(ES12.4E2)', advance = 'no') chgval(j,3,low)
      write(fo4,'(ES12.4E2)', advance = 'no') chgval(j,4,low)
   enddo
endif

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tt + tl
write (*,*)
write (*,'(A,F6.2,A)') '   Time : ', tt, ' min'

write (fo1,*)
write (fo2,*)
write (fo3,*)
write (fo4,*)
write (*,*)

close(fo1)
close(fo2)
close(fo3)
close(fo4)

deallocate (at_num, coord, chgval, bufchgval)

contains
subroutine lvec_refin(lvec)
implicit none

integer, parameter :: fp = selected_real_kind(15)
double precision, intent(inout) :: lvec(3,3)
double precision :: llvec(3), angle(3), dot(3)
integer :: i, j

do i = 1, 3
   llvec(i) = sqrt(lvec(i,1)**2 + lvec(i,2)**2 + lvec(i,3)**2)
enddo

dot(1) = lvec(2,1)*lvec(3,1) + lvec(2,2)*lvec(3,2) + lvec(2,3)*lvec(3,3)
dot(2) = lvec(1,1)*lvec(3,1) + lvec(1,2)*lvec(3,2) + lvec(1,3)*lvec(3,3)
dot(3) = lvec(1,1)*lvec(2,1) + lvec(1,2)*lvec(2,2) + lvec(1,3)*lvec(2,3)

angle(1) = acos(dot(1)/(llvec(2)*llvec(3)))
angle(2) = acos(dot(2)/(llvec(3)*llvec(1)))
angle(3) = acos(dot(3)/(llvec(1)*llvec(1)))

lvec(1,1) = llvec(1)
lvec(1,2) = 0.0_fp
lvec(1,3) = 0.0_fp

lvec(2,1) = llvec(2)*cos(angle(3))
lvec(2,2) = sqrt(llvec(2)**2 - lvec(2,1)**2)
lvec(2,3) = 0.0_fp

lvec(3,1) = llvec(3)*cos(angle(2))
lvec(3,2) = (cos(angle(1))*llvec(2)*llvec(3) - lvec(2,1)*lvec(3,1))/lvec(2,2)
lvec(3,3) = sqrt(llvec(3)**2 - lvec(3,1)**2 - lvec(3,2)**2)

do i = 1, 3
   do j = 1, 3
      if (abs(lvec(i,j)) .lt. 0.000000000001_fp ) then
         lvec(i,j) = 0.0_fp
      endif
   enddo
enddo

end subroutine
end program
