program chgsplit

implicit none

integer, parameter :: fp = selected_real_kind(15), &
&                     fc1 = 10, fc2 = 11, fo = 12

character :: fname1*50, fname2*50, bc100*100, dis*100, &
&            at_name*100, co_typ*100, bc200*200, fout*50
logical :: isfile
double precision :: lc, lvec(3,3,2), tl, tt, fac1, fac2 !, vol
integer :: i, j, k, l, sp_num, at_tnum, tpoint, grid(3,2), gridref(3,2), &
&          chk, t1, t2, trate, col, low, rem, nblock 

integer, allocatable :: at_num(:,:)
character, allocatable :: at_names(:,:)*2
double precision, allocatable :: chgval1(:), chgval2(:), bufchgval(:)

k = iargc()

fac1 = 1.0_fp
fac2 = 1.0_fp
if ((k .lt. 2) .or. (k .gt. 4)) then 
   write (*,*) " !!! ERROR : Invalid input format !!! "
   write (*,*) "             vasp_chgsum [file1] [file2] (fac1) (fac2)"
   write (*,*)
   write (*,*) "             CHGCAR_sum = [file1]*(fac1) + [file2]*(fac2) "
   write (*,*) "             Default fac1 and fac2 : 1.0000"
   write (*,*)
   stop
elseif (k .eq. 3) then
   call getarg(3,bc100)
   read (bc100,*) fac1
   fac2 = fac1
elseif (k .eq. 4) then
   call getarg(3, bc100)
   read (bc100,*) fac1
   call getarg(4, bc100)
   read (bc100,*) fac2
endif

call getarg(1, fname1)
call getarg(2, fname2)
call system_clock(t1, trate)

inquire (file = fname1, exist = isfile)
if (.not. isfile) then
   write (*,*) " !!! ERROR : '",fname1(1:20),"' does not exist !!!"
   stop
endif

inquire (file = fname2, exist = isfile)
if (.not. isfile) then
   write (*,*) " !!! ERROR : '",fname2(1:20),"' does not exist !!!"
   stop
endif

write (*,*)
write (*,*) "  'CHGCAR_sum' will be written "
write (*,'(A,A20,A,F6.3)') "    Charge density 1 : ",fname1," * ",fac1
write (*,'(A,A20,A,F6.3)') "    Charge density 2 : ",fname2," * ",fac2

open (fc1, file = fname1)
open (fc2, file = fname2)
open (fo, file = 'CHGCAR_sum', action = 'write')

read(fc1,'(A)') dis            ! file description
read(fc2,*)
write (fo, '(A)') dis

read(fc1,'(A)') bc100            ! lattice constant
read(fc2,*)
write (fo, '(A)') '1.000000000000'

read(bc100,*) lc

do i = 1,3
   read(fc1,'(A)') bc100         ! lattice vectors
   read(bc100,*) lvec(i,:,1)
   read(fc2,'(A)') bc100
   read(bc100,*) lvec(i,:,2)
   do j = 1, 3
      if (lvec(i,j,1) .ne. lvec(i,j,2)) then
         write (*,*) " !!! ERROR : Lattice does not match !!!"
         stop
      endif
   enddo
enddo

lvec(:,:,1) = lvec(:,:,1)*lc

call lvec_refin(lvec)
do i = 1, 3
   write (fo, '(3F12.6)') lvec(i,:,1)
enddo

!vol = lvec(1,1)*lvec(2,2)*lvec(3,3)

read(fc1,'(A)') at_name         ! atom names
read(fc2,'(A)') bc100
write (fo, '(A)') at_name

sp_num = 0
do i = 1, 99                    ! # of atom spices
   if ((at_name(i:i) .ne. ' ') .and. (at_name(i+1:i+1) .eq. ' ')) then
      sp_num = sp_num + 1
   endif
enddo

allocate (at_num(sp_num,2), at_names(sp_num,2))

read (bc100, *) at_names(:,2)
read (at_name, *) at_names(:,1)

do i = 1, sp_num
   if (at_names(i,1) .ne. at_names(i,2)) then
      write (*,*) " !!! ERROR : Atom types does not match !!! "
      stop
   endif
enddo

read (fc1,'(A)') bc100           ! # of atoms
read (bc100,*) at_num(:,1)
read (fc2,'(A)') bc100           ! # of atoms
read (bc100,*) at_num(:,2)
write (fo, '(A)') bc100

at_tnum = 0
do i = 1, sp_num
   at_tnum = at_tnum + at_num(i,1)
   if (at_num(i,1) .ne. at_num(i,2)) then
      write (*,*) " !!! ERROR : Atom numbers does not match !!! "
      stop
   endif
enddo

read(fc1,'(A)') bc100            ! cartesian or direct
read(fc2,*)
write (fo, '(A)') bc100

do i = 1, at_tnum
   read(fc1,'(A)') bc100         ! coordination of atom
   read(fc2,*)
   write(fo, '(A)') bc100
enddo

read(fc1,'(A)') bc100            ! empty line
read(fc2,*)
write(fo, '(A)') bc100

read(fc1,'(A)') bc100            ! grid
read(bc100, *) gridref(:,1)
read(fc2,'(A)') bc100
read(bc100, *) gridref(:,2)

nblock = 1

tpoint = gridref(1,1)*gridref(2,1)*gridref(3,1)        ! total data number
i      = gridref(1,2)*gridref(2,2)*gridref(3,2)

if (tpoint .ne. i) then
   write (*,*) " !!! ERROR : Number of grid point does not match !!!"
   write (*,*) "        ",fname1(1:15)," : ",tpoint
   write (*,*) "        ",fname2(1:15)," : ",i
   stop
endif
write(fo, '(A)') bc100

write (*,*)
write (*,*) " Total grid points : ",tpoint

read(fc1,'(A)') bc200            ! read first line of file1
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

allocate (chgval1(col), chgval2(col))

if (rem .eq. 0) then
   rem = col
endif

allocate (bufchgval(rem))

read(bc200,*) chgval1(:)             ! put first data
read(fc2,'(A)') bc200
read(bc200,*) chgval2(:)

do k = 1, col
   write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
enddo
write (fo, *)

write (*,*)
write (*,*) " Processing data block 1 ..."
write (*,*) "  0   |   20    |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "

j = 2
do i = 2, low-1
   read(fc1,'(A)') bc200
   read(bc200,*) chgval1(:)          ! data read,  THE FIRST BLOCK
   read(fc2,'(A)') bc200
   read(bc200,*) chgval2(:)
   do k = 1, col
      write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
   enddo
   write (fo, *)
!   write (*,'(I12,A,I12,4ES15.6E3)') i,'/',low, chgval(i,1:4,nblock)
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol
write (*,*)

read(fc1,'(A)') bc200                   ! last line of file 1
read(bc200,*) bufchgval(:)

chgval1(1:rem) = bufchgval(:)!/vol

read(fc2,'(A)') bc200                   ! last line of file 2
read(bc200,*) bufchgval(:)

chgval2(1:rem) = bufchgval(:)!/vol

do k = 1, rem
   write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
enddo
write (fo, *)

write (fo, *)

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
tt = tl
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

check1: do
101 continue
    read(fc1,'(A)', end = 100) bc200
    read(fc2,*, end = 100)
    j = 0
    do i = 1,199
       if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then  ! skip while read grid values
       read(bc200, *, err = 101) grid(:,1)
       if ((gridref(1,1) .eq. grid(1,1)) .and. (gridref(2,1) .eq. grid(2,1)) &
&         .and. (gridref(3,1) .eq. grid(3,1))) then
          write (fo, '(A)') bc200
          nblock = 2
          exit check1
       endif
    endif
end do check1
call system_clock(t1, trate)

write (*,*)
write (*,*) " Processing data block 2 ..."
write (*,*) "  0   |   20    |   40    |    60   |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc1,'(A)') bc200
   read(bc200,*) chgval1(:)          ! data read, THE SECOND BLOCK
   read(fc2,'(A)') bc200
   read(bc200,*) chgval2(:)          ! data read, THE SECOND BLOCK
   do k = 1, col
      write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
   enddo
   write (fo, *)
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc1,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval1(1:rem) = bufchgval(:)!/vol

read(fc2,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval2(1:rem) = bufchgval(:)!/vol

do k = 1, rem
   write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
enddo
write (fo, *)
write (fo, *)

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

check2: do
102 continue
    read(fc1,'(A)', end = 100) bc200
    read(fc2, *, end = 100)
    j = 0
!    write (*,*) 'here?'
    do i = 1,199
       if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then  ! skip while read grid values
       read(bc200, *, err = 102) grid(:,1)
       if ((gridref(1,1) .eq. grid(1,1)) .and. (gridref(2,1) .eq. grid(2,1)) &
&         .and. (gridref(3,1) .eq. grid(3,1))) then
          write (fo, '(A)') bc200
          nblock = 3
          exit check2
       endif
    endif
end do check2
call system_clock(t1,trate)

write (*,*)
write (*,*) " Processing data block 3 ..."
write (*,*) "  0   |   20    |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc1,'(A)') bc200
   read(bc200,*) chgval1(:)          ! data read, THE THIRD BLOCK
   read(fc2,'(A)') bc200
   read(bc200,*) chgval2(:)          ! data read, THE THIRD BLOCK
   do k = 1, col
      write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
   enddo
   write (fo, *)
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)
!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc1,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval1(1:rem) = bufchgval(:)!/vol

read(fc2,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval2(1:rem) = bufchgval(:)!/vol

do k = 1, rem
   write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
enddo
write (fo, *)

write (fo, *)

call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

check3: do
103 continue
    read(fc1,'(A)', end = 100) bc200
    read(fc2,*,end = 100)
    j = 0
    do i = 1,199
       if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
          j = j + 1
       endif
    enddo
    if (j .eq. 3) then  ! skip while read grid values
       read(bc200, *, err = 103) grid(:,1)
       read(fc2, *, end = 100)
       if ((gridref(1,1) .eq. grid(1,1)) .and. (gridref(2,1) .eq. grid(2,1)) &
&         .and. (gridref(3,1) .eq. grid(3,1))) then
          write (fo, '(A)') bc200
          nblock = 4
          exit check3
       endif
    endif
end do check3

call system_clock(t1, trate)

write (*,*)
write (*,*) " Processing data block 4 ..."
write (*,*) "  0   |    20   |   40    |   60    |   80    |  100 % "
write (*,'(A)', advance = 'no') "   "
j = 2

do i = 1, low-1
   read(fc1,'(A)') bc200
   read(bc200,*) chgval1(:)          ! data read, THE FOURTH BLOCK
   read(fc2,'(A)') bc200
   read(bc200,*) chgval2(:)          ! data read, THE FOURTH BLOCK
   do k = 1, col
      write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
   enddo
   write (fo, *)
   if (int(i*100/(low-1)) .eq. j) then
      write (*, '(A1)', advance = 'no') "*"
      j = j + 2
   endif
enddo
write (*,*)

!chgval(1:low-1,:,nblock) = chgval(1:low-1,:,nblock)/vol

read(fc1,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval1(1:rem) = bufchgval(:)!/vol

read(fc2,'(A)') bc200
read(bc200,*) bufchgval(:)

chgval2(1:rem) = bufchgval(:)!/vol

do k = 1, rem
   write (fo, '(ES12.4E2)', advance = 'no') chgval1(k)*fac1 + chgval2(k)*fac2
enddo
write (fo, *)

write (fo, *)


call system_clock(t2, trate)
tl = float(t2-t1)/(float(trate)*60.0_fp)
write (*,'(A,F6.2,A)') '   Time : ', tl, ' min'

100 continue
close(fc1)
close(fc2)
close(fo)

deallocate (at_names, at_num, chgval1, chgval2, bufchgval)

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

dot(1) = dot_product(lvec(2,:), lvec(3,:))
dot(2) = dot_product(lvec(1,:), lvec(3,:))
dot(3) = dot_product(lvec(1,:), lvec(2,:))

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
