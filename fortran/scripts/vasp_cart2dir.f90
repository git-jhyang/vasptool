program vasp_cart2dir

implicit none

integer, parameter :: fp = selected_real_kind(15), fi = 10, fo = 11

integer :: i, j, sp_num, at_tnum
integer, allocatable :: at_num(:)

character :: fn1*50, fn2*50, bc200*200, selc(3)*1
character, allocatable :: at_name(:)*2

double precision :: lcont, lvec(3,3), lenvec(3), ang(3), coord(3)
logical :: isfile, isselc, iscart

call getarg(1, fn1)

inquire (file = fn1, exist = isfile)
if (.not. isfile) then
   if (fn1(1:2) .eq. '  ') then
      write (*,*) " !!! ERROR : Invalid input format !!! "
      write (*,*) "             vasp_cart2dir [filename]"
   else
      write (*,*) " !!! ERROR : File does not exist !!!"
   endif
   stop
endif

fn2 = fn1

do i = 1, 50
   if ( fn1(i:i) .eq. ' ' ) then
      fn2(i:i+4) = '_dir '
      exit
   endif
enddo

open (fi, file = fn1, status = 'old')
open (fo, file = fn2, action = 'write')

read (fi, '(A)') bc200
write (fo, '(A30)') bc200

read (fi, '(A)') bc200
read (bc200, *) lcont

write (fo, '(F22.15)') 1.0_fp

do i = 1, 3
   read (fi, '(A)') bc200
   read (bc200,*) lvec(:,i)
   lvec(:,i) = lvec(:,i)*lcont
enddo

do i = 1, 3
   lenvec(i) = sqrt(lvec(1,i)**2 + lvec(2,i)**2 + lvec(3,i)**2)
enddo

call angle(lvec(:,2), lvec(:,3), ang(1))
call angle(lvec(:,3), lvec(:,1), ang(2))
call angle(lvec(:,1), lvec(:,2), ang(3))

lvec(1,1) = lenvec(1)
lvec(2,1) = 0.0_fp
lvec(3,1) = 0.0_fp

lvec(1,2) = lenvec(2)*cos(ang(3))
lvec(2,2) = sqrt(lenvec(2)**2 - lvec(1,2)**2)
lvec(3,2) = 0.0_fp

lvec(1,3) = lenvec(3)*cos(ang(2))
lvec(2,3) = (cos(ang(1))*lenvec(2)*lenvec(3) - lvec(1,2)*lvec(1,3))/lvec(2,2)
lvec(3,3) = sqrt(lenvec(3)**2 - lvec(1,3)**2 - lvec(2,3)**2)

do i = 1, 3
   do j = 1, 3
      if (abs(lvec(j,i)) .lt. 0.000000000001_fp ) then
         lvec(j,i) = 0.0_fp
      endif
   enddo
   write (fo, '(3F22.15)') lvec(:,i)
enddo

read (fi, '(A)') bc200

sp_num = 0
do i = 1, 199
   if ((bc200(i:i) .ne. ' ') .and. (bc200(i+1:i+1) .eq. ' ')) then
      sp_num = sp_num + 1
   endif
enddo

allocate (at_num(sp_num), at_name(sp_num))
read (bc200, *) at_name

do i = 1, sp_num
   write (fo, '(A5)', advance = 'no') at_name(i)
enddo
write (fo, *)

read (fi, '(A)') bc200
read (bc200, *) at_num

at_tnum = 0
do i = 1, sp_num
   write (fo, '(I5)', advance = 'no') at_num(i)
   at_tnum = at_tnum + at_num(i)
enddo
write (fo, *)

deallocate (at_num, at_name)

read (fi, '(A)') bc200

do i = 1, 199
   if (bc200(i:i) .ne. ' ') then
      exit
   endif
enddo

isselc = .false.
if ((bc200(i:i) .eq. 'S') .or. (bc200(i:i) .eq. 's')) then
   isselc = .true.
   write (fo, '(A)') "Selective Dynamics"
   read (fi, '(A)') bc200
endif

do i = 1, 199
   if (bc200(i:i) .ne. ' ') then
      exit
   endif
enddo

if ((bc200(i:i) .eq. 'C') .or. (bc200(i:i) .eq. 'c')) then
   iscart = .true.
else if ((bc200(i:i) .eq. 'D') .or. (bc200(i:i) .eq. 'd')) then
   write (*,*) " !!! WRANRING : Your 'POSCAR' is already 'Direct' !!!"
   iscart = .false.
endif
write (fo, '(A)') "Direct"

do i = 1, at_tnum
   read (fi, '(A)') bc200
   if ( isselc ) then
      read (bc200, *) coord, selc
   else
      read (bc200, *) coord
   endif
   if ( iscart .and. isselc) then
      call cart_to_dir(lvec, coord)
      write (fo, '(3F22.15, 3A5)') coord, selc
   elseif ( iscart .and. (.not. isselc )) then
      call cart_to_dir(lvec, coord)
      write (fo, '(3F22.15)') coord
   else if ( isselc ) then
      write (fo, '(3F22.15, 3A5)') coord, selc
   else
      write (fo, '(3F22.15)') coord
   endif
enddo

contains

subroutine cart_to_dir(lvec, coord)
implicit none

integer, parameter :: fp = SELECTED_REAL_KIND(16)
double precision, intent(in) :: lvec(3,3)
double precision, intent(inout) :: coord(3)

double precision :: br, conv_mat(3,3), lenvec(3), ang(3), vol, coord_new(3)
integer :: i

do i = 1, 3
   lenvec(i) = sqrt(lvec(1,i)**2 + lvec(2,i)**2 + lvec(3,i)**2)
end do

call angle(lvec(:,2), lvec(:,3), ang(1))
call angle(lvec(:,3), lvec(:,1), ang(2))
call angle(lvec(:,3), lvec(:,2), ang(3))

! http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm

br = cos(ang(1))*cos(ang(2))*cos(ang(3))
br = sqrt(1.0_fp - cos(ang(1))**2 - cos(ang(2))**2 - cos(ang(3))**2 + 2*br)
vol = lenvec(1)*lenvec(2)*lenvec(3)*br

conv_mat(1,1) = 1.0_fp/lenvec(1)
conv_mat(1,2) = -cos(ang(3))/(lenvec(1)*sin(ang(3)))
conv_mat(1,3) = lenvec(2)*lenvec(3)*((cos(ang(3))*(cos(ang(1))-cos(ang(2))*cos(ang(3)))/sin(ang(3)))-cos(ang(2))*sin(ang(3)))/vol

conv_mat(2,1) = 0.0_fp
conv_mat(2,2) = 1.0_fp/(lenvec(2)*sin(ang(3)))
conv_mat(2,3) = -lenvec(1)*lenvec(3)*(cos(ang(1)) - cos(ang(2))*cos(ang(3)))/(sin(ang(3))*vol)

conv_mat(3,1) = 0.0_fp
conv_mat(3,2) = 0.0_fp
conv_mat(3,3) = lenvec(1)*lenvec(2)*sin(ang(3))/vol

coord_new(:) = 0.0_fp

do i = 1, 3
   coord_new(1) = coord_new(1) + conv_mat(1,i)*coord(i)
   coord_new(2) = coord_new(2) + conv_mat(2,i)*coord(i)
   coord_new(3) = coord_new(3) + conv_mat(3,i)*coord(i)
enddo

do i = 1, 3
   do while ( coord_new(i) .lt. 0.0_fp )
      coord_new(i) = coord_new(i) + 1.0_fp
   enddo
   do while ( coord_new(i) .ge. 1.0_fp )
      coord_new(i) = coord_new(i) - 1.0_fp
   enddo
enddo

coord(:) = coord_new(:)

end subroutine

subroutine angle(vec1, vec2, ang)

implicit none

integer, parameter :: fp = SELECTED_REAL_KIND(16)
double precision, intent(in) :: vec1(3), vec2(3)
double precision, intent(out) :: ang

double precision :: br, len1, len2, vecs(3,2)
double precision, parameter :: pi = 3.1415926535897932_fp

integer :: i

len1 = sqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
len2 = sqrt(vec2(1)**2 + vec2(2)**2 + vec2(3)**2)

if (len1 .ne. 0) then
   vecs(:,1) = vec1(:)/len1
endif

if (len2 .ne. 0) then
   vecs(:,2) = vec2(:)/len2
endif

br = 0.0_fp
do i = 1, 3
   br = br + vecs(i,1)*vecs(i,2)
enddo

!write (*,*) br
if ( int(abs(br)*1000000.0)/1000000 .eq. 0 ) then
   ang = pi*0.5_fp
elseif ( int(abs(br)*1000000.0)/1000000 .eq. 1 ) then
   ang = 0.0_fp
else
   ang = acos(br)
endif

end subroutine angle

end program
