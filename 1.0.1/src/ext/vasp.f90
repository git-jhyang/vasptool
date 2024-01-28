module vasp

use constants,  only :  dp, LINELEN, CHKEOF, CHKOK
use util,       only :  cart_to_dir

implicit none

type, public :: poscar
   character(LEN=LINELEN)                        :: Des
   real(kind=dp),    dimension(3,3)              :: LattVec
   integer                                       :: nSpec = 0
   integer                                       :: nAtom
   integer,          dimension(:),   allocatable :: nAtomSpec
   character(LEN=2), dimension(:),   allocatable :: NameAtoms
   logical                                       :: isCart = .false.
   logical                                       :: isSelc = .false.
   real(kind=dp),    dimension(:,:), allocatable :: coo
   logical,          dimension(:,:), allocatable :: selc
end type poscar

type, public :: chgcar
   integer,       dimension(3)                  :: nGrid
   real(kind=dp), dimension(:,:,:), allocatable :: chg
   real(kind=dp)                                :: eTime
end type chgcar

contains

!-----------------------------------------------------------------------------!

 subroutine chg_get(u, DatChg)
   implicit none
   integer,      intent(in)    :: u
   type(chgcar), intent(inout) :: DatChg

   character(LEN=LINELEN) :: line
   integer                :: chk, nLine, nCol, nPnt, nRem
   integer                :: ii, ix1, ix2, ix3, ix4, iy1, iy2, iz1, iz2
   integer                :: t1, t2, trate
   logical                :: isOpen

   nCol = 0

   call system_clock(t1, trate)
 
   read (u,'(A)') line
    do ii = 1, LINELEN-1
       if ((line(ii:ii) .ne. ' ') .and. (line(ii+1:ii+1) .eq. ' ')) then
          nCol = nCol + 1
       endif
    enddo

   nPnt  = DatChg%nGrid(1)*DatChg%nGrid(2)*DatChg%nGrid(3)
   nLine = ceiling(float(nPnt)/float(nCol))

   ix1 = 1
   ix3 = 1
   iy1 = 1
   iz1 = 1
   iz2 = 1
   do ii = 2, nLine
      ix2 = ix1 + nCol - 1
      if (ix2 .le. DatChg%nGrid(1)) then
         read (line, *) DatChg%chg(ix1:ix2,iy1,iz1)
      else
         ix4 = ix2 - DatChg%nGrid(1) 
         ix2 = DatChg%nGrid(1)
         iy2 = iy1 + 1
         if (iy1 .eq. DatChg%nGrid(2)) then
            iy2 = 1
            iz2 = iz1 + 1
         endif
         read (line, *) DatChg%chg(ix1:ix2,iy1,iz1), DatChg%chg(ix3:ix4,iy2,iz2)
         ix2 = ix4
         iy1 = iy2
         iz1 = iz2
      endif
      ix1 = ix2 + 1
      if (ix2 .eq. DatChg%nGrid(1)) then
         ix1 = 1
         iy1 = iy1 + 1
         if (iy1 .gt. DatChg%nGrid(2)) then
            iy1 = 1
            iz1 = iz1 + 1
            iz2 = iz1
         endif
      endif
      read (u,'(A)') line
   enddo
   read (line, *) DatChg%chg(ix1:DatChg%nGrid(1),iy1,iz1)

   call system_clock(t2, trate)
   DatChg%eTime = float(t2-t1)/float(trate)

 end subroutine chg_get

!-----------------------------------------------------------------------------!

 subroutine chg_write(u, DatChg, long)
   implicit none
   integer,           intent(in)    :: u
   type(chgcar),      intent(inout) :: DatChg
   logical, optional, intent(in)    :: long

   integer :: t1, t2, trate
   integer :: ix, iy, iz, cnt, cnt_max
   logical :: isLong

   isLong = .false.
   if (present(long)) isLong = long

   cnt_max = 10
   if (isLong) cnt_max = 5

   call system_clock(t1, trate)

   write (u,'(3I6)') DatChg%nGrid
   if (      isLong) write (u,201) (((DatChg%chg(ix,iy,iz), ix=1,DatChg%nGrid(1)), &
                        iy=1,DatChg%nGrid(2)), iz=1,DatChg%nGrid(3))
   if (.not. isLong) write (u,200) (((DatChg%chg(ix,iy,iz), ix=1,DatChg%nGrid(1)), &
                        iy=1,DatChg%nGrid(2)), iz=1,DatChg%nGrid(3))

200 format(10E13.5E2)
201 format(5E20.12E2)

   call system_clock(t2, trate)
   DatChg%eTime = float(t2-t1)/float(trate)

 end subroutine chg_write

!-----------------------------------------------------------------------------!

 subroutine chg_sum(DatChg1, DatChg2, DatChgOut, Mult)
   implicit none
   type(chgcar),                          intent(in)    :: DatChg1
   type(chgcar),                          intent(in)    :: DatChg2
   type(chgcar),                          intent(inout) :: DatChgOut
   real(kind=dp), dimension(2), optional, intent(in)    :: Mult

   real(kind=dp), dimension(2) :: M
   integer                     :: ix, iy, iz
   integer                     :: t1, t2, trate

   call system_clock(t1, trate)

   M(:) = 1d0

   if (present(Mult)) M = Mult

   do iz = 1, DatChgOut%nGrid(3)
   do iy = 1, DatChgOut%nGrid(2)
   do ix = 1, DatChgOut%nGrid(1)
      DatChgOut%chg(ix,iy,iz) = M(1)*DatChg1%chg(ix,iy,iz) + M(2)*DatChg2%chg(ix,iy,iz) 
   enddo
   enddo
   enddo

   call system_clock(t2, trate)
   DatChgOut%eTime = float(t2-t1)/float(trate)

 end subroutine chg_sum

!-----------------------------------------------------------------------------!

 subroutine poscar_read(u, DatPos)
   implicit none
   integer,      intent(in)  :: u
   type(poscar), intent(out) :: DatPos

   integer                :: ii, chk
   logical                :: iserr, iscart
   character(LEN=LINELEN) :: line
   real(kind=dp)          :: LattCont

   chk = CHKOK

   read (u,'(A)', iostat=chk, end=100) DatPos%Des ! file comment
   read (u,'(A)', iostat=chk, end=100) line
    read (line, *, err=100) LattCont

   do ii = 1, 3
      read (u,'(A)', iostat=chk, end=100) line
       read (line, *, err=100) DatPos%LattVec(:,ii)
       DatPos%LattVec(:,ii) = DatPos%LattVec(:,ii)*LattCont
   enddo

   read (u,'(A)', iostat=chk, end=100) line ! name of atoms, vasp 5 style
    do ii = 1, LINELEN-1
       if ((line(ii:ii) .ne. ' ') .and. (line(ii+1:ii+1) .eq. ' ')) &
                DatPos%nSpec = DatPos%nSpec + 1 
    enddo
    allocate (DatPos%NameAtoms(DatPos%nSpec), DatPos%nAtomSpec(DatPos%nSpec))
    read (line, *, err=100) DatPos%NameAtoms

   read(u,'(A)', iostat=chk, end=100) line  ! number of atoms of each species
    read (line, *, err=100) DatPos%nAtomSpec
    DatPos%nAtom = sum(DatPos%nAtomSpec)
    allocate (DatPos%coo(3,DatPos%nAtom))

   read(u,'(A)', iostat=chk, end=100) line  ! read coordination type
    if (line(1:1) .eq. 'S' .or. line(1:1) .eq. 's') then
       DatPos%isSelc = .true.
       read (u, '(A)', iostat=chk, end=100) line
       allocate (DatPos%selc(3,DatPos%nAtom))
    endif
    if (line(1:1) .eq. 'C' .or. line(1:1) .eq. 'c') DatPos%isCart = .true.

   do ii = 1, DatPos%nAtom      ! read coordination
      read(u,'(A)', iostat=chk, end=100) line
       if (      DatPos%isSelc) read(line, *, err=100) DatPos%coo(:,ii), DatPos%selc(:,ii)
       if (.not. DatPos%isSelc) read(line, *, err=100) DatPos%coo(:,ii), DatPos%selc(:,ii)
   enddo

   return

100 continue
   if (chk .le. CHKEOF) write (*,*) 'touched EoF during read POSCAR'
   if (chk .gt. CHKEOF) write (*,*) 'error while read POSCAR'

 end subroutine poscar_read

!-----------------------------------------------------------------------------!

 subroutine poscar_write(u, DatPos)
   implicit none
   integer,      intent(in) :: u
   type(poscar), intent(in) :: DatPos

   integer                     :: ii
   real(kind=dp), dimension(3) :: coo

   write (u,'(A)') trim(adjustl(DatPos%Des))
   write (u,'(A)') '1.00000000000'
   do ii = 1, 3
      write (u, '(3F22.15)') DatPos%LattVec(:,ii)
   enddo
   do ii = 1, DatPos%nSpec
      write (u,'(A4)',advance='no') DatPos%NameAtoms(ii)
   enddo
   write (u,*)
   do ii = 1, DatPos%nSpec
      write (u,'(I4)',advance='no') DatPos%nAtomSpec(ii)
   enddo
   write (u,*)
   write (u,'(A)') 'Direct'
   do ii = 1, DatPos%nAtom
      coo = DatPos%coo(:,ii)
      if (DatPos%isCart) call cart_to_dir(DatPos%LattVec, coo)
      write (u,'(3F22.15)') coo
   enddo

 end subroutine poscar_write

!-----------------------------------------------------------------------------!

end module vasp

