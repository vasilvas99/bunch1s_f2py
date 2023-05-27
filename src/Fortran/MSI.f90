program MS_driver
   implicit none

   character*8 dir_prefix
   common/drxx/dir_prefix
   real*8 xk
   common/txx/xk

   integer N, M
   integer num_steps
   integer ipasses, iwrite, i
   real*8 bdef
   real*8, dimension(:, :), allocatable :: trajectories
   real*8, dimension(:), allocatable :: time

   character(*), parameter::input_path = "step_trajectories.dat"
   character(*), parameter::meta_input_path = "meta_step_trajectories.dat"
   dir_prefix = "ms1_out/"
   call setup_output_dir(dir_prefix)
   call setup_MS_I

   call read_trajectory_metadata(meta_input_path, N, M, bdef)
   allocate (trajectories(N, M - 1))
   allocate (time(N))
   num_steps = M - 1
   call read_trajectory(input_path, N, M, time, trajectories)
   ipasses = 1
   iwrite = 1
   ! run the actual MSI
   do i = 1, N
      xk = time(i)
      call MS_I(trajectories(i, :), num_steps, bdef, ipasses, iwrite)
   end do
   deallocate (trajectories)
   deallocate (time)
end program MS_driver

subroutine setup_output_dir(output_dir_path)
   character*(*), intent(in)::output_dir_path
   logical dir_exists
   inquire (FILE=output_dir_path, EXIST=dir_exists)
   if (.NOT. dir_exists) then
      call system('mkdir '//output_dir_path) ! should be okay in both windows and linux
   end if
end subroutine

subroutine read_trajectory_metadata(meta_file_path, N, M, bdef)
   character*(*), intent(in)::meta_file_path
   integer, intent(out)::N
   integer, intent(out)::M
   real*8, intent(out)::bdef

   open (unit=15, file=meta_file_path, action='read')
      read (15, *) ! skip header
      read (15, *) N, M, bdef
   close (15)
   print *, "Rows to read: ", N, " Columns to Read: ", M, "Bunch def: ", bdef
end subroutine

subroutine read_trajectory(file_path, N, M, time, trajectories)
   character*(*), intent(in)::file_path
   integer, intent(in)::N
   integer, intent(in)::M
   real*8, intent(out)::time(N)
   real*8, intent(out)::trajectories(N, M - 1)

   integer i, j

   open (16, file=file_path, action="read")
      read (16, *) ! skip header
      do i = 1, N
         read (16, *) time(i), trajectories(i, :)
      end do
   close (16)
   print *, "Trajectories loaded"
end subroutine

! Here be dragons
subroutine setup_MS_I
   character*12 w_h, tn_tm, st_tm, rut
   common/fnlx/w_h, tn_tm, st_tm, rut
   character*21 tname
   character*8 alngdr
   real*8 wo(99), ho(99)
   real*8 mindo(99)
   common/orxx/wo, ho, mindo
   common/drxx/alngdr
   print *, "Output dir: ", alngdr
   do iord = 1, 99
      wo(iord) = 0.0d0
      mindo(iord) = 0.0d0
      ho(iord) = 0.0d0
   end do
   w_h = '!!w-h__.dat'
   tn_tm = '!tn-tm__.dat'
   rut = '!ru-tm__.dat'
   st_tm = '!st-tm__.dat'
   tname = alngdr//w_h
   open (10, file=tname)
   write (10, *)
   close (10)
   tname = alngdr//rut
   open (10, file=tname)
   write (10, *)
   close (10)
   tname = alngdr//tn_tm
   open (10, file=tname)
   write (10, *)
   close (10)
   tname = alngdr//st_tm
   open (10, file=tname)
   write (10, *)
   close (10)
   return
end subroutine

subroutine MS_I(Y, M, bdef, ipasses, iwrite)
   ! All the arguments exchanged with MS_I called previously lstat are on INPUT:
!   Y is an 1-D array and contains the coordinates of the
!   M steps
!   bdef is a key concept - it defines a distance between two neighboring steps
!   ipasses gives the number of passes through the same time interval already done while
!   iwrite gives the number of the time step during the current pass
!
!     the way the program writes in files is changed to have the time as a X axis.
!
   integer M
   Real*8 Y(*), DY(M) ! Here I am leaving DM - the distance between the steps to be only an internal array
!        Bunch / Terrace Definitions:
   real*8 bdef
!        common /bdxx/ bdef
   integer ipasses
!        common /ipxx/ ipasses
   integer iwrite
!        common /iwxx/ iwrite
!
   real*8 wo(99), ho(99)
   real*8 mindo(99)
   real*8 w, h
   real*8 mind
   common/orxx/wo, ho, mindo
   real*8 stime
   real*8 dyi, maxbd, mintd, maxtd
!
   real*8 avtd, tw, ntd
   integer Tno, i
   real*8 one, zero
!     h , w, mind
   data one, zero/1.0d0, 0.0d0/
   character*12 w_h, tn_tm, st_tm, rut
   common/fnlx/w_h, tn_tm, st_tm, rut
   character*21 tname
   real*8 xk
   common/txx/xk
   character*8 alngdr
   common/drxx/alngdr
!
!                        ...initializations...
   w = zero
   h = zero
   tw = zero
   ntd = zero
   mind = bdef
   maxbd = zero
   mintd = M*bdef
   maxtd = zero
!

!          Compute the interstep distances, the minimal one and
!          some other quantities that characterize the bunching:
   do i = 2, M
      dyi = Y(i) - Y(i - 1)
      dy(i) = dyi
      if (dyi .le. bdef) then
!     counting the number of distances in bunches ...
         h = h + one
!                    ...and the total bunch width
         w = w + dyi

         if (dyi .lt. mind) mind = dyi
         if (dyi .gt. maxbd) maxbd = dyi
      else
         tw = tw + dyi
         ntd = ntd + one
         if (dyi .lt. mintd) mintd = dyi
         if (dyi .gt. maxtd) maxtd = dyi
      end if
   end do
   dyi = Y(1) - Y(M) + M*bdef
   dy(1) = dyi
   if (dyi .le. bdef) then
      w = w + dyi
      h = h + one
      if (dyi .lt. mind) mind = dyi
      if (dyi .gt. maxbd) maxbd = dyi
   else
      tw = tw + dyi
      ntd = ntd + one
      if (dyi .lt. mintd) mintd = dyi
      if (dyi .gt. maxtd) maxtd = dyi
   end if
   avtd = tw/ntd
   if (maxbd .gt. bdef) then
      write (*, *) 'Alert:maxbd!'
      stop
   end if
   if (mintd .lt. bdef) then
      write (*, *) 'Alert:mintd!'
      stop
   end if
!

!*********************************************************************************************
!      Counting the terraces and their width:

   Tno = 0
   if (dy(1) .gt. bdef) Tno = 1
   do i = 2, M
      dyi = dy(i)
      if (dyi .gt. bdef) then
         if (dy(i - 1) .lt. bdef) Tno = Tno + 1
      end if
   end do
!     check if the terrace in the beginning spannes the boundaries
   if (dy(1) .gt. bdef .and. dy(M) .gt. bdef) Tno = Tno - 1
!
   if (Tno .gt. 0) then
      stime = dfloat(M)/dfloat(Tno)
      tw = tw/dfloat(Tno)
      h = h/dfloat(Tno)
      w = w/dfloat(Tno)
      h = ((ipasses - 1)*ho(iwrite) + h)/float(ipasses)
      ho(iwrite) = h
      w = ((ipasses - 1)*wo(iwrite) + w)/float(ipasses)
      wo(iwrite) = w
      mind = ((ipasses - 1)*mindo(iwrite) + mind)/float(ipasses)
      mindo(iwrite) = mind
!****
      tname = alngdr//w_h
      open (10, file=tname, access='append')
      write (10, 3077) xk, h + 1, w, w/h, mind          ! Time, N, W, <lb>, lmin_g
      close (10)
      tname = alngdr//st_tm
      open (11, file=tname, access='append')
!        write(11,*)
      write (11, 3077) xk, stime, tw, avtd, tw/avtd
      close (11)
      tname = alngdr//tn_tm
      open (11, file=tname, access='append')
!        write(11,*)
      write (11, *) xk, tno
      close (11)
   else
      print *, 'terrace alert!'
   end if
3075 format(4(f16.8, 1x), e12.4)
3077 format(e15.8, 1x, 4(f12.4, 1x))
!
   return
end subroutine
