!07/26/2017:
!swapped fort.300 to flux.dat
!write nt to actual time in fs

implicit none

integer  :: i,nt
integer  :: errcode
real*8   :: dummy_r1,dummy_r2,N_sum,dt
real*8,allocatable   :: time(:),flux(:) 


  !Open file dp_file and count the number of points it contains
  open(1,file='flux.dat',action='read',status='old',iostat=errcode)
  if (errcode/=0) then
    !write(*,*) 'Error: cannot access file ','fort.300'
   ! write(*,*) 'Error code: ', errcode
   ! write(*,*) 'Error: cannot access file ','flux.dat'
    stop
  endif
  nt=-1
  do while (errcode==0)
    read(1,*,iostat=errcode) dummy_r1,dummy_r2
    nt=nt+1 
  enddo
  !write(6,*) nt
  write(6,*) nt*0.007
  rewind(1)

  allocate(time(nt),flux(nt))
  do i=1,nt
    read(1,*)time(i),flux(i)
  enddo
  close(1)
  dt=time(3)-time(2)

  N_sum=0.d0
  write(400,*)0.d0,N_sum
  do i=2,nt
    N_sum=N_sum+0.5d0*(flux(i-1)+flux(i))*dt
    write(400,*)time(i),N_sum
  enddo

end
















