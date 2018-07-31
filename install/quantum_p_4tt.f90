MODULE quantum

!------------------------------------------------------------------------------------------
!06/06/17:
!Implemented cylinder run_type = 3:
!   Added block for run_type 3 in subroutine init_background_density
!   Added a check_pos_cyl function
!
!06/08/17:
!Changes to SUBROUTINE time_propagator:
!   Deleted tv0
!
!06/09/17:
!Bugs fixed in check_pos_cyl and subroutine init_background_density
!
!06/13/17:
!Charges flipped f=2.d0*hbar/ElectronMass
!
!06/28/17:
!Added factors of 2 to energy(e)
!
!07/07/17:
!Added phase_shift parameter in subroutine time
!phase_shift input will be needed in td_laser.inp
!------------------------------------------------------------------------------------------
  
  real*8,parameter        :: gammaU=-0.1423d0,beta1U=1.0529d0,beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real*8,parameter        :: CU=0.002d0,DU=-0.0116d0
  real*8,parameter        :: E2=14.39964415d0 !1/(4*pi*eps_0), where eps_0=0.005526349868 eV/(Angstrom*V^2)
  real*8,parameter        :: H2M=3.80998174348d0 !(hbar^2)/(2*m_e) in eV*Angstom^2
  real*8,parameter        :: a_B=0.52917720859d0 !Bohr radius
  real*8,parameter        :: hbar=0.658211899d0  !Planck constant over 2*pi in eV*fs
  real*8,parameter        :: Ry=13.60569193d0    !Rydberg constant in eV
  real*8,parameter        :: pi=3.141592653589793d0  
  complex*16,parameter    :: zi=cmplx(0.d0,1.d0)   
  real*8,parameter        :: c_squared=8987551.787368176d0 !in (Angstroms/fs)^2
  real*8,parameter        :: c_light=sqrt(c_squared)
  real*8,parameter        :: ElectronMass=0.5d0*hbar*hbar/H2M !mass of electron in eV*fs^2/Angstrom^2
  
! Run type -- shape of jellium (sphere = 1, diode = 2, cylinder = 3)
  integer           :: run_type

! atom index
  integer           :: atom_index

! Number of lattice points 
  integer           :: N_L(3), N_L_points
   
! Lattice index
  integer,allocatable :: Lattice(:,:),Lattice_inv(:,:,:)
! grid points  
  real*8,allocatable :: grid_point(:,:)  
! grid spacing
  real*8 :: grid_step(3),grid_volume,cell_volume
     
  real*8,allocatable :: V_POT(:),wf(:,:,:),V_exchange(:),rho_bg(:)
  real*8,allocatable :: phi(:),L_phi(:),VH(:),V_ext(:)
  real*8,allocatable :: density(:),density_old(:),H_Phi(:)
  real*8,allocatable :: current_density(:,:)
  
! order of finite  difference
  integer,parameter           :: N_d=4
! max L in multipole expansion
  integer,parameter           :: L_max=4
  real*8,parameter  :: small=1.d-50
    

! Number of orbitals
  integer                     :: N_orbitals
! single_particle energy
  real*8,allocatable          :: sp_energy(:),Psi(:,:)

  integer                     :: N_iteration=4
! Energies 
  real*8                      :: E_hartree,E_exchange
  integer                     :: N_scf_iter=100, N_scf_it_iter=50

! electric field
  double precision              :: E_0,omega,tt,laser_wavelength
  double precision              :: laser_pulse_shift1,laser_pulse_width1
  double precision              :: E_field,A_vec(3)
  logical                       :: dielectric
  double precision              :: kick_strength
  double precision              :: phase_shift
! for time propagation
  complex*16                    :: time_step
  complex*16,allocatable        :: Psi_c(:,:),Psi_g(:,:),Phi_g(:),Phi_r(:),grad_w(:,:),Psi_so(:,:)
! G_vectors
  double precision,allocatable  :: G_vectors(:,:),g_length(:)

  real*8                        :: center(3,2),radius

  integer,parameter             :: energy_file=1001, e_field_file=1002, vec_pot_file=1003, flux_file=1004, dp_file=1005

CONTAINS

SUBROUTINE Exchange_Correlation_potential
!LDA parametrization by Predew and Zunger
!https://journals.aps.org/prb/abstract/10.1103/PhysRevB.23.5048
!http://www.fisica.uniud.it/~giannozz/Corsi/MetNum/LectureNotes/metnum-cap2.pdf
! Equation 2.26
  implicit none
  integer :: i
  real*8  :: rho,c1,c2,rs,rssq,rsln,V_xc
  real*8, parameter :: tol = 1.0d-10
  c1=3.d0/(4.d0*pi)
  c2=4.d0/3.d0*0.4582d0
  do i=1,N_L_points
    rho=Density(i)+1.d-10
    if (rho < tol) then
      rho = tol
    end if
    rs=(c1/rho)**(1.d0/3.d0)
    rs=rs/a_B
    V_xc=-c2/rs
    if(rs>1.d0) then
      rssq=sqrt(rs)
      V_xc=V_xc+gammaU*(1.d0+7.d0/6.d0*beta1U*rssq+4.d0/3.d0*beta2U*rs)/(1.d0+beta1U*rssq+beta2U*rs)**2.d0
    else
      rsln=log(rs)
      V_xc=V_xc+AU*rsln+(BU-AU/3.d0)+2.d0/3.d0*CU*rs*rsln+(2.d0*DU-CU)/3.d0*rs
    endif
    V_Exchange(i)=2.d0*Ry*V_xc
  enddo
END SUBROUTINE Exchange_correlation_potential

SUBROUTINE Exchange_Correlation_Energy
!LDA parametrization by Predew and Zunger
!https://journals.aps.org/prb/abstract/10.1103/PhysRevB.23.5048
  implicit none
  integer :: i
  real*8  :: the_sum,eps_xc,rssq,rho,rsln,rs
  real*8, parameter :: tol=1.0d-10
  the_sum=0.d0
  do i=1,N_L_points
    rho=Density(i)
    if (rho < tol) then
      rho = tol
    end if
    rs=(3.d0/(4.d0*Pi*rho))**(1.d0/3.d0)/a_B
    eps_xc=-.4582d0/rs
    if(rs>1.d0) then
      rssq=sqrt(rs)
      eps_xc=eps_xc+gammaU/(1.d0+beta1U*rssq+beta2U*rs)
    else
      rsln=log(rs)
      eps_xc=eps_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
    endif
    eps_xc=2.d0*Ry*eps_xc
    the_sum=the_sum+density(i)*(eps_xc-V_exchange(i))
  enddo
  E_exchange=the_sum*grid_volume
!  write(6,*)e_exchange
END SUBROUTINE Exchange_correlation_Energy


 SUBROUTINE periodic
!     impose periodic boundary condition
    implicit none
    integer                                 :: k,N(3)
    N=N_L
    do k = 0 ,N_d - 1
      wf(N(1)+k,0:N(2)-1,0:N(3)-1)=wf(k,0:N(2)-1,0:N(3)-1)
      wf(0:N(1)-1,N(2)+k,0:N(3)-1)=wf(0:N(1)-1,k,0:N(3)-1)
      wf(0:N(1)-1,0:N(2)-1,N(3)+k)=wf(0:N(1)-1,0:N(2)-1,k)
    end do     
    do k = - N_d, -1
      wf(k,0:N(2)-1,0:N(3)-1)=wf(k+N(1),0:N(2)-1,0:N(3)-1)
      wf(0:N(1)-1,k,0:N(3)-1)=wf(0:N(1)-1,k+N(2),0:N(3)-1)
      wf(0:N(1)-1,0:N(2)-1,k)=wf(0:N(1)-1,0:N(2)-1,k+N(3))
    end do     
 END SUBROUTINE periodic


 SUBROUTINE solve_poisson
 !solve the Poisson equation
 !equations 13.36 -- 13.37 in Compuational Nanoscience
    implicit none
    integer                                        :: i1,i2,i3,i,mode
    real*8                                         :: f,ff,su1,su2,su3,su       
    complex*16, dimension(:,:,:),allocatable       :: rho
    logical                                        :: inv
    allocate(rho(N_L(1),N_L(2),N_L(3)))
    ff=4.*Pi*E2
    su1=0.d0;su2=0.d0;su3=0.d0      
    do i=1,N_L_points
      i1=Lattice(1,i)+1
      i2=Lattice(2,i)+1
      i3=Lattice(3,i)+1
      rho(i1,i2,i3)=Density(i)
    end do
!     fourier transformed density
    mode=1
    call cfftw(rho,N_L(1),N_L(2),N_L(3),mode)
    do i=1,N_L_points
      i1=Lattice(1,i)+1
      i2=Lattice(2,i)+1
      i3=Lattice(3,i)+1
      if(G_length(i).gt.1.d-10) then
        f=ff/G_length(i)
        rho(i1,i2,i3)=f*rho(i1,i2,i3)
      else
!         G=0
        rho(i1,i2,i3)=0.  
      end if
    end do
!     transformation to real space 
    mode=-1
    call cfftw(rho,N_L(1),N_L(2),N_L(3),mode)
    do i=1,N_L_points
      i1=Lattice(1,i)+1
      i2=Lattice(2,i)+1
      i3=Lattice(3,i)+1
      VH(i)=rho(i1,i2,i3)
    end do
    deallocate(rho)
 END SUBROUTINE solve_poisson

SUBROUTINE init_lattice
  integer :: k1,k2,k3,num,i,k
  ! Setup the lattice bookkeeping
  num=0
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      do k3=1,N_L(3)
        num=num+1 
        Lattice(1,num)=k1-1
        Lattice(2,num)=k2-1
        Lattice(3,num)=k3-1
        Lattice_inv(k1,k2,k3)=num
        grid_point(1,num)=(k1-1)*grid_step(1)
        grid_point(2,num)=(k2-1)*grid_step(2)
        grid_point(3,num)=(k3-1)*grid_step(3)
      enddo
    enddo
  enddo
END SUBROUTINE init_lattice

SUBROUTINE cfftw(chd,n1,n2,n3,mode)
    ! Fourier Transform
!     Interface with FFTW3 library
    implicit none
!     Input/Output variables:
    integer, intent(in) :: n1,n2,n3,mode
    complex*16, intent(inout) :: chd(n1,n2,n3)
!     Work variables:
    integer            :: dirn,flg,len,succ
    integer(kind=8)    :: plan
    complex*16         :: fac
    character (len=10) :: file_name
!     parameters from fftw3.f
    integer, parameter :: fftw_forward = -1
    integer, parameter :: fftw_backward = +1
    integer, parameter :: fftw_estimate = 64
!     ---------------------------------------------------------------
    len = n1*n2*n3 
    if (mode .gt. 0) then
       file_name = 'forwardfft'
       dirn = fftw_forward
    else
       file_name = 'backwrdfft'
       dirn = fftw_backward
    endif
    call dfftw_plan_dft_3d(plan, n1, n2, n3,chd,chd,dirn,fftw_estimate)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    if (mode .gt. 0) then
       fac = 1.d0/len
       chd=fac*chd
    endif
END SUBROUTINE cfftw


SUBROUTINE init_background_density
  ! Setup the confining potential. 
  ! For jellium systems this is the Coulomb attraction to a homogenous positive background charge 
  ! This background charge must have the same density as the atoms being imitated 
  integer :: i,k1,k2,k3,num
  real*8  :: r0(3),x2,d(n_l(1))
  real*8  :: tip_angle,diode_L,diode_R,V_tot1,V_tot2,dR
  real*8  :: atom_rho
  real*8,parameter  :: atom_rs(3)=(/0d0,0d0,1.719826059d0/)  !Agstrom, = 3.25 Bohr  
  real*8,parameter  :: Veps=1d-7
  real*8 :: R_test,dR1,diff1,diff2,diff3
  logical :: quit=.false.
  real*8  :: burn, cyl_L, cyl_R
  ! real*8  ::

  rho_bg=0.d0
  atom_rho=3d0/(4*pi*atom_rs(atom_index)**3)

! sphere
  if(run_type==1) then
    radius=(N_orbitals*2)**(1d0/3)*atom_rs(atom_index)
    write(6,*) 'radius',radius
    r0=(/0.d0,0.d0,0.d0/)
    do i=1,N_L_points
      r0(1)=grid_point(1,i)-center(1,1)
      r0(2)=grid_point(2,i)-center(2,1)
      r0(3)=grid_point(3,i)-center(3,1)
      if(sqrt(r0(1)**2+r0(2)**2+r0(3)**2)<radius) rho_bg(i)=atom_rho
    enddo
!    write(6,*) (4d0/3)*pi*radius**3*atom_rho
  endif

! special diode shape
  if(run_type==2) then
    open(1,file='diode.inp')
    read(1,*) tip_angle
    read(1,*) diode_L
!    read(1,*) diode_R
    close(1)
    !solve for diode radius
    dR=0.1
    R_test=diode_L/3
    DO ii=1,10000000 
      do i=1,1000000000
        diode_R=R_test
        diff2=abs(pi*diode_R**2*(-2*diode_R/tan(tip_angle*pi/360d0)/3+diode_L) - N_orbitals*2/atom_rho)
        diode_R=R_test+dR
        diff3=abs(pi*diode_R**2*(-2*diode_R/tan(tip_angle*pi/360d0)/3+diode_L) - N_orbitals*2/atom_rho)
        diode_R=R_test-dR
        diff1=abs(pi*diode_R**2*(-2*diode_R/tan(tip_angle*pi/360d0)/3+diode_L) - N_orbitals*2/atom_rho)
        if(diff3.lt.diff2) R_test=R_test+dR
        if(diff1.lt.diff2) R_test=R_test-dR
        if(diff1.gt.diff2.and.diff3.gt.diff2) then
          dR=dR*0.1
          exit
        endif
        if(diff1==diff2.and.diff3==diff2) then
          quit=.true.
          exit
        endif
      enddo
      if(quit) exit
    ENDDO
    write(6,*) 'diode radius: ',diode_R
    V_tot1=pi*diode_R**2*(-2*diode_R/tan(tip_angle*pi/360d0)/3+diode_L)
    V_tot2=N_orbitals*2/atom_rho
!    write(6,*) atom_rho
!    write(6,*) V_tot1,V_tot2
    do i=1,N_L_points
      rho_bg(i)=0.d0
      r0(1)=grid_point(1,i)-center(1,1)
      r0(2)=grid_point(2,i)-center(2,1)
      r0(3)=grid_point(3,i)-center(3,1)
      if(check_pos(r0,diode_R,diode_L,tip_angle)) rho_bg(i)=atom_rho
      if(check_pos(r0,diode_R,diode_L,tip_angle))write(400,'(A,3F12.7)') 'C ',grid_point(1,i),grid_point(2,i),grid_point(3,i)
    enddo
  endif

! symmetrical cylinder shape
    if(run_type==3) then
        open(1,file='diode.inp')
        read(1,*) burn
        read(1,*) cyl_L
        close(1)
    ! solve for cyl_R
        cyl_R = sqrt((4d0/3d0)*(N_orbitals*2d0)*(atom_rs(atom_index)**3d0)/cyl_L)
        write(6,*) 'N_orbitals: ', N_orbitals
        write(6,*) 'atom_rs: ', atom_rs
        write(6,*) 'atom_index ', atom_index
        write(6,*) 'cyl_L: ', cyl_L
        write(6,*) 'cylinder radius: ', cyl_R
    ! init background density
        do i=1,N_L_points
            rho_bg(i)=0.d0
            r0(1)=grid_point(1,i)-center(1,1)
            r0(2)=grid_point(2,i)-center(2,1)
            r0(3)=grid_point(3,i)-center(3,1)
            if(check_pos_cyl(r0, cyl_R, cyl_L)) rho_bg(i) = atom_rho
           ! if(check_pos_cyl(r0, cyl_R, cyl_L)) write(400,'(A,3F12.7)') 'C ',grid_point(1,i),grid_point(2,i),grid_point(3,i)
        enddo
    endif

  call save_bov_dat_real(rho_bg,'bg_dens',"./",'density')
  write(6,*) 'sum(rho_bg)*grid_volume', sum(rho_bg)*grid_volume
  density=rho_bg
  call solve_poisson
  V_ext=-VH
 
  d=0.d0
  do i=1,n_l_points
    d(lattice(1,i)+1)=d(lattice(1,i)+1)+v_ext(i)
  end do
  do i=1,n_L(1)
    write(30,*)-0.5d0*grid_step(1)*n_l(1)+i*grid_step(1),d(i)
    write(31,*)-0.5d0*grid_step(1)*n_l(1)+i*grid_step(1),V_ext(Lattice_inv(i,n_L(2)/2,n_L(3)/2))
    write(32,*)-0.5d0*grid_step(1)*n_l(1)+i*grid_step(1),V_ext(Lattice_inv(i,n_L(2)/2,n_L(3)))
  end do
  write(25,*)  
END SUBROUTINE init_background_density

FUNCTION check_pos(r,R0,Ltot,theta)
  ! checks if a point r is within a sharp tipped cylinder
  ! R0, radius of cylinder 
  ! L0, length of cylinder [centered at origin, along z-axiz]
  ! theta, angle of tip [degrees]
  real*8  :: r(3),R0,L0,theta
  real*8  :: R1,L1,rr,Ltot
  logical :: check_pos
  real*8,parameter :: pi=4.d0*atan(1.d0)
    check_pos=.false.
    rr=r(1)**2+r(2)**2
    L1=R0/tan(theta*pi/360d0)
    L0=Ltot-L1
    if(sqrt(rr).le.R0) then
      if(r(3).ge.(-Ltot*0.5d0)) then
        if(r(3).le.(Ltot*0.5d0)) then
          if(r(3).ge.((L0-L1)*0.5d0)) then !inside the tip
            R1=R0-(r(3)-(L0-L1)*0.5d0)*tan(theta*pi/360d0)   
            if(sqrt(rr).le.R1) then   
              check_pos=.true.
            endif
          else !in the body
            check_pos=.true.
          endif
        endif
      endif
    endif
END FUNCTION check_pos

FUNCTION check_pos_cyl(p, R, L_cyl)
    ! checks if a point p is in the symmetrical cylinder
    ! p: point with x,y,z coordinates
    ! R: radius of cylinder
    ! L: length of cylinder along z-axis and centered at origin
    real*8, intent(in)  :: p(3), R, L_cyl
    real*8              :: z1, z2, dL !L_box,x0,y0
    logical             :: check_pos_cyl
    ! calculate constants
!    x0    = 0.5*(N_L(1)*grid_step(1))
!    y0    = 0.5*(N_L(2)*grid_step(2))
!    L_box = N_L(3)*grid_step(3)
!    z1    = 0.5*(L_box - L_cyl)
!    z2    = L_box - z1
    z1 = -0.5*L_cyl
    z2 = -z1
    dL    = sqrt(p(1)**2 + p(2)**2)
    ! check coordinates
    check_pos_cyl = .false.
    if(p(3).ge.z1 .and. p(3).le.z2) then
        if(dL.le.R) then
            check_pos_cyl = .true.
        endif
    endif
END FUNCTION check_pos_cyl

SUBROUTINE laplace_operator
  ! 4th order finite difference representation of the Laplacian operator
  ! sections 2.2 and 7.1 of Computational Nanoscience
  real*8,parameter :: cN0=-205.d0/72.d0,cN1=8.d0/5.d0
  real*8,parameter :: cN2=-1.d0/5.d0,cN3=8.d0/315.d0, cN4=-1.d0/560.d0
  integer          :: i,i1,i2,i3
  real*8           :: K_x,K_y,K_z
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    wf(i1,i2,i3)=Phi(i)
  enddo
  call periodic
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    K_x=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1+1,i2,i3)+wf(i1-1,i2,i3))+& 
             cN2*(wf(i1+2,i2,i3)+wf(i1-2,i2,i3))+& 
             cN3*(wf(i1+3,i2,i3)+wf(i1-3,i2,i3))+& 
             cN4*(wf(i1+4,i2,i3)+wf(i1-4,i2,i3)))/Grid_Step(1)**2
    K_y=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2+1,i3)+wf(i1,i2-1,i3))+& 
             cN2*(wf(i1,i2+2,i3)+wf(i1,i2-2,i3))+& 
             cN3*(wf(i1,i2+3,i3)+wf(i1,i2-3,i3))+& 
             cN4*(wf(i1,i2+4,i3)+wf(i1,i2-4,i3)))/Grid_Step(2)**2
    K_z=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2,i3+1)+wf(i1,i2,i3-1))+& 
             cN2*(wf(i1,i2,i3+2)+wf(i1,i2,i3-2))+& 
             cN3*(wf(i1,i2,i3+3)+wf(i1,i2,i3-3))+& 
             cN4*(wf(i1,i2,i3+4)+wf(i1,i2,i3-4)))/Grid_Step(3)**2
    L_Phi(i)=K_x+K_y+K_z
  enddo
END SUBROUTINE laplace_operator

SUBROUTINE Hamiltonian_wavefn
  ! Calculate H|Phi>
  call laplace_operator
  H_Phi=(V_ext+VH+V_exchange)*Phi-h2m*L_Phi
END SUBROUTINE Hamiltonian_wavefn

SUBROUTINE conjugate_gradient
  ! Conjugate gradient minimization to diagonalize the Hamiltonian
  ! section 18.2 of Computational Nanoscience
  integer :: orbital,iteration,i
  real*8  :: Phi_H_Phi,gamma,beta_Phi,beta_beta,beta_H_Phi
  real*8  :: Phi0_Phi0,Phi0_H_Phi0,delta,A,B,C,omega,overlap
  real*8, allocatable :: alpha(:),beta(:),Phi0(:),H_Phi0(:)
  allocate(alpha(N_L_points),beta(N_L_points),Phi0(N_L_points),H_Phi0(N_L_points))
  do orbital=1,N_orbitals
	  Phi=Psi(:,orbital)
	  Phi0=Phi
	  call Hamiltonian_wavefn
	  H_Phi0=H_Phi
	  Phi0_H_Phi0=sum(Phi*H_Phi)*Grid_Volume
	  Phi0_Phi0=1
	  delta=Phi0_H_Phi0
    do iteration=1,N_iteration
		  alpha=2*(H_Phi0-delta*Phi0)
		  do i=1,orbital-1
			  alpha=alpha-Psi(:,i)*sum(Psi(:,i)*alpha)*Grid_Volume
		  enddo
		  gamma=sum(alpha*alpha)*Grid_Volume
		  if(iteration==1) beta=-alpha
		  if(iteration>1) beta=-alpha+gamma/overlap*beta
		  overlap=gamma
		  beta_Phi=sum(Phi0*beta)*Grid_Volume
		  beta_beta=sum(beta*beta)*Grid_Volume
		  beta_H_Phi=sum(H_Phi0*beta)*Grid_Volume

		  Phi=beta
		  call Hamiltonian_wavefn
		  alpha=H_Phi
		  Phi_H_Phi=sum(Phi*H_Phi)*Grid_Volume
		  A = Phi_H_Phi*beta_Phi  - beta_H_Phi*beta_beta
		  B = Phi_H_Phi*Phi0_Phi0   - Phi0_H_Phi0*beta_beta
		  C = beta_H_Phi*Phi0_Phi0 - Phi0_H_Phi0*beta_Phi
		  omega=(-B+sqrt(B*B-4*A*C))/(2*A)
		  Phi0   = Phi0   +omega*Phi
		  H_Phi0 = H_Phi0 +omega*H_Phi
		  Phi0_Phi0   =sum(Phi0*Phi0)*Grid_Volume
		  Phi0_H_Phi0 =sum(Phi0*H_Phi0)*Grid_Volume
		  delta=Phi0_H_Phi0/Phi0_Phi0
	  enddo
	  Phi0_Phi0=sum(Phi0*Phi0)*Grid_Volume
	  Psi(:,orbital)=Phi0/sqrt(Phi0_Phi0)
  enddo
  deallocate(alpha,beta,Phi0,H_Phi0)
END SUBROUTINE conjugate_gradient


SUBROUTINE orthogonalization
  ! Gram-Schmidt orthogonalization
  ! http://mathworld.wolfram.com/Gram-SchmidtOrthonormalization.html
  integer    :: orbital,i
  real*8     :: s
  complex*16 :: overlap
  do orbital=1,N_orbitals
    do i=1,orbital-1
      overlap=sum(Psi(:,i)*Psi(:,orbital))*Grid_Volume
      Psi(:,orbital)=Psi(:,orbital)-overlap*Psi(:,i)
    enddo
    s=sum(abs(Psi(:,orbital))**2)*Grid_Volume
    Psi(:,orbital)=Psi(:,orbital)/sqrt(s)
  enddo
END SUBROUTINE orthogonalization


SUBROUTINE orthogonalization_c
  ! Gram-Schmidt orthogonalization for complex functions in momentum space
  ! http://mathworld.wolfram.com/Gram-SchmidtOrthonormalization.html
  integer    :: orbital,i
  real*8     :: s
  complex*16 :: overlap
  do orbital=1,N_orbitals
    do i=1,orbital-1
      overlap=sum(conjg(Psi_g(:,i))*Psi_g(:,orbital))
      Psi_g(:,orbital)=Psi_g(:,orbital)-overlap*Psi_g(:,i)
    enddo
    s=sum(abs(Psi_g(:,orbital))**2)
    Psi_g(:,orbital)=Psi_g(:,orbital)/sqrt(s)
  enddo
END SUBROUTINE orthogonalization_c


SUBROUTINE Hamiltonian_density
  ! Calculate the density dependent part of the Kohn-Sham Hamiltonian
  ! Hartree potential 
  call solve_poisson
  ! Exchange-Correlation potential
  call Exchange_Correlation_potential
END SUBROUTINE Hamiltonian_density


SUBROUTINE Hartree_Energy
  ! Calculate the Hartree-Energy
  ! http://www.fisica.uniud.it/~giannozz/Corsi/MetNum/LectureNotes/metnum-cap2.pdf
  ! Equation 2.12
  E_Hartree=-0.5d0*dot_product(density,VH)*grid_volume
END SUBROUTINE Hartree_Energy


SUBROUTINE total_energy(iteration)
  ! Calculate the total energy
  ! http://www.home.uni-osnabrueck.de/apostnik/Lectures/DFT-3.pdf
  ! Equation 3.17
  integer :: i,orbital,iteration
  real*8  :: E_sp,E_total                 
  E_sp=0.d0
  call Hamiltonian_density
    do orbital=1,N_orbitals
      Phi=Psi(:,orbital)
      call Hamiltonian_wavefn
      Sp_Energy(orbital)=sum(Phi*H_Phi)*Grid_Volume
      E_sp=E_sp+2.0*Sp_energy(orbital)
    enddo
  call Exchange_Correlation_energy
  call Hartree_Energy
  write(6,*)'E_sp    ',E_sp
  write(6,*)'E_ex    ',E_exchange
  write(6,*)'E_H     ',E_hartree
  E_total=E_sp+E_Exchange+E_Hartree
  write(6,*)'E_total ',E_total
  write(100,*) E_total
!  if(iteration>0) write(2,'(i10,4ES16.8)')iteration,E_sp,E_exchange,E_hartree,E_total
END SUBROUTINE total_energy


SUBROUTINE energy(e)
  ! Calculate the total energy from a wave function in momentum space
  ! http://www.home.uni-osnabrueck.de/apostnik/Lectures/DFT-3.pdf
  ! Equation 3.17
  ! Equation 13.68 in Compuational Nanoscience (add in the vector potential though)
  implicit none
  integer               :: i,j
  real*8                :: p,t,e
  call Hamiltonian_density
  e=0.d0
  do j=1,N_orbitals
    Phi_g(:)=Psi_g(:,j)
    do i=1,N_L_points
      t=h2m*((g_vectors(1,i)+A_vec(1)/(hbar))**2+(g_vectors(2,i)+A_vec(2)/(hbar))**2+(g_vectors(3,i)+A_vec(3)/(hbar))**2)
!      e=e+t*phi_g(i)*conjg(phi_g(i))*cell_volume
      e=e+2.d0*t*phi_g(i)*conjg(phi_g(i))
    end do
    call fft_to_r_space
    do i=1,N_L_points
      p=V_ext(i)+VH(i)+V_exchange(i)
      e=e+2.d0*p*phi_r(i)*conjg(phi_r(i))*grid_volume
    end do
  end do
  call Exchange_Correlation_energy
  call Hartree_Energy
  write(6,*)'E_sp    ',e
  write(6,*)'E_ex    ',E_exchange
  write(6,*)'E_H     ',E_hartree
  e=e+E_Exchange+E_Hartree
  write(6,*)'E_total ',e
  write(100,*) e
END SUBROUTINE energy


SUBROUTINE calculate_density(c1,c2)
  ! Calculate the density using linear mixing (mixing helps settle into the ground state)
  ! http://www.fisica.uniud.it/~giannozz/Corsi/MetNum/LectureNotes/metnum-cap2.pdf
  ! Equation 2.5
  integer :: orbital,i,k,grid
  real*8  :: c1,c2
  Density_old=Density; Density=0.d0
	  do orbital=1,N_orbitals
	  Phi(:)=Psi(:,orbital)
	    do i=1,N_L_points
	      Density(i)=Density(i)+2.d0*Phi(i)*Phi(i)
	    enddo
	  enddo
  Density=c1*Density_old+c2*Density
END SUBROUTINE Calculate_Density


SUBROUTINE calculate_density_c
implicit none
  ! Calculate the real-space density from a complex wave function in momentum-space
  ! http://www.fisica.uniud.it/~giannozz/Corsi/MetNum/LectureNotes/metnum-cap2.pdf
  ! Equation 2.5
  integer :: orbital,i
  Density=0.d0
  do orbital=1,N_orbitals
    Phi_g(:)=Psi_g(:,orbital)
    call fft_to_r_space
    do i=1,N_L_points
      Density(i)=Density(i)+2.d0*Conjg(Phi_r(i))*Phi_r(i)
    enddo
  enddo
END SUBROUTINE Calculate_Density_c


SUBROUTINE calculate_current
	! Calculate the real-space prob current using wave function in momentum-space
  ! https://en.wikipedia.org/wiki/Probability_current
	integer    :: orbital,i,k,grid,m,ig
  real*8     :: f
  complex*16,allocatable :: wlap_p(:,:)
  allocate(wlap_p(3,N_L_points))
  current_density=0.d0
  do orbital=1,n_orbitals
    !Compute the gradient of an orbital
    do m=1,3
      Phi_g(:)=zi*g_vectors(m,:)*Psi_g(:,orbital)
      call fft_to_r_space
      wlap_p(m,:)=Phi_r(:)
    enddo
! flip the charge
    f=2.d0*hbar/ElectronMass  !electron has negative charge, i.e. -1
    Phi_g(:)=Psi_g(:,orbital)
    call fft_to_r_space
    do grid=1,N_L_points
      do m=1,3
        current_density(m,grid)=current_density(m,grid)+f*imag(conjg(Phi_r(grid))*wlap_p(m,grid))
        current_density(m,grid)=current_density(m,grid)+f*conjg(Phi_r(grid))*Phi_r(grid)*A_vec(m)/(hbar)
      enddo
    enddo
  enddo
  deallocate(wlap_p)
END SUBROUTINE calculate_current


SUBROUTINE calculate_dipole_moment(rho,x)
  !calculate dipole moment (z-direction)
  integer :: j
  real*8  :: rho(*)
  double precision :: x
  x=0.0d0
  do j=1,N_l_points
    x=x+rho(j)*Grid_Point(3,j)
  enddo
  x=x*grid_volume
END SUBROUTINE calculate_dipole_moment


SUBROUTINE calculate_flux(prob_current,value)
    ! current flux at the z-axis boundary
    integer  :: k1,k2,k3
    real*8   :: prob_current(3,*)
    double precision  :: value
    value=0.d0
    k3=1
    do k2=1,N_L(2)
      do k1=1,N_L(1)
        value=value+prob_current(3,Lattice_inv(k1,k2,k3))*grid_step(1)*grid_step(2)
      enddo
    enddo
END SUBROUTINE calculate_flux


SUBROUTINE check_tunneling(potential,energy,Tunnel_prob)
  !calculate an approximate tunneling probabilty
  !assume only two turning points
  integer :: i
  real*8  :: potential(*),energy,Tunnel_prob
  real*8  :: turning_point(2),d,pot1d(N_L(3))
  real*8  :: k,kappa,k2_kappa2_sum_2,k_kappa_prod_2
  logical :: below=.true.,below_prev=.true.
  !find turning points
  turning_point=0.d0
  do k3=1,N_L(3)
    pot1d(k3)=potential(Lattice_inv(N_L(1)/2,N_L(2)/2,k3))
  enddo
  if(energy.lt.pot1d(1)) below_prev=.true.
  do k3=1,N_L(3)
    if(energy.lt.pot1d(k3)) below=.true.
    if(energy.gt.pot1d(k3)) below=.false.
    if(.not.below_prev.and.below) then
      if(turning_point(1).ne.(0.d0)) write(6,*) 'more than one turning point'
      turning_point(1)=grid_point(3,ind)
    endif
    if(.not.below.and.below_prev) then
      if(turning_point(2).ne.(0.d0)) write(6,*) 'more than one turning point'
      turning_point(2)=grid_point(3,ind)
    endif
    below_prev=below
  enddo
  d=turning_point(1)-turning_point(2)
  if(d.lt.(0.d0)) d=turning_point(1)+grid_point(3,N_L(3))+grid_step(3)-turning_point(2)
  V0=maxval(pot1d)
  kappa=sqrt(2d0*ElectronMass*(V0-energy))/hbar
  k=sqrt(2d0*ElectronMass*energy)/hbar
  k2_kappa2_sum_2=(k**2+kappa**2)**2
  k_kappa_prod_2=(2d0*k*kappa)**2
  Tunnel_prob=k_kappa_prod_2/(k2_kappa2_sum_2*sinh(kappa*d)+k_kappa_prod_2)
END SUBROUTINE check_tunneling


FUNCTION Keldysh_parameter(potential,energy)
  !calculate the Keldysh parameter
  !http://www.iuac.res.in/atmol/~safvan/safvan_thesis/node32.html
  integer :: i
  real*8  :: potential(*),energy, Keldysh_parameter
  real*8  :: pot1d(N_L(3))
  do k3=1,N_L(3)
    pot1d(k3)=potential(Lattice_inv(N_L(1)/2,N_L(2)/2,k3))
  enddo
  V0=maxval(pot1d)
  Keldysh_parameter=omega*sqrt(2d0*(V0-energy)*electronmass)/E_0
END FUNCTION Keldysh_parameter


SUBROUTINE setup_g_vectors
  ! initialize g-space (momentum-space) vectors 
  ! equations 2.35 and 2.36 in Computational Nanoscience
  implicit none
  integer                              :: i1,i2,i3,k1,k2,k3,k,n1,n2,n3
  double precision                     :: x,y,z,box(3)
  n1=N_L(1)
  n2=N_L(2)
  n3=N_L(3)
  box(1)=N_L(1)*grid_step(1)
  box(2)=N_L(2)*grid_step(2)
  box(3)=N_L(3)*grid_step(3)
  allocate(G_vectors(3,N_L_points),g_length(N_L_points)) 
! G_vectors
  k=0
  do i1=0,N1-1
    if(i1.le.N1/2) then
      x=2.d0*pi*i1/box(1)
    else
      x=2.d0*pi*(i1-N1)/box(1)
    endif          
    do i2=0,N2-1
      if(i2.le.N2/2) then
        y=2.d0*pi*i2/box(2)
      else
        y=2.d0*pi*(i2-N2)/box(2)
      endif          
      do i3=0,N3-1
        if(i3.le.N3/2) then
          z=2.d0*pi*i3/box(3)
        else
          z=2.d0*pi*(i3-N3)/box(3)
        endif          
        k=k+1
        g_vectors(1,k)=x
        g_vectors(2,k)=y
        g_vectors(3,k)=z
        g_length(k)=x**2+y**2+z**2
      end do
    end do
  end do
END SUBROUTINE setup_g_vectors

SUBROUTINE fft_to_g_space
  ! fourier transform a wave fucntion in real-space (phi_r) to momentum-space (phi_g)
  implicit none
  integer                                        :: i1,i2,i3,i,mode
  complex*16, dimension(:,:,:),allocatable       :: rho
  allocate(rho(N_L(1),N_L(2),N_L(3)))   
  do i=1,N_L_points
    i1=Lattice(1,i)+1
    i2=Lattice(2,i)+1
    i3=Lattice(3,i)+1
    rho(i1,i2,i3)=phi_r(i)
  end do
  ! fourier transformation 
  mode=1
  call cfftw(rho,N_L(1),N_L(2),N_L(3),mode)
   do i=1,N_L_points
     i1=Lattice(1,i)+1
     i2=Lattice(2,i)+1
     i3=Lattice(3,i)+1
     phi_g(i)=rho(i1,i2,i3)
   end do
  Phi_g=Phi_g*(Cell_Volume)**(0.5d0*mode)
  deallocate(rho)
END SUBROUTINE fft_to_g_space

SUBROUTINE fft_to_r_space
  ! fourier transform a wave fucntion in momentum-space (phi_g) to real-space (phi_r)
  implicit none
  integer                                        :: i1,i2,i3,i,mode
  complex*16, dimension(:,:,:),allocatable       :: rho
  allocate(rho(N_L(1),N_L(2),N_L(3)))  
  do i=1,N_L_points
    i1=Lattice(1,i)+1
    i2=Lattice(2,i)+1
    i3=Lattice(3,i)+1
    rho(i1,i2,i3)=phi_g(i)
  end do
  !  transformation to real space 
  mode=-1
  call cfftw(rho,N_L(1),N_L(2),N_L(3),mode)
   do i=1,N_L_points
     i1=Lattice(1,i)+1
     i2=Lattice(2,i)+1
     i3=Lattice(3,i)+1
     Phi_r(i)=rho(i1,i2,i3)
   end do
  Phi_r=Phi_r*(Cell_Volume)**(0.5d0*mode)
  deallocate(rho)
END SUBROUTINE fft_to_r_space


SUBROUTINE imaginary_time
  integer          :: i,k
  double precision :: e,d(N_L(3))
  complex*16 :: fac
  allocate(Psi_c(N_L_points,n_orbitals),Psi_g(N_L_points,n_orbitals), &
&   Phi_g(N_L_Points),Phi_r(N_L_Points),grad_w(N_L_points,3), Psi_so(N_L_points,n_orbitals))
  Psi_c=Psi
  ! Fourier Transform to g-space (momentum-space)
  do i=1,n_orbitals
    phi_r(:)=psi_c(:,i)
    call fft_to_g_space
    Psi_g(:,i)=phi_g(:)
  end do
  Psi_so=Psi_g
  time_step=-zi*0.002d0
  do k=1,N_scf_it_iter
    write(6,*)k
    call time_propagator
    call orthogonalization_c
    call calculate_density_c
    call energy(e)
  enddo  
END SUBROUTINE imaginary_time


SUBROUTINE time_propagator
  ! propagate the wave function from psi(t) to psi(t+dt) via Volkov propagation
  ! https://journals.aps.org/pra/abstract/10.1103/PhysRevA.95.013414
  ! equation 19
!!! small corrections still needed. See notes
  implicit none
  integer                :: i,j
  real*8                 :: p,t
  complex*16,allocatable :: tv(:), tv0(:)
  allocate(tv(N_L_points))
  do i=1,N_L_points
    t=h2m*((g_vectors(1,i)+A_vec(1)/(hbar))**2+(g_vectors(2,i)+A_vec(2)/(hbar))**2+(g_vectors(3,i)+A_vec(3)/(hbar))**2)
    tv(i)=exp(-0.5d0*zi*time_step*t/hbar)
  end do
  do j=1,N_orbitals
    Phi_g(:)=Psi_so(:,j)
    !half of kinetic applied in momentum-space
    phi_g=phi_g*tv
    !potential exponential applied in real-space
    call fft_to_r_space
    do i=1,N_L_points
      p=V_ext(i)+VH(i)+V_exchange(i)
      phi_r(i)=phi_r(i)*exp(-zi*time_step*p/hbar)
    end do
    call fft_to_g_space
    !half of kinetic applied in momentum-space
    phi_g=phi_g*tv
    Psi_so(:,j)=Phi_g(:)
    Psi_g(:,j)=Phi_g(:)
  end do
  deallocate(tv)
END SUBROUTINE time_propagator


SUBROUTINE time
! Propagate the wave function via TDDFT
! Options available for outputing the time-dependent dipole moment after a delta kick
! uses previously found Psi_c
implicit none
  integer  :: it,k,i,j,N_time,itt,k1,k2,k3
  double precision   :: t,e,dp,ds,flux,z,tun_prob
  real*8             :: f,su,dt
  real*8,allocatable :: V_tot(:),J_z(:)
  allocate(V_tot(N_L_points),J_z(N_L_points))
  call init_files
  ! read in laser parameters
  open(1,file='td_laser.inp')
    read(1,*) dt
    time_step=dt
    read(1,*) N_time
    read(1,*) E_0
    read(1,*) laser_wavelength ![nanometer] --> frequency[rad/fs]
    omega=2.d0*pi*sqrt(c_squared)/(10.d0*laser_wavelength)
    tt=8.d0/omega
    read(1,*) laser_pulse_shift1
    read(1,*) laser_pulse_width1
    read(1,*) dielectric
    read(1,*) kick_strength
    read(1,*) phase_shift  
    f=kick_strength
  close(1)
  IF(dielectric) then
    ! delta kick for optical response (z-direction). calculation of dielectric function
    call calculate_density(0.d0,1.d0)
    do k=1,N_orbitals
      do i=1,N_l_points
        z=Grid_Point(3,i)
        Psi_c(i,k)=Psi(i,k)*exp(zi*f*z)
      end do
    end do
    write(dp_file,*) 10000,real(time_step),f
  ENDIF
  !initial dipole moment
  call calculate_dipole_moment(density,ds)

  Psi_so=psi_g
  call calculate_density_c
  call calculate_current
  call calculate_dipole_moment(density,ds)
  call calculate_flux(current_density,flux)
  call energy(e)
  write(energy_file,*)0d0,e
  write(dp_file,*)0d0,0.d0
  write(flux_file,*) 0d0,flux
  write(6,*)e,sum(density)*grid_volume
  V_tot=V_ext+VH+V_exchange
  call check_tunneling(V_tot,e,tun_prob)
  write(500,*)0.d0,tun_prob
  write(6,*) 'Keldysh Parameter: ',Keldysh_parameter(V_tot,sp_energy(N_orbitals))
  write(600,*) Keldysh_parameter(V_tot,sp_energy(N_orbitals))
  J_z(:)=current_density(3,:)
!  call save_bov_dat_real(density,'dens',"./density/",'density',0,0.d0)
!  call save_1D_proj_real(density,'dens1D',"./density1Dp/",0,0.d0)
!  call save_cut_real(density,'dens2D',"./density2D/",0,0.d0) !takes a lot of memory
!  call save_1D_cut_real(density,'dens1D',"./density1D/",0,0.d0)
!  call save_1D_cut_real(V_tot,'pot1D',"./potential1D/",0,0.d0)
!  call save_1D_cut_real(J_z,'curr1D',"./current1D/",0,0.d0)
!  call save_1D_proj_real(J_z,'flux',"./flux/",0,0.d0)

  itt=0
  A_vec=0.d0
  do it=0,N_time
    t=real(time_step)*it
    write(6,*)
    write(6,*) t
    ! define vector potential in the z-direction
    if(.not.dielectric) then
      if(t<3.d0*tt) then
        E_field=E_0*sin(pi*t/(6.d0*tt))*sin(omega*t+phase_shift*2.d0*pi)
      else
        E_field=E_0*sin(omega*t+phase_shift*2.d0*pi)
      endif
      if(laser_pulse_width1.gt.(0.d0)) E_field=E_field*exp(-((t-laser_pulse_shift1)/laser_pulse_width1)**2)
      A_vec(3)=A_vec(3)-E_field*real(time_step)
    endif
    write(e_field_file,*)t,E_field
    write(vec_pot_file,*)t,a_vec(3)

    call time_propagator

    call calculate_density_c
    call calculate_current
    call calculate_dipole_moment(density,dp)
    call calculate_flux(current_density,flux)
    call energy(e)
    write(energy_file,*)t,e
    write(dp_file,*)t,dp-ds
    write(flux_file,*) t,flux
    write(6,*)e,sum(density)*grid_volume
    V_tot=V_ext+VH+V_exchange
    call check_tunneling(V_tot,e,tun_prob)
    write(500,*)t,tun_prob

    if(mod(it+1,1)==0) then
      itt=itt+1
      J_z(:)=current_density(3,:)
!      call save_bov_dat_real(density,'dens',"./density/",'density',itt,t)
!      call save_1D_proj_real(density,'dens1D',"./density1Dp/",itt,t)
!      call save_cut_real(density,'dens2D',"./density2D/",itt,t)
!      call save_1D_cut_real(density,'dens1D',"./density1D/",itt,t)
!      call save_1D_cut_real(V_tot,'pot1D',"./potential1D/",itt,t)
!      call save_1D_cut_real(J_z,'curr1D',"./current1D/",itt,t)
!      call save_1D_proj_real(J_z,'flux',"./flux/",itt,t)
    endif
  end do

  close(dp_file)
  close(energy_file)
  close(e_field_file) 
  close(vec_pot_file) 
  close(flux_file)
  
END SUBROUTINE time


subroutine save_bov_dat_real(real8grid_array,fname_prefix,folder,quantity_name,findex,ftime)
!This is a general subroutine to save scalar real data defined on a 
!grid (given in linear array real8grid_array) into a file. The data is 
!saved in BOV format, i.e. two files are created: <fname_prefix>.bov and 
!<fname_prefix>.dat. Such data can be then visualized with VisIt software. 
  real(8)                     :: real8grid_array(*)
  character(len=*)            :: fname_prefix,folder
  character(len=*), optional  :: quantity_name
  integer, optional           :: findex ! index (may range from 0 to 99999) added
                                        ! between <fname_prefix> and ".bov" and ".dat"
  real(8), optional           :: ftime  ! time which the BOV file corresponds to
  !Local variables
  integer                     :: k1,k2,k3,ind,pl
  character(255)              :: fname_bov,fname_dat,quantity_to_write
  real(8)                     :: time_to_write
  character(5)                :: ch5
  integer,parameter           :: gen_purpose_bovfile=3, gen_purpose_datfile=4
  real(4),allocatable         :: rearranged_arr_real4(:)
  if (present(quantity_name)) then
    quantity_to_write=quantity_name
  else
    quantity_to_write='myquantity'
  endif 
  pl=len_trim(fname_prefix)
  if (present(findex)) then
    write(ch5,'(i5.5)') findex
    write(fname_bov,'(a,a,a)') fname_prefix(1:pl),ch5,'.bov'
    write(fname_dat,'(a,a,a)') fname_prefix(1:pl),ch5,'.dat'
  else
    write(fname_bov,'(a,a)') fname_prefix(1:pl),'.bov'
    write(fname_dat,'(a,a)') fname_prefix(1:pl),'.dat' 
  endif 
  if (present(ftime)) then
    time_to_write=ftime
  else
    time_to_write=0.0d0
  endif 
  !Store actual data in binary form (3D array needs to have its dimension order changed, 
  !linearized, and converted to single precision)
  allocate(rearranged_arr_real4(N_L(1)*N_L(2)*N_L(3)))
  ind=0
  do k3=1,N_L(3)
    do k2=1,N_L(2)
      do k1=1,N_L(1)
        ind=ind+1
        if (Lattice_inv(k1,k2,k3)>0) then
          rearranged_arr_real4(ind)=real8grid_array(Lattice_inv(k1,k2,k3))
        else
          rearranged_arr_real4(ind)=0.0
        endif
      enddo
    enddo
  enddo
  open(gen_purpose_datfile,file=trim(adjustl(folder))//trim(adjustl(fname_dat)),form='unformatted',status='replace')
  write(gen_purpose_datfile) rearranged_arr_real4
  close(gen_purpose_datfile)
  !Now create a header file corresponding to the above datafile. The header file
  !is in ascii text form
  open(gen_purpose_bovfile,file=trim(adjustl(folder))//trim(adjustl(fname_bov)),status='replace')
  write(gen_purpose_bovfile,'(1x,a,1x,es23.15)') 'TIME:',time_to_write
  write(gen_purpose_bovfile,'(1x,a,1x,a)') 'DATA_FILE:',trim(adjustl(fname_dat))
  write(gen_purpose_bovfile,'(1x,a,3(1x,i6))') 'DATA_SIZE:',N_L(1),N_L(2),N_L(3)
  write(gen_purpose_bovfile,'(1x,a)') 'DATA_FORMAT:  FLOAT'
  write(gen_purpose_bovfile,'(1x,a,1x,a)') 'VARIABLE:',trim(adjustl(quantity_to_write))
  write(gen_purpose_bovfile,'(1x,a)') 'DATA_ENDIAN:  LITTLE'
  write(gen_purpose_bovfile,'(1x,a)') 'CENTERING:  zonal'
  write(gen_purpose_bovfile,'(1x,a,3(1x,es23.15))') 'BRICK_ORIGIN:', &
      grid_point(1,1),grid_point(2,1),grid_point(3,1)
  write(gen_purpose_bovfile,'(1x,a,3(1x,es23.15))') 'BRICK_SIZE:',(N_L(1)-1)*grid_step(1),(N_L(2)-1)*grid_step(2),(N_L(3)-1)*grid_step(3)
  write(gen_purpose_bovfile,'(1x,a)') 'BYTE_OFFSET: 4'
  close(gen_purpose_bovfile)
  deallocate(rearranged_arr_real4)
end subroutine save_bov_dat_real


subroutine save_1D_proj_real(real8grid_array,fname_prefix,folder,findex,ftime)
  !print the values of a 3D real function only along the z-axis. I.e. values along (0,0,z)
  real(8)                     :: real8grid_array(*)
  character(len=*)            :: fname_prefix,folder
  integer, optional           :: findex ! index (may range from 0 to 99999) added
                                        ! between <fname_prefix> and ".bov" and ".dat"
  real(8), optional           :: ftime  ! time which the BOV file corresponds to
  !Local variables
  integer                     :: k1,k2,k3,ind,pl
  character(255)              :: fname_bov,fname_dat,quantity_to_write
  real(8)                     :: time_to_write
  character(5)                :: ch5
  integer,parameter           :: gen_purpose_datfile=4
  double precision            :: linedata(N_L(3))
  pl=len_trim(fname_prefix)
  if (present(findex)) then
    write(ch5,'(i5.5)') findex
    write(fname_dat,'(a,a,a)') fname_prefix(1:pl),ch5,'.dat'
  else
    write(fname_dat,'(a,a)') fname_prefix(1:pl),'.dat' 
  endif 
  if (present(ftime)) then
    time_to_write=ftime
  else
    time_to_write=0.0d0
  endif 
  linedata=0.d0
  do i=1,n_l_points
    linedata(lattice(3,i)+1)=linedata(lattice(3,i)+1)+real8grid_array(i)*grid_step(1)*grid_step(2)
  end do
  open(gen_purpose_datfile,file=trim(adjustl(folder))//trim(adjustl(fname_dat)),status='replace')
  do i=1,n_L(3)
    write(gen_purpose_datfile,*)(i-1)*grid_step(3),linedata(i)
  end do
  close(gen_purpose_datfile)
end subroutine save_1D_proj_real


subroutine save_1D_cut_real(real8grid_array,fname_prefix,folder,findex,ftime)
  !print the values of a 3D real function projected along the z-axis. 
  !I.e. every value along z is the integrated values along the x-y plane of the 3D function at that point in z
  real(8)                     :: real8grid_array(*)
  character(len=*)            :: fname_prefix, folder
  integer, optional           :: findex ! index (may range from 0 to 99999) added
                                        ! between <fname_prefix> and ".bov" and ".dat"
  real(8), optional           :: ftime  ! time which the BOV file corresponds to
  !Local variables
  integer                     :: k1,k2,k3,ind,pl
  character(255)              :: fname_bov,fname_dat,quantity_to_write
  real(8)                     :: time_to_write
  character(5)                :: ch5
  integer,parameter           :: gen_purpose_datfile=4
  double precision            :: linedata(N_L(3))
  pl=len_trim(fname_prefix)
  if (present(findex)) then
    write(ch5,'(i5.5)') findex
    write(fname_dat,'(a,a,a)') fname_prefix(1:pl),ch5,'.dat'
  else
    write(fname_dat,'(a,a)') fname_prefix(1:pl),'.dat' 
  endif 
  if (present(ftime)) then
    time_to_write=ftime
  else
    time_to_write=0.0d0
  endif 
  linedata=0.d0
  do k3=1,n_l(3)
    k1=N_L(1)/2
    k2=N_L(2)/2
    ind=Lattice_inv(k1,k2,k3)
    linedata(k3)=real8grid_array(ind)
  end do
  open(gen_purpose_datfile,file=trim(adjustl(folder))//trim(adjustl(fname_dat)),status='replace')
  do i=1,n_L(3)
    write(gen_purpose_datfile,*)(i-1)*grid_step(3),linedata(i)
  end do
  close(gen_purpose_datfile)
end subroutine save_1D_cut_real


subroutine save_cut_real(real8grid_array,fname_prefix,folder,findex,ftime)
  real(8)                     :: real8grid_array(*)
  character(len=*)            :: fname_prefix, folder
  integer, optional           :: findex ! index (may range from 0 to 99999) added
                                        ! between <fname_prefix> and ".bov" and ".dat"
  real(8), optional           :: ftime  ! time which the BOV file corresponds to
  !Local variables
  integer                     :: k1,k2,k3,ind,pl
  character(255)              :: fname_bov,fname_dat,quantity_to_write
  real(8)                     :: time_to_write
  character(5)                :: ch5
  integer,parameter           :: gen_purpose_datfile=4
  double precision            :: data2D(N_L(2),N_L(3))
  pl=len_trim(fname_prefix)
  if (present(findex)) then
    write(ch5,'(i5.5)') findex
    write(fname_dat,'(a,a,a)') fname_prefix(1:pl),ch5,'.dat'
  else
    write(fname_dat,'(a,a)') fname_prefix(1:pl),'.dat' 
  endif 
  if (present(ftime)) then
    time_to_write=ftime
  else
    time_to_write=0.0d0
  endif 
  data2D=0.d0
  do k3=1,n_l(3)
    do k2=1,n_l(2)
      k1=N_L(1)/2
      ind=Lattice_inv(k1,k2,k3)
      data2D(k2,k3)=real8grid_array(ind)
    end do
  end do
  open(gen_purpose_datfile,file=trim(adjustl(folder))//trim(adjustl(fname_dat)),status='replace')
  do k3=1,n_L(3)
    do k2=1,n_l(2)
      write(gen_purpose_datfile,*)(k3-1)*grid_step(3),(k2-1)*grid_step(2),data2D(k2,k3)
    enddo
  end do
  close(gen_purpose_datfile)
end subroutine save_cut_real


subroutine init_files
  call system("mkdir density")
  call system("rm ./density/dens*")
  call system("mkdir potential1D")
  call system("rm ./potential1D/pot1D*")
  call system("mkdir density1D")
  call system("rm ./density1D/dens1D*")
  call system("mkdir density1Dp")
  call system("rm ./density1Dp/dens1D*")
  call system("mkdir current1D")
  call system("rm ./current1D/curr1D*")
  call system("mkdir flux")
  call system("rm ./flux/flux*")
  call system("mkdir density2D")
  call system("rm ./density2D/dens2D*")
  open(energy_file,file='energy.dat')
  open(e_field_file,file='electric_field_z.dat') 
  open(vec_pot_file,file='vector_potential_z.dat') 
  open(flux_file,file='flux.dat')
  open(dp_file,file='dp.dat')
end subroutine init_files



END MODULE quantum




  USE quantum
  implicit none
  integer            :: i,k,i1,i2,i3,num_args,kk,orbital,slice1,slice2,slice3
  real*8             :: x,y,z,x0,y0,z0,r,the_sum
  double precision,allocatable   :: d(:)
  integer            :: N_atoms
  integer,parameter  :: val_elec(3)=(/1,2,1/)
  character*255      :: shape_type
  
  open(1,file='dft.inp')
    read(1,*) atom_index,N_atoms
    if(mod(N_atoms*val_elec(atom_index),2).ne.0) then
      write(6,*) 'must use even number of electrons'
      write(6,*) 'current number: ',N_atoms*val_elec(atom_index)
      write(6,*) 'valence electrons per atom: ',val_elec(atom_index)
      stop
    endif
    N_orbitals=N_atoms*val_elec(atom_index)/2
    read(1,*) shape_type
    run_type=0
    if(trim(adjustl(shape_type))=="sphere") run_type=1
    if(trim(adjustl(shape_type))=="diode") run_type=2
    if(trim(adjustl(shape_type))=="cylinder") run_type=3
    if(run_type==0) then
      write(6,*) 'unknown run_type'
      stop
    endif
    read(1,*) (N_L(i),i=1,3)
    read(1,*) (grid_step(i),i=1,3)
    read(1,*) N_scf_iter
    read(1,*) N_scf_it_iter
  close(1)
  center(:,1)=grid_step*N_L(:)/2

  
  N_L_points=product(N_L)
  grid_volume=product(grid_step)
  cell_volume=grid_volume*N_L_points
  
  allocate(d(N_L(1)))

  allocate(sp_energy(N_orbitals),Psi(N_L_points,N_orbitals),rho_bg(N_L_points))
  allocate(Lattice(3,N_L_points),Lattice_inv(N_L(1),N_L(2),N_L(3)),grid_point(3,N_L_Points))
  allocate(V_POT(N_L_points))
  allocate(density(N_L_points),density_old(N_L_points))
  allocate(phi(N_L_points),L_phi(N_L_points),VH(N_L_points),V_ext(N_L_points))
  allocate(V_exchange(N_L_points),H_Phi(N_L_Points))
  allocate(wf(-N_d:N_L(1)+N_d,-N_d:N_L(2)+N_d,-N_d:N_L(3)+N_d))  
  allocate(current_density(3,N_L_points))
    
  ! Setup the lattice and initial guess for wavefunctions  
  call init_lattice
  call setup_g_vectors
  do k=1,N_orbitals
    call random_number(x0)
    call random_number(y0)
    call random_number(z0)
    do i=1,N_L_points
      x=grid_point(1,i)-2*(0.5d0-x0)-center(1,1)
      y=grid_point(2,i)-2*(0.5d0-y0)-center(2,1)
      z=grid_point(3,i)-2*(0.5d0-z0)-center(3,1)
      Psi(i,k)=exp(-0.5d0*(x**2+y**2+z**2))
    enddo
  enddo
  
  ! Setup the background density
  call init_background_density
  
  ! Orthogonalize the initial wavefunctions, and use them to calculate the initial density and energy
  call orthogonalization
  call calculate_density(0.d0,1.d0)

  write(6,*)'?',sum(density)*grid_volume
  call total_energy(0)

  ! Use the conjugate gradient method to diagonalize the Hamiltonian
  do k=1,N_scf_iter
    write(6,*)k
    call conjugate_gradient
    call orthogonalization
    call calculate_density(0.5d0,0.5d0)
    call total_energy(k)
  enddo  
  write(6,*) 'DONE WITH CONJ GRAD GROUND STATE'
  ! Now Perform imaginary time energy minimization
  call imaginary_time
  write(6,*) 'DONE WITH IMAG TIME GROUND STATE'


  ! Now Perform TDDFT
  call time
          
  deallocate(V_POT,density,density_old)
  deallocate(V_exchange,H_Phi,rho_bg)
  deallocate(phi,L_phi,VH,V_ext,wf)
  deallocate(sp_energy,Psi,lattice,lattice_inv,grid_point)

END 

