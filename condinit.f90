!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
use amr_parameters
use hydro_parameters
implicit none
integer ::nn                            ! Number of cells
real(dp)::dx                            ! Cell size
real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
!================================================================
! This routine generates initial conditions for RAMSES.
! Positions are in user units:
! x(i,1:3) are in [0,boxlen]**ndim.
! U is the conservative variable vector. Conventions are here:
! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
! Q is the primitive variable vector. Conventions are here:
! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
! If nvar >= ndim+3, remaining variables are treated as passive
! scalars in the hydro solver.
! U(:,:) and Q(:,:) are in user units.
!================================================================
integer  :: id, iu, iv, iw, ip, iz, ic
real(dp),dimension(1:nvector,1:nvar),save::q          ! Primitive variables
real(dp) :: D_1,P_1,U_1,V_1,W_1
integer  :: i, j, k, ivar, nmode



id = 1
iu = 2
iv = 3
iw = 4
ip = 5



!Initialize primitive variables to condtions defined in namelist
D_1    = d_region(1)                
P_1    = p_region(1)              
U_1    = 0    
V_1    = 0                          
W_1    = 0                            
do i=1,nn 
q(i,id) = D_1
q(i,iu) = U_1
q(i,iv) = V_1
q(i,iw) = W_1
q(i,ip) = P_1
enddo


! Convert primitive to conservative variables
! density -> density
u(1:nn,1)=q(1:nn,1)
! velocity -> momentum
u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5d0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5d0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5d0*q(1:nn,1)*q(1:nn,4)**2
#endif
! thermal pressure -> total fluid energy
u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)






end subroutine condinit
