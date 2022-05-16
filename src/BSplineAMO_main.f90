!> Computes the eigenvalues of a number of states of the quantum
!! harmonic oscillator, and compares them with the analytically
!! known result, \$ E_n = n + 1/2 \$ (when the parameters are chosen
!! such that hbar and the angular frequency are one.
program main 

  use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, OUTPUT_UNIT
  use ModuleBSpline

  implicit none

  !.. Bspline parameters and variables
  integer        , parameter   :: BS_NNODS = 501
  integer        , parameter   :: BS_ORDER =  12
  real(kind(1d0)), parameter   :: BS_GRMIN = -20.d0
  real(kind(1d0)), parameter   :: BS_GRMAX =  20.d0
  real(kind(1d0)), parameter   :: BS_INTER =  BS_GRMAX - BS_GRMIN
  real(kind(1d0))              :: BS_GRID(BS_NNODS)
  type(ClassBSpline)           :: BSpline

  !.. Hamiltonian, overlap, and spectral-decomposition arrays
  real(kind(1d0)), allocatable :: Hmat(:,:), Smat(:,:)
  integer                      :: iNode, nEn, info

  !.. Initializes the BSpline set
  do iNode=1,BS_NNODS
     BS_GRID(iNode) = BS_GRMIN + &
          BS_INTER * dble(iNode-1) / dble(BS_NNODS-1)
  enddo
  call BSpline%Init( &
       BS_NNODS       , &
       BS_ORDER       , &
       BS_GRID        , &
       info           )
  if(info/=0) error stop

  !.. Due to the continuity requirement on the wavefunction,
  !   the first and last BSplines, which do not vanish at the
  !   boundary of the interval, cannot be used
  nEn = BSpline%GetNBsplines() - 2 

  call FillMatrices()

  call GeneralizedEigenproblem()
  
contains

  subroutine FillMatrices()
    integer         :: iBs1, iBs2
    real(kind(1d0)) :: Overlap, KineticEnergy, PotentialEnergy
    real(kind(1d0)) :: parvec(1)
    procedure(D2DFUN), pointer :: fPtrUni, fPtrPow
    fPtrUni => Unity
    fPtrPow => Power
    allocate(Smat(BS_ORDER,nEn))
    Smat=0.d0
    allocate(Hmat,source=Smat)
    parvec(1)=2.d0
    do iBs2 = 2, nEn + 1
       do iBs1 = max(2, iBs2-BS_ORDER+1), iBs2
          Overlap         = BSpline%Integral(fPtrUni,iBs1,iBs2)
          KineticEnergy   = BSpline%Integral(fPtrUni,iBs1,iBs2,1,1) / 2.d0
          !PotentialEnergy = Bspline%Integral(fPtrPow,iBs1,iBs2,parvec=parvec) / 2.d0
          Smat(iBs1+BS_ORDER-iBs2,iBs2-1) = Overlap      
          Hmat(iBs1+BS_ORDER-iBs2,iBs2-1) = KineticEnergy !+ PotentialEnergy
       enddo
    enddo
  end subroutine FillMatrices


  subroutine GeneralizedEigenproblem()
    real(kind(1d0)), allocatable :: Eval(:), Evec(:,:), Work(:)
    real(kind(1d0)), parameter   :: ERROR_THRESHOLD = 1.d-10 
    real(kind(1d0)), parameter   :: PI = 4.d0*atan(1.d0) 
    real(kind(1d0))              :: EigenvalueError, referenceE
    integer                      :: iEn, info
    allocate(work(3*nEn),Eval(nEn),Evec(1,1))
    work=0.d0
    Eval=0.d0
    Evec=0.d0
    call DSBGV( 'N', 'U', nEn, BS_ORDER-1, BS_ORDER-1, &
         Hmat, BS_ORDER, Smat, BS_ORDER, Eval, Evec, 1, work, info )
    if( info /= 0 )then
       write(ERROR_UNIT,"(a,i0)") "DSBGV info : ", info 
       error stop
    endif
    deallocate(work)
    do iEn = 1, nEn
       !EigenvalueError = Eval(iEn) - (iEn-0.5d0)
       !if( EigenvalueError > ERROR_THRESHOLD )exit
       !write(OUTPUT_UNIT,"(i4,x,e24.16)") iEn, EigenvalueError
       referenceE = 0.5 * ( PI / BS_INTER * iEn )**2
       write(OUTPUT_UNIT,"(i4,*(x,e24.16))") iEn, Eval(iEn), referenceE 
    enddo
    write(OUTPUT_UNIT,"(a,i0)") "Number of Accurate Eigenvalues : ",iEn-1
  end subroutine GeneralizedEigenproblem

  
  !.. Lambdas would be so useful ...
  Pure real(kind(1d0)) function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y=1.d0
  end function Unity

  Pure real(kind(1d0)) function Power(x,parvec) result(res)
    real(kind(1d0)), intent(in) :: x
    real(kind(1d0)), optional, intent(in) :: parvec(*)
    res = x**parvec(1)
  end function Power
  
end program main

