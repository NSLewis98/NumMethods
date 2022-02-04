module ModuleRK
    implicit none
    private
    !TODO: GO through all added routines to ensure they work as intended
    !      Start building the solver. Also may need to modify Ks since it might be an array of d-dim vecs
    abstract interface
        real(kind(1d0)) function func(x, y)
            real(kind(1d0)), optional, intent(in) :: x
            real(kind(1d0)), optional, intent(in) :: y(*)
        end function func
    end interface 
    public func

    real(kind(1d0)), parameter :: MIN_STEP_SIZE = 1d-8 ! Is arbitary

    type, public :: ClassRK
        private
        real(kind(1d0)), dimension(:)  , allocatable :: weights, nodes, Ks
        real(kind(1d0)), dimension(:,:), allocatable :: RkMatrix 
        real(kind(1d0))                              :: step_size
        integer                                      :: order
        
        contains
        !.. All of these procedures have been created. They need to be tested
        procedure, public :: Init        => ClassRK_Init
        procedure, public :: GetKs       => ClassRK_GetKs
        procedure, public :: GetPars     => ClassRK_GetPars
        procedure, public :: GetNodes    => ClassRK_GetNodes
        procedure, public :: GetOrder    => ClassRK_GetOrder
        procedure, public :: GetWeights  => ClassRK_GetWeights
        procedure, public :: GetRkMatrix => ClassRK_GetRkMatrix
        procedure, public :: GetStepSize => ClassRK_GetStepSize
        
        procedure, public :: SetNodes    => ClassRK_SetNodes    
        procedure, public :: SetWeights  => ClassRK_SetWeights
        procedure, public :: SetRkMatrix => ClassRK_SetRkMatrix
        procedure, public :: SetStepSize => ClassRK_SetStepSize
        end type ClassRk
    contains
    ! Works
    subroutine ClassRK_Init(self, order, weights, nodes, RkMatrix, step_size, verbose)
        class(ClassRK) ,            intent(inout) :: self
        real(kind(1d0)),            intent(in)    :: weights(:), nodes(:), RkMatrix(:,:)
        real(kind(1d0)),            intent(in)    :: step_size
        integer        ,            intent(in)    :: order
        logical        , optional,  intent(in)    :: verbose

        if( (size(weights).ne.order) .or.(Checkweights(weights,invLogic=.TRUE.))) then
             if(present(verbose)) then
                if(verbose) write(*,*) "ERROR: weightsnot correct size or does not meet cond"
             endif
            return
        elseif( ( (size(RkMatrix, 1) + size(RkMatrix, 2)) / 2) .ne.order ) then
            if(present(verbose)) then
                if(verbose) write(*,*) "ERROR: RkMatrix not correct size"
             endif
            return
        elseif( (size(nodes).ne.order).or.(Checknodes(nodes, RkMatrix, invLogic=.TRUE.)  ) ) then
            if(present(verbose)) then
                if(verbose) write(*,*) "ERROR: nodes not correct size or does not meet cond"
             endif
            return
        else
            if(present(verbose)) then
                if(verbose) write(*,*) "Runge-Kutta Object Succesfully Initialized..."
             endif
            
            if(allocated(self%nodes))    deallocate(nodes)
            if(allocated(self%weights))  deallocate(weights)
            if(allocated(self%RkMatrix)) deallocate(RkMatrix)

            allocate(self%nodes, source = nodes)
            allocate(self%weights, source = weights)
            allocate(self%RkMatrix, source = RkMatrix)

            self%step_size = step_size
            self%order     = order
        endif
    end subroutine ClassRK_Init
!--------------------------------------------   START OF SUBROUTINES THAT NEED CHECKING    -------------------------------------------------------------   
    ! Check if Works
    subroutine ClassRK_GetPars(self, order, weights, nodes, RkMatrix, step_size)
        class(ClassRK) ,            intent(in)     :: self
        real(kind(1d0)),            intent(out)    :: weights(:), nodes(:), RkMatrix(:,:)
        real(kind(1d0)),            intent(out)    :: step_size
        integer        ,            intent(out)    :: order
        logical        , optional,  intent(out)    :: verbose

        if(allocated(self%RkMatrix)) RkMatrix = self%RkMatrix
        if(allocated(self%weights))  weights  = self%weights
        if(allocated(self%nodes))    nodes    = self%nodes

        step_size = self%step_size
        order     = self%order
        
    end subroutine ClassRK_GetPars

    !.. Add logic when self%weights is not allocated
    subroutine ClassRK_GetNodes(self, nodes)
        class(ClassRK) ,              intent(in)    :: self
        integer        , allocatable, intent(out)   :: nodes(:)
        if(allocated(self%nodes)) allocate(nodes, source=self%nodes)
    end subroutine ClassRK_GetNodes

    subroutine ClassRK_GetOrder(self, order)
        class(ClassRK) , intent(in)     :: self
        integer        , intent(out)    :: order
        order = self%order
    end subroutine ClassRK_GetOrder

    !.. Add logic when self%weights is not allocated
    subroutine ClassRK_GetWeights(self, weights)
        class(ClassRK) ,              intent(in)     :: self
        integer        , allocatable, intent(out)    :: weights(:)
        if(allocated(self%weights)) allocate(weights, source=self%weights)
    end subroutine ClassRK_GetWeights
    
    subroutine ClassRK_GetRkMatrix(self, RkMatrix)
        class(ClassRK) ,              intent(in)  :: self
        real(kind(1d0)), allocatable, intent(out) :: RkMatrix(:,:)

        if( allocated(self%RkMatrix)) then 
            allocate(RkMatrix, source = self%RkMatrix)
        endif
    end subroutine ClassRK_GetRkMatrix

    subroutine ClassRK_GetStepSize(self, step_size)
        class(ClassRK) , intent(in)  :: self
        real(kind(1d0)), intent(out) :: step_size
        step_size = self%step_size
    end subroutine ClassRK_GetStepSize

    subroutine ClassRK_GetKs(self, Ks)
        class(ClassRK) ,              intent(in)  :: self
        real(kind(1d0)), allocatable, intent(out) :: ks(:,:)

        if( allocated(self%ks)) then 
            allocate(ks, source = self%ks)
        endif
    end subroutine ClassRK_GetKs

    subroutine ClassRK_SetNodes(self, nodes, verbose)
        class(ClassRK) ,            intent(inout) :: self
        real(kind(1d0)),            intent(in)    :: nodes(:)
        logical, optional,          intent(in)    :: verbose

        if(size(nodes).eq.self%order ) then
            if(allocated(self%nodes)) deallocate(self%nodes)
            allocate(self%nodes, source=nodes)
            return
        else 
            if(present(verbose)) then
                if(verbose) write(*,*) "Error SetNodes: size(nodes).ne. current order"
            endif
        endif
    end subroutine ClassRK_SetNodes

    subroutine ClassRK_SetWeights(self, weights)
        class(ClassRK) ,            intent(inout) :: self
        real(kind(1d0)),            intent(in)    :: weights(:)
        logical, optional,          intent(in)    :: verbose

        if(size(weights).eq.self%order ) then
            if(allocated(self%weights)) deallocate(self%weights)
            allocate(self%weights, source=weights)
            return
        else 
            if(present(verbose)) then
                if(verbose) write(*,*) "Error SetWeigts: size(weights).ne. current order"
            endif
        endif
    end subroutine ClassRK_SetWeights

    subroutine ClassRK_SetRkMatrix(self, RkMatrix)
        class(ClassRK) ,            intent(inout) :: self
        real(kind(1d0)),            intent(in)    :: RkMatrix(:,:)
        logical, optional,          intent(in)    :: verbose

        if(SUM(SHAPE(RkMatrix))/2.eq.self%order ) then
            if(allocated(self%RkMatrix)) deallocate(self%RkMatrix)
            allocate(self%RkMatrix, source=RkMatrix)
            return
        else 
            if(present(verbose)) then
                if(verbose) write(*,*) "Error SetRkMatrix: SHAPE(RkMatrix).ne.(order, order)"
            endif
        endif
    end subroutine ClassRK_SetRkMatrix

    subroutine ClassRK_SetStepSize(self, step_size, verbose)
        class(ClassRK) ,            intent(inout) :: self
        real(kind(1d0)),            intent(in)    :: step_size
        logical, optional,          intent(in)    :: verbose

        if(step_size.lt.MIN_STEP_SIZE)) then
            if(present(verbose)) then
                if(verbose) write(*,*) "ERROR SetStepSize: step_size.lt.MIN_THRESHOLD=1d-8"
            endif
            return
        else
            self%step_size = step_size
        endif
    end subroutine ClassRK_SetStepSize

    !.. Checks if weights meet the following req:
    !.. \sum_{i} \weights_i = 1
    logical function Checkweights(weights, invLogic) result(res) !Works
        real(kind(1d0)),          intent(in) :: weights(:)
        logical        ,optional, intent(in) :: invLogic
        real(kind(1d0)) :: summation

        summation = SUM(weights)
        write(*,*) "summation", summation
        if( (summation.gt.(1d0 - 1d-5) ).and.( summation.lt.(1d0 + 1d-5) )) then 
            res = .TRUE.
        else
            res = .FALSE.
        endif
        if(present(invLogic)) res = XOR(res, invLogic)
    end function Checkweights
! ----------------------------------------  END ROUTINES THAT NEED CHECKING   -------------------------------------------------------------------------------
    !.. Check if \nodes meets the condition of:
    !... \nodes_s = \sum_{j} \RkMatrix_{sj} 
    logical function Checknodes(nodes, RkMatrix, invLogic) result(res) !Works
        real(kind(1d0)),           intent(in) :: nodes(:), RkMatrix(:,:)
        logical        , optional, intent(in) :: invLogic

        real(kind(1d0)), parameter :: THRESHOLD = 1d-5
        integer                    :: row

        res = .TRUE.
        do row = 1, size(RkMatrix, 1)
            if( ABS(SUM(RkMatrix(row, :)) - nodes(row) ).gt.THRESHOLD ) then
                res = .FALSE.
                exit
            endif
        enddo 
        if(present(invLogic)) res = XOR(res, invLogic)
    end function Checknodes
    
end module ModuleRK
