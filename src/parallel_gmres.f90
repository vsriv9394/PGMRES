module pgmres

    use mpi

    implicit none    

    contains

  
        
        
        
        
        
        
        
  ! FUNCTION NORM =====================================================================

    function norm(v)

        implicit none

        real(kind=8)                           :: norm, norm_local
        real(kind=8), dimension(:), intent(in) :: v
        
        integer                                :: i_err=0, root_process=0

        norm_local = sum(v**2)

        call MPI_Reduce(norm_local,          &
                        norm,                &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        root_process,        &
                        mpi_comm_world,      &
                        i_err)

        norm = sqrt(norm)

    end function norm

  
    
    
    
    
    
    
    
  ! FUNCTION DOT_PROD =================================================================

    function dot_prod(v,w)
        
        implicit none

        real(kind=8)                           :: dot_prod, dot_prod_local
        real(kind=8), dimension(:), intent(in) :: v, w
        
        integer                                :: i_err=0, root_process=0

        dot_prod_local = dot_product(v,w)

        call MPI_Reduce(dot_prod_local,      &
                        dot_prod,            &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        root_process,        &
                        MPI_COMM_WORLD,      &
                        i_err)

    end function dot_prod

  
    
    
    
    
    
    
    
  ! SUBROUTINE GMRES ==================================================================

    subroutine gmres(n, LinearOperator, b, x, m, tol, verbose)

        implicit none

        interface
            subroutine LinearOperator(n,v_in,v_out)
                integer,                       intent(in)  :: n
                real(kind=8),    dimension(n), intent(in)  :: v_in
                real(kind=8),    dimension(n), intent(out) :: v_out
            end subroutine LinearOperator
        end interface

      ! Subroutine arguments -------------------------------------------------------

        real(kind=8),    dimension(n),   intent(in)    :: b
        real(kind=8),    dimension(n),   intent(inout) :: x
        real(kind=8),                    intent(in)    :: tol
        integer,                         intent(in)    :: n, m, verbose

      ! Local arguments ------------------------------------------------------------

        integer                        :: k
        integer                        :: this_proc, ierr
        real(kind=8)                   :: b_norm, res_norm, error
        real(kind=8), dimension(n)     :: res, x_temp
        real(kind=8), dimension(m+1)   :: sn, cs
        real(kind=8), dimension(n,m+1) :: Q
        real(kind=8), dimension(m+1,m) :: H
        real(kind=8), dimension(m+1)   :: beta

      ! Subroutine content ---------------------------------------------------------
        
        call mpi_comm_rank(mpi_comm_world, this_proc, ierr)

        b_norm = norm(b)

        res = 0
        x_temp = 0
        sn = 0
        cs = 0
        Q = 0
        H = 0
        beta = 0

        call LinearOperator(n,x,x_temp)
        res = b - x_temp
        res_norm = norm(res)

        error   = res_norm / b_norm
        beta(1) = res_norm

        Q(:,1)  = res / res_norm

        do k=1,m

            call arnoldi(LinearOperator, Q, H, k, n, m)
            call apply_givens_rotation(H, cs, sn, k)
            
            beta(k+1) = -sn(k)*beta(k)
            beta(k)   =  cs(k)*beta(k)

            error = abs(beta(k+1))/b_norm
            
            if (this_proc==0 .and. verbose==1) then
            
                write(*,"(A,I9,A,E14.6)") "                           Iteration    ", k, "    Error    ", error

            end if

            if (error<tol) then
                
                exit

            end if

        end do

        call back_substitute(H(1:k,1:k), beta(1:k))
        x = x + matmul(Q(:,1:k), beta(1:k))

    end subroutine gmres

  
    
    
    
    
    
    
    
  ! SUBROUTINE ARNOLDI ================================================================

    subroutine arnoldi(LinearOperator, Q, H, k, n, m)

        implicit none

        interface
            subroutine LinearOperator(n,v_in,v_out)
                integer,                       intent(in)  :: n
                real(kind=8),    dimension(n), intent(in)  :: v_in
                real(kind=8),    dimension(n), intent(out) :: v_out
            end subroutine LinearOperator
        end interface

      ! Subroutine arguments -------------------------------------------------------

        integer,                           intent(in)    :: k, n, m
        real(kind=8),    dimension(n,m),   intent(inout) :: Q
        real(kind=8),    dimension(m+1,m), intent(inout) :: H

      ! Local arguments ------------------------------------------------------------

        real(kind=8), dimension(n) :: x_temp
        integer                    :: i

      ! Subroutine content ---------------------------------------------------------
        
        call LinearOperator(n,Q(:,k),Q(:,k+1))

        do i=1,k

            H(i, k)  = dot_prod(Q(:,k+1), Q(:,i))
            Q(:,k+1) = Q(:,k+1) - H(i,k) * Q(:,i)

        end do

        H(k+1, k) = norm(Q(:,k+1))
        Q(:,k+1)  = Q(:,k+1) / H(k+1, k)

    end subroutine arnoldi

  
    
    
  
    
    
    
    
  ! SUBROUTINE APPLY_GIVENS_ROTATION ====================================================

    subroutine apply_givens_rotation(H, cs, sn, k)

        implicit none

      ! Subroutine arguments -------------------------------------------------------

        real(kind=8),    dimension(:,:), intent(inout)   :: H
        integer,                         intent(in)      :: k
        real(kind=8),    dimension(:),   intent(inout)   :: cs, sn

      ! Local arguments ------------------------------------------------------------

        real(kind=8)    :: temp
        integer         :: i

      ! Subroutine content ---------------------------------------------------------

        do i=1,k-1

            temp     =  cs(i)*H(i,k) + sn(i)*H(i+1,k)
            H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
            H(i,k)   =  temp

        end do

        if (H(k,k)==0) then

            cs(k) = 0.0D0
            sn(k) = 1.0D0

        else

            temp  = sqrt(H(k,k)**2 + H(k+1,k)**2)
            cs(k) = abs(H(k,k)) / temp
            sn(k) = cs(k) * H(k+1,k) / H(k,k)

        end if

        H(k,k)   = cs(k)*H(k,k) + sn(k)*H(k+1,k)
        H(k+1,k) = 0.0D0

    end subroutine apply_givens_rotation

  
    
    
    
    
    
    
    
  ! SUBROUTINE BACK_SUBSTITUTION =========================================================

    subroutine back_substitute(H, beta)

        implicit none
        
        real(kind=8),    dimension(:,:), intent(in)    :: H
        real(kind=8),    dimension(:),   intent(inout) :: beta

        integer :: i, k

        k = size(beta)

        beta(k) = beta(k)/H(k,k)

        do i=k-1,1,-1

            beta(i) = (beta(i) - dot_product(H(i,i+1:k),beta(i+1:k)))/H(i,i)

        end do

    end subroutine back_substitute

end module pgmres
