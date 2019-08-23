program main

      call parallel_solve()

end program main

subroutine parallel_solve()

      use pgmres
      implicit none

      real(kind=8), dimension(5)  :: b
      real(kind=8), dimension(5)  :: x
      integer                     :: m=5, n=5
      integer                     :: i_err

      interface
          subroutine Laplacian(n,v_in,v_out)
              integer,                       intent(in)  :: n
              real(kind=8),    dimension(n), intent(in)  :: v_in
              real(kind=8),    dimension(n), intent(out) :: v_out
          end subroutine Laplacian
      end interface

      call mpi_init(i_err)

      b(1) = 1.0D0
      b(2) = 2.0D0
      b(3) = 4.0D0
      b(4) = 5.0D0
      b(5) = 3.0D0

      x = 1.0D0

      write(*,*)
      call gmres(n, Laplacian, b, x, m, 1.0D-7, 1)
      write(*,*)
      write(*,*) "Solution:"
      write(*,*) x

      call mpi_finalize(i_err)

end subroutine parallel_solve

subroutine Laplacian(n, v_in, v_out)
      
      implicit none
      
      integer,                       intent(in)    :: n
      real(kind=8),    dimension(n), intent(in)    :: v_in
      real(kind=8),    dimension(n), intent(out)   :: v_out

      v_out(1) = 4*v_in(1) + 2*v_in(2) + 7*v_in(3) + 1*v_in(4) + 0*v_in(5)
      v_out(2) = 1*v_in(1) + 3*v_in(2) + 0*v_in(3) + 5*v_in(4) + 7*v_in(5)
      v_out(3) = 9*v_in(1) + 5*v_in(2) + 8*v_in(3) + 2*v_in(4) + 5*v_in(5)
      v_out(4) = 0*v_in(1) + 7*v_in(2) + 7*v_in(3) + 3*v_in(4) + 7*v_in(5)
      v_out(5) = 2*v_in(1) + 1*v_in(2) + 5*v_in(3) + 7*v_in(4) + 8*v_in(5)

end subroutine Laplacian
