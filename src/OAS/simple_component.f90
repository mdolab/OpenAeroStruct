subroutine mult(nx, ny, x, y)

  implicit none

  integer, intent(in) :: nx, ny
  real*8, intent(in) :: x(nx)
  real*8, intent(out) :: y(ny)

  integer :: i, j

  y(:) = 0.

  do j=1,ny
    do i=1,nx
      y(j) = y(j) + x(i)**2
    end do
  end do

end subroutine mult
