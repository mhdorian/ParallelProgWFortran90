program sparse

  use sparseMatrix
  !$ use omp_lib
  implicit none
  
  integer, parameter :: n = 10
  integer :: A(n,n)
  integer :: vector(n), r(n), s(n)
  integer :: i,j
  real(8) :: t1, t2, tTotal
  type(sparseMatrix_t) :: Asparse
  
  A = reshape( (/ ((0, i=1,n), j=1,n) /), shape = (/n,n/) )
  A(2,3) = 5
  vector = (/(1, i=1,n)/)
  !--------------------------------------------------
  
  print *, "Here is your matrix:"
  call printMatrix(A,n)
  print *, "Here is your vector:"
  print *, (vector(i), i=1,size(vector))
  
  print *, "converting to sparse matrix..."
  Asparse = convertMatrix(A, n, 1)
  print *, "Here is your sparse matrix:"
  call printSparseMatrix(Asparse)
  
  !------------------------------------------------
  
  print *, "normal matrix-vector mult..."
  
  call cpu_time(t1)
  !$ t1 = omp_get_wtime()
  r = matrixVectorMult(A, vector, size(vector))
  call cpu_time(t2)
  tTotal = t2-t1
  !$ tTotal = omp_get_wtime() - t1
  
  print *, "Result: ", (r(i), i=1,size(r))
  print *, "Time taken to compute (s): ", tTotal 
  
  !---------------------------------------------------
  
  print *, "sparse matrix-vector mult..."
  
  call cpu_time(t1)
  !$ t1 = omp_get_wtime()
  s = spMatrixVectorMult(Asparse, vector, size(vector))
  call cpu_time(t2)
  tTotal = t2-t1
  !$ tTotal = omp_get_wtime() - t1
  
  print *, "Result: ", (s(i), i=1,size(s))
  print *, "Time taken to compute (s): ", tTotal 
  
  print *, "Done."
  
  
!---------------------------------------------------------------------------
 contains
 !---print normal matrix(2D array)----------------------------------------
  subroutine printMatrix(M, sizeM)
    integer, intent(in) :: sizeM
    integer, intent(in) ::  M(sizeM, sizeM)
    integer :: i,j 
    
    do i = 1,sizeM
      print *, (M(i,j), j=1, sizeM)
    end do
  
  end subroutine printMatrix
 !-------------------------------------------
 
 !---brute force Matrix-Vector mult----------------------------------------
  function matrixVectorMult(M, v, sizeN) result (w)
    integer, intent(in) :: sizeN
    integer, dimension(:), intent(in) :: v(sizeN)
    integer, dimension(:,:), intent(in) :: M(sizeN, sizeN)
    integer, dimension(:) :: w(sizeN) 
    integer :: i,j
    
    w = (/(0, i=1,sizeN)/)
    !$omp parallel do   
    do i=1, sizeN
      do j=1, sizeN
	w(i) = w(i) + ( M(i,j) * v(j) )
      enddo
    end do
   !omp end parallel do 
  
  !returns w
  end function matrixVectorMult
 !-------------------------------------------

end program sparse