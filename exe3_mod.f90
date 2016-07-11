module sparseMatrix

  use omp_lib
  implicit none 

  type sparseMatrix_t
    !integer :: n !size/dimention of matrix
    integer, dimension(:), allocatable :: val, col, row 
  end type sparseMatrix_t
  
 contains
 
 !---convert to sparse matrix--------------------------------
  function convertMatrix(M, msize, numVals) result (spMatrix)
    integer, intent(in) :: msize, numVals
    integer, dimension(:,:), intent(in) :: M(msize, msize)
    
    type(sparseMatrix_t) :: spMatrix
    integer :: i,j, valIndex
  
    !spMatrix%n = msize
    allocate(spMatrix%val(1: numVals))
    allocate(spMatrix%col(1: numVals))
  
    allocate(spMatrix%row(1:msize + 1))
    spMatrix%row(msize+1) = numVals + 1
  
    valIndex = 1
    do i = 1, msize
      !store index where row changes
      spMatrix%row(i) = valIndex
      do j = 1, msize
	if (M(i,j) /= 0) then
	  !store value and colum number 
	  spMatrix%val(valIndex) = M(i,j)
	  spMatrix%col(valIndex) = j
	
	  !increment valIndex
	  valIndex = valIndex + 1
	end if 
      enddo 
    
    end do
  
  !returns spMatrix
  end function convertMatrix
  !---end convert to sparse matrix--------------------------------
  
  !---vector-sparse matrix mult--------------------------------
  function spMatrixVectorMult(spM, v, sizeN) result (w)
    integer, intent(in) :: sizeN
    type(sparseMatrix_t), intent(in) :: spM
    integer, dimension(:), intent(in) :: v(sizeN)
    
    integer, dimension(:) :: w(sizeN)
    integer :: i,j
    
    w = (/(0, i=1,sizeN)/)
    
    !omp parallel do
    do i = 1,sizeN
      do j = spM%row(i), (spM%row(i+1) - 1)
	w(i) = w(i) + spM%val(j) * v(spM%col(j))
      enddo
    end do
    !omp end parallel do
  !returns w
  end function spMatrixVectorMult
  !---end vector-sparse matrix mult--------------------------------
  
  
  !---print sparseMatrix_t----------------------------------------
  subroutine printSparseMatrix(M_t)
    type(sparseMatrix_t) :: M_t
    integer :: i
    
    print *, "Val: ", (M_t%val(i), i=1, size(M_t%val))
    print *, "Col: ", (M_t%col(i), i=1, size(M_t%col))
    print *, "Row: ", (M_t%row(i), i=1, size(M_t%row))
  end subroutine printSparseMatrix
  !---end print sparseMatrix_t----------------------------------------  

end module sparseMatrix

