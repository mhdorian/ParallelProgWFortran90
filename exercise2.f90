!mdorian Übung 6 exercise 2
program Monte

  !$use omp lib
  implicit none

  integer :: dimensions, Radius !num of dimensions, radius of circle
  integer :: Nin, N, i,j, k !number of points "darts" within target volume, 
                         !Total darts to be thrown, iteration counters 
  real(4) :: Avol !Acirc !approximate volume of multi-dimensional sphere
                      !approximate area of 2d circle used in part A
  real(4) :: denom, pi, V !exact values of pi and sphere 
  real, dimension(:), allocatable :: dart !array holding coordinate points of dart
  
  pi = 3.14159265
  
  Radius = 1 !unit circle
  dimensions = 2 !set dimentions (4d for now)
  N = 10000 !set number of darts
  !Acirc = 0
  Avol = 0
  Nin = 0 
  
  
  !seed random generator
  call random_seed()
  !initialize dart
  allocate(dart(1:dimensions))
  
  !open write file
  !open(unit= 10, file= "montecarlo4d.dat")
  
  !---begin paralell do--------------------------
  !$omp parallel do private(dart, i) shared(Radius) &
  !$omp reduction(+: Nin)
  do i = 1, N
    
    !generate dart
    do j = 1, dimensions
      call random_number(dart(j))
    enddo
    
    !check if within target
    if (Radius >= sum(dart**2)) then 
      Nin = Nin + 1 !if in volume, increment Nin
    end if
    
    !calulation here when collecting individual plots.
    !moved to outside of loop for finial approximation.
    !Avol = (2**dimensions) * (real(Nin)/i)
    !Acirc = 4 * (real(Nin)/i)
    
    !plot(i, Acirc)
    !write(10, *) i,Avol 
    
  end do 
  !$omp end parallel do
  
  Avol = (2**dimensions) * (real(Nin)/N)
  
  !check against exact volume Note: does not check for bad dimension input
  if (mod(dimensions, 2) == 0) then !even case
    k = dimensions / 2
    denom = 1
    do i= 1, k
      denom = denom * i
    enddo
    V = (pi**k)/denom
  else !odd case
    k = (dimensions - 1) / 2
    denom = 1
    do i = 1,k
      denom = denom * ( ((2*k)-1) * ((2*k)+1) )
    enddo
    V = ((2**(k+1)) * (pi**k)) / denom 
  end if
  
  
  
  print *, "Approximate volume: ", Avol
  print *, "Exact volume: ", V
  !close write file
  !close(10)
  

end program Monte