! ****************************************************
!
!          Regularized Choleski Algorithm
!
!  This program calculates the regularized solutions
!  for a given linear -  possibly singular - system Ax=b
!  where the matrix A is positive semidefinite.
! 
!  Developed by: Sopasakis Pantelis & Patrinos Panos
!  in Fortran 90
!
!  National Technical University of Athens
!  School of Chemical Engineering
!  Automatic Control Laboratory
!  version 1.1.2 (pre-alpha)
!
! ****************************************************
program chol

! Definitions:
! n: the dimension of the problem. dimA=dimY=n
implicit none
integer:: n, flag=0, p
parameter (n=1000)
real, allocatable, dimension(:,:):: A,Y
integer, allocatable, dimension(:):: JJ,b,BB,z,X
real:: sum1=0, sum2=0, sum3=0, sum4=0, sum5=0,sum6=0, W=0, G=0, start, finish
integer:: i,j,k,l

! Allocation of A, Y and JJ
allocate(A(n,n), Y(n,n))
allocate(JJ(n), b(n))
allocate(BB(n))

! *************************************** 
!  Define the matrix A and the vector b:
! ***************************************
!$OMP PARALLEL
!$OMP DO
diagonal: do i=1,n
	   A(i,i)=100+i
	   b(i)=2+i
	  enddo diagonal

!$OMP END DO


!$OMP DO
subHyperDiagonal: do i=1,n-1
		    A(i,i+1)=10+i
		    A(i+1,i)=11+i
		  enddo subHyperDiagonal
!$OMP END DO

! **************************************
! Step 1:
! Calculate Y(i,1)
! **************************************
CALL cpu_time(start)
!Step-1:
flag=1
step_1_If: if (A(1,1).le.0) then
  	    JJ(flag)=1
  	    flag=flag+1 
	    Y(1,1)=1
             ! ----> The following loop (inner_1_1) should be deleted...
	     !$OMP DO
	     inner_1_1: do i=2,n
  		         Y(i,1)=0
 	                enddo inner_1_1
	     !$OMP END DO
	     else
  		Y(1,1)=sqrt(A(1,1))
             !$OMP DO
	        inner_1_2: do i=2,n
   		  	     Y(i,1)=A(i,1)/Y(1,1)
  			   enddo inner_1_2
	     !$OMP END DO
	  endif step_1_If


! **************************************
! Step 2:
! Calculate Y(i,j) for i=1..n, j=1..n-1
! **************************************

!$OMP DO
step2loop: do k=2,n-1
  sum1=0
  do j=1,k-1
   sum1=sum1+Y(j,j)**2
  enddo
  W=A(k,k)-sum1
  
  if (W.le.0) then
    JJ(flag)=k
    flag=flag+1
    Y(k,k)=1
    do i=1,k-1
     Y(k,i)=0
    enddo
    do i=k+1,n
     Y(k,i)=0
    enddo
  else
    Y(k,k)=sqrt(W)
    do i=k+1,n
	sum2=0;
        do l=1,k-1
	  sum2=sum2+Y(k,l)*Y(i,l) 
	enddo
     Y(i,k)=A(i,k)-(  sum2/Y(k,k)  )
    enddo
  endif
enddo step2loop
!$OMP END DO

! **************************************
! Step 3:
! Calculate Y(n,n)
! **************************************
sum3=0
!$OMP DO
step3loop1: 	do j=1,n-1
 		 sum3=sum3+Y(n,j)**2
		enddo step3loop1
!$OMP END DO
G=A(n,n)-sum3;
if (G.le.0) then
  JJ(flag)=k;
  flag = flag+1
  Y(n,n)=1
  !$OMP DO
  do i=1,n-1
   Y(n,i)=0
  enddo
  !$OMP END DO
else
  Y(n,n)=sqrt(G);
endif
!$OMP END PARALLEL


! *****************************************************
! Step 4:
! Calculate BB(k) = b(k) - SUM_[i in JJ] {A(k,i)*b(i)}
! *****************************************************
!$OMP PARALLEL
flag = flag - 1
if (flag.gt.0) then
	!$OMP DO
	do k=1,n
	   if (isIn(JJ,k)) then
		BB(k)=b(k)
	   else
		sum4=0
		do p=1,flag
 		  sum4=sum4+A(k,JJ(p))*b(JJ(p))
	        enddo

	      	BB(k)=b(k)-sum4
	   endif	
	enddo
	!$OMP END DO		

endif
deallocate(A)

! **************************************
! Step 5:
! Calculate z(k), k=1..n
! **************************************
allocate(z(n))
z(1)=BB(1)/Y(1,1)
!$OMP DO

do k=2,n
	sum5=0
	do j=1,k-1
		sum5=sum5+Y(k,j)*z(j)
	enddo
	z(k)=(BB(k)-sum5)/Y(k,k)
enddo
!$OMP END DO


! **************************************
! Step 6:
! Calculate x(k), k=1..n
! **************************************
allocate(X(n))
X(n)=z(n)/Y(n,n)
!$OMP DO
do k=n-1,-1,1
	sum6=0
	do j=k+1,n
		sum6=sum6+Y(j,k)*z(j)
	enddo
	X(k)=(z(k)-sum6)/Y(k,k)
enddo
!$OMP END DO

CALL cpu_time(finish)
!$OMP END PARALLEL
print *,"TIME=",finish-start


!***************************************
! Finalization:
! Deallocate Matrices
!***************************************
deallocate (z,Y,b,BB,JJ,X)


! **************************************
!
!	Subroutines
!
! **************************************
CONTAINS


		! Checks if a given integer belongs to an array
		logical function isIn(J,n) result (belongs)
		 integer::n
		 integer, allocatable, dimension(:):: J
			
			do i=1,size(J)
			  if (n.eq.J(i)) then
				belongs=1
				return
			  endif
			enddo 
			
		end function 


end program chol
