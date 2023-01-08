!  benchinverse.f90 
!
!  FUNCTIONS: Another version like David Frank  Dave_Frank@hotmail.com TEST_FPU
!  benchinverse - Benchmarking for solving inverse matrix.
!  Solver 1: BLAS LAPACK -> Need MKL or Openblas or ACML
!  Solver 2: DSS Intel -> Need MKL
!  Solver 3: Pardiso Intel -> Need MKL
!  Solver 4a: SuperLU   -> Need BLAS from MKL or OpenBlas
!  Solver 4b: SuperLUMT -> Need BLAS from MKL or OpenBlas with Openmp
!  Solver 4c: SuperLU per 1 column solve  -> Need BLAS from MKL or OpenBlas with Openmp
!  Solver 5 : UMFPACK still developing
!  still developing

!  Depedency for intel DSS, BLAS, Pardiso 
!  Install Intel MKL in OS, and compile it using /Qmkl in windows or -Qmkl in linux/mac
!  
!  Depedency for SuperLu 
!  

!****************************************************************************
!
!  PROGRAM: benchinverse
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program benchinverse
    
    implicit none
    
    !   variables input directass
    character*12       :: sizeString
    integer (kind = 4) :: stat,n
    character*1    ::    button
    
    ! Variables
    integer (kind = 4) :: bigsize,bigsize2
    real(kind = 8), dimension(:), allocatable :: mat(:,:),mat2(:,:) ! random numM2DtoCSRbers to invert and check 2d
    real(kind = 8) errors
    character*3 :: inp1,inp2
    !var cases
    integer :: mycase
    logical :: sparse,symmetry
    
    
    !var filing time
    integer (kind = 4) t(8)

!   sparse matrix generator
    integer (kind = 4) :: m,o,c,index,nn, nz, fejlm,nz2,nz1
    integer (kind = 4) , ALLOCATABLE :: snr(:), rnr(:)
    REAL( kind = 8), ALLOCATABLE :: val(:)
    REAL( kind = 8) :: alphax,u

    !var DNS2CSR
    integer ( kind = 4 ) :: ncol
    integer ( kind = 4 ) :: ndns
    integer ( kind = 4 ) :: nrow

    real ( kind = 8 ), allocatable :: a(:)
    integer ( kind = 4 ) :: i
    integer ( kind = 4 ), allocatable ::  ia(:) !(nrow+1)
    integer ( kind = 4 ) :: ierr
    integer ( kind = 4 ) :: j
    integer ( kind = 4 ), allocatable :: ja(:)
    integer ( kind = 4 ) ::next
    integer ( kind = 4 ) :: nzmax
    
    !var CSR2CSC
    real ( kind = 8 ), allocatable :: ao(:)
    integer ( kind = 4 ), allocatable ::  iao(:) !(nrow+1)
    integer ( kind = 4 ), allocatable :: jao(:)
    
    !   variables for benchmarking
    real ( kind = 8 ) :: avg_err, dt(20),avg_err2,pass(3,2)
    integer ( kind = 4 ) :: clock1, clock2,clock3, rate
    logical   :: benchmarkmode


    ! Body of DENSE2HB
    bigsize = 22 !default
    
    n = command_argument_count()
101 if (n > 0) then
        
        call get_command_argument(1, sizeString)
        read(sizeString,*) bigsize
        allocate (mat(bigsize,bigsize),stat=stat )
        if (stat/=0) then 
            stop 'Cannot allocate memory for the main, Argument 1 is size matrix ' 
        endif    
        if (bigsize.lt.22)  then 
            bigsize = 22
            WRITE (*,*)  "Min sparse matrix 22x22 due sparsekit matrix generator"
        endif
        if (bigsize.gt.9000000)  then 
            bigsize = 9000000
            WRITE (*,*)  "Min sparse matrix 9000000 due sparsekit matrix generator"
        endif
        
        !default
        symmetry = .false.
        sparse = .false.
        
102     if (n.ge.2) then
            call get_command_argument(2, sizeString)
            read(sizeString,1000) inp1
            if ((inp1.eq.'SYM').or.(inp1.eq.'sym')) then
                symmetry = .true.
                WRITE (*,*)  "Creating symmetric matrix"
            else
                symmetry = .false.
                WRITE (*,*)  "Creating unsymmetric matrix"
            endif
        endif

103     if (n.ge.3) then
            call get_command_argument(3, sizeString)
            read(sizeString,1000) inp2
            if ((inp2.eq.'SPA').or.(inp2.eq.'spa')) then
                sparse = .true.
                WRITE (*,*)  "Creating sparse matrix"
            else
                sparse = .false.
                WRITE (*,*)  "Creating dense matrix"
            endif
        endif    
    else
        allocate (mat(bigsize,bigsize),stat=stat )
        !default
        symmetry = .false.
        sparse = .false.
        if (stat/=0) then 
            write (*,*) 'Cannot allocate memory for the main' 
            write (*,*) 'Argument 1 is size matrix'
            write (*,*) 'Argument 2 is sym/unsym, Argument 3 is sparse/dense'
            stop 'Please Check'  
        endif    
        
    endif
    
1000 format(a3)
        
    ! Setting zeros for main matrix
    mat = 0.d0
    n = bigsize


    if (.not.sparse) mycase = 1
    if (sparse) mycase = 2
    
	if (mycase.eq.1) then
	    CALL DATE_AND_TIME (values = t)
        CALL RANDOM_SEED()              ! set seed to random number based on time
        CALL RANDOM_NUMBER(mat)         ! fill pool with random data ( 0. -> 1. )
        if (symmetry) then 
            call makesym (mat,bigsize) !from lower
        endif
    endif
    
    if (mycase.eq.2) then
        mat = 0.d0
        ! m and 21 <n < 9000001
        !m row n col c sparsity n~n-10 random 10 < C < N-10   random (1 < INDEX < N-C-8) alpha = 1
        !NN INDEX*M+109 < NN < 9000001 must be specified. nz nonzero, avalues, snr=col rnr=row
        !call matrf2 ( m, n, c, index, alpha, nn, nz, a, snr, rnr, fejlm )
        !for c
        
        CALL RANDOM_NUMBER(u) ! m=10 and o=(n-10)
        m = 10
        o = n-10
        c = m + FLOOR((o+1-m)*u)  ! We want to choose one from m-n+1 integers
        c = max(11,c)
        c = min(11,n-10)
        
        !for index
        CALL RANDOM_NUMBER(u)
        m = 1
        o = n-c-8
        index = m + FLOOR((o+1-m)*u)  ! We want to choose one from m-n+1 integers
        index = max(2,index)
        index = min(2,n-c-8)
        alphax = 1
        
        !for nn
        CALL RANDOM_NUMBER(u)
        m = index*n+109
        o = n*n
        nn = m + FLOOR((o+1-m)*u)  ! We want to choose one from m-n+1 integers
        nn = max(INDEX*M+109+1,nn)
        nn = min(nn,n*n)
        
        allocate (val(nn))
        allocate (snr(nn))
        allocate (rnr(nn))
        val = 0.d0
        call matrf2 (n,n,c,index,alphax,nn,nz,val,snr,rnr,fejlm)
        do m=1,nz
            if ((snr(m).gt.n).or.(snr(m).le.0)) stop "ERROR SNR"
            if ((rnr(m).gt.n).or.(snr(m).le.0)) stop "ERROR RNR"
        enddo    
        
        if (fejlm.ne.0) then
            write (*,*) "ERROR in SPARSE MATRIX GEN; ERROR CODE ", fejlm
            !if (fejlm.eq.0) write (*,*) "0, indicates that the call is successful."
            if (fejlm.eq.1) write (*,*) "1, N is out of range."
            if (fejlm.eq.2) write (*,*) "2, M is out of range."
            if (fejlm.eq.3) write (*,*) "3, C is out of range."
            if (fejlm.eq.4) write (*,*) "4, INDEX is out of range."
            if (fejlm.eq.5) write (*,*) "5, NN is out of range."
            if (fejlm.eq.7) write (*,*) "7, ALPHA is out of range."
            stop
        endif
        mat = 0.d0
        do o=1,nz
            mat(rnr(o),snr(o))=val(o)
        enddo    
        deallocate(val,snr,rnr)
        if (symmetry) then 
            call makesym (mat,bigsize) !from lower

        endif
    endif
    
    call matrix_check_sym (mat,bigsize, symmetry, sparse,nz)
    nz1 = nz
    open(10,file='matrix.dat',form='unformatted')
    write(10) bigsize,mat
    close(10)


    !test accuracy if needed from writing matrix to file
    allocate(mat2(bigsize,bigsize))
    errors = 0.d0
    open(10,file='matrix.dat',form='unformatted')
    read(10) bigsize2,mat2
    close(10)
    do i=1,bigsize
        do j=1,bigsize
            if (abs(mat2(i,j)-mat(i,j)) > 0.d0)  errors = abs(mat2(i,j)-mat(i,j)) + errors    
        enddo
    enddo    
    deallocate(mat2)
    print *, ' Error read = ', errors 
    
     !CSR MATRIX
     ! nrow  = bigsize
     ! ncol  = bigsize
     ! nzmax = nz
     ! dns = mat  
     ! ndns = nrow = bigsize ! real DNS(NDNS,NCOL), an NROW by NCOL dense matrix
     ! output a ! values of a(nz)
     ! output ja !ja(nz) 
     ! output ia  ia(NROW+1)
     ! ierror must be 0 or definitely error
     ! call dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
     allocate(a(nz),ja(nz),ia(bigsize+1))
     a = 0.d0
     ja = 0
     ia = 0
     call dnscsr ( bigsize, bigsize, nz, mat, bigsize , a, ja, ia, ierr )   
     deallocate(mat)
     
     if (ierr.ne.0) stop ' dnscsr ERROR'
     open (11,file='matrix_csr.dat',form='unformatted')
     write(11) nz,bigsize,a,ja,ia
     close(11)
     
    !print matrix csr?
    !write (*,'(A)',advance='no') 'Enter a value of ns '
    !read (*,*) ns
    write(*,'(A)',advance='no')  '  Print the CSR Matrix? ' 
    read(*,*)  button
    if ((button.eq.'y').or.(button.eq.'Y')) then 
        open (13,file='matrix_csr.txt',form='formatted')
        call dump ( bigsize, a, ja, ia, 13 )
        close(13)
    endif   
    
    allocate(ao(nz),jao(nz),iao(bigsize+1))        
    call csrcsc ( bigsize, 1, 1, a, ja, ia, ao, jao, iao )
    deallocate(a,ja,ia)
     open (14,file='matrix_csc.dat',form='unformatted')
     write(14) nz,bigsize,ao,jao,iao
     close(14)
    deallocate(ao,jao,iao)
     

    benchmarkmode = .true.
    n = bigsize
    open(10,file='matrix.dat',form='unformatted')
    read(10) bigsize2
    close(10)
 

    if (bigsize2.eq.bigsize) then
        allocate(mat2(bigsize,bigsize))
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize2,mat2
        close(10)
    else
        write (*,*) bigsize2
        stop ' Size matrix is mismatch - BLAS'
    
    endif    
    
    CALL SYSTEM_CLOCK (clock1,rate)  ! get benchmark (n) start time     
    write(*,*) '  STAGE 1 - Inverse'
    call matrix_check_sym (mat2,bigsize, symmetry, sparse,nz)
    pass(1,1) = mat2(1,1)
    pass(1,2) = mat2(n,n)

    CALL BLASLAPACK_inverse(mat2,n) !solution as input and output
    CALL SYSTEM_CLOCK (clock2,rate)  ! get benchmark (n) start time
    pass(2,1) = mat2(1,1)
    pass(2,2) = mat2(n,n)	

    dt(1) = (clock2-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)" Test1-1: DGETRFI- ", 1," (",n,"x",n,") inverts in ", dt(1), ' seconds'


    write(*,*) '  STAGE 2 - Reinverse'
    call matrix_check_sym (mat2,bigsize, symmetry, sparse,nz)

    CALL BLASLAPACK_inverse(mat2,n) !solution as input and output
    CALL SYSTEM_CLOCK (clock3,rate)  ! get benchmark (n) start time
    pass(3,1) = mat2(1,1)
    pass(3,2) = mat2(n,n)	
    dt(1) = (clock3-clock2)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)" Test1-2: DGETRFI- ", 1," (",n,"x",n,") inverts in ", dt(1), ' seconds'
     
    write(*,*) ' PASS0 :', pass(1,1),pass(1,2)
    write(*,*) ' PASS1 :', pass(2,1),pass(2,2)
    write(*,*) ' PASS2 :', pass(3,1),pass(3,2)

    
    if (stat/=0) stop 'Cannot allocate memory for input working matrix'    
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize2
    close(10)
    
    if (bigsize2.eq.bigsize) then
        allocate (mat(n,n),stat=stat  )
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize,mat
        close(10)
    else
        write (*,*) bigsize2
        stop ' Size matrix is mismatch - BLAS 2'
    endif    
  
	avg_err = SUM(ABS(mat-mat2))/(n*n)    ! invert err.
    dt(1) = (clock3-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)"  Test1: DGETRFI- ", 2," (",n,"x",n,") inverts in ", dt(1), ' seconds  Err=', avg_err
42  FORMAT (A,I4,A,I6,A,I6,A,F7.3,A,F14.12)
    deallocate(mat,mat2)
    
    ! Part4 SuperLU - CSC
    ! Set the problem to be solved for DSS
    allocate (mat2(n,n)) !
    
    if (stat/=0) stop 'Cannot allocate memory for input working matrix'    
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize2
    close(10)
    
    if (bigsize2.eq.bigsize) then
        allocate (mat(n,n),stat=stat  )
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize,mat2
        close(10)
    else
        write (*,*) bigsize2
        stop ' Size matrix is mismatch - Superlu 1'
    endif 
    
    write(*,*) '  STAGE 1 - Inverse'    
    CALL matrix_check_sym (mat2,n, symmetry, sparse,nz)
    allocate(ao(nz),jao(nz),iao(bigsize+1))
    
    pass(1,1) = mat2(1,1)
    pass(1,2) = mat2(n,n)
	
!   2D matrix

    CALL SYSTEM_CLOCK (clock1,rate)  ! get benchmark (n) start time    
    CALL Superlucall2d (mat2,n,nz,ao,jao,iao) !solution as input and output

    deallocate(ao,jao,iao)
    pass(2,1) = mat2(1,1)
    pass(2,2) = mat2(n,n)	
    
    CALL SYSTEM_CLOCK (clock2,rate) 
    dt(4) = (clock2-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec
    WRITE (*,42)" Test4-1: SUPERLU2- ", 1," (",n,"x",n,") inverts in ", dt(4), ' seconds' 


    write(*,*) '  STAGE 2 - Reinverse'
    CALL matrix_check_sym (mat2,n, symmetry, sparse,nz) !get nz
    allocate(ao(nz),jao(nz),iao(bigsize+1))
    
    CALL Superlucall2d (mat2,n,nz,ao,jao,iao) !solution as input and output
	deallocate(ao,jao,iao)
    
    pass(3,1) = mat2(1,1)
    pass(3,2) = mat2(n,n)	
    
    CALL SYSTEM_CLOCK (clock3,rate) 
    dt(4) = (clock3-clock2)/DBLE(rate)  ! get benchmark (n) elapsed sec
    WRITE (*,42)" Test4-2: SUPERLU2- ", 1," (",n,"x",n,") inverts in ", dt(4), ' seconds' 


    write(*,*) ' PASS0 :', pass(1,1),pass(1,2)
    write(*,*) ' PASS1 :', pass(2,1),pass(2,2)
    write(*,*) ' PASS2 :', pass(3,1),pass(3,2)

    if (stat/=0) stop 'Cannot allocate memory for input working matrix'    
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize2
    close(10)
    
    if (bigsize2.eq.bigsize) then
        allocate (mat(n,n),stat=stat  )
        open(10,file='matrix.dat',form='unformatted')
        read(10) bigsize,mat
        close(10)
    else
	   write (*,*) bigsize2
        stop ' Size matrix is mismatch - Superlu 2'
    endif    
        
    avg_err = SUM(ABS(mat2-mat))/(n*n)    ! invert err.
    
    CALL SYSTEM_CLOCK (clock3,rate)
    dt(4) = (clock3-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)"  Test: SUPERLU2- ", 2," (",n,"x",n,") inverts in ", dt(4), ' seconds  Err=', avg_err  
    deallocate(mat,mat2)
    
    CONTAINS

! --------------------------------------------------------------------
    SUBROUTINE Superlucall2d  (a,n,ncc,acc,icc,ccc)     ! Invert matrix by SuperLU v53 method
! --------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER   :: n,ncc
    REAL(KIND=8):: a(n,n)
    integer :: m

    real ( kind = 8 ),     ALLOCATABLE ::acc (:),b(:,:)

    !real ( kind = 8 ) b(n)
    !real ( kind = 8 ) b2(n)
    logical cccadd
    integer nonzero,jika,isama
    integer ( kind = 4 ) , dimension ( n + 1 ) :: ccc 

    integer ( kind = 4 ) factors(8)
    integer ( kind = 4 ) i,j,ncc2,ccc2,jnext,k,l

    integer ( kind = 4 ) , ALLOCATABLE ::icc (:)

  
    integer ( kind = 4 ) info
    integer ( kind = 4 ) iopt
    integer ( kind = 4 ) ldb
    integer ( kind = 4 ) nrhs
  
    real ( kind = 8 ) dts(20)

! Convert to CCC matrix
    ncc2 = 0
    ccc2 = 0
    nonzero = 0
    jika = 0
    
    !allocate (icc(ncc),acc(ncc))
        
    isama = 0
    !DIR$ VECTOR ALIGNED
    do i=1,n
      
      !DIR$ VECTOR ALIGNED 
      do j = 1,n
      if (a(j,i).ne.0.d0) then

          ncc2 = ncc2 + 1
            acc(ncc2) = a(j,i)
            nonzero = nonzero + 1
            if (jika.eq.0) then 
                ccc(1) = i
                ccc2 = 1
                isama = i
                jika = jika + 1
            else
                jika = jika + 1
                if (isama.ne.i) then
                    ccc2 = ccc2+1
                    ccc(ccc2) = jika 
                endif    
                isama = i
            endif    
            
              

            icc(ncc2) = j
     
            
        endif  
      enddo  
     enddo 
    ccc(ccc2+1)=ncc+1  
 
    nrhs = n
    ldb = n
    !allocate (b(n,n))
    a = 0.d0
    do i = 1, n
        a(i,i) = 1.d0
    end do

!
!  Factor the matrix.
!
    iopt = 1
    call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc,ccc, a, ldb, factors, info )

    if ( info /= 0 ) then
    write ( *, '(a)' ) '  Factorization failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
    end if

!
!  Solve the factored system.
!
    iopt = 2
    call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc,ccc, a, ldb, factors, info )

    if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Backsolve failed'
    write ( *, '(a,i4)' ) '  INFO = ', info
    stop 1
    end if
  
!
!  Free memory.
!

    iopt = 3
    call c_fortran_dgssv ( iopt, n, ncc, nrhs, acc, icc, ccc, a, ldb, factors, info )
  
    END SUBROUTINE Superlucall2D
         
    subroutine BLASLAPACK_inverse(a,n)
!   variables for BLAS
    real*8, dimension(:), allocatable :: a(:,:) ! random numbers to invert and check 2d
    integer :: n,nzero
    real*8, dimension(:), allocatable :: work(:) 
    INTEGER, dimension(:), allocatable ::  ipiv(:)
    INTEGER :: info, lwork, ILAENV
    !  LAPACK BLAS
    
    !WRITE(*,*) "BLAS it with DGTRI DGTRF"
    
    lwork = n * ILAENV( 1, 'DGETRI', ' ', n, -1, -1, -1 )
    ALLOCATE ( work(lwork) )
    allocate (ipiv(n))
    CALL DGETRF( n, n, a, n, ipiv, info )
    CALL DGETRI( n, a, n, ipiv, work, lwork, info )
    DEALLOCATE ( work,ipiv )
    end
    
    end program benchinverse

!-------------------------------------------------------------------------------
    
    subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
!*****************************************************************************80
!
!! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
!
!  Discussion:
!
!    This is essentially a transposition operation.  
!
!    It is NOT an in-place algorithm.
!
!    This routine transposes a matrix stored in a, ja, ia format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:

!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, indicates whether or not to fill the values of the
!    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use
!                call 99998 (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
!        for any other normal usage, enter ipos=1.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
!    Compressed Sparse Column format.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iao(n+1)
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
!
!  Compute lengths of rows of A'.
!
  iao(1:n+1) = 0

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k) + 1
      iao(j) = iao(j) + 1
    end do
  end do
!
!  Compute pointers from lengths.
!
  iao(1) = ipos
  do i = 1, n
    iao(i+1) = iao(i) + iao(i+1)
  end do
!
!  Do the actual copying.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if ( job == 1 ) then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next + 1
    end do
  end do
!
!  Reshift IAO and leave.
!
  do i = n, 1, -1
    iao(i+1) = iao(i)
  end do
  iao(1) = ipos

  return
    end
    
!------------------------------------------------------------------------------
    subroutine csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the dimension of the DNS array.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) ndns

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dns(ndns,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nrow
  
  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end
    
!-------------------------------------------------------------------------------
    subroutine dump ( n, a, ja, ia, iout )

!*****************************************************************************80
!
!! DUMP writes the matrix to a file.
!
!  Discussion:
!
!    This routine writes the matrix to a file, one row at a time in a nice 
!    readable format.  This is a simple routine which is useful for debugging.
!
!    The output unit iout will have written in it the matrix in
!    one of two possible formats (depending on the max number of
!    elements per row. the values are output with only two digits
!    of accuracy (D9.2).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer ( kind = 4 ) IOUT, the FORTRAN output unit number.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) maxr
!
!  Select mode horizontal or vertical.
!
  maxr = 0
  do i = 1, n
    maxr = max ( maxr, ia(i+1) - ia(i) )
  end do

  if ( maxr <= 8 ) then
!
!  Able to print one row across line.
!
    do i = 1, n
      write(iout,100) i
      k1 = ia(i)
      k2 = ia(i+1) - 1
      write (iout,101) ja(k1:k2)
      write (iout,102) a(k1:k2)
    end do

  else
!
!  Unable to print one row acros line.  Do three items at a time acros line.
!
     do i = 1, n
        write(iout,200) i
        k1 = ia(i)
        k2 = ia(i+1) - 1
        write (iout,201) (ja(k),a(k), k = k1, k2)
    end do

  end if

 100  format(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  format(' col:',8(i5,6h     :))
 102  format(' val:',8(E9.2,2h :) )
 200  format(1h ,31(1h-),' row',i3,1x,31(1h-),/ &
         3('  columns :   values   *') )
 201  format(3(1h ,i5,6h    : ,D9.2,3h  *) )
  return
end
    
    
!-------------------------------------------------------------------------------   
    subroutine makesym(a,n) !from lower
    implicit none
    real*8  a(n,n) ! random numbers to invert and check 2d  
    integer n,i,j,k
       
    !create symmetry if needed, here symmetry is input
    do i=1,n
        !    k=k-1
        do j=i,n  
         a(j,i) = a(i,j)    
        end do
    enddo

 
    end
!-------------------------------------------------------------------------------  
    subroutine dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
!*****************************************************************************80
!
!! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine does not check whether an element is small.  It considers 
!    that A(I,J) is zero only if it is exactly equal to zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!    
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NZMAX, the maximum number of nonzero elements 
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0 means normal return;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) ndns
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dns(ndns,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nzmax

  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow

    do j = 1, ncol

      if ( dns(i,j) /= 0.0D+00 ) then

        if ( nzmax < next ) then
          ierr = i
          return
        end if

        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1

      end if

    end do

    ia(i+1) = next

  end do

  return
    end

!-------------------------------------------------------------------------------
    subroutine matrf2 ( m, n, c, index, alpha, nn, nz, a, snr, rnr, fejlm )

!*****************************************************************************80
!
!! MATRF2 generates sparse (rectangular or square) matrices.
!
!  Discussion:
!
!    The dimensions of the matrix and the average number of nonzero
!    elements per row can be specified by the user. Moreover, the user
!    can also change the sparsity pattern and the condition number of the
!    matrix. The non-zero elements of the desired matrix will be
!    accumulated (in an arbitrary order) in the first NZ positions of
!    array A. The column and the row numbers of the non-zero element
!    stored in A(I), I = 1,...,NZ, will be found in SNR(I) and RNR(I),
!    respectively. The matrix generated by this routine is of the
!    class F(M,N,C,R,ALPHA) (see reference).
!
!    If A is the sparse matrix of type F(M,N,C,R,ALPHA), then
!
!      min |A(i,j)| = 1/ALPHA,
!
!      max |A(i,j)| = max ( INDEX*N - N, 10*ALPHA ).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski
!
!  Reference:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski,
!    A testing Scheme for Subroutines Solving Large Linear Problems,
!    Computers and Chemistry, 
!    Volume 5, Number 2-3, pages 91-100, 1981.
!
!  Parameters:
!
!   INPUT PARAMETERS
!
!   M    - integer ( kind = 4 ). The number of rows in the desired matrix.
!          N < M+1 < 9000001 must be specified.
!
!
!   N    - integer ( kind = 4 ). The number of columns in the desired matrix.
!          21 < N < 9000001 must be specified.
!
!   C    - integer ( kind = 4 ). The sparsity pattern can be changed by means of this
!          parameter.  10 < C < N-10  must be specified.
!
!   INDEX - integer ( kind = 4 ).  The average number of non-zero elements per row in
!           the matrix will be equal to INDEX.
!           1 < INDEX < N-C-8 must be specified.
!
!   ALPHA - Real. The condition number of the matrix can be changed
!           BY THIS PARAMETER. ALPHA > 0.0 MUST BE SPECIFIED.
!           If ALPHA is approximately equal to 1.0 then the generated
!           matrix is well-conditioned. Large values of ALPHA will
!           usually produce ill-conditioned matrices. Note that no
!           round-off errors during the computations in this routine
!           are made if ALPHA = 2**I (where I is an arbitrary integer ( kind = 4 )
!           which produces numbers in the machine range).
!
!   Input, integer ( kind = 4 ) NN, the length of arrays A, RNR, and SNR.
!   INDEX*M+109 < NN < 9000001 must be specified.
!
!   Output, integer ( kind = 4 ) NZ, the number of nonzero elements in the matrix.
!
!   Output, real A(NN), the nonzero elements of the matrix,
!   accumulated in the first NZ locations of array A.
!
!   Output, integer ( kind = 4 ) SNR(NN), the column number of the non-zero element
!   kept in A(I), I = 1,...NZ.
!
!   Output, integer ( kind = 4 ) RNR(NN), the row number of the non-zero element
!   kept in A(I).
!
!   Output, integer ( kind = 4 ) FEJLM, error indicator.
!   0, indicates that the call is successful.
!   1, N is out of range.
!   2, M is out of range.
!   3, C is out of range.
!   4, INDEX is out of range.
!   5, NN is out of range.
!   7, ALPHA is out of range.
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(nn)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha1
  integer ( kind = 4 ) c
  integer ( kind = 4 ) fejlm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) nz1
  integer ( kind = 4 ) rnr(nn)
  integer ( kind = 4 ) rr1
  integer ( kind = 4 ) rr2
  integer ( kind = 4 ) rr3
  integer ( kind = 4 ) snr(nn)

  m1 = m
  fejlm = 0
  nz1 = index * m + 110
  k = 1
  alpha1 = alpha
  index1 = index - 1
!
!  Check the parameters.
!
  if ( n < 22 ) then
    fejlm = 1
    return
  end if

  if ( 9000000 < n ) then
    fejlm = 1
    return
  end if

  if ( n < m ) then
    fejlm = 2
    return
  end if

  if ( 9000000 < m ) then
    fejlm = 2
    return
  end if

  if ( c < 11 ) then
    fejlm = 3
    return
  end if

  if ( n-c < 11 ) then
    fejlm = 3
    return
  end if

  if ( index < 1 ) then
    fejlm = 4
    return
  end if

  if ( n - c - index < 9 ) then
    fejlm = 4
    return
  end if

  if ( nn < nz1 ) then
    fejlm = 5
    write(*,*) 'nn < nz1'
    return
  end if

  !if ( 9000000 < nn ) then
  !  write(*,*) '9000000 < nn'
  !  fejlm = 5
  !  return
  !end if

  if ( alpha <= 0.0D+00 ) then
    fejlm = 6
    return
  end if
!
!  End of the error check.  Begin to generate the non-zero elements of
!  the required matrix.
!
    a(1:n) = 1.0D+00

    do i = 1, n
    snr(i) = i
    end do

    do i = 1, n
    rnr(i) = i
    end do

    nz = n
    j1 = 1

    do j = 1, index1

    j1 = -j1

    do i = 1, n

      a(nz+i) = real ( j1 * j * i, kind = 8 )

      if ( i + c + j - 1 <= n ) then
        snr(nz+i) = i + c + j - 1
      end if

      if ( i + c + j - 1 > n ) then
        snr(nz+i) = c + i + j - 1 - n
      end if

      rnr(nz+i) = i

    end do

    nz = nz + n

    end do

    rr1 = 10
    rr2 = nz
    rr3 = 1

    do

    do i = 1, rr1
      a(rr2+i) = alpha * real ( i, kind = 8 )
      snr(rr2+i) = n - rr1 + i
      rnr(rr2+i) = rr3
    end do

    if ( rr1 == 1 ) then
      exit
    end if

    rr2 = rr2 + rr1
    rr1 = rr1 - 1
    rr3 = rr3 + 1

    end do

    nz = nz + 55

    do

    m1 = m1 - n
    alpha = 1.0D+00 / alpha

    if ( m1 <= 0 ) then
      exit
    end if

    n2 = k * n

    if ( n <= m1 ) then
      m2 = n
    end if

    if ( m1 < n ) then
      m2 = m1
    end if

    do i = 1, m2
      a(nz+i) = alpha * real ( k + 1, kind = 8 )
      snr(nz+i) = i
      rnr(nz+i) = n2 + i
    end do

    nz = nz + m2

    j1 = 1

    do j = 1, index1
 
      j1 = -j1

      do i = 1, m2

        a(nz+i) = alpha * real ( j * j1, kind = 8 ) &
          * ( real ( ( k + 1 ) * i, kind = 8 ) + 1.0D+00 )

        if ( i + c + j - 1 <= n ) then
          snr(nz+i) = i + c + j - 1
        end if

        if ( n < i + c + j - 1 ) then
          snr(nz+i) = c + i + j - 1 - n
        end if

        rnr(nz+i) = n2 + i

      end do
      nz = nz + m2
    end do

    k = k + 1
    end do


    alpha = 1.0D+00 / alpha1
    rr1 = 1
    rr2 = nz

    do

    do i = 1, rr1
      a(rr2+i) = alpha * real ( rr1 + 1 - i, kind = 8 )
      snr(rr2+i) = i
      rnr(rr2+i) = n - 10 + rr1
    end do

    if ( rr1 == 10 ) then
      exit
    end if

    rr2 = rr2 + rr1
    rr1 = rr1 + 1

    end do

    nz = nz + 55
    alpha = alpha1

    return
    end
!=====
    subroutine matrix_check_sym (a,n, symmetry, sparse,nzero)
    ! for check sparse and symmetry
    implicit none
    real(kind = 8)   :: a(n,n) ! random numbers to invert and check 2d
    integer(kind =4) :: n,t(8),nzero
    logical sparse,symmetry
    integer(kind =4) :: i,j,k
    real(kind = 8) sparseness

    ! Check sparse
    nzero = 0
    sparse = .false.
    do i=1,n
        do j=1,n  
         if (a(i,j).ne.0.d0) then
             nzero = nzero+1
         endif    
        enddo
    enddo
    if (nzero.eq.n*n) then 
        sparse = .false.
        write (*,*) " Matrix is full banded "
    else
        sparse = .true.
        sparseness = ((n*n) - NZERO)
        sparseness = sparseness / (n*n)
        write (*,*) " Matrix has zero values    : ", n*n - NZERO        
        write (*,*) " Matrix has nonzero values : ", NZERO
        write (*,*) " Matrix has total values   : ", n*n
        write (*,*) " Matrix has sparse about   : ", sparseness



    endif    
     
    ! Check symmetry
    k = 0
    symmetry = .false.
    do i=1,n
        do j=1,n  
         if (a(i,j).eq.a(j,i)) then
             k=k+1
             symmetry = .true. !symmetry  as output check
        else
             symmetry = .false. !symmetry  as output check
             exit
             exit
         endif    
        enddo
    enddo
	if (symmetry) then 
 	write (*,*) " Matrix is symmetric"
	else
	write (*,*) " Matrix is unsymmetric"
	endif
    
end subroutine
 
!   Intel DSS MKL
            

    