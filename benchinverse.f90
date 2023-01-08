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
    include 'mkl_dss.f90' ! Include the standard DSS "header file."

    program benchinverse
    
    implicit none
    include 'mkl_pardiso.fi' !for PARADISO
    
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
    
    ! VAR DSS
    !integer :: error
    
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
    
    !hbwrite
    character ( len = 8 ) :: key != 'RUA_32'
    character ( len = 3 ) :: mxtype != 'RUA' OR !PUA 
    !integer ( kind = 4 ) :: ncol != 32
    !integer ( kind = 4 ) :: nrow != 32
    character ( len = 80 ) :: output_file! = 'rua_32_header.txt'
    integer ( kind = 4 ) :: output_unit
    character ( len = 72 ) :: title != '1Real unsymmetric assembled matrix based on IBM32'

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
            stop 'Cannot allocate memory for the main, Argument 1 is size matrix, Argument 2 is sym/unsym, Argument 3 is sparse/dense ' 
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
    
 !   Part2 DSS MKL
    ! Set the problem to be solved for DSS
    n = bigsize
     open (11,file='matrix_csr.dat',form='unformatted')
     read(11) nz2,bigsize2
     close(11) 
     
     if ((nz2.eq.nz1).and.(bigsize2.eq.bigsize)) then
        allocate(a(nz2),ja(nz2),ia(bigsize+1))
        a = 0.d0
        ja = 0
        ia = 0
        open (11,file='matrix_csr.dat',form='unformatted')
        read(11) nz2,bigsize,a,ja,ia
        close(11)
        nz = nz2
     else
         write (*,*) bigsize2,nz2
         stop '  Matrix is mismatch - DSS1'
     endif
     
    allocate (mat2 (n,n)) !as solution
    mat2 = 0.d0
        
!   2D matrix
    CALL SYSTEM_CLOCK (clock1,rate)  ! get benchmark (n) start time
    write(*,*) '  STAGE 1 - Inverse'

    !CALL DSSMKL_inverse2D (n,columns,rowindex,values,nzero,symmetry,solution,benchmarkmode) !solution as input and output
	 CALL DSSMKL_inverse2D (n,ja,ia,a,nz,symmetry,mat2,benchmarkmode) !solution as input and output
	
    pass(2,1) = mat2(1,1)
    pass(2,2) = mat2(n,n)	
    
    !nzero change here
    deallocate(a,ja,ia)

    CALL SYSTEM_CLOCK (clock2,rate) 
    dt(2) = (clock2-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.
    WRITE (*,42)" Test2-1: MKLDSS - ", 1," (",n,"x",n,") inverts in ", dt(2), ' seconds' 

    write(*,*) '  STAGE 2 - ReInverse'
    call matrix_check_sym (mat2,bigsize, symmetry, sparse,nz) !get nz
    allocate(a(nz),ja(nz),ia(bigsize+1))
     a = 0.d0
     ja = 0
     ia = 0
     call dnscsr ( bigsize, bigsize, nz, mat2, bigsize , a, ja, ia, ierr )   
     mat2 = 0.d0
    
    CALL DSSMKL_inverse2D (n,ja,ia,a,nz,symmetry,mat2,benchmarkmode) !solution as input and output
	deallocate(a,ja,ia)
    
    pass(3,1) = mat2(1,1)
    pass(3,2) = mat2(n,n)	

    CALL SYSTEM_CLOCK (clock3,rate) 
    dt(2) = (clock3-clock2)/DBLE(rate)  ! get benchmark (n) elapsed sec.
    WRITE (*,42)" Test2-2: MKLDSS - ", 1," (",n,"x",n,") inverts in ", dt(2), ' seconds' 

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
        stop ' Size matrix is mismatch - DSS2'
    endif    
    
    avg_err = SUM(ABS(mat-mat2))/(n*n)    ! invert err

    CALL SYSTEM_CLOCK (clock3,rate)
    dt(2) = (clock3-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)"  Test2: MKLDSS - ", 2," (",n,"x",n,") inverts in ", dt(2), ' seconds  Err=', avg_err  
    deallocate(mat,mat2)    

! Part3 PARDISO
    ! Set the problem to be solved for DSS
     n = bigsize
     open (11,file='matrix_csr.dat',form='unformatted')
     read(11) nz2,bigsize2
     close(11) 
     
     if ((nz2.eq.nz1).and.(bigsize2.eq.bigsize)) then
        allocate(a(nz2),ja(nz2),ia(bigsize+1))
        a = 0.d0
        ja = 0
        ia = 0
        open (11,file='matrix_csr.dat',form='unformatted')
        read(11) nz2,bigsize,a,ja,ia
        close(11)
        nz = nz2
     else
         write (*,*) bigsize2,nz2
         stop '  Matrix is mismatch - PARDISO 1'
     endif
     
    allocate (mat2 (n,n)) !as solution
    mat2 = 0.d0

    !pass(1,1) same 
    !pass(1,2) same	

    CALL SYSTEM_CLOCK (clock1,rate)  ! get benchmark (n) start time
    write(*,*) '  STAGE 1 - Inverse'

!   2D matrix
    !CALL PARDISOMKL_inverse2D (n,columns,rowindex,values,nzero,symmetry,solution,benchmarkmode) !solution as input and output
    CALL PARDISOMKL_inverse2D (n,ja,ia,a,nz,symmetry,mat2,benchmarkmode) !solution as input and output

    pass(2,1) = mat2(1,1)
    pass(2,2) = mat2(n,n)	
    
    CALL SYSTEM_CLOCK (clock2,rate) 
    dt(3) = (clock2-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec
    WRITE (*,42)" Test3-1: PARDISO- ", 1," (",n,"x",n,") inverts in ", dt(3), ' seconds'  

    
    !nzero change here
    deallocate(a,ja,ia)
    write(*,*) '  STAGE 2 - ReInverse'

    call matrix_check_sym (mat2,bigsize, symmetry, sparse,nz) !get nz
    allocate(a(nz),ja(nz),ia(bigsize+1))
     a = 0.d0
     ja = 0
     ia = 0
     call dnscsr ( bigsize, bigsize, nz, mat2, bigsize , a, ja, ia, ierr )   
     mat2 = 0.d0

    CALL PARDISOMKL_inverse2D (n,ja,ia,a,nz,symmetry,mat2,benchmarkmode) !solution as input and output
	deallocate (a,ja,ia)
    pass(3,1) = mat2(1,1)
    pass(3,2) = mat2(n,n)	
    
    
    CALL SYSTEM_CLOCK (clock3,rate) 
    dt(3) = (clock3-clock2)/DBLE(rate)  ! get benchmark (n) elapsed sec
    WRITE (*,42)" Test3-2: PARDISO- ", 1," (",n,"x",n,") inverts in ", dt(3), ' seconds'  

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
        stop ' Size matrix is mismatch - PARDISO 2'
    endif    
    
    avg_err = SUM(ABS(mat-mat2))/(n*n)    ! invert err

    CALL SYSTEM_CLOCK (clock2,rate)
    dt(3) = (clock3-clock1)/DBLE(rate)  ! get benchmark (n) elapsed sec.   
    WRITE (*,42)"  Test3: PARDISO- ", 2," (",n,"x",n,") inverts in ", dt(3), ' seconds  Err=', avg_err  
    deallocate(mat2,mat)
    
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
        
    
! This Pardiso intel MKL
    
    subroutine PARDiSOMKL_inverse2D(n,columns,rowindex,values,nzero,symmetry,solution,benchmarkmode)
    REAL(KIND=8), dimension(:), allocatable :: solution(:,:) ! random numbers to invert and check 2d
    integer(kind=4) :: n,nzero
    logical symmetry,benchmarkmode
! Define the data arrays and the solution and rhs vectors.
    INTEGER(kind=4), ALLOCATABLE :: columns( : )
    INTEGER (kind=4), ALLOCATABLE :: rowIndex( : )
    REAL(KIND=8), ALLOCATABLE :: values( : )
    
    INTEGER(kind=4) :: error
    INTEGER(kind=4) :: nCols
    INTEGER(kind=4) :: nRhs
    INTEGER(kind=4) :: nRows
    REAL(KIND=8), ALLOCATABLE :: rhs( :,: )

    !   pardiso var
    TYPE(MKL_PARDISO_HANDLE) pt(64)
!.. All other variables
    INTEGER (kind=4) :: maxfct, mnum, mtype, phase, msglvl
    INTEGER (kind=4) :: iparm(64),perm(1)
    INTEGER (kind=4) :: idum(1)
    REAL(kind=8)  :: ddum(1)    
    
!   pt Solver internal data address pointer 0 Must be initialized with zeros and never be modified later 
!.. Initialize the internal solver memory pointer. This is only necessary for the FIRST call of the PARDISO solver.
        DO i = 1, 64
            pt(i)%DUMMY = 0
        END DO

!.. Setup PARDISO control parameter
!..
        DO i = 1, 64
           iparm(i) = 0
        END DO

        !iparm[64] This array is used to pass various parameters to IntelÂ® oneAPI Math Kernel Library
        ! PARDISO and to return some useful information after execution of the solver (see pardiso iparm Parameter for more details) 
        iparm(1) = 1 ! no solver default
        !If iparm[0] =0 IntelÂ® oneAPI Math Kernel Library PARDISO fillsiparm [1] through iparm [63] with default values and uses them. 
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(3) = 1 ! numbers of processors
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! solution on the first n components of x
        iparm(7) = 0 ! not in use
        iparm(8) = 9 ! numbers of iterative refinement steps
        iparm(9) = 0 ! not in use
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(12) = 0 ! not in use
        iparm(13) = 1 ! maximum weighted matching algorithm is ON
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(15) = 0 ! not in use
        iparm(16) = 0 ! not in use
        iparm(17) = 0 ! not in use
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations

        error = 0 ! initialize error flag
         msglvl = 1
	   if (benchmarkmode) msglvl = 0
	   ! print statistical information 1=on 0=off
        !n  !Number of equations in the sparse linear system A*X= B 	n>0         
        
        maxfct = 1 ! Maximal number of factors in memory >0 Generally used value is 1 
        mnum = 1 !The number of matrix (from 1 to maxfct) to solve; Generally used value is 1                     
        
        if (.not.symmetry ) mtype = 11 ! real unsymmetric
	    if (symmetry ) mtype = 1 ! real symmetric
	
!Matrix type 1 Real and structurally symmetric
!            2 Real and symmetric positive definite
!           -2 Real and symmetric indefinite
!            3 Complex and structurally symmetric
!            4 Complex and Hermitian positive definite
!           -4 Complex and Hermitian indefinite
!            6 Complex and symmetric matrix
!           11 Real and nonsymmetric matrix
!           13 Complex and nonsymmetric matrix 

        
! Set the problem to be solved.
    
    nRows = n
    nCols = n
    nRhs = n
        !nrhs
	        ! Number of right-hand sides that need to be solved for >=0
	        ! Generally used value is 1
            ! To obtain better IntelÂ® oneAPI Math Kernel Library PARDISO performance, during the numerical factorization phase you can provide the maximum number of right-hand sides, which can be used further during the solving phase. 
  
    
    !rowindex = ia [n]
              !rowIndex array in CSR3 format >=0 gives the index of the element in array a
              !that contains the first non-zero element from row i of A
              ! The last element ia (n) is taken to be equal to the number of non-zero elements in A
              ! Note: iparm [34] indicates whether row/column indexing starts from 1 or 0. 
    iparm(34) = 1
    

    !ja = columns
              ! columns array in CSR3 format >=0
	          ! The column indices for each row of A must be sorted in increasing order. For structurally symmetric matrices zero diagonal elements must be stored in a
              ! and ja. Zero diagonal elements should be stored for symmetric matrices, although they are not required. For symmetric matrices, the solver needs only the upper triangular part of the system.
              !Note: iparm [34] indicates whether row/column indexing starts from 1 or 0. 
    

              !A = Contains the non-zero elements of the coefficient matrix A*
              !The size of a is the same as that of ja , and the coefficient matrix can be either real or complex. 
              !The matrix must be stored in the 3-array variation of compressed sparse row (CSR3) format 
              ! with increasing values of ja for each row
    
    perm(1) = 0
        !perm[n]
	          !Holds the permutation vector of size n , specifies elements used for computing a partial solution, or specifies differing values of the input matrices for low rank update >=0 
    
    
    !rhs and solution matrix allocated already
    ! rhs or B
        ! b[n*nrhs]
	    ! Right-hand side vectors 
        ! On entry, contains the right-hand side vector/matrix B , which is placed contiguously in memory. The b[i+k*n]
        ! element must hold the i-th component of k-th right-hand side vector. Note that b is only accessed in the solution phase.
        ! On output, the array is replaced with the solution if iparm [5] =1. 
    
    !solution or x
        ! x [n*nrhs] Solution vectors 
        ! On output, if iparm [5] =0, contains solution vector/matrix X which is placed contiguously in memory. The 
        ! x[i+k*n] element must hold the i-th component of k-th solution vector. Note that x is only accessed in the solution phase. 
     
!.. Reordering and Symbolic Factorization, This step also allocates
!   all memory that is necessary for the factorization
        ! phase = Controls the execution of the solver For iparm[35]> 0
        ! phases 331, 332, and 333 perform a different decomposition. See the phase parameter of pardiso for details.
	    ! 11 Analysis
	    ! 12 Analysis, numerical factorization
        ! 13 Analysis, numerical factorization, solve
        ! 22 Numerical factorization
        ! 23 Numerical factorization, solve
        ! 33 Solve, iterative refinement
        ! 331 phase =33, but only forward substitution
        ! 332 phase =33, but only diagonal substitution
        ! 333 phase =33, but only backward substitution
        ! 0 Release internal memory for L and U of the matrix number mnum
        ! -1 Release all internal memory for all matrices         
    
    !error
        !	0 No error
        !   -1 Input inconsistent
        ! -2 Not enough memory
        ! -3 Reordering problem
        ! -4 Zero pivot, numerical factorization or iterative refinement problem
        ! -5 Unclassified (internal) error
        ! -6 Reordering failed (matrix types 11 and 13 only)
        ! -7 Diagonal matrix is singular
        ! -8 32-bit integer overflow problem
        ! -9 Not enough memory for OOC
        ! -10 Problems with opening OOC temporary files
        ! -11 Read/write problems with the OOC data file     
    
        phase = 11 ! only reordering and symbolic factorization        
        !    pardiso (pt, maxfct, mnum, mtype, phase,n,a     ,ia       ,ja       ,perm, nrhs, iparm, msglvl, b   , x   , error)
        CALL pardiso (pt, maxfct, mnum, mtype, phase,n,values, rowindex, columns, idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Reordering completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
        !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
        !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!.. Factorization.
        phase = 22 ! only factorization
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowindex, columns,idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Factorization completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
        
!.. Back substitution and iterative refinement
        iparm(8) = 2 ! max numbers of iterative refinement steps
        phase = 33 ! only factorization
        solution = 0.d0
        
        
!   matrix identity back to 1
        ALLOCATE(RHS(n,n))
        RHS = 0.d0 
        DO I=1,n
            RHS(I,I) = 1.d0
        ENDDO   
        
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowindex, columns, idum, nrhs, iparm, msglvl, rhs( : , 1 ), solution( : , 1 ), error)
        
!.. Termination and release of memory
        phase = -1 ! release internal memory
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum,  idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
        if (ALLOCATED(RHS)) DEALLOCATE(RHS)
    end subroutine

    
   subroutine DSSMKL_inverse2D (n,columns,rowindex,values,nzero,symmetry,solution,benchmarkmode)
      use mkl_dss
      implicit none
      real*8, dimension(:), allocatable :: solution(:,:) ! random numbers to invert and check 2d
      integer :: n,nzero
      logical symmetry,benchmarkmode
      integer :: error  
  ! Define the data arrays and the solution and rhs vectors.
      INTEGER, ALLOCATABLE :: columns( : )
      INTEGER :: nCols
      INTEGER :: nRhs
      INTEGER :: nRows
      REAL(KIND=8), ALLOCATABLE :: rhs( :,: )
      INTEGER, ALLOCATABLE :: rowIndex( : )
      REAL(KIND=8), ALLOCATABLE :: values( : )
      
      TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
      REAL(KIND=8),ALLOCATABLE::statOUt( : )
      CHARACTER*15 statIn
      INTEGER perm(1)
      INTEGER, PARAMETER :: bufLen = 20
      INTEGER buff(bufLen)
      
      nRows = n
      nCols = n 
      nRhs = n   
      perm(1) = 0
      solution = 0.d0
  
  !    WRITE(*,*) "Intel Direct Sparse Solver"
  ! Initialize the solver.
      error = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
      IF (error /= MKL_DSS_SUCCESS) GOTO 999
  ! Define the non-zero structure of the matrix.
      if (.not.symmetry) then
      error = DSS_DEFINE_STRUCTURE(handle, MKL_DSS_NON_SYMMETRIC, rowIndex, nRows, &
                                   nCols, columns, nzero)
      else
      error = DSS_DEFINE_STRUCTURE(handle, MKL_DSS_SYMMETRIC_STRUCTURE, rowIndex, nRows, &
                                   nCols, columns, nzero)
      endif
      IF (error /= MKL_DSS_SUCCESS) GOTO 999
  ! Reorder the matrix.
      error = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
      IF (error /= MKL_DSS_SUCCESS) GOTO 999
  ! Factor the matrix.
      error = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, values)
      IF (error /= MKL_DSS_SUCCESS) GOTO 999    

      ALLOCATE(rhs(nRows,NCols))
      !   matrix identity back to 1
      RHS = 0.0d0    
      DO I=1,n
       RHS(I,I) = 1.0d0
      ENDDO 
      
      ! Allocate the solution vector and solve the problem.
      error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, rhs( : , 1 ), nRhs, solution( : , 1 ))
      DEALLOCATE(rhs)
      
      IF (error /= MKL_DSS_SUCCESS) GOTO 999

  ! Print Out the determinant of the matrix (no statistics for a diagonal matrix)
      ALLOCATE(statOut(5))

      if (.not.benchmarkmode) then    
        statIn = 'Flops'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
        WRITE(*,"('Number of floating point operations  : '(5F10.3))") (statOut(1))
  	
   	    statIn = 'Peakmem'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
        WRITE(*,"('Total peak memory in kilobytes       : '(5F10.3))") (statOut(1))
  	
   	    statIn = 'Factormem'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
        WRITE(*,"('Permanent memory in kilobytes        : '(5F10.3))") (statOut(1))
  	
   	    statIn = 'Solvemem'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
        WRITE(*,"('Total double precision memory        : '(5F10.3))") (statOut(1))
  	
        IF(nRows < nzero) THEN
         statIn = 'determinant'
         call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
         error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
         IF (error /= MKL_DSS_SUCCESS) GOTO 999
         WRITE(*,"('pow of determinant : '(5F10.3))") (statOut(1))
         WRITE(*,"('base of determinant: '(5F10.3))") (statOut(2))
         WRITE(*,"('Determinant        : '(5F10.3))") ((10**statOut(1))*statOut(2))
      	END IF

  	endif
  	
  ! Deallocate solver storage and various local arrays.
      error = DSS_DELETE(handle, MKL_DSS_DEFAULTS)
      IF (error /= MKL_DSS_SUCCESS) GOTO 999
      IF (ALLOCATED(statOut)) DEALLOCATE(statOut)
  
      GOTO 1000
  ! Print an error message and exit
  999 WRITE(*,*) "Solver returned error code", error
      STOP 1
  1000 CONTINUE
      
      end subroutine
   
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
subroutine prtmt ( nrow, ncol, a, ja, ia, rhs, guesol, title, key, type, &
  ifmt, job, iounit )

!*****************************************************************************80
!
!! PRTMT writes a matrix in Harwell-Boeing format into a file.
!
!  Discussion:
!
!    This routine assumes that the matrix is stored in CSC format
!    (Compressed Sparse Column format).
!    There is some limited functionality for right hand sides.
!
!    This code attempts to pack as many elements as possible per
!    80-character line.
!
!    This code attempts to avoid as much as possible to put
!    blanks in the formats that are written in the 4-line header
!    This is done for purely esthetical reasons since blanks
!    are ignored in format descriptors.
!
!    sparse formats for right hand sides and guesses not supported.
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
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NCOL+1), the matrix in CSC
!    Compressed Sparse Column format.
!
!    Input, real RHS(*), contains the right hand sides and optionally
!    the associated initial guesses and/or exact solutions
!    in this order.  See also GUESOL for details.   RHS will
!    be used only if 2 < JOB.  Only full storage for the right hand 
!    sides is supported.
!
! guesol = a 2-character string indicating whether an initial guess
!          (1-st character) and / or the exact solution (2-nd)
!          character) is provided with the right hand side.
!         if the first character of guesol is 'G' it means that an
!          an intial guess is provided for each right hand sides.
!          These are assumed to be appended to the right hand sides in
!          the array rhs.
!         if the second character of guesol is 'X' it means that an
!          exact solution is provided for each right hand side.
!          These are assumed to be appended to the right hand sides
!          and the initial guesses (if any) in the array rhs.
!
! title  = character*71 = title of matrix test ( character a*71 ).
! key    = character*8  = key of matrix
! type   = charatcer*3  = type of matrix.
!
! ifmt       = integer ( kind = 4 ) specifying the format chosen for the real values
!         to be output (i.e., for a, and for rhs-guess-sol if
!          applicable). the meaning of ifmt is as follows.
!        * if (ifmt < 100) then the E descriptor is used,
!           format Ed.m, in which the length (m) of the mantissa is
!           precisely the integer ( kind = 4 ) ifmt (and d = ifmt+6)
!        * if (ifmt > 100) then prtmt will use the
!           F- descriptor (format Fd.m) in which the length of the
!           mantissa (m) is the integer ( kind = 4 ) mod(ifmt,100) and the length
!           of the integer ( kind = 4 ) part is k = ifmt/100 (and d = k+m+2)
!          Thus  ifmt= 4   means  E10.4  +.xxxxD+ee    while
!                ifmt=104  means  F7.4   +x.xxxx
!                ifmt=205  means  F9.5   +xx.xxxxx
!          Note: formats for ja, and ia are internally computed.
!
! job       = integer ( kind = 4 ) to indicate whether matrix values and
!         a right hand side is available to be written
!          job = 1   write srtucture only, i.e., the arrays ja and ia.
!          job = 2   write matrix including values, i.e., a, ja, ia
!          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
!         job = nrhs+2 write matrix and nrhs successive right hand sides
!         Note that there cannot be any right hand side if the matrix
!         has no values. Also the initial guess and exact solutions when
!          provided are for each right hand side. For example if nrhs=2
!          and guesol='GX' there are 6 vectors to write.
!
!
! iounit = logical unit number where to write the matrix into.
!
! on return:
!
! the matrix a, ja, ia will be written in output unit iounit
! in the Harwell-Boeing format. Noe of the inputs is modofied.
!
  implicit none

  integer ( kind = 4 ) ncol

  real ( kind = 8 ) a(*)
  character ( len = 2 ) guesol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(ncol+1)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ifmt
  integer ( kind = 4 ) ihead
  integer ( kind = 4 ) indcrd
  character ( len = 16 ) indfmt
  integer ( kind = 4 ) iounit
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) job
  character ( len = 8 ) key
  integer ( kind = 4 ) len
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nperli
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) ptrcrd
  character ( len = 16 ) ptrfmt
  real ( kind = 8 ) rhs(*)
  integer ( kind = 4 ) rhscrd
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer ( kind = 4 ) totcrd
  character ( len = 3 ) type
  integer ( kind = 4 ) valcrd
  character ( len = 20 ) valfmt
!
!  Compute pointer format.
!
  nnz = ia(ncol+1) - 1
  len = int ( log10 ( 0.1D+00 + real ( nnz + 1, kind = 8 ) ) ) + 1
  nperli = 80 / len
  ptrcrd = ncol / nperli + 1

  if ( 9 < len ) then
     assign 101 to ix
  else
     assign 100 to ix
  end if

  write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
!
!  Compute the ROW index format.
!
  len = int ( log10 ( 0.1D+00 + real ( nrow, kind = 8 ) ) ) + 1
  nperli = min ( 80 / len, nnz )
  indcrd = ( nnz - 1 ) / nperli + 1
  write (indfmt,100) nperli,len
!
!  Compute values and RHS format (using the same for both).
!
  valcrd = 0
  rhscrd = 0
!
!  Skip this part if no values provided.
!
  if ( job <= 1 ) then
    go to 20
  end if

  if ( 100 <= ifmt ) then
     ihead = ifmt / 100
     ifmt = ifmt - 100 * ihead
     len = ihead + ifmt + 2
     nperli = 80 / len

     if ( len <= 9 ) then
        assign 102 to ix
     elseif ( ifmt <= 9 ) then
        assign 103 to ix
     else
        assign 104 to ix
     end if

     write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )

  else
     len = ifmt + 6
     nperli = 80 / len
!
!  Try to minimize the blanks in the format strings.
!
     if ( nperli <= 9 ) then
      if ( len <= 9 ) then
         assign 105 to ix
      else if ( ifmt <= 9 ) then
         assign 106 to ix
      else
         assign 107 to ix
      end if
   else
      if ( len <= 9 ) then
         assign 108 to ix
      else if ( ifmt <= 9 ) then
         assign 109 to ix
      else
           assign 110 to ix
        end if
     end if

     write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hE,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hE,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hE,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hE,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hE,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hE,i2,1h.,i2,1h) )

  end if
  valcrd = ( nnz - 1 ) / nperli + 1
  nrhs = job - 2

  if ( 1 <= nrhs ) then
     i = ( nrhs * nrow - 1 ) / nperli + 1
     rhscrd = i
     if ( guesol(1:1) == 'G' ) then
       rhscrd = rhscrd + i
     end if
     if ( guesol(2:2) == 'X' ) then
       rhscrd = rhscrd + i
     end if
     rhstyp = 'F' // guesol
  end if

 20   continue

  totcrd = ptrcrd + indcrd + valcrd + rhscrd
!
!  Write four line or five line header.
!
  write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd, &
       rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt

  if ( 1 <= nrhs ) then
    write (iounit,11) rhstyp, nrhs
  end if

 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i4)

  write(iounit,ptrfmt) ia(1:ncol+1)
  write(iounit,indfmt) ja(1:nnz)

  if ( job <= 1 ) then
    return
  end if

  write(iounit,valfmt) (a(i), i = 1, nnz)
  if ( job <= 2 ) then
    return
  end if
  len = nrow * nrhs
  next = 1
  iend = len
  write(iounit,valfmt) (rhs(i), i = next, iend)
!
!  Write initial guesses if available
!
  if ( guesol(1:1) == 'G' ) then
     next = next + len
     iend = iend + len
     write(iounit,valfmt) (rhs(i), i = next, iend)
  end if
!
!  Write exact solutions if available.
!
  if ( guesol(2:2) == 'X' ) then
     next = next + len
     iend = iend + len
     write(iounit,valfmt) (rhs(i), i = next, iend)
  end if

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
            

    