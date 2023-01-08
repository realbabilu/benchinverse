# benchinverse
Benchmark for Solving Dense Matrix or Sparse Matrix
Version 0.0.1 MKL and nonMKL version int32 or integer*4 version fortran
uploaded 08/01/2023
This will inverse the matrix and reinverse matrix again to check the different from original matrix

how to compile?
1. install intel mkl and intel compiler for MKL version, OR

   intell gcc equation.com for non MKL, need fortran compiler, c compiler, blas-lapack (openblas,acml, etc), or mys or cygwin or linux compiler.
   gcc equation.com will need libgcc_s_seh-1.dll libgfortran-5.dll libopenblas.dll libquadmath-0.dll libwinpthread-1.dll
   
2. download superlu first in  other  folder
3. extract SRC in superlu to separate folder.  also c_fortran_dgssv.c  c_fortran_sgssv.c  c_fortran_cgssv.c  c_fortran_zgssv.c 
   in FORTRAN folder must be copied in the same folder. 
   build the static in that folder
    
    example for windows intel Oneapi (c and Fortran) via Visual OneAPI command prompt x64
    
    icl *.c /c  /QaxCORE-AVX2 -DF77_CALL_C=UPCASE /O3

    example for windows MSVC via Visual OneAPI command prompt x64. Upcase by intel fortran MSVC (history: MS fortran -> Digital Compaq Fortran -> Intel Fortran) 
    
    icl *.c /c  /QaxCORE-AVX2 -DF77_CALL_C=UPCASE /O3
    lib *.obj /OUT:libsuperlu.lib 
    
    example for windows gcc equation.com (use batch file instead)
    
    gcc -c -Ofast -ffast-math -march:core-avx2 *.c
    ar.exe qc libsuperlu.a  *.o
    ranlib.exe libsuperlu.a
    
4   copy the lib (a or lib) to the folder of benchinverse.f90 
5.  Create exe
     in windows intel oneapi
     
     ifort benchinverse.f90 superlu.lib -DF77_CALL_C=UPCASE  /Qmkl /Ox /Qpar /QaxCORE-AVX2  

    in windows gcc equation.com with openblas
   
    gfortran benchinverse_nonMKL.f90  libsuperlu.a libopenblas.dll.a -o benchinverse_nonMKL.exe

6. how to use it?
    no argument --> default 22x22 random 
    1st argument --> size n x n
    2nd argument  -->  sym or unsym
    3rd argument  -->  sparse or full banded
    
    22x22 full banded unsymetric example 
    
    benchinverse.exe       
    
    2000x2000 full banded unsymetric problem example
    
    benchinverse.exe 2000 
    
    3000x3000 full banded symetric problem example
  
    benchinverse.exe 3000  sym 
    
    4000x4000 full banded unsymetric problem example
    
    benchinverse.exe 4000  unsym 
    
    5000x5000 symetric sparse problem example
    
    benchinverse.exe 5000  sym  sparse  
    
    6000x6000 unsymetric sparse problem example
    
    benchinverse.exe 6000  unsym  sparse 

Note : 
1. UPCASE calling is needed in intel fortran windows -DF77_CALL_C=UPCASE
 2. ADD_ calling is using for others like gfortran -DF77_CALL_C=ADD_
 3. MKL need Intel compiler, even MSVC will gave error in DSS.
 4. Benchinverse_nonMKL will not use DSS and PARDISO
