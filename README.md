# benchinverse
Benchmark for Solving Dense Matrix or Sparse Matrix
Version 0.0.1 MKL version int32 or integer*4 version fortran
uploaded 08/01/2023

how to compile?
1. install intel mkl for MKL version
2. download superlu first in  other  folder
3. extract SRC in superlu to separate folder.  also c_fortran_dgssv.c  c_fortran_sgssv.c  c_fortran_cgssv.c  c_fortran_zgssv.c 
   in FORTRAN folder must be copied in the same folder. 
   build the static in that folder
    
    example for windows intel Oneapi (c and Fortran) via Visual OneAPI command prompt x64
    
    icl *.c /c  /QaxCORE-AVX2 -DF77_CALL_C=UPCASE /O3

    example for windows MSVC via Visual OneAPI command prompt x64. Upcase by intel fortran MSVC (history: MS fortran -> Digital Compaq Fortran -> Intel Fortran) 
    
    icl *.c /c  /QaxCORE-AVX2 -DF77_CALL_C=UPCASE /O3
    lib *.obj /OUT:libsuperlu.lib 
 
    example for windows gcc equation.com 
    
    gcc -c -Ofast -ffast-math -march=core-avx2 *.c  (use batch file instead)
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

    benchinverse.exe           --> 22x22 full banded unsymetric
    benchinverse.exe 2000  --> 2000x2000 full banded unsymetric problem
    benchinverse.exe 2000  sym --> 2000x2000 full banded symetric problem
    benchinverse.exe 2000  unsym --> 2000x2000 full banded unsymetric problem
    benchinverse.exe 2000  sym  sparse --> 2000x2000 full banded symetric sparse problem
    benchinverse.exe 2000  unsym  sparse --> 2000x2000 full banded unsymetric sparse problem

Note : 
1. UPCASE calling is needed in intel fortran windows -DF77_CALL_C=UPCASE
 2. ADD_ calling is using for others like gfortran -DF77_CALL_C=ADD_
 3. MKL need Intel compiler, even MSVC will gave error in DSS.
 4. Benchinverse_nonMKL will not use DSS and PARDISO
