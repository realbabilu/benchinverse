cd superlu
gcc_win_libsuperlu_create.bat
cd ..
copy superlu\libsuperlu.a
gfortran  libsuperlu.a libopenblas.dll.a -o benchinverse_nonMKL.exe 