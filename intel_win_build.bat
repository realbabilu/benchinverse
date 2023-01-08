cd superlu
intelwin_libsuperlu_create.bat
cd ..
copy superlu\libsuperlu.lib
ifort benchinverse.f90 libsuperlu.lib -DF77_CALL_C=UPCASE  /Qmkl /Ox /Qpar