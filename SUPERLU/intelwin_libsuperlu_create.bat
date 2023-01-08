
icl *.c /c /Qmkl /QaxCORE-AVX2 -DF77_CALL_C=UPCASE /Ox /Qpar
lib *.obj /OUT:libsuperlu.lib
del *.obj