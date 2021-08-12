AC_INIT([Polly], 0.1.0)

## Use pkg-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: pkg-config(1)
AC_PATH_PROG([PKG_CONFIG], [pkg-config])
## If pkg-config was found, let's use it
if test "${PKG_CONFIG}" != ""; then
    # If pkg-config has ZLIB information, proceed
    if test "${PKG_CONFIG} --exists 'zlib' --print-errors" != ""; then
        ZLIB_CFLAGS=`${PKG_CONFIG} --cflags zlib`
        ZLIB_LIBS=`${PKG_CONFIG} --libs zlib`
    else
        AC_MSG_ERROR([ZLIB not found, is zlib installed?])
    fi
    # If pkg-config has BZIP2 information, proceed
    if test "${PKG_CONFIG} --exists 'bzip2' --print-errors" != ""; then
        BZIP2_CFLAGS=`${PKG_CONFIG} --cflags bzip2`
        BZIP2_LIBS=`${PKG_CONFIG} --libs bzip2`
    else
        AC_MSG_ERROR([BZIP2 not found, is bzip2 installed?])
    fi
    # If pkg-config has HTS information, proceed
    if test "${PKG_CONFIG} --exists 'htslib' --print-errors" != ""; then
        HTS_CFLAGS=`${PKG_CONFIG} --cflags htslib`
        HTS_LIBS=`${PKG_CONFIG} --libs htslib`
    else
        AC_MSG_ERROR([htslib not found, is htslib installed?])
    fi
else
    AC_MSG_ERROR([pkg-config not found, is PKG installed?])
fi

# Now substitute these variables in Makevars.in to create src/Makevars
AC_SUBST(ZLIB_CFLAGS)
AC_SUBST(ZLIB_LIBS)
AC_SUBST(BZIP2_CFLAGS)
AC_SUBST(BZIP2_LIBS)
AC_SUBST(HTS_CFLAGS)
AC_SUBST(HTS_LIBS)

AC_OUTPUT(src/Makevars)
