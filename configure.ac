AC_INIT([mgsa], 0.99.0)

dnl ******************************
dnl openmp checking
dnl ******************************

AX_OPENMP

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)

dnl and do substitution in the src/Makevars.in
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
