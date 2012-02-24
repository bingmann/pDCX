AC_PREREQ(2.67)
AC_INIT(pDCX, 0.0.0)
AC_CONFIG_SRCDIR([src/dc.cpp])
AC_CONFIG_AUX_DIR(acscripts)

AM_INIT_AUTOMAKE(foreign)

# set debug info flag if no optimization flags are set.

if test "$CFLAGS" == ""; then
    CFLAGS="-g"
fi
if test "$CXXFLAGS" == ""; then
    CXXFLAGS="-g"
fi

# enable full optimization by configure switch

AC_ARG_ENABLE(optimize,
	AS_HELP_STRING([--enable-optimize],
		       [Build with full optimization @<:@default=no@:>@]),
        [ case "${enableval}" in
          yes)
                CFLAGS="$CFLAGS -O3 -fomit-frame-pointer -march=native";
                CXXFLAGS="$CXXFLAGS -O3 -fomit-frame-pointer -march=native";
          	;;
          no) ;;
          *) AC_MSG_ERROR(bad value ${enableval} for --enable-optimize) ;;
          esac ],
        [ optimize=false ])

# enable compilation with assertions by configure switch

AC_ARG_ENABLE(ndebug,
	AS_HELP_STRING([--enable-ndebug],
		       [Build without assertions @<:@default=no@:>@]),
        [ case "${enableval}" in
          yes)
                CFLAGS="$CFLAGS -DNDEBUG";
                CXXFLAGS="$CXXFLAGS -DNDEBUG";
          	;;
          no) ;;
          *) AC_MSG_ERROR(bad value ${enableval} for --enable-ndebug) ;;
          esac ],
        [ ndebug=false ])

# find MPI and configure build environment

AX_PROG_CXX_MPI(true,[],[ AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.]) ])
AC_LANG([C++])

# Checks for programs.

AC_SYS_LARGEFILE

# construct Makefiles

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
