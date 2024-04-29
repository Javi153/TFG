#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
AC_PREREQ(2.59)
AC_INIT([CPSP], [2.4.6], [], [cpsp], [http://cpsp.informatik.uni-freiburg.de])
AC_CONFIG_AUX_DIR([.])
AC_CANONICAL_HOST


# backwards compatibility for autoconf >= 2.64
# PACKAGE_URL should be the fifth argument of AC_INIT
m4_ifndef([AC_PACKAGE_URL],
	[m4_define([AC_PACKAGE_URL], [http://www.bioinf.uni-freiburg.de/~mmann/])
	AC_DEFINE(PACKAGE_URL, "AC_PACKAGE_URL", [Package URL])
	AC_SUBST(PACKAGE_URL, AC_PACKAGE_URL)])

# special linker flags for cygwin etc. under MS Windows
lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac


# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS


AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (def=no)]
	),
  	debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE(DEBUG)
	CPPFLAGS="$CPPFLAGS -g -Wall -Werror -Wno-uninitialized -O0 "
else
	AC_DEFINE(NDEBUG)
	CPPFLAGS="$CPPFLAGS -O3 "
fi

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([cpsp], [misc/doxygen/doxygen.config], [doc])


AM_INIT_AUTOMAKE

AC_LANG([C++])
AC_DEFINE([__cplusplus])


###############################################################################
# BEGIN STATIC LINKING CHECK
###############################################################################

AC_ARG_ENABLE(static-linking,
    AC_HELP_STRING(
    	[--enable-static-linking],
	    [enables static linking of the executables (default=disabled)]
	),
    [enable_static_linking=$enableval],
    [enable_static_linking=no]
)

if test "$enable_static_linking" = "yes"; then
	LDFLAGS="-static $LDFLAGS"
fi
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT($enable_static_linking)


###############################################################################
# END STATIC LINKING CHECK
###############################################################################



# BIOINFORMATIC UTILITY LIBRARY - BIU
LIBS="-lbiu $LIBS"

# BIU library path support, if not installed in usual directories
AC_ARG_WITH(BIU,
	AC_HELP_STRING(
	    [--with-BIU=prefix],
	    [try this for the BIU-library prefix install directory]
	  ),
	  BIUPATHSET=1,
	  BIUPATHSET=0
)
if test  $BIUPATHSET = 1 ; then
		CPPFLAGS="$CPPFLAGS -I$with_BIU/include"
		LDFLAGS="-L$with_BIU/lib $LDFLAGS"
fi

# GECODE LIBRARIES
LIBS="-lgecodeminimodel -lgecodeint -lgecodesearch -lgecodekernel $LIBS"

# GECODE library path support, if not installed in usual directories
AC_ARG_WITH(gecode,
	AC_HELP_STRING(
	    [--with-gecode=prefix],
	    [try this for the Gecode-library prefix install directory]
	  ),
	  GECODEPATHSET=1,
	  GECODEPATHSET=0
)
if test  $GECODEPATHSET = 1 ; then
		CPPFLAGS="$CPPFLAGS -I$with_gecode/include"
		LDFLAGS="-L$with_gecode/lib $LDFLAGS"
fi



AC_SUBST(GLOBALCPPFLAGS, $CPPFLAGS )
AC_SUBST(GLOBALLDFLAGS, $LDFLAGS )

AC_LANG_CPLUSPLUS

# Checks for libraries and header files.


# check for stl headers
AC_CHECK_HEADERS([vector])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# check for biu headers
AC_MSG_CHECKING([for the BIU C++ library])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <biu/VirtualList.hh>]],[[return 0;]] )],
     [ 
       AC_MSG_RESULT([yes])
       BIUNOTFOUND=0;
     ],
     [
       AC_MSG_RESULT([no])
       BIUNOTFOUND=1;
     ]
)

# check for gecode headers
AC_MSG_CHECKING([for the Gecode C++ library])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <gecode/kernel.hh>]],[[return 0;]] )],
     [ 
       AC_MSG_RESULT([yes])
       GECODENOTFOUND=0;
     ],
     [
       AC_MSG_RESULT([no])
       GECODENOTFOUND=1;
     ]
)



DEPENDENCYNOTFOUND=0;

# error output if BIU-headers not found
if test $BIUNOTFOUND = 1; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE(the BIU C++ library is required for building CPSP-tools.)
	AC_MSG_NOTICE(it can be obtained from http://www.bioinf.uni-freiburg.de/SW/BIU/.)
	AC_MSG_NOTICE()
 	if test $BIUPATHSET = 1; then
		AC_MSG_NOTICE(cant find BIU in given path '$with_BIU'.)
 	else
		AC_MSG_NOTICE(if BIU is installed elsewhere use --with-BIU=prefix.)
 	fi
 	DEPENDENCYNOTFOUND=1;
fi

# error output if Gecode-headers not found
if test $GECODENOTFOUND = 1; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE(the Gecode C++ library is required for building CPSP-tools.)
	AC_MSG_NOTICE(it can be obtained from http://www.gecode.org.)
	AC_MSG_NOTICE()
 	if test $GECODEPATHSET = 1; then
		AC_MSG_NOTICE(cant find Gecode in given path '$with_gecode'.)
 	else
		AC_MSG_NOTICE(if Gecode is installed elsewhere use --with-gecode=prefix.)
 	fi
 	DEPENDENCYNOTFOUND=1;
fi



# error ABORT if one of the libraries was not found
if test $DEPENDENCYNOTFOUND = 1; then
	AC_MSG_NOTICE()
	AC_MSG_ERROR(library not found! Try './configure --help'.)
fi

AC_CONFIG_FILES([ Makefile src/doxygen.hh ])
AC_CONFIG_FILES([ src/Makefile src/cpsp/Makefile ])
AC_CONFIG_FILES([ src/cpsp/gecode/Makefile ])
AC_CONFIG_FILES([ src/cpsp/SideChain/Makefile ])
#AC_CONFIG_FILES([ src/cpsp/ilog/Makefile ])
AC_CONFIG_FILES([ src/bin/version.hh ])
AC_CONFIG_FILES([ src/bin/Makefile ])
AC_CONFIG_FILES([ misc/Makefile misc/doxygen/Makefile ])
AC_CONFIG_FILES([ tests/Makefile ])

AC_OUTPUT
