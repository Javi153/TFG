/* gecode/support/config.hpp.  Generated from config.hpp.in by configure.  */
/* gecode/support/config.hpp.in.  Generated from configure.ac by autoheader.  */

/* Whether to build with default memory allocator */
#define GECODE_ALLOCATOR /**/

/* Whether to include audit code */
/* #undef GECODE_AUDIT */

/* User-defined prefix of dll names */
#define GECODE_DLL_USERPREFIX ""

/* User-defined suffix of dll names */
#define GECODE_DLL_USERSUFFIX ""

/* Supported version of FlatZinc */
#define GECODE_FLATZINC_VERSION "1.6"

/* Whether gcc understands visibility attributes */
#define GECODE_GCC_HAS_CLASS_VISIBILITY /**/

/* whether __builtin_ffsll is available */
/* #undef GECODE_HAS_BUILTIN_FFSLL */

/* whether __builtin_popcountll is available */
/* #undef GECODE_HAS_BUILTIN_POPCOUNTLL */

/* Whether counting-based search support available */
/* #undef GECODE_HAS_CBS */

/* Whether CPProfiler support available */
/* #undef GECODE_HAS_CPPROFILER */

/* Whether to build FLOAT variables */
#define GECODE_HAS_FLOAT_VARS /**/

/* Whether Gist is available */
/* #undef GECODE_HAS_GIST */

/* Whether GNU hash_map is available */
#define GECODE_HAS_GNU_HASH_MAP /**/

/* Whether to build INT variables */
#define GECODE_HAS_INT_VARS /**/

/* Whether MPFR is available */
/* #undef GECODE_HAS_MPFR */

/* Whether we have mtrace for memory leak debugging */
/* #undef GECODE_HAS_MTRACE */

/* Whether Qt is available */
/* #undef GECODE_HAS_QT */

/* Whether to build SET variables */
#define GECODE_HAS_SET_VARS /**/

/* Whether unistd.h is available */
#define GECODE_HAS_UNISTD_H 1

/* Whether unordered_map is available */
#define GECODE_HAS_UNORDERED_MAP /**/

/* Gecode version */
#define GECODE_LIBRARY_VERSION "6-2-0"

/* Heap memory alignment */
/* #undef GECODE_MEMORY_ALIGNMENT */

/* How to check allocation size */
/* #undef GECODE_MSIZE */

/* Whether to track peak heap size */
/* #undef GECODE_PEAKHEAP */

/* Whether we need malloc.h */
/* #undef GECODE_PEAKHEAP_MALLOC_H */

/* Whether we need malloc/malloc.h */
/* #undef GECODE_PEAKHEAP_MALLOC_MALLOC_H */

/* Whether we are compiling static libraries */
/* #undef GECODE_STATIC_LIBS */

/* Whether we have Mac OS threads */
/* #undef GECODE_THREADS_OSX */

/* Whether we have Mac OS threads (new version) */
/* #undef GECODE_THREADS_OSX_UNFAIR */

/* Whether we have posix threads */
#define GECODE_THREADS_PTHREADS 1

/* Whether we have posix spinlocks */
#define GECODE_THREADS_PTHREADS_SPINLOCK 1

/* Whether we have windows threads */
/* #undef GECODE_THREADS_WINDOWS */

/* Use clock() for time-measurement */
/* #undef GECODE_USE_CLOCK */

/* Use gettimeofday for time-measurement */
#define GECODE_USE_GETTIMEOFDAY 1

/* Gecode version */
#define GECODE_VERSION "6.2.0"

/* Gecode version */
#define GECODE_VERSION_NUMBER 600200

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the <gmp.h> header file. */
/* #undef HAVE_GMP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <malloc.h> header file. */
/* #undef HAVE_MALLOC_H */

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the <mpfr.h> header file. */
/* #undef HAVE_MPFR_H */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "users@gecode.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "GECODE"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "GECODE 6.2.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "gecode"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "6.2.0"

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* How to tell the compiler to really, really inline */
#define forceinline inline
