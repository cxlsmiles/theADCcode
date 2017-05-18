/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the `efence' library (-lefence). */
/* #undef HAVE_LIBEFENCE */

/* Define to 1 if you have the `ma' library (-lma). */
/* #undef HAVE_LIBMA */

/* Define to 1 if you have the `molcas' library (-lmolcas). */
/* #undef HAVE_LIBMOLCAS */

/* Major version number */
#define MAJOR 0

/* Minor version number */
#define MINOR 12

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "victor.vysotskiy@pci.uni-heidelberg.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "PHIS"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "PHIS 0.12.-0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "PHIS"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.12.-0"

/* Patch level */
#define PATCH -0

/* should the GUK backend be included */
#define WITH_GUK 

/* should the MOLCAS backend be included */
/* #undef WITH_MOLCAS */

/* configure back-end stub */
/* #undef WITH_STUB */

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
