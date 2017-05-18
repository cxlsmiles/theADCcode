AC_INIT(PHIS, 0.12.-0, [victor.vysotskiy@pci.uni-heidelberg.de], [PHIS])
AC_CONFIG_HEADER(include/config.h include/phis.h)

dnl ====================== configure flags =============================
AC_ARG_ENABLE( debug,
[
  --enable-debug          enable compiler flags for debugging],
[CPPFLAGS="$CPPFLAGS -DDEBUG";
 AC_CHECK_LIB(efence, main)]
)

AC_ARG_WITH( molcas,
[
  --with-molcas           enable MOLCAS support],
)

AC_ARG_WITH( guk,
[
  --with-guk              enable GAMESS-UK support],
)

AC_ARG_WITH( stub,
[
  --with-stub             enable back-end stub (for debugging)],
)

AC_ARG_WITH( i8,
[
  --with-i8               treat integer as integer(8)],
)

AC_ARG_WITH( large_arrays,
[
  --with-large_arrays     support for large arrays (>2Gb)],
)

dnl ===================	system specific checks =========================
AC_CANONICAL_TARGET
dnl ===================	check for compilers  ===========================	
AC_PROG_F77([pgf90 g95 ifort g95 gfortran f90], Fortran 90)
AC_PROG_CC
dnl ===================	check for paths  ===============================
AC_ARG_VAR(PERL,[Path to the Perl interpreter, should be version 5.5 or later.]) 
AC_PATH_PROG(PERL,$PERL perl)
AC_PATH_PROG(CPP,cpp,cpp)
AC_PATH_PROG(RM,rm,rm)
dnl ===================	Check for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_BIGENDIAN
AC_C_INLINE
AC_F77_WRAPPERS
dnl ===================	target & compiler dependent settings ===========
dnl ===================	DEC options  ===================================
if test "$target_vendor" = "dec"; then
  if test "$F77" = "f90"; then
     FFLAGS="-i8 -fpe0 -fpconstant -ieee"
    if test -n "$enable_debug"; then
     FFLAGS="$FFLAGS -O0 -g -check overflow -check underflow";
    else
     FFLAGS="$FFLAGS -O2 -tune host -arch host";  
    fi 
  fi
  if test "$GCC" = "yes"; then
    CFLAGS="$CFLAGS -mieee";
  fi
fi
FC=`echo $F77 | awk --field-separator=/ '{print $NF}'`
CCO=`echo $CC | awk --field-separator=/ '{print $NF}'`
dnl =================== X86_64 options =================================
if test  $target_cpu = "x86_64"; then
   if test "$FC" = "pgf90" -o "$FC" = "pgf77" ; then
dnl =================== settings for PGI compiler ======================
      FFLAGS="-tp k8-64";
      if test -n "$enable_debug"; then 
         FFLAGS="$FFLAGS -g -O0 -Mchkfpstk -Mchkstk  -Ktrap=align,denorm,divz,fp,inexact,inv,ovf,unf -Minform=inform -Kieee";
      else
         FFLAGS="$FFLAGS -fastsse -Mipa=fast "   
      fi
      if test -n "$with_i8"; then
         FFLAGS="$FFLAGS -i8";
      fi 
      if test -n "$with_large_arrays"; then
         FFLAGS="$FFLAGS -mcmodel=medium -Mquad -Msmartalloc -Msmart  -Mnomovnt -Mnontemporal";
      fi
   fi
   if test "$FC" = "ifort"  ; then
dnl =================== settings for PGI compiler ======================
      FFLAGS="-axT -xT";
      if test -n "$enable_debug"; then
         FFLAGS="$FFLAGS -g -O0 -mp -no-prec-div -debug=all";
      else
         FFLAGS="$FFLAGS -O3 -ip -pad "
      fi
      if test -n "$with_i8"; then
         FFLAGS="$FFLAGS -i8";
      fi
      if test -n "$with_large_arrays"; then
         FFLAGS="$FFLAGS -mcmodel=medium";
      fi
   fi

   if test "$CCO" = "icc" ; then 
      if test -n "$enable_debug"; then
         CFLAGS="-O0 -g -mp -no-prec-div -debug=all";
      else  
         CFLAGS="-axT -xT -O3 -ip";
      fi
      if test -n "$with_large_arrays"; then
         CFLAGS="$CFLAGS -mcmodel=medium";
      fi
   fi

   if test "$CCO" = "pgcc" ; then 
      if test -n "$enable_debug"; then
         CFLAGS="-O0 -g -Mchkfpstk -Mchkstk  -Kieee";
      else  
         CFLAGS="-fastsse -Mipa=fast";
      fi
      if test -n "$with_large_arrays"; then
         CFLAGS="$CFLAGS -mcmodel=medium -Mquad -Msmartalloc -Msmart  -Mnomovnt ";
      fi
   fi
dnl ===================  settings for GFortran compiler =================== 
if test "$F77" = "gfortran" ; then
dnl settings for GNU_G95 compiler
   FFLAGS="-fno-second-underscore";

  if test -n "$enable_debug"; then 
     FFLAGS="$FFLAGS -g -O0 -Wall -pedantic -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -mieee-fp -ftrace=full";
  else
     FFLAGS="$FFLAGS -O2 -funroll-loops -fomit-frame-pointer  -mtune=opteron -mfpmath=sse,387 -m64 -march=opteron"   
  fi

  if test -n "$with_i8"; then
     FFLAGS="$FFLAGS -fdefault-integer-8";
  fi 
  
  if test -n "$with-large_arrays"; then
     FFLAGS="$FFLAGS -mcmodel=medium";
  fi
fi

dnl settings for G95 compiler
if test "$FC" = "g95" ; then
dnl settings for GNU_G95 compiler
     FFLAGS="-fno-second-underscore -mieee-fp -ftrace=full";

  if test -n "$enable_debug"; then 
     FFLAGS="$FFLAGS -g -O0 -Wall -pedantic -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal";
  else
     FFLAGS="$FFLAGS -O2 -funroll-loops -fomit-frame-pointer  -mtune=opteron -mfpmath=sse,387 -m64 -march=opteron "
  fi

  if test -n "$with_i8"; then
     FFLAGS="$FFLAGS -i8";
  fi 
  
  if test -n "$with-large_arrays"; then
     FFLAGS="$FFLAGS -mcmodel=medium";
  fi

fi
dnl end of X86_64 options
fi
dnl ================= end system specific checks =======================

dnl ================= create MOLCAS's include files ====================
    PPFLAGS=`fgrep -iw CPPFLAGS $MOLCAS/Symbols | awk -F= '{print $2}' | tr "'" " "`
    PPFLAGS=${PPFLAGS//-I\$\{INCDIR\}/""};
    PPFLAGS=${PPFLAGS//-I\$\{GAINC\}/""};
    PPFLAGS=${PPFLAGS//-I\./""};
    cpp $PPFLAGS $MOLCAS/src/Include/itmax.fh   | egrep -vi '^(#|\*|Celse|Cend|Cstart|c;)' >  include/info.inc ;
    cpp $PPFLAGS $MOLCAS/src/Include/info.fh    | egrep -vi '^(#|\*|Celse|Cend|Cstart|c;)' >> include/info.inc ;
    cpp $PPFLAGS $MOLCAS/src/Include/WrkSpc.inc | egrep -vi '^(#|\*|Celse|Cend|Cstart|c;)' >> include/info.inc ;
dnl ============== back-end specific configurations ====================

dnl --- MOLCAS ---
if test -n "$with_molcas"; then
  AC_DEFINE(WITH_MOLCAS, ,[should the MOLCAS backend be included])

  if test -z "$MOLCAS"; then
    AC_MSG_ERROR('environment variable MOLCAS not set')
  else  
       VERSION=`grep "Data Molcas_version" $MOLCAS/src/Include/Molcas_version.inc  | awk '{print $5}' | tr "}'/" " "`
       VERSION=${VERSION%% *}
    if test -z "$VERSION"; then
       AC_MSG_ERROR('MOLCAS VERSION is not detected');
    elif test "$VERSION" == "6.2" ; then
       Vers="6_2";
       Status="FALSE";
    elif test "$VERSION" == "6.4" ; then
       Vers="6_4";
       Status="TRUE";
    else 
       AC_MSG_ERROR('The version $VERSION is not supported yet');
    fi    
  fi

dnl FIXME: Maybe there is a machine independent way?
  AC_SUBST(RtoI)
  AC_SUBST(ItoB)
  AC_SUBST(Vers)
  AC_SUBST(Vendor)
  AC_SUBST(Status)
  AC_SUBST(DYN_FLAGS)
  if test $target_vendor = "dec"; then
    RtoI=1; ItoB=8;
  else
    RtoI=2; ItoB=4;
  fi
dnl VVP _D_I8_ 
  if test -n "$with_i8"; then
    RtoI=1; ItoB=8;
  fi 

  if test "$FC" = "pgf90" -o "$FC" = "pgf77" ; then
     Vendor="pgi";
     DYN_FLAGS="-shared -Bstatic_pgi"
  elif test "$FC" = "ifort"; then
     Vendor="intel";
     DYN_FLAGS="-i-static -shared -static-intel"
  elif test "$FC" = "gfortran" -o "$FC" = "g95"; then
     Vendor="gnu";
     DYN_FLAGS="-shared  -shared-libgcc  -symbolic"
  fi

dnl This is pretty non-standard, maybe I should use autoconf to
dnl include the relevant files??
  FFLAGS="$FFLAGS"

dnl FIXME: Not sure if this is the way to do it...
  MA=`ls -d $MOLCAS/g/lib/*` 
  LDFLAGS=" "
  tmp_path=`find $MOLCAS/g/lib/ -name 'libma*'`;
  tmp_path=`dirname $tmp_path` 
dnl  if test -n "$tmp_path"; then
dnl    LDFLAGS="$LDFLAGS -L$tmp_path";
dnl    AC_CHECK_LIB(ma, MA_alloc_get,
dnl    , dnl default here
dnl      AC_MSG_ERROR(cannot link MOLCAS memory manager) 
dnl    )
dnl  fi
  AC_CHECK_LIB(molcas, main)
  subdirs="$subdirs molcas"
  ac_output="$ac_output molcas/Makefile  include/molcas.fh Makefile"
fi;

dnl --- GUK ---
if test -n "$with_guk"; then
  AC_DEFINE(WITH_GUK, ,[should the GUK backend be included])
  if test "$FC" = "pgf90" -o "$FC" = "pgf77" ; then
     Vendor="pgi";
     DYN_FLAGS="-shared -Bstatic_pgi"
  elif test "$FC" = "ifort"; then
     Vendor="intel";
     DYN_FLAGS="-shared"
  elif test "$FC" = "gfortran" -o "$FC" = "g95"; then
     Vendor="gnu";
     DYN_FLAGS="-shared  -shared-libgcc  -symbolic"
  fi
  subdirs="$subdirs guk"
  ac_output="$ac_output guk/Makefile"
fi;

dnl --- STUB ---
if test -n "$with_stub"; then
  AC_DEFINE(WITH_STUB, ,[configure back-end stub])
  subdirs="$subdirs stub"
  ac_output="$ac_output stub/Makefile"
fi;

dnl ============ end back-end specific configurations ==================

dnl ===========================	output =================================

dnl --- make the version info available in the config header ---
AC_PROG_AWK
MAJOR=`echo $PACKAGE_VERSION | $AWK -F'.' '{ print $1 }'`
MINOR=`echo $PACKAGE_VERSION | $AWK -F'.' '{ print $2 }'`
PATCH=`echo $PACKAGE_VERSION | $AWK -F'.' '{ print $3 }'`
AC_DEFINE_UNQUOTED(MAJOR, $MAJOR, [Major version number])
AC_DEFINE_UNQUOTED(MINOR, $MINOR, [Minor version number])
AC_DEFINE_UNQUOTED(PATCH, $PATCH, [Patch level])

LIBS="../lib/libphis.a $LIBS $MOLCAS/lib/libmolcas.a $MA/libma.a"

subdirs="$subdirs phis examples"

AC_OUTPUT($ac_output phis/Makefile examples/Makefile Makefile)

dnl ====================================================================
