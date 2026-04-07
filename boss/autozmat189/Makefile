########################################################
###  User-Configurable Options (1-6)
###  (This Makefile is configured for UNIX system.
###   On DOS/Windows, try to use DJGPP or CYGWIN)
########################################################

###-----------------------------------------------------
### 1. C compiler and linker
CC = cc
SYS_CFLAGS = -g
AR = ar
ARFLAGS = cr
###-----------------------------------------------------

###-----------------------------------------------------
### 2. Whether you want to use 'signal' UNIX system call
###    Comment out if you are on non-UNIX system or
###    are not using DJGPP/CYGWIN.
SIGNAL_FLAGS = -DHAVE_SIGNAL
###-----------------------------------------------------

###-----------------------------------------------------
### 3. Atomtype library

### If you don't have atomtype library, use below.
ATOMTYPE_FLAGS =
MY_LIBS = -lmol

### If you have Ruhong Zhou's atomtype library, use below.
###-----------------------------------------------------

###-----------------------------------------------------
### 4. 'tar' options
###    (Used for making an archive for distribution.
###     Normal users can skip this part.)
#SGI# TAR_ARCHIVE = tar Lcvf
#UNIX# TAR_ARCHIVE = tar Lcvf
### On Linux, use below
TAR_ARCHIVE = tar hcvf
###-----------------------------------------------------

###-----------------------------------------------------
### 5. Shell, awk, etc.
SHELL = /bin/csh
AWK = awk
MV = /bin/mv
RM = /bin/rm
###-----------------------------------------------------

###-----------------------------------------------------
### 6. PWD
### If PWD environment variable is screwy on your system,
### set PWD to current directory manually.
# PWD = ${HOME}/autozmat
###-----------------------------------------------------

########################################################
### Done! Do not edit below.
########################################################
TARGET = autozmat
VERSION = 1.89

SYS_LIBS = -lm
ROOT_DIR = ${PWD}
LOCAL_INST_LIB_DIR = ${ROOT_DIR}/lib
LOCAL_INST_INC_DIR = ${ROOT_DIR}/include
LOCAL_INST_BIN_DIR = ${ROOT_DIR}/bin

UPPER_INST_BIN_DIR = ${PWD}/../../bin
UPPER_INST_INC_DIR = ${PWD}/../../include
UPPER_INST_LIB_DIR = ${PWD}/../../lib

LIBS = ${MY_LIBS} ${SYS_LIBS}
CFLAGS = ${SYS_CFLAGS} ${SIGNAL_FLAGS} ${ATOMTYPE_FLAGS}
TARGET_VERSION = ${TARGET}${VERSION}

all:
	(cd libmol; \
		make CC="$(CC)" CFLAGS="$(CFLAGS)" \
		INST_LIB_DIR="${LOCAL_INST_LIB_DIR}" \
		INST_INC_DIR="${LOCAL_INST_INC_DIR}" \
		INST_BIN_DIR="${LOCAL_INST_BIN_DIR}" \
		install)

	-(cd libatomtype; \
		if ("${ATOMTYPE_FLAGS}" != "") \
			make CC="$(CC)" \
			CFLAGS="$(CFLAGS)" \
			INST_LIB_DIR="${LOCAL_INST_LIB_DIR}" \
			INST_INC_DIR="${LOCAL_INST_INC_DIR}" \
			INST_BIN_DIR="${LOCAL_INST_BIN_DIR}" \
			AR="$(AR)" ARFLAGS="$(ARFLAGS)" install \
	)

	(cd src; make remove_target; \
		make CC="$(CC)" \
		CFLAGS="$(CFLAGS)" \
		INST_LIB_DIR="${LOCAL_INST_LIB_DIR}" \
		INST_INC_DIR="${LOCAL_INST_INC_DIR}" \
		INST_BIN_DIR="${LOCAL_INST_BIN_DIR}" \
		VERSION="$(VERSION)" \
		LIBS="$(LIBS)")

noatomtype:
	make ATOMTYPE_FLAGS="" \
	CC="$(CC)" CFLAGS="$(SYS_CFLAGS) ${SIGNAL_FLAGS}" \
	LIBS="-lmol ${SYS_LIBS}" \
	all

dist: dist_atomtype
dist_atomtype:
	(make clean && \
		cd dist && \
		make ARCHIVE="${TARGET_VERSION}-atomtype" \
		TAR_ARCHIVE="${TAR_ARCHIVE}"\
		$@ \
	)

dist_noatomtype:
	(make clean && \
		cd dist && \
		make ARCHIVE="${TARGET_VERSION}-noatomtype" \
		TAR_ARCHIVE="${TAR_ARCHIVE}"\
		$@ \
	)

clean:
	(cd libmol; make inc; make clean;)
	(cd src; make clean;)
	if ("${ATOMTYPE_FLAGS}" != "") make clean_atomtype;
	(cd sample; make clean;)

clean_installed:
	(cd libmol; make inc; make clean_installed;)
	(cd src; make clean_installed;)
	if ("${ATOMTYPE_FLAGS}" != "") make clean_atomtype_installed;

clean_atomtype:
	-(cd libatomtype && make clean;)

clean_atomtype_installed:
	-(cd libatomtype && make clean_installed;)

clean_dist:
	(cd dist && make clean)

clean_installed_upper: clean_installed clean_dist
	(cd dist && make $@)

install_upper: download all
	(cd dist; \
		if ("${HAVE_ATOMTYPE}" != "") \
			make \
			UPPER_INST_BIN_DIR="${UPPER_INST_BIN_DIR}" \
			UPPER_INST_INC_DIR="${UPPER_INST_INC_DIR}" \
			UPPER_INST_LIB_DIR="${UPPER_INST_LIB_DIR}" \
			$@; \
		\
		if ("${HAVE_ATOMTYPE}" == "") \
			make \
			UPPER_INST_BIN_DIR="${UPPER_INST_BIN_DIR}" \
			UPPER_INST_INC_DIR="${UPPER_INST_INC_DIR}" \
			UPPER_INST_LIB_DIR="${UPPER_INST_LIB_DIR}" \
			install_upper_bin; \
	)

install:
	@ csh INSTALL

test:
	-(cd sample && make ATOMTYPE_OPT=${ATOMTYPE_OPT})

down: download
download:
	(cd libmol && make LIB_MOL_DIR="${LIB_MOL_DIR}" download)

################################################################
### make Makefile for a given operating system
################################################################
makefile:
	sed 's/^#${OS_VERSION}# //g' Makefile.tpl | sed 's/^#${HAVE_ATOMTYPE}# //g' > Makefile.tmp
	if ("`diff Makefile Makefile.tmp |& head -1`" != "") ${MV} Makefile.tmp Makefile
	-@ ${RM} -f Makefile.tmp

