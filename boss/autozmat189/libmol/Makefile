CC = cc
CFLAGS = -g

INST_LIB_DIR = ../../lib
INST_INC_DIR = ../../include

INC_DIR = ../include
LOCAL_INC_DIR = ../local_include

# authentic libmol directory (leave it unknown)
LIB_MOL_DIR =

all:
	(cd src; make CC="$(CC)" CFLAGS="$(CFLAGS)" \
		INST_LIB_DIR="${INST_LIB_DIR}" \
		INST_INC_DIR="${INST_INC_DIR}" \
		INC_DIR="${INC_DIR}" \
		LOCAL_INC_DIR="${LOCAL_INC_DIR}";)

install:
	(cd src; make CC="$(CC)" CFLAGS="$(CFLAGS)" \
		INST_LIB_DIR="${INST_LIB_DIR}" \
		INST_INC_DIR="${INST_INC_DIR}" \
		INC_DIR="${INC_DIR}" \
		LOCAL_INC_DIR="${LOCAL_INC_DIR}" install;)

inc:
	(cd local_include; make)

clean:
	(cd src; make clean)

clean_installed:
	(cd src; make clean_installed)

download:
	csh Download.csh ${LIB_MOL_DIR}

