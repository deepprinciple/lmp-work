AUTOZMAT = ../bin/autozmat
AUTHENTIC_DIR = Authentic
SCRATCH_DIR = Scratch
ATOMTYPE_OPT =

#----------- DO NOT EDIT BELOW --------------
SHELL = /bin/csh
OPTIONS = AUTOZMAT=${AUTOZMAT} AUTHENTIC_DIR=${AUTHENTIC_DIR} \
	SCRATCH_DIR=${SCRATCH_DIR}

test:
	sh TEST.sh ${OPTIONS} INSTALL=FALSE ATOMTYPE_OPT=${ATOMTYPE_OPT}

install: clean
	sh TEST.sh ${OPTIONS} INSTALL=TRUE ATOMTYPE_OPT=
	sh TEST.sh ${OPTIONS} INSTALL=TRUE ATOMTYPE_OPT=-a

clean:
	-/bin/rm -f core ${SCRATCH_DIR}/*

clean_installed: clean
	-/bin/rm -f ${AUTHENTIC_DIR}/atomtyping/* ${AUTHENTIC_DIR}/noatomtyping/*

