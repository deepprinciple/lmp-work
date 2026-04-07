#CC = gcc
CC = cc

INST_LIB_DIR = ../../lib
INST_INC_DIR = ../../include
INC_DEST = ${INST_INC_DIR}/mol

INC_DIR = ../include
LOCAL_INC_DIR = ../local_include

CFLAGS = -g
COMFLAGS = ${CFLAGS} -I${INC_DIR} -I${LOCAL_INC_DIR} -I/usr/include

LDFLAGS = 
AR = ar
ARFLAGS = cr
LIBMOL = libmol.a

CP = /bin/cp

### following three lines are need for 'link' only.
SHELL = /bin/csh
SRCDIR = ${HOME}/XChemEdit/libmol/src
LINK_FILE_SUFFIX = .c

OBJ = \
	analgeom.o \
	autozmat.o \
	bossinput.o \
	bosspar.o \
	cartesian.o \
	elemdata.o \
	error.o \
	filepdb.o \
	gausszmat.o \
	hash.o \
	hydrogen.o \
	jaguarzmat.o \
	list.o \
	mdl.o \
	mmodel.o \
	mol.o \
	mopaczmat.o \
	pdb.o \
	pepz.o \
	strscan.o \
	tripos.o \
	utils.o \
	zmatrix.o \
	connect.o

LOCAL_INC_OBJ = \
	$(LOCAL_INC_DIR)/hash.h \
	$(LOCAL_INC_DIR)/pdb.h \
	$(LOCAL_INC_DIR)/strscan.h \
	$(LOCAL_INC_DIR)/list.h \
	$(LOCAL_INC_DIR)/orbital.h \
	$(LOCAL_INC_DIR)/mol.h \
	$(LOCAL_INC_DIR)/elemdata.h \
	$(LOCAL_INC_DIR)/error.h \
	$(LOCAL_INC_DIR)/zmatrix.h \
	$(LOCAL_INC_DIR)/autozmat.h \
	$(LOCAL_INC_DIR)/gausszmat.h \
	$(LOCAL_INC_DIR)/jaguar.h \
	$(LOCAL_INC_DIR)/mopaczmat.h \
	$(LOCAL_INC_DIR)/boss.h \
	$(LOCAL_INC_DIR)/bossinput.h \
	$(LOCAL_INC_DIR)/bosspar.h \
	$(LOCAL_INC_DIR)/mmodel.h \
	$(LOCAL_INC_DIR)/cartesian.h \
	$(LOCAL_INC_DIR)/utils.h \
	$(LOCAL_INC_DIR)/hydrogen.h \
	$(LOCAL_INC_DIR)/mdl.h \
	$(LOCAL_INC_DIR)/pepz.h \
	$(LOCAL_INC_DIR)/analgeom.h \
	$(LOCAL_INC_DIR)/tripos.h \
	$(LOCAL_INC_DIR)/connect.h

INST_INC = ${INC_DEST}/molecule.h

all: $(LIBMOL)

$(LIBMOL): ${INC_DIR}/molecule.h ${OBJ}
	${AR} ${ARFLAGS} ${LIBMOL} ${OBJ}

.c.o :
	$(CC) $(COMFLAGS) -c -o $@ $<

clean:
	-rm -f core a.out ${LIBMOL}
	-rm -f *.o

clean_installed: clean
	-rm -f ${INC_DEST}/*
	-rm -f ${INST_LIB_DIR}/${LIBMOL}

##### installation #####
install: ${INST_LIB_DIR}/${LIBMOL} ${INST_INC}

${INST_LIB_DIR}/${LIBMOL}: ${LIBMOL}
	${CP} ${LIBMOL} ${INST_LIB_DIR}
${INC_DEST}/molecule.h: ${INC_DIR}/molecule.h
	${CP} ${INC_DIR}/molecule.h ${INC_DEST}

${INC_DIR}/molecule.h: $(LOCAL_INC_OBJ)
	(cd ${LOCAL_INC_DIR}; make;)

#------------------------------------------------------
# local use only
#------------------------------------------------------
link:
	echo "${OBJ}" |awk '\
	BEGIN {\
		suffix = ".c";\
		if ($LINK_FILE_SUFFIX != "") suffix = $LINK_FILE_SUFFIX;\
	}\
	{\
		for(i=1;i<=NF;i++) {\
			base = $$i;\
			sub(/\.o/, "", base);\
			orig_file = sprintf("${SRCDIR}/%s.C", base);\
			link_file = sprintf("%s%s", base, suffix);\
			cmd = sprintf("/bin/rm -f %s", link_file);\
			system(cmd);\
			cmd = sprintf("ln -s %s %s", orig_file, link_file);\
			system(cmd);\
		}\
	}'

upload:
	@echo 'Execute the following on the command line (in csh):'
	@echo 'foreach file (*.c)'
	@echo '	/bin/cp -i $$file ${SRCDIR}/$$file'
	@echo 'end'

