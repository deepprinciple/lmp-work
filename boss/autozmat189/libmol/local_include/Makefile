SHELL = /bin/csh
#INCDIR = ../../../../libmol/local_include
INCDIR = ${HOME}/XChemEdit/libmol/local_include

TARGET = ../include/molecule.h
OBJ = \
	hash.h \
	pdb.h \
	strscan.h \
	list.h \
	orbital.h \
	mol.h \
	elemdata.h \
	error.h \
	zmatrix.h \
	autozmat.h \
	gausszmat.h \
	jaguar.h \
	mopaczmat.h \
	boss.h \
	bossinput.h \
	bosspar.h \
	mmodel.h \
	cartesian.h \
	utils.h \
	hydrogen.h \
	mdl.h \
	pepz.h \
	analgeom.h \
	tripos.h \
	connect.h

OBJ_EXTRA = macros.h

all: $(OBJ)
	/bin/rm -f ${TARGET}
	touch ${TARGET}
	echo "#ifndef _MOLECULE_H_"   >> ${TARGET}
	echo "#define _MOLECULE_H_"   >> ${TARGET}
	echo ""                       >> ${TARGET}
	echo "#ifdef HAVE_PLOT"       >> ${TARGET}
	echo '#include "plot.h"'      >> ${TARGET}
	echo "#endif /* HAVE_PLOT */" >> ${TARGET}
	cat ${OBJ}                    >> ${TARGET}
	echo "#endif /* _MOLECULE_H_ */" >> ${TARGET}

link:
	echo "${OBJ} ${OBJ_EXTRA}" |awk '{\
		for(i=1;i<=NF;i++) {\
			cmd = sprintf("/bin/rm -f %s", $$i);\
			system(cmd);\
			cmd = sprintf("ln -s ${INCDIR}/%s %s", $$i, $$i);\
			system(cmd);\
		}\
	}'

upload:
	@echo 'Execute the following on the command line:'
	@echo 'foreach file (*.h)'
	@echo ' /bin/cp -i $$file ${INCDIR}/$$file'
	@echo 'end'

