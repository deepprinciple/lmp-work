#-------------------------------------------------------#
### SGI and generic UNIX ###
#TAR_ARCHIVE = tar Lcvf

### LINUX ###
TAR_ARCHIVE = tar hcvf
#-------------------------------------------------------#

SHELL = /bin/sh
CP = /bin/cp
RM = /bin/rm -f

TARGET = autozmat
VERSION = 1.86
ARCHIVE = ${TARGET}${VERSION}

### for installation
UPPER_INST_BIN_DIR = ../../../bin
UPPER_INST_INC_DIR = ../../../include
UPPER_INST_LIB_DIR = ../../../lib

LOCAL_LIB_DIR = ../lib/opls
DEST_LIB_DIR = ${UPPER_INST_LIB_DIR}/opls
DEST_LIB_FILES = \
	$(DEST_LIB_DIR)/FG_hash.tbl \
	$(DEST_LIB_DIR)/FG_simil.tbl \
	$(DEST_LIB_DIR)/element.lis \
	$(DEST_LIB_DIR)/oplsaa.type \
	$(DEST_LIB_DIR)/oplsaa_simil.tbl \
	$(DEST_LIB_DIR)/paramstd.dat \
	$(DEST_LIB_DIR)/rings.lis

LOCAL_BIN_DIR = ../bin
DEST_BIN_DIR = ${UPPER_INST_BIN_DIR}
DEST_BIN_FILES = ${DEST_BIN_DIR}/${TARGET}

### --- distribution directories and files
DIST_DIRS = lib dist
DIST_LINKS = 00README Copyright INSTALL bin dist/Makefile docs \
	include lib/opls lib/boss libmol sample src utils

###########################################################################
#####     DO NOT EDIT BELOW!
###########################################################################

### default values -------------
DIRS = ${DIST_DIRS}
LINKS = ${DIST_LINKS}
###-----------------------------

dist: dist_atomtype
dist_atomtype:
	make DIRS="${DIST_DIRS}" LINKS="${DIST_LINKS} libatomtype" make_tar_dir \
	ARCHIVE="${ARCHIVE}" TAR_ARCHIVE="${TAR_ARCHIVE}"

dist_noatomtype:
	make DIRS="${DIST_DIRS}" LINKS="${DIST_LINKS}" \
	ATOMTYPE="no" \
	make_tar_dir \
	ARCHIVE="${ARCHIVE}" TAR_ARCHIVE="${TAR_ARCHIVE}"

make_tar_dir:
	@ if [ -r ${ARCHIVE} ]; then true; else mkdir ${ARCHIVE}; fi

	@ for file in ${DIRS}; do \
		if [ -r ${ARCHIVE}/$$file ]; then ${RM} -r ${ARCHIVE}/$$file ; fi ;\
		mkdir ${ARCHIVE}/$$file ;\
	  done
	@ for file in ${LINKS}; do \
		if [ -r ${ARCHIVE}/$$file ]; then ${RM} -r ${ARCHIVE}/$$file ; fi ;\
		if [ -L ${ARCHIVE}/$$file ]; then ${RM} -r ${ARCHIVE}/$$file ; fi ;\
		if [ `dirname $$file` = "." ]; then ln -s ../../$$file ${ARCHIVE}/$$file ; \
		else ln -s ../../../$$file ${ARCHIVE}/$$file ; fi ;\
	  done

	# prepare Makefile
	${RM} ${ARCHIVE}/Makefile;
	@ if [ "${ATOMTYPE}" = "no" ]; then \
		sed '/^ATOMTYPE_FLAGS.*DHAVE_ATOMTYPE/d' ../Makefile |\
		sed '/^MY_LIBS.*latomtype/d' > ${ARCHIVE}/Makefile;\
	else \
		cat ../Makefile > ${ARCHIVE}/Makefile;\
	fi

	@ if [ -r ${ARCHIVE}.tar ]; then ${RM} ${ARCHIVE}.tar ; fi
	@ if [ -r ${ARCHIVE}.tar.gz ]; then ${RM} ${ARCHIVE}.tar.gz ; fi
	@ ${TAR_ARCHIVE} ${ARCHIVE}.tar ${ARCHIVE}
	@ gzip ${ARCHIVE}.tar

clean:
	-${RM} -r ${TARGET}*

clean_installed_upper:
	-${RM} ${DEST_BIN_FILES} ${DEST_LIB_FILES}

install_upper: install_upper_bin install_upper_lib
install_upper_bin: ${DEST_BIN_FILES}
install_upper_lib: ${DEST_LIB_FILES}

${DEST_BIN_DIR}/${TARGET}: ${LOCAL_BIN_DIR}/${TARGET}
	${CP} ${LOCAL_BIN_DIR}/${TARGET} $@
	-strip $@

$(DEST_LIB_DIR)/FG_hash.tbl: $(LOCAL_LIB_DIR)/FG_hash.tbl
	$(CP) $(LOCAL_LIB_DIR)/FG_hash.tbl $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/FG_simil.tbl: $(LOCAL_LIB_DIR)/FG_simil.tbl
	$(CP) $(LOCAL_LIB_DIR)/FG_simil.tbl $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/element.lis: $(LOCAL_LIB_DIR)/element.lis
	$(CP) $(LOCAL_LIB_DIR)/element.lis $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/oplsaa.type: $(LOCAL_LIB_DIR)/oplsaa.type
	$(CP) $(LOCAL_LIB_DIR)/oplsaa.type $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/oplsaa_simil.tbl: $(LOCAL_LIB_DIR)/oplsaa_simil.tbl
	$(CP) $(LOCAL_LIB_DIR)/oplsaa_simil.tbl $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/paramstd.dat: $(LOCAL_LIB_DIR)/paramstd.dat
	$(CP) $(LOCAL_LIB_DIR)/paramstd.dat $(DEST_LIB_DIR)
$(DEST_LIB_DIR)/rings.lis: $(LOCAL_LIB_DIR)/rings.lis
	$(CP) $(LOCAL_LIB_DIR)/rings.lis $(DEST_LIB_DIR)

