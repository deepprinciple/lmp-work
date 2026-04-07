CC = cc
CXX = ${CC}

INST_INC_DIR = ../include
INST_LIB_DIR = ../lib
INST_BIN_DIR = ../bin

CFLAGS = -g
COMFLAGS = ${CFLAGS} -I${INST_INC_DIR} -DDEBUG
VERSION = 1.87
LDFLAGS = -L${INST_LIB_DIR}
LIBS = -lmol -lm
MV = /bin/mv

OBJ = autozmat.o get_opt.o error.o oplstype.o
TARGET = autozmat

$(TARGET):	$(OBJ)
	$(CC) $(COMFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LIBS)
	${MV} $(TARGET) ${INST_BIN_DIR}

.c.o :
	$(CC) $(COMFLAGS) -DVERSION=\"$(VERSION)\" -c -o $@ $<

.c :
	$(CC) $(COMFLAGS) -DVERSION=\"$(VERSION)\" $< -o $@ $(LDFLAGS) $(LIBS)

remove_target:
	/bin/rm -f $(TARGET)

clean:
	-@ rm -f *.o $(TARGET) $(CLEAN) core *~

clean_installed: clean
	-@ rm -f $(INST_BIN_DIR)/$(TARGET)
	-@ rm -f $(INST_BIN_DIR)/$(TARGET)*

#########################################################
autozmat.o: autozmat.h oplstype.h
oplstype.o: autozmat.h oplstype.h

