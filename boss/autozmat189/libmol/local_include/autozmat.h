#ifndef __AUTOZMAT_H_
#define __AUTOZMAT_H_

extern int	AssignZmatProp (ZmatPtr zmatlist);
extern ZmatPtr	NewZmatFromAtom (AtomPtr atom, int z1, int z2, int z3);

/* masks for flags */
#define Z_ADD_BGN_DUMMIES	0x1
#define Z_PUT_SERNO_INFO	0x2

extern ZmatPtr	CreateAutoZmatFromMol (MolPtr mol, int flags);

#endif
