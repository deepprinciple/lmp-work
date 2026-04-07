#ifndef __MMODEL_H_
#define __MMODEL_H_

#define MMOD_MAX_BONDS	6

extern int	MModelToAn (int type);
extern int	GetMModelAtomType (AtomPtr atom);
extern int	GetMModelBondType (int bondtype);
extern MolPtr	FLoadFullMModel (FILE *fp, int natoms, char *title);
extern MolPtr	FLoadMModel (FILE *fp);
extern void	FPrintMModel (FILE *fp, MolPtr m);
#endif
