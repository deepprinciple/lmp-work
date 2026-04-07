#ifndef __TRIPOS_H_
#define __TRIPOS_H_

extern MolPtr	FLoadTriposMol (FILE *fp);
extern MolPtr	FLoadTriposMol2 (FILE *fp);
extern void	FPrintTriposMol (FILE *fp, MolPtr m);
extern void	FPrintTriposMol2 (FILE *fp, MolPtr m);

#endif
