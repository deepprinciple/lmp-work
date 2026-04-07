#ifndef __MDL_H_
#define __MDL_H_

extern MolPtr	FLoadMdlMol (FILE *);
extern MolPtr	FLoadMdlSd (FILE *);
extern void	FPrintMdlMol (FILE *, MolPtr);
extern void	FPrintMdlSd (FILE *, MolPtr);

#endif
