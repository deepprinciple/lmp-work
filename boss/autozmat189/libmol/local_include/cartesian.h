#ifndef __CARTESIAN_H_
#define __CARTESIAN_H_

#define MINDTOOL_BREAK_PERCENT	0x1

extern MolPtr	FLoadZmatCartesian (FILE *fp, ZmatPtr *zmatlist_ret);
extern MolPtr	FLoadCartesian (FILE *fp);
extern MolPtr	FLoadMindtool (FILE *fp, int flags);
extern void	FPrintCartesian (FILE *fp, ChainPtr chainlist);
extern void	FPrintMindtool (FILE *fp, ChainPtr chainlist);

#endif
