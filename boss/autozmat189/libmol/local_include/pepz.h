#ifndef __PEPZ_H_
#define __PEPZ_H_

/* flags for FPrintPepz() */
#define PEPZ_ATOMTYPE_FROM_ZMAT		0x1

extern void	FPrintPepz (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist, BossZmatPtr boss, ListPtr atomtypelist, int flags);

#endif
