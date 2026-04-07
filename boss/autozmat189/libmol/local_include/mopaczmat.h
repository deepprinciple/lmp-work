#ifndef __MOPACZMAT_H_
#define __MOPACZMAT_H_
extern MopacZmatPtr	FLoadMopacZmat (FILE *fp);
extern char	*GetMopacRoute(void);
extern void	SetMopacRoute(char *);
extern void	SPrintMopacZmat(char **, MopacZmatPtr);
extern void	FPrintMopacZmat(FILE *, MopacZmatPtr);
#endif
