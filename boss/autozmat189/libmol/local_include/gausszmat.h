#ifndef __GAUSSZMAT_H_
#define __GAUSSZMAT_H_
extern void	FreeGaussZmat (GaussZmatPtr *gauss);
extern GaussZmatPtr	CopyGaussZmat (GaussZmatPtr oldgauss);
extern GaussZmatPtr	FLoadGaussZmat (FILE *fp);
extern MolPtr	CreateMolFromGauss (GaussZmatPtr gauss);
extern GaussZmatPtr	NewGaussZmat (void);
extern GaussZmatPtr	CreateGaussZmatFromMol (MolPtr mol);
extern GaussZmatPtr	CreateAutoGaussZmatFromMol (MolPtr mol);

extern char	*GetGaussRoute(void);
extern void	SetGaussRoute(char *);
extern void	SPrintGaussRoute (char **str_sum, GaussZmatPtr gauss);
extern void	SPrintGaussCartZmat (char **str_sum, GaussZmatPtr gauss);
extern void	SPrintGaussZmat (char **str_sum, GaussZmatPtr gauss);
extern void	FPrintGaussCartZmat (FILE *fp, GaussZmatPtr gauss);
extern void	FPrintGaussZmat (FILE *fp, GaussZmatPtr gauss);
#endif
