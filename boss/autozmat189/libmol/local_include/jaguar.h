#ifndef __JAGUAR_H_
#define __JAGUAR_H_

#define L_JAGUAR_TYPE	data1
#define L_JAGUAR_LINE	str1

/* JaguarPtr->flags */
#define JAGUAR_CARTESIAN	0x1

typedef struct {
	unsigned long	flags;
	ListPtr	header;
	ZmatPtr	zmatlist;
	ListPtr	gen;
	ListPtr	gvb;
	ListPtr	lmp2;
	ListPtr	atomic;
	ListPtr	vdw;
	ListPtr	hess;
	ListPtr	guess;
	ListPtr	pointch;
	ListPtr	efields;
	ListPtr	ham;
	ListPtr	orbman;
	ListPtr	echo;
	ListPtr	path;
	} JaguarZmat, *JaguarZmatPtr;

extern void	FreeJaguarZmat (JaguarZmatPtr *jaguar);
extern JaguarZmatPtr	FLoadJaguarZmat (FILE *fp);
extern void	SPrintJaguarZmat (char **str_sum, JaguarZmatPtr j);
extern void	FPrintJaguarZmat (FILE *fp, JaguarZmatPtr j);
extern GaussZmatPtr	CreateGaussZmatFromJaguar (JaguarZmatPtr j);
extern JaguarZmatPtr	CreateJaguarZmatFromGauss (GaussZmatPtr g);
extern JaguarZmatPtr	CreateJaguarZmatFromMol (MolPtr mol);
extern JaguarZmatPtr	CopyJaguarZmat (JaguarZmatPtr oldjaguar);

#endif
