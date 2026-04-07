#ifndef __ORBITAL_H_
#define __ORBITAL_H_

typedef struct _Cont {
	float	*point;
	float	*prop;
	int	npoints;
	int	color;
	struct _Cont	*next;
	} Cont, *ContPtr;

/* attribute types for SetOrbitalAttributes() and SetDefaultOrbitalAttributes() */
enum	{ORB_ATTR_CUTOFF, ORB_ATTR_SCALEFAC};
#define DEF_ORB_CUTOFF		0.075
#define DEF_ORB_SCALEFAC	1.0
#define DEF_DENSITY_CUTOFF	0.1

/* masks for Orbital->flags */
#define ORB_DENSITY		0x1
#define ORB_POTENTIAL		0x2
#define ORB_NO_NEGATIVE_CONT	0x4
#define ORB_AUTO_CONT		0x8
#define ORB_GAUSS_CUBE		0x10
#define ORB_FORMATTED_CUBE	0x20
#define ORB_SIGN_BY_LINEWIDTH	0x40
#define ORB_SIGN_BY_LINECOLOR	0x80

#define ORB_SIGN_FLAGS	(ORB_SIGN_BY_LINEWIDTH|ORB_SIGN_BY_LINECOLOR)

typedef struct {
	unsigned long	flags;
	int	refno;		/* reference number (filename) */
	int	property_refno;	/* property grid (e.g., potential) */
	int	curorb;
	int	norbs;
	int	nbasis;
	int	ngrids;
	int	mo_first, mo_last;
	int	charge;
	int	multiplicity;
	char	basis[32];
	char	title[84];

	float	scalefactor;
	float	cutoff;
	float	zaxis[3];	/* +Z-axis, used for keeping track of rotation transformation
					in monatomic case. */

	float	prop_min, prop_max;	/* in case, contour has property values */

	FILE	*fp;

	float	*coeff;		/* MO coefficients */
	float	*zeta;		/* zeta values for SEMI empirical MO's */

	Cont	*cont;		/* orbital contours */

	float	*coord;		/* atom cooridnates */
	int	*an;		/* atomic numbers */
	int	natoms;		/* number of atoms */
	} Orbital, *OrbitalPtr;

extern void	FreeOrbital(Orbital **);
extern Orbital	*ReadPsi1(char *);
extern Orbital	*FReadPsi1(FILE *);
extern Cont	*HideContourLines(FILE *, float (*)[4], float *, float *, float, int);
extern void	PrintContour(Cont *);
extern Orbital	*FLoadOrbital (FILE *fp);
extern Orbital	*CopyOrbital(Orbital *);
extern void	FreeCont (Cont **cont);
extern ContPtr	CopyCont (ContPtr oldcont);
extern void	RotateCont (Cont *cont, float (*rotmat)[4], float *center, float *disp);
extern int	ExecutePSI (OrbitalPtr);
extern int	ExecutePSICON (OrbitalPtr);
extern void	FPrintPsi1 (FILE *fp, Orbital *orb);
extern Orbital	*FLoadCubeDensity (FILE *fp);
extern void	SetOrbitalAttributes(int attr, float val);
extern void	SetDefaultOrbitalAttributes(int attr);
extern void	SetPsi1Path (char *);
extern void	SetPsiconPath (char *);

#endif
