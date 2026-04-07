#ifndef __ZMATRIX_H_
#define __ZMATRIX_H_

/* Z-Matrix types */
enum	{Z_TYPE_GAUSS, Z_TYPE_JAGUAR, Z_TYPE_BOSS};

#define Z_GAUSS_DEFAULT_TITLE	"GAUSSIAN Z-MATRIX"
#define Z_MOPAC_DEFAULT_TITLE	"MOPAC Z-MATRIX"
#define Z_BOSS_DEFAULT_TITLE	"BOSS Z-MATRIX"
#define Z_MCPRO_DEFAULT_TITLE	"MCPRO Z-MATRIX"

/* Zmat->flags */
#define Z_VAR_LENGTH	0x1
#define Z_VAR_ANGLE	0x2
#define Z_VAR_DIHED	0x4
#define Z_CARTESIAN	0x8
#define Z_CHIRAL_PLUS	0x10
#define Z_CHIRAL_MINUS	0x20
#define Z_NEG_DIHED	0x40
#define Z_POS_DIHED	0x80
#define Z_GRAD_LENGTH	0x100
#define Z_GRAD_ANGLE	0x200
#define Z_GRAD_DIHED	0x400
#define Z_COUNTER_POISE	0x800

#define Z_MARKED	0x1000000
#define Z_MARKED_LENGTH	0x2000000
#define Z_MARKED_ANGLE	0x4000000
#define Z_MARKED_DIHED	0x8000000

/* array size */
#define Z_LABEL_LEN	4	/* max length of atom label in BOSS Z-matrix */
#define Z_RES_LEN	3	/* max length of residue name in BOSS Z-matrix */
#define Z_NAME_LEN	32	/* max length of atom and variable names (as in GAUSSIAN Z-matrix) */

typedef struct _Zmat {
	unsigned long	flags;
	int	an;			/* atomic number */
	char	label[Z_LABEL_LEN+1];	/* atomic symbol or label */
	char	resname[Z_RES_LEN+1];	/* residue name */
	int	resseq;			/* residue sequence */

	int	itype;			/* initial atom type */
	int	ftype;			/* final   atom type */

	char	atomname[4][Z_NAME_LEN+1];	/* atom identifiers */
	char	varname[4][Z_NAME_LEN+1];	/* variable identifiers */

	int	zdef[3];		/* reference atoms */
	double	zval[3];		/* length, angle, dihedral */
	double	cart[3];		/* Cartesian coord */

	int	atom;			/* atom serial number */

	/* charge, sigma, and epsilon are used to store user-supplied values */
	double	charge;			/* charge */
	double	sigma;			/* sigma */
	double	epsilon;		/* epsilon */

	struct _Zmat	*next;
	struct _Zmat	*prev;
	} Zmat, *ZmatPtr;

/* GaussZmat->flags */
#define GAUSS_CARTESIAN		0x1
#define GAUSS_MOPAC_ZMAT	0x2

typedef struct {
	unsigned long	flags;
	char	percent[1024];
	char	route[1024];
	char	title[256];
	int	charge, multiplicity;
	ZmatPtr	zmatlist;

	char	*data;
	} GaussZmat, *GaussZmatPtr;

#define MopacZmat	GaussZmat
#define MopacZmatPtr	GaussZmatPtr

#define ForEachZmat(zlist,z)	for(z=zlist;z;z=z->next)

/* zmatrix.C */
#define UpdateZmatCartWithVal		ZmatToCartesian
#define UpdateZmatListCartWithVal	ZmatListToCartesian

extern void	EnterZmat (ZmatPtr newptr, ZmatPtr *top);
extern ZmatPtr	NewZmat (void);
extern ZmatPtr	EnterNewZmat (ZmatPtr *top);
extern ZmatPtr	CopyZmat (ZmatPtr oldlist);
extern ZmatPtr	AllocZmat (int n);
extern ZmatPtr	ReallocZmat (ZmatPtr top, int nitems);
extern void	PrintZmat (ZmatPtr top);
extern ZmatPtr	SearchZmat (ZmatPtr top, int i);
extern ZmatPtr	SearchZmatByAtomSerno (ZmatPtr top, int atom_serno);
extern int	SearchZmatByPointer (ZmatPtr top, ZmatPtr zmat);

extern int	FindZmatByName (ZmatPtr zmatlist, char *name);
extern ZmatPtr	DeleteZmat (ZmatPtr del, ZmatPtr *top);
extern void	InsertZmat (ZmatPtr *top, ZmatPtr node, ZmatPtr newlist);
extern void	FreeZmat (ZmatPtr *top);
extern void	FreeZmatEntry (ZmatPtr *entry);
extern ListPtr	CopyZmatList (ListPtr oldzmatlist);
extern ZmatPtr	LinkZmatList (ListPtr zmatlist);
extern void	BreakZmatList (ListPtr zmatlist);
extern void	FreeZmatList (ListPtr *zmatlist);
extern void	PrintZmatList (ListPtr zmatlist);
extern int	CheckInternalCoord (int z0, int z1, int z2, int z3, double len, double theta, double phi);
extern void	InternalToCartesian (double *a1, double *a2, double len, double *a3, double theta, double *a4, double phi);
extern int	ZmatToCartesian (ZmatPtr zmatlist);
extern int	ZmatListToCartesian (ListPtr zmatlist);
extern char	*GetSymbolFromZmatCenter (char *str);
extern int	ParseZmatCenter (char *str, ZmatPtr zmatlist, ZmatPtr zmat);
extern int	ParseZmatRef (char *str, ZmatPtr zmatlist, ZmatPtr zmat, int n);
extern MolPtr	CreateMolFromZmat (ZmatPtr zmatlist);
extern int	IsRedundantIntZmatVar (ZmatPtr zmatlist, int nzmat, ZmatPtr cur_zmat, int vartype);
extern int	CountIntZmatVar (ZmatPtr zmatlist);
extern int	IsRedundantCartZmatVar (ZmatPtr zmatlist, int nzmat, ZmatPtr cur_zmat, int n);
extern int	CountCartZmatVar (ZmatPtr zmatlist);

extern int	ConvertZmat (int src_zmat_type, unsigned long src_zmat, int dst_zmat_type, unsigned long *dst_zmat);
extern ZmatPtr	GetZmatOrderedByAtomSerno (ZmatPtr zmatlist);

extern void	UpdateZmatValWithCart (ZmatPtr zmat);
extern int	UpdateZmatWithMolCoord (ZmatPtr zmatlist, MolPtr mol);
extern int	UpdateZmatListWithMolCoord (ListPtr zmatlist, MolPtr mol);

/* flags for UpdateMolWithZmat and UpdateMolWithZmatList */
#define Z_UPDATE_ATOM_ATTR	0x1
extern int	UpdateMolWithZmat (ZmatPtr zmatlist, MolPtr mol, int flags);
extern int	UpdateMolWithZmatList (ListPtr zmatlist, MolPtr mol, int flags);

extern int	IsSequentialZmat (ZmatPtr zmatlist);

extern int	SetTemplateZmat (unsigned long zmat, int zmat_type);
extern unsigned long	ApplyTemplateZmat (MolPtr mol, int zmat_type, int user_flags);
extern int	GetTemplateZmatType (void);

#endif
