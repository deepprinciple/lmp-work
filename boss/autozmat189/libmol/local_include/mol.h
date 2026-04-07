#ifndef __MOL_H_
#define __MOL_H_

/* atom flags */
#define A_HETATM	0x1
#define A_BREAK		0x2
#define A_HIDDEN	0x4
#define A_DOMAIN	0x8

#define A_MARKED	0x1000000
#define A_SELECTED	0x2000000
#define A_PICKED	0x4000000
#define A_NEW		0x8000000

/* atom->label2D_align */
enum	{A_LABEL_LEFT, A_LABEL_CENTER, A_LABEL_RIGHT, A_LABEL_TOP, A_LABEL_BOTTOM, A_LABEL_POS_MAX};

/* bond flags */
#define B_BOLD		0x1000
#define B_HASHED	0x2000
#define B_OPEN_WEDGED	0x4000
#define B_FILL_WEDGED	0x8000
#define B_HASH_WEDGED	0x10000
#define B_WAVED		0x20000
#define B_ARROWED	0x40000

#define B_INTERSECT	0x100000
#define B_BACKWARD	0x200000
#define B_POS_DIR	0x400000
#define B_NEG_DIR	0x800000

#define B_MARKED	0x1000000
#define B_SELECTED	0x2000000
#define B_PICKED	0x4000000

#define B_BONDSHAPE	(B_BOLD|B_HASHED|B_OPEN_WEDGED|B_FILL_WEDGED|B_HASH_WEDGED|B_WAVED|B_ARROWED)
#define BONDSHAPE(flags)	(flags & B_BONDSHAPE)

/* bond types */
#define B_PARTIAL	0
#define B_SINGLE	1
#define B_DOUBLE	2
#define B_TRIPLE	3
#define B_AROMA		4
#define B_DBLAROMA	5
#define B_HBOND		7

/* residue flags */
#define R_HETATM	0x1
#define R_HIDDEN	0x2
#define R_MARKED	0x1000000
#define R_SELECTED	0x2000000

/* residue structure flags */
#define R_HELIX		0x1
#define R_SHEET		0x2
#define R_TURN		0x4
#define R_BREAK		0x100
#define R_SSINIT	0x200
#define R_SSTERM	0x400

#define SSTYPES		(R_HELIX|R_SHEET|R_TURN)

/* SS flags */
#define S_HIDDEN	0x1
#define S_MARKED	0x1000000
#define S_SELECTED	0x2000000

/* chain flags */
#define C_HIDDEN	0x1
#define C_MARKED	0x1000000
#define C_SELECTED	0x2000000

#define MAXNEIGHBOR	10

/*
typedef struct {
	float	x, y, z;
	} Vec;
*/

typedef struct {
	float	x1, y1, z1, x2, y2, z2;
	} BBox3;

typedef struct {
	float	x1, y1, x2, y2;
	} BBox;

typedef struct {
	unsigned long	flags;	/* label flags */
	int	precision;	/* label precision */
	int	fontsize;	/* label font size */
	int	fontstyle;	/* label font style */
	} AtomLabel;

typedef struct {
	unsigned long	mflags;
	unsigned long	viewtype;
	AtomLabel	label;
	} Domain, *DomainPtr;

struct _Bond;
struct _Residue;
struct _Chain;

typedef struct _ATOM {
	unsigned long	flags;
	int	an;		/* atomic number */
	double	x, y, z;	/* coord */
	long	data;		/* scratch space */

	int	refno;
	int	serno;
	float	tempfac;
	char	altloc;

	int	col;		/* color index */
	float	label1, label2;
	int	usertag;

	int	type;
	double	charge;
	double	sigma;
	double	epsilon;
	char	ambertype[3];

	int	nneighbors;
	struct _ATOM	*nbatom[MAXNEIGHBOR];
	struct _Bond	*nbbond[MAXNEIGHBOR];

	struct _Residue	*residue;	/* parent residue */

	/* for 2D drawing */
	float	px, py;		/* projected points */
	int	label2D_align;	/* alignment of 2D atom label */
	char	label2D[32];	/* 2D atom label */

	/* domain properties */
	Domain	domain;

	struct _ATOM	*next;
	struct _ATOM	*prev;
	} ATOM, *AtomPtr;

/*
	vp1               vn2
	 | +             - |
	atom1 ----------- atom2
	 | -             + |
	vn1               vp2
*/

typedef struct _Bond {
	unsigned long	flags;
	int	type;			/* bond type */
	ATOM	*atom1, *atom2;

	/* for 2D drawing */
	float	vp1[2], vn2[2], vp2[2], vn1[2];	/* projected vertices for drawing bond */

	/* for H-bond */
	ATOM	*CA1, *CA2;	/* source and destination alpha carbons */
	float	energy;
	int	offset;

	struct _Bond	*next;
	struct _Bond	*prev;
	} Bond, *BondPtr;

typedef struct _Residue {
	unsigned long	flags;		/* general flags */
	unsigned long	ss_flags;	/* structure flags */
	ATOM		*atom;		/* list of atoms */
	ATOM		*backbone;	/* backbone atom */
	float		normal[3];	/* normal vector */
	float		smooth[3];	/* smoothed backbone */
	int		seqno;		/* sequence number */
	int		refno;		/* residue index */

	BBox		bbox;		/* bounding box */

	struct _Chain	*chain;		/* parent chain */

	struct _Residue	*next;
	struct _Residue	*prev;
	} Residue, *ResiduePtr;

typedef struct _Chain {
	unsigned long	flags;
	Residue		*residue;	/* list of residues */
	Bond		*bond;		/* list of back bonds */
	char		id;		/* chain identifier */
	
	BBox		bbox;		/* bounding box */

	struct _Chain	*next;
	struct _Chain	*prev;
	} Chain, *ChainPtr;

/* secondary structure */
typedef struct _SS {
	unsigned long	flags;

	int	type;		/* structure type */
	char	chain_id;	/* chain ID */
	int	init;		/* sequence number of initial residue (for immediate use only) */
	int	term;		/* sequence number of terminal residue (for immediate use only) */
	Residue	*init_res;	/* initial residue */
	Residue	*term_res;	/* terminal residue */

	BBox	bbox;		/* bounding box */

	struct _SS	*next;
	} SS, *SSPtr;

enum	ParamType {PARAM_OPLS, PARAM_MMODEL, PARAM_MM2, PARAM_MM3};
typedef struct _Mol {
	unsigned long	flags;
	int		param_type;	/* parameter type (OPLS, MacroModel, MM2, etc) in atom->type */
	int		header;
	Bond		*ssbond;	/* list of S-S bonds */
	Bond		*hbond;		/* list of H-bonds */
	Chain		*chain;		/* list of chains */
	Bond		*bond;		/* list of bonds */
	SS		*ss;		/* list of secondary structures */

	BBox		bbox;		/* boudning box in world coord */
	BBox3		bbox3;		/* 3D bbox in atom coord */

	struct _Mol	*next;
	struct _Mol	*prev;
	} Mol, *MolPtr;

#define ForEachMol(mlist, m)	for(m=mlist;m;m=m->next)
#define ForEachChain(clist, c)  for(c=clist;c;c=c->next)
#define ForEachRes(rlist, r)    for(r=rlist;r;r=r->next)
#define ForEachAtom(alist, a)   for(a=alist;a;a=a->next)
#define ForEachBond(blist, b)   for(b=blist;b;b=b->next)
#define ForEachChainRes(clist, c, r)    ForEachChain(clist,c) ForEachRes(c->residue,r)
#define ForEachChainResAtom(clist, c, r, a)     ForEachChain(clist,c) ForEachRes(c->residue,r) ForEachAtom(r->atom,a)
#define ForEachResAtom(rlist, r, a)     ForEachRes(rlist,r) ForEachAtom(r->atom,a)
#define ForEachSS(sslist, ss)	for(ss=sslist;ss;ss=ss->next)

#define UnmarkAllAtomsInChains(clist,c,r,a)	ForEachChainRes(clist, c, r) a->flags &= ~A_MARKED
#define UnmarkAllAtomsInResidues(rlist,r,a)	ForEachResAtom(rlist, r, a) a->flags &= ~A_MARKED
#define UnmarkAllAtoms(alist,a)	ForEachAtom(alist,a) a->flags &= ~A_MARKED

#ifndef SQUARE
#define SQUARE(x)	(x)*(x)
#endif

#define AROMATICBOND(b)		(b==B_AROMA||b==B_DBLAROMA)
#define DOUBLEBOND(b)		(b==B_AROMA||b==B_DOUBLE)
#define TRIPLEBOND(b)		(b==B_DBLAROMA||b==B_TRIPLE)
#define ATOMSQDIST(a1,a2)	(SQUARE(a1->x-a2->x)+SQUARE(a1->y-a2->y)+SQUARE(a1->z-a2->z))
#define ATOMDIST(a1,a2)		sqrt(SQUARE(a1->x-a2->x)+SQUARE(a1->y-a2->y)+SQUARE(a1->z-a2->z))

/***************************************************/
/* pdb.C */
extern MolPtr	FLoadPdbFile(FILE *fp);
extern void	FPrintPdbMol(FILE *fp, MolPtr mlist);
extern ResiduePtr	CreatePdbResidue (ChainPtr *, ChainPtr *, ResiduePtr *, PdbResidue *, int);
extern AtomPtr	CreatePdbAtom (ChainPtr *, ChainPtr *, ResiduePtr *, AtomPtr *, PdbAtom *, BondPtr *, int hetatm);

/* orbital.C */
extern MolPtr	CreateMolFromOrbital(Orbital *orb);

/****************************************
 ******    mol.C    *********************
 ****************************************/

extern void	EnterChain (ChainPtr New, ChainPtr *top);
extern ChainPtr	NewChain (void);
extern ChainPtr	EnterNewChain (ChainPtr *top);
extern ChainPtr	DeleteChain (ChainPtr del, ChainPtr *top);
extern void	FreeChainEntry (ChainPtr *entry);
extern void	FreeChain (ChainPtr *top);
extern void	EnterResidue (ResiduePtr New, ResiduePtr *top);
extern ResiduePtr	NewResidue (void);
extern ResiduePtr	EnterNewResidue (ResiduePtr *top);
extern ResiduePtr	DeleteResidue (ResiduePtr del, ResiduePtr *top);
extern void	FreeResidueEntry (ResiduePtr *entry);
extern void	FreeResidue (ResiduePtr *top);
extern void	EnterMol (MolPtr New, MolPtr *top);
extern MolPtr	NewMol (void);
extern MolPtr	EnterNewMol (MolPtr *top);
extern MolPtr	DeleteMol (MolPtr del, MolPtr *top);
extern void	FreeMolEntry (MolPtr *entry);
extern void	FreeMol (MolPtr *top);
extern void	EnterAtom (AtomPtr New, AtomPtr *top);
extern AtomPtr	NewAtom (void);
extern AtomPtr	EnterNewAtom (AtomPtr *top);
extern AtomPtr	DeleteAtom (AtomPtr del, AtomPtr *top);
extern void	FreeAtomEntry (AtomPtr *entry);
extern void	FreeAtom (AtomPtr *top);
extern int	SearchAtomInChainByPointer (ChainPtr clist, AtomPtr atom);
extern AtomPtr	SearchAtomInChainByNumber (ChainPtr clist, int num);
extern AtomPtr	SearchAtomInChainBySerno (ChainPtr clist, int serno);
extern AtomPtr	NextAtom (AtomPtr atom);
extern AtomPtr	FirstAtom (ChainPtr chain);
extern void	InitAtomData (AtomPtr atom, ResiduePtr res, int an, double x, double y, double z);
extern int	CountResidueInChains (ChainPtr clist);
extern int	CountAtomInChains (ChainPtr clist);
extern int	CountAtomInResidues (ResiduePtr rlist);
extern int	CountAtom (AtomPtr atomlist);

extern void	EnterBond (BondPtr New, BondPtr *top);
extern BondPtr	NewBond (void);
extern BondPtr	EnterNewBond (BondPtr *top);
extern BondPtr	DeleteBond (BondPtr del, BondPtr *top);
extern void	FreeBondEntry (BondPtr *entry);
extern void	FreeBond (BondPtr *top);
extern int	CountBond (BondPtr bondlist);

extern SS	*NewSS (void);
extern void	EnterSS (SS *New, SS **top);
extern SS	*EnterNewSS (SS **top);
extern SS	*DeleteSS (SS *del, SS **top);
extern void	RemoveSS (SS *del, SS **top);
extern void	FreeSS (SS **top);
extern BondPtr	GetSernoBond (BondPtr bondlist, int serno1, int serno2);
extern BondPtr	GetAtomBonded (BondPtr bondlist, AtomPtr atom1, AtomPtr atom2);
extern BondPtr	GetBondBetweenAtoms (AtomPtr atom1, AtomPtr atom2);
extern BondPtr	GetNeighborBond (AtomPtr atom1, AtomPtr atom2);
extern AtomPtr	GetSernoAtom (ChainPtr chainlist, int serno);
extern int	ProcessBond (MolPtr m, int serno1, int serno2);
extern void	CalcBonds (MolPtr m);
extern void	CalcBondsWithinResidues (MolPtr m);
extern void	CalcBondsBetweenResidues (MolPtr m);
extern void	CalcBondsBetweenTwoResidues (Residue *res1, Residue *res2, BondPtr *bondlist_ret);
extern int	TestBonded (AtomPtr atom1, AtomPtr atom2);
extern void	SetMolBBox3 (MolPtr m);

extern int	DetermineStruct (MolPtr m);
extern int	CalcHbonds (MolPtr m);
extern int	FindAlphaHelix (MolPtr m, int pitch, int flag);

extern AtomPtr	FindCA (Residue *res);
extern Chain	*FindChainByID (Chain *clist, char id);
extern void	GetNeighbor (MolPtr mol);
extern SS	*MakeSSList (Chain *chainlist);
extern SS	*SearchSSList (SS *sslist, Residue *res);
extern void	FPrintSSEntry (FILE *fp, SS *ss);
extern void	FPrintSS (FILE *fp, SS *sslist);
extern void	FPrintAtomEntry (FILE *fp, AtomPtr atom);
extern void	FPrintAtom (FILE *fp, AtomPtr atomlist);

extern ChainPtr	CopyChainEntry(ChainPtr);
extern ChainPtr	CopyChain(ChainPtr);
extern ResiduePtr	CopyResidueEntry(ResiduePtr);
extern ResiduePtr	CopyResidue(ResiduePtr);
extern MolPtr	CopyMol(MolPtr);
extern AtomPtr	CopyAtomEntry(AtomPtr);
extern AtomPtr	CopyAtom(AtomPtr);
extern BondPtr	CopyBondEntry(BondPtr);
extern BondPtr	CopyBond(BondPtr);
extern SSPtr	CopySS(SSPtr, ChainPtr, ChainPtr);

extern void	ClearBBox (BBox *bbox);
extern void	ZeroBBox (BBox *bbox);
extern void	SetBBox (BBox *bbox, float x1, float y1, float x2, float y2);
extern void	UpdateBBox (BBox *bbox, float x1, float y1, float x2, float y2);
extern void	SubtBBox (BBox *bbox, float x, float y);
extern void	ClearBBox3 (BBox3 *bbox);
extern void	SetBBox3 (BBox3 *bbox, float x1, float y1, float z1, float x2, float y2, float z2);
extern void	UpdateBBox3 (BBox3 *bbox, float x1, float y1, float z1, float x2, float y2, float z2);
extern void	UpdateBBox3Point (BBox3 *bbox, float x, float y, float z);

extern float	BondType2Order (int type);
extern float	GuessBondLength(AtomPtr, AtomPtr, int);
extern int	GuessBondType(AtomPtr, AtomPtr);
extern char	*GetFormula (ChainPtr clist);

extern void	FPrintMol (FILE *fp, MolPtr mol);
extern MolPtr	FLoadMol(FILE *fp, int version);

extern void	AddAtomListToMol (MolPtr mol, ListPtr atomlist, int prepend);

extern AtomPtr	CreateDummy(double *, double *, double *);
extern ListPtr	MakeAllDummyAtoms(ChainPtr clist);
extern ListPtr	MakeDummyAtAtom (AtomPtr atom, int serno);
extern AtomPtr	MakeDummyAtomsAtCenter (ChainPtr clist, int ndummies, ListPtr *dummylist);
extern int	PrependDummyAtomsToMol (MolPtr mol, int ndummies);

extern void	RemoveMarkedAtomsFromMol (MolPtr mol);
extern void	RemoveDummiesFromMol (MolPtr mol);
extern void	FixAtomSerno (MolPtr mol);

extern void	FixConnectionTable (MolPtr mol);

#endif
