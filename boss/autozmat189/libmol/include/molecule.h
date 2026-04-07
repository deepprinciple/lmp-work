#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#ifdef HAVE_PLOT
#include "plot.h"
#endif /* HAVE_PLOT */
#ifndef __HASH_H_
#define __HASH_H_

extern int	LookUpStrTable(char *str);
extern void	CleanUpStrTable(void);
extern char	*GetStrValue(int n);

#endif
#ifndef __PDB_H_
#define __PDB_H_

#define PDB_ILEN	3	/* id */
#define PDB_DLEN	9	/* date */
#define PDB_RLEN	3	/* residue   GLY, ALA, ... */
#define PDB_ALEN	4	/* atom name C12, H24, ... */
#define PDB_PLEN	4	/* pdb name */

#define PDB_REMARK_LEN	59

#define PDB_ISIZE	(PDB_ILEN + 1)
#define PDB_RSIZE	(PDB_RLEN + 1)
#define PDB_ASIZE	(PDB_ALEN + 1)
#define PDB_PSIZE	(PDB_PLEN + 1)
#define PDB_DSIZE	(PDB_DLEN + 1)

#define PDB_REMARK_SIZE	(PDB_REMARK_LEN+1)

/* PDB record type in alphabetical order */
enum	{PDB_UNK=0,  PDB_ANISOU, PDB_ATOM,   PDB_AUTHOR, PDB_COMPND, PDB_CONECT, PDB_CRYST1, PDB_END,
	 PDB_ENDMDL, PDB_EXPDTA, PDB_FORMUL, PDB_FTNOTE, PDB_HEADER, PDB_HELIX,  PDB_HET,    PDB_HETATM,
	 PDB_JRNL,   PDB_MASTER, PDB_MODEL,  PDB_MTRIX,  PDB_OBSLTE, PDB_ORIGX,  PDB_REMARK, PDB_REVDAT,
	 PDB_SCALE,  PDB_SEQRES, PDB_SHEET,  PDB_SIGATM, PDB_SIGUIJ, PDB_SITE,   PDB_SOURCE, PDB_SPRSDE,
	 PDB_SSBOND, PDB_TER,    PDB_TURN,   PDB_TVECT,  PDB_MAX};

typedef struct	{	/* residue */
	char	name[PDB_RSIZE];		/* residue name */
	char	chain_id;			/* chain id */
	int	seq;				/* sequence */
	char	insert;				/* insertion code */
	} PdbResidue;

typedef struct {	/* strand register */
	char		name[PDB_ASIZE];	/* atom name */
	PdbResidue	residue;
	} PdbStrandReg;

typedef struct	{
	char	text[81];
	} PdbUnk;

typedef struct	{
	int		serial;
	char		name[PDB_ASIZE];
	char		altloc;
	PdbResidue	residue;
	int		u[6];			/* anisotropic temp factors */
	} PdbAnisou;
typedef PdbAnisou	PdbSiguij;

typedef struct	{
	int		serial;			/* atom serial # */
	char		name[PDB_ASIZE];	/* atom name */
	char		altloc;			/* alternate location indicator */
	PdbResidue	residue;		/* residue */
	double		x, y, z;		/* orthogonal coord in Angstrom */
	float		occupancy;		/* occupance */
	float		tempfac;		/* temperature scale factor */
	int		footnote;		/* footnote number */
	} PdbAtom;

typedef struct	{
	char	cont;
	char	text[61];
	} PdbAuthor;
typedef PdbAuthor	PdbCompnd;
typedef PdbAuthor	PdbExpdta;
typedef PdbAuthor	PdbSource;
typedef PdbAuthor	PdbJrnl	;

typedef struct	{
	int	serial;			/* serial # */
	int	covalent[4];		/* covalent bond connectivity */
	int	donor[2];		/* "serial" acts as h-bond donor */
	int	salt_minus;		/* "serial" has excess - charge */
	int	acceptor[2];		/* "serial" acts as h-bond acceptor */
	int	salt_plus;		/* "serial" has excess + charge */
	} PdbConect;

typedef struct	{
	float	a, b, c;		/* a, b, c in angstrom */
	float	alpha, beta, gamma;	/* alpha, beta, gamma in degree */
	char	spacegrp[12];		/* space group (left justified) */
	int	z;
	} PdbCryst1;

typedef struct	{
	int	component;		/* component number */
	char	hetid[PDB_RSIZE];	/* non-standard group (HET) id */
	int	cont;			/* continuation # */
	char	exclude;		/* exclude from M.W. calc */
	char	formula[52];		/* formula of non-std grp */
	} PdbFormul;

typedef struct	{
	int	num;
	char	text[60];
	} PdbFtnote;
typedef PdbFtnote	PdbRemark;

typedef struct	{
	char	funcclass[41];		/* functional classification */
	char	date[PDB_DSIZE];	/* date of deposition */
	char	id[PDB_PSIZE];		/* identification code */
	} PdbHeader;

typedef struct	{
	int		serial;		/* serial # (helix #) */
	char		id[PDB_ISIZE];	/* helix id (right justified) */
	PdbResidue	initial;	/* initial residue of helix */
	PdbResidue	terminal;	/* terminal residue of helix */
	int		helixclass;	/* class of helix */
	char		comment[31];	/* comment */
	} PdbHelix;

typedef struct	{
	PdbResidue	residue;	/* non-std grp (heterogen) residue */
	int		natoms;		/* # of atoms in non-std grp */
	char		text[41];	/* text */
	} PdbHet;
typedef PdbAtom	PdbHetatm;
typedef PdbAtom	PdbSigatm;

typedef struct	{
	int	nremark;		/* # of REMARK */
	int	nftnote;		/* # of FTNOTE */
	int	nhet;			/* # of HET */
	int	nhelix;			/* # of HELIX */
	int	nsheet;			/* # of SHEET */
	int	nturn;			/* # of TURN */
	int	nsite;			/* # of NSITE */
	int	ntransform;		/* # of (ORIGX + SCALE + MTRIX) */
	int	natom;			/* # of (ATOM + HETATM) */
	int	nter;			/* # of TER */
	int	nconect;		/* # of CONECT */
	int	nseqres;		/* # of SEQRES */
	} PdbMaster;

typedef struct	{
	int	row;
	int	serial;
	float	m1, m2, m3, v;
	int	given;
	} PdbMtrix;

typedef struct	{
	int	cont;
	char	date[PDB_DSIZE];	/* date this entry was replaced */
	char	old_id[PDB_PSIZE];	/* obsolete id */
	char	new_id[8][PDB_PSIZE];	/* new id's */
	} PdbObslte;

typedef struct	{
	int	row;
	float	o1, o2, o3, t;
	} PdbOrigx;

typedef struct	{
	int	mod;			/* modification # */
	int	cont;			/* continuation field */
	char	date[PDB_DSIZE];	/* date */
	char	id[8];			/* id name used for correction */
	char	modtype;		/* modification type */
	char	correct[31];		/* record types that were corrected */
	} PdbRevdat;

typedef struct	{
	int	row;
	float	s1, s2, s3, u;
	} PdbScale;

typedef struct	{
	int	serial;			/* serial # of SEQRES record for current chain */
	char	chain_id;		/* chain id */
	int	nresidues;		/* # of residues in this chain */
	char	names[13][PDB_RSIZE];	/* residue name */
	} PdbSeqres;

typedef struct	{
	int		strand;		/* strand number */
	char		id[PDB_ISIZE];	/* strand id */
	int		nstrands;	/* number of strands */
	PdbResidue	initial;	/* initial residue */
	PdbResidue	terminal;	/* terminal residue */
	int		sense;		/* sense of this strand */
	PdbStrandReg	current;	/* in current strand registration */
	PdbStrandReg	prev;		/* in previous strand registration */
	} PdbSheet;

typedef struct	{
	int		seq;		/* sequence # */
	char		id[PDB_ISIZE];	/* site id (right justified) */
	int		nresidues;	/* # of residues comprising site */
	PdbResidue	residue1;	/* 1st residue */
	PdbResidue	residue2;
	PdbResidue	residue3;
	PdbResidue	residue4;
	} PdbSite;

typedef struct	{
	int	cont;			/* continuation field */
	char	date[PDB_DSIZE];	/* date this entry superseded the old */
	char	id[PDB_PSIZE];		/* id of this entry */
	char	old_id[8][PDB_PSIZE];	/* id's of entries which are being replaced */
	} PdbSprsde;

typedef struct	{
	int		seq;
	PdbResidue	residue1;
	PdbResidue	residue2;
	char		comment[31];
	} PdbSsbond;

typedef struct	{
	int		serial;
	PdbResidue	residue;
	} PdbTer;

typedef struct	{
	int		seq;			/* sequence # (turn #) */
	char		id[PDB_ISIZE];		/* turn id */
	PdbResidue	residue1;		/* residue i */
	PdbResidue	residue2;		/* residue i + 3 */
	char		comment[31];
	} PdbTurn;

typedef struct	{
	int	serial;				/* serial # */
	float	t1, t2, t3;			/* translation vector */
	char	comment[31];			/* comment */
	} PdbTvect;

typedef struct	{
	int	serial;				/* model serial # */
	} PdbModel;

typedef union	{
	PdbUnk		unk;
	PdbAnisou	anisou;
	PdbAtom		atom;
	PdbAuthor	author;
	PdbCompnd	compnd;
	PdbConect	conect;
	PdbCryst1	cryst1;
	PdbFormul	formul;
	PdbFtnote	ftnote;
	PdbHeader	header;
	PdbHelix	helix;
	PdbHet		het;
	PdbHetatm	hetatm;
	PdbJrnl		jrnl;
	PdbMaster	master;
	PdbMtrix	mtrix;
	PdbObslte	obslte;
	PdbOrigx	origx;
	PdbRemark	remark;
	PdbRevdat	revdat;
	PdbScale	scale;
	PdbSeqres	seqres;
	PdbSheet	sheet;
	PdbSigatm	sigatm;
	PdbSiguij	siguij;
	PdbSite		site;
	PdbSource	source;
	PdbSprsde	sprsde;
	PdbSsbond	ssbond;
	PdbTer		ter;
	PdbTurn		turn;
	PdbTvect	tvect;
	PdbModel	model;
	} PdbRec, *PdbRecPtr;

typedef struct _PdbRecord {
	int		type;		/* PDB record type (PDB_ATOM, PDB_HETATM, etc) */
	PdbRec		rec;		/* union that holds record data */
	struct _PdbRecord	*next;	/* pointer to next PdbRecord structure */
	} PdbRecord, *PdbRecordPtr;

extern PdbRecord	readpdbrecord (FILE *);
extern int	getpdbrecordtype (char *);
extern int	scanpdbrecord (PdbRecord *, char *, char *);
extern char	*getpdbrecordscanfmt (int type);
extern char	*getpdbrecordprintfmt (int type);
extern void	fprintpdbrecord (PdbRecord *, FILE *);

#endif	/* PDB_H */
#ifndef __STRSCAN_H_
#define __STRSCAN_H_

#define DATA_INT	'i'
#define DATA_FLOAT	'f'
#define DATA_STR	's'

extern char	*strscan (char *buf, char *fmt, ...);
extern int	filescan (FILE *fp, char *fmt, ...);
extern int	IsAlphabet (char *str);
extern int	IsAlphaNumeric (char *str, int sign_flag);
extern int	IsInt (char *thestr);
extern int	IsNumber (char *thestr);
extern int	IsFloat (char *str);
extern int	GetDataType (char *thestr);
#endif
#ifndef __LIST_H_
#define __LIST_H

/* masks for QsortList */
#define M_DATA_1	0x1
#define M_DATA_2	0x2
#define M_DATA_3	0x4
#define M_DATA_4	0x8
#define M_DATA_5	0x10
#define M_DATA_6	0x20
#define M_DATA_7	0x40
#define M_DATA_8	0x80

/**************************************
* aliases for List structure elements *
***************************************/
/* general */
#define L_DATA_1	data1
#define L_DATA_2	data2
#define L_DATA_3	data3
#define L_DATA_4	data4
#define L_DATA_5	data5
#define L_DATA_6	data6
#define L_DATA_7	data7
#define L_DATA_8	data8
#define L_VAL_1		val1
#define L_VAL_2		val2
#define L_VAL_3		val3
#define L_VAL_4		val4
#define L_VAL_5		val5
#define L_VAL_6		val6
#define L_VAL_7		val7
#define L_VAL_8		val8
#define L_STR_1		str1
#define L_STR_2		str2
#define L_STR_3		str3
#define L_STR_4		str4
#define L_STR_5		str5
#define L_STR_6		str6
#define L_STR_7		str7
#define L_STR_8		str8

/* Z-matrix */
#define L_ZMAT		data1

/* boss atom types */
#define L_TYPE		data1
#define L_AN		data2
#define L_CHARGE	val1
#define L_SIGMA		val2
#define L_EPSILON	val3
#define L_AII		val4
#define L_CII		val5
#define L_AMBER		str1
#define L_COMMENT	str8

/* boss torsional parameters (L_TYPE = data1, L_COMMENT = str8) */
#define L_N1		data2
#define L_N2		data3
#define L_N3		data4
#define L_V0		val1
#define L_V1		val2
#define L_V2		val3
#define L_V3		val4

#define L_P1		val5
#define L_P2		val6
#define L_P3		val7
#define L_ETORSION	val8

/* in comments of torsion parameters, [A1-A2-A3-A4] represents the
corresponding atom types */
#define L_A1		data5
#define L_A2		data6
#define L_A3		data7
#define L_A4		data8

/* boss bond stretching and angle bending (for strbnd file) */
#define L_K		val1
#define L_EQUIL		val2
#define L_AMBER1	str1
#define L_AMBER2	str2
#define L_AMBER3	str3
#define L_AMBER4	str4

#define L_ESTRETCH	val8
#define L_EBEND		val8

/* 1,4 and nonboneded interactions */
#define L_Q1		val1
#define L_Q2		val2
#define L_AII1		val3
#define L_AII2		val4
#define L_CII1		val5
#define L_CII2		val6
#define L_ENB_ELEC	val7
#define L_ENB_VDW	val8

/* for calculation of energy (store pointers to atoms) */
#define L_ATOMPTR1	data5
#define L_ATOMPTR2	data6
#define L_ATOMPTR3	data7
#define L_ATOMPTR4	data8

/* geometry variations */
#define L_CENTER	data1
#define L_PARAM		data2
#define L_VAL		val1

/* variable bonds  (L_CENTER = data1) */
/* variable angles (L_CENTER = data1) */

/* additional bonds, bond angles,  */
#define L_ATOM1	data1
#define L_ATOM2	data2
#define L_ATOM3	data3
#define L_ATOM4	data4

/* harmonic constraints (L_ATOM1 = data1, L_ATOM2 = data2) */
#define L_RIJ0	val1
#define L_K0	val2
#define L_K1	val3
#define L_K2	val4

/* variable dihedral angles (L_CENTER = data1) */
#define L_ITYPE	data5
#define L_FTYPE	data6
#define L_IMPROPER	data8
#define L_RANGE	val1

/* additional dihedral angles
L_ATOM1 = data1, L_ATOM2 = data2, L_ATOM3 = data3, L_ATOM4 = data4,
L_ITYPE = data5, L_FTYPE = data6, L_IMPROPER = data8 */

/* domain definitions */
#define L_FIRST_DOM1	data1
#define L_LAST_DOM1	data2
#define L_FIRST_DOM2	data3
#define L_LAST_DOM2	data4

/* conformational search (L_CENTER = data1, L_PARAM = data2) */
#define L_RMIN	val1
#define L_RMAX	val2

/* local heating */
#define L_RESIDUE	data1

/* Gaussian Z-matrix variables list (var=val) */
#define L_VARNAME	str1
#define L_VARVAL	val1
#define L_VARFLAGS	data1

/* data structure for scratch space. */
#define SL	128	/* string length */
typedef struct _List {
	size_t	data1,    data2,    data3,    data4,    data5,    data6,    data7,    data8;
	float	val1,     val2,     val3,     val4,     val5,     val6,     val7,     val8;
	char	str1[SL], str2[SL], str3[SL], str4[SL], str5[SL], str6[SL], str7[SL], str8[SL];
	struct _List	*next;
	} List, *ListPtr;


#define ForEachList(list,l)	for(l=list;l;l=l->next)

#define NEW_LIST(type) \
	(type *) calloc(1, sizeof(type));

#define APPEND_LIST(p,tmp,top) \
	if (!*top) *top = p; \
	else { \
		for(tmp=(*top);tmp->next;tmp=tmp->next); \
		tmp->next = p; \
	}

#define ENTER_LIST	APPEND_LIST

#define ENTER_NEW_LIST(p,tmp,top,type)	\
	if ((p = NEW_LIST(type))) APPEND_LIST(p,tmp,top)

#define FREE_LIST(p,tmp,top) \
	if (!top || !*top) return; \
	p = *top; \
	while (p) {tmp = p->next; free((void *)p); p = tmp;} \
	*top = NULL;

#define PRINT_LIST(p,top) \
	PRINT("List:\n"); \
	ForEachList(top,p) PRINT(" %d\n", p);

#define FOR_EACH_LIST(top,p) \
	for(p=top;p;p=p->next)

#define DELETE_LIST(p,tmp,top) \
	if (p && tmp && top && *top) { \
		if (p == *top) *top = p->next; \
		else { \
			for(tmp=(*top);tmp->next!=p;tmp=tmp->next); \
			tmp->next = p->next; \
		} \
	}


extern ListPtr	NewList ();
extern ListPtr	EnterNewList (ListPtr *top);
extern void	EnterList (ListPtr newptr, ListPtr *top);
extern void	InsertList (ListPtr *top, ListPtr node, ListPtr newlist);
extern ListPtr	CopyList (ListPtr oldlist);
extern void	PrintList (ListPtr top);
extern void	DeleteList (ListPtr del, ListPtr *top);
extern void	RemoveList (ListPtr del, ListPtr *top);
extern void	FreeList (ListPtr *top);
extern void	FreeListEntry (ListPtr *list);
extern ListPtr	SearchList (ListPtr top, int i);
extern ListPtr	FindListData (ListPtr top, int data);
extern int	FindListDataEntry (ListPtr top, int data);

extern int	CompareList (ListPtr list1, ListPtr list2, unsigned int flags);
extern void	SwapList (ListPtr listi, ListPtr listj);
extern void	QsortList (ListPtr toplist, int left, int right, int mode);
extern int	CountList (ListPtr top);

#endif


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
#ifndef __ELEM_DATA_H_
#define __ELEM_DATA_H_

typedef struct {
	char	symbol[4];	/* atomic symbol */
	char	name[13];	/* atomic name */
	float	weight;		/* atomic weight */
	float	weight_isotope;	/* atomic weight of the most abundant isotope */
	float	bsrad;		/* Bragg-Slater radius */
	float	vdwrad;		/* VDW radius */
	int	col;		/* CPK color */
	} ElemData;

typedef struct {
	char	*name;
	int	refno;
	} RefTab;

/* color indices */
enum	{
	CPK_LIGHT_GREY, CPK_SKY_BLUE,  CPK_RED,        CPK_YELLOW,
	CPK_BLACK,      CPK_PINK,      CPK_GOLDEN_ROD, CPK_BLUE,
	CPK_ORANGE,     CPK_DARK_GREY, CPK_BROWN,      CPK_PURPLE,
	CPK_DEEP_PINK,  CPK_GREEN,     CPK_FIRE_BRICK, CPK_MID_GREEN,
	CPK_MAX};

#define MAXELEM	104

enum	{ELEM_X=0, ELEM_H, ELEM_HE, ELEM_LI, ELEM_BE, ELEM_B, ELEM_C, ELEM_N,	/*  0 -  7 */
	ELEM_O, ELEM_F, ELEM_NE, ELEM_NA, ELEM_MG, ELEM_AL, ELEM_SI, ELEM_P,	/*  8 - 15 */
	ELEM_S, ELEM_CL,							/* 16 - 17 */
	ELEM_LP=98, ELEM_GH=99};

#define ELEM_DATA {\
	{"X ", "Dummy       ",   0.00000,   0.00000, 0.35, 1.20, CPK_GOLDEN_ROD},	/* 0 */ \
	{"H ", "Hydrogen    ",   1.00790,   1.00783, 0.35, 1.20, CPK_LIGHT_GREY}, \
	{"He", "Helium      ",   4.00260,   4.00260, 0.50, 1.40, CPK_YELLOW}, \
	{"Li", "Lithium     ",   6.94100,   7.01600, 1.45, 1.82, CPK_GOLDEN_ROD}, \
	{"Be", "Beryllium   ",   9.01218,   9.01218, 1.05, 1.40, CPK_GREEN}, \
	{"B ", "Boron       ",  10.81000,  11.00931, 0.85, 1.80, CPK_ORANGE}, \
	{"C ", "Carbon      ",  12.01100,  12.00000, 0.70, 1.70, CPK_BLACK}, \
	{"N ", "Nitrogen    ",  14.00670,  14.00307, 0.65, 1.55, CPK_BLUE}, \
	{"O ", "Oxygen      ",  15.99940,  15.99491, 0.60, 1.52, CPK_RED}, \
	{"F ", "Fluorine    ",  18.99840,  18.99840, 0.55, 1.47, CPK_GREEN}, \
	{"Ne", "Neon        ",  20.17900,  19.99244, 0.65, 1.54, CPK_ORANGE},	/* 10 */ \
	{"Na", "Sodium      ",  22.98977,  22.98977, 1.80, 2.27, CPK_ORANGE}, \
	{"Mg", "Magnesium   ",  24.30500,  23.98504, 1.50, 1.73, CPK_SKY_BLUE}, \
	{"Al", "Aluminium   ",  26.98154,  26.98154, 1.25, 2.40, CPK_SKY_BLUE}, \
	{"Si", "Silicon     ",  28.08550,  27.97693, 1.10, 2.10, CPK_RED}, \
	{"P ", "Phosphorus  ",  30.97376,  30.97376, 1.00, 1.80, CPK_PINK}, \
	{"S ", "Sulphur     ",  32.06000,  31.97207, 1.00, 1.80, CPK_ORANGE}, \
	{"Cl", "Chlorine    ",  35.45300,  34.96885, 0.96, 1.75, CPK_PINK}, \
	{"Ar", "Argon       ",  39.94800,  39.96240, 0.95, 1.88, CPK_YELLOW}, \
	{"K ", "Potassium   ",  39.09830,  38.96371, 2.20, 2.75, CPK_RED}, \
	{"Ca", "Calcium     ",  40.08000,  39.96259, 1.80, 2.00, CPK_RED},	/* 20 */ \
	{"Sc", "Scandium    ",  44.95590,  44.95592, 1.60, 2.30, CPK_GREEN}, \
	{"Ti", "Titanium    ",  47.90000,  45.94800, 1.40, 2.00, CPK_GREEN}, \
	{"V ", "Vanadium    ",  50.94150,  50.94400, 1.35, 1.90, CPK_GREEN}, \
	{"Cr", "Chromium    ",  51.99600,  51.94050, 1.40, 1.80, CPK_GREEN}, \
	{"Mn", "Manganese   ",  54.93800,  54.93810, 1.40, 1.75, CPK_GREEN}, \
	{"Fe", "Iron        ",  55.84700,  55.93490, 1.40, 1.70, CPK_GREEN}, \
	{"Co", "Cobalt      ",  58.93320,  58.93320, 1.35, 1.70, CPK_RED}, \
	{"Ni", "Nickel      ",  58.70000,  57.93530, 1.35, 1.63, CPK_GREEN}, \
	{"Cu", "Copper      ",  63.54600,  62.92980, 1.16, 1.40, CPK_RED}, \
	{"Zn", "Zinc        ",  65.38000,  63.92910, 1.35, 1.39, CPK_RED},	/* 30 */ \
	{"Ga", "Gallium     ",  69.73500,  68.92570, 1.35, 1.87, CPK_RED}, \
	{"Ge", "Germanium   ",  72.59000,  73.92190, 1.30, 2.20, CPK_RED}, \
	{"As", "Arsenic     ",  74.92160,  74.92160, 1.25, 1.85, CPK_GREEN}, \
	{"Se", "Selenium    ",  78.96000,  79.91650, 1.15, 1.90, CPK_MID_GREEN}, \
	{"Br", "Bromine     ",  79.90400,  79.91830, 1.15, 1.85, CPK_YELLOW}, \
	{"Kr", "Krypton     ",  83.80000,  83.91200, 1.15, 2.02, CPK_BLUE}, \
	{"Rb", "Rubidium    ",  85.46780,  84.91170, 1.10, 2.85, CPK_SKY_BLUE}, \
	{"Sr", "Strontium   ",  87.62000,  87.90560, 2.35, 2.10, CPK_RED}, \
	{"Y ", "Yttrium     ",  88.90590,  88.90540, 2.00, 2.20, CPK_YELLOW}, \
	{"Zr", "Zirconium   ",  91.22000,  89.90430, 1.80, 2.10, CPK_SKY_BLUE},	/* 40 */ \
	{"Nb", "Niobium     ",  92.90640,  92.90640, 1.55, 2.00, CPK_YELLOW}, \
	{"Mo", "Molybdenum  ",  95.94000,  95.94000, 1.45, 1.90, CPK_YELLOW}, \
	{"Tc", "Technetium  ",  98.90620,  98.90620, 1.45, 1.80, CPK_RED}, \
	{"Ru", "Ruthenium   ", 101.17000, 101.17000, 1.35, 1.75, CPK_RED}, \
	{"Rh", "Rhodium     ", 102.90550, 102.90550, 1.30, 1.70, CPK_GREEN}, \
	{"Pd", "Palladium   ", 106.40000, 106.40000, 1.35, 1.63, CPK_GREEN}, \
	{"Ag", "Silver      ", 107.86800, 107.86800, 1.40, 1.72, CPK_LIGHT_GREY}, \
	{"Cd", "Cadmium     ", 112.41000, 112.41000, 1.60, 1.58, CPK_GREEN}, \
	{"In", "Indium      ", 114.82000, 114.82000, 1.55, 1.93, CPK_GREEN}, \
	{"Sn", "Tin         ", 118.69000, 118.69000, 1.55, 2.40, CPK_RED},	/* 50 */ \
	{"Sb", "Antimony    ", 121.75000, 121.75000, 1.45, 2.10, CPK_RED}, \
	{"Te", "Tellurium   ", 127.60000, 127.60000, 1.45, 2.06, CPK_RED}, \
	{"I ", "Iodine      ", 126.90450, 126.90040, 1.40, 1.98, CPK_PURPLE}, \
	{"Xe", "Xenon       ", 131.30000, 131.30000, 1.40, 2.16, CPK_YELLOW}, \
	{"Cs", "Cesium      ", 132.90540, 132.90540, 1.50, 2.50, CPK_RED}, \
	{"Ba", "Barium      ", 137.33000, 137.33000, 1.50, 1.50, CPK_RED}, \
	{"La", "Lanthanum   ", 138.90550, 138.90550, 1.50, 1.50, CPK_RED}, \
	{"Ce", "Cerium      ", 140.12000, 140.12000, 1.50, 1.50, CPK_RED}, \
	{"Pr", "Praseodymium", 140.90770, 140.90770, 1.50, 1.50, CPK_RED}, \
	{"Nd", "Neodymium   ", 144.24001, 144.24001, 1.50, 1.50, CPK_RED},	/* 60 */ \
	{"Pm", "Promethium  ", 145.00000, 145.00000, 1.50, 1.50, CPK_RED}, \
	{"Sm", "Samarium    ", 150.39999, 150.39999, 1.50, 1.50, CPK_RED}, \
	{"Eu", "Europium    ", 151.96001, 151.96001, 1.50, 1.50, CPK_RED}, \
	{"Gd", "Gadolinium  ", 157.25000, 157.25000, 1.50, 1.50, CPK_RED}, \
	{"Tb", "Terbium     ", 158.92540, 158.92540, 1.50, 1.50, CPK_RED}, \
	{"Dy", "Dysprosium  ", 162.50000, 162.50000, 1.50, 1.50, CPK_RED}, \
	{"Ho", "Holmium     ", 164.90340, 164.90340, 1.50, 1.50, CPK_RED}, \
	{"Er", "Erbium      ", 167.25999, 167.25999, 1.50, 1.50, CPK_RED}, \
	{"Tm", "Thulium     ", 168.93420, 168.93420, 1.50, 1.50, CPK_RED}, \
	{"Yb", "Ytterium    ", 173.03999, 173.03999, 1.50, 1.50, CPK_RED},	/* 70 */ \
	{"Lu", "Luterium    ", 174.96700, 174.96700, 1.50, 1.50, CPK_RED}, \
	{"Hf", "Hafnium     ", 178.49001, 178.49001, 1.50, 1.50, CPK_RED}, \
	{"Ta", "Tantalum    ", 180.94791, 180.94791, 1.50, 1.50, CPK_RED}, \
	{"W ", "Tungsten    ", 183.85001, 183.85001, 1.50, 1.50, CPK_GREEN}, \
	{"Re", "Rhenium     ", 186.20700, 186.20700, 1.50, 1.50, CPK_GREEN}, \
	{"Os", "Osmium      ", 190.20000, 190.20000, 1.50, 1.50, CPK_GREEN}, \
	{"Ir", "Iridium     ", 192.22000, 192.22000, 1.50, 1.50, CPK_GREEN}, \
	{"Pt", "Platinum    ", 195.09000, 195.09000, 1.50, 1.50, CPK_LIGHT_GREY}, \
	{"Au", "Gold        ", 196.96651, 196.96651, 1.44, 1.44, CPK_GOLDEN_ROD}, \
	{"Hg", "Mercury     ", 200.59000, 200.59000, 1.50, 1.50, CPK_LIGHT_GREY},	/* 80 */ \
	{"Tl", "Thallium    ", 204.37000, 204.37000, 1.50, 1.50, CPK_RED}, \
	{"Pb", "Lead        ", 207.20000, 207.20000, 1.50, 1.50, CPK_RED}, \
	{"Bi", "Bismuth     ", 208.98039, 208.98039, 1.50, 1.50, CPK_RED}, \
	{"Po", "Polonium    ", 209.00000, 209.00000, 1.50, 1.50, CPK_RED}, \
	{"At", "Astatine    ", 210.00000, 210.00000, 1.50, 1.50, CPK_RED}, \
	{"Rn", "Radon       ", 222.00000, 222.00000, 1.50, 1.50, CPK_RED}, \
	{"Fr", "Francium    ", 223.00000, 223.00000, 1.50, 1.50, CPK_RED}, \
	{"Ra", "Radium      ", 226.02541, 226.02541, 1.50, 1.50, CPK_RED}, \
	{"Ac", "Actinium    ", 227.00000, 227.00000, 1.50, 1.50, CPK_RED}, \
	{"Th", "Thorium     ", 232.03810, 232.03810, 1.50, 1.50, CPK_RED},	/* 90 */ \
	{"Pa", "Protactinium", 231.03590, 231.03590, 1.50, 1.50, CPK_RED}, \
	{"U ", "Uranium     ", 238.02901, 238.02901, 1.50, 1.50, CPK_RED}, \
	{"Np", "Neptunium   ", 237.04820, 237.04820, 1.50, 1.50, CPK_RED}, \
	{"Pu", "Plutonium   ", 244.00000, 244.00000, 1.50, 1.50, CPK_RED}, \
	{"Am", "Americium   ", 243.00000, 243.00000, 1.50, 1.50, CPK_RED}, \
	{"Cm", "Curium      ", 247.00000, 247.00000, 1.50, 1.50, CPK_RED}, \
	{"Bk", "Berkelium   ", 247.00000, 247.00000, 0.35, 1.50, CPK_RED}, \
	{"LP", "LonePair    ", 251.00000, 251.00000, 0.35, 1.20, CPK_GOLDEN_ROD},	/* 98 = lone pair */ \
	{"Gh", "Ghost       ", 254.00000, 254.00000, 0.35, 1.20, CPK_LIGHT_GREY},	/* 99 = ghost */ \
	{"Fm", "Fermium     ", 257.00000, 257.00000, 1.00, 1.00, CPK_RED},	/* 100 */ \
	{"Md", "Mendelevium ", 258.00000, 258.00000, 1.00, 1.00, CPK_RED}, \
	{"No", "Nobelium    ", 259.00000, 259.00000, 1.00, 1.00, CPK_RED}, \
	{"Lr", "Lawrencium  ", 260.00000, 260.00000, 1.00, 1.00, CPK_RED}}

/* Originally, atomic numbers 98 and 99 are Cf and Es, but
since they're rarely found in organic and biochemistry,
their slots are used for lone pair and ghost atom.

{"Cf", "Californium ", 251.00000, 251.00000, 0.35, 1.20},
{"Es", "Einsteinium ", 254.00000, 254.00000, 0.35, 1.20},
{"Fm", "Fermium     ", 257.00000, 257.00000, 1.00, 1.00},
{"Md", "Mendelevium ", 258.00000, 258.00000, 1.00, 1.00},
{"No", "Nobelium    ", 259.00000, 259.00000, 1.00, 1.00},
{"Lr", "Lawrencium  ", 260.00000, 260.00000, 1.00, 1.00}}
*/


/* ordering of residues from RasMol */
enum	{RES_ALA, RES_GLY, RES_LEU, RES_SER, RES_VAL,	/* essential amino acids (RES_ALA~RES_HYP) */
	 RES_THR, RES_LYS, RES_ASP, RES_ILE, RES_ASN,
	 RES_GLU, RES_PRO, RES_ARG, RES_PHE, RES_GLN,
	 RES_TYR, RES_HIS, RES_CYS, RES_MET, RES_TRP,
	 RES_ASX, RES_GLX, RES_PCA, RES_HYP,
	 RES_HIP,					/* non standard amino acids */
	 RES_ACE, RES_AME, 				/* peptide caps */

	 RES_A,   RES_C,   RES_G,   RES_T,		/* DNA nucleotides (RES_A~RES_T) */
	 RES_U,   RES_PU,  RES_I,   RES_1MA, RES_5MC,	/* RNA nucleotides (RES_U~RES_PSU) */
	 RES_OMC, RES_1MG, RES_2MG, RES_M2G, RES_7MG,
	 RES_OMG, RES_YG,  RES_H2U, RES_5MU, RES_PSU,
	 RES_UNK, RES_FOR, RES_HOH, RES_DOD,	/* misc */
	 RES_SO4, RES_PO4, RES_NAD, RES_COA,
	 RES_MAX};

/* RefTab format */
#define RES_REF_TAB {     \
	{"ALA", RES_ALA}, \
	{"GLY", RES_GLY}, \
	{"LEU", RES_LEU}, \
	{"SER", RES_SER}, \
	{"VAL", RES_VAL}, \
	{"THR", RES_THR}, \
	{"LYS", RES_LYS}, \
	{"ASP", RES_ASP}, \
	{"ILE", RES_ILE}, \
	{"ASN", RES_ASN}, \
	{"GLU", RES_GLU}, \
	{"PRO", RES_PRO}, \
	{"ARG", RES_ARG}, \
	{"PHE", RES_PHE}, \
	{"GLN", RES_GLN}, \
	{"TYR", RES_TYR}, \
	{"HIS", RES_HIS}, \
	{"CYS", RES_CYS}, \
	{"MET", RES_MET}, \
	{"TRP", RES_TRP}, \
	{"ASX", RES_ASX}, \
	{"GLX", RES_GLX}, \
	{"PCA", RES_PCA}, \
	{"HYP", RES_HYP}, \
	                  \
	{"HIP", RES_HIP}, \
	                  \
	{"ACE", RES_ACE}, \
	{"AME", RES_AME}, \
	                  \
	{"  A", RES_A  }, \
	{"  C", RES_C  }, \
	{"  G", RES_G  }, \
	{"  T", RES_T  }, \
	{"  U", RES_U  }, \
	{" PU", RES_PU }, \
	{"  I", RES_I  }, \
	{"1MA", RES_1MA}, \
	{"5MC", RES_5MC}, \
	{"OMC", RES_OMC}, \
	{"1MG", RES_1MG}, \
	{"2MG", RES_2MG}, \
	{"M2G", RES_M2G}, \
	{"7MG", RES_7MG}, \
	{"OMG", RES_OMG}, \
	{" YG", RES_YG }, \
	{"H2U", RES_H2U}, \
	{"5MU", RES_5MU}, \
	{"PSU", RES_PSU}, \
	{"UNK", RES_UNK}, \
	{"FOR", RES_FOR}, \
	{"HOH", RES_HOH}, \
	{"DOD", RES_DOD}, \
	{"SO4", RES_SO4}, \
	{"PO4", RES_PO4}, \
	{"NAD", RES_NAD}, \
	{"COA", RES_COA}, \
                          \
	{"ADE", RES_A  }, \
	{"CPR", RES_PRO}, \
	{"CSH", RES_CYS}, \
	{"CSM", RES_CYS}, \
	{"CYH", RES_CYS}, \
	{"CYT", RES_C  }, \
	{"D2O", RES_DOD}, \
	{"GUA", RES_G  }, \
	{"H2O", RES_HOH}, \
	{"SOL", RES_HOH}, \
	{"SUL", RES_SO4}, \
	{"THY", RES_T  }, \
	{"TIP", RES_HOH}, \
	{"TRY", RES_TRP}, \
	{"URI", RES_U  }, \
	{"WAT", RES_HOH}}

/* BS = BackSpace, A after digit = Asterisk, '*' */

enum	{ATM_N,   ATM_CA,  ATM_C,   ATM_O,	/* Amino acid backbone */
	 ATM_CBS, ATM_OT,  ATM_S,   ATM_P,	/* Shapely Amino Backbone */
	 ATM_O1P, ATM_O2P, ATM_O5A, ATM_C5A,
	 ATM_C4A, ATM_O4A, ATM_C3A, ATM_O3A,
	 ATM_C2A, ATM_O2A, ATM_C1A,		/* Nucleic Acid Backbone  */
	 ATM_CA2, ATM_SG,  ATM_N1,  ATM_N2,
	 ATM_N3,  ATM_N4,  ATM_N6,  ATM_O2,
	 ATM_O4,  ATM_O6,			/* Nucleic Acid H-Bonding */
	 ATM_UNK,				/* Misc */
	 ATM_MAX};

#define ATM_REF_TAB { \
	{" N  ", ATM_N  }, \
	{" CA ", ATM_CA }, \
	{" C  ", ATM_C  }, \
	{" O  ", ATM_O  }, \
	{" C\\ ",ATM_CBS}, \
	{" OT ", ATM_OT }, \
	{" S  ", ATM_S  }, \
	{" P  ", ATM_P  }, \
	{" O1P", ATM_O1P}, \
	{" O2P", ATM_O2P}, \
	{" O5*", ATM_O5A}, \
	{" C5*", ATM_C5A}, \
	{" C4*", ATM_C4A}, \
	{" O4*", ATM_O4A}, \
	{" C3*", ATM_C3A}, \
	{" O3*", ATM_O3A}, \
	{" C2*", ATM_C2A}, \
	{" O2*", ATM_O2A}, \
	{" C1*", ATM_C1A}, \
	{" CA2", ATM_CA2}, \
	{" SG ", ATM_SG }, \
	{" N1 ", ATM_N1 }, \
	{" N2 ", ATM_N2 }, \
	{" N3 ", ATM_N3 }, \
	{" N4 ", ATM_N4 }, \
	{" N6 ", ATM_N6 }, \
	{" O2 ", ATM_O2 }, \
	{" O4 ", ATM_O4 }, \
	{" O6 ", ATM_O6 }, \
	{"    ", ATM_UNK}, \
                           \
	{" OP1", ATM_O1P}, \
	{"1OP ", ATM_O1P}, \
	{" OP2", ATM_O2P}, \
	{"2OP ", ATM_O2P}}

#define IsAmino(x)       ((x)<=RES_AME)
#define IsAminoNucleo(x) ((x)<=RES_PSU)
#define IsNucleo(x)      (((x)>=RES_A) && ((x)<=RES_PSU))
#define IsProtein(x)     (((x)<=RES_HYP) || (((x)>=RES_UNK) && ((x)<=RES_FOR)))
#define IsDNA(x)         (((x)>=RES_A) && ((x)<=RES_T))
#define IsSolvent(x)     (((x)>=RES_HOH) && ((x)<=RES_PO4))
#define IsWater(x)       (((x)==RES_HOH) || ((x)==RES_DOD))
#define IsIon(x)         (((x)==RES_SO4) || ((x)==RES_PO4))

#define IsPyrimidine(x)  (IsCytosine(x) || IsThymine(x))
#define IsPurine(x)      (IsAdenine(x) || IsGuanine(x))
#define IsRNA(x)         (IsNucleo(x) && !IsThymine(x))
#define NucleicCompl(x)  ((x)^3)


#define IsProline(x)     ((x)==RES_PRO)
#define IsCysteine(x)    ((x)==RES_CYS)
#define IsAdenine(x)     ((x)==RES_A)
#define IsCytosine(x)    ((x)==RES_C)
#define IsGuanine(x)     ((x)==RES_G)
#define IsThymine(x)     ((x)==RES_T)
#define IsNADGroup(x)    ((x)==RES_NAD)
#define IsCOAGroup(x)    ((x)==RES_COA)

#define IsAlphaCarbon(x)     ((x)==ATM_CA)
#define IsSugarPhosphate(x)  ((x)==ATM_P)
#define IsAminoBackbone(x)   ((x)<=ATM_O)
#define IsShapelyBackbone(x) ((x)<=ATM_P)
#define IsNucleicBackbone(x) (((x)>=ATM_P) && ((x)<=ATM_C1A))
#define IsShapelySpecial(x)  ((x)==ATM_CA2)
#define IsCysteineSulphur(x) ((x)==ATM_SG)

#define IsProteinAlphaCarbon(res, atm)  (IsProtein(res) && IsAlphaCarbon(atm))
#define IsNucleoSugarPhosphate(res,atm)	(IsNucleo(res) && IsSugarPhosphate(atm))

#define IsMetal(an)	(an==3 || an==4 || an==11 || an==12 || an==13 || (an>=19 && an<=32) || \
	((an>=37 && an <= 103) && !(an==52 || an==53 || an==54 || an==86)))
#define IsDummy(an)	(an<0||an>MAXELEM||an==ELEM_GH||an==ELEM_LP)

#define AN2ATOMREFNO(an)	GetSimpleAtomRefno(GetElemSymbol(an))
#define AN2ATOMNAME(an)		GetAtomName(AN2ATOMREFNO(an))

extern ElemData	*GetElemData(int an);
extern char	*GetElemSymbol(int an);
extern int	GetElemNumber(char *sym);
extern float	GetElemBSRadius(int an);
extern float	GetElemVdwRadius(int an);
extern float	GetElemWeight(int an);
extern float	GetElemIsotopeWeight(int an);
extern int	GetElemCPKColor(int an);

extern void	SetElemData(int an, ElemData *edata);
extern void	SetElemBSRadius(int an, float r);
extern void	SetElemVdwRadius(int an, float r);
extern void	SetElemWeight(int an, float w);
extern void	SetElemIsotopeWeight(int an, float w);
extern void	SetElemCPKColor(int an, int col);

extern int	TestMetalElem(int an);
extern short	*GetElemColor();
extern void	SetElemColor();

extern char	*GetAtomName(int refno);
extern char	*GetPdbAtomName (int atm_refno, int res_refno);
extern char	*MakePdbAtomName (char *str, int an);
extern char	*GetResName(int refno);
extern int	GetPdbAtomRefno(char *name, int res_refno);
extern int	GetSimpleAtomRefno(char *name);
extern int	GetAtomRefno(char *name);
extern int	GetResRefno(char *resname);
extern int	GetPdbAtomicNumber(int atm_refno, int res_refno);
extern Residue	*FindResBySeqno(Residue *reslist, int seqno);
extern AtomPtr	FindAtomInRes (Residue *res, int refno);
#endif
#ifndef __MOL_ERROR_H_
#define __MOL_ERROR_H_

/* error definitions */
enum {MERR_DEFAULT,
	MERR_NOFILE,
	MERR_MEM,
	MERR_EOF,
	MERR_FORMAT,

	MERR_INVALID_FLOAT,
	MERR_INVALID_NUMBER,
	MERR_UNDEFINED_CENTER,
	MERR_UNKNOWN_ATOMSYMBOL,
	MERR_MULTI_DEF_CENTER,
	MERR_ILLEGAL_CENTERNAME,

	MERR_INVALID_VAR_ASSIGN,
	MERR_VAR_NOT_FOUND,
	MERR_PARFORMAT,
	MERR_OPENFILE,
	MERR_RANGE,

	MERR_INTERNAL_COORD,
	MERR_INVALID_CENTER,
	MERR_INVALID_BONDLEN,
	MERR_INVALID_BONDANG,
	MERR_MAX};

extern int	merrno;

extern void	setmolerrorfunc (int (*func)(char *, ...));
extern int	setmolerror (char *format, ...);
extern int	catmolerror (char *format, ...);
extern void	clearmolerror (void);
extern void	setmolerrortext (char *format, ...);
extern int	setmolerrorno (int ierr);
extern int	getmolerrorno (void);
extern void	printmolerror (void);
extern void	setmoldebuglevel (int);
extern int	getmoldebuglevel (void);

#endif
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
#ifndef __AUTOZMAT_H_
#define __AUTOZMAT_H_

extern int	AssignZmatProp (ZmatPtr zmatlist);
extern ZmatPtr	NewZmatFromAtom (AtomPtr atom, int z1, int z2, int z3);

/* masks for flags */
#define Z_ADD_BGN_DUMMIES	0x1
#define Z_PUT_SERNO_INFO	0x2

extern ZmatPtr	CreateAutoZmatFromMol (MolPtr mol, int flags);

#endif
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
#ifndef __MOPACZMAT_H_
#define __MOPACZMAT_H_
extern MopacZmatPtr	FLoadMopacZmat (FILE *fp);
extern char	*GetMopacRoute(void);
extern void	SetMopacRoute(char *);
extern void	SPrintMopacZmat(char **, MopacZmatPtr);
extern void	FPrintMopacZmat(FILE *, MopacZmatPtr);
#endif
#ifndef __BOSS_H_
#define __BOSS_H_

/*************** BOSS Calculation Types **********************/
enum	{BOSS_CALC_OPT, BOSS_CALC_OPT_NMODE, BOSS_CALC_NMODE,
	BOSS_CALC_SP, BOSS_CALC_DIHED, BOSS_CALC_CSEARCH,
	BOSS_CALC_MC_GAS, BOSS_CALC_MC_LIQ, BOSS_CALC_MC_LIQ_CLUSTER,
	BOSS_CALC_MC_SOL, BOSS_CALC_MC_SOL_CLUSTER, BOSS_CALC_TYPE_MAX};

#define IS_CALC_OPT(t)		(t==BOSS_CALC_OPT||t==BOSS_CALC_OPT_NMODE||t==BOSS_CALC_SP||t==BOSS_CALC_DIHED)
#define IS_CALC_CSEARCH(t)	(t==BOSS_CALC_CSEARCH)
#define IS_CALC_FREQ(t)		(t==BOSS_CALC_NMODE||t==BOSS_CALC_OPT_NMODE)
#define IS_CALC_MC_GAS(t)	(t==BOSS_CALC_MC_GAS)
#define IS_CALC_MC_LIQ(t)	(t==BOSS_CALC_MC_LIQ||t==BOSS_CALC_MC_LIQ_CLUSTER)
#define IS_CALC_MC_SOL(t)	(t==BOSS_CALC_MC_SOL||t==BOSS_CALC_MC_SOL_CLUSTER)
#define IS_CALC_MC_GAS_SOL(t)	(IS_CALC_MC_GAS(t)||IS_CALC_MC_SOL(t))
#define IS_CALC_MC_LIQ_SOL(t)	(IS_CALC_MC_LIQ(t)||IS_CALC_MC_SOL(t))
#define IS_CALC_MC(t)		(IS_CALC_MC_GAS(t)||IS_CALC_MC_LIQ(t)||IS_CALC_MC_SOL(t))
#define IS_CALC_MC_CLUSTER(t)	(t==BOSS_CALC_MC_LIQ_CLUSTER||t==BOSS_CALC_MC_SOL_CLUSTER)

#define NEED_OPTIMIZER(t)	(t==BOSS_CALC_OPT||t==BOSS_CALC_OPT_NMODE||t==BOSS_CALC_DIHED||t==BOSS_CALC_CSEARCH)
#define NEED_SOLVENT(t)		(IS_CALC_MC_LIQ_SOL(t))

/**************** BOSS/MCPRO Z-matrix components ****************/
enum    {BOSS_GEOMVAR, BOSS_VARBOND, BOSS_ADDBOND, BOSS_HARMONIC,
	BOSS_VARANGLE, BOSS_ADDANGLE, BOSS_VARDIHED, BOSS_ADDDIHED,
	BOSS_DOMAIN, BOSS_CONFORM, BOSS_LOCALHEAT, BOSS_EXCLUDE};

/**************** BOSS versions ****************/
enum	{BOSS_V33=33, BOSS_V34, BOSS_V35, BOSS_V36, BOSS_V37, BOSS_V38, BOSS_V39, BOSS_V40,
	BOSS_V41, BOSS_V42, BOSS_V43,
	BOSS_VMAX};

enum	{
	P_SVMOD1=1,
	P_SLFMT,

	P_QMNAME, P_CMNAME, P_QMSCAL, P_ICH0, P_ICH1, P_ICH2,	/* BOSS4.0 QM */

	P_SOLOR,
	P_ICAPAT, P_CAPRAD, P_CAPFK,
	P_OPTIM,  P_FTOL,   P_NEWZM,  P_NSAEVL, P_NSATR,  P_SATEMP, P_SART,
	P_ICSOPT, P_CSETOL, P_CSRMS,  P_CSRTOL, P_CSREUP, P_CSAELO, P_CSAEUP,
	P_NMODE,  P_NSYM,   P_NCONF,  P_FSCALE, P_FRQCUT,
	P_NMOL,   P_IBOX,   P_BCUT,
	P_SVMOD2, P_NMOL2,
	P_NCENT1, P_NCENT2, P_NCENTS,
	P_NROTA1, P_NROTA2,
	P_ICUT,   P_ICUTAT,
	P_IRECNT, P_INDSOL, P_IZLONG, P_MAXOVL, P_NOXX, P_NOSS, P_NOBNDV, P_NOANGV,
	P_NRDF,   P_NRDFS,
	P_NRDFS1, P_NRDFS2,
	P_NRDFA1, P_NRDFA2,
	P_RDLMIN, P_RDLINC,
	P_EPRMIN, P_EPRINC,
	P_EDMIN,  P_EDINC,
	P_ESSMIN, P_ESSINC,
	P_EBSMIN, P_EBSINC,
	P_NVCHG,  P_NSCHG,  P_MAXVAR,
	P_NCONSV, P_NBUSE,  P_NOCUT,  P_NOSMTH,
	P_VDEL,   P_WKC,
	P_RDEL,   P_ADEL,   P_RADSOL,
	P_RDELS1, P_ADELS1,
	P_RDELS2, P_ADELS2,
	P_LHTSOL, P_TLHT,
	P_RCUT,   P_SCUT,   P_CUTNB,
	P_T,      P_P,      P_TORCUT,
	P_DIELEC, P_SCL14C, P_SCL14L, P_DIELRF,
	P_PLTFMT, P_ISOLEC,
	P_NOBACK, P_IBDPDB, P_NSCHG2,	/* MCPRO-only parameters */
	P_MAXPLUS1};

typedef struct {
	char	*name;		/* parameter name */
	char	*defval;	/* default parameter value */
	char	*val;		/* parameter value */
	char	*scanfmt;	/* scan format */
	char	*printfmt;	/* print format */
	int	size;		/* char val[size] */
	int	id;		/* parameter ID */
	int	version;	/* version-specific parameters */
	char	*comment;	/* comment */
	} BossParEntry;

#define BP_MCPRO	0x1	/* same as B_MCPRO */
typedef struct {
	int	flags;
	int	version;
	int	maxentry;	/* MAXBOSSPARENTRY or MAXMCPROPARENTRY */
	BossParEntry	*par;
	} BossPar, *BossParPtr;

/* flags for BossZmatPtr->flags */
#define B_MCPRO	BP_MCPRO	/* 0x1 */
#define B_QM			0x2
#define B_SOLVENT_ZMAT		0x4
/* 0x8 */
/* 0x10 */
/* 0x20 */
/* 0x40 */
/* 0x80 */
#define B_AUTO_ADDBOND		0x8000
#define B_AUTO_ADDANGLE		0x10000
#define B_AUTO_ADDDIHED		0x20000
#define B_AUTO_CSEARCH		0x40000
#define B_VARDIHED_EXPLICIT	0x80000

#define B_NO_REWIND_PARFILE	0x100000	/* temporary flag */
#define B_PARAM_IN_ZMAT		0x200000	/* temporary flag */

#define B_USER_CHARGE		0x400000
#define B_USER_SIGMA		0x800000
#define B_USER_EPSILON		0x1000000

#define B_AUTO_ADD	(B_AUTO_ADDBOND|B_AUTO_ADDANGLE|B_AUTO_ADDDIHED)

typedef struct {
	unsigned long	flags;
	char	title[256];
	int	version;

	ListPtr	zmatlist;	/* Z-matrix */

	ListPtr	atomtypes;	/* atom types */
	ListPtr	torsiontypes;	/* torsion types */

	ListPtr	geomvar;	/* geometry variations */
	ListPtr	varbond;	/* variable bonds */
	ListPtr	addbond;	/* additional bonds */
	ListPtr	harmonic;	/* harmonic constraints */
	ListPtr	varangle;	/* variable angles */
	ListPtr	addangle;	/* additional angles */
	ListPtr	vardihed;	/* variable dihedrals */
	ListPtr	adddihed;	/* additional dihedrals */
	ListPtr	domain;		/* domain definitions */
	ListPtr	conform;	/* internal coords to sample in conformational search (BOSS only) */
	ListPtr	exclude;	/* excluded atom list (MCPRO only) */
	ListPtr	localheat;	/* local heating (BOSS only) */

	BossPar	*param;	/* list of parameters */
	} BossZmat, *BossZmatPtr;

/************************************************************
 * Data structure for atomtype and torsion parameters.
 * These parameters are usually read from a BOSS/MCPRO
 * parameter file or "oplsaa.par".
 ************************************************************/

#define OPLS_TYPE_ATOM		0x1
#define OPLS_TYPE_TORSION	0x2
#define OPLS_TYPE_STRETCH	0x4
#define OPLS_TYPE_BENDING	0x8

typedef struct {
	ListPtr	atomtypelist;
	ListPtr	torsionlist;
	} OplsTypeTorsion;

/* data structure for stretching and bending parameters.
 * These are read from "oplsaa.sb"
 */
typedef struct {
	ListPtr	stretchlist;
	ListPtr	bendinglist;
	} OplsStretchBend;

/************************************************
 ***   BOSS ADD ALL Z-MATRIX VARIABLES
 ************************************************/
/* masks for add all */
#define M_BOSS_ZMAT_VAR_BOND		0x1
#define M_BOSS_ZMAT_VAR_ANGLE		0x2
#define M_BOSS_ZMAT_VAR_DIHED		0x4
#define M_BOSS_ZMAT_ADD_BOND		0x8
#define M_BOSS_ZMAT_ADD_ANGLE		0x10
#define M_BOSS_ZMAT_ADD_DIHED		0x20

#define M_BOSS_ZMAT_AUTO_ADD_BOND	0x40
#define M_BOSS_ZMAT_AUTO_ADD_ANGLE	0x80
#define M_BOSS_ZMAT_AUTO_ADD_DIHED	0x100
#define M_BOSS_ZMAT_AUTO_CSEARCH	0x200

#define M_BOSS_ZMAT_CSEARCH		0x400

#define M_BOSS_ZMAT_VAR_DIHED_EXPLICIT	0x800

#define M_BOSS_ZMAT_VAR		(M_BOSS_ZMAT_VAR_BOND|M_BOSS_ZMAT_VAR_ANGLE|M_BOSS_ZMAT_VAR_DIHED)
#define M_BOSS_ZMAT_ADD		(M_BOSS_ZMAT_ADD_BOND|M_BOSS_ZMAT_ADD_ANGLE|M_BOSS_ZMAT_ADD_DIHED)
#define M_BOSS_ZMAT_VARADD	(M_BOSS_ZMAT_VAR|M_BOSS_ZMAT_ADD)

#define M_BOSS_ZMAT_AUTO_ADD	(M_BOSS_ZMAT_AUTO_ADD_BOND|M_BOSS_ZMAT_AUTO_ADD_ANGLE|M_BOSS_ZMAT_AUTO_ADD_DIHED)
#define M_BOSS_ZMAT_AUTO	(M_BOSS_ZMAT_AUTO_ADD|M_BOSS_ZMAT_AUTO_CSEARCH)

#define M_BOSS_ZMAT_DEFAULT_BOND_ANGLE	(M_BOSS_ZMAT_VAR_BOND|M_BOSS_ZMAT_ADD_BOND|M_BOSS_ZMAT_VAR_ANGLE|\
				M_BOSS_ZMAT_AUTO_ADD_ANGLE)
/* implicit variable and additional dihedrals */
#define M_BOSS_ZMAT_DEFAULT	(M_BOSS_ZMAT_DEFAULT_BOND_ANGLE|\
				M_BOSS_ZMAT_VAR_DIHED|M_BOSS_ZMAT_AUTO_ADD_DIHED)
/* explicit variable dihedrals */
#define M_BOSS_ZMAT_DEFAULT2	(M_BOSS_ZMAT_DEFAULT_BOND_ANGLE|\
				M_BOSS_ZMAT_VAR_DIHED|M_BOSS_ZMAT_AUTO_ADD_DIHED|M_BOSS_ZMAT_VAR_DIHED_EXPLICIT)
/* explicit additional dihedrals */
#define M_BOSS_ZMAT_DEFAULT3	(M_BOSS_ZMAT_DEFAULT_BOND_ANGLE|\
				M_BOSS_ZMAT_VAR_DIHED|M_BOSS_ZMAT_ADD_DIHED)
/* explicit variable and additional dihedrals */
#define M_BOSS_ZMAT_DEFAULT4	(M_BOSS_ZMAT_DEFAULT_BOND_ANGLE|\
				M_BOSS_ZMAT_VAR_DIHED|M_BOSS_ZMAT_ADD_DIHED|M_BOSS_ZMAT_VAR_DIHED_EXPLICIT)

/* options for solvent Z-matrix */
#define M_BOSS_ZMAT_SOLVENT	(M_BOSS_ZMAT_VARADD|M_BOSS_ZMAT_VAR_DIHED_EXPLICIT)

typedef struct {
	int	var_dihed_itype;
	int	var_dihed_ftype;
	float	var_dihed_range;

	float	var_bond_rmin;
	float	var_bond_rmax;
	float	var_angle_rmin;
	float	var_angle_rmax;
	float	var_dihed_rmin;
	float	var_dihed_rmax;
	} BossZmatAddAllOpt;

#endif

#ifndef __BOSSINPUT_H_
#define __BOSSINPUT_H_

extern ListPtr	FLoadBossList (FILE *fp, char *scanfmt[], int nvar);
extern void	FreeBossZmat (BossZmatPtr *bp);
extern void	FPrintBossZmat (FILE *fp, BossZmatPtr boss);
extern void	SPrintBossZmat (char **str_sum, BossZmatPtr boss);
extern void	FPrintBossZmatTypes (FILE *fp, BossZmatPtr boss);
extern void	SPrintBossZmatTypes (char **str_sum, BossZmatPtr boss);
extern void	FPrintBossZmatWithTypes (FILE *fp, BossZmatPtr bos);
extern void	SPrintBossZmatWithTypes (char **str_sum, BossZmatPtr bos);
extern void	FixBossZmatAtomType (BossZmatPtr boss);
extern void	UpdateBossZmatTypes (BossZmatPtr boss, OplsTypeTorsion *opls);

extern BossZmatPtr	FLoadBossZmatInput (FILE *fp, int flags);
extern BossZmatPtr	FLoadBossZmat (FILE *fp, FILE *parfp, int flags);
extern BossZmatPtr	CopyBossZmat (BossZmatPtr oldboss);
extern BossZmatPtr	CreateBossZmatFromMol (MolPtr mol, int flags, int autozmat_flags);
extern BossZmatPtr	CreateBossZmatFromZmat (ZmatPtr zmatlist, int flags);
extern BossZmatPtr	CreateBossZmatFromGauss (GaussZmatPtr gauss);
extern BossZmatPtr	CreateBossZmatFromJaguar (JaguarZmatPtr jaguar);

extern MolPtr	CreateMolFromBossZmat (BossZmatPtr boss);
extern MolPtr	FLoadBossPdbInput (FILE *fp);
extern MolPtr	FLoadBossMindtoolInput (FILE *fp);

extern int	UpdateBossGeomVarWithMolCoord (BossZmatPtr boss, MolPtr mol);
extern ZmatPtr	GetZmatFromBossGeomVar (BossZmatPtr boss);
extern ListPtr	GetZmatListFromBossGeomVar (BossZmatPtr boss);
extern int	EvalBossGeomVar (BossZmatPtr boss, MolPtr mol);

extern GaussZmatPtr	CreateGaussZmatFromBoss (BossZmatPtr boss);
extern JaguarZmatPtr	CreateJaguarZmatFromBoss (BossZmatPtr boss);
extern void	SetImproperTorsionFlags (ListPtr dihedlist, ChainPtr chainlist, BondPtr bondlist);

extern int	AddAllBossZmatOptions (BossZmatPtr boss, MolPtr mol, unsigned long flags, BossZmatAddAllOpt *addall);
extern BossZmatAddAllOpt	*GetAddAllBossZmatOptions (void);

#endif
#ifndef __BOSSPAR_H_
#define __BOSSPAR_H_
extern BossPar	*CreateBossPar (int flags);
extern BossPar	*CopyBossPar(BossPar *);
extern void	FreeBossPar (BossPar **bp);
extern BossParEntry	*FindBossParEntry (BossPar *bp, int param_id);
extern char	*GetBossParDefaultValue (BossPar *bp, int param_id);
extern char	*GetBossParValue (BossPar *bp, int param_id);
extern void	SetBossParValue (BossPar *bp, int param_id, char *val);

extern int	FGuessBossParVersion (FILE *fp);
extern BossPar	*FLoadBossPar (FILE *fp, int version, int flags);
extern void	FPrintBossPar (FILE *fp, BossPar *bp, int version);
extern ListPtr	FLoadBossAtomTypes (FILE *fp);
extern ListPtr	FLoadBossAtomTypesFromZmat (FILE *fp);
extern void	FPrintBossOneAtomType (FILE *fp, ListPtr atomtype);
extern void	FPrintBossAtomTypes (FILE *fp, ListPtr atomtypes, int version);
extern ListPtr	FLoadBossTorsion (FILE *fp);
extern void	FPrintBossOneTorsionType (FILE *fp, ListPtr torsiontype);
extern void	FPrintBossTorsionTypes (FILE *fp, ListPtr torsiontypes, int version);
extern ListPtr	FindBossAtomType (ListPtr atomtypelist, int type);
extern void	AssignBossPar (ListPtr zmatlist, ListPtr atomtypelist);

extern void	FPrintBossParFile (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes, int version);
extern void	FPrintBossParFileWithZmatTypes (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes,
	ListPtr zmat_atomtypes, ListPtr zmat_torsiontypes, int version);

extern int	FLoadBossParFile (FILE *fp, int flags, BossParPtr *bp_ret, ListPtr *at_ret, ListPtr *tt_ret, int *version_ret);

extern OplsTypeTorsion	*FLoadOplsTypeTorsion (FILE *fp);
extern OplsTypeTorsion	*FLoadDefaultOplsTypeTorsion (FILE *fp);
extern OplsTypeTorsion	*LoadOplsTypeTorsion (char *filename);
extern OplsTypeTorsion	*LoadDefaultOplsTypeTorsion (char *filename);
extern OplsTypeTorsion	*GetDefaultOplsTypeTorsion (void);
extern void	FreeOplsTypeTorsion (OplsTypeTorsion **p);

extern ListPtr	FLoadOplsStretch (FILE *fp);
extern ListPtr	FLoadOplsBending (FILE *fp, int rewind_flag);
extern OplsStretchBend	*FLoadOplsStretchBend (FILE *fp);
extern OplsStretchBend	*FLoadDefaultOplsStretchBend (FILE *fp);
extern OplsStretchBend	*LoadOplsStretchBend (char *filename);
extern OplsStretchBend	*LoadDefaultOplsStretchBend (char *filename);
extern OplsStretchBend	*GetDefaultOplsStretchBend (void);
extern void	FreeOplsStretchBend (OplsStretchBend **p);

extern void	AssignDefaultBossParByCalcType (BossParPtr bp, int calc_type);
extern char	*GetBossCalcTypeText (int calc_type);

#endif
#ifndef __MMODEL_H_
#define __MMODEL_H_

#define MMOD_MAX_BONDS	6

extern int	MModelToAn (int type);
extern int	GetMModelAtomType (AtomPtr atom);
extern int	GetMModelBondType (int bondtype);
extern MolPtr	FLoadFullMModel (FILE *fp, int natoms, char *title);
extern MolPtr	FLoadMModel (FILE *fp);
extern void	FPrintMModel (FILE *fp, MolPtr m);
#endif
#ifndef __CARTESIAN_H_
#define __CARTESIAN_H_

#define MINDTOOL_BREAK_PERCENT	0x1

extern MolPtr	FLoadZmatCartesian (FILE *fp, ZmatPtr *zmatlist_ret);
extern MolPtr	FLoadCartesian (FILE *fp);
extern MolPtr	FLoadMindtool (FILE *fp, int flags);
extern void	FPrintCartesian (FILE *fp, ChainPtr chainlist);
extern void	FPrintMindtool (FILE *fp, ChainPtr chainlist);

#endif
#ifndef __UTILS_H_
#define __UTILS_H_

extern void	TrimHead (char *str, char *delim);
extern void	TrimTail (char *str, char *delim);
extern void	TrimStr (char *str, char *delim);
extern void	AddtoString (char **str_sum, char *str);
extern int	IsBlank (char *str);
extern void	StringUpperCase (char *str);
extern void	StringLowerCase (char *str);

extern void	FPrintString (void *dest, char *str);
extern void	SPrintString (void *dest, char *str);
extern void	BgnPrintString (void);
extern void	ContinuePrintString (char *str);
extern char	*EndPrintString (void);

#endif
#ifndef __HYDROGEN_H_
#define __HYDROGEN_H_

enum	{HYBRID_UNK, HYBRID_SP, HYBRID_SP2, HYBRID_SP3, HYBRID_SP3_PLUS};

enum    {
	H_MAX_SUM_ANG, H_MAX_C1_ANG, H_MAX_C2_ANG,
	H_MAX_CC2_LEN, H_MAX_CC3_LEN, H_MIN_AMIDE_ANG, H_MAX_AMIDE_LEN,
	H_MAX_CN2_LEN, H_MAX_CN3_LEN, H_MAX_CO2_LEN, H_MAX_CS2_LEN,
	H_PREF_MAX};

#define HYDROGEN_PREF_TEXT	{\
	"Max sum of three Csp3 angles",\
	"Max Csp3 (X-C-Y) angle",\
	"Max Csp2 (X=C-Y) angle",\
	"Max C double bond length",\
	"Max C triple bond length",\
	"Min amide angle",\
	"Max amide bond length",\
	"Max C=N bond length",\
	"Max C-N triple bond length",\
	"Max C=O bond length",\
	"Max C=S bond length"}

#define HYDROGEN_HELP_TEXT	{\
	"Csp3: max sum of three Csp3 angles (\260)",\
	"Csp3: max Csp3 (X-C-Y) angle (\260)",\
	"Csp2: max Csp2 (X=C-Y) angle (\260)",\
	"Csp2: max C double bond length (\305)",\
	"Csp: max C triple bond length (\305)",\
	"Namide: min amide angle (\260)",\
	"Namide: max amide bond length (\260)",\
	"Nsp2: max C-N double bond length (\305)",\
	"Nsp: max C-N triple bond length (\305)",\
	"Osp2: max C=O bond length (\305)",\
	"Ssp2: max C=S bond length (\305)"}

extern double	MAX_SUM_ANG;
extern double	MAX_C1_ANG;
extern double	MAX_C2_ANG;
extern double	MAX_CC2_LEN;
extern double	MAX_CC3_LEN;
extern double	MIN_AMIDE_ANG;
extern double	MAX_CN2_LEN;
extern double	MAX_AMIDE_LEN;
extern double	MAX_CN3_LEN;
extern double	MAX_CO2_LEN;
extern double	MAX_CS2_LEN;

extern int	AddHydrogenToMol (MolPtr mol);
extern AtomPtr	AddHydrogenToAtom (AtomPtr theatom);
extern int	AddHydrogenToMolByNumber (MolPtr mol, int nC, int nN, int nO, int nS);
extern AtomPtr	AddHydrogenToAtomByNumber (AtomPtr theatom, int nhydrogens);
extern void	SetAddHydrogenPref (int id, float val);
extern float	GetAddHydrogenPref (int id, float *val);

extern int	GetHybridization (AtomPtr atom);

#endif
#ifndef __MDL_H_
#define __MDL_H_

extern MolPtr	FLoadMdlMol (FILE *);
extern MolPtr	FLoadMdlSd (FILE *);
extern void	FPrintMdlMol (FILE *, MolPtr);
extern void	FPrintMdlSd (FILE *, MolPtr);

#endif
#ifndef __PEPZ_H_
#define __PEPZ_H_

/* flags for FPrintPepz() */
#define PEPZ_ATOMTYPE_FROM_ZMAT		0x1

extern void	FPrintPepz (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist, BossZmatPtr boss, ListPtr atomtypelist, int flags);

#endif
#ifndef __ANALGEOM_H_
#define __ANALGEOM_H_
extern double	GetAtomDistance (AtomPtr atom1, AtomPtr atom2);
extern double	GetAtomSQDist (AtomPtr atom1, AtomPtr atom2);
extern double	GetAngle (double *v1, double *v2, double *v3);
extern double	GetAtomAngle (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3);
extern double	GetDihedral (double *v1, double *v2, double *v3, double *v4);
extern double	GetAtomDihedral (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, AtomPtr atom4);
extern void	MakeAtomVectors (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *v1, float *v2);
extern int	GetAtomNormal (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *norm);
extern int	IsAtomLinear (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomCis (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomTrans (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomPerp (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsLinear (double *v1, double *v2, double *v3, float epsilon);
#endif
#ifndef __TRIPOS_H_
#define __TRIPOS_H_

extern MolPtr	FLoadTriposMol (FILE *fp);
extern MolPtr	FLoadTriposMol2 (FILE *fp);
extern void	FPrintTriposMol (FILE *fp, MolPtr m);
extern void	FPrintTriposMol2 (FILE *fp, MolPtr m);

#endif
#ifndef __CONNECT_H_
#define __CONNECT_H_

extern ListPtr	LookUpBondList (ListPtr bondlist, int a1, int a2);
extern ListPtr	LookUpAngleList (ListPtr anglelist, int a1, int a2, int a3);
extern ListPtr	LookUpDihedList (ListPtr dihedlist, int a1, int a2, int a3, int a4);
extern ListPtr	LookUpImproperList (ListPtr improperlist, ListPtr dihedlist, int a1, int a2, int a3, int a4);
extern ListPtr	MakeBondList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeAngleList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeDihedList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeImproperList (ChainPtr clist, ZmatPtr zmatlist, ListPtr dihedlist);
extern int	MakeConnect (ChainPtr clist, ZmatPtr zmatlist, ListPtr *listbnd, ListPtr *listang, ListPtr *listdih, ListPtr *listimp);
extern ListPtr	GetAddBondList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	GetAddAngleList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	GetAddDihedList (ChainPtr clist, ZmatPtr zmatlist);

#endif
#endif /* _MOLECULE_H_ */
