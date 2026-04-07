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
