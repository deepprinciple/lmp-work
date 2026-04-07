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

