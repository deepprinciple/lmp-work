#ifndef __AUTOZMAT_H__
#define __AUTOZMAT_H__

#ifndef PATH_MAX
#define PATH_MAX	1024
#endif

/* EGCS-MINGW32 */
#if defined(__i386__) && !defined(linux)
#define MSDOS
#endif

#ifndef MSDOS
#define DIRSEP	'/'
#else
#define DIRSEP	'\\'
#endif

/* macros */
#define BZERO(s, n)		memset((char *)s, (int)'\0', n)
#define BCOPY(s1, s2, n)	memmove((char *)s2, (char *)s1, n)
#define ABS(x)			(((x) > 0) ? (x) : -(x))

/* masks for autozmat_flags */
#define M_VERBOSE		0x1	/* set if -V option is specified */
#define M_VERBOSE_ATOMTYPE	0x2	/* causes atomtyping to be verbose */
#define M_REMOVE_DUMMY		0x4	/* set if -x option is specified */
#define M_NO_BGN_DUMMY		0x8	/* When creating BOSS z-matrix, don't add dummies at the beginning */
#define M_NO_ATOM_TYPING	0x10	/* instructs not to perform atom typing */
#define M_NO_TORSION_TYPING	0x20	/* instructs not to perform torsion typing */
#define M_ABORT_ON_ATOMTYPE_ERR	0x40	/* abort if error occurs during atomtyping */
#define M_BREAK_MIND		0x80	/* break atoms delimited by '%' into separate molecules */
#define M_USER_CHARGE		0x100	/* in OPLSAA atomtyping, replace charges by user-supplied ones */
#define M_USER_SIGMA		0x200	/* in OPLSAA atomtyping, replace sigmas by user-supplied ones */
#define M_USER_EPSILON		0x400	/* in OPLSAA atomtyping, replace epsilons by user-supplied ones */
#define M_REPORT_MW		0x800	/* reports the molecular weight */
#define M_QUIET			0x1000	/* be quiet */

#define M_ASSIGN_XYZ		0x10000	/* temporary use only (no equivalent cmd-line option) */

/*---------- convenience flags ------------*/
#define M_NO_TYPING			(M_NO_ATOM_TYPING|M_NO_TORSION_TYPING)
#define DO_ATOM_TYPING(flags)		!(flags & M_NO_ATOM_TYPING)
#define DO_TORSION_TYPING(flags)	!(flags & M_NO_TORSION_TYPING)
#define DO_TYPING(flags)		(DO_ATOM_TYPING(flags)||DO_TORSION_TYPING(flags))
#define M_REPORT			(M_REPORT_MW)

/* file format defines */
#define NEED_PARAM_FILE(f)		(f==FILE_BOSS_ZMAT||f==FILE_MCPRO_ZMAT||f==FILE_BOSS_SOLVENT_ZMAT)

#define IS_BOSSZMAT_FORMAT(f)		(f==FILE_BOSS_ZMAT||f==FILE_MCPRO_ZMAT||f==FILE_BOSS_QM_ZMAT||\
					f==FILE_MCPRO_QM_ZMAT||f==FILE_BOSS_SOLVENT_ZMAT)

#define IS_QM_BOSSZMAT_FORMAT(f)	(f==FILE_BOSS_QM_ZMAT||f==FILE_MCPRO_QM_ZMAT)
#define IS_ZMAT_FORMAT(f)		(IS_BOSSZMAT_FORMAT(f)||f==FILE_GAUSS_ZMAT||\
					f==FILE_MOPAC_ZMAT||f==FILE_JAGUAR_ZMAT)
#define NEED_CONTAB(f)			(f==FILE_MMODEL||f==FILE_MDL||f==FILE_TRIPOS_MOL||\
					f==FILE_TRIPOS_MOL2||f==FILE_PDB)

/* masks */
#define M_ADD_H_TO_C		0x1
#define M_ADD_H_TO_N		0x2
#define M_ADD_H_TO_O		0x4
#define M_ADD_H_TO_S		0x8

#define M_ADD_H		(M_ADD_H_TO_C|M_ADD_H_TO_N|M_ADD_H_TO_O|M_ADD_H_TO_S)

/* file formats */
enum	{
	FILE_BOSS_QM_ZMAT, FILE_BOSS_ZMAT, FILE_BOSS_SOLVENT_ZMAT,
	FILE_MCPRO_QM_ZMAT, FILE_MCPRO_ZMAT,
	FILE_GAUSS_ZMAT, FILE_JAGUAR_ZMAT, FILE_MOPAC_ZMAT,
	FILE_CARTESIAN, FILE_PDB,
	FILE_MMODEL, FILE_MDL, FILE_MINDTOOL, FILE_PEPZ,
	FILE_TRIPOS_MOL, FILE_TRIPOS_MOL2,
	FILE_MAX};

typedef struct {
	char	*text;
	char	*desc;
	} FileFormat;

#define FILE_FORMAT	{\
	{"qm_boss",  "BOSS Z-matrix for QM/MM calcs"},\
	{"boss",     "BOSS Z-matrix"},\
	{"solv_boss","BOSS Z-matrix for solvent"},\
	{"qm_mcpro", "MCPRO Z-matrix for QM/MM calcs"},\
	{"mcpro",    "MCPRO Z-matrix"},\
	{"gauss",    "GAUSSIAN Z-matrix"},\
	{"jaguar",   "JAGUAR Z-matrix"},\
	{"mopac",    "MOPAC Z-matrix"},\
	{"xyz",      "Cartesian coordinates"},\
	{"pdb",      "Protein DataBase format"},\
	{"mmodel",   "MacroModel"},\
	{"mdl",      "MDL MOL or SD format"},\
	{"mind",     "Mindtool (WL Jorgensen)"},\
	{"pepz",     "PEPZ (WL Jorgensen)"},\
	{"mol",      "Tripos Sybyl MOL"},\
	{"mol2",     "Tripos Sybyl MOL2"}\
	}

/* default input, output, and template Z-matrix file formats */
#define DEFAULT_IFMT	FILE_PDB
#define DEFAULT_OFMT	FILE_BOSS_QM_ZMAT
#define DEFAULT_TFMT	FILE_GAUSS_ZMAT

/******** -z option **********/
typedef struct {
	int	flags;	/* flag(s) */
	char	*text;	/* text in the command line argument */
	char	*desc;	/* description */
	} ZOption;

#define ZOPTION_TAB	{\
	{M_BOSS_ZMAT_VAR_BOND,  "varbond",  "inserts variable bonds"},\
	{M_BOSS_ZMAT_ADD_BOND,  "addbond",  "inserts additional bonds"},\
	{M_BOSS_ZMAT_VAR_ANGLE, "varangle", "inserts variable angles"},\
	{M_BOSS_ZMAT_ADD_ANGLE, "addangle", "inserts additional angles"},\
	{M_BOSS_ZMAT_VAR_DIHED, "vardihed", "inserts variable dihedrals"},\
	{M_BOSS_ZMAT_VAR_DIHED|M_BOSS_ZMAT_VAR_DIHED_EXPLICIT,\
		"vardihed_explicit", "inserts variable dihedrals with explicit torsion types"},\
	{M_BOSS_ZMAT_ADD_DIHED, "adddihed", "inserts additional dihedrals"},\
	{M_BOSS_ZMAT_VARADD,    "varadd",   "inserts all of above"},\
\
	{M_BOSS_ZMAT_AUTO_ADD_BOND,  "auto_addbond",  "inserts AUTO for additional bonds"},\
	{M_BOSS_ZMAT_AUTO_ADD_ANGLE, "auto_addangle", "inserts AUTO for additional angles"},\
	{M_BOSS_ZMAT_AUTO_ADD_DIHED, "auto_adddihed", "inserts AUTO for additional dihedrals"},\
	{M_BOSS_ZMAT_AUTO_ADD,       "auto_add",      "inserts AUTO for additional bonds, angles, and dihedrals"},\
\
	{M_BOSS_ZMAT_AUTO_CSEARCH,   "auto_search",  "inserts AUTO for conformational searching"},\
	{M_BOSS_ZMAT_CSEARCH, "search",  "inserts conformational search parameters"},\
\
	{M_BOSS_ZMAT_DEFAULT, "default",\
		"equivalent to 'varbond,addbond,varangle,auto_addangle,vardihed,auto_adddihed'"},\
	{M_BOSS_ZMAT_DEFAULT2, "default2",\
		"equivalent to 'varbond,addbond,varangle,auto_addangle,vardihed_explicit,auto_adddihed'"},\
	{M_BOSS_ZMAT_DEFAULT3, "default3",\
		"equivalent to 'varbond,addbond,varangle,auto_addangle,vardihed,adddihed'"},\
	{M_BOSS_ZMAT_DEFAULT4, "default4",\
		"equivalent to 'varbond,addbond,varangle,auto_addangle,vardihed_explicit,adddihed'"},\
\
	{M_BOSS_ZMAT_DEFAULT2, "default_var", "equivalent to 'default2'"},\
	{M_BOSS_ZMAT_DEFAULT3, "default_add", "equivalent to 'default3'"},\
	{M_BOSS_ZMAT_DEFAULT4, "default_varadd", "equivalent to 'default4'"},\
	{M_BOSS_ZMAT_SOLVENT, "solvent", "'equivalent to varbond,addbond,varangle,addangle,vardihed_explicit,adddihed'"}\
	}

/* general options (-g option) */
#define OPTION_TAB	{\
	{M_VERBOSE,              "verbose",              "prints verbose info"},\
	{M_VERBOSE_ATOMTYPE,     "verbose_atomtype",     "prints a more verbose information in atomtyping"},\
	{M_REMOVE_DUMMY,         "remove_dummy",         "removes dummy atoms before printing"},\
	{M_NO_BGN_DUMMY,         "no_bgn_dummy",         "don't add dummies at the beginning in BOSS Z-matrix."},\
	{M_NO_ATOM_TYPING,       "no_atom_typing",       "instructs not to assign atom types"},\
	{M_NO_TORSION_TYPING,    "no_torsion_typing",    "instructs not to assign torsion types"},\
	{M_ABORT_ON_ATOMTYPE_ERR,"abort_on_atomtype_err","aborts if error occurs during atomtyping"},\
	{M_BREAK_MIND,           "break_mind",           "breaks atoms delimited by '%' into separate molecules"},\
	{M_USER_CHARGE,          "user_charge",          "in OPLSAA atomtyping, use user-supplied charges"},\
	{M_USER_SIGMA,           "user_sigma",           "in OPLSAA atomtyping, use user-supplied sigmas"},\
	{M_USER_EPSILON,         "user_epsilon",         "in OPLSAA atomtyping, use user-supplied epsilons"},\
	{M_REPORT_MW,            "report_mw",            "reports the molecular weight"},\
	{M_QUIET,                "quiet",                "be quiet!"},\
	{0,                      "pepz=",                "sets the starting atom type in PEPZ, .e.g, '-g pepz=500'"}\
	}

/* external variables */
extern int	autozmat_flags;
extern char	default_opls_par_file[];

/* functions */
extern void	init_mol_info(MolPtr, char *, int, int);

#endif	/* __AUTOZMAT_H__ */
