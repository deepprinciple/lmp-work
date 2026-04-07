#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <ctype.h>

#include "macros.h"
#include "molecule.h"

/*++++++++++++++++++++++++++++++++++++++++++++*/
static OplsTypeTorsion	*opls_typetorsion=NULL;
static OplsStretchBend	*opls_stretchbend=NULL;
/*++++++++++++++++++++++++++++++++++++++++++++*/

static char	_f[8][8];
static char	*scanfmt[8] = {_f[0],_f[1],_f[2],_f[3],_f[4],_f[5],_f[6],_f[7]};

/* Note
 * BossParEntry->version = 0 by default
 *    this is marked by, e.g, BOSS_V38, if it is introduced at version 3.8
 */

/* this is template BOSS parameter entries */
static BossParEntry	bossparentry[] = {
	{0, 0, 0, 0, 0, 0, 0, 0,
		"SET-UP AND POTENTIAL FUNCTION PARAMETERS:"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"SVMOD1 (SOLVENT MODEL) - TIP3P, TIP4P, CH3OH, DMSO"},
	{"SVMOD1", "TIP4P",  NULL, "%5c", "%-5s", 5, P_SVMOD1, 0,
		"  (1A5)    THF, CCL4, ARGON, CH3CN, MEOME, C3H8, MECL2, CHCL3, OTHER or NONE"},

/* BOSS 4.0 QM variables */
	{0, 0, 0, 0, 0, 0, 0, BOSS_V40,
		"QM Parameters (2A4,F5.2,3I3)"},
	{"QMNAME", "NONE", NULL, "%4c", "%4s",   4, P_QMNAME, BOSS_V40, 0},
	{"CMNAME", "CM1 ", NULL, "%4c", "%4s",   4, P_CMNAME, BOSS_V40, 0},
	{"QMSCAL", "1.20", NULL, "%5c", "%5.2f", 5, P_QMSCAL, BOSS_V40, 0},
	{"ICH0  ", " 0 ",  NULL, "%3c", "%3d",   3, P_ICH0,   BOSS_V40, 0},
	{"ICH1  ", " 0 ",  NULL, "%3c", "%3d",   3, P_ICH1,   BOSS_V40, 0},
	{"ICH2  ", " 0 ",  NULL, "%3c", "%3d",   3, P_ICH2,   BOSS_V40, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ORIGIN OF INITIAL SOLUTE COORDINATES - FROM Z-MATRIX, PDB FILE, MIND FILE,"},
	{"SLFMT ", "ZMAT", NULL, "%5c", "%-5s", 5, P_SLFMT, 0,
		"    (1A5)    IN FILE, ZMATRIX & IN FILE (ZMAT, PDB, MIND, IN, ZIN)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ORIGIN OF INITIAL SOLVENT COORDINATES - FROM STORED BOXES, IN FILE or NONE"},
	{"SOLOR ", "BOXES",NULL, "%5c", "%-5s", 5, P_SOLOR, 0,
		"    (1A5)                               (BOXES, IN, NONE)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"FOR SOLVENT CAP - ICAPAT (I4), CAP RADIUS ABOUT ICAPAT, RESTORING FORCE CONST (2F10.4)"},
	{"ICAPAT", "0", NULL, "%4c",  "%04d",    4, P_ICAPAT, 0, 0},
	{"CAPRAD", "0", NULL, "%10c", "%10.4f", 10, P_CAPRAD, 0, 0},
	{"CAPFK ", "",  NULL, "%10c", "%10.4f", 10, P_CAPFK, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"OPTIMIZER (A6),TOLERANCE(F10.4),NEWZM(A5),NSAEVL,NSATR(2I8),SATEMP,SART(2F10.4)"},
	{"OPTIM ", "SIMPLX",  NULL, "%6c",  "%-6s",    6, P_OPTIM,  0, 0},
	{"FTOL  ", "0.001",   NULL, "%10c", "%9.6f ", 10, P_FTOL,   0, 0},
	{"NEWZM ", "NEWZM",   NULL, "%7c",  "%-7s",    7, P_NEWZM,  0, 0},
	{"NSAEVL", "1000000", NULL, "%8c",  "%8d",     8, P_NSAEVL, BOSS_V36, 0},
	{"NSATR ", "5",       NULL, "%8c",  "%8d",     8, P_NSATR,  BOSS_V36, 0},
	{"SATEMP", "5.0",     NULL, "%10c", "%10.4f", 10, P_SATEMP, BOSS_V36, 0},
	{"SART  ", "0.5",     NULL, "%10c", "%10.4f", 10, P_SART,   BOSS_V36, 0},

	{0, 0, 0, 0, 0, 0, 0, BOSS_V36,
		"CONFORM SEARCH OPTION(8I1),CSETOL,CSRMS,CSRTOL,CSREUP,CSAELO,CSAEUP(6F10.4)"},
	{"ICSOPT", "00000000",   NULL, "%8c",  "%-8s",    8, P_ICSOPT, BOSS_V36, 0},
	{"CSETOL", "      0.01", NULL, "%10c", "%10.4f", 10, P_CSETOL, BOSS_V36, 0},
	{"CSRMS ", "      0.25", NULL, "%10c", "%10.4f", 10, P_CSRMS,  BOSS_V36, 0},
	{"CSRTOL", "      5.00", NULL, "%10c", "%10.4f", 10, P_CSRTOL, BOSS_V36, 0},
	{"CSREUP", "     100.0", NULL, "%10c", "%10.4f", 10, P_CSREUP, BOSS_V36, 0},
	{"CSAELO", " -1.00E+06", NULL, "%10c", "%10.3E", 10, P_CSAELO, BOSS_V36, 0},
	{"CSAEUP", "  1.00E+06", NULL, "%10c", "%10.3E", 10, P_CSAEUP, BOSS_V36, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"FOR NORMAL MODE CALC (NMODE, NSYM, NCONF, FSCALE, FRQCUT)"},
	{"NMODE ", "0",     NULL, "%3c",  "%3d",      3, P_NMODE,  0, 0},
	{"NSYM  ", "1",     NULL, "%3c",  "%3d",      3, P_NSYM,   0, 0},
	{"NCONF ", "1",     NULL, "%3c",  "%3d",      3, P_NCONF,  0, 0},
	{"FSCALE", "1.0",   NULL, "%10c", "% 10.4f", 10, P_FSCALE, 0, 0},
	{"FRQCUT", "500.0", NULL, "%10c", "% 10.4f", 10, P_FRQCUT, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NMOL, IBOX, BCUT (2I4,F6.2) - # OF SOLVENT MOLECULES (9999 FOR AUTO),"},

	{"NMOL  ", "260", NULL, "%4c", "%4d",   4, P_NMOL, 0, 0},
	{"IBOX  ", "267", NULL, "%4c", "%4d",   4, P_IBOX, 0, 0},
	{"BCUT  ", "0.0", NULL, "%6c", "%6.2f", 6, P_BCUT, 0,
		"    STORED BOX SIZE, & BOX CUTOFF FOR AUTO (DEF = 7.0)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"2ND SOLVENT (SVMOD2, choices as above) & # OF 2ND SOLVENT MOLECULES (NMOL2)"},
	{"SVMOD2", "NONE", NULL, "%5c", "%-5s", 5, P_SVMOD2, 0, 0},
	{"NMOL2 ", "0",    NULL, "%4c", "%4d",  4, P_NMOL2,  0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NCENT1, NCENT2 (THE ATOMS USED TO DEFINE THE CENTER OF THE SOLUTE(S)) (I3)"},
	{"NCENT1", "1", NULL, "%3c", "%3d", 3, P_NCENT1, 0, 0},
	{"NCENT2", "1", NULL, "%3c", "%3d", 3, P_NCENT2, 0, 0},
	{"NCENTS", "0", NULL, "%3c", "%3d", 3, P_NCENTS, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NROTA1, NROTA2  (THE ATOMS SOLUTES 1 & 2 ARE ROTATED ABOUT - NOTE:"},
	{"NROTA1", "1", NULL, "%3c", "%3d", 3, P_NROTA1, 0, 0},
	{"NROTA2", "1", NULL, "%3c", "%3d", 3, P_NROTA2, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ICUT (CUTOFF METHOD), ICUTAT (SOLUTE ATOMS USED FOR ICUT = 4 - MAX 9)"},
	{"ICUT  ", "2", NULL, "%3c",  "%3d",   3, P_ICUT,   0, 0},
	{"ICUTAT", "",  NULL, "%45c", "%-45s", 45, P_ICUTAT, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"IRECNT(RECENTR) INDSOL(TRANS X's INDP) IZLONG MAXOVL NOXX NOSS NOBNDV NOANGV"},
	{"IRECNT", "1", NULL, "%3c", "%3d", 3, P_IRECNT, 0, 0},
	{"INDSOL", "0", NULL, "%3c", "%3d", 3, P_INDSOL, 0, 0},
	{"IZLONG", "0", NULL, "%3c", "%3d", 3, P_IZLONG, 0, 0},
	{"MAXOVL", "0", NULL, "%3c", "%3d", 3, P_MAXOVL, 0, 0},
	{"NOXX  ", "0", NULL, "%3c", "%3d", 3, P_NOXX,   BOSS_V38, 0},
	{"NOSS  ", "0", NULL, "%3c", "%3d", 3, P_NOSS,   BOSS_V38, 0},
	{"NOBNDV", "0", NULL, "%3c", "%3d", 3, P_NOBNDV, BOSS_V40, 0},
	{"NOANGV", "0", NULL, "%3c", "%3d", 3, P_NOANGV, BOSS_V40, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDF (# OF SOLVENT-SOLVENT RDFS), NRDFS (# OF SOLUTE-SOLVENT RDFS)"},
	{"NRDF  ", "0", NULL, "%3c", "%3d", 3, P_NRDF, 0, 0},
	{"NRDFS ", "0", NULL, "%3c", "%3d", 3, P_NRDFS, 0,
		"             MAXIMUM # IS 3 AND 12"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDFS1 (SOLVENT ATOMS FOR SOLVENT-SOLVENT RDFS IN ASCENDING ORDER)"},
	{"NRDFS1", " ", NULL, "%36c", "%-36s", 36, P_NRDFS1, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDFS2 (SOLVENT ATOMS FOR SOLVENT-SOLVENT RDFS IN ASCENDING ORDER)"},
	{"NRDFS2", " ", NULL, "%36c", "%-36s", 36, P_NRDFS2, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDFA1 (SOLUTE ATOMS FOR SOLUTE-SOLVENT RDFS IN ASCENDING ORDER)"},
	{"NRDFA1", " ", NULL, "%36c", "%-36s", 36, P_NRDFA1, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDFA2 (SOLVENT ATOMS FOR SOLUTE-SOLVENT RDFS IN ASCENDING ORDER)"},
	{"NRDFA2", " ", NULL, "%36c", "%-36s", 36, P_NRDFA2, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDLMIN, RDLINC (MINIMUM AND INCREMENT FOR THE RDFS IN ANGSTROMS)"},
	{"RDLMIN", "2.45", NULL, "%10c", "% 10.4f", 10, P_RDLMIN, 0, 0},
	{"RDLINC", "0.10", NULL, "%10c", "% 10.4f", 10, P_RDLINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EPRMIN, EPRINC (MIN AND INC FOR SOLVENT-SOLVENT ENERGY PAIR DIST.)"},
	{"EPRMIN", "-10.25", NULL, "%10c", "% 10.4f", 10, P_EPRMIN, 0, 0},
	{"EPRINC", "0.5",    NULL, "%10c", "% 10.4f", 10, P_EPRINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EDMIN, EDINC   (MIN AND INC FOR SOLVENT-SOLVENT ENERGY DISTRIBUTION)"},
	{"EDMIN ", "-40.5", NULL, "%10c", "% 10.4f", 10, P_EDMIN, 0, 0},
	{"EDINC ", "1.0",   NULL, "%10c", "% 10.4f", 10, P_EDINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ESSMIN, ESSINC (MIN AND INC FOR SOLUTE-SOLVENT  ENERGY PAIR DIST.)"},
	{"ESSMIN", "-12.15", NULL, "%10c", "% 10.4f", 10, P_ESSMIN, 0, 0},
	{"ESSINC", "0.3",    NULL, "%10c", "% 10.4f", 10, P_ESSINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EBSMIN, EBSINC (MIN AND INC FOR SOLUTE-SOLVENT ENERGY DISTRIBUTION)"},
	{"EBSMIN", "-100.0", NULL, "%10c", "% 10.4f", 10, P_EBSMIN, 0, 0},
	{"EBSINC", "2.0",    NULL, "%10c", "% 10.4f", 10, P_EBSINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NVCHG, NSCHG (FREQ OF VOLUME AND SOLUTE MOVES) SET NVCHG=999999 FOR NVT"},
	{"NVCHG ", "0", NULL, "%6c", "%6d", 6, P_NVCHG, 0, 0},
	{"NSCHG ", "5", NULL, "%6c", "%6d", 6, P_NSCHG, 0, 0},
	{"MAXVAR", "0", NULL, "%6c", "%6d", 6, P_MAXVAR, 0,
		"       MAXVAR (3I6) WILL BE COMPUTED FROM 0 DEFAULT VALUES"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NCONSV (FREQ OF COORD OUPUT), NBUSE=1 (USE SOLVENT NB LIST), NOCUT, NOSMTH"},
	{"NCONSV", "999999", NULL, "%6c", "%6d", 6, P_NCONSV, 0, 0},
	{"NBUSE ", "0",      NULL, "%6c", "%6d", 6, P_NBUSE,  0, 0},
	{"NOCUT ", "0",      NULL, "%6c", "%6d", 6, P_NOCUT,  BOSS_V40, 0},
	{"NOSMTH", "0",      NULL, "%6c", "%6d", 6, P_NOSMTH, BOSS_V40, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"VDEL (SIZE OF VOLUME MOVES), WKC (FOR PREF. SAMPLING) - BOTH DEFAULT 0"},
	{"VDEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_VDEL, 0, 0},
	{"WKC   ", "0", NULL, "%10c", "% 10.4f", 10, P_WKC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDEL, ADEL (RANGES FOR SOLVENT TRANSLATIONS AND ROTATIONS), RADIUS FOR SASA"},
	{"RDEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_RDEL,   0, 0},
	{"ADEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_ADEL,   0, 0},
	{"RADSOL", "0", NULL, "%10c", "% 10.4f", 10, P_RADSOL, BOSS_V40, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDELS1, ADELS1 (RANGES FOR TRANSLATIONS AND ROTATIONS OF SOLUTE 1)"},
	{"RDELS1", "0.08", NULL, "%10c", "% 10.4f", 10, P_RDELS1, 0, 0},
	{"ADELS1", "8.0",  NULL, "%10c", "% 10.4f", 10, P_ADELS1, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDELS2, ADELS2 (RANGES FOR TRANSLATIONS AND ROTATIONS OF SOLUTE 2)"},
	{"RDELS2", "0.08", NULL, "%10c", "% 10.4f", 10, P_RDELS2, 0, 0},
	{"ADELS2", "8.0",  NULL, "%10c", "% 10.4f", 10, P_ADELS2, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, BOSS_V38,
		"SOLUTE NUMBER (I4) and TEMPERATURE (F10.4) FOR LOCAL HEATING"},
	{"LHTSOL", "0",      NULL, "%4c",  "%4d",      4, P_LHTSOL, BOSS_V38, 0},
	{"TLHT  ", "1000.0", NULL, "%10c", "% 10.4f", 10, P_TLHT,   BOSS_V38, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RCUT, SCUT, CUTNB (CUTOFFS FOR SOLVENT-SOLVENT, SOLUTE-SOLVENT & INTRASOLUTE"},
	{"RCUT  ", "9.0",  NULL, "%10c", "% 10.4f", 10, P_RCUT, 0, 0},
	{"SCUT  ", "9.0",  NULL, "%10c", "% 10.4f", 10, P_SCUT, 0, 0},
	{"CUTNB ", "100.0",NULL, "%10c", "% 10.4f", 10, P_CUTNB, 0,
		"                NONBONDED INTERACTIONS"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"T, P (TEMPERATURE IN DEGREES C, AND PRESSURE IN ATMOSPHERES), TORCUT"},
	{"T     ", "25.0", NULL, "%10c", "% 10.4f", 10, P_T, 0, 0},
	{"P     ", "1.0",  NULL, "%10c", "% 10.4f", 10, P_P, 0, 0},
	{"TORCUT", "0.0",  NULL, "%10c", "% 10.4f", 10, P_TORCUT, BOSS_V36, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"DIELEC CONST (FOR VACUUM/CONTINUUM SIMULATIONS) & SCL14C, SCL14L, e (RXN FIELD)"},
	{"DIELEC", "1.0", NULL, "%10c", "%10.4f", 10, P_DIELEC, 0, 0},
	{"SCL14C", "2.0", NULL, "%10c", "%10.4f", 10, P_SCL14C, 0, 0},
	{"SCL14L", "2.0", NULL, "%10c", "%10.4f", 10, P_SCL14L, 0, 0},
	{"DIELRF", "1.0", NULL, "%10c", "%10.4f", 10, P_DIELRF, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"PLT FILE FORMAT (MIND,PDB,PDB2,PDBB), ISOLEC(SOLUTE FOR ENERGY COMPONENTS)(1A5,I4)"},
	{"PLTFMT", "PDB",  NULL, "%5c", "%-5s", 5, P_PLTFMT, 0, 0},
	{"ISOLEC", "0001", NULL, "%4c", "%04d", 4, P_ISOLEC, BOSS_V36, 0},
	};

#define MAXBOSSPARENTRY	(sizeof(bossparentry)/sizeof(BossParEntry))

static BossParEntry	mcproparentry[] = {
	{0, 0, 0, 0, 0, 0, 0, 0,
		"SET-UP AND POTENTIAL FUNCTION PARAMETERS:"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"SVMOD1 (SOLVENT MODEL) - TIP3P, TIP4P, CH3OH, DMSO"},
	{"SVMOD1", "TIP4P",  NULL, "%5c", "%-5s", 5, P_SVMOD1, 0,
		"  (1A5)    THF, CCL4, ARGON, CH3CN, MEOME, C3H8, MECL2, CHCL3, OTHER or NONE"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"ORIGIN OF INITIAL SOLUTE COORDINATES - FROM Z-MATRIX, PDB FILE, MIND FILE,"},
	{"SLFMT ", "BOXES",  NULL, "%5c", "%-5s", 5, P_SLFMT, 0,
		"    (1A5)    IN FILE, ZMATRIX & IN FILE (ZMAT, PDB, MIND, IN, ZIN)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ORIGIN OF INITIAL SOLVENT COORDINATES - FROM STORED BOXES, IN FILE or NONE"},
	{"SOLOR ", "ZMAT",   NULL, "%5c", "%-5s", 5, P_SOLOR, 0,
		"    (1A5)                               (BOXES, IN, NONE)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"FOR SOLVENT CAP - ICAPAT (I4), CAP RADIUS ABOUT ICAPAT, RESTORING FORCE CONST (2F10.4)"},
	{"ICAPAT", "0", NULL, "%4c",  "%04d",    4, P_ICAPAT, 0, 0},
	{"CAPRAD", "0", NULL, "%10c", "%10.4f", 10, P_CAPRAD, 0, 0},
	{"CAPFK ", "",  NULL, "%10c", "%10.4f", 10, P_CAPFK, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"OPTIMIZER (A6),TOLERANCE(F10.4),NEWZM(A5),NSAEVL,NSATR(2I8),SATEMP,SART(2F10.4)"},
	{"OPTIM ", "SIMPLX", NULL, "%6c",  "%-6s",    6, P_OPTIM, 0, 0},
	{"FTOL  ", "0.001",  NULL, "%10c", "%9.6f ", 10, P_FTOL,   0, 0},
	{"NEWZM ", "NEWZM",  NULL, "%7c",  "%-7s",    7, P_NEWZM,  0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NMOL, IBOX, BCUT (2I4,F6.2) - # OF SOLVENT MOLECULES (9999 FOR AUTO),"},

	{"NMOL  ", "260", NULL, "%4c", "%4d",   4, P_NMOL, 0, 0},
	{"IBOX  ", "267", NULL, "%4c", "%4d",   4, P_IBOX, 0, 0},
	{"BCUT  ", "0.0", NULL, "%6c", "%6.2f", 6, P_BCUT, 0,
		"    STORED BOX SIZE, & BOX CUTOFF FOR AUTO (DEF = 7.0)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"2ND SOLVENT (SVMOD2, choices as above) & # OF 2ND SOLVENT MOLECULES (NMOL2)"},
	{"SVMOD2", "NONE", NULL, "%5c", "%-5s", 5, P_SVMOD2, 0, 0},
	{"NMOL2 ", "0",    NULL, "%4c", "%4d",  4, P_NMOL2,  0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NCENT1, NCENT2 (THE ATOMS USED TO DEFINE THE CENTER OF THE SOLUTE(S))"},
	{"NCENT1", "1", NULL, "%4c", "%4d", 4, P_NCENT1, 0, 0},
	{"NCENT2", "1", NULL, "%4c", "%4d", 4, P_NCENT2, 0,
		"        INTEGERS IN I3 UNLESS NOTED"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NROTA1, NROTA2  (THE ATOMS SOLUTES 1 & 2 ARE ROTATED ABOUT - NOTE:"},
	{"NROTA1", "1", NULL, "%4c", "%4d", 4, P_NROTA1, 0, 0},
	{"NROTA2", "1", NULL, "%4c", "%4d", 4, P_NROTA2, 0,
		"           ROTATIONS DO NOT CHANGE THEIR RELATIVE POSITIONS)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NOBACK = 1 IF THE BACKBONE FOR SOLUTE 1 IS FIXED"},
	{"NOBACK", "1", NULL, "%4c", "%4d", 4, P_NOBACK, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"IRECNT(RECENTER X) INDSOL(TRANSLATE X's INDEP.) IZLONG(ORIENT X ALONG Z)"},
	{"IRECNT", "1", NULL, "%4c", "%4d", 4, P_IRECNT, 0, 0},
	{"INDSOL", "1", NULL, "%4c", "%4d", 4, P_INDSOL, 0, 0},
	{"IZLONG", "1", NULL, "%4c", "%4d", 4, P_IZLONG, 0, 0},
	{"MAXOVL", "0", NULL, "%4c", "%4d", 4, P_MAXOVL, 0,
		"             MAXOVL = 1(MAX OVERLAY SOLUTES) (4I4)"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NRDF (# OF SOLVENT-SOLVENT RDFS) - MAXIMUM # IS 3"},
	{"NRDF  ", "0", NULL, "%4c", "%4d", 4, P_NRDF, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"IBDPDB = 1 TO AVOID FALSE BONDS BETWEEN ATOMS IN DIFFERENT RESIDUES"},
	{"IBDPDB", "1", NULL, "%4c", "%4d", 4, P_IBDPDB, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"The following line is blank for compatibility with BOSS"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		""},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDLMIN, RDLINC (MINIMUM AND INCREMENT FOR THE RDFS IN ANGSTROMS)"},
	{"RDLMIN", "2.45", NULL, "%10c", "% 10.4f", 10, P_RDLMIN, 0, 0},
	{"RDLINC", "0.10", NULL, "%10c", "% 10.4f", 10, P_RDLINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EPRMIN, EPRINC (MIN AND INC FOR SOLVENT-SOLVENT ENERGY PAIR DIST.)"},
	{"EPRMIN", "-10.25", NULL, "%10c", "% 10.4f", 10, P_EPRMIN, 0, 0},
	{"EPRINC", "0.50",   NULL, "%10c", "% 10.4f", 10, P_EPRINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EDMIN, EDINC   (MIN AND INC FOR SOLVENT-SOLVENT ENERGY DISTRIBUTION)"},
	{"EDMIN ", "-40.50", NULL, "%10c", "% 10.4f", 10, P_EDMIN, 0, 0},
	{"EDINC ", "1.00",   NULL, "%10c", "% 10.4f", 10, P_EDINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"ESSMIN, ESSINC (MIN AND INC FOR SOLUTE-SOLVENT  ENERGY PAIR DIST.)"},
	{"ESSMIN", "-12.15", NULL, "%10c", "% 10.4f", 10, P_ESSMIN, 0, 0},
	{"ESSINC", "0.30",   NULL, "%10c", "% 10.4f", 10, P_ESSINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"EBSMIN, EBSINC (MIN AND INC FOR SOLUTE-SOLVENT ENERGY DISTRIBUTION)"},
	{"EBSMIN", "-100.0", NULL, "%10c", "% 10.4f", 10, P_EBSMIN, 0, 0},
	{"EBSINC", "2.00",   NULL, "%10c", "% 10.4f", 10, P_EBSINC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NVCHG, NSCHG, NSCHG2 (FREQUENCIES OF VOLUME, SOLUTE & 2ND SOLUTE MOVES),"},
	{"NVCHG ", "0",  NULL, "%6c", "%6d", 6, P_NVCHG, 0, 0},
	{"NSCHG ", "15", NULL, "%6c", "%6d", 6, P_NSCHG, 0, 0},
	{"NSCHG2", "0",  NULL, "%6c", "%6d", 6, P_NSCHG2, 0,
		"       (3I6) WILL BE COMPUTED FROM 0 DEFAULT VALUES"},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"NCONSV       (FREQUENCY OF WRITING ALL COORDINATES TO THE SAVE FILE - I6)"},
	{"NCONSV", "999999", NULL, "%6c", "%6d", 6, P_NCONSV, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"VDEL (SIZE OF VOLUME MOVES), WKC (FOR PREF. SAMPLING) - BOTH DEFAULT 0"},
	{"VDEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_VDEL, 0, 0},
	{"WKC   ", "0", NULL, "%10c", "% 10.4f", 10, P_WKC, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDEL, ADEL (RANGES FOR SOLVENT TRANSLATIONS AND ROTATIONS) - BOTH DEFAULT 0"},
	{"RDEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_RDEL, 0, 0},
	{"ADEL  ", "0", NULL, "%10c", "% 10.4f", 10, P_ADEL, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDELS1, ADELS1 (RANGES FOR TRANSLATIONS AND ROTATIONS OF SOLUTE 1)"},
	{"RDELS1", "0.08", NULL, "%10c", "% 10.4f", 10, P_RDELS1, 0, 0},
	{"ADELS1", "8.0",  NULL, "%10c", "% 10.4f", 10, P_ADELS1, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RDELS2, ADELS2 (RANGES FOR TRANSLATIONS AND ROTATIONS OF SOLUTE 2)"},
	{"RDELS2", "0.08", NULL, "%10c", "% 10.4f", 10, P_RDELS2, 0, 0},
	{"ADELS2", "8.0",  NULL, "%10c", "% 10.4f", 10, P_ADELS2, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"RCUT, SCUT, CUTNB (CUTOFFS FOR SOLVENT-SOLVENT, SOLUTE-SOLVENT & INTRASOLUTE"},
	{"RCUT  ", "9.0",  NULL, "%10c", "% 10.4f", 10, P_RCUT, 0, 0},
	{"SCUT  ", "9.0",  NULL, "%10c", "% 10.4f", 10, P_SCUT, 0, 0},
	{"CUTNB ", "100.0",NULL, "%10c", "% 10.4f", 10, P_CUTNB, 0,
		"                NONBONDED INTERACTIONS"},
	{0, 0, 0, 0, 0, 0, 0, 0,
		"T, P (TEMPERATURE IN DEGREES C, AND PRESSURE IN ATMOSPHERES)"},
	{"T     ", "25.0", NULL, "%10c", "% 10.4f", 10, P_T, 0, 0},
	{"P     ", "1.0",  NULL, "%10c", "% 10.4f", 10, P_P, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"DIELECTRIC CONSTANT (ONLY FOR VACUUM & CONTINUUM SIMULATIONS) & SCL14C, SCL14L"},
	{"DIELEC", "1.0", NULL, "%10c", "%10.4f", 10, P_DIELEC, 0, 0},
	{"SCL14C", "2.0", NULL, "%10c", "%10.4f", 10, P_SCL14C, 0, 0},
	{"SCL14L", "2.0", NULL, "%10c", "%10.4f", 10, P_SCL14L, 0, 0},

	{0, 0, 0, 0, 0, 0, 0, 0,
		"PLT FILE FORMAT (MIND,PDB,PDB2,PDBB), ISOLEC(SOLUTE FOR ENERGY COMPONENTS)(1A5,I4)"},
	{"PLTFMT", "PDB",  NULL, "%5c", "%-5s", 5, P_PLTFMT, 0, 0},
	{"ISOLEC", "0001", NULL, "%4c", "%04d", 4, P_ISOLEC, 0, 0},
	};

#define MAXMCPROPARENTRY	(sizeof(mcproparentry)/sizeof(BossParEntry))

BossPar	*CreateBossPar (int flags)
{
	int	i, maxentry, parsize;
	BossPar	*bp;
	BossParEntry	*p;

	maxentry = (flags & B_MCPRO) ? MAXMCPROPARENTRY : MAXBOSSPARENTRY;
	p = (flags & B_MCPRO) ? mcproparentry : bossparentry;
	parsize = sizeof(BossParEntry)*maxentry;

	if (!(bp = (BossPar *)calloc(1, sizeof(BossPar))) ||
	    !(bp->par = (BossParEntry *)calloc(1, parsize))) {
		if (bp) free((char *)bp);
		return NULL;
	}
	BCOPY(p, bp->par, parsize);
	bp->flags = flags;
	bp->maxentry = maxentry;
	bp->version = 0;

	for(i=0,p=bp->par;i<maxentry;i++,p++) {
		if (p->size > 0) {
			p->val = (char *)calloc(1, p->size+1);
			strncpy(p->val, p->defval, p->size);
		}
	}
	return bp;
}

void	FreeBossPar (BossPar **bp)
{
	int	i;
	BossParEntry	*p;

	if (!bp || !*bp) return;
	for(i=0,p=(*bp)->par;i<(*bp)->maxentry;i++,p++) {
		if (p->val) free(p->val);
	}
	*bp = NULL;
}

BossPar	*CopyBossPar (BossPar *oldbp)
{
	int	i;
	BossPar	*bp;
	BossParEntry	*p, *op;

	if (!oldbp) return NULL;
	if (!(bp = CreateBossPar(oldbp->flags))) return NULL;
	p = bp->par;
	BCOPY(oldbp, bp, sizeof(BossPar));
	bp->par = p;

	for(i=0,p=bp->par,op=oldbp->par;i<bp->maxentry;i++,p++,op++) {
		if (p->val && op->val) BCOPY(op->val, p->val, p->size+1);
	}

	return bp;
}

BossParEntry	*FindBossParEntry (BossPar *bp, int param_id)
{
	int	i;
	BossParEntry	*p;

	if (!bp) return NULL;
	for(i=0,p=bp->par;i<bp->maxentry;i++,p++) {
		if (p->id == param_id) return p;
	}

	return NULL;
}

char	*GetBossParDefaultValue (BossPar *bp, int param_id)
{
	BossParEntry	*p;
	if (!(p = FindBossParEntry(bp, param_id))) return NULL;
	return p->defval;
}

char	*GetBossParValue (BossPar *bp, int param_id)
{
	BossParEntry	*p;
	if (!(p = FindBossParEntry(bp, param_id))) return NULL;
	return p->val;
}

void	SetBossParValue (BossPar *bp, int param_id, char *val)
{
	BossParEntry	*p;
	if (!(p = FindBossParEntry(bp, param_id)) || p->size <= 0) return;
	strncpy(p->val, val ? val : "", p->size);
	p->val[p->size] = '\0';
}

int	FGuessBossParVersion (FILE *fp)
{
	int	line, version;
	fpos_t	pos;
	char	card[128];

	if (!fp) return 0;
	if (fp == stdin) return BOSS_VMAX-1;
	rewind(fp);
	fgetpos(fp, &pos);
	for(line=0;READLINE(fp, card);line++) {
		if (strncmp(card, "TYP ", 4) == 0 || strncmp(card, " TYP ", 5) == 0) break;
	}
	fsetpos(fp, &pos);

	switch (line) {
	case 66: version = BOSS_V35; break;
	case 68: version = BOSS_V36; break;
	case 70: version = BOSS_V38; break;
	case 72:
	default:
		version = BOSS_VMAX-1;
		break;
	}
	return version;
}

BossPar	*FLoadBossPar (FILE *fp, int version, int flags)
{
	int	i;
	char	str[256], *c;
	BossPar	*bp;
	BossParEntry	*p;

	if (!fp || !(bp = CreateBossPar(flags))) return NULL;
	bp->version = version;
	for(i=0,p=bp->par;i<bp->maxentry;) {
		if (!p->name || !p->scanfmt) {
			i++;
			p++;
			if (p->version == 0 || version >= p->version) {
				if (!READLINE(fp,str)) break;;
			}
			continue;
		}
		if (p->version == 0 || version >= p->version) {
			if (!READLINE(fp, str)) break;
			if ((c = strrchr(str, '\n'))) *c = '\0';	/* chop newline char */
			c = str;
			for(;i<bp->maxentry;i++,p++) {
				if (!p->name || !p->scanfmt) break;
				c = strscan(c, p->scanfmt, p->val);
			}
		} else {
			for(;i<bp->maxentry;i++,p++) {
				if (!p->name || !p->scanfmt) break;
			}
		}
	}
	return bp;
}

void	FPrintBossPar (FILE *fp, BossPar *bp, int version)
{
	int	i, j, k, n;
	float	x;
	BossParEntry	*p;

	if (!fp || !bp) return;
	if (bp->flags & B_MCPRO) fprintf(fp, "MCPRO%d ", version);
	else fprintf(fp, "BOSS%d ", version);
	for(i=0,p=bp->par;i<bp->maxentry;) {
		if (!p->name || !p->scanfmt) {
			if (p->version == 0 || version >= p->version) {
				fprintf(fp, "%s\n", p->comment);
			}
			i++;
			p++;
			continue;
		}
		if (p->version == 0 || version >= p->version) {
			for(j=0;i<bp->maxentry;i++,p++,j++) {
				if (!p->name || !p->scanfmt || !p->printfmt) break;

				if (strchr(p->printfmt, 'd')) {
					sscanf(p->val, "%d", &n);
					fprintf(fp, p->printfmt, n);
				} else if (strchr(p->printfmt, 'f') ||
						strchr(p->printfmt, 'g') ||
						strchr(p->printfmt, 'G') ||
						strchr(p->printfmt, 'e') ||
						strchr(p->printfmt, 'E')) {
					sscanf(p->val, "%f", &x);
					fprintf(fp, p->printfmt, x);
				} else {
					for(k=0;k<p->size;k++) p->val[k] = toupper(p->val[k]);
					fprintf(fp, p->printfmt, p->val);
				}

				if (p->comment) fprintf(fp, "%s", p->comment);
			}
			if (j!=i) fprintf(fp, "\n");
		} else {
			for(j=0;i<bp->maxentry;i++,p++,j++) {
				if (!p->name || !p->scanfmt || !p->printfmt) break;
			}
		}
	}
}

/* fp = file pointer
   from_zmat = 1 --- parameters are being read from Z-matrix file.
                     Do not wind file. start reading atom types right away.
                     Atom types at the end of Z-matrix may have 4 digits.
               0 --- find a line beginning with 'TYP'
*/
static ListPtr	_FLoadBossAtomTypes (FILE *fp, int from_zmat)
{
	int	i;
	char	str[256], type[256], *c, *val[8];
	ListPtr	list, nextlist, atomtypelist=NULL;

	if (!fp) return NULL;
	clearmolerror();

	strcpy(scanfmt[0], "%3c% ");
	strcpy(scanfmt[1], "%2c% ");
	strcpy(scanfmt[2], "%2c% ");
	strcpy(scanfmt[3], "%10c");
	strcpy(scanfmt[4], "%10c");
	strcpy(scanfmt[5], "%10c");
	strcpy(scanfmt[6], "%64c");

	/* find out the beginning of the atom types */
	if (!from_zmat) {
		while (1) {
			if (!READLINE(fp, str)) {
	setmolerrortext("Can't find line containing \"TYP AN AT   CHARGE     SIGMA    EPSILON\"");
	catmolerror("                        or \" TYP AN AT   CHARGE     SIGMA    EPSILON\"");
				return (ListPtr)setmolerrorno(MERR_PARFORMAT);
			}

			if (strncmp(str, "TYP ", 4) == 0) break;
			if (strncmp(str, " TYP ", 5) == 0) {
			        strcpy(scanfmt[0], "%4c% ");
				break;
			}
			if (strstr(str, "EndBossPar:") || strstr(str, "EndCoord") ||
				strstr(str, "EndBossZmat")) return NULL;
		}
		/* read atom types */
		atomtypelist = FLoadBossList(fp, scanfmt, 7); 
	} else {
		if (!READLINE(fp, str)) return NULL;
		sscanf(str, "%4c", type);
		if (strncmp(type, " 800", 4) == 0 ||
		    strncmp(type, "0800", 4) == 0 ||
		    strncmp(type, "9500", 4) == 0) {	/* 4-digit atom types */
			strcpy(scanfmt[0], "%4c% ");
		}

		if (!(list = NewList())) {
			setmolerrorno(MERR_MEM);
			return NULL;
		}
		if ((c = strrchr(str, '\n'))) *c = '\0';
		c = str;
		val[0] = list->str1;
		val[1] = list->str2;
		val[2] = list->str3;
		val[3] = list->str4;
		val[4] = list->str5;
		val[5] = list->str6;
		val[6] = list->str7;
		val[7] = list->str8;
		for(i=0;i<7;i++) {
			c = strscan(c, scanfmt[i], val[i]);
			if (!*c) break;
		}

		/* read the rest of atom types */
		atomtypelist = FLoadBossList(fp, scanfmt, 7); 
		if (list) EnterList(list, &atomtypelist);
	}

	for(list=atomtypelist;list;list=nextlist) {
		nextlist = list->next;
		if (!sscanf(list->str1, "%d", &list->L_TYPE) ||
		    !sscanf(list->str2, "%d", &list->L_AN) ||
		    !sscanf(list->str4, "%f", &list->L_CHARGE) ||
		    !sscanf(list->str5, "%f", &list->L_SIGMA) ||
		    !sscanf(list->str6, "%f", &list->L_EPSILON)) {
			DeleteList(list, &atomtypelist);
			list->next = NULL;
			FreeList(&list);
			continue;
		}
		/* this is a heck...
		 * Make sure L_AMBER != str3 and L_COMMENT != str7
		 */
		strcpy(list->L_AMBER, strncmp(list->str3, "   ", 3) == 0 ? "??" : list->str3);
		strcpy(list->L_COMMENT, list->str7 + strspn(list->str7, " \t"));
	}
	if (getmolerrorno() != -1) {
		FreeList(&atomtypelist);
		return NULL;
	}
	return atomtypelist;
}

ListPtr	FLoadBossAtomTypes (FILE *fp)
{
	return _FLoadBossAtomTypes(fp, 0);
}

ListPtr	FLoadBossAtomTypesFromZmat (FILE *fp)
{
	return _FLoadBossAtomTypes(fp, 1);
}

void	FPrintBossOneAtomType (FILE *fp, ListPtr atomtype)
{
	char	amber[4];

	if (!fp || !atomtype) return;
	if (atomtype->L_TYPE <= 0) return;
	if (strlen(atomtype->L_AMBER) < 3) sprintf(amber, "%-2s ", atomtype->L_AMBER);
	else strncpy(amber, atomtype->L_AMBER, 3);
	amber[3] = '\0';

	if (strncmp(amber, "  ", 2) == 0 || strncmp(amber, "??", 2) == 0) {
		sprintf(amber, "%2s ", GetElemSymbol(atomtype->L_AN));
	}

	fprintf(fp, "%04d %02d %3s% 10.6f %9.6f %9.6f     %s\n",
		atomtype->L_TYPE, atomtype->L_AN, amber, atomtype->L_CHARGE,
		atomtype->L_SIGMA, atomtype->L_EPSILON, atomtype->L_COMMENT);

}

void	FPrintBossAtomTypes (FILE *fp, ListPtr atomtypes, int version)
{
	ListPtr	list;

	fprintf(fp, "  OPLS PARAMETERS FOR ORGANIC AND BIOCHEMICAL SYSTEMS\n");
	fprintf(fp, " TYP AN AT   CHARGE     SIGMA    EPSILON\n");
	ForEachList(atomtypes,list) {
		FPrintBossOneAtomType(fp, list);
	}
}

static void	process_torsion_comment (ListPtr list)
{
	char	str[256], *p, *tok, *atype[4];
	int	i, n, type[4];
	if (!list) return;
	str[0] = '\0';
	sscanf(list->L_COMMENT, "%*[^<]<%[^>\n]", str);
	if (!str[0] || (strchr(str, '[') && strchr(str, ']'))) return;

	p = str;
	for(n=0;(tok = strtok(p, "-"))!=(char *)NULL && n<4;n++) {
		p = (char *)NULL;
		atype[n] = tok;
	}
	if (n < 4) return;
	for(i=0;i<4;i++) {
		if (sscanf(atype[i], "%d", &type[i]) != 1) type[i] = 0;
	}
	list->L_A1 = type[0];
	list->L_A2 = type[1];
	list->L_A3 = type[2];
	list->L_A4 = type[3];
}

/* 
 * FLoadBossTorsion --- load OPLS torsion types from parameter file
 * fp = file pointer
 */

ListPtr	FLoadBossTorsion (FILE *fp)
{
	char	str[64];
	int	torsion_order_flag=1;
	ListPtr	list, nextlist, torsionlist=NULL;

	/* skip to the beginning of torsion types */
	while (1) {
		if (!READLINE(fp, str)) {
	setmolerrortext("Can't find a line containing \"Type   V1         V2    f   V3    f   V4    f\"");
			return (ListPtr)setmolerrorno(MERR_PARFORMAT);
		}
		if (strncasecmp(str, "Type", 4) == 0) {
			sscanf(str, "%*8c%d", &torsion_order_flag);
			break;
		}
		if (strstr(str, "EndCoord") || strstr(str, "EndBossPar") ||
			strstr(str, "EndBossZmat")) return NULL;
	}

	/* read torsion types */
	strcpy(scanfmt[0], "%15c");	/* "%3c%2 %8c% %c" */
	strcpy(scanfmt[1], "%8c% ");
	strcpy(scanfmt[2], "%c");
	strcpy(scanfmt[3], "%8c% ");
	strcpy(scanfmt[4], "%c");
	strcpy(scanfmt[5], "%8c% ");
	strcpy(scanfmt[6], "%c");
	strcpy(scanfmt[7], "%80c");

	torsionlist = FLoadBossList(fp, scanfmt, 8);

	for(list=torsionlist;list;list=nextlist) {
		nextlist = list->next;

		list->str1[3] = list->str1[13] = '\0';

		if (torsion_order_flag == 0) {
			if (!sscanf(list->str1,   "%d", &list->L_TYPE) ||
			    !sscanf(list->str1+5, "%f", &list->L_V1) ||
			    !sscanf(list->str1+14,"%d", &list->L_N1) ||
			    !sscanf(list->str2,   "%f", &list->L_V2) ||
			    !sscanf(list->str3,   "%d", &list->L_N2) ||
			    !sscanf(list->str4,   "%f", &list->L_V3) ||
			    !sscanf(list->str5,   "%d", &list->L_N3) ||
			    !sscanf(list->str6,   "%f", &list->L_V0)) {
				/* if L_COMMENT != str8, we should copy str8 to L_COMMENT too */
				DeleteList(list, &torsionlist);
				list->next = NULL;
				FreeList(&list);
				continue;
			}
		} else {
			if (!sscanf(list->str1,   "%d", &list->L_TYPE) ||
			    !sscanf(list->str1+5, "%f", &list->L_V0) ||
			    !sscanf(list->str2,   "%f", &list->L_V1) ||
			    !sscanf(list->str3,   "%d", &list->L_N1) ||
			    !sscanf(list->str4,   "%f", &list->L_V2) ||
			    !sscanf(list->str5,   "%d", &list->L_N2) ||
			    !sscanf(list->str6,   "%f", &list->L_V3) ||
			    !sscanf(list->str7,   "%d", &list->L_N3)) {
				/* if L_COMMENT != str8, we should copy str8 to L_COMMENT too */
				DeleteList(list, &torsionlist);
				list->next = NULL;
				FreeList(&list);
				continue;
			}
		}
		list->L_P1 = (float)list->L_N1 * M_PI_2;
		list->L_P2 = (float)list->L_N2 * M_PI_2;
		list->L_P3 = (float)list->L_N3 * M_PI_2;

		strcpy(list->str1, list->str8);	/* temporary place to hold str8 */
		strcpy(list->L_COMMENT, list->str1 + strspn(list->str1, " \t"));

		/* extract [A1-A2-A3-A4] part from the comment, if any */
		process_torsion_comment(list);
	}
	if (getmolerrorno() != -1) {
		FreeList(&torsionlist);
		return NULL;
	}
	return torsionlist;
}

void	FPrintBossOneTorsionType (FILE *fp, ListPtr list)
{
	if (!fp || !list) return;

	fprintf(fp, "%03d  % 8.4f  % 8.4f %d% 8.4f %d% 8.4f %d   %s\n",
		list->L_TYPE, list->L_V0, list->L_V1, list->L_N1, list->L_V2, list->L_N2,
		list->L_V3, list->L_N3, list->L_COMMENT);
}

void	FPrintBossTorsionTypes (FILE *fp, ListPtr torsiontypes, int version)
{
	ListPtr	list;

	fprintf(fp, "      More charge and L-J parameters should be added above this line.\n");
	fprintf(fp, "      The following are the Fourier coefficients for dihedral angles.\n");
	fprintf(fp, "Type   V1         V2    f   V3    f   V4    f      Description\n");

	ForEachList(torsiontypes,list) {
		FPrintBossOneTorsionType(fp, list);
	}
}

ListPtr	FindBossAtomType (ListPtr atomtypelist, int type)
{
	ListPtr	list;
	ForEachList(atomtypelist,list) {
		if (list->L_TYPE == type) return list;
	}
	return NULL;
}

void	AssignBossPar (ListPtr zmatlist, ListPtr atomtypelist)
{
	ZmatPtr	zmat;
	ListPtr	list;
	int	unknown_type=0, nerr=0;

	ForEachZmat(LinkZmatList(zmatlist), zmat) {
		if (zmat->itype == -1) zmat->an = -1;
		else if ((list = FindBossAtomType(atomtypelist, zmat->itype))) zmat->an = list->L_AN;
		else {
			if (zmat->itype > 0 && zmat->itype <= 53) {
				zmat->an = zmat->itype;
				continue;
			}

			zmat->an = 2;
			if (nerr > 10) continue;
			if (nerr == 10) {
				setmolerror("More errors... Check parameter file if atom types exist.");
				printmolerror();
			} else {
				setmolerror("Atom type %d not found in parameter file. Was set to He.",
					zmat->itype);
				printmolerror();
			}
			if (nerr == 0) unknown_type = zmat->itype;
			nerr++;
		}
	}
	if (nerr) {
		setmolerror("WARNING: OPLS atom types (%d etc) not found. Ignored.", unknown_type);
	}
	BreakZmatList(zmatlist);
}

void	FPrintBossParFile (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes, int version)
{
	if (!fp || !bp) return;
	FPrintBossPar(fp, bp, version);
	if (!atomtypes) return;
	FPrintBossAtomTypes(fp, atomtypes, version);
	if (torsiontypes) FPrintBossTorsionTypes(fp, torsiontypes, version);
}

void	FPrintBossParFileWithZmatTypes (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes,
	ListPtr zmat_atomtypes, ListPtr zmat_torsiontypes, int version)
{
	ListPtr	list;

	if (!fp || !bp) return;
	FPrintBossPar(fp, bp, version);
	if (!atomtypes && !zmat_atomtypes) return;

	fprintf(fp, "  OPLS PARAMETERS FOR ORGANIC AND BIOCHEMICAL SYSTEMS\n");
	fprintf(fp, " TYP AN AT   CHARGE     SIGMA    EPSILON\n");

	ForEachList(atomtypes,list) {
		if (FindListData(zmat_atomtypes, list->L_TYPE)) {
			continue;
		}
		FPrintBossOneAtomType(fp,list);
	}
	ForEachList(zmat_atomtypes,list) {
		FPrintBossOneAtomType(fp,list);
	}

	if (!torsiontypes && !zmat_torsiontypes) return;

	fprintf(fp, "      More charge and L-J parameters should be added above this line.\n");
	fprintf(fp, "      The following are the Fourier coefficients for dihedral angles.\n");
	fprintf(fp, "Type   V1         V2    f   V3    f   V4    f      Description\n");

	ForEachList(torsiontypes,list) {
		if (FindListData(zmat_torsiontypes, list->L_TYPE)) {
			continue;
		}
		FPrintBossOneTorsionType(fp,list);
	}
	ForEachList(zmat_torsiontypes,list) {
		FPrintBossOneTorsionType(fp,list);
	}
}

int	FLoadBossParFile (FILE *fp, int flags, BossParPtr *bp_ret, ListPtr *at_ret, ListPtr *tt_ret, int *version_ret)
{
	if (*version_ret < 0 || *version_ret >= BOSS_VMAX) *version_ret = FGuessBossParVersion(fp);
	*bp_ret = FLoadBossPar(fp, *version_ret, flags);
	(*at_ret = FLoadBossAtomTypes(fp)) && (*tt_ret = FLoadBossTorsion(fp));

	return (*bp_ret && *at_ret && *tt_ret);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Below are used for loading "oplsaa.par" (or BOSS input parameter file).
 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
OplsTypeTorsion	*FLoadOplsTypeTorsion (FILE *fp)
{
	OplsTypeTorsion	*opls;

	clearmolerror();
	if (!(opls = (OplsTypeTorsion *)calloc(1, sizeof(OplsTypeTorsion)))) return NULL;
	opls->atomtypelist = FLoadBossAtomTypes(fp);
	opls->torsionlist = FLoadBossTorsion(fp);

	return opls;
}

OplsTypeTorsion	*FLoadDefaultOplsTypeTorsion (FILE *fp)
{
	OplsTypeTorsion	*opls;

	if (!(opls = FLoadOplsTypeTorsion(fp))) return NULL;
	if (opls_typetorsion) FreeOplsTypeTorsion(&opls_typetorsion);
	opls_typetorsion = opls;
	return opls;
}

static OplsTypeTorsion	*_LoadOplsTypeTorsion (char *filename, int default_flag)
{
	FILE	*fp;
	OplsTypeTorsion	*p;

	clearmolerror();
	if (!filename || !*filename) return NULL;
	if (!(fp = fopen(filename, "r"))) {
		setmolerror("Can't open parameter file %s.", filename);
		return NULL;
	}

	if (default_flag) p = FLoadDefaultOplsTypeTorsion(fp);
	else p = FLoadOplsTypeTorsion(fp);

	fclose(fp);
	return p;
}

OplsTypeTorsion	*LoadOplsTypeTorsion (char *filename)
{
	return _LoadOplsTypeTorsion(filename, 0);
}

OplsTypeTorsion	*LoadDefaultOplsTypeTorsion (char *filename)
{
	return _LoadOplsTypeTorsion(filename, 1);
}

OplsTypeTorsion	*GetDefaultOplsTypeTorsion (void)
{
	return opls_typetorsion;
}

void	FreeOplsTypeTorsion (OplsTypeTorsion **p)
{
	OplsTypeTorsion	*q;

	if (!p || !(q = *p)) return;
	if (q->atomtypelist) FreeList(&q->atomtypelist);
	if (q->torsionlist) FreeList(&q->torsionlist);
	free(*p);
	*p = NULL;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Below are used for loading "oplsaa.sb".
 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
ListPtr	FLoadOplsStretch (FILE *fp)
{
	ListPtr	list, stretchlist=NULL;
	char	*p, type1[8], type2[8], comment[128], str[256];
	float	Rk, Req;

	if (!fp) return NULL;
	while (READLINE(fp, str)) {
		if (strncmp(str, "   ", 3) == 0) break;	/* begining of bond angle bending */
		if (strchr(" \t\n", str[0])) break;
		if (strncmp(str, "***", 3) == 0) continue;
		strscan(str, "%2c% %2c%10f%10f%80c", type1, type2, &Rk, &Req, comment);
		if (!type1[0] || strncmp(type1, "  ", 2) == 0) continue;
		if (!(list = EnterNewList(&stretchlist))) return stretchlist;
		list->L_K = Rk;
		list->L_EQUIL = Req;
		strncpy(list->L_AMBER1, type1, 2);
		strncpy(list->L_AMBER2, type2, 2);
		strncpy(list->L_COMMENT, comment, 80);
		list->L_AMBER1[2] = '\0';
		list->L_AMBER2[2] = '\0';
	}
	return stretchlist;
}

ListPtr	FLoadOplsBending (FILE *fp, int rewind_flag)
{
	ListPtr	list, bendinglist=NULL;
	char	type1[8], type2[8], type3[8], comment[128], str[256];
	float	Tk, Teq;

	if (!fp) return NULL;
	if (rewind_flag) {
		while (1) {
			if (!READLINE(fp, str)) {
				setmolerror("Error in stretch-bend file.");
				return NULL;
			}
			if (strncmp(str, "   ", 3) == 0) break;	/* begining of bond angle bending */
		}
	}

	while (READLINE(fp, str)) {
		if (strncmp(str, "***", 3) == 0) continue;
		strscan(str, "%2c% %2c% %2c%2 %10f%10f%80c", type1, type2, type3, &Tk, &Teq, comment);
		if (!type1[0] || strncmp(type1, "  ", 2) == 0) continue;
		if (!(list = EnterNewList(&bendinglist))) return bendinglist;
		list->L_K = Tk;
		list->L_EQUIL = Teq;
		strncpy(list->L_AMBER1, type1, 2);
		strncpy(list->L_AMBER2, type2, 2);
		strncpy(list->L_AMBER3, type3, 2);
		strncpy(list->L_COMMENT, comment, 80);
		list->L_AMBER1[2] = '\0';
		list->L_AMBER2[2] = '\0';
		list->L_AMBER3[2] = '\0';
	}
	return bendinglist;
}

static OplsStretchBend	*_FLoadOplsStretchBend (FILE *fp, int default_flag)
{
	OplsStretchBend	*opls;

	clearmolerror();
	if (!(opls = (OplsStretchBend *)calloc(1, sizeof(OplsStretchBend)))) return NULL;
	opls->stretchlist = FLoadOplsStretch(fp);
	opls->bendinglist = FLoadOplsBending(fp, 0);

	if (default_flag) {
		if (opls_stretchbend) FreeOplsStretchBend(&opls_stretchbend);
		opls_stretchbend = opls;
	}

	return opls;
}

OplsStretchBend	*FLoadOplsStretchBend (FILE *fp)
{
	return _FLoadOplsStretchBend(fp, 0);
}

OplsStretchBend	*FLoadDefaultOplsStretchBend (FILE *fp)
{
	return _FLoadOplsStretchBend(fp, 1);
}

static OplsStretchBend	*_LoadOplsStretchBend (char *filename, int default_flag)
{
	FILE	*fp;
	OplsStretchBend	*p;

	clearmolerror();
	if (!filename || !*filename) return NULL;
	if (!(fp = fopen(filename, "r"))) {
		setmolerror("Can't open stretch bend file %s.", filename);
		return NULL;
	}

	if (default_flag) p = FLoadDefaultOplsStretchBend(fp);
	else p = FLoadOplsStretchBend(fp);

	fclose(fp);
	return p;
}

OplsStretchBend	*LoadOplsStretchBend (char *filename)
{
	return _LoadOplsStretchBend(filename, 0);
}

OplsStretchBend	*LoadDefaultOplsStretchBend (char *filename)
{
	return _LoadOplsStretchBend(filename, 1);
}

OplsStretchBend	*GetDefaultOplsStretchBend (void)
{
	return opls_stretchbend;
}

void	FreeOplsStretchBend (OplsStretchBend **p)
{
	OplsStretchBend	*q;

	if (!p || !(q = *p)) return;
	if (q->stretchlist) FreeList(&q->stretchlist);
	if (q->bendinglist) FreeList(&q->bendinglist);
	free(*p);
	*p = NULL;
}

/**********************************************************************
 ***** Default BOSS Parameter Generation Based on Calculation Type ****
 **********************************************************************/
/*
 * NPT ---> NVCHG = 0 (use default value)
 * NVT ---> NVCHG = 999999 (no volume move)
 */
void	AssignDefaultBossParByCalcType (BossParPtr bp, int calc_type)
{
	float	x;
	BossParEntry	*p;

	if (!bp || !(p = bp->par)) return;

	SetBossParValue(bp, P_NOSMTH, GetBossParDefaultValue(bp,P_NOSMTH));
	SetBossParValue(bp, P_NCENTS, GetBossParDefaultValue(bp,P_NCENTS));

	if (IS_CALC_MC_LIQ(calc_type)||IS_CALC_MC_SOL(calc_type)) {
		if (strcmp(GetBossParValue(bp,P_SVMOD1), "OTHER") == 0) {
			SetBossParValue(bp,P_NCENTS,"1");
		}
		if (sscanf(GetBossParValue(bp,P_DIELRF), "%f", &x) == 1 &&
			ABS(x-1.0) < 0.001) SetBossParValue(bp,P_NOSMTH,"1");
	}
}

static char	*boss_calc_type_text[] = {
	"opt", "opt", "freq", "sp",
	"dihed", "search",
	"mc_gas", "mc_liq",
	"mc_liq_cluster", "mc_sol",
	"mc_sol_cluster"};

char	*GetBossCalcTypeText (int calc_type)
{
	if (calc_type < 0 || calc_type >= BOSS_CALC_TYPE_MAX) calc_type = 0;
	return boss_calc_type_text[calc_type];
}

