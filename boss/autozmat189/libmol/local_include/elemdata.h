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
