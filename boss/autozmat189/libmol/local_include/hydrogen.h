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
