#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

/*********************************************************************************************
                         MacroModel atom types
----------------------------------------------------------------------------------------------
No. Symbol Description                     Equivalencies
                                        MM3     Charmm     Amber
----------------------------------------------------------------------------------------------
  1  C1    Csp                           4
  2  C2    Csp2                         2,3     CT,C                C,C*,CA,CB,CC,CM,CF,CG,CN
  3  C3    Csp3                          1      CT                  CT
  4  CA    United Atom CH  - sp3                CH1E                CH
  5  CB         "      CH2 - sp3                CH2E                C2
  6  CC         "      CH3 - sp3                CH3E                C3
  7  CD         "      CH  - sp2                CR1E                CD,CE,CJ,CP
  8  CE         "      CH2 - sp2
  9  CF         "      CH  - sp
 10  CM    Carbanion
 11  CP    Carbocation                   30
 12  CR    Carbon free radical           29
 13
 14  C0    Any carbon
 15  O2    Oxygen - double bond          7      O                   O,O2
 16  O3    
 17  OA    United Atom OH                       OH1E
 18  OM    O- (Alkoxide, carboxyalte)    47                         OC
 19  OW    United Atom H2O - Water              OH2E
 20  OP    Oxonium (sp2) =O(+)-
 21  OQ    Oxonium (sp3) R3O+ 
 22
 23  O0    Any oxygen
 24  N1    Nitrogen - sp                 10
 25  N2    N - sp2                        9     N,NR,NP,NH1,NH2     N,NB,NC,N*,N2
 26  N3    N - sp3                        8     NH3                 NT
 27  NA    United Atom NH  - sp3
 28  NB         "      NH2 - sp3                NH1E
 29
 30  ND    United Atom NH2 - sp2                NH2E
 31  N4    N+  - sp2                     46     NC2                 NA
 32  N5    N+  - sp3                     39                         N3
 33  NE    United Atom NH+  - sp3
 34  NF    United Atom NH2+ - sp3
 35  NG    United Atom NH3+ - sp3               NH3E
 36  NH    United Atom NH+  - sp2
 37  NI    United Atom NH2+ - sp2               NC2E
 38
 39
 40  N0    Any nitrogen
 41  H1    H - electroneut (e.g., C,S)  5,44    HA                  HC,HS
 42  H2    H - O (neut)                 21,24   H,HC                HO
 43  H3    H - N (neut)                  23     H,HC                H,H2
 44  H4    H - Cation                    48                         H3
 45  H5    H - Anion
 46
 47
 48  H0    Any Hydrogen
 49  S1    Sulfur                    15,17,18   S                   S,SH
 50  SA    United Atom SH                       SH1E
 51  SM    Thiolate anion
 52  S0    Any Sulfur                           S
 53  P0    Phosphorus                     25                        P
 54  B2    Boron (sp2)                    26
 55  B3    Boron (sp3)                    27
 56  F0    Fluorine                       11
 57  Cl    Chlorine                       12
 58  Br    Bromine                        13
 59  I0    Iodine                         14
 60  Si    Silicon                        19
 61  Du    Dummy atom for FEP
 62  Z0    Special Atom to be defined
 63  Lp    Lone pair
 64  00    Any atom                        0     *
 65  Li    Lithium
 66  Na    Sodium
 67  K     Potasium
----------------------------------------------------------------------------------------------
*/

#define MAX_MMODEL_TYPES	64
#define LONEPAIR		98

static int	mmodeltype[MAX_MMODEL_TYPES+1] = {0,
	 6,  6,  6,  6,  6,	6, 6, 6, 6, 6,	6, 6, 6, 6,		/*  1-14 C */
	 8,  8,  8,  8,  8,	8, 8, 8, 8,				/* 15-23 O */
	 7,  7,  7,  7,  7,	7, 7, 7, 7, 7,	7, 7, 7, 7, 7,	7, 7,	/* 24-40 N */
	 1,  1,  1,  1,  1,	1, 1, 1,				/* 41-48 H */
	16, 16, 16, 16,		/* 49-52 S */
	15,			/* 53    P */
	 5,  5, 		/* 54,55 B */
	 9,			/* 56    F */
	17,			/* 57   Cl */
	35,			/* 58   Br */
	53,			/* 59    I */
	14,			/* 60   Si */
	-1,			/* 61   DU */
	 0,
	LONEPAIR,		/* 63   LP */
	 0
	 };

/* return atomic number for a macromodel type "type" */
int	MModelToAn (int type)
{
	if (type > MAX_MMODEL_TYPES || type < 0) type = 0;
	return mmodeltype[type];
}

/* returns a mamcomodel atomtype based on atomic number (an),
the number of neighbors (nb), attached atom, etc.
Uses 62 rather than 64 for unknown
atom types as described in macromodel doc to avoid bugs in macromodel.
Input:
	atom = atom pointer
*/
int	GetMModelAtomType (AtomPtr atom)
{
	int	type;
	int	an, nb;
	AtomPtr	nbatom;

	if (!atom) return 0;
	an = atom->an;
	nb = atom->nneighbors;
	nbatom = atom->nbatom[0];

	switch (an) {
	case 0:
	case 1:
		type = 48;	/* any hydrogen */
		if (!nbatom) break;
		switch (nbatom->an) {
		case 6: type = 41; break;
		case 8: type = 42; break;
		case 7: type = 43; break;
		}
		break;
	case 6:
		switch (nb) {
		case 2: type = 1; break;
		case 3: type = 2; break;
		case 4: type = 3; break;
		default: type = 14;
		}
		break;
	case 7:
		switch (nb) {
		case 1: type = 24; break;
		case 2: type = 25; break;
		case 3: type = 26; break;
		case 4: type = 32; break;
		default: type = 40;
		}
		break;
	case 8:
		switch (nb) {
		case 1: type = 15; break;
		case 2: type = 16; break;
		default: type = 23;
		}
		break;
	case 9: type = 56; break;
	case 14: type = 60; break;
	case 15: type = 53; break;
	case 16:
		switch (nb) {
		case 2: type = 49; break;
		default: type = 52; break;
		}
		break;
	case 17: type = 57; break;
	case 35: type = 58; break;
	case 53: type = 59; break;
	case 98: type = 63; break;	/* LP */
	default: type = 61;		/* dummy */
	}
	return type;
}

/* MacroModel bond types
----------------------------------------
No.	Symbol	Description
----------------------------------------
 0	  .	zero-order bond
 1	  -	single bond
 2	  =	double bond
 3	  %	triple bond
 4	  *	undefined (excluding 0)
----------------------------------------
*/

int	GetMModelBondType (int bondtype)
{
	switch (bondtype) {
	case B_HBOND:
	case B_PARTIAL:
	case B_SINGLE: return 1;
	case B_AROMA:
	case B_DOUBLE: return 2;
	case B_DBLAROMA:
	case B_TRIPLE: return 3;
	default: return 0;
	}
}

/* MacroModel file format 
Header: (1X,I5,2X,A70) natoms, title
Atom Entries: (1X,I3,6(1X,I5,1X,I1),1X,3(F11.6,1X),I5,2A1,I4,2F9.5,1X,A4,1X,A4)
	I3	atomtype
	6(1X,I5,1X,I1)
		I5	sequence number of attached atom
		I1	bond order
	F11.6	x coord
	F11.6	y coord
	F11.6	z coord
	- optional -
	I5	residue number (used only with biopolymer)
	A1	single-char MacorModel residue code
	A1	single-char residue chain name
	I4	atom color
	F9.5	charge for standard charge/charge electrostatics
	F9.5	charge for charge/multipole electrostatistics
	A4	ascii residue name (PDB residue name)
	A4	ascii atom name (PDB atom name)

E.g.,
  41     5 1     0 0     0 0     0 0     0 0     0 0   -2.034824    -.118289     .807813     00A   0   .00000   .00000 RESI ATOM
*/

static char	*fmt = " %3d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %5d %1d %11lf %11lf %11lf";
static char	*opt_fmt = " %5d%c%c%4d%9lf%9lf %4c %4c";

/* read macromodel format and assign atomlist and bondlist.
*/
static AtomPtr	readmmodelline (char *str, Chain **top_chn, Chain **cur_chn, Residue **cur_res, AtomPtr *conn_atom, BondPtr *bondlist, int serial)
{
	int	an, n, j, hetatm, atom_type, atom_color, res_seq, bondorder[MMOD_MAX_BONDS], nb[MMOD_MAX_BONDS];
	double	x, y, z, cc_charge, cm_charge;
	char	res_id, chain_id, res_name[5], atom_name[5];
	AtomPtr	atom;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	str = strscan(str, fmt, &atom_type,
		&nb[0], &bondorder[0], &nb[1], &bondorder[1], &nb[2], &bondorder[2],
		&nb[3], &bondorder[3], &nb[4], &bondorder[4], &nb[5], &bondorder[5],
		&x, &y, &z);
	res_name[0] = '\0';
	atom_name[0] = '\0';
	if (str) strscan(str, opt_fmt,
		&res_seq, &res_id, &chain_id, &atom_color, &cc_charge, &cm_charge, res_name, atom_name);

	an = MModelToAn(atom_type);
	patom->x = x;
	patom->y = y;
	patom->z = z;
	patom->serial = serial;

	pres->seq = res_seq;
	if (atom_name[0] && !IsBlank(atom_name)) strncpy(patom->name, atom_name, 4);
	else strcpy(patom->name, AN2ATOMNAME(an));

	if (res_name[0] && !IsBlank(res_name)) strncpy(pres->name, res_name, 3);
	else strcpy(pres->name, "UNK");
	hetatm = 1;

	atom = CreatePdbAtom(top_chn, cur_chn, cur_res, conn_atom, patom, bondlist, hetatm);
	atom->an = an;
	atom->type = atom_type;
	atom->charge = cc_charge;

	n = 0;
	for(j=0;j<MMOD_MAX_BONDS;j++) {
		if (nb[j] <= 0) continue;
		if (nb[j] == serial) continue;
		/* nbatom and nbbond arrays store serial number and bondtype, temporarily */
		atom->nbatom[n] = (AtomPtr)nb[j];
		atom->nbbond[n] = (BondPtr)bondorder[n];
		n++;
	}
	atom->nneighbors = n;

	return atom;
}

static int	get_mmodel_header (char *str, int *natoms, char title[])
{
	char	*c;

	if (!natoms) return 0;
	*natoms = 0;
	if ((c = strrchr(str, '\n'))) *c = '\0';
	strscan(str, " %5d  %70s", natoms, title);
	setmolerrortext(str);
	if (*natoms == 0) return setmolerrorno(MERR_FORMAT);
	return 1;
}

static int	read_mmodel_header (FILE *fp, int *natoms, char title[])
{
	char	str[256];
	if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
	return get_mmodel_header(str, natoms, title);
}

MolPtr	FLoadFullMModel (FILE *fp, int natoms, char *title)
{
	MolPtr	m=NULL;
	Chain	*chain;
	Residue	*res;
	AtomPtr	atom, conn_atom;
	BondPtr	bond;
	int	i, n, a;
	char	*c, str[256];

	clearmolerror();
	if (!fp || natoms <= 0) return NULL;

	if (!(m = NewMol())) return NULL;
	m->header = LookUpStrTable(title);
	m->param_type = PARAM_MMODEL;

	m->chain = chain = NULL;
	m->bond = NULL;
	res = NULL;
	conn_atom = NULL;

	for(i=0;i<natoms;i++) {
		if (!READLINE(fp, str)) break;
		if ((c = strrchr(str, '\n'))) *c = '\0';
		if (!(atom = readmmodelline(str, &m->chain, &chain, &res, &conn_atom, &m->bond, i+1))) break;
	}
	natoms = i;

	if (!m->chain || getmolerrorno() != -1) {
		FreeMol(&m);
		setmolerrortext(str);
		return (MolPtr)setmolerrorno(MERR_FORMAT);
	}

	ForEachChainResAtom(m->chain, chain, res, atom) {
		int	neighbors[MAXNEIGHBOR+1], bondtypes[MAXNEIGHBOR+1];

		for(i=0,n=0;i<atom->nneighbors;i++) {
			a = (int)atom->nbatom[i];
			if (a <= 0 || a > natoms || a == atom->serno) continue;
			neighbors[n] = a;
			bondtypes[n] = (int)atom->nbbond[i];
			n++;
		}
		atom->nneighbors = n;
		for(i=0;i<n;i++) {
			atom->nbatom[i] = (AtomPtr)neighbors[i];
			atom->nbbond[i] = (BondPtr)bondtypes[i];
		}
	}

	m->bond = NULL;
	ForEachChainResAtom(m->chain, chain, res, atom) {
		for(i=0;i<atom->nneighbors;i++) {
			a = (int)atom->nbatom[i];
			/* Check if bond already exists */
			if (GetSernoBond(m->bond, atom->serno, a)) continue;
			if (!(bond = EnterNewBond(&m->bond))) {
				setmolerrorno(MERR_MEM);
				break;
			}
			bond->type = (int)atom->nbbond[i];
			bond->atom1 = atom;
			bond->atom2 = GetSernoAtom(m->chain, a);
		}
		if (getmolerrorno() != -1) break;
	}
	/* clear nbatom and nbbond fields --- they will be filled when GetNeighbor is called */
	ForEachChainResAtom(m->chain, chain, res, atom) {
		for(i=0;i<MAXNEIGHBOR;i++) {
			atom->nbatom[i] = NULL;
			atom->nbbond[i] = NULL;
		}
	}

	if (getmolerrorno() != -1) {
		FreeMol(&m);
		return NULL;
	}
	SetMolBBox3(m);

	return m;
}

MolPtr	FLoadMModel (FILE *fp)
{
	MolPtr	mlist=NULL, m=NULL, orig_mol=NULL;
	AtomPtr	atom;
	int	n, nf, natoms, orig_natoms, nmols, got_header;
	char	str[256], title[256];
	double	x, y, z;

	nmols = 0;
	clearmolerror();
	if (!fp) return NULL;

	if (!read_mmodel_header(fp, &natoms, title)) return NULL;
	if (natoms <= 0) return (MolPtr)setmolerrorno(MERR_FORMAT);
	clearmolerror();

	if (!(mlist = FLoadFullMModel(fp, natoms, title))) return NULL;
	orig_natoms = CountAtomInChains(mlist->chain);
	orig_mol = mlist;
	nmols++;

	/* check if there are more molecules */
	got_header = 0;
	while (1) {
		if (!got_header && !read_mmodel_header(fp, &natoms, title)) {
			clearmolerror();	/* could be EOF --- this is not an error */
			break;
		}
		got_header = 0;
		clearmolerror();

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "Number of Atoms = %d, Title = %s\n", natoms, title);
}

		if (natoms > 0) {
			if (!(m = FLoadFullMModel(fp, natoms, title))) {
				catmolerror("Failed to read the connection table No. %d.", nmols+1);
				break;
			}
			nmols++;
			orig_natoms = CountAtomInChains(m->chain);
			orig_mol = m;
			EnterMol(m, &mlist);
		} else if (natoms < 0 && -natoms == orig_natoms) {
			natoms = -natoms;

			/* NOTE: CopyMol copies a single molecule */
			if (!(m = CopyMol(orig_mol))) break;
			
			/* read a compressed connection table */
			while (READLINE(fp, str)) {
				nf = sscanf(str, "%d %lf %lf %lf", &n, &x, &y, &z);
				if (n <= 0 || n > natoms) {
					get_mmodel_header(str, &natoms, title);
					got_header = 1;
					break;
				}
				if (nf != 4 || !(atom = SearchAtomInChainByNumber(m->chain,n))) {

if (getmoldebuglevel() > 0) {
					setmolerrortext(str);
					setmolerrorno(MERR_FORMAT);
					printmolerror();
}

					FreeMol(&m);
					return mlist;
				}
				atom->x = x;
				atom->y = y;
				atom->z = z;
				if (n == natoms) break;
			}

			m->header = LookUpStrTable(title);
			EnterMol(m, &mlist);
		} else {
			setmolerrortext("Number of atoms in compressed connection table is different from the parent.");
			setmolerrorno(MERR_FORMAT);
			break;
		}
	}
	return mlist;
}

void	FPrintOneMModel (FILE *fp, MolPtr m)
{
	Chain	*c;
	Residue	*r;
	AtomPtr	atom, nbatom;
	BondPtr	nbbond;
	int	i, natoms, nconn, conn[MMOD_MAX_BONDS], type[MMOD_MAX_BONDS];
	char	chain_id, *title;

	if (!m) return;
	natoms = CountAtomInChains(m->chain);

	/* title */
	title = GetStrValue(m->header);
	fprintf(fp, "%6d  %s\n", natoms, title ? title : "");

	ForEachChainResAtom(m->chain, c, r, atom) {
		if (m->param_type != PARAM_MMODEL) {
			fprintf(fp, " %3d", GetMModelAtomType(atom));
		} else fprintf(fp, " %3d", atom->type);

		for(i=0;i<MMOD_MAX_BONDS;i++) conn[i] = type[i] = 0;
		for(i=0,nconn=0;i<atom->nneighbors && nconn < MMOD_MAX_BONDS;i++) {
			nbatom = atom->nbatom[i];
			nbbond = atom->nbbond[i];
			conn[nconn] = nbatom->serno;
			type[nconn] = GetMModelBondType(nbbond->type);
			nconn++;
		}

		for(i=0;i<MMOD_MAX_BONDS;i++) fprintf(fp, " %5d %d", conn[i], type[i]);
		fprintf(fp, " % 11.6f % 11.6f % 11.6f", atom->x, atom->y, atom->z);

		/* optional part */
		if (c->id == '\0') chain_id = ' ';
		else chain_id = c->id;
		fprintf(fp, " %5d%c%c%4d% 9.5f% 9.5f %4s %4s",
				0, '0', chain_id, 0, atom->charge, atom->charge, GetResName(r->refno), GetAtomName(atom->refno));

		fprintf(fp, "\n");
	}
}

void	FPrintMModel (FILE *fp, MolPtr mlist)
{
	MolPtr	m;

	if (!mlist) return;
	ForEachMol(mlist,m) FPrintOneMModel(fp, m);
}

