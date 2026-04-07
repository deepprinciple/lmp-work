#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

#define GetElemNumberFromType(t)	((SearchAtomTypeTabByType(t))->an)
#define GetElemNumberFromString(s)	((SearchAtomTypeTabByString(s))->an)
#define GetMolTypeFromAtom(a)		((SearchAtomTypeTabByAtom(a))->mol)
#define GetMol2TypeFromAtom(a)		((SearchAtomTypeTabByAtom(a))->mol2)

typedef struct {
	int	id;
	int	an;
	char	mol2[6];
	int	mol;
	char	*comment;
	} AtomTypeTab;

enum	{
	MOL_DUM, MOL_C3, MOL_C2, MOL_CAR, MOL_C1,
	MOL_N4, MOL_N3, MOL_N2, MOL_N1, MOL_NAR, MOL_NAM, MOL_NPL3,
	MOL_O3, MOL_O2, MOL_OCO2, MOL_S3, MOL_S2, MOL_SO, MOL_SO2,
	MOL_P3, MOL_H, MOL_BR, MOL_CL, MOL_F, MOL_I, MOL_LP,
	MOL_SI, MOL_NA, MOL_K, MOL_CA, MOL_LI, MOL_AL,
	MOL_TYPE_MAX};

/* Sybyl MOL and MOL2 atom types */
static AtomTypeTab	atom_type_tab[] = {
	{MOL_DUM, -1, "Du   ",  26, "Dummy"},

	{MOL_C3,   6, "C.3  ",   1, "Csp3"},
	{MOL_C2,   6, "C.2  ",   2, "Csp2"},
	{MOL_CAR,  6, "C.ar ",   3, "C aromatic"},
	{MOL_C1,   6, "C.1  ",   4, "Csp"},

	{MOL_N4,   7, "N.4  ",  31, "Nsp3 +charged"},
	{MOL_N3,   7, "N.3  ",   5, "Nsp3"},
	{MOL_N2,   7, "N.2  ",   6, "Nsp2"},
	{MOL_N1,   7, "N.1  ",   7, "Nsp"},
	{MOL_NAR,  7, "N.ar ",  11, "N aromatic"},
	{MOL_NAM,  7, "N.am ",  28, "N amide"},
	{MOL_NPL3, 7, "N.pl3",  19, "N trigonal planar"},

	{MOL_O3,   8, "O.3  ",   8, "Osp3"},
	{MOL_O2,   8, "O.2  ",   9, "Osp2"},
	{MOL_OCO2, 8, "O.CO2",   9, "O in CO2"},

	{MOL_S3,  16, "S.3  ",  10, "Ssp3"},
	{MOL_S2,  16, "S.2  ",  18, "Ssp2"},
	{MOL_SO,  16, "S.o  ",  29, "S sulfoxide"},
	{MOL_SO2, 16, "S.o2 ",  30, "S sulfone"},

	{MOL_P3,  15, "P.3  ",  12, "P"},
	{MOL_H,    1, "H    ",  13, "H"},
	{MOL_BR,  35, "Br   ",  14, "Br"},
	{MOL_CL,  17, "Cl   ",  15, "Cl"},
	{MOL_F,    9, "F    ",  16, "F"},
	{MOL_I,   53, "I    ",  17, "I"},
	{MOL_LP,  98, "LP   ",  20, "Lone pair"},
	{MOL_SI,  14, "Si   ",  27, "Si"},
	{MOL_NA,  11, "Na   ",  21, "Na"},
	{MOL_K,   19, "K    ",  22, "K"},
	{MOL_CA,  20, "Ca   ",  23, "Ca"},
	{MOL_LI,   3, "Li   ",  24, "Li"},
	{MOL_AL,  13, "Al   ",  25, "Al"}
	};

#define ATOM_TYPE_SIZE	(sizeof(atom_type_tab)/sizeof(AtomTypeTab))

static AtomTypeTab	*SearchAtomTypeTabByID (int id)
{
	int	i;
	AtomTypeTab	*p;

	for(i=0,p=atom_type_tab;i<ATOM_TYPE_SIZE;i++,p++) {
		if (p->id == id) return p;
	}
	return atom_type_tab;
}

static AtomTypeTab	*SearchAtomTypeTabByType (int type)
{
	int	i;
	AtomTypeTab	*p;

	for(i=0,p=atom_type_tab;i<ATOM_TYPE_SIZE;i++,p++) {
		if (p->mol == type) return p;
	}
	return atom_type_tab;
}

/* an = atomic number
   nb = number of neighbors
 */
static AtomTypeTab	*SearchAtomTypeTabByElemNumber (int an, int nb)
{
	int	i, id=(-1);
	AtomTypeTab	*p;

	switch (an) {
	case 1: id = MOL_H;  break;
	case 3: id = MOL_LI; break;
	case 6:
		switch (nb) {
		case 2: id = MOL_C1; break;
		case 3: id = MOL_C2; break;
		case 4: id = MOL_C3; break;
		default: id = MOL_C3;
		}
		break;
	case 7:
		switch (nb) {
		case 1: id = MOL_N1; break;
		case 2: id = MOL_N2; break;
		case 3: id = MOL_N3; break;
		case 4: id = MOL_N4; break;
		default: id = MOL_N3; break;
		}
		break;
	case 8:
		switch (nb) {
		case 1: id = MOL_O2; break;
		case 2: id = MOL_O3; break;
		default: id = MOL_O3;
		}
		break;
	case 16:
		switch (nb) {
		case 1: id = MOL_S2; break;
		case 2: id = MOL_S3; break;
		default: id = MOL_S3; break;
		}
		break;

	case  9: id = MOL_F;  break;
	case 11: id = MOL_NA; break;
	case 13: id = MOL_AL; break;
	case 14: id = MOL_SI; break;
	case 15: id = MOL_P3; break;
	case 17: id = MOL_CL; break;
	case 19: id = MOL_K;  break;
	case 20: id = MOL_CA; break;

	case 35: id = MOL_BR; break;
	case 53: id = MOL_I;  break;
	}

	if (id == -1) {
		for(i=0,p=atom_type_tab;i<ATOM_TYPE_SIZE;i++,p++) {
			if (p->an == an) id = p->id;
		}
	}

	return SearchAtomTypeTabByID(id);
}

static AtomTypeTab	*SearchAtomTypeTabByAtom (AtomPtr atom)
{
	int	i, hybrid, an, nb, id=(-1);
	AtomTypeTab	*p;

	if (!atom) return atom_type_tab;
	hybrid = GetHybridization(atom);
	an = atom->an;
	nb = atom->nneighbors;

	switch (an) {
	case 1: id = MOL_H;  break;
	case 3: id = MOL_LI; break;
	case 6:
		switch (hybrid) {
		case HYBRID_SP:  id = MOL_C1; break;
		case HYBRID_SP2: id = MOL_C2; break;
		case HYBRID_SP3: id = MOL_C3; break;
		default: id = MOL_C3;
		}
		break;
	case 7:
		switch (hybrid) {
		case HYBRID_SP:       id = MOL_N1; break;
		case HYBRID_SP2:      id = MOL_N2; break;
		case HYBRID_SP3:      id = MOL_N3; break;
		case HYBRID_SP3_PLUS: id = MOL_N4; break;
		default: id = MOL_N3; break;
		}
		break;
	case 8:
		switch (hybrid) {
		case HYBRID_SP2: id = MOL_O2; break;
		case HYBRID_SP3: id = MOL_O3; break;
		default: id = MOL_O3;
		}
		break;
	case 16:
		switch (hybrid) {
		case HYBRID_SP2: id = MOL_S2; break;
		case HYBRID_SP3: id = MOL_S3; break;
		default: id = MOL_S3; break;
		}
		break;

	case  9: id = MOL_F;  break;
	case 11: id = MOL_NA; break;
	case 13: id = MOL_AL; break;
	case 14: id = MOL_SI; break;
	case 15: id = MOL_P3; break;
	case 17: id = MOL_CL; break;
	case 19: id = MOL_K;  break;
	case 20: id = MOL_CA; break;

	case 35: id = MOL_BR; break;
	case 53: id = MOL_I;  break;
	}

	if (id == -1) {
		for(i=0,p=atom_type_tab;i<ATOM_TYPE_SIZE;i++,p++) {
			if (p->an == an) id = p->id;
		}
	}

	return SearchAtomTypeTabByID(id);
}

static AtomTypeTab	*SearchAtomTypeTabByString (char *type)
{
	int	i, n;
	char	s[6];
	AtomTypeTab	*p;

	if (!type || !*type) return atom_type_tab;
	strncpy(s, type, 5);
	s[5] = '\0';
	TrimStr(s, " \t");
	n = strlen(s);
	for(i=n;i<5;i++) s[i] = ' ';

	for(i=0,p=atom_type_tab;i<ATOM_TYPE_SIZE;i++,p++) {
		if (strcasecmp(p->mol2, s) == 0) return p;
	}
	return atom_type_tab;
}

/*************************************************
 *****    Load Tripos MOL format file       ******
 *************************************************/

MolPtr	FLoadTriposMol (FILE *fp)
{
	char	str[256], title[81], sym[5], s1[16], s2[16], s3[16], s4[16];
	int	i, serno, ierr=0, iend=0, natoms, nbonds, ia, ja, type;
	ChainPtr	chain;
	ResiduePtr	res;
	AtomPtr	atom, atom1, atom2, conn_atom;
	BondPtr	bond;
	MolPtr	m=NULL;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	clearmolerror();
	if (!fp) return NULL;

	if (!READLINE(fp, str)) return NULL;	/* natoms */
	s1[0] = title[0] = '\0';
	strscan(str, "%4c %3c%10 %40c", s1, s2, title);
	if (sscanf(s1, "%d", &natoms) != 1 || natoms <= 0) return NULL;

	/* create a molecule */
	if (!(m = NewMol())) return NULL;
	m->header = LookUpStrTable(title);

	/* initialization */
	chain = NULL;
	res = NULL;
	conn_atom = NULL;
	strcpy(patom->name, " ");
	strcpy(pres->name, "UNK");
	pres->seq = 0;
	ierr = 0;
	iend = 0;

	nbonds = 0;
	for(i=0,serno=0;i<natoms;i++) {
		if (!READLINE(fp, str)) {
			iend = 1;
			break;
		}
		if (strncasecmp(str+5, "MOL", 3) == 0) {
			strscan(str, "%4c", s1);
			if (sscanf(s1, "%d", &nbonds) != 1 || nbonds <= 0) iend = 1;
			break;
		}

		strscan(str, "%4 %4c%9c%9c%9c%4c", s1, s2, s3, s4, sym);

		if (!(IsInt(s1) && IsNumber(s2) && IsNumber(s3) && IsNumber(s4))) {
			ierr = 1;
			setmolerrortext(str);
			setmolerrorno(MERR_FORMAT);
			break;
		}

		serno++;
		patom->serial = serno;
		strcpy(patom->name, sym);
		sscanf(s1, "%d", &type);
		sscanf(s2, "%lf", &patom->x);
		sscanf(s3, "%lf", &patom->y);
		sscanf(s4, "%lf", &patom->z);

		atom = CreatePdbAtom(&m->chain, &chain, &res, &conn_atom, patom, &m->bond, 1);
		atom->an = GetElemNumberFromType(type);
		atom->type = type;
		atom->refno = GetAtomRefno(sym);
	}

	if (ierr) {
		FreeMol(&m);
		return NULL;
	}

	if (!iend && nbonds == 0) {
		if (READLINE(fp, str) && strncasecmp(str+5, "MOL", 3) == 0) {
			strscan(str, "%4c", s1);
			if (sscanf(s1, "%d", &nbonds) != 1 || nbonds <= 0) iend = 1;
		}
	}

	if (iend || nbonds <= 0) {
		GetNeighbor(m);
		CalcBondsWithinResidues(m);
		FixConnectionTable(m);
		return m;
	}

	/* read connectivity */
	for(i=0;i<nbonds;i++) {
		if (!READLINE(fp, str) || strncasecmp(str+5, "MOL", 3) == 0) {
			iend = 1;
			break;
		}

		strscan(str, "%4 %4c%4c%9 %4c", s1, s2, s3);

		if (!(IsInt(s1) && IsInt(s2))) {
			ierr = 1;
			setmolerrortext(str);
			setmolerrorno(MERR_FORMAT);
			break;
		}

		sscanf(s1, "%d", &ia);
		sscanf(s2, "%d", &ja);
		sscanf(s3, "%d", &type);

		if (ia <= 0 || ia > serno || ja <= 0 || ja > serno) continue;
		switch (type) {
		case 1: type = B_SINGLE; break;
		case 2: type = B_DOUBLE; break;
		case 3: type = B_TRIPLE; break;

		case 4: type = B_AROMA;  break;	/* amide */
		case 5: type = B_AROMA;  break;	/* aromatic */

		default:type = B_SINGLE; break;
		}

		atom1 = GetSernoAtom(m->chain, ia);
		atom2 = GetSernoAtom(m->chain, ja);
		if (!atom1 || !atom2) continue;

		bond = EnterNewBond(&m->bond);
		bond->atom1 = atom1;
		bond->atom2 = atom2;
		bond->type = type;
	}

	if (iend || !m->bond) {
		GetNeighbor(m);
		CalcBondsWithinResidues(m);
		FixConnectionTable(m);
	}

	return m;
}

/************************************************
 ******  Sybyl Mol2 Format
 ******
 ************************************************/

static int read_molecule_record (FILE *fp, char *title, int *natoms, int *nbonds, int *charge_type)
{
	char	*p, str[256];

	if (!fp) return 0;

	/* mol_name */
	if (!READLINE(fp, str)) return 0;
	strcpy(title, str);
	if ((p = strchr(title, '\n'))) *p = '\0';

	/* num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
	 *
	 * num_atoms (integer) = the number of atoms in the molecule. 
	 * num_bonds (integer) = the number of bonds in the molecule. 
	 * num_subst (integer) = the number of substructures in the molecule. 
	 * num_feat (integer) = the number of features in the molecule. 
	 * num_sets (integer) = the number of sets in the molecule. 
	 */
	if (!READLINE(fp, str)) return 0;
	sscanf(str, "%d %d", natoms, nbonds);

	/* mol_type (SMALL, PROTEIN, NUCLEIC_ACID, SACCHARIDE) */
	if (!READLINE(fp, str)) return 0;

	/* charge_type
	 *
	 * NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN,
	 * GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_CHARGES, USER_CHARGES 
	 */
	if (!READLINE(fp, str)) return 0;
	p = str + strspn(str, " \t");
	if (strncasecmp(p, "NO_CHARGES", 10) == 0) *charge_type = 0;
	else *charge_type = 1;

	return 1;
}

/*
@<tripos>atom

The data records associated with this RTI consist of a single data line. This data line contains all
the information necessary to reconstruct one atom contained within the molecule. The atom ID
numbers associated with the atoms in the molecule will be assigned sequentially when the mol2
file is read into SYBYL.

Format:

atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]

Fields:

    atom_id (integer) = the ID number of the atom at the time the mol2 file was created. This
    is provided for reference only and is not used when the mol2 file is read into SYBYL. 
    atom_name (string) = the name of the atom. 
    x (real) = the x coordinate of the atom. 
    y (real) = the y coordinate of the atom. 
    z (real) = the z coordinate of the atom. 
    atom_type (string) = the SYBYL atom type for the atom. 
    subst_id (integer) = the ID number of the substructure containing the atom. 
    subst_name (string) = the name of the substructure containing the atom. 
    charge (real) = the charge associated with the atom. 
    status_bit (string) = the internal SYBYL status bits associated with the atom. These
    should never be set by the user. Valid status bits are DSPMOD, TYPECOL, CAP,
    BACKBONE, DICT, ESSENTIAL, WATER and DIRECT. 

Example:

1 CA -0.149 0.299 0.000 C.3 1 ALA1 0.000 BACKBONE|DICT|DIRECT
1 CA -0.149 0.299 0.000 C.3

In the first example the atom has ID number 1. It is named CA and is located at (-0.149, 0.299,
0.000). Its atom type is C.3. It belongs to the substructure with ID 1 which is named ALA1. The
charge associated with the atom is 0.000 and the SYBYL status bits associated with the atom are
BACKBONE, DICT, and DIRECT. Example two is the minimal information necessary for MOL2
to create an atom. 

*/

static int	read_atom_record (FILE *fp, MolPtr m, int natoms, int charge_type)
{
	char	str[256], type[16], name[16], subst_name[64];
	double	x, y, z, charge;
	int	i, subst_id;

	ChainPtr	chain;
	ResiduePtr	res;
	AtomPtr	atom, conn_atom;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	clearmolerror();
	if (!fp || !m || natoms <= 0) return 0;

	/* initialization */
	chain = NULL;
	res = NULL;
	conn_atom = NULL;
	strcpy(patom->name, " ");
	strcpy(pres->name, "UNK");
	pres->seq = 0;

	for(i=0;i<natoms;i++) {
		if (!READLINE(fp, str)) return 0;
		charge = 0;
		if (sscanf(str, "%*s %s %lf %lf %lf %s %d %s %lf",
			name, &x, &y, &z, type, &subst_id, subst_name, &charge) < 5) return 0;

		pres->seq = subst_id;
		strncpy(pres->name, subst_name, PDB_RLEN);
		pres->name[PDB_RLEN] = '\0';

		patom->serial = i+1;
		strcpy(patom->name, name);

		patom->x = x;
		patom->y = y;
		patom->z = z;

		atom = CreatePdbAtom(&m->chain, &chain, &res, &conn_atom, patom, &m->bond, 1);
		atom->an = GetElemNumberFromString(type);
		if (charge_type != 0) atom->charge = charge;
		atom->refno = GetAtomRefno(name);
	}

	return 1;
}

static int	read_bond_record (FILE *fp, MolPtr m, int natoms, int nbonds)
{
	char	str[256], type[8];
	int	i, n, ia, ja;
	AtomPtr	atom1, atom2;
	BondPtr	bond;

	if (!fp || !m || !m->chain) return 0;

	for(i=0;i<nbonds;i++) {
		if (!READLINE(fp, str)) return 0;
		if (sscanf(str, "%*d %d %d %s", &ia, &ja, type) < 2 ||
			ia <= 0 || ja <= 0) continue;

		if (strncasecmp(type, "ar", 2) == 0 || strncasecmp(type, "am", 2) == 0) n = B_AROMA;
		else if (sscanf(type, "%d", &n) != 1 || n <= 0 || n > 3) n = 1;

		if (ia <= 0 || ia > natoms || ja <= 0 || ja > natoms) continue;

		atom1 = GetSernoAtom(m->chain, ia);
		atom2 = GetSernoAtom(m->chain, ja);
		if (!atom1 || !atom2) continue;

		bond = EnterNewBond(&m->bond);
		bond->atom1 = atom1;
		bond->atom2 = atom2;
		bond->type = n;
	}

	if (!m->bond) {
		GetNeighbor(m);
		CalcBondsWithinResidues(m);
		FixConnectionTable(m);
	}

	return 1;
}

MolPtr	FLoadTriposMol2 (FILE *fp)
{
	MolPtr	m, mnext, mlist=NULL;
	char	str[256], title[256], *p, *q;
	int	n, natoms, nbonds, charge_type;

	while (READLINE(fp, str)) {
		p = str + strspn(str, " \t");
		if (strncasecmp(p, "@<TRIPOS>", 9) != 0) continue;
		q = p + 9;
		p = q + strspn(q, " \t");

		if (strncasecmp(q, "MOLECULE", 8) == 0) {
			if (!read_molecule_record(fp, title, &natoms, &nbonds, &charge_type)) continue;
			m = EnterNewMol(&mlist);
			if (title[0]) m->header = LookUpStrTable(title);
			else m->header = LookUpStrTable("");
		}
		else if (strncasecmp(q, "ATOM", 4) == 0 && m && !m->chain) read_atom_record(fp, m, natoms, charge_type);
		else if (strncasecmp(q, "BOND", 4) == 0 && m && !m->bond)  read_bond_record(fp, m, natoms, nbonds);
	}

	/* remove molecules with no atoms */
	for(m=mlist,n=1;m;m=mnext,n++) {
		mnext = m->next;
		if (!m->chain) {
			catmolerror("Molecule #%d (title=%s) had a format error.\n", n, GetStrValue(m->header));
			DeleteMol(m, &mlist);
			m->prev = m->next = NULL;
			FreeMol(&m);
		}
	}
	if (getmolerrorno() != -1) printmolerror();

	return mlist;
}

/*****************************************************
 *****   Write molecule in Tripos MOL format   *******
 *****************************************************/

void	FPrintTriposMol (FILE *fp, MolPtr m)
{
	char	title[128], *p;
	int	n, ia, ja, type, natoms, nbonds;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	BondPtr	b;

	clearmolerror();
	if (!fp || !m) return;

	/* count number of atoms and bonds */
	natoms = CountAtomInChains(m->chain);

	nbonds = 0;
	ForEachBond(m->bond,b) nbonds++;

	if (natoms >= 10000) {
		setmolerror("Warning: The number of atoms %d exceeds the maximum 9999.", natoms);
		natoms = 9999;
	}
	if (nbonds >= 10000) {
		setmolerror("\nWarning: The number of bonds %d exceeds the maximum 9999.", nbonds);
		nbonds = 9999;
	}

	/* print header */
	if (m->header >= 0) p = GetStrValue(m->header);
	title[0] = '\0';
	if (p && *p) {
		strcpy(title, p);
		if ((p = strchr(title, '\n'))) *p = '\0';
	}
	if (!title[0]) strcpy(title, "Tripos Molecule");

	fprintf(fp, "%4d MOL          %-40s\n", natoms, title);

	/* print atom records */
	n = 0;
	ForEachChainResAtom(m->chain,c,r,a) {
		n++;
		type = GetMolTypeFromAtom(a);
		fprintf(fp, "%4d%4d% 9.4f% 9.4f% 9.4f%-4s\n",
			n, type, a->x, a->y, a->z, GetAtomName(a->refno));
	}

	/* print bond records */
	fprintf(fp, "%4d MOL\n", nbonds);
	n = 0;
	ForEachBond(m->bond,b) {
		n++;
		ia = SearchAtomInChainByPointer(m->chain, b->atom1);
		ja = SearchAtomInChainByPointer(m->chain, b->atom2);
		if (ia > natoms || ja > natoms) continue;

		switch (b->type) {
		case B_SINGLE: type = 1; break;
		case B_DOUBLE: type = 2; break;
		case B_TRIPLE: type = 3; break;
		case B_AROMA:  type = 5; break;
		default:       type = 1; break;
		}
		fprintf(fp, "%4d%4d%4d        %4d\n", n, ia, ja, type);
	}

	fprintf(fp, "%4d MOL\n", 0);
}

/*****************************************************
 *****   Write molecule in Tripos MOL2 format  *******
 *****************************************************/

void	FPrintTriposMol2 (FILE *fp, MolPtr m)
{
	char	title[128], *type, *p, *q;
	int	n, ia, ja, natoms, nbonds;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	BondPtr	b;

	clearmolerror();
	if (!fp || !m) return;

	/* count number of atoms and bonds */
	natoms = CountAtomInChains(m->chain);

	nbonds = 0;
	ForEachBond(m->bond,b) nbonds++;

	/* print molecule record */
	if (m->header >= 0) p = GetStrValue(m->header);
	else p = NULL;
	title[0] = '\0';
	if (p && *p) {
		strcpy(title, p);
		if ((p = strchr(title, '\n'))) *p = '\0';
	}
	if (!title[0]) strcpy(title, "Tripos Molecule");

	fprintf(fp, "@<TRIPOS>MOLECULE\n");
	fprintf(fp, "%s\n", title);
	fprintf(fp, "%d %d\n", natoms, nbonds);
	fprintf(fp, "SMALL\n");
	fprintf(fp, "USER_CHARGES\n");
	fprintf(fp, "\n\n");	/* two blank lines */

	/* print atom records */
	fprintf(fp, "@<TRIPOS>ATOM\n");
	n = 0;
	ForEachChainResAtom(m->chain,c,r,a) {
		n++;
		type = GetMol2TypeFromAtom(a);
		if (!(p = GetAtomName(a->refno))) p = GetElemSymbol(a->an);
		if (!(q = GetResName(a->residue->refno)) || strspn(q, " \t") == strlen(q)) q = "UNK";
		fprintf(fp, "%4d %5s % 10.5f % 10.5f % 10.5f %5s %4d %8s % 10.5f\n",
			n, p, a->x, a->y, a->z, type,
			a->residue->seqno, q, a->charge);
	}

	/* print bond records */
	fprintf(fp, "@<TRIPOS>BOND\n");
	n = 0;
	ForEachBond(m->bond,b) {
		n++;
		ia = SearchAtomInChainByPointer(m->chain, b->atom1);
		ja = SearchAtomInChainByPointer(m->chain, b->atom2);

		switch (b->type) {
		case B_SINGLE: type = "1";  break;
		case B_DOUBLE: type = "2";  break;
		case B_TRIPLE: type = "3";  break;
		case B_AROMA:  type = "AR"; break;
		default:       type = "1";  break;
		}
		fprintf(fp, "%4d %4d %4d %4s\n", n, ia, ja, type);
	}
}

