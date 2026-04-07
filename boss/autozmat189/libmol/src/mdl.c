#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

MolPtr	FLoadMdlMol (FILE *fp)
{
	char	str[256], title[81], sym[4], s1[16], s2[16], s3[16], *p;
	int	i, an, serno, ierr=0, iend=0, natoms, nbonds, ia, ja, type;
	ChainPtr	chain;
	ResiduePtr	res;
	AtomPtr	atom, atom1, atom2, conn_atom;
	BondPtr	bond;
	MolPtr	m=NULL;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	clearmolerror();
	if (!fp) return NULL;

	if (!READLINE(fp, str)) return NULL;	/* title */
	sscanf(str, "%80c", title);
	title[80] = '\0';
	if ((p = strchr(title, '\n'))) *p = '\0';

	if (!READLINE(fp, str)) return NULL;	/* program details */
	if (!READLINE(fp, str)) return NULL;	/* comments */
	if (!READLINE(fp, str)) return NULL;	/* natoms and nbonds */
	strscan(str, "%3d%3d", &natoms, &nbonds);

	if (natoms <= 0) return NULL;
	if (natoms >= 1000) natoms = 999;

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

	for(i=0,serno=0;i<natoms;i++) {
		if (!READLINE(fp, str)) {
			iend = 1;
			break;
		}
		if (strncmp(str, "$$$$", 4) == 0 || strncmp(str, "M  END", 6) == 0) {
			iend = 1;
			break;
		}

		strscan(str, "%10c%10c%10c %3s", s1, s2, s3, sym);
		if (strcmp(sym, "A") == 0 || strcmp(sym, "Q") == 0) an = -1;
		else an = GetElemNumber(sym);

		if (!(IsNumber(s1) && IsNumber(s2) && IsNumber(s3))) {
			ierr = 1;
			setmolerrortext(str);
			setmolerrorno(MERR_FORMAT);
			break;
		}

		serno++;
		patom->serial = serno;
		strcpy(patom->name, AN2ATOMNAME(an));
		sscanf(s1, "%lf", &patom->x);
		sscanf(s2, "%lf", &patom->y);
		sscanf(s3, "%lf", &patom->z);
		atom = CreatePdbAtom(&m->chain, &chain, &res, &conn_atom, patom, &m->bond, 1);
		atom->an = an;
	}

	if (ierr) {
		FreeMol(&m);
		return NULL;
	}

	if (iend || nbonds <= 0) {
		CalcBondsWithinResidues(m);
		return m;
	}

	/* read connectivity */
	for(i=0;i<nbonds;i++) {
		if (!READLINE(fp, str)) {
			iend = 1;
			break;
		}
		if (strncmp(str, "$$$$", 4) == 0 || strncmp(str, "M  END", 6) == 0) {
			iend = 1;
			break;
		}

		strscan(str, "%3c%3c%3c", s1, s2, s3);

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
		case 4: type = B_AROMA;  break;
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
		CalcBondsWithinResidues(m);
	}

	return m;
}

MolPtr	FLoadMdlSd (FILE *fp)
{
	MolPtr	m, mlist=NULL;
	char	str[256];
	int	ieof=0;

	while ((m = FLoadMdlMol(fp))) {
		EnterMol(m, &mlist);
		while (1) {
			if (!READLINE(fp, str)) {
				ieof = 1;
				break;
			}
			if (strncmp(str, "$$$$", 4) == 0) break;
		}
		if (ieof) break;
	}

	return mlist;
}

void	FPrintMdlMol (FILE *fp, MolPtr m)
{
	char	title[128], *p;
	int	n, an, ia, ja, type, natoms, nbonds;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	BondPtr	b;

	clearmolerror();
	if (!fp || !m) return;

	/* count number of atoms and bonds */
	natoms = CountAtomInChains(m->chain);
	if (natoms == 0) return;
	if (natoms >= 1000) {
		setmolerror("Warning: The number of atoms %d exceeds the maximum 999.", natoms);
		natoms = 999;
	}

	nbonds = 0;
	ForEachBond(m->bond,b) nbonds++;
	if (nbonds >= 1000) {
		catmolerror("\nWarning: The number of bonds %d exceeds the maximum 999.", nbonds);
		nbonds = 999;
	}

	/* print header */
	if (m->header >= 0) p = GetStrValue(m->header);
	title[0] = '\0';
	if (p && *p) {
		strcpy(title, p);
		if ((p = strchr(title, '\n'))) *p = '\0';
	}
	if (!title[0]) strcpy(title, "MDL MOL");

	fprintf(fp, "%s\n", title);
	fprintf(fp, "WJ BOSS4.0          3D\n\n");
	fprintf(fp, "%3d%3d  0  0  0  0  0  0  0  0  1 V2000\n", natoms, nbonds);

	/* atom coordinates */
	n = 0;
	ForEachChainResAtom(m->chain,c,r,a) {
		n++;
		if (n > natoms) break;
		if (a->an == 0) an = 0;
		p = GetElemSymbol(a->an);
		if (strcmp(p, "X ") == 0) p = "Q ";
		fprintf(fp, "% 10.4f% 10.4f% 10.4f %2s  0  0  0  0  0  0  0  0  0  0  0  0\n",
			a->x, a->y, a->z, p);
	}

	/* connectivity */
	n = 0;
	ForEachBond(m->bond,b) {
		ia = SearchAtomInChainByPointer(m->chain, b->atom1);
		ja = SearchAtomInChainByPointer(m->chain, b->atom2);
		if (ia > natoms || ja > natoms) continue;

		n++;
		if (n > 999) break;
		switch (b->type) {
		case B_SINGLE: type = 1; break;
		case B_DOUBLE: type = 2; break;
		case B_TRIPLE: type = 3; break;
		case B_AROMA:  type = 4; break;
		default:       type = 1; break;
		}
		fprintf(fp, "%3d%3d%3d%3d\n", ia, ja, type, 0);
	}

	/* print footer */
	fprintf(fp, "M  END\n");
	fprintf(fp, "\n");
	fprintf(fp, "$$$$\n");
}

void	FPrintMdlSd (FILE *fp, MolPtr mlist)
{
	MolPtr	m;

	ForEachMol(mlist,m) {
		FPrintMdlMol(fp, m);
	}
}

