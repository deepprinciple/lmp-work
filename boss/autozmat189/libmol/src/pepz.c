#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

static ListPtr	dihed_type_list=NULL;

static void	fprint_pepz_dihedral_list (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist,
	ListPtr dihedlist, ListPtr atomtypelist, int add_dihed_flags)
{
	ListPtr	list, list2;
	ZmatPtr	z1, z2, z3, z4;
	char	*amber, *amber1, *amber2, *amber3, type[64];
	AtomPtr	atom1, atom2, atom3, atom4;
	int	improper, a1, a2, a3, a4;

	ForEachList(dihedlist,list) {
		if (!add_dihed_flags) {	/* variable dihedrals */
			a1 = (int)list->L_CENTER;
			z1 = SearchZmat(zmatlist,a1);
			if (!z1) continue;
			a2 = z1->zdef[0];
			a3 = z1->zdef[1];
			a4 = z1->zdef[2];
		} else {		/* additional dihedrals */
			a1 = list->L_ATOM1;
			a2 = list->L_ATOM2;
			a3 = list->L_ATOM3;
			a4 = list->L_ATOM4;
			z1 = SearchZmat(zmatlist,a1);
		}
		z2 = SearchZmat(zmatlist,a2);
		z3 = SearchZmat(zmatlist,a3);
		z4 = SearchZmat(zmatlist,a4);
		if (!z1 || !z2 || !z3 || !z4) continue;

		atom1 = SearchAtomInChainByNumber(chainlist, z1->atom);
		atom2 = SearchAtomInChainByNumber(chainlist, z2->atom);
		atom3 = SearchAtomInChainByNumber(chainlist, z3->atom);
		atom4 = SearchAtomInChainByNumber(chainlist, z4->atom);

		improper = 0;
		if (atom2 && atom4) {
			/* if atom2 and atom4 are bonded, it's an improper dihedral */
			if (GetBondBetweenAtoms(atom2, atom4)) improper = 1;
		}

		/* ambertype */
		amber = amber1 = amber2 = amber3 = NULL;
		ForEachList(atomtypelist,list2) {
			if (list2->L_TYPE == z1->itype)  amber = list2->L_AMBER;
			if (list2->L_TYPE == z2->itype) amber1 = list2->L_AMBER;
			if (list2->L_TYPE == z3->itype) amber2 = list2->L_AMBER;
			if (list2->L_TYPE == z4->itype) amber3 = list2->L_AMBER;
			if (amber && amber1 && amber2 && amber3) break;
		}
		strcpy(type, (improper || list->L_ITYPE == 0) ? "i":"2");

		ForEachList(dihed_type_list,list2) {
			if (   (list2->L_ATOM1 == z1->itype &&
				list2->L_ATOM2 == z2->itype &&
				list2->L_ATOM3 == z3->itype &&
				list2->L_ATOM4 == z4->itype) ||
			       (list2->L_ATOM1 == z4->itype &&
				list2->L_ATOM2 == z3->itype &&
				list2->L_ATOM3 == z2->itype &&
				list2->L_ATOM4 == z1->itype)) {
				if (strcmp(list2->L_COMMENT, type) == 0) break;
			}
		}

		if (!list2) {
			list2 = EnterNewList(&dihed_type_list);
			list2->L_ATOM1 = z1->itype;
			list2->L_ATOM2 = z2->itype;
			list2->L_ATOM3 = z3->itype;
			list2->L_ATOM4 = z4->itype;
			strcpy(list2->L_COMMENT, type);
		} else continue;
		fprintf(fp, "%4d-%4d-%4d-%4d | %-2s-%-2s-%-2s-%-2s | %04d |   %s |\n",
			z1->itype, z2->itype, z3->itype, z4->itype,
			amber  ? amber  : "",
			amber1 ? amber1 : "",
			amber2 ? amber2 : "",
			amber3 ? amber3 : "",
			list->L_ITYPE,
			type);
	}
}

static void	fprint_pepz_dihedral (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist,
	BossZmatPtr boss, ListPtr atomtypelist)
{
	ZmatPtr	z1,z2;
	ListPtr	list;

	if (!boss) return;

	ForEachList(boss->addbond,list) {
		z1 = SearchZmat(zmatlist,list->L_ATOM1);
		z2 = SearchZmat(zmatlist,list->L_ATOM2);
		if (!z1 || !z2) continue;
		fprintf(fp, "%-3s%5d %-3s%5d\n", z1->atomname[0], 0, z2->atomname[0], 0);
	}

	/* dihedral info */
	fprintf(fp, "#\n");
	fprintf(fp, "# List of Dihedrals: format(4(i4,1x),2x,4(a2,1x),2x,i4,2x,i4)\n");
	fprintf(fp, "#\n");
	fprintf(fp, "#   BOSS atom type  | AMBER Atype | dihe | 1,4 |\n");
	fprintf(fp, "# i    j    k    l  |  i  j  k  l |number|scale|\n");
	fprintf(fp, "#----------------------------------------------|\n");

	dihed_type_list = NULL;
	fprint_pepz_dihedral_list(fp, chainlist, zmatlist, boss->vardihed, atomtypelist, 0);
	fprint_pepz_dihedral_list(fp, chainlist, zmatlist, boss->adddihed, atomtypelist, 0);
	fprintf(fp, "\n");

	if (dihed_type_list) FreeList(&dihed_type_list);
}

static void	fprint_atomtypes (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist, ListPtr atomtypelist, int flags)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	atom;
	int	type, ntypes=0;
	ListPtr	list, newlist, tmp_atomtype_list=NULL;
	ZmatPtr	zmat;

	fprintf(fp, "#\n");
	fprintf(fp, "# List of OPLS Atom Types:\n");
	fprintf(fp, "#\n");
	fprintf(fp, " TYP AN AT   CHARGE     SIGMA    EPSILON\n");

	ForEachChainResAtom(chainlist,c,r,atom) {
		if (atom->an < 0) continue;	/* dummy atom */
		zmat = NULL;
		if ((flags & PEPZ_ATOMTYPE_FROM_ZMAT) &&
		    (zmat = SearchZmatByAtomSerno(zmatlist, atom->serno))) {
			type = zmat->itype;
		} else type = atom->type;

		if (type <= 0 || type > 9999) continue;
		if (SearchList(tmp_atomtype_list, type)) continue;

		if (!(list = SearchList(atomtypelist, type))) {
			newlist = EnterNewList(&tmp_atomtype_list);
			newlist->L_AN = atom->an == 0 ? 1 : atom->an;
			newlist->L_TYPE = type;
			newlist->L_CHARGE = atom->charge;
			newlist->L_SIGMA = atom->sigma;
			newlist->L_EPSILON = atom->epsilon;
			strcpy(newlist->L_AMBER, atom->ambertype);
		} else {
			newlist = NewList();
			BCOPY(list, newlist, sizeof(List));
			newlist->next = NULL;
			EnterList(newlist, &tmp_atomtype_list);
		}
		ntypes++;
	}

	if (!tmp_atomtype_list) return;

	/* sort the atom types in ascending order.
	 * M_DATA_1 = data1 = L_TYPE
	 */
	QsortList(tmp_atomtype_list, 1, ntypes, M_DATA_1);

	/* print out */
	ForEachList(tmp_atomtype_list, list) {
		fprintf(fp, "%04d %02d %-2s % 10.6f %9.6f %9.6f     %s\n",
			list->L_TYPE, list->L_AN, list->L_AMBER, list->L_CHARGE,
			list->L_SIGMA, list->L_EPSILON, list->L_COMMENT);
	}

	/* free the list */
	FreeList(&tmp_atomtype_list);
}

void	FPrintPepz (FILE *fp, ChainPtr chainlist, ZmatPtr zmatlist, BossZmatPtr boss, ListPtr atomtypelist, int flags)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	atom;
	ZmatPtr	z,z1,z2,z3;
	ListPtr	list;
	int	natoms=0, nzmat=0, naddbonds=0;
	char	tree, sym[8];

	if (!fp || !chainlist || !zmatlist) return;
	ForEachChainResAtom(chainlist,c,r,atom) natoms++;
	if (natoms == 0) return;

	/* set atom names */
	nzmat = 0;
	ForEachZmat(zmatlist,z) {
		nzmat++;
		strcpy(sym, (z->an<=0 || z->an>=98) ? "X" : GetElemSymbol(z->an));
		if (sym[1] == ' ') sym[1] = '\0';

		sprintf(z->atomname[0], "%s%d", sym, nzmat);
		z->atomname[0][3] = '\0';
		z1 = SearchZmat(zmatlist, z->zdef[0]);
		z2 = SearchZmat(zmatlist, z->zdef[1]);
		z3 = SearchZmat(zmatlist, z->zdef[2]);
		if (z1) strcpy(z->atomname[1], z1->atomname[0]);
		if (z2) strcpy(z->atomname[2], z2->atomname[0]);
		if (z3) strcpy(z->atomname[3], z3->atomname[0]);

		switch (nzmat) {
		case 1:
			strcpy(z->atomname[1], "DU3");
			strcpy(z->atomname[2], "DU2");
			strcpy(z->atomname[3], "DU1");
			break;
		case 2:
			strcpy(z->atomname[2], "DU3");
			strcpy(z->atomname[3], "DU2");
			break;
		case 3:
			strcpy(z->atomname[3], "DU3");
			break;
		}
	}

	/* count the number of additional bonds */
	naddbonds = 0;
	if (boss) {
		ForEachList(boss->addbond,list) naddbonds++;
	}

	/* print out first line in PEPZ */
	fprintf(fp, "RES%4d%4d%4d%4d%4d %-10s %3s %1s\n",
		natoms, naddbonds, 0, 0, 0, "OPLS", "ALL", "N");

	/* print out Z-matrix */
	nzmat = 0;
	ForEachZmat(zmatlist,z) {
		nzmat++;
		tree = 'M';

		if ((atom = SearchAtomInChainBySerno(chainlist,z->atom))) {
			if (nzmat > 3 && atom->nneighbors <= 1) tree = 'E';
		}

		fprintf(fp, "  %-3s  %03d %-3s %3d % 10.3f  %-3s %3d  % 10.2f %-3s %3d  % 10.1f %3s%c\n",
			z->atomname[0], z->itype,
			z->atomname[1], nzmat <= 1 ? -1 : 0, nzmat <= 1 ?   1.0 : z->zval[0],
			z->atomname[2], nzmat <= 2 ? -1 : 0, nzmat <= 2 ?  90.0 : z->zval[1],
			z->atomname[3], nzmat <= 3 ? -1 : 0, nzmat <= 3 ? 180.0 : z->zval[2],
			"   ", tree);
	}

	/* print out torsion types */
	if (boss) {
		fprint_pepz_dihedral(fp, chainlist, zmatlist, boss, atomtypelist);
	}

	/* print out atom types */
	fprint_atomtypes(fp, chainlist, zmatlist, atomtypelist, flags);
}

