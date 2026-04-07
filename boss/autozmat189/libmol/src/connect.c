#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

/***************** GENERATION OF CONNECTIVITY ******************************/
ListPtr	LookUpBondList (ListPtr bondlist, int a1, int a2)
{
	ListPtr	list;
	ForEachList(bondlist,list) {
		if (list->L_ATOM1 == a1 && list->L_ATOM2 == a2) return list;
		if (list->L_ATOM1 == a2 && list->L_ATOM2 == a1) return list;
	}
	return NULL;
}

ListPtr	LookUpAngleList (ListPtr anglelist, int a1, int a2, int a3)
{
	ListPtr	list;
	ForEachList(anglelist,list) {
		if (list->L_ATOM1 == a1 && list->L_ATOM2 == a2 && list->L_ATOM3 == a3) return list;
		if (list->L_ATOM1 == a3 && list->L_ATOM2 == a2 && list->L_ATOM3 == a1) return list;
	}
	return NULL;
}

ListPtr	LookUpDihedList (ListPtr dihedlist, int a1, int a2, int a3, int a4)
{
	ListPtr	list;

	ForEachList(dihedlist,list) {
		if ((list->L_ATOM1 == a1 && list->L_ATOM2 == a2 && list->L_ATOM3 == a3 && list->L_ATOM4 == a4) ||
		    (list->L_ATOM1 == a4 && list->L_ATOM2 == a3 && list->L_ATOM3 == a2 && list->L_ATOM4 == a1)) {
			return list;
		}
	}
	return NULL;
}

static int	lookup_imp (int ia, int ja, int ka, int la, int a1, int a2, int a3, int a4)
{
	if ((a1 == ia && a2 == ja && a3 == ka && a4 == la) ||
	    (a1 == ia && a2 == ja && a3 == la && a4 == ka) ||
	    (a1 == ka && a2 == ja && a3 == ia && a4 == la) ||
	    (a1 == ka && a2 == ja && a3 == la && a4 == ia) ||
	    (a1 == la && a2 == ja && a3 == ia && a4 == ka) ||
	    (a1 == la && a2 == ja && a3 == ka && a4 == ia)) return 1;
	else return 0;
}

ListPtr	LookUpImproperList (ListPtr improperlist, ListPtr dihedlist, int a1, int a2, int a3, int a4)
{
	ListPtr	list, implist;
	int	ia, ja, ka, la;

	ForEachList(improperlist,list) {
		ia = list->L_ATOM1;
		ja = list->L_ATOM2;
		ka = list->L_ATOM3;
		la = list->L_ATOM4;
		if ((implist = LookUpDihedList(dihedlist, ia, ja, ka, la))) return implist;
		if (lookup_imp(ia, ja, ka, la, a1, a2, a3, a4)) return list;
		if (lookup_imp(ia, ja, ka, la, a4, a3, a2, a1)) return list;
	}
	return NULL;
}

ListPtr	MakeBondList (ChainPtr clist, ZmatPtr zmatlist)
{
	ChainPtr	c;
	ResiduePtr	r;
	register AtomPtr	atom1, atom2;
	register int	i;
	ListPtr	bondlist, list;
	ZmatPtr	zmat1, zmat2;
	int	z1, z2;

	if (!clist) return NULL;

	bondlist = NULL;
	ForEachChainResAtom(clist,c,r,atom1) {
		if (atom1->an < 0) continue;
		zmat1 = SearchZmatByAtomSerno(zmatlist, atom1->serno);
		z1 = SearchZmatByPointer(zmatlist, zmat1);

		for(i=0;i<atom1->nneighbors;i++) {
			atom2 = atom1->nbatom[i];
			if (atom2->an < 0) continue;
			if (atom1 == atom2) continue;
			zmat2 = SearchZmatByAtomSerno(zmatlist, atom2->serno);
			z2 = SearchZmatByPointer(zmatlist, zmat2);

			if (LookUpBondList(bondlist, z1, z2)) continue;
			if (!(list = EnterNewList(&bondlist))) return bondlist;

			list->L_ATOM1 = z1;
			list->L_ATOM2 = z2;
			list->L_ATOMPTR1 = (size_t) atom1;
			list->L_ATOMPTR2 = (size_t) atom2;
		}
	}

	return bondlist;
}

ListPtr	MakeAngleList (ChainPtr clist, ZmatPtr zmatlist)
{
	ChainPtr	c;
	ResiduePtr	r;

	register AtomPtr	atom1, atom2, atom3;
	register int		i, j;
	ListPtr	anglist, list;

	ZmatPtr	zmat1, zmat2, zmat3;
	int	z1, z2, z3;

	if (!clist) return NULL;

	anglist = NULL;
	ForEachChainResAtom(clist,c,r,atom1) {
		if (atom1->an < 0) continue;
		zmat1 = SearchZmatByAtomSerno(zmatlist, atom1->serno);
		z1 = SearchZmatByPointer(zmatlist, zmat1);

		for(i=0;i<atom1->nneighbors;i++) {
			atom2 = atom1->nbatom[i];
			if (atom2->an < 0) continue;
			if (atom1 == atom2) continue;
			zmat2 = SearchZmatByAtomSerno(zmatlist, atom2->serno);
			z2 = SearchZmatByPointer(zmatlist, zmat2);

			for(j=0;j<atom2->nneighbors;j++) {
				atom3 = atom2->nbatom[j];
				if (atom3->an < 0) continue;
				if (atom3 == atom1 || atom3 == atom2) continue;
				zmat3 = SearchZmatByAtomSerno(zmatlist, atom3->serno);
				z3 = SearchZmatByPointer(zmatlist, zmat3);

				if (LookUpAngleList(anglist, z1, z2, z3)) continue;
				if (!(list = EnterNewList(&anglist))) return anglist;

				list->L_ATOM1 = z1;
				list->L_ATOM2 = z2;
				list->L_ATOM3 = z3;
				list->L_ATOMPTR1 = (size_t) atom1;
				list->L_ATOMPTR2 = (size_t) atom2;
				list->L_ATOMPTR3 = (size_t) atom3;
			}
		}
	}
	return anglist;
}

ListPtr	MakeDihedList (ChainPtr clist, ZmatPtr zmatlist)
{
	ChainPtr	c;
	ResiduePtr	r;
	register AtomPtr	atom1, atom2, atom3, atom4;
	register int		i, j, k;
	ListPtr	dihedlist, list;

	ZmatPtr	zmat1, zmat2, zmat3, zmat4;
	int	z1, z2, z3, z4;

	if (!clist) return NULL;

	dihedlist = NULL;
	ForEachChainResAtom(clist,c,r,atom1) {
		if (atom1->an < 0) continue;
		zmat1 = SearchZmatByAtomSerno(zmatlist, atom1->serno);
		z1 = SearchZmatByPointer(zmatlist, zmat1);

		for(i=0;i<atom1->nneighbors;i++) {
			atom2 = atom1->nbatom[i];
			if (atom2->an < 0) continue;
			if (atom1 == atom2) continue;
			zmat2 = SearchZmatByAtomSerno(zmatlist, atom2->serno);
			z2 = SearchZmatByPointer(zmatlist, zmat2);

			for(j=0;j<atom2->nneighbors;j++) {
				atom3 = atom2->nbatom[j];
				if (atom3->an < 0) continue;
				if (atom3 == atom1 || atom3 == atom2) continue;
				zmat3 = SearchZmatByAtomSerno(zmatlist, atom3->serno);
				z3 = SearchZmatByPointer(zmatlist, zmat3);

				for(k=0;k<atom3->nneighbors;k++) {
					atom4 = atom3->nbatom[k];
					if (atom4->an < 0) continue;
					if (atom4 == atom1 || atom4 == atom2 || atom4 == atom3) continue;
					zmat4 = SearchZmatByAtomSerno(zmatlist, atom4->serno);
					z4 = SearchZmatByPointer(zmatlist, zmat4);

					if (LookUpDihedList(dihedlist, z1, z2, z3, z4)) continue;
					if (!(list = EnterNewList(&dihedlist))) return dihedlist;

					list->L_ATOM1 = z1;
					list->L_ATOM2 = z2;
					list->L_ATOM3 = z3;
					list->L_ATOM4 = z4;
					list->L_ATOMPTR1 = (size_t) atom1;
					list->L_ATOMPTR2 = (size_t) atom2;
					list->L_ATOMPTR3 = (size_t) atom3;
					list->L_ATOMPTR4 = (size_t) atom4;
				}
			}
		}
	}

	return dihedlist;
}

ListPtr	MakeImproperList (ChainPtr clist, ZmatPtr zmatlist, ListPtr dihedlist)
{
	ChainPtr	c;
	ResiduePtr	r;
	register AtomPtr	atom1, atom2, atom3, atom4;
	register int		i, k;
	ListPtr	improperlist, list;

	ZmatPtr	zmat1, zmat2, zmat3, zmat4;
	int	z1, z2, z3, z4;

	if (!clist) return NULL;

	improperlist = NULL;
	ForEachChainResAtom(clist,c,r,atom1) {
		if (atom1->an < 0) continue;
		zmat1 = SearchZmatByAtomSerno(zmatlist, atom1->serno);
		z1 = SearchZmatByPointer(zmatlist, zmat1);

		for(i=0;i<atom1->nneighbors;i++) {
			atom2 = atom1->nbatom[i];
			if (atom2->an < 0) continue;
			if (atom1 == atom2) continue;
			zmat2 = SearchZmatByAtomSerno(zmatlist, atom2->serno);
			z2 = SearchZmatByPointer(zmatlist, zmat2);

			atom3 = atom2->nbatom[0];
			if (atom3->an < 0) continue;
			if (atom3 == atom1 || atom3 == atom2) continue;
			zmat3 = SearchZmatByAtomSerno(zmatlist, atom3->serno);
			z3 = SearchZmatByPointer(zmatlist, zmat3);

			for(k=1;k<atom2->nneighbors;k++) {
				atom4 = atom2->nbatom[k];
				if (atom4->an < 0) continue;
				if (atom4 == atom1 || atom4 == atom2 || atom4 == atom3) continue;
				zmat4 = SearchZmatByAtomSerno(zmatlist, atom4->serno);
				z4 = SearchZmatByPointer(zmatlist, zmat4);

				if (LookUpImproperList(improperlist, dihedlist, z1, z2, z3, z4)) continue;
				if (!(list = EnterNewList(&improperlist))) return improperlist;

				list->L_ATOM1 = z1;
				list->L_ATOM2 = z2;
				list->L_ATOM3 = z3;
				list->L_ATOM4 = z4;
				list->L_ATOMPTR1 = (size_t) atom1;
				list->L_ATOMPTR2 = (size_t) atom2;
				list->L_ATOMPTR3 = (size_t) atom3;
				list->L_ATOMPTR4 = (size_t) atom4;
			}
		}
	}

	return improperlist;
}

/* assume 'clist' has correct neighbor information
(atom->nbatom, atom->nbbond, atom->nneighbors).
Based on this info, construct bond, angle, dihedral and improper torsion lists. */

int	MakeConnect (ChainPtr clist, ZmatPtr zmatlist, ListPtr *listbnd, ListPtr *listang, ListPtr *listdih, ListPtr *listimp)
{
	if (!clist || !zmatlist) return 0;
	if (listbnd) *listbnd = MakeBondList(clist, zmatlist);
	if (listang) *listang = MakeAngleList(clist, zmatlist);
	if (listdih) *listdih = MakeDihedList(clist, zmatlist);
	if (listimp) *listimp = MakeImproperList(clist, zmatlist, *listdih);

	return 1;
}

/*******************************************************************
* functions for getting BOSS addtional bonds, angles and dihedrals *
********************************************************************/

static int	iszmatbond (ZmatPtr zmatlist, int z1, int z2)
{
	ZmatPtr	zmat;
	int	a1, a2;

	a1 = 0;
	ForEachZmat(zmatlist,zmat) {
		a1++;
		a2 = zmat->zdef[0];
		if ((a1 == z1 && a2 == z2) || (a1 == z2 && a2 == z1)) return 1;
	}
	return 0;
}

static int	iszmatangle (ZmatPtr zmatlist, int z1, int z2, int z3)
{
	ZmatPtr	zmat;
	int	a1, a2, a3;

	a1 = 0;
	ForEachZmat(zmatlist,zmat) {
		a1++;
		a2 = zmat->zdef[0];
		a3 = zmat->zdef[1];
		if ((a1 == z1 && a2 == z2 && a3 == z3) || (a1 == z3 && a2 == z2 && a3 == z1)) return 1;
	}
	return 0;
}

static int	iszmatdihedral (ZmatPtr zmatlist, int z1, int z2, int z3, int z4)
{
	ZmatPtr	zmat;
	int	a1, a2, a3, a4;

	a1 = 0;
	ForEachZmat(zmatlist,zmat) {
		a1++;
		a2 = zmat->zdef[0];
		a3 = zmat->zdef[1];
		a4 = zmat->zdef[2];
		if ((a1 == z1 && a2 == z2 && a3 == z3 && a4 == z4) ||
		    (a1 == z4 && a2 == z3 && a3 == z2 && a4 == z1)) return 1;

	}
	return 0;
}

static int	iszmatimproper (ZmatPtr zmatlist, int z1, int z2, int z3, int z4)
{
	ZmatPtr	zmat;
	int	a1, a2, a3, a4;

	a1 = 0;
	ForEachZmat(zmatlist,zmat) {
		a1++;
		a2 = zmat->zdef[0];
		a3 = zmat->zdef[1];
		a4 = zmat->zdef[2];
		if (lookup_imp(z1, z2, z3, z4, a1, a2, a3, a4)) return 1;
		if (lookup_imp(z1, z2, z3, z4, a4, a3, a2, a1)) return 1;
	}
	return 0;
}

ListPtr	GetAddBondList (ChainPtr clist, ZmatPtr zmatlist)
{
	ListPtr	list, next, bondlist;

	if (!(bondlist = MakeBondList(clist, zmatlist))) return NULL;
	for(list=bondlist;list;) {
		next = list->next;
		if (iszmatbond(zmatlist, list->L_ATOM1, list->L_ATOM2)) {
			DeleteList(list, &bondlist);
		}
		list = next;
	}
	return bondlist;
}

ListPtr	GetAddAngleList (ChainPtr clist, ZmatPtr zmatlist)
{
	ListPtr	list, next, anglelist;

	if (!(anglelist = MakeAngleList(clist, zmatlist))) return NULL;
	for(list=anglelist;list;) {
		next = list->next;
		if (iszmatangle(zmatlist, list->L_ATOM1, list->L_ATOM2, list->L_ATOM3)) {
			DeleteList(list, &anglelist);
		}
		list = next;
	}
	return anglelist;
}

ListPtr	GetAddDihedList (ChainPtr clist, ZmatPtr zmatlist)
{
	ListPtr	list, next, dihedlist, improperlist;

	dihedlist = MakeDihedList(clist, zmatlist);
	improperlist = MakeImproperList(clist, zmatlist, dihedlist);
	if (!dihedlist && !improperlist) return NULL;

	for(list=dihedlist;list;) {
		next = list->next;
		list->L_IMPROPER = 0;
		if (iszmatdihedral(zmatlist, list->L_ATOM1, list->L_ATOM2, list->L_ATOM3, list->L_ATOM4)) {
			DeleteList(list, &dihedlist);
		}
		list = next;
	}
	for(list=improperlist;list;) {
		next = list->next;
		list->L_IMPROPER = 1;
		if (iszmatdihedral(zmatlist, list->L_ATOM1, list->L_ATOM2, list->L_ATOM3, list->L_ATOM4) ||
		    (iszmatimproper(zmatlist, list->L_ATOM1, list->L_ATOM2, list->L_ATOM3, list->L_ATOM4))) {
			DeleteList(list, &improperlist);
		}
		list = next;
	}
	EnterList(improperlist, &dihedlist);
	return dihedlist;
}


