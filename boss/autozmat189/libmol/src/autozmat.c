#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

/* find the closest atom from the isolated atom */
static AtomPtr	find_closest_atom (ChainPtr clist, AtomPtr isolated_atom, int flags)
{
	double	mindist, dist;
	AtomPtr	minatom, a;
	Residue	*r;
	Chain	*c;

	mindist = -1.0E+6;
	ForEachChainResAtom(clist,c,r,a) {
		if (flags && !(a->flags & flags)) continue;
		if (isolated_atom == a) continue;
		dist = GetAtomSQDist(isolated_atom, a);
		if (mindist <= 0 || mindist > dist) {
			mindist = dist;
			minatom = a;
		}
	}
	return minatom;
}

/* return 1 if atom1 has a higher priority than atom2.
   else return 0.

   Priority Criteria:
    1. number of neighbor atoms (i.e., substituents)
    2. atomic number
 */
static int	has_higher_priority (AtomPtr atom1, AtomPtr atom2)
{
	int	i, n1, n2, an1, an2;
	AtomPtr	nbatom;

	n1 = atom1->nneighbors;
	n2 = atom2->nneighbors;
	an1 = atom1->an;
	an2 = atom2->an;

	if (n1 != n2) return (n1 > n2);		/* number of substituents */
	if (an1 != an2) return (an1 > an2);	/* atomic number */

	/* sum number of substituents and atomic numbers of neighbor atoms */
	for(i=0,n1=an1=0;i<atom1->nneighbors;i++) {
		nbatom = atom1->nbatom[i];
		n1 += nbatom->nneighbors;
		an1 += nbatom->an;
	}
	for(i=0,n2=an2=0;i<atom2->nneighbors;i++) {
		nbatom = atom2->nbatom[i];
		n2 += nbatom->nneighbors;
		an2 += nbatom->an;
	}

	if (n1 != n2) return (n1 > n2);		/* sum of number of substituents */
	if (an1 != an2) return (an1 > an2);	/* sum of atomic number */

	/* OK, give up */
	return 1;
}

static void	sort_neighbor (AtomPtr nbatoms[], BondPtr nbbonds[], int nneighbors)
{
	int	i, j;
	AtomPtr	atom;
	BondPtr	bond;

	/* bubble sort neighbor atoms according to a priority determined
	 * by preset criteria.
	 */
	for(i=0;i<nneighbors-1;i++) {
		for(j=i+1;j<nneighbors;j++) {
			if (!has_higher_priority(nbatoms[i], nbatoms[j])) {
				atom = nbatoms[i];
				nbatoms[i] = nbatoms[j];
				nbatoms[j] = atom;
				bond = nbbonds[i];
				nbbonds[i] = nbbonds[j];
				nbbonds[j] = bond;
			}
		}
	}
}

static int	search_zmat_flags (ZmatPtr zmatlist, AtomPtr atom, ZmatPtr *zmat_ret)
{
	unsigned long	a = (unsigned long)atom;
	int	n=0;
	ZmatPtr	zmat;

	if (zmat_ret) *zmat_ret = NULL;
	ForEachZmat(zmatlist,zmat) {
		n++;
		if (zmat->flags == a) {
			if (zmat_ret) *zmat_ret = zmat;
			return n;
		}
	}
	return 0;
}

static int	is_valid_zdef (int z0, int z1, int z2, int z3)
{
	return !(z0 <= 0 || (z0 == 2 && z1 <= 0) || (z0 == 3 && (z1 <= 0 || z2 <= 0))
		|| (z0 > 3 && (z1 <= 0 || z2 <= 0 || z3 <= 0)));
}

static int	pick_zmat_connection (ZmatPtr zmatlist, AtomPtr atom, int nzmat, int z1, int *z2, int *z3)
{
	int	i, z;
	ZmatPtr	zmat1, zmat2, zmat3;
	AtomPtr	atom1, atom2, atom3, nbatom, atom2_new;

	if (!zmatlist || !atom) return 0;

	zmat1 = SearchZmat(zmatlist,  z1);
	zmat2 = SearchZmat(zmatlist, *z2);
	zmat3 = SearchZmat(zmatlist, *z3);

	if (!zmat1 || !zmat2 || !zmat3) return 0;

	atom1 = (AtomPtr)zmat1->flags;
	atom2 = (AtomPtr)zmat2->flags;
	atom3 = (AtomPtr)zmat3->flags;

	if (!atom1 || !atom2 || !atom3) return 0;
	if (atom2->an > 0 && atom3->an > 0) return 0;

	/* CASE 1: atom--atom1--atom2(=dummy)--atom3 */

	if (atom2->an <= 0) {
		atom2_new = NULL;
		for(i=0;i<atom1->nneighbors;i++) {
			nbatom = atom1->nbatom[i];
			if (nbatom->an <= 0 || !(nbatom->flags & A_MARKED)
				|| nbatom == atom2 || nbatom == atom) {
				continue;
			} else if (!(z = search_zmat_flags(zmatlist, nbatom, NULL)) || z >= nzmat) {
				continue;
			} else if (!IsAtomLinear(atom, atom1, nbatom, 5.0)) {
				atom2_new = nbatom;
				break;
			}
		}
		if (!atom2_new) {
			/* atom--atom1--dummy is the only way to define Z-matrix */
			return 0;
		}

		/* update atom2, zmat2, and *z2 */
		atom2 = atom2_new;
		for(zmat2=zmatlist,(*z2)=1;zmat2;zmat2=zmat2->next,(*z2)++) {
			if (zmat2->flags == (unsigned long)atom2) break;
		}
		if (!zmat2) return 0;
	}

	/* now we have atom--atom1--atom2--??
	 * So, we need to find atom3.
	 */
	atom3 = NULL;
	zmat3 = NULL;
	*z3 = 0;

	for(i=0;i<atom2->nneighbors;i++) {
		nbatom = atom2->nbatom[i];
		if (nbatom->an <= 0 || !(nbatom->flags & A_MARKED) ||
			nbatom == atom1 || nbatom == atom) {
			continue;
		} else if (!(z = search_zmat_flags(zmatlist, nbatom, NULL)) || z >= nzmat) {
			continue;
		} else if (!IsAtomLinear(atom, atom1, nbatom, 5.0) &&
			!IsAtomLinear(atom1,atom2,nbatom,5.0)) {
			atom3 = nbatom;
			break;
		}
	}

	if (!atom3) {
		/* look for an improper torsion
		 * atom--atom1--atom2
		 *         |
		 *       atom3
		 */
		for(i=0;i<atom1->nneighbors;i++) {
			nbatom = atom1->nbatom[i];
			if (nbatom->an <= 0 || !(nbatom->flags & A_MARKED)
				|| nbatom == atom2 || nbatom == atom) {
				continue;
			} else if (!(z = search_zmat_flags(zmatlist, nbatom, NULL)) || z >= nzmat) {
				continue;
			} else if (!IsAtomLinear(atom, atom1, nbatom, 5.0) &&
				!IsAtomLinear(atom1,atom2,nbatom,5.0)) {
					atom3 = nbatom;
				break;
			}
		}
	}
	if (!atom3) return 0;	/* can't help any more */

	for(zmat3=zmatlist,(*z3)=1;zmat3;zmat3=zmat3->next,(*z3)++) {
		if (zmat3->flags == (unsigned long)atom3) break;
	}
	if (!zmat3) return 0;

	return 1;
}

static int	autozmat (ChainPtr clist, AtomPtr theatom, ZmatPtr thezmat, int z0, int *nzmat, ZmatPtr *zmatlist)
{
	int	i, j, z1, z2, z3, z1_new, z2_new, zref, nzmat_old;
	AtomPtr	atom, a1, a2;
	AtomPtr	nbatom[MAXNEIGHBOR];
	ZmatPtr	zmat, zmat1, zmat2;

	if (!zmatlist) return 0;
	if (!*zmatlist) {
		z0 = 1;
		*nzmat = 1;
		thezmat = NewZmatFromAtom(theatom, 0, 0, 0);
		thezmat->flags = (unsigned long)theatom;	/* store atom pointer temporarily */
		EnterZmat(thezmat, zmatlist);
		theatom->flags |= A_MARKED;

	}

	z1 = thezmat->zdef[0];
	z2 = thezmat->zdef[1];
	z3 = thezmat->zdef[2];

	if (!is_valid_zdef(z0, z1, z2, z3)) {
		return setmolerror("autozmat ERROR: z = %d %d %d %d", z0, z1, z2, z3);
	}

	/* make Z-matrix definitions for atoms attached to 'theatom' */
	for(i=0;i<theatom->nneighbors;i++) nbatom[i] = theatom->nbatom[i];

	/* this takes cares of improper torsions */
	if (*nzmat >= 3) zref = *nzmat+1;
	else zref = 0;
	nzmat_old = *nzmat;

	for(i=0;i<theatom->nneighbors;i++) {
		atom = nbatom[i];
		if (atom->flags & A_MARKED) continue;

		/* this takes of improper torsions */
		if (*nzmat > nzmat_old && zref > 0) z2 = zref;

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "Original Z-matrix: N = %d Z0=%d, Z1=%d Z2=%d\n", *nzmat+1, z0, z1, z2);
}

		z1_new = z1;
		z2_new = z2;
		if (*nzmat >= 5 && pick_zmat_connection(*zmatlist, atom, *nzmat+1, z0, &z1_new, &z2_new)) {
			z1 = z1_new;
			z2 = z2_new;
			if (zref > 0) zref = *nzmat+1;

if (getmoldebuglevel() > 0) {
	fprintf(stderr, " Re-Picked Z-matrix: N=%d Z0=%d, Z1=%d Z2=%d\n", *nzmat+1, z0, z1, z2);
}

		}

		if (!(zmat = NewZmatFromAtom(atom, z0, z1, z2))) return 0;

		zmat->flags = (unsigned long)atom;	/* store atom pointer temporarily */
		EnterZmat(zmat, zmatlist);
		atom->flags |= A_MARKED;
		*nzmat = *nzmat+1;

		if (*nzmat == 2) {
			zmat->zdef[0] = 1;
			zmat->zdef[1] = 0;
			zmat->zdef[2] = 0;
		} else if (*nzmat == 3) {
			zmat->zdef[1] = 3 - zmat->zdef[0];
			zmat->zdef[2] = 0;
		} else if (*nzmat > 3) {
			if (zmat->zdef[1] <= 0) {
				zmat->zdef[1] = 2;
				zmat->zdef[2] = 3;
			} else if (zmat->zdef[2] <= 0) {
				zmat->zdef[2] = 3;
			}
		}

/*
fprintf(stderr, "Zmat: %s%d, %d-%d-%d-%d\n",
GetElemSymbol((int)atom->an), atom->serno,
*nzmat, zmat->zdef[0], zmat->zdef[1], zmat->zdef[2]);
*/
		if (*nzmat <= 2) continue;

		/* test linearity */
		zmat1 = SearchZmat(*zmatlist, zmat->zdef[0]);
		zmat2 = SearchZmat(*zmatlist, zmat->zdef[1]);
		/*
		if (IsLinear(zmat->cart, zmat1->cart, zmat2->cart, 5.0))
		*/
		if (IsAtomLinear((AtomPtr)zmat->flags, (AtomPtr)zmat1->flags, (AtomPtr)zmat2->flags, 5.0)) {
/*
fprintf(stderr, "THREE ATOMS LIE IN A STRAIGHT LINE.\n");
*/
			if (i != theatom->nneighbors-1) {	/* undo */
/*
fprintf(stderr, "   SWAPPING %d and %d\n", i,i+1);
*/
				nbatom[i] = nbatom[i+1];
				nbatom[i+1] = atom;
				DeleteZmat(zmat, zmatlist);
				atom->flags &= ~A_MARKED;
				*nzmat = *nzmat-1;
				zmat = NULL;
				i--;
				continue;
			}

			/* check the original neighbor atom list to see if zmat1 has neighbors
			 * other than zmat and zmat2
			 */
			a2 = (AtomPtr)zmat2->flags;
			for(j=0;j<theatom->nneighbors;j++) {
				a1 = theatom->nbatom[j];
				if ((a1->flags & A_MARKED) && a1 != atom && a1 != a2) {
					zmat->zdef[2] = zmat->zdef[1];
					zmat->zdef[1] = search_zmat_flags(*zmatlist, a1, NULL);
/*
fprintf(stderr, "   zdef[2] = %d zdef[1] = %d\n",
zmat->zdef[2], zmat->zdef[1]);
*/
					break;
				}
			}

			if (j == theatom->nneighbors) {	/* linearity */
/*
fprintf(stderr, "LINEARITY ERROR:\n");
*/
				return setmolerror("autozmat: THREE ATOMS LIE IN A STRAIGHT LINE.");
			}
		}
	}

	theatom->flags |= A_PICKED;
	for(i=0;i<theatom->nneighbors;i++) {
		atom = theatom->nbatom[i];
		if (atom->flags & A_PICKED) continue;
		j = search_zmat_flags(*zmatlist, atom, &zmat);
		if (!autozmat(clist, atom, zmat, j, nzmat, zmatlist)) return 0;
	}

	return 1;
}

/***********************************************************
 **********  Public Functions  *****************************
 ***********************************************************/

ZmatPtr	NewZmatFromAtom (AtomPtr atom, int z1, int z2, int z3)
{
	ZmatPtr	zmat;
	char	*p, label[Z_LABEL_LEN+1];

	if (!(zmat = NewZmat())) return NULL;
	zmat->flags = (unsigned long)atom;
	zmat->zdef[0] = z1;
	zmat->zdef[1] = z2;
	zmat->zdef[2] = z3;
	zmat->cart[0] = atom->x;
	zmat->cart[1] = atom->y;
	zmat->cart[2] = atom->z;
	zmat->an = atom->an;
	zmat->atom = atom->serno;
	zmat->itype = atom->type;

	zmat->charge = atom->charge;
	zmat->sigma = atom->sigma;
	zmat->epsilon = atom->epsilon;

	if ((p = GetAtomName(atom->refno))) {
		strncpy(label, p+strspn(p," "), Z_LABEL_LEN);
		label[Z_LABEL_LEN] = '\0';
		if (label[0] && strcmp(label, GetElemSymbol(atom->an)) != 0) {
			strcpy(zmat->label, label);
		}
	}

	zmat->resseq = atom->residue->seqno;
	if ((p = GetResName(atom->residue->refno)) && strcmp(p, "UNK") != 0) {
		strncpy(zmat->resname, p, Z_RES_LEN);
		zmat->resname[Z_RES_LEN] = '\0';
	} else strcpy(zmat->resname, "UNK");

	return zmat;
}

/*
 * AssignZmatProp --- assign Zmat fields based on atomic number (an), 
 * internal coords (zdef), and Cartesian coords (cart).
 *
 *  Input: an, zdef, cart
 *  Output: label, atomname, flags (Z_VAR_LENGTH, Z_VAR_ANGLE, Z_VAR_DIHED),
	zval, varname,
*/

int	AssignZmatProp (ZmatPtr zmatlist)
{
	int	nzmat=0;
	ZmatPtr	zmat, zmat1, zmat2, zmat3;

	ForEachZmat(zmatlist,zmat) {
		nzmat++;
		if (!zmat->label[0]) strcpy(zmat->label, GetElemSymbol(zmat->an));
		strcpy(zmat->atomname[0], GetElemSymbol(zmat->an));
		zmat->flags &= ~(Z_VAR_LENGTH | Z_VAR_ANGLE | Z_VAR_DIHED);
		zmat->zval[0] = zmat->zval[1] = zmat->zval[2] = 0.0;

		if (nzmat <= 1) continue;
		if (zmat->an >= 0) zmat->flags |= Z_VAR_LENGTH;
		if (!(zmat1 = SearchZmat(zmatlist, zmat->zdef[0]))) return 0;
		zmat->zval[0] = DIST(zmat->cart, zmat1->cart);
		sprintf(zmat->varname[0], "R%0d", nzmat);

		if (nzmat <= 2) continue;
		if (zmat->an >= 0) zmat->flags |= Z_VAR_ANGLE;
		if (!(zmat2 = SearchZmat(zmatlist, zmat->zdef[1]))) return 0;
		zmat->zval[1] = GetAngle(zmat->cart, zmat1->cart, zmat2->cart);
		sprintf(zmat->varname[1], "A%0d", nzmat);

		if (nzmat <= 3) continue;
		if (zmat->an >= 0) zmat->flags |= Z_VAR_DIHED;
		if (!(zmat3 = SearchZmat(zmatlist, zmat->zdef[2]))) return 0;
		zmat->zval[2] = GetDihedral(zmat->cart, zmat1->cart, zmat2->cart, zmat3->cart);
		sprintf(zmat->varname[2], "D%0d", nzmat);
	}

	return 1;
}

static void	add_atoms_from_list (MolPtr mol, ListPtr atomlist)
{
	AtomPtr	a1, a2;
	BondPtr	b;
	ListPtr	list;

	ForEachList(atomlist,list) {
		/* a1 = an atom where a dummy has been added.
		   a2 = dummy attached to a1
		*/
		if (!(a1 = (AtomPtr)list->L_ATOM1) || !(a2 = (AtomPtr)list->L_ATOM2)) continue;
		EnterAtom(a2, &a1->residue->atom);
		a2->residue = a1->residue;
		b = EnterNewBond(&mol->bond);
		b->atom1 = a2;
		b->atom2 = a1;
		b->type = B_SINGLE;

		a2->flags |= A_NEW;
	}
	GetNeighbor(mol);
}

/* CreateAutoZmatFromMol (MolPtr mol, int flags)
   --- automatically generates a Z-matrix based only on
       the connectivity information.

   --- needs mol->chain and mol->bond

   Once this function is returned, one should update
   natoms, nbonds, segments, etc of
   the molecule.

 Masks for flags
 Z_ADD_BGN_DUMMIES --- add beginning two dummies at the center atom
 */

ZmatPtr	CreateAutoZmatFromMol (MolPtr mol, int flags)
{
	ChainPtr	c, clist;
	ResiduePtr	r;
	AtomPtr	a, a1, a2, first_atom=NULL;
	BondPtr	bond;
	ZmatPtr	zmat, zmatlist=NULL;
	int	i, nzmat;
	ListPtr	list, dummylist=NULL;

	if (!mol || !(clist = mol->chain)) return NULL;

	/* Look for isolated atoms. If found, link them with
	 * the nearest atom.
	 */
	ForEachChainResAtom(clist,c,r,a1) {
		a1->flags &= ~(A_MARKED|A_PICKED);
		if (a1->nneighbors == 0) {
			if (!(a2 = find_closest_atom(clist,a1,0))) {
				setmolerror("Can't find the closest atom from %s%d.",
					GetElemSymbol((int)a1->an), a1->serno);
					return 0;
			}
			/* make a bond between a1 and a2 */
			bond = EnterNewBond(&mol->bond);
			bond->atom1 = a1;
			bond->atom2 = a2;
			bond->type = B_PARTIAL;
			a1->nbatom[a1->nneighbors] = a2;
			a2->nbatom[a2->nneighbors] = a1;
			a1->nbbond[a1->nneighbors] = bond;
			a2->nbbond[a2->nneighbors] = bond;
			a1->nneighbors++;
			a2->nneighbors++;
		}
	}

	/* add dummy atoms if necessary */
	dummylist = MakeAllDummyAtoms(clist);
	add_atoms_from_list(mol, dummylist);
	if (dummylist) FreeList(&dummylist);

	first_atom = NULL;
	if (flags & Z_ADD_BGN_DUMMIES) {
		/* add two dummies at a central atom */
		first_atom = MakeDummyAtomsAtCenter(clist, 2, &list);
		if (list) {
			add_atoms_from_list(mol, list);
			FreeList(&list);
		}
	}

	/* sort neighbor atoms by atomic number */
	ForEachChainResAtom(clist,c,r,a) {
		sort_neighbor(a->nbatom, a->nbbond, a->nneighbors);
	}

	/* Find out first heavy or dummy atom */
	if (!first_atom) {
		for(a=FirstAtom(clist);a;a=NextAtom(a)) {
			if (a->an == -1 || a->nneighbors > 1) break;
		}
		first_atom = a ? a : FirstAtom(clist);
	}

	/* start with first_atom */
	nzmat = 0;
	zmatlist = NULL;
	ForEachChainResAtom(clist,c,r,a) a->flags &= ~(A_MARKED|A_PICKED);
	if (!autozmat(clist, first_atom, NULL, 1, &nzmat, &zmatlist)) return NULL;

	/* make sure all the atoms have been processed */
	ForEachChainResAtom(clist,c,r,a1) {
		if (a1->flags & A_MARKED) continue;
		if (!(a2 = find_closest_atom(clist,a1,A_MARKED))) {
			setmolerror("Can't find the closest atom from %s%d.",
			GetElemSymbol((int)a1->an), a1->serno);
			FreeZmat(&zmatlist);
			return 0;
		}

		i = 0;
		ForEachZmat(zmatlist,zmat) {
			i++;
			if (zmat->flags == (unsigned long)a2) break;
		}

		if (!zmat || !autozmat(clist, a1, zmat, i, &nzmat, &zmatlist)) {
			setmolerror("Error in autozmat.");
			FreeZmat(&zmatlist);
			return NULL;
		}
	}

/* print out serial number info if requested */
if (getmoldebuglevel() > 0 || (flags & Z_PUT_SERNO_INFO)) {
	char	label[8];
	int	old_serno, new_serno;

	fprintf(stderr, "Changes in serial numbers during the file format conversion:\n");
	fprintf(stderr, "AtomName OldSerNo NewSerNo      Comment\n");
	fprintf(stderr, "-------------------------------------------------\n");

	old_serno = 0;
	ForEachChainResAtom(clist,c,r,a) {
		if (!(a->flags & A_NEW)) old_serno++;

		new_serno = 0;
		ForEachZmat(zmatlist,zmat) {
			new_serno++;
			if (zmat->flags == (unsigned long)a) break;
		}

		strncpy(label, GetAtomName(a->refno), 7);
		label[7] = '\0';
		TrimStr(label, " ");

		fprintf(stderr, "%7s  %6d   %6d   %s\n",
			label,
			(a->flags & A_NEW) ? 0 : old_serno,
			new_serno,
			(a->flags & A_NEW) ? "(Newly created dummy)" : "");

		a->flags &= ~A_NEW;
	}
	fprintf(stderr, "-------------------------------------------------\n");
}


	/* we have used zmat->flags for storing atom pointer.
	 * So, clear it here. This must be done before calling AssignZmatProp().
	 */
	ForEachZmat(zmatlist,zmat) zmat->flags = 0;

	AssignZmatProp(zmatlist);

	/* if atom->type is not PARAM_OPLS type, zero it */
	if (mol->param_type != PARAM_OPLS) {
		ForEachZmat(zmatlist,zmat) zmat->itype = zmat->ftype = 0;
	}
	
	/* if zmat->itype = 0, assign atomic numbers */
	ForEachZmat(zmatlist,zmat) {
		if (zmat->itype == 0) {
			if (zmat->an == ELEM_LP) zmat->itype = 99;
			else if (zmat->an == ELEM_GH) zmat->itype = 100;
			else zmat->itype = zmat->an;
		}
	}

	/* Remember to update natoms, nbonds, segments, bbox, etc. */

	return zmatlist;
}

