#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "macros.h"
#include "molecule.h"

static int	current_version=0;

int	TestBonded (AtomPtr atom1, AtomPtr atom2);

/* chain */
void	EnterChain (ChainPtr New, ChainPtr *top)
{
	ChainPtr	tmp;

	if (!New || !top) return;
	if (!*top) {
		New->prev = NULL;
		*top = New;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;
		/* append New item */
		New->prev = tmp;
		tmp->next = New;
	}
}

ChainPtr	NewChain (void)
{
	ChainPtr	New;

	if (!(New = (Chain *) malloc(sizeof(Chain)))) return NULL;
	BZERO(New, sizeof(Chain));
	return New;
}

ChainPtr	EnterNewChain (ChainPtr *top)
{
	ChainPtr	New;
	if (!(New = NewChain())) return NULL;
	EnterChain(New, top);
	return New;
}

ChainPtr	DeleteChain (ChainPtr del, ChainPtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;
	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

void	FreeChainEntry (ChainPtr *entry)
{
	ChainPtr	p;
	if (!entry || !*entry) return;
	p = *entry;

	if (p->residue) FreeResidue(&p->residue);
	if (p->bond) FreeBond(&p->bond);

	free((char *)p);
	*entry = NULL;
}

void	FreeChain (ChainPtr *top)
{
	ChainPtr	p, next;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		FreeChainEntry(&p);
		p = next;
	}
	*top = NULL;
}

ChainPtr	CopyChainEntry (ChainPtr oldc)
{
	ChainPtr	c;
	ResiduePtr	r;

	if (!oldc) return NULL;
	if (!(c = NewChain())) return NULL;
	BCOPY(oldc, c, sizeof(Chain));
	c->residue = CopyResidue(oldc->residue);
	ForEachRes(c->residue,r) r->chain = c;
	c->bond = NULL;	/* ??? */
	c->prev = c->next = NULL;

	return c;
}

ChainPtr	CopyChain (ChainPtr oldclist)
{
	ChainPtr	clist=NULL,c,oldc;

	if (!oldclist) return NULL;
	ForEachChain(oldclist,oldc) {
		if (!(c = CopyChainEntry(oldc))) continue;
		EnterChain(c, &clist);
	}

	return clist;
}

/* residue */

void	EnterResidue (ResiduePtr New, ResiduePtr *top)
{
	ResiduePtr	tmp;

	if (!New || !top) return;
	if (!*top) {
		New->prev = NULL;
		*top = New;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;
		/* append New item */
		New->prev = tmp;
		tmp->next = New;
	}
}

ResiduePtr	NewResidue (void)
{
	ResiduePtr	New;

	if (!(New = (Residue *) malloc(sizeof(Residue)))) return NULL;
	BZERO(New, sizeof(Residue));
	return New;
}

ResiduePtr	EnterNewResidue (ResiduePtr *top)
{
	ResiduePtr	New;
	if (!(New = NewResidue())) return NULL;
	EnterResidue(New, top);
	return New;
}

ResiduePtr	DeleteResidue (ResiduePtr del, ResiduePtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;
	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

void	FreeResidueEntry (ResiduePtr *entry)
{
	ResiduePtr	p;
	if (!entry || !*entry) return;
	p = *entry;

	if (p->atom) FreeAtom(&p->atom);

	free((char *)p);
	*entry = NULL;
}

void	FreeResidue (ResiduePtr *top)
{
	ResiduePtr	p, next;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		FreeResidueEntry(&p);
		p = next;
	}
	*top = NULL;
}

ResiduePtr	CopyResidueEntry (ResiduePtr oldr)
{
	ResiduePtr	r;
	AtomPtr	a, olda;

	if (!oldr) return NULL;
	if (!(r = NewResidue())) return NULL;
	BCOPY(oldr, r, sizeof(Residue));
	r->atom = CopyAtom(oldr->atom);
	ForEachAtom(r->atom,a) a->residue = r;

	r->backbone = NULL;
	for(olda=oldr->atom,a=r->atom;olda && a; olda=olda->next,a=a->next) {
		if (olda == oldr->backbone) {
			r->backbone = a;
			break;
		}
	}

	r->chain = NULL;
	r->prev = r->next = NULL;

	return r;
}

ResiduePtr	CopyResidue (ResiduePtr oldrlist)
{
	ResiduePtr	rlist=NULL, r, oldr;

	if (!oldrlist) return NULL;
	ForEachRes(oldrlist,oldr) {
		if (!(r = CopyResidueEntry(oldr))) continue;
		EnterResidue(r, &rlist);
	}

	return rlist;
}

/* molecule */
void	EnterMol (MolPtr New, MolPtr *top)
{
	MolPtr	tmp;

	if (!New || !top) return;
	if (!*top) {
		New->prev = NULL;
		*top = New;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;
		/* append New item */
		New->prev = tmp;
		tmp->next = New;
	}
}

MolPtr	NewMol (void)
{
	MolPtr	New;

	if (!(New = (Mol *) malloc(sizeof(Mol)))) return NULL;
	BZERO(New, sizeof(Mol));
	New->header = -1;
	return New;
}

MolPtr	EnterNewMol (MolPtr *top)
{
	MolPtr	New;
	if (!(New = NewMol())) return NULL;
	EnterMol(New, top);
	return New;
}

MolPtr	DeleteMol (MolPtr del, MolPtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;
	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

void	FreeMolEntry (MolPtr *entry)
{
	MolPtr	p;
	if (!entry || !*entry) return;
	p = *entry;

	if (p->ssbond) FreeBond(&p->ssbond);
	if (p->hbond)  FreeBond(&p->hbond);
	if (p->chain)  FreeChain(&p->chain);
	if (p->bond)   FreeBond(&p->bond);
	if (p->ss)     FreeSS(&p->ss);

	free((char *)p);
	*entry = NULL;
}

void	FreeMol (MolPtr *top)
{
	MolPtr	p, next;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		FreeMolEntry(&p);
		p = next;
	}
	*top = NULL;
}

MolPtr	CopyMol (MolPtr oldm)
{
	MolPtr	newm=NULL;
	ChainPtr	oldc, c;
	ResiduePtr	oldr, r;
	AtomPtr		olda, a;
	BondPtr		oldb, b;

	if (!oldm || !(newm = NewMol())) return NULL;
	BCOPY(oldm, newm, sizeof(Mol));
	newm->next = newm->prev = NULL;

	newm->ssbond = CopyBond(oldm->ssbond);
	newm->hbond = CopyBond(oldm->hbond);
	newm->chain = CopyChain(oldm->chain);
	newm->ss = CopySS(oldm->ss, oldm->chain, newm->chain);

	/* copy bonds */
	newm->bond = CopyBond(oldm->bond);
	for(oldc=oldm->chain,c=newm->chain;oldc && c;oldc=oldc->next,c=c->next) {
		for(oldr=oldc->residue,r=c->residue;oldr && r;oldr=oldr->next,r=r->next) {
			for(olda=oldr->atom,a=r->atom;olda && a;olda=olda->next,a=a->next) {
				if (olda == oldr->backbone) r->backbone = a;
				for(oldb=oldm->bond,b=newm->bond;oldb && b;oldb=oldb->next,b=b->next) {
					if (olda == oldb->atom1) b->atom1 = a;
					if (olda == oldb->atom2) b->atom2 = a;
				}
			}
		}
	}
	GetNeighbor(newm);	/* update atom->nneighbors, nbatom[], nbbond[] */

	return newm;
}

/* atom */
void	EnterAtom (AtomPtr New, AtomPtr *top)
{
	AtomPtr	tmp;

	if (!New || !top) return;
	if (!*top) {
		New->prev = NULL;
		*top = New;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;
		/* append New item */
		New->prev = tmp;
		tmp->next = New;
	}
}

AtomPtr	NewAtom (void)
{
	AtomPtr	New;

	if (!(New = (AtomPtr) malloc(sizeof(ATOM)))) return NULL;
	BZERO(New, sizeof(ATOM));
	return New;
}

AtomPtr	EnterNewAtom (AtomPtr *top)
{
	AtomPtr	New;
	if (!(New = NewAtom())) return NULL;
	EnterAtom(New, top);
	return New;
}

AtomPtr	DeleteAtom (AtomPtr del, AtomPtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;
	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

void	FreeAtomEntry (AtomPtr *entry)
{
	AtomPtr	p;
	if (!entry || !*entry) return;
	p = *entry;
	free((char *)p);
	*entry = NULL;
}

void	FreeAtom (AtomPtr *top)
{
	AtomPtr	p, next;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		FreeAtomEntry(&p);
		p = next;
	}
	*top = NULL;
}

int	SearchAtomInChainByPointer (ChainPtr clist, AtomPtr atom)
{
	int	n=0;
	AtomPtr	a;
	Residue	*r;
	Chain	*c;

	n = 0;
	ForEachChainResAtom(clist,c,r,a) {
		n++;
		if (a == atom) return n;
	}
	return 0;
}

AtomPtr	SearchAtomInChainByNumber (ChainPtr clist, int num)
{
	int	n=0;
	AtomPtr	a;
	Residue	*r;
	Chain	*c;

	n = 0;
	ForEachChainResAtom(clist,c,r,a) {
		n++;
		if (n == num) return a;
	}
	return NULL;
}

AtomPtr	SearchAtomInChainBySerno (ChainPtr clist, int serno)
{
	AtomPtr	a;
	Residue	*r;
	Chain	*c;

	ForEachChainResAtom(clist,c,r,a) {
		if (a->serno == serno) return a;
	}
	return NULL;
}

AtomPtr	NextAtom (AtomPtr atom)
{
	Chain	*c;
	Residue	*r;

	if (!atom) return NULL;
	if (atom->next) return atom->next;
	if (!(r = atom->residue)) return NULL;
	if (r->next) return r->next->atom;
	if (!(c = r->chain)) return NULL;
	if (c->next && c->next->residue) return c->next->residue->atom;
	return NULL;
}

AtomPtr	FirstAtom (ChainPtr chain)
{
	Residue	*r;
	AtomPtr	a;

	if (!chain || !(r = chain->residue) || !(a = r->atom)) return NULL;
	return a;
}

void	InitAtomData (AtomPtr atom, ResiduePtr res, int an, double x, double y, double z)
{
	BZERO(atom, sizeof(ATOM));
	atom->an = (char)an;
	atom->x = x;
	atom->y = y;
	atom->z = z;
	atom->residue = res;
	atom->refno = AN2ATOMREFNO(an);
	atom->col = GetElemCPKColor(atom->an);
}

/* note that nbatom and nbbond fields are not filled here. */
AtomPtr	CopyAtomEntry (AtomPtr olda)
{
	AtomPtr	a;

	if (!olda) return NULL;
	if (!(a = NewAtom())) return NULL;
	BCOPY(olda, a, sizeof(ATOM));
	a->prev = a->next = NULL;
	a->residue = NULL;
	a->nneighbors = 0;
	BZERO(a->nbatom, sizeof(a->nbatom));
	BZERO(a->nbbond, sizeof(a->nbbond));

	return a;
}

AtomPtr	CopyAtom (AtomPtr oldatomlist)
{
	AtomPtr	atomlist=NULL, oldatom, atom;

	if (!oldatomlist) return NULL;
	ForEachAtom(oldatomlist, oldatom) {
		if (!(atom = CopyAtomEntry(oldatom))) continue;
		EnterAtom(atom, &atomlist);
	}

	return atomlist;
}

int	CountResidueInChains (ChainPtr clist)
{
	register int	nres=0;
	register ChainPtr	c;
	register ResiduePtr	r;

	ForEachChainRes(clist,c,r) nres++;
	return nres;
}

int	CountAtomInChains (ChainPtr clist)
{
	register int	natoms=0;
	register ChainPtr	c;
	register ResiduePtr	r;
	register AtomPtr	a;

	ForEachChainResAtom(clist,c,r,a) natoms++;
	return natoms;
}

int	CountAtomInResidues (ResiduePtr rlist)
{
	register int	natoms=0;
	register ResiduePtr	r;
	register AtomPtr	a;

	ForEachResAtom(rlist,r,a) natoms++;
	return natoms;
}

int	CountAtom (AtomPtr atomlist)
{
	register int	natoms=0;
	register AtomPtr	a;

	ForEachAtom(atomlist,a) natoms++;
	return natoms;
}

/* bond */
void	EnterBond (BondPtr New, BondPtr *top)
{
	BondPtr	tmp;

	if (!New || !top) return;
	if (!*top) {
		New->prev = NULL;
		*top = New;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;
		/* append New item */
		New->prev = tmp;
		tmp->next = New;
	}
}

BondPtr	NewBond (void)
{
	BondPtr	New;

	if (!(New = (Bond *) malloc(sizeof(Bond)))) return NULL;
	BZERO(New, sizeof(Bond));
	return New;
}

BondPtr	EnterNewBond (BondPtr *top)
{
	BondPtr	New;
	if (!(New = NewBond())) return NULL;
	EnterBond(New, top);
	return New;
}

BondPtr	DeleteBond (BondPtr del, BondPtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;
	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

void	FreeBondEntry (BondPtr *entry)
{
	BondPtr	p;
	if (!entry || !*entry) return;
	p = *entry;
	free((char *)p);
	*entry = NULL;
}

void	FreeBond (BondPtr *top)
{
	BondPtr	p, next;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		FreeBondEntry(&p);
		p = next;
	}
	*top = NULL;
}

/* note that atom1 and atom2 fields are not filled here. */
BondPtr	CopyBondEntry (BondPtr oldbond)
{
	BondPtr	bond;

	if (!oldbond) return NULL;
	if (!(bond = NewBond())) return NULL;
	BCOPY(oldbond, bond, sizeof(Bond));
	bond->atom1 = bond->atom2 = NULL;
	bond->prev = bond->next = NULL;

	return bond;
}

BondPtr	CopyBond (BondPtr oldbondlist)
{
	BondPtr	bondlist=NULL, oldbond, bond;

	if (!oldbondlist) return NULL;
	ForEachBond(oldbondlist, oldbond) {
		if (!(bond = CopyBondEntry(oldbond))) continue;
		EnterBond(bond, &bondlist);
	}

	return bondlist;
}

int	CountBond (BondPtr bondlist)
{
	register int	nbonds=0;
	register BondPtr	b;

	ForEachBond(bondlist,b) nbonds++;
	return nbonds;
}

/****************************************
*** Singly Linked list handling routines
*****************************************/

SS	*NewSS (void)
{
	SS	*New;
	if (!(New = (SS *) malloc(sizeof(SS)))) return NULL;
	BZERO(New, sizeof(SS));
	return New;
}

void	EnterSS (SS *New, SS **top)
{
	SS	*tmp;
	if (!New || !top) return;
	if (!*top) *top = New;	/* first item */
	else {	/* go to the end and append the New item */
		for(tmp=(*top);tmp->next;tmp=tmp->next) ;
		tmp->next = New;
	}
}

SS	*EnterNewSS (SS **top)
{
	SS	*New;
	if (!(New = NewSS())) return NULL;
	EnterSS(New, top);
	return New;
}

SS	*CopySS (SS *oldsslist, Chain *oldclist, Chain *clist)
{
	SSPtr	sslist=NULL, oldss, ss;
	Chain	*c, *oldc;
	Residue	*r, *oldr;

	if (!oldsslist) return NULL;
	ForEachSS(oldsslist,oldss) {
		if (!(ss = NewSS())) break;
		BCOPY(oldss, ss, sizeof(SS));
		ss->next = NULL;
		EnterSS(ss, &sslist);
	}

	for(oldc=oldclist,c=clist;oldc && c;oldc=oldc->next,c=c->next) {
		for(oldr=oldc->residue,r=c->residue;oldr && r;oldr=oldr->next,r=r->next) {
			for(oldss=oldsslist,ss=sslist;oldss && ss;oldss=oldss->next,ss=ss->next) {
				if (oldr == oldss->init_res) ss->init_res = r;
				if (oldr == oldss->term_res) ss->term_res = r;
			}
		}
	}

	return sslist;
}

SS	*DeleteSS (SS *del, SS **top)
{
	SS	*ss;

	if (!top || !*top) return NULL;
	if (!del) return *top;
	if (del == *top) *top = del->next;
	else {
		for(ss=(*top);ss->next != del;ss=ss->next);
		ss->next = del->next;
	}
	return *top;
}

void	RemoveSS (SS *del, SS **top)
{
	if (!del) return;
	DeleteSS(del, top);
	del->next = NULL;
	FreeSS(&del);
}

void	FreeSS (SS **top)
{
	SS	*next, *p;
	if (!top || !*top) return;
	p = *top;
	while (p) {
		next = p->next;
		free((char *)p);
		p = next;
	}
	*top = NULL;
}










/* check if serial numbers serno1 and serno2 are bonded.
If yes, return that bond. Otherwise return NULL.
*/
BondPtr	GetSernoBond (BondPtr bondlist, int serno1, int serno2)
{
	register BondPtr	bond;

	for(bond=bondlist;bond;bond=bond->next) {
		if ((bond->atom1->serno == serno1 && bond->atom2->serno == serno2) ||
		    (bond->atom1->serno == serno2 && bond->atom2->serno == serno1)) return bond;
	}
	return NULL;
}

BondPtr	GetAtomBonded (BondPtr bondlist, AtomPtr atom1, AtomPtr atom2)
{
	register BondPtr	bond;

	ForEachBond(bondlist, bond) {
		if ((bond->atom1 == atom1 && bond->atom2 == atom2) ||
		    (bond->atom1 == atom2 && bond->atom2 == atom1)) return bond;
	}
	return NULL;
}

BondPtr	GetBondBetweenAtoms (AtomPtr atom1, AtomPtr atom2)
{
	int	i;

	if (!atom1 || !atom2) return NULL;
	for(i=0;i<atom1->nneighbors;i++) {
		if (atom1->nbatom[i] == atom2) return atom1->nbbond[i];
	}
	return NULL;
}

BondPtr	GetNeighborBond (AtomPtr atom1, AtomPtr atom2)
{
	int	i;

	if (!atom1 || !atom2) return NULL;
	for(i=0;i<atom1->nneighbors;i++) {
		if (atom1->nbatom[i] == atom2) return atom1->nbbond[i];
	}
	return NULL;
}

AtomPtr	GetSernoAtom (ChainPtr chainlist, int serno)
{
	register ChainPtr	c;
	register ResiduePtr	r;
	register AtomPtr	a;

	ForEachChainResAtom(chainlist,c,r,a) if (a->serno == serno) return a;
	return NULL;
}

int	ProcessBond (MolPtr m, int serno1, int serno2)
/* serno1 and serno2 are bonded */
{
	BondPtr	bond;

	/* check if the bond already exists. If yes, increase the bond order */
	if ((bond = GetSernoBond(m->bond, serno1, serno2))) {
		bond->type++;
		return 1;
	}
	/* create a new bond */
	bond = EnterNewBond(&m->bond);
	bond->atom1 = GetSernoAtom(m->chain, serno1);
	bond->atom2 = GetSernoAtom(m->chain, serno2);
	bond->type = B_SINGLE;

	return 1;
}

void	CalcBondsWithinResidues (MolPtr m)
{
	register int	ntest;
	register Chain	*chain;
	register Residue	*res;
	register AtomPtr	atom1, atom2;
	register BondPtr	bond;
	BondPtr	bondlist=NULL;

	ntest = 0;
	ForEachChainRes(m->chain, chain, res) {
		for(atom1=res->atom;atom1 && atom1->next;atom1=atom1->next) {
			for(atom2=atom1->next;atom2;atom2=atom2->next) {
				if (GetBondBetweenAtoms(atom1, atom2)) continue;
				ntest++;
				if (TestBonded(atom1, atom2)) {
					bond = EnterNewBond(&bondlist);
					bond->atom1 = atom1;
					bond->atom2 = atom2;
					bond->type = B_SINGLE;
				}
			}
		}
	}
	if (bondlist) EnterBond(bondlist, &m->bond);
}

/* It is assumed that neighbor information (atom->nneighbors, atom->nbatom,
   atom->nbbond) is not available yet. Use m->bond to test if two atoms are bonded.
   Newly found bonds are appeneded to m->bond.
*/
void	CalcBondsBetweenResidues (MolPtr m)
{
	register int	ntest;
	register Chain	*chain;
	register Residue	*res1, *res2;
	register AtomPtr	atom1, atom2;
	register BondPtr	bond;
	BondPtr	bondlist=NULL;

	ntest = 0;
	ForEachChain(m->chain,chain) {
		for(res1=chain->residue;res1 && res1->next;res1=res1->next) {
			for(res2=res1->next;res2;res2=res2->next) {
				for(atom1=res1->atom;atom1;atom1=atom1->next) {
					for(atom2=res2->atom;atom2;atom2=atom2->next) {
						if (GetBondBetweenAtoms(atom1, atom2)) continue;
						ntest++;
						if (TestBonded(atom1, atom2)) {
							bond = EnterNewBond(&bondlist);
							bond->atom1 = atom1;
							bond->atom2 = atom2;
							bond->type = B_SINGLE;
						}
					}
				}
			}
		}
	}
	if (bondlist) EnterBond(bondlist, &m->bond);
}

/* It is assumed that neighbor information (atom->nneighbors, atom->nbatom,
   atom->nbbond) is not available yet. Use bondlist_ret to test if two atoms are bonded.
   Newly found bonds are appeneded to bondlist_ret.
*/
void	CalcBondsBetweenTwoResidues (Residue *res1, Residue *res2, BondPtr *bondlist_ret)
{
	register AtomPtr	atom1, atom2;
	register BondPtr	bond;
	BondPtr	bondlist=NULL;

	if (!res1 || !res2 || !bondlist_ret) return;
	for(atom1=res1->atom;atom1;atom1=atom1->next) {
		for(atom2=res2->atom;atom2;atom2=atom2->next) {
			if (*bondlist_ret && GetAtomBonded(*bondlist_ret, atom1, atom2)) continue;
			if (TestBonded(atom1, atom2)) {
				bond = EnterNewBond(&bondlist);
				bond->atom1 = atom1;
				bond->atom2 = atom2;
				bond->type = B_SINGLE;
			}
		}
	}
	if (bondlist) EnterBond(bondlist, bondlist_ret);
}

#define MinBondDist	(0.4*0.4)
static float	metal_bond    = 0.8;
static float	nonmetal_bond = 1.5;

int	TestBonded (AtomPtr atom1, AtomPtr atom2)
{
	float	r1, r2, r;
	float	dist;

	r1 = GetElemBSRadius(atom1->an);
	r2 = GetElemBSRadius(atom2->an);
	r1 *= TestMetalElem(atom1->an) ? metal_bond : nonmetal_bond;
	r2 *= TestMetalElem(atom2->an) ? metal_bond : nonmetal_bond;

	r = SQUARE(r1+r2);
	dist = SQUARE(atom1->x-atom2->x) + SQUARE(atom1->y-atom2->y) + SQUARE(atom1->z-atom2->z);

	if (IsDummy(atom1->an) || IsDummy(atom2->an)) {
		return (dist < 1.01);
	} else {
		return (dist >= MinBondDist && dist <= r);
	}
}

void	SetMolBBox3 (MolPtr m)
{
	register float	xmin, ymin, zmin, xmax, ymax, zmax;
	register ChainPtr	chain;
	register ResiduePtr	res;
	register AtomPtr	atom;

	if (!m || !m->chain) return;
	xmin = ymin = zmin =  1000000;
	xmax = ymax = zmax = -1000000;
	ForEachChainResAtom(m->chain, chain, res, atom) {
		xmin = MIN(xmin, atom->x);
		ymin = MIN(ymin, atom->y);
		zmin = MIN(zmin, atom->z);
		xmax = MAX(xmax, atom->x);
		ymax = MAX(ymax, atom->y);
		zmax = MAX(zmax, atom->z);
	}

	SetBBox3(&m->bbox3, xmin, ymin, zmin, xmax, ymax, zmax);
}

AtomPtr	FindCA (Residue *res)
{
	register AtomPtr	atom;

	if (!res) return NULL;
	ForEachAtom(res->atom, atom) if (atom->refno == ATM_CA) return atom;
	return NULL;
}

Chain	*FindChainByID (Chain *clist, char id)
{	
	Chain	*c;
	ForEachChain(clist,c) if (c->id == id) return c;
	return NULL;
}

/* given a bond list (mol->bond), neighbor info (atom->nneighbors,
 * atom->nbatom, and atom->nbbond) is rebuilt.
 */
void	GetNeighbor (MolPtr mol)
{
	MolPtr  m;
	BondPtr b;
	AtomPtr a1, a2, a;
	Residue	*r;
	Chain	*c;
	int	i;

	ForEachMol(mol, m) {
		ForEachChainResAtom(m->chain, c, r, a) {
			a->nneighbors = 0;
			for(i=0;i<MAXNEIGHBOR;i++) {
				a->nbatom[i] = NULL;
				a->nbbond[i] = NULL;
			}
		}
		ForEachBond(m->bond, b) {
			a1 = b->atom1;
			a2 = b->atom2;
			if (a1->nneighbors >= MAXNEIGHBOR || a2->nneighbors >= MAXNEIGHBOR) continue;
			a1->nbatom[a1->nneighbors] = a2;
			a2->nbatom[a2->nneighbors] = a1;
			a1->nbbond[a1->nneighbors] = b;
			a2->nbbond[a2->nneighbors] = b;
			a1->nneighbors++;
			a2->nneighbors++;
		}
	}
}

static SS	*_MakeSSList (Residue *init, Residue *term, int type, char chain_id)
{
	Residue	*r;
	SS	*ss, *sslist;

	if (!init || !term || init == term) return NULL;

	/* skip no-backbone residues forwards */
	for(;!init->backbone && init && init != term->next;init=init->next);

	if (init == term) return NULL;

	/* skip no-backbone residues backwards */
	for(;!term->backbone && term && term != init->prev;term=term->prev);

	if (init == term) return NULL;

	/* check if there is an R_SSTERM between init and term */
	sslist = NULL;
	for(r=init;r && r != term->next;r=r->next) {
		if (!r->backbone) {
/*
fprintf(stderr, "%s%d BACKBONE=0\n", GetResName(r->refno), r->seqno);
*/
			return sslist;
		}

		if (r == term || (r != init && (r->ss_flags & R_SSTERM))) {
			if ((ss = EnterNewSS(&sslist))) {
				ss->init_res = init;
				ss->term_res = r;
				ss->init = init->seqno;
				ss->term = r->seqno;
				ss->type = type;
				ss->chain_id = chain_id;
			}
			init = r;
		}
	}
	return sslist;
}

SS	*MakeSSList (Chain *chainlist)
{
	ChainPtr	chain;
	ResiduePtr	res, init, term, reslast;
	AtomPtr		atom, atom_prev;
	SSPtr	ss, sslist;
	int	type, old_type;

	if (!chainlist || !chainlist->residue) return NULL;

	/* check if the residue sequence nubmers are in increasing order */
	ForEachChainRes(chainlist,chain,res) {
		atom_prev = NULL;
		ForEachAtom(res->atom,atom) {
			if (atom_prev && atom->residue->seqno < atom_prev->residue->seqno) {
				return NULL;
			}
			atom_prev = atom;
		}
	}

	sslist = NULL;
	for(chain=chainlist;chain;chain=chain->next) {
		for(res=chain->residue;res && res->next && !(res->next->flags & R_HETATM);res=res->next);
		reslast = res;
		init = chain->residue;
		old_type = init->ss_flags & SSTYPES;
		for(res=chain->residue;res;res=res->next) {
			if (res->flags & R_HETATM) break;
			if (!res->backbone) {
				ForEachAtom(res->atom,atom) {
					if (IsAlphaCarbon(atom->refno)) {
						res->backbone = atom;
						break;
					}
				}
			}
			if ((type = (res->ss_flags & SSTYPES)) != old_type || !res->backbone
				|| (res->ss_flags & R_BREAK) || res == reslast) {
				if (init && (init != res)) {
					if (res->ss_flags & R_BREAK) term = res->prev;
					else if (res == reslast) term = res;
					else term = old_type ? (res->prev ? res->prev : res) : res;

					if ((ss = _MakeSSList(init, term, old_type, chain->id))) {
						EnterSS(ss, &sslist);
					}

					if (type || (res->ss_flags & R_BREAK)) init = res;
					else init = res->prev;
				} else {
					if (!res->backbone) init = res->next;
					else if ((!init && res->backbone) || (res->ss_flags & R_BREAK)) init = res;
				}

				old_type = (init && !res->backbone) ? (init->ss_flags & SSTYPES) : type;
			}
		}
	}

/*
FPrintSS(stderr, sslist);
*/
	return sslist;
}

SS	*SearchSSList (SS *sslist, Residue *res)
{
	Residue	*r;
	SS	*ss;

	for(ss=sslist;ss;ss=ss->next) {
		for(r=ss->init_res;r && r != ss->term_res->next;r=r->next) {
			if (r == res) return ss;
		}
	}
	return NULL;
}

void	FPrintSSEntry (FILE *fp, SS *ss)
{
	fprintf(fp, "SS type = %d init = %s %d term = %s %d\n",
		ss->type, GetResName(ss->init_res->refno), ss->init_res->seqno,
		GetResName(ss->term_res->refno), ss->term_res->seqno);
}

void	FPrintSS (FILE *fp, SS *sslist)
{
	SS	*ss;
	fprintf(fp, "List of Secondary Structures: (0=Coil, 1=Helix, 2=Sheet, 4=Turn)\n");
	for(ss=sslist;ss;ss=ss->next) FPrintSSEntry(fp, ss);
}

void	FPrintAtomEntry (FILE *fp, AtomPtr atom)
{
	int	i;

	fprintf(fp, "%s %d % 9.6f % 9.6f % 9.6f",
		GetAtomName(atom->refno), atom->serno, atom->x, atom->y, atom->z);
	for(i=0;i<atom->nneighbors;i++) {
		fprintf(fp, " % 4d", atom->nbatom[i]->serno);
	}
	fprintf(fp, "\n");
}

void	FPrintAtom (FILE *fp, AtomPtr atomlist)
{
	AtomPtr	atom;
	fprintf(fp, "List of Atoms:\n");
	ForEachAtom(atomlist, atom) FPrintAtomEntry(fp, atom);
}

/* bounding box handling routines */
void	ClearBBox (BBox *bbox)
{
	bbox->x1 = 1.0E+3;
	bbox->y1 = 1.0E+3;
	bbox->x2 = 1.0E-3;
	bbox->y2 = 1.0E-3;
}

void	ZeroBBox (BBox *bbox)
{
	bbox->x1 = 0.0;
	bbox->y1 = 0.0;
	bbox->x2 = 0.0;
	bbox->y2 = 0.0;
}

void	SetBBox (BBox *bbox, float x1, float y1, float x2, float y2)
{
	bbox->x1 = x1;
	bbox->y1 = y1;
	bbox->x2 = x2;
	bbox->y2 = y2;
}

void	UpdateBBox (BBox *bbox, float x1, float y1, float x2, float y2)
{
	bbox->x1 = MIN(bbox->x1, x1);
	bbox->y1 = MIN(bbox->y1, y1);
	bbox->x2 = MAX(bbox->x2, x2);
	bbox->y2 = MAX(bbox->y2, y2);
}

void	SubtBBox (BBox *bbox, float x, float y)
{
	bbox->x1 -= x;
	bbox->y1 -= y;
	bbox->x2 -= x;
	bbox->y2 -= y;
}

void	ClearBBox3 (BBox3 *bbox)
{
	bbox->x1 = 1.0E+6;
	bbox->y1 = 1.0E+6;
	bbox->z1 = 1.0E+6;
	bbox->x2 = 1.0E-6;
	bbox->y2 = 1.0E-6;
	bbox->z2 = 1.0E-6;
}

void	SetBBox3 (BBox3 *bbox, float x1, float y1, float z1, float x2, float y2, float z2)
{
	bbox->x1 = x1;
	bbox->y1 = y1;
	bbox->z1 = z1;
	bbox->x2 = x2;
	bbox->y2 = y2;
	bbox->z2 = z2;
}

void	UpdateBBox3 (BBox3 *bbox, float x1, float y1, float z1, float x2, float y2, float z2)
{
	bbox->x1 = MIN(bbox->x1, x1);
	bbox->y1 = MIN(bbox->y1, y1);
	bbox->z1 = MIN(bbox->z1, z1);
	bbox->x2 = MAX(bbox->x2, x2);
	bbox->y2 = MAX(bbox->y2, y2);
	bbox->z2 = MAX(bbox->z2, z2);
}

void	UpdateBBox3Point (BBox3 *bbox, float x, float y, float z)
{
	bbox->x1 = MIN(bbox->x1, x);
	bbox->y1 = MIN(bbox->y1, y);
	bbox->z1 = MIN(bbox->z1, z);
	bbox->x2 = MAX(bbox->x2, x);
	bbox->y2 = MAX(bbox->y2, y);
	bbox->z2 = MAX(bbox->z2, z);
}

#define MULTIPLEBOND(t)	(t==B_AROMA||t==B_DOUBLE||t==B_DBLAROMA||t==B_TRIPLE)

/* Estimate bondlength between atom1 and atom2 based on the
 * number of their neighbors. atom1 and atom2 do not have to be bonded.
 */
float	GuessBondLength (AtomPtr atom1, AtomPtr atom2, int bondtype)
{
	int	an1, an2, n1, n2;
	float	r1, r2;

	if (!atom1 || !atom2) return 0.0;
	if (bondtype <= 0) bondtype = B_SINGLE;
	an1 = (int)atom1->an;
	an2 = (int)atom2->an;
	n1 = atom1->nneighbors;
	n2 = atom2->nneighbors;
	if (an1 > an2) {
		ISWAP(an1, an2)
		ISWAP(n1, n2)
	}

	switch (an1) {
	case 0: return 1.0;
	case 1:
		switch (an2) {
		case 1: break;	/* H-H */
		case 3: break;	/* H-Li */
		case 6:
			switch (n2) {
			case 1: break;
			case 2:	break; /* sp  (H-CR)  */
			case 3:	break; /* sp2 (H-CR2) */
			case 4:	break; /* sp3 (H-CR3) */
			}
			break;
		case 7:
			switch (n2) {
			case 1: break;
			case 2: break;	/* sp2 (H-NR) */
			case 3: return 1.08;	/* sp3 (H-NR2) */
			case 4: break;	/* sp3 (H-NR3+) */
			}
			break;
		case 8:
			switch (n2) {
			case 1: break;
			case 2: break;	/* sp3 (H-OR) */
			}
			break;
		case 9:
		case 11:
		case 14:
		case 15:
		case 16:
		case 17:
		case 19:
		case 35:
		case 53:
			break;
		}
		break;
	case 3:
		switch (an2) {
		case 3:
		case 6:
		case 7:
		case 8:
		case 9:
			break;
		}
		break;
	case 6:
		switch (an2) {
		case 6:
			switch (bondtype) {
			case B_AROMA:	return 1.400;	/* C:C */
			case B_DOUBLE:	return 1.335;	/* C=C */
			case B_DBLAROMA:
			case B_TRIPLE:	return 1.204;	/* C%C */
			}
			return 1.540;	/* C-C */
		case 7:
			switch (bondtype) {
			case B_AROMA:
			case B_DOUBLE:	return 1.270;	/* C=N */
			case B_DBLAROMA:
			case B_TRIPLE:	return 1.158;	/* C%N */
			}
			return 1.470;	/* C-N */
		case 8:
			return MULTIPLEBOND(bondtype) ? 1.220 : 1.430;
		case 14:	/* Si */
			break;
		case 15:	/* P */
			return 1.830;
		case 16:	/* S */
			return MULTIPLEBOND(bondtype) ? 1.710 : 1.817;
		}
		break;
	case 7:
		switch (an2) {
		case 7:
			switch (bondtype) {
			case B_AROMA:	return 1.330;
			case B_DOUBLE:	return 1.346;	/* N=N */
			case B_DBLAROMA:
			case B_TRIPLE:	return 1.200;	/* N%N */
			}
			return 1.500;	/* N-N */
		case 8:
			return 1.4;
		case 15:
			break;
		case 16:
			return 1.625;
		}
		break;

	case 8:
		switch (an2) {
		case 8:
			return 1.48;
		case 15:
			return MULTIPLEBOND(bondtype) ? 1.49 : 1.6;
		case 16:
			return MULTIPLEBOND(bondtype) ? 1.45 : 1.5;
		}
		break;
	case 15:
		break;
	case 16:
		switch (an2) {
		case 16:
			return MULTIPLEBOND(bondtype) ? 2.0 : 2.03;
		}
		break;
	}

	/* default bond length (sum of Bragg-Slater radii of two atoms */
	r1 = GetElemBSRadius(an1);
	r2 = GetElemBSRadius(an2);

	return (r1+r2);
}

float	BondType2Order (int type)
{
	switch (type) {
	case B_SINGLE: return 1.0;
	case B_DOUBLE: return 2.0;
	case B_TRIPLE: return 3.0;
	case B_PARTIAL:  return 0.5;
	case B_AROMA:    return 1.5;
	case B_DBLAROMA: return 2.5;
	}
	return 0.0;
}

/* try to get bond order between atom1 and atom2, based on the number of 
substituents on atom1 and atom2. Before calling this function, neighbor
information should be set correctly. */

int	GuessBondType (AtomPtr atom1, AtomPtr atom2)
{
	int	an1, an2, n1, n2;

	if (!atom1 || !atom2) return B_SINGLE;
	an1 = (int)atom1->an;
	an2 = (int)atom2->an;
	n1 = atom1->nneighbors;
	n2 = atom2->nneighbors;

	if (an1 > an2) {
		ISWAP(an1, an2)
		ISWAP(n1, n2)
	}

	switch (an1) {
	case 0:		/* unfilled valency */
	case 1:		/* H  */
	case 3:		/* Li */
	case 9:		/* F  */
	case 11:	/* Na */
	case 17:	/* Cl */
	case 19:	/* K  */
	case 35:	/* Br */
	case 53:	/* I  */
		return B_SINGLE;
	}

	switch (an2) {
	case 0:		/* unfilled valency */
	case 1:		/* H  */
	case 3:		/* Li */
	case 9:		/* F  */
	case 11:	/* Na */
	case 17:	/* Cl */
	case 19:	/* K  */
	case 35:	/* Br */
	case 53:	/* I  */
		return B_SINGLE;
	}

	switch (an1) {
	case 6:
		switch (an2) {
		case 6:
			switch (n1) {
			case 1:
			case 2:
				switch (n2) {
				case 1:
				case 2: return B_TRIPLE;
				case 3: return B_DOUBLE;
				}
				break;
			case 3:
				switch (n2) {
				case 1:
				case 2:
				case 3: return B_DOUBLE;
				}
				break;
			}
			return B_SINGLE;
		case 7:
			switch (n1) {
			case 1:
			case 2:
				switch (n2) {
				case 1: return B_TRIPLE;
				case 2: return B_DOUBLE;
				}
				break;
			case 3:
				switch (n2) {
				case 1:
				case 2: return B_DOUBLE;
				}
				break;
			}
			return B_SINGLE;
		case 8:
			switch (n1) {
			case 1:
			case 2:
			case 3:
				if (n2 == 1) return B_DOUBLE;
				break;
			}
			return B_SINGLE;
		
		case 15:
			switch (n1) {
			case 1:
			case 2:
			case 3:
				switch (n2) {
				case 1:
				case 2:
				case 3:
				case 4: return B_DOUBLE;
				}
				break;
			}
			return B_SINGLE;
		
		case 16:
			switch (n1) {
			case 1:
			case 2:
			case 3:
				if (n2 <= 4) return B_DOUBLE;
				break;
			}
			return B_SINGLE;
		}
		break;
	case 7:
		switch (an2) {
		case 7:
			switch (n1) {
			case 1:
				switch (n2) {
				case 1: return B_TRIPLE;
				case 2: return B_DOUBLE;
				}
				break;
			case 2:
				switch (n2) {
				case 1:
				case 2: return B_DOUBLE;
				}
				break;
			}
			return B_SINGLE;
		case 8:
			switch (n1) {
			case 1:
			case 2:
				if (n2 == 1) return B_DOUBLE;
				break;
			}
			return B_SINGLE;
		
		case 15:
			switch (n1) {
			case 1:
			case 2: if (n2 < 5) return B_DOUBLE;
				break;
			}
			return B_SINGLE;
		
		case 16:
			switch (n1) {
			case 1:
			case 2: if (n2 <= 4) return B_DOUBLE;
			}
			return B_SINGLE;
		}
		break;
	case 8:
		switch (an2) {
		case 8:
			if (n1 == 1 && n2 == 1) return B_DOUBLE;
			return B_SINGLE;
		case 15:
		case 16:
			if (n1 == 1 && n2 < 5) return B_DOUBLE;
			return B_SINGLE;
		}
		break;
	case 15:
		switch (an2) {
		case 15:
		case 16:
			if (n1 < 5 && n2 < 5) return B_DOUBLE;
			return B_SINGLE;
		}
		break;
	case 16:
		switch (an2) {
		case 16:
			if (n1 < 5 && n2 < 5) return B_DOUBLE;
			return B_SINGLE;
		}
		break;
	}

	return B_SINGLE;
}

static char	chem_formula[256];
char	*GetFormula (ChainPtr clist)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr		a;
	int	an, natoms[104];
	char	sym[4];

	BZERO(natoms, sizeof(natoms));
	ForEachChainResAtom(clist,c,r,a) {
		an = (int)a->an;
		an = MAX(MIN(an, 103), 0);
		natoms[an]++;
	}
	BZERO(chem_formula, sizeof(chem_formula));
	for(an=0;an<=103;an++) {
		if (natoms[an] == 0) continue;
		strcpy(sym, GetElemSymbol(an));
		if (sym[1] == ' ') sym[1] = '\0';
		sprintf(chem_formula, "%s%s%d", chem_formula, sym, natoms[an]);
	}
	if (chem_formula[0]) return chem_formula;
	else return NULL;
}

void	FPrintMol (FILE *fp, MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr		a;
	BondPtr		b;
	SSPtr		ss;
	long	flags;
	int	backbone_serno;
	char	ambertype[3];

	if (!fp || !mol || !mol->chain) return;

	fprintf(fp, "BgnMol:\n");
	fprintf(fp, "BgnAtm:\n");
	ForEachChain(mol->chain,c) {
		flags = c->flags;
		flags &= ~(C_MARKED|C_SELECTED);
		fprintf(fp, "BgnChn: %ld %d\n", flags, (int)c->id);
		ForEachRes(c->residue,r) {
			flags = r->flags;
			flags &= ~(R_MARKED|R_SELECTED);
			backbone_serno = r->backbone ? r->backbone->serno : 0;

			fprintf(fp, "BgnRes: %ld %ld %d %d %s\n",
				flags, r->ss_flags, r->seqno, backbone_serno, GetResName(r->refno));

			ForEachAtom(r->atom,a) {
				flags = a->flags;
				flags &= ~(A_MARKED|A_SELECTED|A_PICKED|A_NEW);

				/* old version (up to 0.80) */
				/*
				fprintf(fp, "%3d % 9.6f % 9.6f % 9.6f %d [%s] "\
					"%6.2f %d %9.6f %9.6f %9.6f "\
					"%d %d %.3f %.3f %d %s\n",
					(int)a->an, a->x, a->y, a->z, a->serno, GetAtomName(a->refno),
					a->tempfac, (int)a->altloc, a->label1, a->label2, a->charge,
					flags, a->col,
					a->px, a->py, a->label2D_align, a->label2D[0]
					);
				*/

				/* version 0.82 */
				strcpy(ambertype, a->ambertype);
				if (ambertype[0] == '?') ambertype[1] = '?';
				if (!isascii(ambertype[0]) || !isascii(ambertype[1])) strcpy(ambertype, "??");
				if (!ambertype[1]) ambertype[1] = ' ';
				ambertype[2] = '\0';

				fprintf(fp,
					"%3d %.6f %.6f %.6f %d [%4s] %.2f %d "\
					"%d %.6f %.6f %.6f [%2s] "\
					"%ld %d %ld %.6f %.6f "\
					"%.3f %.3f %d %s %s "\
					"%ld %ld "\
					"%ld %d %d %d\n",

					a->an, a->x, a->y, a->z, a->serno,
						MakePdbAtomName(GetAtomName(a->refno), a->an),
						a->tempfac, a->altloc,
					a->type, a->charge, a->sigma, a->epsilon, ambertype,
					flags, a->col, a->data, a->label1, a->label2,
					a->px, a->py, a->label2D_align,
						a->label2D[0] ? a->label2D : "*",
						a->usertag > 0 ? GetStrValue(a->usertag) : "*",

						/* domain properties */
						a->domain.mflags,
						a->domain.viewtype,

						/* domain label properties */
						a->domain.label.flags,
						a->domain.label.precision,
						a->domain.label.fontsize,
						a->domain.label.fontstyle
					);
			}
			fprintf(fp, "EndRes:\n");
		}
		fprintf(fp, "EndChn:\n");
	}
	fprintf(fp, "EndAtm:\n");

	fprintf(fp, "BgnBnd:\n");
	ForEachBond(mol->bond,b) {
		flags = b->flags;
		flags &= ~B_MARKED;
		fprintf(fp, "%d %d %d %ld\n", b->atom1->serno, b->atom2->serno, b->type, flags);
	}
	fprintf(fp, "EndBnd:\n");

	fprintf(fp, "BgnSS:\n");
	ForEachSS(mol->ss,ss) {
		flags = ss->flags;
		flags &= ~(S_MARKED|S_SELECTED);
		fprintf(fp, "SSEntry: %ld %d %d %d %d\n",
			flags, ss->type, (int)ss->chain_id, ss->init_res->seqno, ss->term_res->seqno);
	}
	fprintf(fp, "EndSS:\n");

	fprintf(fp, "EndMol:\n");
}

static ChainPtr	load_atom (FILE *fp)
{
	ChainPtr	clist, c;
	ResiduePtr	rlist, r;
	AtomPtr		alist, a;
	double	x, y, z, charge, sigma, epsilon;
	float	tempfac, label1, label2, px, py;
	int	an, serno, n, type, flags, col, label2D_align, altloc, zmatnum;
	char	card[256], keywd[256];	/* card and keywd should have the same size */
	char	name[64], usertag[64], label2D[32], ambertype[3];
	int	domain_mflags, domain_viewtype, domain_label_flags,
		domain_label_precision, domain_label_fontsize, domain_label_fontstyle;

	clist = c = NULL;
	rlist = r = NULL;
	alist = a = NULL;

	while (READLINE(fp, card)) {
		if (sscanf(card, "%[^:]", keywd) != 1) continue;
		if (strcmp(keywd, "BgnChn") == 0) {
			if (c) {
				c->residue = rlist;
				rlist = NULL;
			}
			c = EnterNewChain(&clist);
			sscanf(card, "%*[^:]: %ld %d", &c->flags, &n);
			c->id = (char)n;
		} else if (strcmp(keywd, "BgnRes") == 0) {
			if (!c) {	/* this is a problem */
				c = EnterNewChain(&clist);
			}
			if (r) {
				r->atom = alist;
				alist = NULL;
			}
			r = EnterNewResidue(&rlist);
			r->chain = c;

			sscanf(card, "%*[^:]: %ld %ld %d %d %s",
				&r->flags, &r->ss_flags, &r->seqno, &serno, name);
			r->refno = GetResRefno(name);
			r->backbone = (AtomPtr)serno;
		} else if (strcmp(keywd, "EndChn") == 0) {
			if (c) {
				c->residue = rlist;
				rlist = NULL;
				c = NULL;
			}
		} else if (strcmp(keywd, "EndRes") == 0) {
			if (r) {
				r->atom = alist;
				alist = NULL;
				r = NULL;
			}
		} else if (strcmp(keywd, "EndAtm") == 0) {
			if (r) r->atom = alist;
			if (c) c->residue = rlist;
			break;
		} else {	/* read atom entry */
			if (!r) {	/* problem... */
				r = EnterNewResidue(&rlist);
			}

			label2D[0] = '\0';
			usertag[0] = '\0';
			ambertype[0] = ambertype[1] = ambertype[2] = '\0';
			zmatnum = 0;
			name[4] = '\0';

			/* domain properties */
			domain_mflags = 0;
			domain_viewtype = 0;
			domain_label_flags = 0;
			domain_label_precision = 0;
			domain_label_fontsize = 0;
			domain_label_fontstyle = 0;

			if (current_version > 0 && current_version <= 800) {	/* up to version 0.800 */
				n = sscanf(card,
					"%d %lf %lf %lf %d [%4c] %f %d "\
					"%f %f %lf "\
					"%d %d %f %f %d %s",
					&an, &x, &y, &z, &serno, name, &tempfac, &altloc,
					&label1, &label2, &charge,
					&flags, &col, &px, &py, &label2D_align, label2D);
			} else if (current_version > 0 && current_version < 900) {	/* earlier than 0.9 */
				n = sscanf(card,
					"%d %lf %lf %lf %d%*[^[][%[^]]] %f %d "\
					"%d %lf %lf %lf%*[^[][%[^]]] "\
					"%d %d %d %f %f "\
					"%f %f %d %s %s",
					&an, &x, &y, &z, &serno, name, &tempfac, &altloc,
					&type, &charge, &sigma, &epsilon, ambertype,
					&flags, &col, &zmatnum, &label1, &label2,
					&px, &py, &label2D_align, label2D, usertag
					);
			} else {
				n = sscanf(card,
					"%d %lf %lf %lf %d%*[^[][%[^]]] %f %d "\
					"%d %lf %lf %lf%*[^[][%[^]]] "\
					"%d %d %d %f %f "\
					"%f %f %d %s %s "\
					"%d %d %d %d %d %d",
					&an, &x, &y, &z, &serno, name, &tempfac, &altloc,
					&type, &charge, &sigma, &epsilon, ambertype,
					&flags, &col, &zmatnum, &label1, &label2,
					&px, &py, &label2D_align, label2D, usertag,
					
					/* domain properties */
					&domain_mflags,
					&domain_viewtype,

					/* domain label properties */
					&domain_label_flags,
					&domain_label_precision,
					&domain_label_fontsize,
					&domain_label_fontstyle
					);
			}

			if (label2D[0] == '\n' || strcmp(label2D, "*") == 0) label2D[0] = '\0';
			if (usertag[0] == '\n' || !usertag[0]) strcpy(usertag, "*");
			ambertype[2] = '\0';
			if (ambertype[0] == '\n' || strcmp(ambertype, "  ") == 0)
				ambertype[0] = ambertype[1] = '?';

			if (n < 4) break;
			if (!(a = EnterNewAtom(&alist))) break;

			a->an = an;
			a->x = x;
			a->y = y;
			a->z = z;
			a->serno = serno;
			a->refno = GetPdbAtomRefno(name, r->refno);
			a->tempfac = tempfac;
			a->altloc = (char)altloc;

			a->type = type;
			a->charge = charge;
			a->sigma = sigma;
			a->epsilon = epsilon;
			strcpy(a->ambertype, ambertype);

			a->flags = flags;
			a->col = col;
			a->data = zmatnum;
			a->label1 = label1;
			a->label2 = label2;

			a->px = px;
			a->py = py;
			a->label2D_align = label2D_align;
			strcpy(a->label2D, label2D);
			a->usertag = LookUpStrTable(usertag);

			a->residue = r;

			/* domain properties */
			domain_viewtype = MAX(1, domain_viewtype);
			a->domain.mflags   = domain_mflags;
			a->domain.viewtype = domain_viewtype;
			a->domain.label.flags     = domain_label_flags;
			a->domain.label.precision = domain_label_precision;
			a->domain.label.fontsize  = domain_label_fontsize;
			a->domain.label.fontstyle = domain_label_fontstyle;
		}
	}

	/* find backbones */
	ForEachChainRes(clist,c,r) {
		if (!r->backbone) continue;
		serno = (int)r->backbone;
		r->backbone = NULL;
		ForEachAtom(r->atom,a) {
			if (a->serno == serno) {
				r->backbone = a;
				break;
			}
		}
	}

	return clist;
}

static BondPtr	load_bond (FILE *fp, ChainPtr clist)
{
	int	type, a1, a2;
	long	flags;
	AtomPtr	atom1, atom2;
	BondPtr	b, blist;
	char	card[128], keywd[128];

	if (!clist) return NULL;
	blist = NULL;
	while (READLINE(fp, card)) {
		if (sscanf(card, "%[^:]", keywd) != 1) continue;
		if (strcmp(keywd, "EndBnd") == 0) break;
		flags = 0;
		if (sscanf(card, "%d %d %d %ld", &a1, &a2, &type, &flags) >= 3 && a1 != a2) {
			if ((atom1 = GetSernoAtom(clist, a1)) &&
			    (atom2 = GetSernoAtom(clist, a2))) {
				b = EnterNewBond(&blist);
				b->flags = flags;
				b->type = type;
				b->atom1 = atom1;
				b->atom2 = atom2;
			}
		}
	}

	return blist;
}

static SSPtr	load_ss (FILE *fp, ChainPtr clist)
{
	SSPtr	sslist=NULL, ss;
	long	flags;
	int	type, init, term, chain_id;
	char	str[256], keywd[64];
	ChainPtr	c;
	ResiduePtr	init_res, term_res;

	while (READLINE(fp, str)) {
		if (sscanf(str, "%[^:]", keywd) != 1) continue;
		if (strcmp(keywd, "EndSS") == 0) break;
		else if (strcmp(keywd, "SSEntry") != 0) continue;
		if (sscanf(str, "%*[^:]:%ld %d %d %d %d", &flags, &type, &chain_id, &init, &term) != 5) continue;

		if (!(c = FindChainByID(clist, chain_id))) continue;
		if (!(init_res = FindResBySeqno(c->residue, init))) continue;
		if (!(term_res = FindResBySeqno(c->residue, term))) continue;

		if (!(ss = EnterNewSS(&sslist))) continue;
		ss->flags = flags;
		ss->type = type;
		ss->chain_id = (char)chain_id;
		ss->init = init;
		ss->term = term;
		ss->init_res = init_res;
		ss->term_res = term_res;
	}

	return sslist;
}

MolPtr	FLoadMol (FILE *fp, int version)
{
	MolPtr	mol;
	int	ierr;
	char	str[256], keywd[64];

	if (!fp) return NULL;
	if (!(mol = NewMol())) return NULL;
	current_version = version;

	ierr = 0;
	while (READLINE(fp, str)) {
		if (sscanf(str, "%[^:]", keywd) != 1) continue;
		if (strcmp(keywd, "BgnAtm") == 0) mol->chain = load_atom(fp);
		else if (strcmp(keywd, "BgnBnd")   == 0) mol->bond = load_bond(fp, mol->chain);
		else if (strcmp(keywd, "BgnSSBnd") == 0) ;
		else if (strcmp(keywd, "BgnHBnd")  == 0) ;
		else if (strcmp(keywd, "BgnSS")    == 0) mol->ss = load_ss(fp, mol->chain);
		else if (strcmp(keywd, "Flags")    == 0) sscanf(str, "%*[^:]:%ld", &mol->flags);
		else if (strcmp(keywd, "Header")   == 0) ;
		else if (strcmp(keywd, "EndMol")   == 0) break;
	}
	GetNeighbor(mol);
	return mol;
}


/*********************************************************
 *****   Determine secondary structure of protein   ******
 *********************************************************/

/* Coupling constant for Electrostatic Energy   */
/* QConst = -332 * 0.42 * 0.2 * 1000.0 [*250.0] */
#define QConst	(-6972000.0)
#define MAX_HBOND_DIST		9.0
#define MIN_HBOND_DIST		0.5
#define MAX_HBOND_SQDIST	(9.0*9.0)
#define MIN_HBOND_SQDIST	(0.5*0.5)
#define MAX_HBOND_ENERGY	1.0E6

static BondPtr	CreateHbond (AtomPtr srcCA, AtomPtr dstCA, AtomPtr src, AtomPtr dst, float energy, int offset)
{
	BondPtr	hb;

	hb = NewBond();

	hb->type = B_HBOND;
	hb->atom1 = src;
	hb->atom2 = dst;
	hb->CA1 = srcCA;
	hb->CA2 = dstCA;
	hb->energy = energy;
	hb->offset = (offset>=-128) && (offset<127) ? (char)offset : 0;

/*
fprintf(stderr, "Creating HBOND: %d-%d\n", src->serno, dst->serno);
*/

	return hb;
}

static float	CalcHbondEnergy (ResiduePtr r, AtomPtr H, AtomPtr N)
{
	double	HO_dist, HC_dist, NC_dist, NO_dist, dist;
	AtomPtr	C, O;
	float	energy;

	if (!(C = FindAtomInRes(r, ATM_C)) || !(O = FindAtomInRes(r, ATM_O))) return 0;
	HO_dist = GetAtomSQDist(H, O);
	if (HO_dist < MIN_HBOND_SQDIST) return 0;
	HO_dist = sqrt(HO_dist);

	HC_dist = GetAtomSQDist(H, C);
	if (HC_dist < MIN_HBOND_SQDIST) return 0;
	HC_dist = sqrt(HC_dist);

	NC_dist = GetAtomSQDist(N, C);
	if (NC_dist < MIN_HBOND_SQDIST) return 0;
	NC_dist = sqrt(NC_dist);

	NO_dist = GetAtomSQDist(N, O);
	if (NO_dist < MIN_HBOND_SQDIST) return 0;
	NO_dist = sqrt(NO_dist);

	dist = (HO_dist - HC_dist + NC_dist - NO_dist);
	if (ABS(dist) < 1.0E-6) energy = MAX_HBOND_ENERGY;
	else energy = QConst/dist;

	return energy;
}

void	CalcProteinHbonds (ChainPtr c1, BondPtr *hblist)
{
	ChainPtr	c2;
	ResiduePtr	r1, r2;
	AtomPtr	CA1, CA2, PC1, PO1, N1, O;
	AtomPtr	best1CA, best2CA, best1, best2;
	BondPtr	hb;
	int	pos1, pos2, offset, off1, off2;
	double	dist, CO_dist;
	float	energy, R1, R2;

	static ATOM	H;

	BZERO(&H, sizeof(ATOM));

	pos1 = 0;
	PC1 = PO1 = NULL;
	ForEachRes(c1->residue,r1) {
		pos1++;

		if (PC1 && PO1) {
			PC1 = FindAtomInRes(r1, ATM_C);
			PO1 = FindAtomInRes(r1, ATM_O);
		} else {
			PC1 = FindAtomInRes(r1, ATM_C);
			PO1 = FindAtomInRes(r1, ATM_O);
			continue;
		}

/*
fprintf(stderr, "PC1 = %d PO1 = %d\n", PC1->serno, PO1->serno);
*/

		if (!IsAmino(r1->refno) || IsProline(r1->refno)) continue;
		if (!(CA1 = FindAtomInRes(r1, ATM_CA)) || !(N1 = FindAtomInRes(r1, ATM_N))) continue;

/*
fprintf(stderr, "CA1 = %d N1 = %d\n", CA1->serno, N1->serno);
*/

		dist = GetAtomSQDist(PC1, PO1);
		if (dist > MAX_HBOND_SQDIST || dist < MIN_HBOND_SQDIST) continue;
		CO_dist = sqrt(dist);
		H.x = N1->x + (PC1->x-PO1->x)/CO_dist;
		H.y = N1->y + (PC1->y-PO1->y)/CO_dist;
		H.z = N1->z + (PC1->z-PO1->z)/CO_dist;

		c2 = c1;
		pos2 = 0;
		R1 = R2 = 0;
		best1 = best2 = best1CA = best2CA = NULL;

		ForEachRes(c2->residue,r2) {
			pos2++;
			if (r2 == r1 || r2->next == r1) continue;
			if (!IsAmino(r2->refno) || !(CA2 = FindAtomInRes(r2,ATM_CA))) continue;

			dist = GetAtomSQDist(CA1, CA2);
			if (dist > MAX_HBOND_SQDIST || dist < MIN_HBOND_SQDIST) continue;

			if (!(O = FindAtomInRes(r2, ATM_O))) continue;
			if ((energy = CalcHbondEnergy(r2, &H, N1))) {
				if (c1 == c2) offset = pos1 - pos2;
				else offset = 0;
				if (energy < R1) {
					best2CA = best1CA;
					best1CA = CA2;
					best2 = best1;
					best1 = O;
					R2 = R1;
					R1 = energy;
					off2 = off1;
					off1 = offset;
				} else if (energy < R2) {
					best2CA = CA2;
					best2 = O;
					R2 = energy;
					off2 = offset;
				}
			}
		}

		if (best1 || best2) {
			if (best2) hb = CreateHbond(CA1, best2CA, N1, best2, R2, off2);
			else hb = CreateHbond(CA1, best1CA, N1, best1, R1, off1);
			if (hb) EnterBond(hb, hblist);
		}
	}

/*
fprintf(stderr, "HBOND LIST:\n");
for(hb=(*hblist);hb;hb=hb->next) {
	fprintf(stderr, "SRC=%d DST=%d\n", hb->atom1->serno, hb->atom2->serno);
}
*/

}

void	CalcNucleicHbonds (ChainPtr clist, ChainPtr c1, BondPtr *hblist)
{
	ChainPtr	c2;
	ResiduePtr	r1, r2, best;
	AtomPtr		CA1, CA2, N1, best1;
	BondPtr	hb;
	double	max, dist;
	int	refno;

	if (!clist || !c1) return;

	ForEachRes(c1->residue,r1) {
		if (!IsPurine(r1->refno)) continue;
		if (!(N1 = FindAtomInRes(r1, ATM_N1))) continue;	/* purine N1 */

		refno = NucleicCompl(r1->refno);
		max = 5.0;
		best = NULL;

		ForEachChain(clist,c2) {
			if (c2 == c1 || !c2->residue || !IsDNA(c2->residue->refno)) continue;
			ForEachRes(c2->residue,r2) {
				if (r2->refno != refno) continue;
				if (!(CA1 = FindAtomInRes(r2,ATM_CA))) continue;
				dist = GetAtomSQDist(CA1, N1);

				best1 = CA1;
				best = r2;
				max = dist;
			}
		}
		if (!best) continue;

		CA1 = FindAtomInRes(r1, ATM_P);
		CA2 = FindAtomInRes(best, ATM_P);

		if ((hb = CreateHbond(CA1, CA2, N1, best1, 0, 0))) EnterBond(hb, hblist);
		if (IsGuanine(r1->refno)) {	/* G-C */
			if ((CA1 = FindAtomInRes(r1, ATM_N2)) &&	/* G; N2 */
			    (CA2 = FindAtomInRes(best, ATM_O2))) {	/* C; O2 */
				if ((hb = CreateHbond(NULL, NULL, CA1, CA2, 0, 0))) EnterBond(hb, hblist);
			}

			if ((CA1 = FindAtomInRes(r1, ATM_O6)) &&	/* G; O6 */
			    (CA2 = FindAtomInRes(best, ATM_N4))) {	/* C; N4 */
				if ((hb = CreateHbond(NULL, NULL, CA1, CA2, 0, 0))) EnterBond(hb, hblist);
			}
		} else {	/* A-T */
			if ((CA1 = FindAtomInRes(r1, ATM_N6)) &&	/* A; N6 */
			    (CA2 = FindAtomInRes(best, ATM_O4))) {	/* T; O4 */
				if ((hb = CreateHbond(NULL, NULL, CA1, CA2, 0, 0))) EnterBond(hb, hblist);
			}
		}
	}
}

int	CalcHbonds (MolPtr m)
{
	ChainPtr	c;
	int	refno;

	if (m->hbond) FreeBond(&m->hbond);	/* free */
	m->hbond = NULL;

	ForEachChain(m->chain, c) {
		if (!c->residue) continue;
		refno = c->residue->refno;
		if (IsProtein(refno)) CalcProteinHbonds(c, &m->hbond);
		else if (IsDNA(refno)) CalcNucleicHbonds(m->chain, c, &m->hbond);
	}

	return 1;
}

int	FindAlphaHelix (MolPtr m, int pitch, int flag)
{
	BondPtr	hb;
	ChainPtr	c;
	ResiduePtr	r, r_first, ptr;
	AtomPtr	srcCA, dstCA;
	int	prev, dist, found, nhelix;

	hb = m->hbond;
	ForEachChain(m->chain, c) {
		if (!(r_first = c->residue)) continue;
		if (!IsProtein(r_first->refno)) {
			if (IsNucleo(r_first->refno)) {
				while (hb && !IsAminoBackbone(hb->atom1->refno)) hb= hb->next;
			}
			continue;
		}

		prev = dist = 0;
		for(r=c->residue;r && hb;r=r->next) {
			if (IsAmino(r->refno) && (srcCA = FindAtomInRes(r, ATM_CA))) {
				if (dist == pitch) {
					found = 0;
					dstCA = FindAtomInRes(r_first, ATM_CA);
					while (hb && hb->CA1 == srcCA) {
						if (hb->CA2 == dstCA) found = 1;
						hb = hb->next;
					}
					if (found) {
						if (prev) {
							if (!(r_first->ss_flags & R_HELIX)) nhelix++;
							for(ptr=r_first;ptr && ptr!=r;ptr=ptr->next)
								ptr->ss_flags |= flag;
						} else prev = 1;
					} else prev = 0;
				} else while (hb && hb->CA1 == srcCA) hb = hb->next;
			} else prev = 0;

			if (r->ss_flags & R_HELIX) {
				r_first = r;
				prev = 0;
				dist = 1;
			} else if (dist == pitch) {
				r_first = r_first->next;
			} else dist++;
		}
	}

	return 1;
}

int	DetermineStruct (MolPtr m)
{
	/* calculate H-bonds */
	CalcHbonds(m);

	/* helix */
/*
FindAlphaHelix(m, 0, 0);
*/

	return 1;
}

/* create a dummy atom at p, where p is

p3---p2
     |
     |
 p---p1

 p-p1 = 1.0A, p-p1-p2 = 90, p-p1-p2-p3 = 0
*/

AtomPtr	CreateDummy (double *p1, double *p2, double *p3)
{
	double	q[3];
	AtomPtr	a;

	InternalToCartesian(q, p1, 1.0, p2, RADIAN(90.0), p3, RADIAN(0.0));
	if (!(a = NewAtom())) return NULL;

	a->an = -1;
	a->refno = AN2ATOMREFNO(a->an);
	a->col = GetElemCPKColor((int)a->an);
	a->residue = NULL;
	a->type = -1;

	a->x = q[0];
	a->y = q[1];
	a->z = q[2];

	return a;
}

/* add a dummy atom if three atoms lie in a near straight line.
   list->L_ATOM1 = existing atom where a dummy is added.
   list->L_ATOM2 = a new dummy
*/
ListPtr	MakeAllDummyAtoms (ChainPtr clist)
{
	double	p[3];
	int	i, natoms, max_serno;
	AtomPtr	a, a1, a2, a3, a4, dummy;
	ChainPtr	c1, c2, c3;
	Residue	*r1, *r2, *r3;
	ListPtr	list, dummylist;

	/* count the number of atoms */
	natoms = 0;
	max_serno = 0;
	ForEachChainResAtom(clist,c1,r1,a1) {
		natoms++;
		max_serno = MAX(max_serno, a1->serno);
	}
	if (natoms < 3) return NULL;
	max_serno = MAX(max_serno, natoms);

	dummylist = NULL;
	ForEachChainResAtom(clist,c1,r1,a1) {
		ForEachChainResAtom(clist,c2,r2,a2) {
			if (a2==a1 || !GetBondBetweenAtoms(a1,a2)) continue;

			/* max dummy an atom can have = 1 */
			ForEachList(dummylist,list) {
				if ((AtomPtr)list->L_ATOM1 == a2) break;
			}
			if (list) continue;	/* a2 is found in dummylist */

			ForEachChainResAtom(clist,c3,r3,a3) {
				if (a3==a1 || a3==a2 || !GetBondBetweenAtoms(a2,a3) ||
					!IsAtomLinear(a1,a2,a3,5.0)) continue;

				/* check if a2 has a bond other than a2-a1 and a2-a3 */
				for(i=0;i<a2->nneighbors;i++) {
					if (a2->nbatom[i] != a1 && a2->nbatom[i] != a3) {
						break;
					}
				}
				if (i != a2->nneighbors) continue;

				/* a1-a2-a3 lie in a straight line and a2 doesn't have any other bond.
				 * check if a1 and a3 have other bonds, and make sure a4-a1-a2 and
				 * a2-a3-a4 are not colinear.
				 */
				for(i=0;i<a3->nneighbors;i++) {
					if (a3->nbatom[i] != a2) {
						a4 = a3->nbatom[i];
						if (!IsAtomLinear(a2,a3,a4,5.0)) {
							a = a3;
							break;
						}
					}
					a4 = NULL;
				}

				/* 4th atom couldn't be found in a3. check a1 */
				if (!a4) {
					for(i=0;i<a1->nneighbors;i++) {
						if (a1->nbatom[i] != a2) {
							a4 = a1->nbatom[i];
							if (!IsAtomLinear(a2,a1,a4,5.0)) {
								a = a1;
								break;
							}
						}
						a4 = NULL;
					}
				}

				if (!a4) {	/* 4 atoms are in a straight line */
					p[0] = a1->x-1.0;
					p[1] = a1->y-1.0;
					p[2] = a1->z-1.0;
					if (!(dummy = CreateDummy(&a2->x, &a1->x, p))) {
						return dummylist;
					}
				} else {
					if (!(dummy = CreateDummy(&a2->x, &a->x, &a4->x))) {
						return dummylist;
					}
				}

				max_serno++;
				dummy->serno = max_serno;
				dummy->residue = a2->residue;

				/* add a dummy to list */
				list = EnterNewList(&dummylist);
				list->L_ATOM1 = (size_t)a2;
				list->L_ATOM2 = (size_t)dummy;
			}
		}
	}

	return dummylist;
}

/* add atoms in atomlist to residues where those atoms belong.
 * mol --- molecule
 * atomlist --- list of atoms to be added
 *   L_ATOM1 = parent atom
 *   L_ATOM2 = new atom to be added to L_ATOM1
 * prepend --- 0 or 1
 *   0 = add atoms at the beginning of residue
 *   1 = add atoms at the end of residue
 */
void	AddAtomListToMol (MolPtr mol, ListPtr atomlist, int prepend)
{
	AtomPtr	a1, a2;
	BondPtr	b;
	ListPtr	list;

	ForEachList(atomlist,list) {
		/* a1 = an atom where a dummy has been added.
		   a2 = dummy attached to a1
		*/
		if (!(a1 = (AtomPtr)list->L_ATOM1) || !(a2 = (AtomPtr)list->L_ATOM2)) continue;
		if (prepend) {
			a2->prev = NULL;
			a2->next = a1->residue->atom;
			a1->residue->atom->prev = a2;
			a1->residue->atom = a2;
		} else {
			EnterAtom(a2, &a1->residue->atom);
		}
		a2->residue = a1->residue;
		b = EnterNewBond(&mol->bond);
		b->atom1 = a2;
		b->atom2 = a1;
		b->type = B_SINGLE;

		a2->flags |= A_NEW;
	}
	GetNeighbor(mol);
}

/* add a dummy at atom with a serial number of seno.
 */
ListPtr	MakeDummyAtAtom (AtomPtr atom, int serno)
{
	AtomPtr	a1, a2, a3, dummy=NULL;
	int	i, j, found;
	double	p[3];
	ListPtr	list;

	if (!atom) return NULL;
	a1 = atom;
	a2 = a3 = NULL;
	for(i=0,found=0;i<a1->nneighbors;i++) {
		a2 = a1->nbatom[i];
		for(j=0;j<a2->nneighbors;j++) {
			a3 = a2->nbatom[j];
			if (a3 == a1) continue;
			if (!IsAtomLinear(a1,a2,a3,5.0)) {
				found = 1;
				break;
			}
		}
		if (found) break;
	}

	if (!found) {
		p[0] = a1->x-1.0;
		p[1] = a1->y-1.0;
		p[2] = a1->z-1.0;
		a2 = atom->nbatom[0];
		if (!a2 || !(dummy = CreateDummy(&a1->x, &a2->x, p))) {
			return NULL;
		}
	} else {
		if (!(dummy = CreateDummy(&a1->x, &a2->x, &a3->x))) {
			return NULL;
		}
	}

	if (!dummy) return NULL;
	dummy->serno = serno;
	dummy->residue = a1->residue;

	/* add a dummy to list */
	list = NewList();
	list->L_ATOM1 = (size_t)a1;
	list->L_ATOM2 = (size_t)dummy;

	return list;
}

static AtomPtr	get_linked_dummy (AtomPtr atom)
{
	int	i;
	AtomPtr	nbatom;

	for(i=0;i<atom->nneighbors;i++) {
		nbatom = atom->nbatom[i];
		if (nbatom->an == -1 && !(nbatom->flags & A_MARKED)) {
			return nbatom;
		}
	}
	return NULL;
}

/* add a specified number of dummy atoms at a central atom.
   list->L_ATOM1 = existing atom where a dummy is added.
   list->L_ATOM2 = a new dummy
   
   returns the first atom where autozmat may start from
*/
AtomPtr	MakeDummyAtomsAtCenter (ChainPtr clist, int ndummies, ListPtr *dummylist_ret)
{
	double	cx, cy, cz, dist, min_dist, heavy_min_dist;
	int	i, natoms, max_serno;
	AtomPtr	a, dummy, dummy2, center_atom, heavy_center_atom;
	Chain	*c;
	Residue	*r;
	ListPtr	list, dummylist=NULL;
	
	*dummylist_ret = NULL;

	/* count the number of atoms */
	natoms = 0;
	max_serno = 0;
	cx = cy = cz = 0;
	ForEachChainResAtom(clist,c,r,a) {
		natoms++;
		max_serno = MAX(max_serno, a->serno);
		cx += a->x;
		cy += a->y;
		cz += a->z;
	}
	if (natoms == 0) return NULL;
	if (natoms < 2) return FirstAtom(clist);
	max_serno = MAX(max_serno, natoms);

	/* find the closest atom to the center */
	cx /= (double)natoms;
	cy /= (double)natoms;
	cz /= (double)natoms;

	center_atom = a = FirstAtom(clist);
	heavy_center_atom = NULL;

	min_dist = SQUARE(a->x-cx) + SQUARE(a->y-cy) + SQUARE(a->z-cx);
	heavy_min_dist = 1.0E+12;

	ForEachChainResAtom(clist,c,r,a) {
		dist = SQUARE(a->x-cx) + SQUARE(a->y-cy) + SQUARE(a->z-cx);
		if (dist < min_dist) {
			min_dist = dist;
			center_atom = a;
		}
		if (a->an > 1 && dist < heavy_min_dist) {
			heavy_min_dist = dist;
			heavy_center_atom = a;
		}
	}

	dummylist = NULL;
	dummy = heavy_center_atom ? heavy_center_atom : center_atom;

	/* check if the center atom already has a dummy */
	if (dummy->an == -1) {
		ndummies--;
		dummy->flags |= A_MARKED;
	}

	while (ndummies > 0) {
		if ((dummy2 = get_linked_dummy(dummy))) {
			ndummies--;
			dummy = dummy2;
			dummy->flags |= A_MARKED;
		} else break;
	}

	for(i=0;i<ndummies;i++) {
		max_serno++;
		if ((list = MakeDummyAtAtom(dummy, max_serno))) EnterList(list, &dummylist);
		else break;
		a = (AtomPtr)list->L_ATOM2;
		a->nneighbors = 1;
		a->nbatom[0] = dummy;
		dummy = a;
	}

	*dummylist_ret = dummylist;

	return dummy;	/* first atom in the Z-matrix */
}

int	PrependDummyAtomsToMol (MolPtr mol, int ndummies)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a, parent_atom=NULL;
	ListPtr	list, dummylist;
	int	i, max_serno;

	if (!mol) return 0;
	if (ndummies <= 0) return 1;	/* nothing to do */

	max_serno = 0;
	ForEachChainResAtom(mol->chain,c,r,a) {
		if (!parent_atom) parent_atom = a;
		max_serno = MAX(max_serno, a->serno);
	}

	if (!parent_atom) return 0;
	dummylist = NULL;
	for(i=0;i<ndummies;i++) {
		max_serno++;
		if ((list = MakeDummyAtAtom(parent_atom, max_serno))) EnterList(list, &dummylist);
		else {
			if (dummylist) FreeList(&dummylist);
			return setmolerror("Not enough memory.");
		}
		a = (AtomPtr)list->L_ATOM2;
		a->nneighbors = 1;
		a->nbatom[0] = parent_atom;
		parent_atom = a;
	}

	AddAtomListToMol(mol, dummylist, 1);
	return 1;
}

/* delete atoms if (atom->flags & A_MARKED) != NULL */

void	RemoveMarkedAtomsFromMol (MolPtr mol)
{
	ChainPtr	c, c1, cnext, delchainlist;
	ResiduePtr	r, r1, rnext, delreslist;
	AtomPtr	a, anext, a1, a2, delatomlist;
	BondPtr	b, bnext;
	int	first=1;

	if (!mol || !mol->chain) return;

	delatomlist = NULL;
	delreslist = NULL;
	delchainlist = NULL;

	ForEachChainRes(mol->chain,c,r) {
		for(a=r->atom,first=1;a;a=anext) {
			anext = a->next;
			if (!(a->flags & A_MARKED)) continue;
			if (a == r->backbone) {
				r->backbone = NULL;
				r->ss_flags |= R_BREAK;
			}

/*
			if (first) {
				first = 0;
				if (!(r->flags & R_HETATM)) {
					r->refno = GetResRefno("UNK");
				}
			}
*/

			DeleteAtom(a, &r->atom);
			a->prev = a->next = NULL;
			EnterAtom(a, &delatomlist);

			ForEachChainResAtom(mol->chain,c1,r1,a1) {
				if (a1->serno > a->serno) a1->serno--;
			}
		}
	}

	/* remove bonds */
	for(b=mol->bond;b;b=bnext) {
		bnext = b->next;
		a1 = b->atom1;
		a2 = b->atom2;
		if ((a1->flags & A_MARKED) || (a2->flags & A_MARKED)) {
			DeleteBond(b, &mol->bond);
			FreeBondEntry(&b);
		}
	}

	/* update neighbor list */
	GetNeighbor(mol);

	/* now delete residues and chains if necessary */
	for(c=mol->chain;c;c=cnext) {
		cnext = c->next;
		for(r=c->residue;r;r=rnext) {
			rnext = r->next;
			/* if there is no atom left, remove the residue too */
			if (!r->atom) {	/* delete the residue */
				if ((r->ss_flags & R_BREAK) && r->next) r->next->ss_flags |= R_BREAK;
				DeleteResidue(r, &c->residue);
				r->prev = r->next = NULL;
				EnterResidue(r, &delreslist);
			}
		}
		if (!c->residue) {	/* delete the chain */
			DeleteChain(c, &mol->chain);
			c->prev = c->next = NULL;
			EnterChain(c, &delchainlist);
		}
	}

	/* free atoms, residues and chains which have been deleted */
	if (delatomlist) FreeAtom(&delatomlist);
	if (delreslist) FreeResidue(&delreslist);
	if (delchainlist) FreeChain(&delchainlist);
}

void	RemoveDummiesFromMol (MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;

	if (!mol || !mol->chain) return;

	ForEachChainResAtom(mol->chain,c,r,a) {
		if (a->an < 0) a->flags |= A_MARKED;
		else a->flags &= ~A_MARKED;
	}
	RemoveMarkedAtomsFromMol(mol);
}

void	FixAtomSerno (MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	int	n;

	if (!mol || !mol->chain) return;
	n = 0;
	ForEachChainResAtom(mol->chain,c,r,a) {
		n++;
		a->serno = n;
	}
}

/* maximum number of bonds each atom can possess.
*/
static int	maxbonds[MAXELEM+1] = {0,
	1, 0, 1, 2, 3, 4, 4, 2, 	/* 1-8 */
	1, 0, 1, 2, 6, 6, 5, 6, 	/* 9-16 */
	1, 0, 1, 2, 6, 6, 6, 6, 	/* 17-24 */
	8, 6, 6, 6, 6, 6, 3, 4, 	/* 25-32 */
	3, 2, 1, 0, 1, 2, 6, 6, 	/* 33-40 */
	6, 6, 6, 6, 6, 6, 6, 6, 	/* 41-48 */
	3, 4, 3, 2, 1, 0, 1, 2, 	/* 49-56 */
	12,6, 6, 6, 6, 6, 6, 6, 	/* 57-64 */
	6, 6, 6, 6, 6, 6, 6, 6, 	/* 65-72 */
	6, 6, 6, 6, 6, 6, 6, 6, 	/* 73-80 */
	3, 4, 3, 2, 1, 0, 1, 2, 	/* 81-88 */
	6, 6, 6, 6, 6, 6, 0, 0, 	/* 89-96 */
	0, 0, 0, 0, 0, 0, 0, 0}; 	/* 97-104 */

static float	get_dist_delta (AtomPtr atom1, AtomPtr atom2)
{
	float	rad1, rad2, dist, ref_dist, del;

	rad1 = GetElemBSRadius(atom1->an);
	rad1 *= TestMetalElem(atom1->an) ? metal_bond : nonmetal_bond;
	rad2 = GetElemBSRadius(atom2->an);
	rad2 *= TestMetalElem(atom2->an) ? metal_bond : nonmetal_bond;
	ref_dist = rad1+rad2;

	dist = ATOMDIST(atom1, atom2);
	del = ABS(dist-ref_dist);

	return del;
}

void	FixConnectionTable (MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a, nbatom, nbatom1, nbatom2;
	BondPtr	nbbond, delbond1, delbond2, delbond3, dist_bond[MAXNEIGHBOR];
	double	rad1, rad2, dist, bond_dist, dist_delta[MAXNEIGHBOR];
	float	minangle, angle;
	int	i, j, n;

if (getmoldebuglevel() > 0) {
fprintf(stderr, "Fixing connection table...\n");
}

	if (!mol) return;
	GetNeighbor(mol);
	ForEachChainResAtom(mol->chain,c,r,a) {
		if (IsDummy(a->an)) continue;
		for(i=0,n=0;i<a->nneighbors;i++) if (!IsDummy(a->nbatom[i]->an)) n++;
		if (n <= maxbonds[a->an]) continue;

		/* remove bad H-H, Cl-Cl contacts, etc */
		for(i=0,j=0;i<a->nneighbors;i++) {
			if (maxbonds[a->an] != 1) continue;
			nbatom = a->nbatom[i];
			nbbond = a->nbbond[i];
			if (IsDummy(nbatom->an) || maxbonds[nbatom->an] != 1) continue;

			/* delete the bond */
			j++;
			DeleteBond(nbbond, &mol->bond);
			FreeBondEntry(&nbbond);

if (getmoldebuglevel() > 0) {
fprintf(stderr, "Deleting a bad contact between monovalent atoms; %d-%d\n", a->serno, nbatom->serno);
}

		}
		n -= j;
		if (j > 0) GetNeighbor(mol);
	}

	/* delete too long or too short contacts */
	ForEachChainResAtom(mol->chain,c,r,a) {
		if (IsDummy(a->an)) continue;
		for(i=0,n=0;i<a->nneighbors;i++) if (!IsDummy(a->nbatom[i]->an)) n++;

		if (n <= maxbonds[a->an]) continue;

		rad1 = GetElemBSRadius(a->an);
		rad1 *= TestMetalElem(a->an) ? metal_bond : nonmetal_bond;

		for(i=0;i<a->nneighbors;i++) {
			nbatom = a->nbatom[i];
			nbbond = a->nbbond[i];

			dist_delta[i] = 0.0;
			dist_bond[i] = nbbond;

			if (IsDummy(nbatom->an)) continue;

			rad2 = GetElemBSRadius(nbatom->an);
			rad2 *= TestMetalElem(nbatom->an) ? metal_bond : nonmetal_bond;
			
			bond_dist = (rad1+rad2);
			dist = ATOMDIST(a,nbatom);

			dist_delta[i] = bond_dist - dist;
			dist_bond[i] = nbbond;
		}

		/* sort dist_delta in descending order */
		for(i=0;i<a->nneighbors-1;i++) {
			for(j=i+1;j<a->nneighbors;j++) {
				if (dist_delta[i] > dist_delta[j]) {
					dist = dist_delta[j];
					dist_delta[j] = dist_delta[i];
					dist_delta[i] = dist;
					
					nbbond = dist_bond[j];
					dist_bond[j] = dist_bond[i];
					dist_bond[i] = nbbond;
				}
			}
		}

		/* mark bonds to delete */
		for(i=0;i<a->nneighbors;i++) {
			nbbond = dist_bond[i];
			if (IsDummy(nbbond->atom1->an) || IsDummy(nbbond->atom2->an)) continue;
			n--;

if (getmoldebuglevel() > 0) {
fprintf(stderr, "Deleting too long or too short contact; %d-%d\n",
	nbbond->atom1->serno, nbbond->atom2->serno);
}
			DeleteBond(nbbond, &mol->bond);
			FreeBondEntry(&nbbond);

			if (n <= maxbonds[a->an]) break;
		}

		/* update the neighbor info */
		GetNeighbor(mol);
	}

	/* check for bad bond angles */
	ForEachChainResAtom(mol->chain,c,r,a) {
		if (IsDummy(a->an)) continue;
		for(i=0,n=0;i<a->nneighbors;i++) if (!IsDummy(a->nbatom[i]->an)) n++;

		minangle = 55.0;
		delbond1 = NULL;
		delbond2 = NULL;
		delbond3 = NULL;

		for(i=0;i<a->nneighbors;i++) {
			nbatom1 = a->nbatom[i];
			if (IsDummy(nbatom1->an)) continue;

			for(j=0;j<a->nneighbors;j++) {
				nbatom2 = a->nbatom[j];
				if (IsDummy(nbatom2->an)) continue;
				if (nbatom1 == nbatom2) continue;

				angle = GetAtomAngle(nbatom1, a, nbatom2);
				if (angle > minangle) continue;
				minangle = angle;
				delbond1 = a->nbbond[i];
				delbond2 = a->nbbond[j];
				delbond3 = GetAtomBonded(mol->bond, nbatom1, nbatom2);
			}
		}

		if (!delbond1) continue;

		dist_delta[0] = get_dist_delta(delbond1->atom1, delbond1->atom2);
		dist_bond[0] = delbond1;
		dist_delta[1] = get_dist_delta(delbond2->atom1, delbond2->atom2);
		dist_bond[1] = delbond2;
		if (delbond3) {
			dist_delta[2] = get_dist_delta(delbond3->atom1, delbond3->atom2);
			dist_bond[2] = delbond3;
		} else {
			dist_delta[2] = -1.0E6;
			dist_bond[2] = NULL;
		}

		/* sort dist_delta in descending order */
		for(i=0;i<2;i++) {
			for(j=i+1;j<3;j++) {
				if (dist_delta[i] > dist_delta[j]) {
					dist = dist_delta[j];
					dist_delta[j] = dist_delta[i];
					dist_delta[i] = dist;
					
					nbbond = dist_bond[j];
					dist_bond[j] = dist_bond[i];
					dist_bond[i] = nbbond;
				}
			}
		}
		
		/* delete the worst bond */
		if (!dist_bond[0]) continue;
if (getmoldebuglevel() > 0) {
fprintf(stderr, "Deleting a bond forming a bad angle; %d-%d\n",
	dist_bond[0]->atom1->serno, dist_bond[0]->atom2->serno);
}
		DeleteBond(dist_bond[0], &mol->bond);
		FreeBondEntry(&dist_bond[0]);

		GetNeighbor(mol);
	}
}

