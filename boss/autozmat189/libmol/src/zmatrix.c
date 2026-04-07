#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "macros.h"
#include "molecule.h"

typedef double	double3[3];

/******************************************************/
/*** Zmat doubly linked list **************************/
/******************************************************/

void	EnterZmat (ZmatPtr newptr, ZmatPtr *top)
{
	ZmatPtr	tmp;

	if (!newptr || !top) return;

	if (!*top) {
		newptr->prev = NULL;
		*top = newptr;
	} else {
		/* go to the last item */
		for(tmp = *top; tmp->next; tmp = tmp->next) ;

		/* append newptr item */
		newptr->prev = tmp;
		tmp->next = newptr;
	}
}

ZmatPtr	NewZmat (void)
{
	ZmatPtr	newptr;

	if (!(newptr = (ZmatPtr) malloc(sizeof(Zmat)))) return NULL;
	BZERO(newptr, sizeof(Zmat));
	return newptr;
}

ZmatPtr	EnterNewZmat (ZmatPtr *top)
{
	ZmatPtr	newptr;

	if (!(newptr = NewZmat())) return NULL;
	EnterZmat(newptr, top);
	return newptr;
}

ZmatPtr	CopyZmat (ZmatPtr oldlist)
{
	ZmatPtr	old, newptr, newlist;

	if (!oldlist) return NULL;

	newlist = NULL;
	for(old=oldlist;old;old=old->next) {
		if (!(newptr = NewZmat())) {
			if (newlist) FreeZmat(&newlist);
			return NULL;
		}
		BCOPY((void *)old, (void *)newptr, (int)(sizeof(Zmat)));
		newptr->next = newptr->prev = NULL;
		EnterZmat(newptr, &newlist);
	}
	return newlist;
}

/* allocate Zmat sructures in a contiguous memory space */
ZmatPtr	AllocZmat (int n)
{
	ZmatPtr	top, p, q;
	int	i;

	if (n <= 0) return NULL;
	if (!(top = (ZmatPtr) malloc(n * sizeof(Zmat)))) return NULL;
	BZERO((void *)top, sizeof(n * sizeof(Zmat)));

	p = &top[0];
	p->prev = NULL;
	if (n == 1) {
		p->next = NULL;
		return p;
	} else p->next = &top[1];

	for(i=1;i<n;i++) {
		q = p->next;
		q->prev = p;
		q->next = (i<n-1) ? &top[i+1] : NULL;
		p = q;
	}

	return &top[0];
}

/*
top;	pointer to Zmat structure which has been previously allocated by AllocZmat 
nitems;	new number of Zmat structures to be allocated 
*/
ZmatPtr	ReallocZmat (ZmatPtr top, int nitems)
{
	int	n, i;
	ZmatPtr	p, q;

	if (!top) return NULL;
	for(n=0, p = top; p; n++, p = p->next) ;

	if (n == 0 || nitems == n) return top;

	if (!(top = (ZmatPtr) realloc(top, nitems * sizeof(Zmat)))) return NULL;

	p = &top[0];
	p->prev = NULL;
	if (nitems == 1) {
		p->next = NULL;
		return p;
	} else p->next = &top[1];

	for(i=1;i<nitems;i++) {
		q = p->next;
		q->prev = p;
		q->next = (i<nitems-1) ? &top[i+1] : NULL;
		p = q;
	}

	return &top[0];
}

/* print out Zmat sructure info */
void	PrintZmat (ZmatPtr top)
{
	while (top) {
		fprintf(stderr, "%3d % 9.6f % 9.6f % 9.6f   %3d %3d %3d % 9.6f % 9.3f % 9.3f\n",
			top->an, top->cart[0], top->cart[1], top->cart[2],
			top->zdef[0], top->zdef[1], top->zdef[2],
			top->zval[0], top->zval[1], top->zval[2]);
		top = top->next;
	}
}

/* return a pointer to the i-th list item.
if there is no match, NULL is returned. */
ZmatPtr	SearchZmat (ZmatPtr top, int i)
{
	int	count = 1;

	if (i <= 0) return NULL;
	while (top) {
		if (i == count) return top;
		top = top->next;
		count++;
	}
	return NULL;
}

/* return a pointer to the zmatrix in which
atom_serno == zmat->atom.
if there is no match, NULL is returned. */
ZmatPtr	SearchZmatByAtomSerno (ZmatPtr top, int atom_serno)
{
	if (atom_serno <= 0) return NULL;
	while (top) {
		if (atom_serno == top->atom) return top;
		top = top->next;
	}
	return NULL;
}

int	SearchZmatByPointer (ZmatPtr top, ZmatPtr zmat)
{
	int	n=0;

	if (!zmat) return 0;
	while (top) {
		n++;
		if (top == zmat) return n;
		top = top->next;
	}
	return 0;
}


int	FindZmatByName (ZmatPtr zmatlist, char *name)
{
	int	count=0;
	ZmatPtr	zmat;

	ForEachZmat(zmatlist,zmat) {
		count++;
		if (strcmp(zmat->atomname[0], name) == 0) return count;
	}
	return 0;
}

/* delete a given item from a doubly linked list.
return top item in the list. */
ZmatPtr	DeleteZmat (ZmatPtr del, ZmatPtr *top)
{
	if (!top) return NULL;
	if (!del) return *top;

	if (del->prev) del->prev->next = del->next;
	else *top = del->next;
	if (del->next) del->next->prev = del->prev;
	return *top;
}

/* newlist is inserted before node */
void	InsertZmat (ZmatPtr *top, ZmatPtr node, ZmatPtr newlist)
{
	ZmatPtr	newlistend, prev;

	if (!top || !node || !newlist) return;
	/* go to the end */
	for(newlistend=newlist;newlistend && newlistend->next;newlistend=newlistend->next) ;
	for(prev=(*top);prev;prev=prev->next) {
		if (prev->next == node) break;
	}
	if (node == (*top)) {
		*top = newlist;
		(*top)->prev = NULL;
		newlistend->next = node;
		node->prev = newlistend;
	} else {
		if (!prev) return;
		prev->next = newlist;
		newlist->prev = prev;
		newlistend->next = node;
		node->prev = newlistend;
	}
}

void	FreeZmat (ZmatPtr *top)
{
	/* free the atom list starting from top */

	ZmatPtr	p, tmp;

	if (!top || !*top) return;
	p = *top;
	while (p) {
		tmp = p->next;
		free((char *)p);
		p = tmp;
	}
	*top = NULL;
}

/* free only one Zmat entry */
void	FreeZmatEntry (ZmatPtr *entry)
{
	if (!entry || !*entry) return;
	free(*entry);
	*entry = NULL;
}

/******************************************************/
/*** Z-matrix list handling routines ******************/
/******************************************************/

ListPtr	CopyZmatList (ListPtr oldzmatlist)
{
	ListPtr	zmatlist=NULL, oldlist, list;

	if (!oldzmatlist) return NULL;
	ForEachList(oldzmatlist, oldlist) {
		list = EnterNewList(&zmatlist);
		list->L_ZMAT = (size_t)CopyZmat((ZmatPtr)oldlist->L_ZMAT);
	}
	return zmatlist;
}

ZmatPtr	LinkZmatList (ListPtr zmatlist)
{
	ListPtr	l;
	ZmatPtr	zmat, newzmatlist=NULL;
	ForEachList(zmatlist,l) {
		if (!(zmat = (ZmatPtr)l->L_ZMAT)) continue;
		EnterZmat(zmat, &newzmatlist);
	}
	return newzmatlist;
}

void	BreakZmatList (ListPtr zmatlist)
{
	ListPtr	l;
	ZmatPtr	zmat;

	ForEachList(zmatlist,l) {
		if (!(zmat = (ZmatPtr)l->L_ZMAT)) continue;
		if (zmat->prev) zmat->prev->next = NULL;
		zmat->prev = NULL;
	}
}

void	FreeZmatList (ListPtr *zmatlist)
{
	ListPtr	l;
	ZmatPtr	zmat;

	ForEachList(*zmatlist, l) {
		if ((zmat = (ZmatPtr)l->L_ZMAT)) FreeZmat(&zmat);
	}
	FreeList(zmatlist);
}

void	PrintZmatList (ListPtr zmatlist)
{
	ListPtr	l;
	ZmatPtr	zmat;

	ForEachList(zmatlist, l) {
		if ((zmat = (ZmatPtr)l->L_ZMAT)) PrintZmat(zmat);
	}
}

/******************************************************/

int	CheckInternalCoord (int z0, int z1, int z2, int z3, double len, double theta, double phi)
{
	return !((z0 > 1 && (z1 >= z0)) ||
	(z0 > 2 && (z2 >= z0 || z1 == z2)) ||
	(z0 > 3 && (z3 >= z0 || z1 == z3 || z2 == z3)) ||
	(z0 > 1 && len < 1.0E-6) ||
	(z0 > 2 && (theta < 0.0 || theta > 180.0)));
}

static void	acrossb (double *a, double *b, double *c)
{
	double	s, x, y, z;
	/* c = a x b
	 * c is returned (normalized)
	 */
	x = a[1] * b[2] - b[1] * a[2];
	y = b[0] * a[2] - a[0] * b[2];
	z = a[0] * b[1] - b[0] * a[1];

	s = x*x + y*y + z*z;
	if (s < 1.0E-12) s = 1.0E-12;	/* to avoid a division by zero */

	s = 1.0/sqrt(s);
	c[0] = x*s;
	c[1] = y*s;
	c[2] = z*s;
}

void	InternalToCartesian (double *a1, double *a2, double len, double *a3, double theta, double *a4, double phi)
{
	int	i;
	double	xp[3], yp[3], zp[3], vca[3], vcb[3];
	double	st, xpd, ypd, zpd;

	/* get the coordinate of an atom, a1, where bond distance between a1 and a2
	 * is "length", bond angle formed by a1, a2, and a3 is "angle", and dihedral
	 * angle formed by a1, a2, a3, and a4 is "dihedral".
	 */

	st  = sin(theta);
	xpd =  len*st*cos(phi);
	ypd =  len*st*sin(phi);
	zpd = -len*cos(theta);

	/* determine transformation (direction cosine) matrix */
	for(i=0;i<3;i++) {
		vca[i] = a4[i] - a2[i];
		vcb[i] = a3[i] - a2[i];
	}
	acrossb(vca, vcb, yp);
	acrossb(vcb, yp, xp);
	acrossb(xp, yp, zp);

	/* transform to original system centered on a3 and recenter relative to the origin */
	a1[0] = xp[0] * xpd + yp[0] * ypd + zp[0] * zpd + a2[0];
	a1[1] = xp[1] * xpd + yp[1] * ypd + zp[1] * zpd + a2[1];
	a1[2] = xp[2] * xpd + yp[2] * ypd + zp[2] * zpd + a2[2];
}

int	ZmatToCartesian (ZmatPtr zmatlist)
{
	ZmatPtr	zmat;
	double3	*cart;
	double	len, theta, phi, y;
	int	n, i, z1, z2, z3;

	if (!zmatlist) return 0;
	for(zmat=zmatlist,n=0;zmat;zmat=zmat->next,n++) {
		for(i=0;i<3;i++) zmat->cart[i] = 0;
	}
	if (n == 1) return 1;

	/* allocate an array to hold Cartesian coords */
	if (!(cart = (double3 *)calloc(n*3, sizeof(double)))) return 0;

	zmat = zmatlist;	/* first Z-matrix */
	zmat = zmat->next;	/* 2nd */
	len = zmat->zval[0];
	zmat->cart[0] = cart[1][0] = len;
	if (n == 2) return 1;

	zmat = zmat->next;	/* 3rd */
	len = zmat->zval[0];
	theta = RADIAN(zmat->zval[1]);
	zmat->cart[1] = cart[2][1] = len * sin(theta);
	y = len * cos(theta);
	zmat->cart[0] = cart[2][0] = zmat->zdef[0] == 1 ? y : cart[1][0] - y;
	if (n == 3) return 1;

	zmat = zmat->next;	/* 4th */
	for(i=3;i<n;i++,zmat=zmat->next) {
		z1 = zmat->zdef[0]-1;
		z2 = zmat->zdef[1]-1;
		z3 = zmat->zdef[2]-1;
		len   = zmat->zval[0];
		theta = RADIAN(zmat->zval[1]);
		phi   = RADIAN(zmat->zval[2]);

		/* internal to Cartesian conversion */
		InternalToCartesian(cart[i], cart[z1], len, cart[z2], theta, cart[z3], phi);
	}

	for(zmat=zmatlist,i=0;zmat;zmat=zmat->next,i++) {
		zmat->cart[0] = cart[i][0];
		zmat->cart[1] = cart[i][1];
		zmat->cart[2] = cart[i][2];
	}
	if (cart) free((char *)cart);

	return 1;
}

int	ZmatListToCartesian (ListPtr zmatlist)
{
	int	ok;
	ZmatPtr	zmat=NULL;

	if (!(zmat = LinkZmatList(zmatlist))) return 0;
	ok = ZmatToCartesian(zmat);
	BreakZmatList(zmatlist);
	return ok;
}

/*
 * return the atomic symbol part of a GAUSSIAN or MOPAC Z-matrix center
 * specification, which should be composed of "atomic symbol + digits".
 */

static char	atom_symbol[128];

char	*GetSymbolFromZmatCenter (char *str)
{
	register short	i;
	register char	*p;
	char	*digit_index;

	atom_symbol[0] = '\0';
	if (!str || !*str || !isalpha(str[0])) return NULL;
	digit_index = NULL;

	for(i=0,p=str;i<strlen(str);i++,p++) {
		if (!digit_index) {
			if (isalpha(*p)) continue;
			else if (isdigit(*p)) digit_index = p;
			else return NULL;
		} else {
			if (!isdigit(*p)) return NULL;
		}
	}

	if (digit_index) {
		strncpy(atom_symbol, str, digit_index-str);
		atom_symbol[digit_index-str] = '\0';
	} else strcpy(atom_symbol, str);

	return atom_symbol;
}

int	ParseZmatCenter (char *str, ZmatPtr zmatlist, ZmatPtr zmat)
{
	char	*sym;

	if (IsInt(str)) {
		sscanf(str, "%d", &zmat->an);
		strcpy(zmat->label, GetElemSymbol(zmat->an));
		return 1;
	}
	if ((sym = GetSymbolFromZmatCenter(str))) {	/* alphabet or alphabet+number */
		if ((zmat->an = GetElemNumber(sym)) == 0) {
			if (strcmp(sym, "-") == 0 || strcasecmp(sym, "X") == 0) {	/* dummy */
				zmat->an = -1;
				strcpy(zmat->label, "X");
			} else if (strcasecmp(sym, "BQ") == 0) {	/* ghost */
				strcpy(zmat->label, "Bq");
			} else if (strcasecmp(sym, "XX") == 0) {	/* MOPAC dummy */
				zmat->an = -1;
			} else {
				return setmolerrorno(MERR_UNKNOWN_ATOMSYMBOL);
			}
		} else {
			if (strcmp(sym, str) != 0 && FindZmatByName(zmatlist, str)) {
				return setmolerrorno(MERR_MULTI_DEF_CENTER);
			}
			strcpy(zmat->label, sym);
		}
	} else return setmolerrorno(MERR_ILLEGAL_CENTERNAME);
	return 1;
}

int	ParseZmatRef (char *str, ZmatPtr zmatlist, ZmatPtr zmat, int n)
{
	if (!str || !zmatlist || !zmat || n < 0 || n > 2) return 0;
	if (IsInt(str)) {
		sscanf(str, "%d", &zmat->zdef[n]);
		return 1;
	}
	if (GetSymbolFromZmatCenter(str)) {	/* atomic symbol + number */
		strcpy(zmat->atomname[n+1], str);
		if (!(zmat->zdef[n] = FindZmatByName(zmatlist, str)))
			return setmolerrorno(MERR_UNDEFINED_CENTER);
		else return 1;
	}
	return setmolerrorno(MERR_FORMAT);
}

MolPtr	CreateMolFromZmat (ZmatPtr zmatlist)
{
	MolPtr	mol;
	AtomPtr	atom, conn_atom=NULL;
	BondPtr	bondlist=NULL;
	ResiduePtr	res=NULL;
	ChainPtr	chain=NULL, chainlist=NULL;
	ZmatPtr	zmat;
	int	serno, hetatm;
	int	chain_id = 'A'-1;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	if (!zmatlist) return NULL;
	patom->name[0] = ' ';
	chainlist = NULL;
	chain = NULL;
	res = NULL;
	serno = 0;
	chain_id++;

	strcpy(patom->name, " ");
	strcpy(pres->name, "UNK");
	pres->chain_id = chain_id;
	pres->seq = 0;
	ForEachZmat(zmatlist,zmat) {
		serno++;
		patom->serial = serno;
		strcpy(patom->name, AN2ATOMNAME(zmat->an));
		if (zmat->resname[0]) {
			strcpy(pres->name, zmat->resname);
			pres->seq = zmat->resseq;
		}

		patom->x = zmat->cart[0];
		patom->y = zmat->cart[1];
		patom->z = zmat->cart[2];
		hetatm = 1;
		atom = CreatePdbAtom(&chainlist, &chain, &res, &conn_atom, patom, &bondlist, hetatm);
		atom->an = zmat->an;
		atom->type = zmat->itype;
		atom->refno = GetAtomRefno(zmat->label[0] ? zmat->label : AN2ATOMNAME(atom->an));
	}
	if (!chainlist) return NULL;
	mol = NewMol();
	mol->chain = chainlist;
	if (bondlist) mol->bond = bondlist;

	GetNeighbor(mol);
	CalcBondsWithinResidues(mol);
	FixConnectionTable(mol);

	/* SS list */
	mol->ss = NULL;

	return mol;
}

/* check for redundant variables */
int	IsRedundantIntZmatVar (ZmatPtr zmatlist, int nzmat, ZmatPtr cur_zmat, int vartype)
{
	int	i, j;
	char	*p, *var;
	ZmatPtr	zmat;

	switch (vartype) {
	case Z_VAR_LENGTH: var = cur_zmat->varname[0]; break;
	case Z_VAR_ANGLE:  var = cur_zmat->varname[1]; break;
	case Z_VAR_DIHED:  var = cur_zmat->varname[2]; break;
	default: return 0;
	}
	if (!*var) return 0;
	if (vartype == Z_VAR_ANGLE) {
		if (strcmp(cur_zmat->varname[0], cur_zmat->varname[1]) == 0) return 1;
	} else if (vartype == Z_VAR_DIHED) {
		if (*var == '-' || *var == '+') var++;
		if (strcmp(var, cur_zmat->varname[0]) == 0 ||
		    strcmp(var, cur_zmat->varname[1]) == 0) return 1;
	}

	for(i=0,zmat=zmatlist;i<nzmat;i++,zmat=zmat->next) {
		for(j=0;j<3;j++) {
			p = &zmat->varname[j][0];
			if (j == 3 && (*p == '-' || *p == '+')) p++;
			if (strcmp(p, var) == 0) return 1;
		}
	}
	return 0;
}

int	CountIntZmatVar (ZmatPtr zmatlist)
{
	ZmatPtr	zmat;
	int	nzmat, nvar;

	/* count the number of variables */
	for(zmat=zmatlist,nzmat=nvar=0;zmat;zmat=zmat->next,nzmat++) {
		if (nzmat == 0) continue;
		if (!IsRedundantIntZmatVar(zmatlist, nzmat, zmat, Z_VAR_LENGTH)) nvar++;
		if (nzmat == 1) continue;
		if (!IsRedundantIntZmatVar(zmatlist, nzmat, zmat, Z_VAR_ANGLE )) nvar++;
		if (nzmat == 2) continue;
		if (!IsRedundantIntZmatVar(zmatlist, nzmat, zmat, Z_VAR_DIHED )) nvar++;
	}
	return nvar;
}

int	IsRedundantCartZmatVar (ZmatPtr zmatlist, int nzmat, ZmatPtr cur_zmat, int n)
{
	int	i, j;
	char	*p, *var;
	ZmatPtr	zmat;

	if (!zmatlist || !cur_zmat) return 0;
	if (n < 0 || n > 2) return 1;
	var = cur_zmat->varname[n];
	if (!*var) return 1;
	if (*var == '-' || *var == '+') var++;

	if (n == 1 || n == 2) {
		p = cur_zmat->varname[0];
		if (*p == '-' || *p == '+') p++;
		if (strcmp(p, var) == 0) return 1;

		if (n == 2) {
			p = cur_zmat->varname[1];
			if (*p == '-' || *p == '+') p++;
			if (strcmp(var, cur_zmat->varname[1]) == 0) return 1;
		}
	}

	for(i=0,zmat=zmatlist;i<nzmat;i++,zmat=zmat->next) {
		for(j=0;j<3;j++) {
			p = zmat->varname[j];
			if (*p == '-' || *p == '+') p++;
			if (strcmp(p, var) == 0) return 1;
		}
	}
	return 0;
}

int	CountCartZmatVar (ZmatPtr zmatlist)
{
	ZmatPtr	zmat;
	int	nzmat, nvar;

	/* count the number of variables */
	for(zmat=zmatlist,nvar=0;zmat;zmat=zmat->next) {
		if (!IsRedundantCartZmatVar(zmatlist, nzmat, zmat, 0)) nvar++;
		if (!IsRedundantCartZmatVar(zmatlist, nzmat, zmat, 1)) nvar++;
		if (!IsRedundantCartZmatVar(zmatlist, nzmat, zmat, 2)) nvar++;
	}
	return nvar;
}

int	ConvertZmat (int src_zmat_type, unsigned long src_zmat, int dst_zmat_type, unsigned long *dst_zmat)
{
	GaussZmatPtr	gauss=NULL;
	JaguarZmatPtr	jaguar=NULL;
	BossZmatPtr	boss=NULL;

	if (!src_zmat || !dst_zmat) return 0;
	*dst_zmat=0;

	switch (src_zmat_type) {
	case Z_TYPE_GAUSS:
		gauss = (GaussZmatPtr)src_zmat;
		switch (dst_zmat_type) {
		case Z_TYPE_GAUSS:
			break;
		case Z_TYPE_JAGUAR:
			jaguar = CreateJaguarZmatFromGauss(gauss);
			*dst_zmat = (unsigned long)jaguar;
			break;
		case Z_TYPE_BOSS:
			boss = CreateBossZmatFromGauss(gauss);
			*dst_zmat = (unsigned long)boss;
			break;
		default:
			return 0;
		}
		break;

	case Z_TYPE_JAGUAR:
		jaguar = (JaguarZmatPtr)src_zmat;
		switch (dst_zmat_type) {
		case Z_TYPE_GAUSS:
			gauss = CreateGaussZmatFromJaguar(jaguar);
			*dst_zmat = (unsigned long)gauss;
			break;
		case Z_TYPE_JAGUAR:
			break;
		case Z_TYPE_BOSS:
			boss = CreateBossZmatFromJaguar(jaguar);
			*dst_zmat = (unsigned long)boss;
			break;
		default:
			return 0;
		}
		break;

	case Z_TYPE_BOSS:
		boss = (BossZmatPtr)src_zmat;
		switch (dst_zmat_type) {
		case Z_TYPE_GAUSS:
			gauss = CreateGaussZmatFromBoss(boss);
			*dst_zmat = (unsigned long)gauss;
			break;
		case Z_TYPE_JAGUAR:
			if ((gauss = CreateGaussZmatFromBoss(boss))) break;
			jaguar = CreateJaguarZmatFromGauss(gauss);
			FreeGaussZmat(&gauss);
			*dst_zmat = (unsigned long)jaguar;
			break;
		case Z_TYPE_BOSS:
			break;
		default:
			return 0;
		}
		break;
	}
	
	if (boss) {
		ZmatPtr	zmat, zmatlist;
		zmatlist = LinkZmatList(boss->zmatlist);
		ForEachZmat(zmatlist,zmat) if(zmat->itype == 0) zmat->itype = zmat->an;
		BreakZmatList(boss->zmatlist);
	}

	return (*dst_zmat ? 1 : 0);
}

/* check if zmat->atom in zmatlist is sequential */
int	IsSequentialZmat (ZmatPtr zmatlist)
{
	ZmatPtr	zmat;

	ForEachZmat(zmatlist,zmat) {
		if (zmat->next && zmat->next->atom != zmat->atom+1) return 0;
	}
	return 1;
}

/* update ZmatPtr->zval with ZmatPtr->cart */
void	UpdateZmatValWithCart (ZmatPtr zmat)
{
	ZmatPtr	z, z1, z2, z3;

	ForEachZmat(zmat,z) {
		z1 = SearchZmat(zmat,z->zdef[0]);
		z2 = SearchZmat(zmat,z->zdef[1]);
		z3 = SearchZmat(zmat,z->zdef[2]);
		if (!z1) continue;
		z->zval[0] = DIST(z->cart, z1->cart);
		if (!z2) continue;
		z->zval[1] = GetAngle(z->cart, z1->cart, z2->cart);
		if (!z3) continue;
		z->zval[2] = GetDihedral(z->cart, z1->cart, z2->cart, z3->cart);
	}
}

ZmatPtr	GetZmatOrderedByAtomSerno (ZmatPtr zmatlist)
{
	ZmatPtr	newzmatlist=NULL;
	ZmatPtr	zmat, newzmat, next, prev;
	int	n, i, ierr=0;

	n = 0;
	ForEachZmat(zmatlist,zmat) n++;

	newzmatlist = NULL;
	for(i=1;i<=n;i++) {
		ForEachZmat(zmatlist,zmat) {
			if (zmat->atom == i) break;
		}
		if (!zmat) {
			ierr = 1;
			break;
		}

		prev = zmat->prev;
		next = zmat->next;
		zmat->prev = zmat->next = NULL;
		newzmat = CopyZmat(zmat);
		zmat->prev = prev;
		zmat->next = next;

		if (!newzmat) {
			ierr = 1;
			break;
		}
		EnterZmat(newzmat, &newzmatlist);
	}
	if (ierr) {
		FreeZmat(&newzmatlist);
		return NULL;
	}
	return newzmatlist;
}

/* update ZmatPtr->zval and cart with coordinates from Mol */
int	UpdateZmatWithMolCoord (ZmatPtr zmatlist, MolPtr mol)
{
	AtomPtr	a;
	ZmatPtr	zmat;

	if (!zmatlist || !mol || !mol->chain) return 0;
	if (IsSequentialZmat(zmatlist)) {
		a = FirstAtom(mol->chain);
		ForEachZmat(zmatlist,zmat) {
			if (!a) return 0;
			zmat->cart[0] = a->x;
			zmat->cart[1] = a->y;
			zmat->cart[2] = a->z;
			a = NextAtom(a);
		}
		UpdateZmatValWithCart(zmatlist);
		return 1;
	}

	ForEachZmat(zmatlist,zmat) {
		if (!(a = SearchAtomInChainBySerno(mol->chain,zmat->atom))) return 0;
		zmat->cart[0] = a->x;
		zmat->cart[1] = a->y;
		zmat->cart[2] = a->z;
	}
	UpdateZmatValWithCart(zmatlist);

	return 1;
}

int	UpdateZmatListWithMolCoord (ListPtr zmatlist, MolPtr mol)
{
	ZmatPtr	zmat;
	int	ret_val;

	if (!zmatlist || !mol || !mol->chain) return 0;
	if (!(zmat = LinkZmatList(zmatlist))) return 0;
	ret_val = UpdateZmatWithMolCoord(zmat, mol);
	BreakZmatList(zmatlist);
	return ret_val;
}

/* flags
 * Z_UPDATE_ATOM_ATTR --- update atom attributes using zmatrix
 * coordinates are always updated.
 */
int	UpdateMolWithZmat (ZmatPtr zmatlist, MolPtr mol, int flags)
{
	AtomPtr	a;
	ZmatPtr	zmat;

	if (!zmatlist || !mol || !mol->chain) return 0;
	if (!ZmatToCartesian(zmatlist)) return 0;
	if (IsSequentialZmat(zmatlist)) {
		a = FirstAtom(mol->chain);
		ForEachZmat(zmatlist,zmat) {
			if (!a) return 0;
			a->x = zmat->cart[0];
			a->y = zmat->cart[1];
			a->z = zmat->cart[2];
			if (flags & Z_UPDATE_ATOM_ATTR) {
				if (a->an != zmat->an) {
					a->an = zmat->an;
					a->col = GetElemCPKColor(a->an);
					a->type = zmat->itype;
				}
			}
			a = NextAtom(a);
		}
		return 1;
	}

	ForEachZmat(zmatlist,zmat) {
		if (!(a = SearchAtomInChainBySerno(mol->chain,zmat->atom))) return 0;
		a->x = zmat->cart[0];
		a->y = zmat->cart[1];
		a->z = zmat->cart[2];
		if (flags & Z_UPDATE_ATOM_ATTR) {
			if (a->an != zmat->an) {
				a->an = zmat->an;
				a->col = GetElemCPKColor(a->an);
				a->type = zmat->itype;
			}
		}
	}

	return 1;
}

int	UpdateMolWithZmatList (ListPtr zmatlist, MolPtr mol, int flags)
{
	ZmatPtr	zmat;
	int	retval;

	if (!(zmat = LinkZmatList(zmatlist))) return 0;
	retval = UpdateMolWithZmat(zmat, mol, flags);
	BreakZmatList(zmatlist);
	return retval;
}

/* ---------------------------------------------------------------------- */
/*  Setting and Applying Template Z-matrix                                */
/* ---------------------------------------------------------------------- */

typedef struct {
	int	type;
	GaussZmatPtr	gauss;
	JaguarZmatPtr	jaguar;
	BossZmatPtr	boss;

	ZmatPtr	zmatlist;
	} ZmatTypeTab;

static ZmatTypeTab	_template_zmat, *template_zmat=(&_template_zmat);
#define IS_VALID_ZMAT_TYPE(t)	(t==Z_TYPE_GAUSS||t==Z_TYPE_JAGUAR||t==Z_TYPE_BOSS)

static void	fill_zmat_with_atom (ZmatPtr zmatlist, ChainPtr chain, int user_flags)
{
	ZmatPtr	zmat;
	AtomPtr	atom;

	/* assign zmat->an, atom, cart, itype, ftype, sigma, epsilon */
	atom = FirstAtom(chain);
	ForEachZmat(zmatlist,zmat) {
		zmat->cart[0] = atom->x;
		zmat->cart[1] = atom->y;
		zmat->cart[2] = atom->z;
		zmat->an = atom->an;
		zmat->atom = atom->serno;

		if (atom->type == 0 && (atom->an == ELEM_LP || atom->an == ELEM_GH)) {
			zmat->itype = zmat->ftype = atom->an == ELEM_LP ? 99 : 100;
		} else {
			zmat->itype = zmat->ftype = atom->type;
		}

		if (user_flags & B_USER_CHARGE) zmat->charge = atom->charge;
		if (user_flags & B_USER_SIGMA) zmat->sigma = atom->sigma;
		if (user_flags & B_USER_EPSILON) zmat->epsilon = atom->epsilon;

		zmat->resname[0] = '\0';
		zmat->resseq = 0;

		atom = NextAtom(atom);
	}

	/* assign label, atomname, flags, zval, and varname */
	AssignZmatProp(zmatlist);
}

/* Create a Z-matrix using template.
 *  Input: chain --- list of chains
 *         title --- title of the molecule
 *         zmat_tab->type --- type of Z-matrix to create
 *  Output:
 *         one of zmat_tab->gauss, zmat_tab->jaguar, and zmat_tab->boss
 */

static int	create_zmat_using_template (ChainPtr chain, char *title, int user_flags, ZmatTypeTab *zmat_tab)
{
	if (!chain || !zmat_tab || !template_zmat->zmatlist) return 0;
	if (!title) title = "Z-Matrix";

	zmat_tab->gauss = NULL;
	zmat_tab->jaguar = NULL;
	zmat_tab->boss = NULL;

	switch (zmat_tab->type) {
	case Z_TYPE_GAUSS:
		switch (template_zmat->type) {
		case Z_TYPE_GAUSS:
			zmat_tab->gauss = CopyGaussZmat(template_zmat->gauss);
			zmat_tab->zmatlist = zmat_tab->gauss->zmatlist;
			break;
		case Z_TYPE_JAGUAR:
		case Z_TYPE_BOSS:
			zmat_tab->zmatlist = CopyZmat(template_zmat->zmatlist);
			zmat_tab->gauss = NewGaussZmat();
			zmat_tab->gauss->zmatlist = zmat_tab->zmatlist;
			break;
		}
		fill_zmat_with_atom(zmat_tab->zmatlist, chain, user_flags);
		strcpy(zmat_tab->gauss->title, title);
		break;

	case Z_TYPE_JAGUAR:
		switch (template_zmat->type) {
		case Z_TYPE_JAGUAR:
			zmat_tab->jaguar = CopyJaguarZmat(template_zmat->jaguar);
			zmat_tab->zmatlist = zmat_tab->jaguar->zmatlist;
			break;
		case Z_TYPE_GAUSS:
		case Z_TYPE_BOSS:
			zmat_tab->zmatlist = CopyZmat(template_zmat->zmatlist);
			zmat_tab->jaguar = (JaguarZmat *)calloc(1, sizeof(JaguarZmat));
			zmat_tab->jaguar->zmatlist = zmat_tab->zmatlist;
			break;
		}
		fill_zmat_with_atom(zmat_tab->zmatlist, chain, user_flags);
		break;

	case Z_TYPE_BOSS:
		switch (template_zmat->type) {
		case Z_TYPE_BOSS:
			BreakZmatList(template_zmat->boss->zmatlist);
			zmat_tab->boss = CopyBossZmat(template_zmat->boss);
			FreeBossPar(&zmat_tab->boss->param);
			/*
			FreeList(&zmat_tab->boss->atomtypes);
			FreeList(&zmat_tab->boss->torsiontypes);
			*/

			template_zmat->zmatlist = LinkZmatList(template_zmat->boss->zmatlist);

			zmat_tab->zmatlist = LinkZmatList(zmat_tab->boss->zmatlist);
			fill_zmat_with_atom(zmat_tab->zmatlist, chain, user_flags);
			BreakZmatList(zmat_tab->boss->zmatlist);
			strcpy(zmat_tab->boss->title, title);
			break;
		case Z_TYPE_GAUSS:
		case Z_TYPE_JAGUAR:
			zmat_tab->zmatlist = CopyZmat(template_zmat->zmatlist);
			zmat_tab->boss = (BossZmat *)calloc(1, sizeof(BossZmat));
			zmat_tab->boss->zmatlist = NewList();
			zmat_tab->boss->zmatlist->L_ZMAT = (size_t)zmat_tab->zmatlist;
			fill_zmat_with_atom(zmat_tab->zmatlist, chain, user_flags);
			strcpy(zmat_tab->boss->title, title);
			break;
		}
		break;
	}

/*
if (zmat_tab->gauss) FPrintGaussZmat(stderr, zmat_tab->gauss);
if (zmat_tab->jaguar) FPrintJaguarZmat(stderr, zmat_tab->jaguar);
if (zmat_tab->boss) FPrintBossZmat(stderr, zmat_tab->boss);
*/

	return 1;
}

/**********************************************
 ***   Public Functions                     ***
 **********************************************/

int	SetTemplateZmat (unsigned long zmat, int zmat_type)
{
	if (!zmat) return 0;
	if (!IS_VALID_ZMAT_TYPE(zmat_type)) return 0;

	/* free template */
	if (template_zmat->gauss) FreeGaussZmat(&template_zmat->gauss);
	if (template_zmat->jaguar) FreeJaguarZmat(&template_zmat->jaguar);
	if (template_zmat->boss) FreeBossZmat(&template_zmat->boss);
	template_zmat->zmatlist = NULL;

	template_zmat->type = zmat_type;
	switch (zmat_type) {
	case Z_TYPE_GAUSS:
		template_zmat->gauss = CopyGaussZmat((GaussZmatPtr)zmat);
		template_zmat->zmatlist = template_zmat->gauss->zmatlist;
		break;
	case Z_TYPE_JAGUAR:
		template_zmat->jaguar = CopyJaguarZmat((JaguarZmatPtr)zmat);
		template_zmat->zmatlist = template_zmat->jaguar->zmatlist;
		break;
	case Z_TYPE_BOSS:
		template_zmat->boss = CopyBossZmat((BossZmatPtr)zmat);
		template_zmat->zmatlist = LinkZmatList(template_zmat->boss->zmatlist);
		break;
	}
	return 1;
}

unsigned long	ApplyTemplateZmat (MolPtr mol, int zmat_type, int user_flags)
{
	int	ierr=0, ndummy_atoms=0, ndummy_zmats=0, natoms=0, nzmats=0;
	AtomPtr	atom;
	ZmatPtr	zmat;
	ZmatTypeTab	zmat_tab;

	if (!template_zmat->zmatlist) {
		return (unsigned long)setmolerror("Template Z-matrix has not been set.");
	}
	if (!mol) return (unsigned long)setmolerror("No molecule to apply the template");
	if (!IS_VALID_ZMAT_TYPE(zmat_type)) {
		return (unsigned long)setmolerror("%d is not a valid Z-matrix type.", zmat_type);
	}

	BZERO((void *)&zmat_tab, sizeof(ZmatTypeTab));
	zmat_tab.type = zmat_type;

	/* count the number of atoms in mol */
	natoms = ndummy_atoms = 0;
	natoms = CountAtomInChains(mol->chain);
	for(atom=FirstAtom(mol->chain);atom;atom=NextAtom(atom)) {
		if (atom->an > 0) break;
		ndummy_atoms;
	}

	/* count the number of internal coords in the template */
	nzmats = ndummy_zmats = 0;
	ForEachZmat(template_zmat->zmatlist,zmat) nzmats++;
	ForEachZmat(template_zmat->zmatlist,zmat) {
		if (zmat->an > 0) break;
		ndummy_zmats++;
	}

	if (natoms != nzmats) {
		ierr = 1;
		if (natoms == nzmats-ndummy_zmats) {
			/* add ndummy_zmats of internal coords at the beginning of mol */
			if (PrependDummyAtomsToMol(mol, ndummy_zmats)) ierr = 0;
		}
		if (ierr) {
			return (unsigned long)setmolerror("ERROR: Number of atoms in the template Z-matrix and selected molecule is different.");
		}
	}

	create_zmat_using_template(mol->chain, GetStrValue(mol->header), user_flags, &zmat_tab);

	/* update zmat->val fields using zmat->cart */
	UpdateZmatValWithCart(zmat_tab.zmatlist);

	switch (zmat_type) {
	case Z_TYPE_GAUSS: return (unsigned long)zmat_tab.gauss;
	case Z_TYPE_JAGUAR: return (unsigned long)zmat_tab.jaguar;
	case Z_TYPE_BOSS: return (unsigned long)zmat_tab.boss;
	}

	return 0;
}

int	GetTemplateZmatType (void)
{
	if (!template_zmat->zmatlist) return -1;
	else return template_zmat->type;
}

