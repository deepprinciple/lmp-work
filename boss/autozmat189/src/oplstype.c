#ifdef HAVE_ATOMTYPE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "atomtype/atomtype.h"
#include "atomtype/parameter.h"
#include "atomtype/buildct.h"
#include "atomtype/bondorder.h"
#include "atomtype/aromRing.h"
#include "atomtype/ringPerc.h"
#include "atomtype/util_Type.h"
#include "atomtype/initMol.h"

#include "mol/molecule.h"
#include "error.h"
#include "autozmat.h"
#include "oplstype.h"

/*-------- System-dependent part -------*/
/* on non-UNIX systems, you may have to undefine HAVE_SIGNAL in Makefile */
#ifdef HAVE_SIGNAL
#include <signal.h>
#include <setjmp.h>
static jmp_buf	jump_buffer;
#endif
static int	idletime = 120;	/* 120 sec */
/*--------------------------------------*/

extern int	print_flag;

#define BGN_TORSION_TYPE	601

#define ISDUMMY(an)	(an<0||an>MAXELEM||an==ELEM_GH||an==ELEM_LP)
#define VERBOSITY(n)	(debug_flag > 0 || print_flag >= n)

static int	*orig_serno=NULL;
static int	debug_flag=0;
static OplsTypeTorsion	*opls_typetorsion=NULL;
static int	boss_user_flags = 0;

#define ATOM_SIGMA(serno)	molecule->atom[serno].sigma
#define ATOM_EPSILON(serno)	molecule->atom[serno].epsilon

#define ATOM_ORI_NAME(serno)	molecule->atom[serno].ori_name
#define RES_NAME(serno)		molecule->atom[serno].resname
#define MOLECULE_NAME		molecule->name
#define MOLECULE_NAME_SIZE	sizeof(molecule->name)
#define RADIAN(x)		((x)*0.0174532925193591373)

static Molecule	_molecule, *molecule=(&_molecule);
static RingSet	_rings, *rings=(&_rings);
static DataBase	_DB, *DB=(&_DB);
static Parameter	_PM, *PM=(&_PM);
static Element	elements[MAX_ELEMENTS];
/***************************************************************/

/* maximum number of bonds each atom can possess.
   based on Ruhong Zhou's atomtype code.
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


static void	print_atomtype_func (char *str)
{
	if (print_flag != -1) print_info(str);
}

static int	error_atomtype_func (char *str)
{
	if (print_flag == -1) print_flag = 0;
	print_atomtype_func("ERROR: ");
	print_atomtype_func(str);

#ifdef __i386__
	print_atomtype_func("\n");
	print_atomtype_func("Try '-a' option to bypass this error.\n");
#endif

#ifdef HAVE_SIGNAL
	raise(SIGUSR1);
#else
	exit(1);
#endif

	return 0;
}

/* add a new atomtype to the list */
static int	enter_new_atomtype (ListPtr *atomtypelist, int type, int an, double q, double s, double e, char *amber)
{
	ListPtr	list;

	if (!atomtypelist) return 0;
	ForEachList(*atomtypelist,list) if (list->L_TYPE == type) break;
	if (list) return 1;	/* already exists */
	if (!(list = EnterNewList(atomtypelist))) {
		return seterror("Not enough memory.");
	}
	list->L_TYPE = type;
	list->L_AN = an;
	list->L_CHARGE = q;
	list->L_SIGMA = s;
	list->L_EPSILON = e;
	strcpy(list->L_AMBER, amber);

	return 1;
}

/*--- use List data structure to handle domains ---*/
#define DOMAIN_NATOMS(list)	list->data1
#define DOMAIN_NBONDS(list)	list->data2
#define DOMAIN_ATOM(list)	list->data3
#define DOMAIN_BOND(list)	list->data4
#define DOMAIN_TITLE(list)	list->str1

#define ATOM_SERNO(atom)	atom->data
/*-------------------------------------------------*/

static int	assign_boss_atom_types (BossZmatPtr boss, ListPtr domain)
{
	AtomPtr	a, *atomlist;
	ZmatPtr	zmat, zmatlist;
	ListPtr	list, atomtypes;
	int	i, n;

	if (!boss || !domain) return 0;
	atomlist = (AtomPtr *)DOMAIN_ATOM(domain);

	atomtypes = NULL;
	for(i=0;i<DOMAIN_NATOMS(domain);i++) {
		a = atomlist[i];
		if (a->type <= 0) continue;
		enter_new_atomtype(&atomtypes,
			a->type, a->an, a->charge, a->sigma, a->epsilon, a->ambertype);
	}

	/* copy atom types */
	zmatlist = LinkZmatList(boss->zmatlist);
	ForEachZmat(zmatlist,zmat) {
		for(i=0;i<DOMAIN_NATOMS(domain);i++) {
			a = atomlist[i];
			if (a->serno == zmat->atom) {
				zmat->itype = zmat->ftype = a->type;
				break;
			}
		}
	}
	BreakZmatList(boss->zmatlist);

	/* append newly found atom types to the main list */
	EnterList(atomtypes, &boss->atomtypes);

	/* sort atom type list by type (L_TYPE) in ascending order */
	n = 0;
	ForEachList(boss->atomtypes,list) n++;
	QsortList(boss->atomtypes, 1, n, M_DATA_1);

	return 1;
}

/* check if newly built torsion types can be found in the
 * default torsion parameters.
 */
static int	find_default_torsion_type (double v[])
{
	ListPtr	torsion;
	float	dv0, dv1, dv2, dv3;

	static int	tried=0;
	static int	type = BGN_TORSION_TYPE;

	if (!opls_typetorsion) {
		if (!tried) tried = 1;
		else return 0;

		opls_typetorsion = LoadDefaultOplsTypeTorsion(default_opls_par_file);
	}
	if (!opls_typetorsion) return 0;

	ForEachList(opls_typetorsion->torsionlist,torsion) {
		dv0 = ABS(v[0] - torsion->L_V0);
		dv1 = ABS(v[1] - torsion->L_V1);
		dv2 = ABS(v[2] - torsion->L_V2);
		dv3 = ABS(v[3] - torsion->L_V3);

		if (dv0 < 1.0E-4 && dv1 < 1.0E-4 && dv2 < 1.0E-4 && dv3 < 1.0E-4) break;
	}
	if (torsion) {	/* match found */
		return (int)torsion->L_TYPE;
	} else {
		type++;
		return (type-1);
	}
}

static int	assign_boss_torsion_types (BossZmatPtr boss, MolPtr mol, ListPtr domain,
	ListPtr properlist, ListPtr improperlist)
{
	double	v[4];
	char	*symbol[4];
	int	atom[4], oplstype[4];
	int	i, j, n, serno;
	ListPtr	list, dihedlist, torsion, torsiontypes;
	AtomPtr	a, *atomlist;

	if (!boss || !domain) return 0;
	atomlist = (AtomPtr *)DOMAIN_ATOM(domain);

	torsiontypes = NULL;
	for(n=0;n<molecule->n_proptors+molecule->n_imprtors;n++) {
		Torsion	*ttype;

		ttype = &molecule->torsion_list[n];
		for(i=0;i<4;i++) {
			oplstype[i] = OPLSTYPE(ttype->atom[i]);
			atom[i] = 0;

			for(j=0;j<DOMAIN_NATOMS(domain);j++) {
				a = atomlist[j];
				if (ATOM_SERNO(a) == ttype->atom[i]) {
					atom[i] = a->serno;	/* original serial number */
					break;
				}
			}

			v[i] = ttype->v[i];
			symbol[i] = ttype->symbol[i];
			if (strlen(symbol[i]) == 0 || strcmp(symbol[i], " ") == 0 ||
				strcmp(symbol[i], "  ") == 0) {
				symbol[i] = "??";
			}
		}

		list = LookUpDihedList(properlist,atom[0],atom[1],atom[2],atom[3]);
		if (!list) {
			list = LookUpImproperList(improperlist,properlist,atom[0],atom[1],atom[2],atom[3]);
		}

		if (!list || !(dihedlist = (ListPtr)list->L_DATA_7)) continue;

		ForEachList(torsiontypes,torsion) {
			if (ABS(torsion->L_V0 - v[0]) < 1.0E-4 && 
			    ABS(torsion->L_V1 - v[1]) < 1.0E-4 &&
			    ABS(torsion->L_V2 - v[2]) < 1.0E-4 &&
			    ABS(torsion->L_V3 - v[3]) < 1.0E-4) break;
		}

		/* new torsion type */
		if (!torsion) {
			if (!(torsion = EnterNewList(&torsiontypes))) break;
			torsion->L_TYPE = find_default_torsion_type(v);
			torsion->L_V0 = v[0];
			torsion->L_V1 = v[1];
			torsion->L_V2 = v[2];
			torsion->L_V3 = v[3];
			sprintf(torsion->L_COMMENT, "%-2s-%-2s-%-2s-%-2s (%03d-%03d-%03d-%03d)",
				symbol[0], symbol[1], symbol[2], symbol[3],
				oplstype[0], oplstype[1], oplstype[2], oplstype[3]);
		}

		dihedlist->L_ITYPE = dihedlist->L_FTYPE = torsion->L_TYPE;
	}

	/* check if all the torsions have been assigned a type */
	ForEachList(boss->vardihed,list) {
		Torsion	_ttype, *ttype=(&_ttype);

		if (list->L_ITYPE > 0) continue;	/* type has been assigned */

		BZERO(ttype, sizeof(Torsion));
		for(i=0;i<4;i++) {
			switch (i) {
			case 0: n = list->L_ATOM1; break;
			case 1: n = list->L_ATOM2; break;
			case 2: n = list->L_ATOM3; break;
			case 3: n = list->L_ATOM4; break;
			}

			for(j=0;j<DOMAIN_NATOMS(domain);j++) {
				a = atomlist[j];
				if (a->serno == n) {
					break;
				}
			}
			if (j == DOMAIN_NATOMS(domain)) continue;

			ttype->atom[i] = n;
			ttype->proper = list->L_IMPROPER ? -1 : 1;
			ttype->an[i] = a->an;

			strcpy(ttype->symbol[i], OPLSNAME(n));
		}

		/* try to match ttype to a torsion type in the parameter file.
		 * match_torsion_type() will fill ttype->v if successful
		 */
		if (match_torsion_type(ttype, PM, molecule, 0) == -1) {	/* no match */
			continue;
		}

		ForEachList(torsiontypes,torsion) {
			if (torsion->L_V0 == ttype->v[0] && 
			    torsion->L_V1 == ttype->v[1] &&
			    torsion->L_V2 == ttype->v[2] &&
			    torsion->L_V3 == ttype->v[3]) break;
		}

		/* new torsion type */
		if (!torsion) {
			if (!(torsion = EnterNewList(&torsiontypes))) break;
			torsion->L_TYPE = find_default_torsion_type(ttype->v);
			torsion->L_V0 = ttype->v[0];
			torsion->L_V1 = ttype->v[1];
			torsion->L_V2 = ttype->v[2];
			torsion->L_V3 = ttype->v[3];

			for(i=0;i<4;i++) {
				oplstype[i] = OPLSTYPE(ttype->atom[i]);
				symbol[i] = ttype->symbol[i];
				if (strlen(symbol[i]) == 0 || strcmp(symbol[i], " ") == 0 ||
					strcmp(symbol[i], "  ") == 0) {
					symbol[i] = "??";
				}
			}

			sprintf(torsion->L_COMMENT, "%-2s-%-2s-%-2s-%-2s (%03d-%03d-%03d-%03d)",
				symbol[0], symbol[1], symbol[2], symbol[3],
				oplstype[0], oplstype[1], oplstype[2], oplstype[3]);

		}

		list->L_ITYPE = list->L_FTYPE = torsion->L_TYPE;
	}

	EnterList(torsiontypes, &boss->torsiontypes);

	/* sort atom type list by type (L_TYPE) in ascending order */
	n = 0;
	ForEachList(boss->torsiontypes,list) n++;
	QsortList(boss->torsiontypes, 1, n, M_DATA_1);

	return 1;
}

static void	cleanup_atomtype (void)
{
	cleanup_rings(rings);
	release_mol(molecule);
	release_DB(DB);
	release_PM(PM);
}

#ifdef HAVE_SIGNAL
static void	sig_user_handler (int sig)
{
	sig = sig;
	longjmp(jump_buffer, 1);
}

static void	sig_alarm_handler (int sig)
{
	sig = sig;
	fprintf(stderr, "Timed out.\n");
	longjmp(jump_buffer, 1);
}
#endif

/*
 * process_atomtype --- driver function for the atomtype library
 *
 * domain --- domain information
 * errmsg --- message to be issued on error condition
 *
 */

static int	process_atomtype (ListPtr domain, char *errmsg)
{
	AtomPtr	*atomlist, atom;
	BondPtr	*bondlist, bond;
	int	i, j, n, serno, an, ia1, ia2;
	char	symbol[4];

#ifdef HAVE_SIGNAL
/* set jump position */
if (setjmp(jump_buffer)) {
	if (!(autozmat_flags & M_QUIET)) print_atomtype_func(errmsg);
	return 0;
}
signal(SIGUSR1, sig_user_handler);
signal(SIGALRM, sig_alarm_handler);
alarm(idletime);
#endif

	if (!domain) return 0;
	atomlist = (AtomPtr *)DOMAIN_ATOM(domain);
	bondlist = (BondPtr *)DOMAIN_BOND(domain);

	/* -------------------------------------------------------------- */
	/* Initialize the atomtype library and pass the connection table. */
	/* -------------------------------------------------------------- */

	cleanup_atomtype();

	/* initialize and allocate molecule data */
	initialize_mol(molecule);
	TOTALATOMS = DOMAIN_NATOMS(domain);
	allocate_mol(molecule);

	/* title */
	n = MOLECULE_NAME_SIZE;
	strncpy(MOLECULE_NAME, DOMAIN_TITLE(domain), n-1);
	MOLECULE_NAME[n-1] = '\0';

	/* put the data into the atomtype main code */
	for(i=0,serno=0;i<DOMAIN_NATOMS(domain);i++) {
		atom = atomlist[i];
		an = atom->an;
		if (ISDUMMY(an)) continue;
		serno++;
		ATOM_SERNO(atom) = serno;

		if (an == ELEM_X) {	/* an=0, growing-point */
			an = 1;	/* becomes hydrogen */
		}
		strcpy(symbol, GetElemSymbol(an));
		if (symbol[1] == ' ') symbol[1] = '\0';

		ATOMICNUMBER(serno) = an;
		strcpy(ATOMNAME(serno), symbol);
		strcpy(ATOM_ORI_NAME(serno), symbol);
		strcpy(RES_NAME(serno), GetResName(atom->residue->refno));

		if (autozmat_flags & M_ASSIGN_XYZ) {
			ATOMCOORD(serno).x = atom->x;
			ATOMCOORD(serno).y = atom->y;
			ATOMCOORD(serno).z = atom->z;
		}
	}

	/* initialize connection table info */
	TOTALBONDS = 0;
	for(i=1;i<=TOTALATOMS;i++) {
		ATOMNCONN(i) = 0;
		for(j=0;j<MAX_BONDS;j++){
			CONNECTION(i,j) = 0;
			BO(i,j) = 0;
		}
	}

	/* put the bond data into the atomtype main code */
	for(i=0,TOTALBONDS=0;i<DOMAIN_NBONDS(domain);i++) {
		bond = bondlist[i];
		if (ISDUMMY(bond->atom1->an) || ISDUMMY(bond->atom2->an)) continue;
		ia1 = ATOM_SERNO(bond->atom1);
		ia2 = ATOM_SERNO(bond->atom2);

		CONNSTART(TOTALBONDS) = ia1;
		CONNEND(TOTALBONDS) = ia2;

		BO(ia1, ATOMNCONN(ia1)) = 0;	/* must be set to zero */
		CONNECTION(ia1, ATOMNCONN(ia1)) = ia2;
		ATOMNCONN(ia1)++;

		BO(ia2, ATOMNCONN(ia2)) = 0;	/* must be set to zero */
		CONNECTION(ia2, ATOMNCONN(ia2)) = ia1;
		ATOMNCONN(ia2)++;

		/* increment the total number of bonds */
		TOTALBONDS++;
	}

if (VERBOSITY(2)) {
	PRINT("Initial connection table:\n");
	print_connection_table(molecule);
}

	/* prepare for setting up connection table */
	check_residue_name(molecule);
	read_element_file(elements);
	pre_treatment(molecule,elements);
	assign_radii(molecule,elements);

if (VERBOSITY(1)) {
	PRINT("TOTALATOMS = %d, TOTALBONDS = %d\n", TOTALATOMS, TOTALBONDS);
	for(i=1;i<=TOTALATOMS;i++) {
		PRINT("   Atom %4d: nconn = %d, max_bonds = %d\n",
			i, molecule->atom[i].nconnects, molecule->atom[i].max_bonds);
	}
}

	/* sort connection table */
	qsort(molecule->connections,TOTALBONDS,sizeof(connect_type), QSORT_PROTO sort_connections);

	/* check connection table */
	check_connection_table(molecule);

if (VERBOSITY(2)) print_connection_table(molecule);

	/* search for Smallest Set of Smallest Rings (SSSR) */
	search_SSSR(molecule, rings);

	/* find aromatic rings then */
	find_aromatic_rings(molecule, rings);

	/* assign bond orders and figure out atom valences */
	assign_bond_order(molecule, rings);

	/* assign OPLSAA atom type */
	assign_atom_type(molecule, rings, DB);

	/* assign parameters for stretching, bending, and torsions */
	if (DO_TORSION_TYPING(autozmat_flags)) {
		assign_parameters(molecule, rings, PM);
if (VERBOSITY(2)) report_parameters(molecule);
	}

	/* -------------------- Done! ---------------------------- */

	/* retrieve atom types */
	for(i=0,n=0;i<DOMAIN_NATOMS(domain);i++) {
		atom = atomlist[i];
		if (ISDUMMY(atom->an)) {	/* dummy atom */
			atom->type = -1;
			continue;
		}

		n++;
		atom->type = OPLSTYPE(n);
		strcpy(atom->ambertype, OPLSNAME(n));

if (!(boss_user_flags & B_USER_CHARGE)) atom->charge = ATOMCHARGE(n);
if (!(boss_user_flags & B_USER_SIGMA)) atom->sigma = ATOM_SIGMA(n);
if (!(boss_user_flags & B_USER_EPSILON)) atom->epsilon = ATOM_EPSILON(n);

	}

	/* -- NOTE ---
	 * The atomtype library must not be cleaned up at this point, because molecule
	 * database is still going to be used in the calling routine.
	 */

	alarm(0U);
	return 1;
}

static ListPtr	enter_new_domain_from_list (ListPtr list, ListPtr *top)
{
	ListPtr	domain, l;
	AtomPtr	atom, *atomlist;
	BondPtr	bond, *bondlist;
	int	i, natoms, nbonds;

	if (!top || !list) return NULL;
	if (!(domain = EnterNewList(top))) return NULL;

	/* count the number of atoms and bonds */
	natoms = nbonds = 0;
	for(l=list;l;l=l->next) {
		atom = (AtomPtr)DOMAIN_ATOM(l);
		if (ISDUMMY(atom->an)) continue;
		natoms++;
		ATOM_SERNO(atom) = natoms;

		for(i=0;i<atom->nneighbors;i++) {
			bond = atom->nbbond[i];
			bond->flags &= ~B_MARKED;
			if (bond->flags & B_SELECTED) continue;
			if (!ISDUMMY(bond->atom1->an) && !ISDUMMY(bond->atom2->an)) nbonds++;
		}
	}

	if (natoms == 0) return 0;
	if (nbonds % 2 != 0) {
		PRINT("ERROR in the number of bonds: %d\n", nbonds);
		return NULL;
	}

	nbonds = nbonds/2;	/* bonds were counted twice */

	DOMAIN_NATOMS(domain) = natoms;
	DOMAIN_NBONDS(domain) = nbonds;
	atomlist = (AtomPtr *)calloc(1, natoms*sizeof(AtomPtr));
	if (nbonds > 0) {
		bondlist = (BondPtr *)calloc(1, nbonds*sizeof(BondPtr));
	} else bondlist = NULL;
	DOMAIN_ATOM(domain) = (size_t)atomlist;
	DOMAIN_BOND(domain) = (size_t)bondlist;

	for(l=list,natoms=0,nbonds=0;l;l=l->next) {
		atom = (AtomPtr)DOMAIN_ATOM(l);
		if (ISDUMMY(atom->an)) continue;
		atomlist[natoms] = atom;
		natoms++;
		for(i=0;i<atom->nneighbors;i++) {
			bond = atom->nbbond[i];
			if (ISDUMMY(bond->atom1->an) || ISDUMMY(bond->atom2->an)) continue;
			if (bond->flags & B_MARKED) continue;
			if (bond->flags & B_SELECTED) continue;
			bond->flags |= B_MARKED;
			bondlist[nbonds] = bond;
			nbonds++;
		}
	}

	return domain;
}

static AtomPtr	search_domain_by_serno (ListPtr domain, int serno)
{
	int	i;
	AtomPtr	atom, *atomlist;

	atomlist = (AtomPtr *)DOMAIN_ATOM(domain);
	for(i=0;i<DOMAIN_NATOMS(domain);i++) {
		atom = atomlist[i];
		if (atom->serno == serno) return atom;
	}
	return NULL;
}

static void	free_domain (ListPtr *top)
{
	ListPtr	p, tmp;

	if (!top || !*top) return;
	p = *top;
	while (p) {
		tmp = p->next;
		if (DOMAIN_ATOM(p)) free((void *)DOMAIN_ATOM(p));
		if (DOMAIN_BOND(p)) free((void *)DOMAIN_BOND(p));
		free((void *)p);
		p = tmp;
	}
	*top = NULL;
}

static void	print_domain (ListPtr domain)
{
	AtomPtr	atom, nbatom, *atomlist;
	BondPtr	nbbond;
	int	i, j, natoms, nconn, conn[MMOD_MAX_BONDS], type[MMOD_MAX_BONDS];

	atomlist = (AtomPtr *)DOMAIN_ATOM(domain);
	natoms = DOMAIN_NATOMS(domain);

	PRINT("%6d  %s\n", natoms, "TITLE");
	for(i=0;i<natoms;i++) {
		atom = atomlist[i];
		for(j=0;j<MMOD_MAX_BONDS;j++) conn[j] = type[j] = 0;
		for(j=0,nconn=0;j<atom->nneighbors && nconn < MMOD_MAX_BONDS;j++) {
			nbatom = atom->nbatom[j];
			nbbond = atom->nbbond[j];
			if (ISDUMMY(nbatom->an)) continue;
			conn[nconn] = ATOM_SERNO(nbatom);
			type[nconn] = GetMModelBondType(nbbond->type);
			nconn++;
		}

		PRINT(" %3d", GetMModelAtomType(atom));
		for(j=0;j<MMOD_MAX_BONDS;j++) PRINT(" %5d %d", conn[j], type[j]);
		PRINT(" % 11.6f % 11.6f % 11.6f\n", atom->x, atom->y, atom->z);
	}
}

/*
 * domain         --- list of atoms which forms a domain
 * start, destiny --- start and destiny of searching
 * atom           --- current atom
 */
#define MAX_DOMAIN_SEARCH_DEPTH	100000
static int	domain_search_depth=0;
static int	_get_domain_list (ListPtr *domain, AtomPtr start, AtomPtr destiny, AtomPtr atom)
{
	AtomPtr	nbatom;
	BondPtr	nbbond;
	ListPtr	list;
	int	i;

	domain_search_depth++;
	if (domain_search_depth > MAX_DOMAIN_SEARCH_DEPTH) {	/* recursion too deep */
		return 0;
	}

	if (atom->flags & A_MARKED) return 0;
	atom->flags |= A_MARKED;
	if (!(list = EnterNewList(domain))) return 0;
	DOMAIN_ATOM(list) = (size_t)atom;

	for(i=0;i<atom->nneighbors;i++) {
		nbatom = atom->nbatom[i];
		nbbond = atom->nbbond[i];
		if (ISDUMMY(nbatom->an)) continue;
		if (nbbond->flags & B_SELECTED) continue;
		if ((nbatom == destiny && atom != start) ||
			(!(nbatom->flags & A_MARKED) && !_get_domain_list(domain, start, destiny, nbatom))) {
			if (*domain) FreeList(domain);
			*domain = (List *)NULL;
			return 0;
		}
	}
	return 1;
}

/* before calling this routine first time, turn off A_MARKED bit in atom->flags.
 */
static List	*get_domain_list (AtomPtr start, AtomPtr destiny)
{
	ListPtr	list=NULL;

	if (!start || start == destiny) return 0;
	if (destiny) destiny->flags |= A_MARKED;
	domain_search_depth = 0;
	_get_domain_list(&list, start, destiny, start);
	return list;
}

static void	prepare_torsion (BossZmatPtr boss, MolPtr mol,
	ListPtr *properlist_ret, ListPtr *improperlist_ret, ListPtr *last_vardihed_ret)
{
	ZmatPtr	zmat,zmatlist;
	ListPtr	list,dihedlist;
	ListPtr	properlist=NULL,improperlist=NULL,lastlist=NULL;

	/* zero return values */
	*properlist_ret = NULL;
	*improperlist_ret = NULL;
	*last_vardihed_ret = NULL;

	/* We need to assign L_ATOM2, L_ATOM3, L_ATOM4 in variable diherals,
	   since LookUpDihedList depends on them.
	*/
	zmatlist = LinkZmatList(boss->zmatlist);
	ForEachList(boss->vardihed,list) {
		if (!(zmat = SearchZmat(zmatlist,list->L_CENTER))) continue;
		list->L_ATOM2 = zmat->zdef[0];
		list->L_ATOM3 = zmat->zdef[1];
		list->L_ATOM4 = zmat->zdef[2];
	}
	BreakZmatList(boss->zmatlist);

	/*** link additional dihedrals to variable dihedrals.
	   we'll break them later.
	***/

	lastlist = NULL;
	for(list=boss->vardihed;list && list->next;list=list->next);
	lastlist = list;
	EnterList(boss->adddihed, &boss->vardihed);

	/* set improper torsion flag */
	SetImproperTorsionFlags(boss->vardihed, mol->chain, mol->bond);

	/* generate proper and improper torsion lists */
	properlist = NULL;
	improperlist = NULL;
	ForEachList(boss->vardihed,list) {
		/* set the torstion type to 0 */
		list->L_ITYPE = list->L_FTYPE = 0;

		dihedlist = NewList();
		BCOPY(list, dihedlist, sizeof(List));
		dihedlist->next = NULL;
		/*
		   L_DATA_7 is free. use it for storing the pointer of this list.
		*/
		dihedlist->L_DATA_7 = (size_t)list;	/* data7 */
		/* default torsion type=0 (undecided) */
		dihedlist->L_ITYPE = dihedlist->L_FTYPE = 0;

		/* append to properlist or improperlist */
		EnterList(dihedlist, list->L_IMPROPER ? &improperlist : &properlist);
	}
	
	*properlist_ret = properlist;
	*improperlist_ret = improperlist;
	*last_vardihed_ret = lastlist;
}

#define ASSIGN_OPLS_ERR_MSG	"\
\n\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
!!! Error occurred while trying to assign atomtypes.                !!!\n\
!!! Atomtypes will not be assigned, but everything else will be OK. !!!\n\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
\n"

int	assign_opls_parameters (BossZmatPtr boss)
{
	MolPtr	mol;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a, *atomlist;
	BondPtr	b;
	ListPtr	list, newlist, domain, domainlist=NULL;
	int	i, type, natoms, ndomains, ierr=0;
	ZmatPtr	zmat, zmatlist;

	ListPtr	properlist=NULL, improperlist=NULL, last_vardihed=NULL;
	ListPtr	saved_atomtypes=NULL, saved_torsiontypes=NULL;
	ListPtr	newatomtypes=NULL;

	set_atomtype_error_func(error_atomtype_func);
	set_atomtype_print_func(print_atomtype_func);

	if (!boss) return 0;
	if (!(mol = CreateMolFromBossZmat(boss))) return 0;
	boss_user_flags = boss->flags & (B_USER_CHARGE|B_USER_SIGMA|B_USER_EPSILON);

	if (DO_ATOM_TYPING(autozmat_flags)) {
		if (boss->atomtypes) saved_atomtypes = boss->atomtypes;
		boss->atomtypes = NULL;
	}
	
	if (DO_TORSION_TYPING(autozmat_flags)) {
		if (boss->torsiontypes) saved_torsiontypes = boss->torsiontypes;
		boss->torsiontypes = NULL;
	}

	init_mol_info(mol, boss->title, 0, 0);
	FixConnectionTable(mol);
	ForEachBond(mol->bond,b) b->flags &= ~B_SELECTED;

	/* if a molecule is too big, process atomtype residue by residue */
	if ((natoms = CountAtomInChains(mol->chain)) > 1000) {
		/*
		PRINT("INFO: Molecule is too big (%d atoms) for atomtyping.\n", natoms);
		PRINT("      Atomtyping will be processed on residue by residue basis.\n");
		*/
		ForEachBond(mol->bond,b) {
			/*
			if (b->atom1->residue != b->atom2->residue) b->flags |= B_SELECTED;
			*/
		}
	}


	/* mark dummy atoms so that they are not included in domains */
	ForEachChainResAtom(mol->chain,c,r,a) {
		a->flags &= ~A_MARKED;
		ATOM_SERNO(a) = -1;
		if (ISDUMMY(a->an)) {
			a->flags |= A_MARKED;
			a->type = -1;
		}
	}

	/* build domains */
	ndomains = 0;
	ForEachChainResAtom(mol->chain,c,r,a) {
		if (a->flags & A_MARKED) continue;

		if (!(list = get_domain_list(a, (AtomPtr)NULL))) break;
		domain = enter_new_domain_from_list(list, &domainlist);
		FreeList(&list);
		ndomains++;

if (VERBOSITY(2)) {
	PRINT("Domain #%d: %d atoms\n", ndomains, DOMAIN_NATOMS(domain));
	print_domain(domain);
}

	}

	/* prepare dihedrals for torsion type assignment (ie, assign_boss_torsion_types) */
	if (DO_TORSION_TYPING(autozmat_flags)) {
		prepare_torsion(boss, mol, &properlist, &improperlist, &last_vardihed);
	}

	/* run atomtype for each domain */
	ndomains = 0;
	ForEachList(domainlist,domain) {
		ndomains++;
		natoms = DOMAIN_NATOMS(domain);

if (VERBOSITY(2)) {
	PRINT("Processing domain #%d (%d atoms)...\n", ndomains, natoms);
}

		if (natoms > MAX_DOMAIN_ATOMS) {
			PRINT("INFO: It may take a while to assign atom types to this domain due to\n");
			PRINT("      its size. If time matters, use '-a' option to skip atomtyping.\n");
		}

		if (!process_atomtype(domain, ASSIGN_OPLS_ERR_MSG)) {
			ierr = 1;
			break;
		}

		/* assign atom types */
		if (DO_ATOM_TYPING(autozmat_flags)) {
			assign_boss_atom_types(boss, domain);
		}

		/* assign torsion types */
		if (natoms > 3 && DO_TORSION_TYPING(autozmat_flags)) {
			assign_boss_torsion_types(boss, mol, domain, properlist, improperlist);
		}
	}

	/***** clean up torsions *****/

	if (domainlist) FreeList(&domainlist);
	if (properlist) FreeList(&properlist);
	if (improperlist) FreeList(&improperlist);

	/* break boss->vardihed into the original state */
	if (DO_TORSION_TYPING(autozmat_flags)) {
		if (last_vardihed) last_vardihed->next = NULL;
		else boss->vardihed = NULL;
	}

	/* clean-up */
	cleanup_atomtype();
	
	if (ierr) {	/* use the original atomtypes and torsiontypes list */
		if (boss->atomtypes) FreeList(&boss->atomtypes);
		if (boss->torsiontypes) FreeList(&boss->torsiontypes);
		boss->atomtypes = saved_atomtypes;
		boss->torsiontypes = saved_torsiontypes;
	} else {
		if (DO_ATOM_TYPING(autozmat_flags)) {
			if (saved_atomtypes) FreeList(&saved_atomtypes);
		}
		if (DO_TORSION_TYPING(autozmat_flags)) {
			if (saved_torsiontypes) FreeList(&saved_torsiontypes);
		}

if (boss_user_flags) {
		/* make a new atom type for each atom */
		newatomtypes = NULL;
		type = 801;
		ForEachChainResAtom(mol->chain,c,r,a) {
			if (a->type < 0) continue;
			ForEachList(boss->atomtypes,list) if (list->L_TYPE == a->type) break;
			if (!list) {	/* hmm... this is strange */
				continue;
			}

			if (!(newlist = EnterNewList(&newatomtypes))) {
				seterror("Not enough memory.");
				ierr =1;
				break;
			}

			newlist->L_TYPE = type;
			newlist->L_AN   = list->L_AN;
			strcpy(newlist->L_AMBER, list->L_AMBER);

			newlist->L_CHARGE  = (boss_user_flags & B_USER_CHARGE)  ? a->charge  : list->L_CHARGE;
			newlist->L_SIGMA   = (boss_user_flags & B_USER_SIGMA)   ? a->sigma   : list->L_SIGMA;
			newlist->L_EPSILON = (boss_user_flags & B_USER_EPSILON) ? a->epsilon : list->L_EPSILON;

			a->type = type;
			zmatlist = LinkZmatList(boss->zmatlist);
			ForEachZmat(zmatlist,zmat) {
				if (a->serno == zmat->atom) {
					zmat->itype = zmat->ftype = a->type;
					break;
				}
			}
			BreakZmatList(boss->zmatlist);

			type++;
		}

		if (newatomtypes) {
			FreeList(&boss->atomtypes);
			boss->atomtypes = newatomtypes;
		}
}

	}

	return !ierr;
}

static void	print_opls_force_fields ()
{
	int	i, j, n;
	int	atom[4];
	float	r0, rk, t0, tk, v[4];

	/* assign parameters for stretching, bending, and torsions */
	assign_parameters(molecule, rings, PM);

	if (VERBOSITY(2)) report_parameters(molecule);

	for(n=0;n<molecule->n_stretchs;n++) {
		atom[0] = molecule->stretch_list[n].atom[0];
		atom[1] = molecule->stretch_list[n].atom[1];
		r0 = molecule->stretch_list[n].r0;
		rk = molecule->stretch_list[n].k;
	}

	for(n=0;n<molecule->n_bendings;n++) {
		for(i=0;i<3;i++) atom[i] = molecule->bending_list[n].atom[i];
		t0 = RADIAN(molecule->bending_list[n].theta0);
		tk = molecule->bending_list[n].k;
	}

	for(n=0;n<molecule->n_proptors;n++) {
		for(i=0;i<4;i++) {
			atom[i] = molecule->torsion_list[n].atom[i];
			v[i] = molecule->torsion_list[n].v[i];
		}
	}

	for(n=0;n<molecule->n_imprtors;n++) {
		for(i=0;i<4;i++) {
			j = n+molecule->n_proptors;
			atom[i] = molecule->torsion_list[j].atom[i];
			v[i] = molecule->torsion_list[j].v[i];
		}
	}
}

#define BUILD_CONN_ERR_MSG "\
\n\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
!!! Error occurred while trying to build connectivity information.  !!!\n\
!!! Assigned bond orders (especially aromatic rings) may be wrong.  !!!\n\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
\n"

int	build_connectivity (MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	atom, atom1, atom2, *atomlist;
	BondPtr	bond, bondlist;
	ListPtr	list, domain, domainlist=NULL;

	int	i, j, k, n, an, ia1, ia2, natoms, ndomains;
	int	ierr = 0, saved_flags=0;

	set_atomtype_error_func(error_atomtype_func);
	set_atomtype_print_func(print_atomtype_func);

	if (!mol) return 0;

	init_mol_info(mol, GetStrValue(mol->header), 0, 0);
	FixConnectionTable(mol);
	ForEachBond(mol->bond,bond) bond->flags &= ~B_SELECTED;

	/* if a molecule is too big, process atomtype residue by residue */
	if ((natoms = CountAtomInChains(mol->chain)) > 1000) {
		/*
		PRINT("INFO: Molecule is too big (%d atoms) for atomtyping.\n", natoms);
		PRINT("      Atomtyping will be processed on residue by residue basis.\n");
		*/
		ForEachBond(mol->bond,bond) {
			/*
			if (bond->atom1->residue != bond->atom2->residue) bond->flags |= B_SELECTED;
			*/
		}
	}

	/* mark dummy atoms so that they are not included in domains */
	ForEachChainResAtom(mol->chain,c,r,atom) {
		atom->flags &= ~A_MARKED;
		ATOM_SERNO(atom) = -1;
		if (ISDUMMY(atom->an)) {
			atom->flags |= A_MARKED;
			atom->type = -1;
		}
	}

	/* build domains */
	n = 0;
	ForEachChainResAtom(mol->chain,c,r,atom) {
		if (atom->flags & A_MARKED) continue;
		if (!(list = get_domain_list(atom, (AtomPtr)NULL))) break;
		domain = enter_new_domain_from_list(list, &domainlist);
		n++;

if (VERBOSITY(2)) {
	PRINT("Domain #%d: %d atoms\n", n, DOMAIN_NATOMS(domain));
	print_domain(domain);
}

		FreeList(&list);
	}

	/* run atomtype for each domain */
	saved_flags = autozmat_flags;
	autozmat_flags = M_ASSIGN_XYZ|M_NO_TORSION_TYPING;
	autozmat_flags &= ~M_NO_ATOM_TYPING;

	bondlist = NULL;
	ndomains = 0;
	ForEachList(domainlist,domain) {
		ndomains++;
		natoms = DOMAIN_NATOMS(domain);

if (VERBOSITY(2)) {
	PRINT("Processing domain #%d (%d atoms)...\n", ndomains, natoms);
}

		if (natoms > MAX_DOMAIN_ATOMS) {
			PRINT("INFO: It may take a while to assign atom types to this domain due to\n");
			PRINT("      its size. If time matters, use '-a' option to skip atomtyping.\n");
		}

		if (!process_atomtype(domain, ASSIGN_OPLS_ERR_MSG)) {
			ierr = 1;
			break;
		}

		/*
		 * Rebuild the connectivity.
		 */
		atomlist = (AtomPtr *)DOMAIN_ATOM(domain);
		ia1 = 0;
		for(i=0;i<DOMAIN_NATOMS(domain);i++) {
			atom1 = atomlist[i];
			if (ISDUMMY(atom1->an)) {	/* dummy atom */
				atom1->type = -1;
				continue;
			}

			ia1++;

			for(j=0;j<ATOMNCONN(ia1);j++) {
				ia2 = CONNECTION(ia1, j);

				atom2 = NULL;
				for(k=0;k<DOMAIN_NATOMS(domain);k++) {
					atom = atomlist[k];
					if (ATOM_SERNO(atom) == ia2) {
						atom2 = atom;
						break;
					}
				}
				if (!atom2) continue;

				/* Check if bond already exists */
				if (GetAtomBonded(bondlist, atom1, atom2)) continue;
				if (!(bond = EnterNewBond(&bondlist))) {
					ierr = 1;
					break;
				}
				bond->type = BO(ia1, j);
				bond->atom1 = atom1;
				bond->atom2 = atom2;
			}
			if (ierr) break;
		}
	}
	
	/* replace the old bond list by new one */
	if (!ierr && bondlist) {
		FreeBond(&mol->bond);
		mol->bond = bondlist;
		GetNeighbor(mol);
		mol->param_type = PARAM_OPLS;
	}

	if (domainlist) FreeList(&domainlist);

	/* clean-up */
	cleanup_atomtype();

	autozmat_flags = saved_flags;
	return !ierr;
}
void	set_idle_time (int tsec) {idletime = tsec;}
void	set_opls_debug_level (int level) {debug_flag = level;}
void	set_atomtype_debug_level (int level) {print_flag = level;}

#else
/* atomtype library is not used.
 * just give dummy definitions.
 */
#include <stdio.h>
#include "mol/molecule.h"
int	assign_opls_parameters (BossZmatPtr boss) {return 1;}
void	set_idle_time (int tsec) {tsec=tsec;}
void	set_opls_debug_level (int level) {level=level;}
void	set_atomtype_debug_level (int level) {level=level;}

#endif	/* HAVE_ATOMTYPE */


