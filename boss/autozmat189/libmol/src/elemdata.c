#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

#define IsValidElemNo(an)	(an >= 0 && an < MAXELEM)
#define BCOPY(s1, s2, n)	memmove((char *)s2, (char *)s1, n)

static ElemData	elemdata[] = ELEM_DATA;

/*#######################################################################*/

void	FPrintElemData (FILE *fp)
{
	int	i;
	ElemData	*e;

	for(i=0,e=elemdata;i<MAXELEM;i++,e++) {
		fprintf(fp, "{\"%-3s\", \"%-12s\", %9.5f, %9.5f, %4.2f, %4.2f}\n",
			e->symbol, e->name, e->weight, e->weight_isotope, e->bsrad, e->vdwrad);
	}
}

ElemData	*GetElemData (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return &elemdata[an];
}

char	*GetElemSymbol (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].symbol;
}

int	GetElemNumber (char *sym)
{
	char	symbol[3];
	register int	i;
	register ElemData	*e;

	if (!sym || !*sym) return 0;

	if (sym[0] == ' ' || isdigit(sym[0])) {
		symbol[0] = sym[1];
		symbol[1] = ' ';
	} else {
		symbol[0] = sym[0];
		symbol[1] = sym[1] ? sym[1] : ' ';
	}
	symbol[2] = '\0';

	if (strcasecmp(symbol, "D ") == 0) return  1;	/* deuterium */
	if (strcasecmp(symbol, "T ") == 0) return  1;	/* tritium */

	if (strcasecmp(symbol, "DU") == 0 || strcasecmp(symbol, "XX") == 0 ||
		strcasecmp(symbol, "X ") == 0 || strcasecmp(symbol, " X") == 0) return -1; /* dummy */

	for(i=0,e=elemdata;i<MAXELEM;i++,e++) if (strcasecmp(symbol, e->symbol) == 0) return i;
	
	/* unknown atomic symbol */

	return 0;
}

float	GetElemBSRadius (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].bsrad;
}

float	GetElemVdwRadius (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].vdwrad;
}

float	GetElemWeight (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].weight;
}

float	GetElemIsotopeWeight (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].weight_isotope;
}

int	GetElemCPKColor (int an)
{
	if (!IsValidElemNo(an)) an = 0;
	return elemdata[an].col;
}

void	SetElemData (int an, ElemData *edata)
{
	if (!IsValidElemNo(an)) return;
	BCOPY(&elemdata[an], edata, sizeof(ElemData));
}

void	SetElemBSRadius (int an, float r)
{
	if (!IsValidElemNo(an)) return;
	elemdata[an].bsrad = r;
}

void	SetElemVdwRadius (int an, float r)
{
	if (!IsValidElemNo(an)) return;
	elemdata[an].vdwrad = r;
}

void	SetElemWeight (int an, float w)
{
	if (!IsValidElemNo(an)) return;
	elemdata[an].weight = w;
}

void	SetElemIsotopeWeight (int an, float w)
{
	if (!IsValidElemNo(an)) return;
	elemdata[an].weight_isotope = w;
}

void	SetElemCPKColor (int an, int col)
{
	if (!IsValidElemNo(an)) return;
	elemdata[an].col = col;
}

int	TestMetalElem (int an)
{
	/* return TRUE if the atom (represented by the atomic number an) is metal */
	return ((an >=  3 && an <=  4) || (an >= 11 && an <= 13) ||
		(an >= 19 && an <= 32) || (an >= 37 && an <= 51) ||
		(an >= 55 && an <= 85) || (an >= 87));
}

/*#######################################################################*/

int	GetPdbAtomicNumber (int atm_refno, int res_refno)
/* atm_refno = reference number of the atom */
/* res_refno = reference number of the residue */
{
	char     sym[3], *name;

	name = GetAtomName(atm_refno);
	if (IsNADGroup(res_refno) || IsCOAGroup(res_refno)) {
		sym[0] = name[1];
		sym[1] = ' ';
	} else {
		sym[0] = name[0];
		sym[1] = name[1];
		if (sym[0] == 'H' && IsProtein(res_refno)) return 1;	/* handle HG, HD, etc */
	}
	sym[2] = '\0';
	return GetElemNumber(sym);
}

char	*GetPdbAtomName (int atm_refno, int res_refno)
/* atm_refno = reference number of the atom */
/* res_refno = reference number of the residue */
{
	static char	name[PDB_ASIZE];

	if (GetPdbAtomicNumber(atm_refno, res_refno) <= 0) {
		strcpy(name, "HE");	/* Helium */
		return name;
	} else {
		return GetAtomName(atm_refno);
	}
}

char	*MakePdbAtomName (char *str, int an)
{
	char	sym[3];
	static char	name[PDB_ASIZE];	/* PDB_ASIZE = PDB_ALEN+1 */

	name[0] = '\0';
	if (!str) return " X";
	strncpy(name, str+strspn(str," "), PDB_ASIZE);
	name[PDB_ASIZE-1] = '\0';
	if (!name[0]) return " X";

	if (isdigit(name[1]) || name[1] == ' ' || !name[1]) {
		sym[0] = ' ';
		sym[1] = name[0];
		sym[2] = '\0';

		name[3] = name[2];
		name[2] = name[1];
		name[1] = sym[1];
		name[0] = sym[0];
	} else {
		sym[0] = name[0];
		sym[1] = name[1];
		sym[2] = '\0';
		if (an > 0 && GetElemNumber(sym) != an) {
			sym[0] = ' ';
			sym[1] = name[0];
			sym[2] = '\0';

			name[3] = name[2];
			name[2] = name[1];
			name[1] = sym[1];
			name[0] = sym[0];
		}
	}

	/* hydrogen */
	if (an == 1 && strcmp(sym, " D") == 0 && strcmp(sym, " T") == 0) {
		return name;
	}

	/* dummy */
	if (strcmp(sym, " X") == 0 || strcmp(sym, "DU") == 0) {
		return name;
		/*
		sprintf(name, "HE");
		*/
	} else if (GetElemNumber(sym) != an) {
		strcpy(sym, GetElemSymbol(an));
		name[0] = (sym[1] == ' ') ? ' ' : sym[0];
		name[1] = (sym[1] == ' ') ? sym[0] : sym[1];
		name[2] = ' ';
		name[3] = ' ';
	}

	return name;
}

typedef struct {
	int	refno;	/* reference number */
	int	n;	/* position in string table */
	} StrTab;

static RefTab	atm_ref_tab[] = ATM_REF_TAB;
static RefTab	res_ref_tab[] = RES_REF_TAB;

static StrTab	*atm_tab = NULL;
static StrTab	*res_tab = NULL;

static int	cur_res_tab_size = 0;
static int	max_res_tab_size = 0;
static int	cur_atm_tab_size = 0;
static int	max_atm_tab_size = 0;

static int	max_res_refno = 0;
static int	max_atm_refno = 0;
static int	res_tab_sorted = 0;
static int	atm_tab_sorted = 0;

static int	str_tab_initialized = 0;

static void	SortStrTab (StrTab *tab, int left, int right)
{
	register int	i, j, n;
	register StrTab	tmp;

	i = left;
	j = right;
	n = tab[(i+j)/2].n;

	do {
		while (tab[i].n < n && i < right) i++;
		while (tab[j].n > n && j > left)  j--;
		if (i<=j) {
			tmp.refno = tab[i].refno;
			tmp.n = tab[i].n;

			tab[i].refno = tab[j].refno;
			tab[i].n = tab[j].n;

			tab[j].refno = tmp.refno;
			tab[j].n = tmp.n;

			i++; j--;
		}
	} while (i<=j);
	if (left < j)  SortStrTab(tab, left, j);
	if (i < right) SortStrTab(tab, i, right);
}

static void	InitStrTab (void)
{
	int	i;
	StrTab	*tab;

	if (!str_tab_initialized) str_tab_initialized = 1; else return;
	cur_res_tab_size = max_res_tab_size = sizeof(res_ref_tab)/sizeof(RefTab);
	cur_atm_tab_size = max_atm_tab_size = sizeof(atm_ref_tab)/sizeof(RefTab);

	res_tab = (StrTab *)malloc(cur_res_tab_size * sizeof(StrTab));
	atm_tab = (StrTab *)malloc(cur_atm_tab_size * sizeof(StrTab));
	if (!res_ref_tab || !atm_ref_tab) return;

	/* store names in the string table */
	for(i=0;i<cur_res_tab_size;i++) {
		res_tab[i].refno = res_ref_tab[i].refno;
		res_tab[i].n = LookUpStrTable(res_ref_tab[i].name);
	}
	for(i=0;i<cur_atm_tab_size;i++) {
		atm_tab[i].refno = atm_ref_tab[i].refno;
		atm_tab[i].n = LookUpStrTable(atm_ref_tab[i].name);
	}

	/* sort tables */
	SortStrTab(res_tab, 0, cur_res_tab_size-1);
	SortStrTab(atm_tab, 0, cur_atm_tab_size-1);

	res_tab_sorted = 1;
	atm_tab_sorted = 1;

	/* find the highest refno */
	for(i=0,tab=res_tab,max_res_refno=0;i<cur_res_tab_size;i++,tab++) {
		if (tab->refno > max_res_refno) max_res_refno = tab->refno;
	}
	for(i=0,tab=atm_tab,max_atm_refno=0;i<cur_atm_tab_size;i++,tab++) {
		if (tab->refno > max_atm_refno) max_atm_refno = tab->refno;
	}
}

#define STEPSIZE	64

static int	FindRefno (char *name, int *tab_sorted, StrTab **tab, int *max_refno, int *cur_tab_size, int *max_tab_size)
{
	register int	n, i, lo, hi, mid, sortflag;
	register StrTab	*t;

	if (!name || !*name) return -1;
	if (!str_tab_initialized) InitStrTab();

	n = LookUpStrTable(name);

	/* see if it's in tab */
	if (*tab_sorted) {	/* binary search */
		t = *tab;
		lo = 0;
		hi = *cur_tab_size-1;
		while (lo <= hi) {
			mid = (lo+hi) >> 1;
			if (n == t[mid].n) return t[mid].refno;
			else if (n < t[mid].n) hi = mid-1;
			else lo = mid+1;
		}
	} else {	/* sequential search */
		for(i=0,t=(*tab);i<*cur_tab_size;i++,t++) if (t->n == n) return t->refno;
	}

	/* The name is not found in the tab. So, increment max_refno */
	*max_refno = *max_refno + 1;

	/* store a new entry in tab */
	if (*cur_tab_size >= *max_tab_size) {
		*max_tab_size += STEPSIZE;
		if (!(*tab = (StrTab *)realloc(*tab, *max_tab_size * sizeof(StrTab)))) return -1;
		sortflag = 1;
	} sortflag = 0;

	t = *tab + *cur_tab_size;
	t->n = n;
	t->refno = *max_refno;
	*cur_tab_size = *cur_tab_size + 1;

	if (sortflag) {
		SortStrTab(*tab, 0, *cur_tab_size-1);
		*tab_sorted = 1;
	}

	return t->refno;
}

int	GetResRefno (char *resname)
{
	int	refno;

	refno = FindRefno(resname, &res_tab_sorted, &res_tab,
		&max_res_refno, &cur_res_tab_size, &max_res_tab_size);
	return (refno == -1 ? RES_UNK : refno);
}

int	GetAtomRefno (char *name)
{
	int	refno;

	refno = FindRefno(name, &atm_tab_sorted, &atm_tab,
		&max_atm_refno, &cur_atm_tab_size, &max_atm_tab_size);
	return (refno == -1 ? ATM_UNK : refno);
}

/*
 * name = atom name 
 * res_refno =  reference number of the residue
 */
int	GetPdbAtomRefno (char *name, int res_refno)
{
	static char	str[5];
	register int	i;

	if (isdigit(name[1])) {
		str[0] = str[2] = str[3] = ' ';
		str[1] = toupper(name[0]);
	} else for(i=0;i<4;i++) str[i] = toupper(name[i]);
	str[4] = '\0';

	if (IsProtein(res_refno) && str[0] == 'H') {
		str[0] = ' ';
		str[1] = 'H';
	} else if (IsNucleo(res_refno) && str[3] == '\'') str[3] = '*';
	return GetAtomRefno(str);
}

int	GetSimpleAtomRefno (char *name)
{
	char	str[5];

	str[2] = str[3] = ' ';
	if (name[1] && name[1] != ' ') {
		str[0] = toupper(name[0]);
		str[1] = toupper(name[1]);
	} else {
		str[0] = ' ';
		str[1] = name[0];
	}
	str[4] = '\0';
	return GetAtomRefno(str);
}

char	*GetResName (int refno)
{
	register int	i;
	register StrTab	*t;

	for(i=0,t=res_tab;i<cur_res_tab_size;i++,t++) if (t->refno == refno) return GetStrValue(t->n);
	return NULL;
}

char	*GetAtomName (int refno)
{
	register int	i;
	register StrTab	*t;

	for(i=0,t=atm_tab;i<cur_atm_tab_size;i++,t++) if (t->refno == refno) return GetStrValue(t->n);
	return NULL;
}

Residue	*FindResBySeqno (Residue *reslist, int seqno)
{
	Residue	*r;
	for(r=reslist;r;r=r->next) if (r->seqno == seqno) return r;
	return NULL;
}

Residue	*FindResByRefno (Residue *reslist, int refno)
{
	Residue	*r;

	for(r=reslist;r;r=r->next) if (r->refno == refno) return r;
	return NULL;
}

AtomPtr	FindAtomInRes (Residue *res, int refno)
{
	register AtomPtr	atom;
	if (!res) return NULL;
	ForEachAtom(res->atom, atom) if (atom->refno == refno) return atom;
	return NULL;
}

