#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

static int	readgaussroute (FILE *fp, GaussZmat *gauss)
{
	int	maxlen;
	char	str[256], *c;

	if (!fp || !gauss) return 0;

	/* % cards */
	gauss->percent[0] = '\0';
	gauss->route[0] = '\0';
	maxlen = sizeof(gauss->percent);
	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		if (str[0] == '$') continue;
		if (str[0] == '%') {
			if (strlen(str)+strlen(gauss->percent) < maxlen) strcat(gauss->percent, str);
		} else {
			strcpy(gauss->route, str);	/* 1st route card */
			break;
		}
	}

	/* remaining route cards */
	maxlen = sizeof(gauss->route);
	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		if (BLANKLINE(str)) break;	/* break at a blank line */
		if (strlen(str)+strlen(gauss->route) < maxlen) strcat(gauss->route, str);
	}

	/* title */
	maxlen = sizeof(gauss->title);
	gauss->title[0] = '\0';
	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		if (BLANKLINE(str)) break;	/* break at a blank line */
		if (strlen(str)+strlen(gauss->title) < maxlen) strcat(gauss->title, str);
	}

	/* trim */
	TrimTail(gauss->percent, " \n");
	TrimTail(gauss->route,   " \n");
	TrimTail(gauss->title,   " \n");

	/* charge and multiplicity */
	if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
	if ((c = strchr(str, '\n'))) *c = '\0';
	if ((c = strchr(str, ','))) *c = ' ';
	if (sscanf(str, "%d %d", &gauss->charge, &gauss->multiplicity) != 2) {
		setmolerrortext(str);
		return setmolerrorno(MERR_FORMAT);
	}

	return 1;
}

static int	readgausszmat (FILE *fp, GaussZmat *gauss, int *eof)
{
	int	nzmat, ncart, count;
	char	str[256], line[256], *c, *p;
	ZmatPtr	zmat;

	if (!fp || !gauss || !eof) return 0;
	*eof = 0;
	gauss->flags &= ~GAUSS_CARTESIAN;

	/* Z-matrix is terminated by a blank line OR
	 * when the first 8 char of the fist string of the line
	 * is "Variable" (case-insensitive)
	 */
	
	gauss->zmatlist = NULL;
	nzmat = 0;
	while (1) {
		if (!READLINE(fp, str)) {
			if (!(gauss->flags & GAUSS_CARTESIAN) && !gauss->zmatlist) setmolerrorno(MERR_EOF);
			break;
		}
		if (strstr(str, "$EndCoord:") || strstr(str, "EndGaussZmat:")) {
			*eof = 1;
			break;
		}
		if ((c = strrchr(str, '\n'))) *c = '\0';
		setmolerrortext(str);
		strcpy(line, str+strspn(str, " \t"));

		/* test if the line is blank, or the first string is "Variable" */
		if (BLANKLINE(line) || strncasecmp(line, "Variables:", 8) == 0) break;	/* End of Z-matrix */

		if (!(zmat = NewZmat())) {
			setmolerrorno(MERR_MEM);
			break;
		}
		nzmat++;
		zmat->atom = nzmat;

		/* break up the line into tokens */
		for(c=line,count=0,ncart=0;(p = strtok(c, " ,\t\n")) && ncart < 3;count++,c=(char *)NULL) {
			switch (count) {
			case 0:	/* center (alphabetic or alphabetic followed by number), or atomic number */
				strcpy(zmat->atomname[0], p);
				ParseZmatCenter(p, gauss->zmatlist, zmat);
				break;
			case 1:	/* may be positive integer, center name, floating-point number, 0, or -1 */
				if (strcmp(p, "0") == 0 || strcmp(p, "-1") == 0) {	/* cartesian flag */
					strcpy(zmat->atomname[1], p);
					sscanf(p, "%d", &zmat->zdef[0]);
					zmat->flags |= Z_CARTESIAN;
					gauss->flags |= GAUSS_CARTESIAN;
				} else if (IsFloat(p)) {
					zmat->flags |= Z_CARTESIAN;
					gauss->flags |= GAUSS_CARTESIAN;
					sscanf(p, "%lf", &zmat->cart[ncart]);
					ncart++;
				} else {	/* integer or symbol+digits */
					ParseZmatRef(p, gauss->zmatlist, zmat, 0);
				}
				break;
			case 2:	/* may be alphanumeric or floating-point number */
				if (zmat->flags & Z_CARTESIAN) {
					if (sscanf(p, "%lf", &zmat->cart[ncart]) != 1) {
						setmolerrorno(MERR_INVALID_FLOAT);
						break;
					} else ncart++;
				} else {
					if (IsAlphaNumeric(p, 0)) strcpy(zmat->varname[0], p);
					else if (IsFloat(p)) sscanf(p, "%lf", &zmat->zval[0]);
					else setmolerrorno(MERR_FORMAT);
				}
				break;
			case 3:	/* may be integer or center name or floating point number */
				if (zmat->flags & Z_CARTESIAN) {
					if (sscanf(p, "%lf", &zmat->cart[ncart]) != 1) {
						setmolerrorno(MERR_INVALID_FLOAT);
						break;
					} else ncart++;
				} else ParseZmatRef(p, gauss->zmatlist, zmat, 1);
				break;
			case 4:	/* must be alphanumberic or floating-point number */
				if (zmat->flags & Z_CARTESIAN) {
					if (sscanf(p, "%lf", &zmat->cart[ncart]) != 1) {
						setmolerrorno(MERR_INVALID_FLOAT);
						break;
					} else ncart++;
				} else {
					if (IsAlphaNumeric(p, 0)) strcpy(zmat->varname[1], p);
					else if (IsFloat(p)) sscanf(p, "%lf", &zmat->zval[1]);
					else setmolerrorno(MERR_FORMAT);
				}
				break;
			case 5:	/* may be integer or center name */
				ParseZmatRef(p, gauss->zmatlist, zmat, 2);
				break;
			case 6:	/* (+/-)alphanumeric or floating-point number */
				if (IsAlphaNumeric(p, 1)) strcpy(zmat->varname[2], p);
				else if (IsFloat(p)) sscanf(p, "%lf", &zmat->zval[2]);
				else setmolerrorno(MERR_FORMAT);
				break;
			case 7:
				if (strcmp(p, "1") == 0) zmat->flags |= Z_CHIRAL_PLUS;
				else if (strcmp(p, "-1") == 0) zmat->flags |= Z_CHIRAL_MINUS;
				break;
			}
			if (ncart == 3) break;
			if (getmolerrorno() != -1) break;
		}

		/* error checking */
		if (getmolerrorno() != -1) {
			FreeZmat(&zmat);
			break;
		} else EnterZmat(zmat, &gauss->zmatlist);

		if (!(zmat->flags & Z_CARTESIAN) && !(gauss->flags & GAUSS_CARTESIAN)) {
			double	len, theta, phi;
			/* just assign legal values if varname is not NULL,
			 * because their values have not been assigned yet
			 */
			if (zmat->varname[0]) len = 1.0;
			if (zmat->varname[1]) theta = 90.0;
			phi = 0.0;

			if (!CheckInternalCoord(nzmat, zmat->zdef[0], zmat->zdef[1], zmat->zdef[2],
				len, theta, phi)) {
				setmolerrorno(MERR_INTERNAL_COORD);
				break;
			}
		}
	}
	if (getmolerrorno() != -1) {
		FreeZmat(&gauss->zmatlist);
		return 0;
	}
	return 1;
}

typedef char	StrName[Z_NAME_LEN+1];
static int	readgaussvar (FILE *fp, GaussZmat *gauss, int *eof)
{
	int	i, nzmat, maxvar, var_section, nvar, zmat_nvar;
	int	*varflag;
	double	*varval;
	StrName	*varname;
	char	*c, str[256];
	ZmatPtr	zmat;

	if (!fp || !gauss || !eof) return 0;
	*eof = 0;
	for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next) nzmat++;
	if (nzmat == 0) return 0;
	zmat_nvar = CountIntZmatVar(gauss->zmatlist);

	maxvar = nzmat * 3;
	if (!(varflag = (int   *)calloc(maxvar, sizeof(int))))    return setmolerrorno(MERR_MEM);
	if (!(varval  = (double*)calloc(maxvar, sizeof(double)))) return setmolerrorno(MERR_MEM);
	if (!(varname = (StrName *)calloc(maxvar, sizeof(StrName))))  return setmolerrorno(MERR_MEM);

	var_section = 1;
	nvar = 0;
	while (1) {
		if (!READLINE(fp, str)) break;
		if ((c = strrchr(str, '\n'))) *c = '\0';
		if (strstr(str, "$EndCoord:") || strstr(str, "EndGaussZmat:")) {
			*eof = 1;
			break;
		}
		setmolerrortext(str);
		if (var_section) {
			if (BLANKLINE(str) || strncmp(str+strspn(str, " \t"), "Constant", 8) == 0) {
				var_section = 0;
				continue;
			}
		} else if (BLANKLINE(str)) break;

		if ((c = strrchr(str, '='))) *c = ' ';
		if (sscanf(str, "%s %lf", varname[nvar], &varval[nvar]) != 2) {
			setmolerrorno(MERR_INVALID_VAR_ASSIGN);
			break;
		}
		varflag[nvar] = var_section;
		nvar++;

		if (nvar >= maxvar || nvar >= zmat_nvar) break;
	}

	/* assign values */
	for(zmat=gauss->zmatlist,nzmat=1;zmat;zmat=zmat->next,nzmat++) {
		if (nzmat == 1) continue;
		if (zmat->varname[0][0]) {	/* length id */
			for(i=0;i<nvar;i++) {
				if (strcmp(varname[i], zmat->varname[0]) == 0) {
					zmat->zval[0] = varval[i];
					if (varflag[i]) zmat->flags |= Z_VAR_LENGTH;
					break;
				}
			}
			if (i == nvar) {
				setmolerrortext(zmat->varname[0]);
				setmolerrorno(MERR_VAR_NOT_FOUND);
				break;
			}
		}
		if (nzmat == 2) continue;
		if (zmat->varname[1][0]) {	/* angle ID */
			for(i=0;i<nvar;i++) {
				if (strcmp(varname[i], zmat->varname[1]) == 0) {
					zmat->zval[1] = varval[i];
					if (varflag[i]) zmat->flags |= Z_VAR_ANGLE;
					break;
				}
			}
			if (i == nvar) {
				setmolerrortext(zmat->varname[1]);
				setmolerrorno(MERR_VAR_NOT_FOUND);
				break;
			}
		}
		if (nzmat == 3) continue;
		if (zmat->varname[2][0]) {	/* dihedral identifier */
			c = zmat->varname[2];
			if (c[0] == '-' || c[0] == '+') c++;
			for(i=0;i<nvar;i++) {
				if (strcmp(varname[i], c) == 0) {
					zmat->zval[2] = varval[i];
					if (zmat->varname[2][0] == '-') zmat->zval[2] *= -1.0;
					if (varflag[i]) zmat->flags |= Z_VAR_DIHED;
					break;
				}
			}
			if (i == nvar) {
				setmolerrortext(zmat->varname[1]);
				setmolerrorno(MERR_VAR_NOT_FOUND);
				break;
			}
		}
	}

	free(varname);
	free(varval);
	free(varflag);

	/* read the remaining Z-matrix */
	if ((getmolerrorno() == -1) && !*eof) {	/* no error has been detected */
		int	nlines=0;

		gauss->data = NULL;
		while (READLINE(fp, str)) {
			nlines++;
			if (nlines > 1000) break;	/* something is wrong */
			if (strstr(str, "$EndCoord:") || strstr(str, "EndGaussZmat:")) break;
			AddtoString(&gauss->data, str);
		}
	}

	return (getmolerrorno() == -1);
}

void	FreeGaussZmat (GaussZmatPtr *gauss)
{
	if (!gauss || !*gauss) return;
	if ((*gauss)->zmatlist) FreeZmat(&(*gauss)->zmatlist);
	if ((*gauss)->data) free((*gauss)->data);
	*gauss = NULL;
}

GaussZmatPtr	CopyGaussZmat (GaussZmatPtr oldgauss)
{
	GaussZmatPtr	gauss;
	int	n;

	if (!oldgauss) return NULL;
	if (!(gauss = NewGaussZmat())) return NULL;
	BCOPY(oldgauss, gauss, sizeof(GaussZmat));
	gauss->zmatlist = CopyZmat(oldgauss->zmatlist);

	if (oldgauss->data && ((n = strlen(oldgauss->data)) > 0)) {
		gauss->data = (char *)malloc(n+1);
		strncpy(gauss->data, oldgauss->data, n);
		gauss->data[n] = '\0';
	}

	return gauss;
}

GaussZmatPtr	FLoadGaussZmat (FILE *fp)
{
	GaussZmatPtr	gauss;
	int	eof=0;

	clearmolerror();
	if (!(gauss = NewGaussZmat())) return NULL;

	/* read %, $, route, title, and charge and multiplicity */
	readgaussroute(fp, gauss);	/* ignore errors */

	if (!readgausszmat(fp, gauss, &eof) || !readgaussvar(fp, gauss, &eof)) {
		FreeGaussZmat(&gauss);
		return NULL;
	}

	/* convert internal coordinates to cartesian */
	if (!(gauss->flags & GAUSS_CARTESIAN)) ZmatToCartesian(gauss->zmatlist);

	return gauss;
}

MolPtr	CreateMolFromGauss (GaussZmatPtr gauss)
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

	if (!gauss || !gauss->zmatlist) return NULL;
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
	ForEachZmat(gauss->zmatlist,zmat) {
		serno++;
		patom->serial = serno;
		strcpy(patom->name, AN2ATOMNAME(zmat->an));
		patom->x = zmat->cart[0];
		patom->y = zmat->cart[1];
		patom->z = zmat->cart[2];
		if (zmat->resname[0]) {
			strcpy(pres->name, zmat->resname);
			pres->seq = zmat->resseq;
		}

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

#define DEF_GAUSS_ROUTE	"#N RHF/6-31G* OPT SCF=DIRECT POP=CHELPG IOP(6/42=3)"
char	gauss_route[512] = DEF_GAUSS_ROUTE;

char	*GetGaussRoute (void)
{
	return gauss_route;
}

void	SetGaussRoute (char *route)
{
	strcpy(gauss_route, route);
}

GaussZmatPtr	NewGaussZmat (void)
{
	GaussZmatPtr	gauss;

	if (!(gauss = (GaussZmatPtr)calloc(1, sizeof(GaussZmat)))) return NULL;
	strcpy(gauss->route, gauss_route);
	strcpy(gauss->title, "GAUSSIAN Z-MATRIX");
	gauss->charge = 0;
	gauss->multiplicity = 1;

	return gauss;
}

GaussZmatPtr	CreateGaussZmatFromMol (MolPtr mol)
{
	GaussZmatPtr	gauss;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	ZmatPtr	zmat;

	if (!mol || !(gauss = NewGaussZmat())) return NULL;
	gauss->zmatlist = NULL;
	gauss->flags |= GAUSS_CARTESIAN;

	ForEachChainResAtom(mol->chain,c,r,a) {
		if (!(zmat = EnterNewZmat(&gauss->zmatlist))) {
			FreeGaussZmat(&gauss);
			return NULL;
		}
		zmat->an = (int)a->an;
		if (zmat->an == 0) zmat->an = 1;
		zmat->cart[0] = a->x;
		zmat->cart[1] = a->y;
		zmat->cart[2] = a->z;
		zmat->atom = a->serno;

		strncpy(zmat->label, GetAtomName(a->refno), Z_LABEL_LEN);
		zmat->label[Z_LABEL_LEN] = '\0';

		if (!zmat->atomname[0][0]) {
			if (zmat->an < 0 || zmat->an >= 98) strcpy(zmat->atomname[0], "X");
			else strcpy(zmat->atomname[0], GetElemSymbol(zmat->an));
		}

		zmat->resseq = a->residue->seqno;
		strcpy(zmat->resname, GetResName(a->residue->refno));

		if (zmat->an >= 0) zmat->flags |= (Z_VAR_LENGTH | Z_VAR_ANGLE | Z_VAR_DIHED);
		else zmat->flags &= ~(Z_VAR_LENGTH | Z_VAR_ANGLE | Z_VAR_DIHED);
	}

	return gauss;
}

GaussZmatPtr	CreateAutoGaussZmatFromMol (MolPtr mol)
{
	GaussZmatPtr	g;
	ZmatPtr	zmatlist;

	if (!(zmatlist = CreateAutoZmatFromMol(mol, 0))) return NULL;
	g = NewGaussZmat();
	g->zmatlist = zmatlist;
	return g;
}

/******************** Z-matrix printing **************************/

static void	GPrintGaussVar (GaussZmatPtr gauss, void *dest, void (*print_func)(void *, char *))
{
	char	str[256];
	ZmatPtr	zmat;
	int	nconst = 0, nzmat;

	if (!gauss) return;

	sprintf(str, "        Variables:\n");
	(*print_func)(dest, str);
	for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
		if (nzmat < 1) continue;
		if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, zmat->flags & Z_VAR_LENGTH)) {
			continue;
		}
		if (zmat->varname[0][0]) {
			if (zmat->flags & Z_VAR_LENGTH) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[0], zmat->zval[0]);
				(*print_func)(dest, str);
			} else nconst++;
		}
	}
	for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
		if (nzmat < 2) continue;
		if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, zmat->flags & Z_VAR_ANGLE)) {
			continue;
		}
		if (zmat->varname[1][0]) {
			if (zmat->flags & Z_VAR_ANGLE) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[1], zmat->zval[1]);
				(*print_func)(dest, str);
			} else nconst++;
		}
	}
	for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
		if (nzmat < 3) continue;
		if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, zmat->flags & Z_VAR_DIHED)) {
			continue;
		}
		if (zmat->varname[2][0]) {
			if (zmat->flags & Z_VAR_DIHED) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[2], zmat->zval[2]);
				(*print_func)(dest, str);
			} else nconst++;
		}
	}

	if (nconst) {
		sprintf(str, "        Constants:\n");
		(*print_func)(dest, str);
		for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
			if (nzmat < 1) continue;
			if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, Z_VAR_LENGTH)) {
				continue;
			}
			if (zmat->varname[0][0] && !(zmat->flags & Z_VAR_LENGTH)) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[0], zmat->zval[0]);
				(*print_func)(dest, str);
			}
		}
		for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
			if (nzmat < 2) continue;
			if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, Z_VAR_ANGLE)) {
				continue;
			}
			if (zmat->varname[1][0] && !(zmat->flags & Z_VAR_ANGLE)) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[1], zmat->zval[1]);
				(*print_func)(dest, str);
			}
		}
		for(zmat=gauss->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
			if (nzmat < 3) continue;
			if (IsRedundantIntZmatVar(gauss->zmatlist, nzmat, zmat, Z_VAR_DIHED)) {
				continue;
			}
			if (zmat->varname[2][0] && !(zmat->flags & Z_VAR_DIHED)) {
				sprintf(str, " %-8s % 13.6f\n", zmat->varname[2], zmat->zval[2]);
				(*print_func)(dest, str);
			}
		}
	}

	if (gauss->data) {
		(*print_func)(dest, gauss->data);
	}
}

static void	GPrintGaussRoute (GaussZmatPtr gauss, void *dest, void (*print_func)(void *, char *))
{
	char	str[256];

	if (gauss->percent[0]) {
		sprintf(str, "%s\n", gauss->percent);
		(*print_func)(dest, str);
	}

	if (!gauss->route[0] || (gauss->flags & GAUSS_MOPAC_ZMAT)) {
		(*print_func)(dest, GetGaussRoute());
	} else (*print_func)(dest, gauss->route);
	(*print_func)(dest, "\n\n");

	if (!gauss->title[0]) {
		sprintf(str, "GAUSSIAN Z-MATRIX (Cartesian Coordinates)\n\n");
		(*print_func)(dest, str);
	} else {
		(*print_func)(dest, gauss->title);
		(*print_func)(dest, "\n\n");
	}
	sprintf(str, "%d %d\n", gauss->charge, gauss->multiplicity);
	(*print_func)(dest, str);
}

static void	GPrintGaussCartZmat (GaussZmatPtr gauss, void *dest, void (*print_func)(void *, char *))
{
	ZmatPtr	zmat;
	char	str[256], *sym;

	GPrintGaussRoute(gauss, dest, print_func);
	for(zmat=gauss->zmatlist;zmat;zmat=zmat->next) {
		if (strncasecmp(zmat->atomname[0], "Du", 2) == 0 ||
		    strncasecmp(zmat->atomname[0], "Gh", 2) == 0)
			sym = "X";
		else sym = (zmat->an == 0) ? "H" : &zmat->atomname[0][0];

		sprintf(str, "%-3s % 12.6f % 12.6f % 12.6f\n",
			sym, zmat->cart[0], zmat->cart[1], zmat->cart[2]);
		(*print_func)(dest, str);
	}
	/* needs a blank line at the end */
	(*print_func)(dest, "\n");
}

static void	GPrintGaussZmat (GaussZmatPtr gauss, void *dest, void (*print_func)(void *, char *))
{
	ZmatPtr	zmat;
	int	i, nzmat, done;
	char	str[256], sym[Z_NAME_LEN], *p;

	if (!gauss) return;
	if (gauss->flags & GAUSS_CARTESIAN) {
		GPrintGaussCartZmat(gauss, dest, print_func);
		return;
	}
	GPrintGaussRoute(gauss, dest, print_func);

	nzmat = 0;
	for(zmat=gauss->zmatlist;zmat;zmat=zmat->next) {
		nzmat++;
		for(i=0,done=0;i<8;i++) {
			str[0] = '\0';
			switch (i) {
			case 0:
				p = &zmat->atomname[0][0];
				if (strncasecmp(p, "Du", 2) == 0 ||
				    strncasecmp(p, "Gh", 2) == 0)
					sprintf(sym, "X%s", p+2);
				else if (zmat->an == 0 && *p == 'X') sprintf(sym, "H%s", p+1);
				else strcpy(sym, p);

				sprintf(str, " %-3s", sym);
				break;
			case 1:
				if (!zmat->atomname[1][0]) {	/* integer representation of center */
					if (nzmat == 1) {done = 1; break;}
					sprintf(str, " %-4d", zmat->zdef[0]);
				} else {
					if (strcmp(zmat->atomname[0], "0") == 0 || strcmp(zmat->atomname[0], "-1") == 0) {
						/* cartesian coord */
						sprintf(str, " %-4d % 9.6f % 9.6f % 9.6f",
							zmat->zdef[0], zmat->cart[0], zmat->cart[1], zmat->cart[2]);
						done = 1;
					} else sprintf(str, " %-4s", zmat->atomname[1]);
				}
				break;
			case 2:
				if (!zmat->varname[0][0]) sprintf(str, "   % 9.6f", zmat->zval[0]);
				else sprintf(str, "   %-9s", zmat->varname[0]);
				if (nzmat == 2) done = 1;
				break;
			case 3:
				if (!zmat->atomname[2][0]) sprintf(str, " %-4d", zmat->zdef[1]);
				else sprintf(str, " %-4s", zmat->atomname[2]);
				break;
			case 4:
				if (!zmat->varname[1][0]) sprintf(str, "   %10.6f", zmat->zval[1]);
				else sprintf(str, "   %-10s", zmat->varname[1]);
				if (nzmat == 3) done = 1;
				break;
			case 5:
				if (!zmat->atomname[3][0]) sprintf(str, " %-4d", zmat->zdef[2]);
				else sprintf(str, " %-4s", zmat->atomname[3]);
				break;
			case 6:
				if (!zmat->varname[2][0]) sprintf(str, "   % 11.6f", zmat->zval[2]);
				else sprintf(str, "   %-11s", zmat->varname[2]);
				break;
			case 7:
				if (zmat->flags & Z_CHIRAL_PLUS) sprintf(str, "  1");
				else if (zmat->flags & Z_CHIRAL_MINUS) sprintf(str, " -1");
				break;
			}
			if (str[0]) (*print_func)(dest, str);
			if (done) break;
		}
		(*print_func)(dest, "\n");
	}

	GPrintGaussVar(gauss, dest, print_func);
}

void	SPrintGaussZmat (char **str_sum, GaussZmatPtr gauss)
{
	if (!str_sum || !gauss) return;
	ContinuePrintString(*str_sum);
	GPrintGaussZmat(gauss, NULL, SPrintString);
	*str_sum = EndPrintString();
}

void	FPrintGaussZmat (FILE *fp, GaussZmatPtr gauss)
{
	if (!fp || !gauss) return;
	GPrintGaussZmat(gauss, (void *)fp, FPrintString);
}

void	FPrintGaussCartZmat (FILE *fp, GaussZmatPtr gauss)
{
	if (!fp || !gauss || !gauss->zmatlist) return;
	GPrintGaussCartZmat(gauss, (void *)fp, FPrintString);
}

void	SPrintGaussCartZmat (char **str_sum, GaussZmatPtr gauss)
{
	if (!str_sum || !gauss) return;
	ContinuePrintString(*str_sum);
	GPrintGaussCartZmat(gauss, NULL, SPrintString);
	*str_sum = EndPrintString();
}

