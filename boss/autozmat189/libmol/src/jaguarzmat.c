#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <ctype.h>

#include "macros.h"
#include "molecule.h"

#define VAR_GRADIENT	0x1
#define VAR_CONSTRAINED	0x2

/* var must be alphanumeric, and val must be a number.
   val#  --- constrained var
   val*  --- compute the gradient of the energy both at the
	     orginal geom and at a geom for which *-marked
	     coord has been changed slightly.
(E.g.)
 var1=val1,var2=val2 var3=val3# var4=val4*#
*/

static void	parse_var_flags (char *str, char *name, unsigned long *flags)
{
	char	*p, *q, buf[64];

	strcpy(buf, str);
	if ((p = strchr(buf, '*')) || (q = strchr(buf, '#'))) {
		if (p) {
			(*flags) |= VAR_GRADIENT;
			if ((q = strchr(p+1, '#'))) (*flags) |= VAR_CONSTRAINED;
			*p = '\0';
		} else {
			(*flags) |= VAR_CONSTRAINED;
			if ((p = strchr(q+1, '*'))) (*flags) |= VAR_GRADIENT;
			*q = '\0';
		}
	}
	strcpy(name, buf);
}

static int	parse_var_name (char *str, unsigned long sign_flags, char *var_name, unsigned long *var_flags)
{
	parse_var_flags(str, var_name, var_flags);
	return IsAlphaNumeric(var_name, sign_flags);
}

static int	parse_var (char *line, ListPtr *varlist)
{
	int	n;
	char	*s, *p, name[64];
	ListPtr	list;
	unsigned long	flags;

/*
fprintf(stderr, "VARLINE=%s\n", line);
*/

	for(s=line,n=0;(p=strtok(s, " =,\t\n"));n++,s=(char *)NULL) {
		if (n%2==0) {
			list = EnterNewList(varlist);
			if (IsAlphaNumeric(p, 0)) strcpy(list->L_VARNAME, p);
			else return setmolerrorno(MERR_INVALID_VAR_ASSIGN);
		} else {
			flags = 0;
			parse_var_flags(p, name, &flags);
			list->L_VARFLAGS |= flags;

			if (IsNumber(name)) sscanf(name, "%f", &list->L_VARVAL);
			else return setmolerrorno(MERR_INVALID_VAR_ASSIGN);
		}
	}

	return 1;
}

static void	print_var (ListPtr varlist)
{
	ListPtr	list;

	ForEachList(varlist, list) {
		fprintf(stderr, "VarName=%s: Val=%f Flags=%d\n",
			list->L_VARNAME, list->L_VARVAL, list->L_VARFLAGS);
	}
}

static int	assign_var (ZmatPtr zmatlist, ListPtr varlist, unsigned long flags)
{
	ZmatPtr	zmat;
	ListPtr	list;
	int	i;
	char	*p;
	double	sign;

	ForEachZmat(zmatlist,zmat) {
		for(i=0;i<3;i++) {
			p = zmat->varname[i];
			if (!*p) continue;
			sign = (*p == '-') ? -1.0 : 1.0;
			if (*p == '-' || *p == '+') p++;

			ForEachList(varlist,list) if (strcmp(list->L_VARNAME, p) == 0) break;
			if (!list) {
				setmolerrortext(zmat->varname[i]);
				setmolerrorno(MERR_VAR_NOT_FOUND);
				break;
			}

			if (flags & JAGUAR_CARTESIAN) zmat->cart[i] = sign * list->L_VARVAL;
			else zmat->zval[i] = sign * list->L_VARVAL;

			if (list->L_VARFLAGS & VAR_CONSTRAINED) {
				switch (i) {
				case 0: zmat->flags &= ~Z_VAR_LENGTH; break;
				case 1: zmat->flags &= ~Z_VAR_ANGLE;  break;
				case 2: zmat->flags &= ~Z_VAR_DIHED;  break;
				}
			}
			if (list->L_VARFLAGS & VAR_GRADIENT) {
				switch (i) {
				case 0: zmat->flags |= Z_GRAD_LENGTH; break;
				case 1: zmat->flags |= Z_GRAD_ANGLE;  break;
				case 2: zmat->flags |= Z_GRAD_DIHED;  break;
				}
			}
		}
		if (getmolerrorno() != -1) break;
	}
	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	parse_jaguar_zmat_center (char *str, ZmatPtr zmatlist, ZmatPtr zmat)
{
	int	i;
	char	*p, symbol[8];

	if (IsInt(str)) return setmolerrorno(MERR_ILLEGAL_CENTERNAME);
	strcpy(zmat->label, str);

	if ((p = strrchr(str, '@')) && p == str+strlen(str)-1) {
		*p = '\0';
		zmat->flags |= Z_COUNTER_POISE;
	}

	strcpy(zmat->atomname[0], str);
	strcpy(symbol, str);
	symbol[7] = '\0';
	for(i=0;i<strlen(symbol) && i<8;i++) {
		if (isdigit(symbol[i])) {
			symbol[i] = '\0';
			break;
		}
	}

	if ((zmat->an = GetElemNumber(symbol)) == 0) {
		if (strcasecmp(symbol, "X") == 0) zmat->an = -1;
		else return setmolerrorno(MERR_UNKNOWN_ATOMSYMBOL);
	} else {
		if (strcmp(symbol, str) != 0 && FindZmatByName(zmatlist, str)) {
			return setmolerrorno(MERR_MULTI_DEF_CENTER);
		}
	}
	return 1;
}

static int	parse_cart_zmat (char *line, ZmatPtr zmatlist, ZmatPtr zmat)
{
	int	count, flags;
	char	*c, *p, var_name[64];
	unsigned long	var_flags;

/*
fprintf(stderr, "CARTZMAT=%s:\n", line);
*/

	/* break up the line into tokens */
	for(c=line,count=0;(p = strtok(c, " ,\t\n"));count++,c=(char *)NULL) {
/*
fprintf(stderr, "TOKEN=%s:\n", p);
*/
		flags = (Z_VAR_LENGTH|Z_VAR_ANGLE|Z_VAR_DIHED);
		switch (count) {
		case 0:	/* atom label (alphanumeric, optionally followed by @), max 4 chars  */
			if (*p == '#') continue;
			if (!parse_jaguar_zmat_center(p, zmatlist, zmat)) break;
			break;

		case 1:	/* alphanumeric or a number, optionally followed by '*' or '#'. */
		case 2:
		case 3:
			var_flags = 0;
			parse_var_flags(p, var_name, &var_flags);

/*
fprintf(stderr, "varname=%s: flags=%d\n", var_name, var_flags);
*/

			if (IsAlphaNumeric(var_name, 1)) strcpy(zmat->varname[count-1], var_name);
			else if (IsNumber(var_name)) sscanf(var_name, "%lf", &zmat->cart[count-1]);
			else {
/*
fprintf(stderr, "ERROR in VARNAME: %s\n", p);
*/

				setmolerrorno(MERR_FORMAT);
				break;
			}

			if (var_flags & VAR_GRADIENT) {
				switch (count) {
				case 1: flags |= Z_GRAD_LENGTH; break;
				case 2: flags |= Z_GRAD_ANGLE;  break;
				case 3: flags |= Z_GRAD_DIHED;  break;
				}
			}
			if (var_flags & VAR_CONSTRAINED) {
				switch (count) {
				case 1: flags &= ~Z_VAR_LENGTH; break;
				case 2: flags &= ~Z_VAR_ANGLE;  break;
				case 3: flags &= ~Z_VAR_DIHED;  break;
				}
			}

			zmat->flags |= flags;
			break;
		}
		if (getmolerrorno() != -1) break;
	}
	if (count < 4 && getmolerrorno() == -1) setmolerrorno(MERR_FORMAT);
	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	parse_zmat (char *line, ZmatPtr zmatlist, ZmatPtr zmat, int nzmat)
{
	int	count;
	char	*c, *p;
	double	len, theta, phi;

	zmat->flags |= (Z_VAR_LENGTH|Z_VAR_ANGLE|Z_VAR_DIHED);

	/* break up the line into tokens */
	for(c=line,count=0;(p = strtok(c, " ,\t\n"));count++,c=(char *)NULL) {
		switch (count) {
		case 0:	/* atom label (alphanumeric, optionally followed by @), max 4 chars  */
			if (*p == '#') continue;
			parse_jaguar_zmat_center(p, zmatlist, zmat);
			break;

		case 1:	/* may be positive integer or center name */
		case 3:
		case 5:
			ParseZmatRef(p, zmatlist, zmat, (int)(count/2));
			break;

		case 2:	/* alphanumeric or a number */
		case 4:	/* alphanumeric or a number */
		case 6:	/* (+/-)alphanumeric or a number */
			if (IsAlphaNumeric(p, count==6)) strcpy(zmat->varname[(int)(count/3)], p);
			else if (IsNumber(p)) sscanf(p, "%lf", &zmat->zval[(int)(count/3)]);
			else setmolerrorno(MERR_FORMAT);
			break;
		}
		if (getmolerrorno() != -1) break;
	}

	/* check internal coord */
	if (getmolerrorno() == -1) {
		/* just assign legal values if varname is not NULL,
		 * because their values have not been assigned yet
		 */
		if (zmat->varname[0]) len = 1.0;
		if (zmat->varname[1]) theta = 90.0;
		phi = 0.0;
		if (!CheckInternalCoord(nzmat, zmat->zdef[0], zmat->zdef[1], zmat->zdef[2], len, theta, phi)) {
			setmolerrorno(MERR_INTERNAL_COORD);
		}
	}

	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	read_jaguar_zmat (FILE *fp, JaguarZmatPtr j, ListPtr *var_def_list)
{
	ZmatPtr	zmat;
	int	count, nzmat=0;
	char	str[256], line[256], *c, *p;

	if (!fp || !j || !var_def_list) return 0;
	clearmolerror();

	j->flags &= ~JAGUAR_CARTESIAN;
	j->zmatlist = NULL;

	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		TrimHead(str, " \t");
		TrimTail(str, " \t\n");
/*
fprintf(stderr, "LINE=%s:\n", str);
*/

		if (str[0] == '#' || !str[0]) continue;
		if (str[0] == '&' || str[0] == '$') break;
		setmolerrortext(str);

		if (strchr(str, '=')) {
			parse_var(str, var_def_list);
			continue;
		}

		if (!j->zmatlist) {
			strcpy(line, str);
			for(c=line,count=0;(p = strtok(c, " ,\t\n"));c=(char *)NULL) count++;
			if (count >= 4) j->flags |= JAGUAR_CARTESIAN;
/*
fprintf(stderr, "zmat type (cartesian?) = %d\n", j->flags & JAGUAR_CARTESIAN);
*/
		}

		nzmat++;
		if (!(zmat = NewZmat())) {
			setmolerrorno(MERR_MEM);
			break;
		}
		zmat->atom = nzmat;

		if (j->flags & JAGUAR_CARTESIAN) parse_cart_zmat(str, j->zmatlist, zmat);
		else parse_zmat(str, j->zmatlist, zmat, nzmat);

		if (getmolerrorno() != -1) {
			FreeZmat(&zmat);
			break;
		} else EnterZmat(zmat, &j->zmatlist);
	}

	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	read_jaguar_zvar (FILE *fp, ListPtr *var_def_list)
{
	char	str[256];

	if (!fp || !var_def_list) return 0;
	clearmolerror();

	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		TrimHead(str, " \t");
		TrimTail(str, " \t\n");

		if (str[0] == '#' || !str[0]) continue;
		if (str[0] == '&' || str[0] == '$') break;
		setmolerrortext(str);

		if (strchr(str, '=')) {
			parse_var(str, var_def_list);
			continue;
		}
		if (getmolerrorno() != -1) break;
	}

	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	read_jaguar_section (FILE *fp, ListPtr *top_list)
{
	char	str[256];
	ListPtr	list;

	if (!fp || !top_list) return 0;
	clearmolerror();
	*top_list = NULL;

	while (1) {
		if (!READLINE(fp, str)) return setmolerrorno(MERR_EOF);
		TrimHead(str, " \t");
		TrimTail(str, " \t\n");

		if (str[0] == '#' || !str[0]) continue;
		if (str[0] == '&' || str[0] == '$') break;

		if (!(list = EnterNewList(top_list))) return setmolerrorno(MERR_MEM);
		strcpy(list->L_JAGUAR_LINE, str);
	}

	return (getmolerrorno() == -1) ? 1 : 0;
}

static int	read_jaguar_input (FILE *fp, JaguarZmatPtr j)
{
	int	header=1;
	char	str[256], keywd[16];
	ListPtr	list, var_def_list=NULL;

	clearmolerror();
	while (1) {
		if (!READLINE(fp, str)) {
			if (!j->zmatlist) return setmolerrorno(MERR_EOF);
			else break;
		}
		if (strstr(str, "$EndCoord:")) break;
		TrimHead(str, " \t");
		TrimTail(str, " \t\n");
		if (str[0] == '#' || !str[0]) continue;
		if (str[0] != '&' && str[0] != '$') continue;

		/* parse main sections */
		keywd[0] = '\0';
		sscanf(str+1, "%s", keywd);
		if (strcmp(keywd, "zmat") == 0) {
/*
fprintf(stderr, "Reading jaguar zmat...\n");
*/
			header = 0;
			read_jaguar_zmat(fp, j, &var_def_list);
		}
		else if (strcmp(keywd, "zvar") == 0) read_jaguar_zvar(fp, &var_def_list);

		else if (strcmp(keywd, "gen")  == 0)    read_jaguar_section(fp, &j->gen);
		else if (strcmp(keywd, "gvb")  == 0)    read_jaguar_section(fp, &j->gvb);
		else if (strcmp(keywd, "lmp2") == 0)    read_jaguar_section(fp, &j->lmp2);
		else if (strcmp(keywd, "atomic") == 0)  read_jaguar_section(fp, &j->atomic);
		else if (strcmp(keywd, "vdw") == 0)     read_jaguar_section(fp, &j->vdw);
		else if (strcmp(keywd, "hess") == 0)    read_jaguar_section(fp, &j->hess);
		else if (strcmp(keywd, "guess") == 0)   read_jaguar_section(fp, &j->guess);
		else if (strcmp(keywd, "pointch") == 0) read_jaguar_section(fp, &j->pointch);
		else if (strcmp(keywd, "efields") == 0) read_jaguar_section(fp, &j->efields);
		else if (strcmp(keywd, "ham") == 0)     read_jaguar_section(fp, &j->ham);
		else if (strcmp(keywd, "orbman") == 0)  read_jaguar_section(fp, &j->orbman);
		else if (strcmp(keywd, "echo") == 0)    read_jaguar_section(fp, &j->echo);
		else if (strcmp(keywd, "path") == 0)    read_jaguar_section(fp, &j->path);
		else if (header) {	/* parse header */
			list = EnterNewList(&j->header);
			strcpy(list->L_JAGUAR_LINE, str);
		}
		if (getmolerrorno() != -1) break;
	}

	if (getmolerrorno() == -1) {
		assign_var(j->zmatlist, var_def_list, j->flags);
		if (var_def_list) FreeList(&var_def_list);
		/* convert internal coordinates to cartesian */
		if (!(j->flags & JAGUAR_CARTESIAN)) ZmatToCartesian(j->zmatlist);
	}

	return (getmolerrorno() == -1) ? 1 : 0;
}

void	FreeJaguarZmat (JaguarZmatPtr *jaguar)
{
	JaguarZmatPtr	j;

	if (!jaguar || !(j = *jaguar)) return;

	if (j->header)  FreeList(&j->header);
	if (j->zmatlist)FreeZmat(&j->zmatlist);
	if (j->gen)     FreeList(&j->gen);
	if (j->gvb)     FreeList(&j->gvb);
	if (j->lmp2)    FreeList(&j->lmp2);
	if (j->atomic)  FreeList(&j->atomic);
	if (j->vdw)     FreeList(&j->vdw);
	if (j->hess)    FreeList(&j->hess);
	if (j->guess)   FreeList(&j->guess);
	if (j->pointch) FreeList(&j->pointch);
	if (j->efields) FreeList(&j->efields);
	if (j->ham)     FreeList(&j->ham);
	if (j->orbman)  FreeList(&j->orbman);
	if (j->echo)    FreeList(&j->echo);
	if (j->path)    FreeList(&j->path);

	free((void *)j);
	*jaguar = NULL;
}

JaguarZmatPtr	FLoadJaguarZmat (FILE *fp)
{
	JaguarZmatPtr	j;

	clearmolerror();
	if (!(j = (JaguarZmat *)calloc(1, sizeof(JaguarZmat)))) return NULL;
	if (!read_jaguar_input (fp, j)) {
		FreeJaguarZmat(&j);
		return NULL;
	}
	return j;
}

/*********** creating Z-matrix ****************/
GaussZmatPtr	CreateGaussZmatFromJaguar (JaguarZmatPtr j)
{
	GaussZmatPtr	g;

	if (!j || !(g = NewGaussZmat())) return NULL;
	if (j->flags & JAGUAR_CARTESIAN) g->flags |= GAUSS_CARTESIAN;
	g->zmatlist = CopyZmat(j->zmatlist);

	return g;
}

JaguarZmatPtr	CreateJaguarZmatFromGauss (GaussZmatPtr g)
{
	JaguarZmatPtr	j;

	if (!g || !(j = (JaguarZmat *)calloc(1, sizeof(JaguarZmat)))) return NULL;
	if (g->flags & GAUSS_CARTESIAN) j->flags |= JAGUAR_CARTESIAN;
	j->zmatlist = CopyZmat(g->zmatlist);

	return j;
}

JaguarZmatPtr	CreateJaguarZmatFromMol (MolPtr mol)
{
	JaguarZmatPtr	j;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr		a;
	ZmatPtr		zmat;
	ListPtr		list;

	if (!mol || !(j = (JaguarZmat *)calloc(1, sizeof(JaguarZmat)))) return NULL;
	j->zmatlist = NULL;
	j->flags |= JAGUAR_CARTESIAN;

	ForEachChainResAtom(mol->chain,c,r,a) {
		if (!(zmat = EnterNewZmat(&j->zmatlist))) {
			FreeJaguarZmat(&j);
			return NULL;
		}
		zmat->an = (int)a->an;
		if (zmat->an == 0) zmat->an = 1;
		zmat->cart[0] = a->x;
		zmat->cart[1] = a->y;
		zmat->cart[2] = a->z;
		zmat->atom = a->serno;
		if (zmat->an < 0 || zmat->an >= 98) strcpy(zmat->label, "X");
		else strcpy(zmat->label, GetElemSymbol(zmat->an));
		strcpy(zmat->atomname[0], zmat->label);
		if (zmat->an >= 0) zmat->flags |= (Z_VAR_LENGTH | Z_VAR_ANGLE | Z_VAR_DIHED);
		else zmat->flags &= ~(Z_VAR_LENGTH | Z_VAR_ANGLE | Z_VAR_DIHED);       
	}

	/* gen section */
	list = EnterNewList(&j->gen);
	strcpy(list->L_JAGUAR_LINE, "igeopt=1");
	list = EnterNewList(&j->gen);
	strcpy(list->L_JAGUAR_LINE, "basis=6-31G*");
	list = EnterNewList(&j->gen);
	strcpy(list->L_JAGUAR_LINE, "ip192=2");

	return j;
}

JaguarZmatPtr	CopyJaguarZmat (JaguarZmatPtr oldjaguar)
{
	JaguarZmatPtr	jaguar, j;

	if (!oldjaguar) return NULL;
	if (!(jaguar = j = (JaguarZmat *)calloc(1, sizeof(JaguarZmat)))) return NULL;
	BCOPY(oldjaguar, jaguar, sizeof(JaguarZmat));

	jaguar->zmatlist = CopyZmat(oldjaguar->zmatlist);

	if (oldjaguar->header)  j->header  = CopyList(oldjaguar->header);
	if (oldjaguar->gen)     j->gen     = CopyList(oldjaguar->gen);
	if (oldjaguar->gvb)     j->gvb     = CopyList(oldjaguar->gvb);
	if (oldjaguar->lmp2)    j->lmp2    = CopyList(oldjaguar->lmp2);
	if (oldjaguar->atomic)  j->atomic  = CopyList(oldjaguar->atomic);
	if (oldjaguar->vdw)     j->vdw     = CopyList(oldjaguar->vdw);
	if (oldjaguar->hess)    j->hess    = CopyList(oldjaguar->hess);
	if (oldjaguar->guess)   j->guess   = CopyList(oldjaguar->guess);
	if (oldjaguar->pointch) j->pointch = CopyList(oldjaguar->pointch);
	if (oldjaguar->efields) j->efields = CopyList(oldjaguar->efields);
	if (oldjaguar->ham)     j->ham     = CopyList(oldjaguar->ham);
	if (oldjaguar->orbman)  j->orbman  = CopyList(oldjaguar->orbman);
	if (oldjaguar->echo)    j->echo    = CopyList(oldjaguar->echo);
	if (oldjaguar->path)    j->path    = CopyList(oldjaguar->path);

	return jaguar;
}

/*********** printing Z-matrix ****************/
static void	GPrintJaguarCartZmat (JaguarZmatPtr j, void *dest, void (*print_func)(void *, char *))
{
	ZmatPtr	zmat;
	char	str[256], *sym;

	if (!j) return;
	ForEachZmat(j->zmatlist, zmat) {
		if (strncasecmp(zmat->atomname[0], "Du", 2) == 0 ||
		    strncasecmp(zmat->atomname[0], "Gh", 2) == 0)
			sym = "X";
		else sym = (zmat->an == 0) ? "H" : &zmat->atomname[0][0];

		sprintf(str, "%-3s % 12.6f % 12.6f % 12.6f\n",
			sym, zmat->cart[0], zmat->cart[1], zmat->cart[2]);
		(*print_func)(dest, str);
	}
}

static void	GPrintJaguarIntZmat (JaguarZmatPtr j, void *dest, void (*print_func)(void *, char *))
{
	ZmatPtr	zmat;
	char	str[256], sym[32], *p;
	int	i, nzmat, done;

	if (!j) return;
	for(zmat=j->zmatlist,nzmat=1;zmat;zmat=zmat->next,nzmat++) {
		done = 0;
		for(i=0;i<7;i++) {
			str[0] = '\0';
			switch (i) {
			case 0:
				p = &zmat->atomname[0][0];
				if (strncasecmp(p, "Du", 2) == 0 ||
				    strncasecmp(p, "Gh", 2) == 0) sprintf(sym, "X%s", p+2);
				else if (zmat->an == 0 && *p == 'X') sprintf(sym, "H%s", p+1);
				else strcpy(sym, p);

				sprintf(str, " %-3s", sym);
				if (nzmat == 1) done = 1;
				break;
			case 1:
				if (!zmat->atomname[1][0]) sprintf(str, " %-4d", zmat->zdef[0]);
				else sprintf(str, " %-4s", zmat->atomname[1]);
				break;
			case 2:
				if (!zmat->varname[0][0]) sprintf(str, "   % 9.6f", zmat->zval[0]);
				else sprintf(str, "   %-9s", zmat->varname[0]);
				if (nzmat == 2) done = 1;
				break;
			case 3:
				if (!zmat->atomname[2][0]) sprintf(str, " %-4d", zmat->zdef[1]);
				else sprintf(str, "   %-4s", zmat->atomname[2]);
				break;
			case 4:
				if (!zmat->varname[1][0]) sprintf(str, "   %10.6f", zmat->zval[1]);
				else sprintf(str, "   %-10s", zmat->varname[1]);
				if (nzmat == 3) done = 1;
				break;
			case 5:
				if (!zmat->atomname[3][0]) sprintf(str, " %-4d", zmat->zdef[2]);
				else sprintf(str, "   %-4s", zmat->atomname[3]);
				break;
			case 6:
				if (!zmat->varname[2][0]) sprintf(str, "   % 11.6f", zmat->zval[2]);
				else sprintf(str, "   %-11s", zmat->varname[2]);
				break;
			}
			if (str[0]) (*print_func)(dest, str);
			if (done) break;
		}
		(*print_func)(dest, "\n");
	}
}

static void	GPrintJaguarSection (ListPtr top_list, char *title, void *dest, void (*print_func)(void *, char *))
{
	char	str[256];
	ListPtr	list;

	if (title) {
		sprintf(str, "&%s\n", title);
		(*print_func)(dest, str);
	}

	ForEachList(top_list, list) {
		sprintf(str, "%s\n", list->L_JAGUAR_LINE);
		(*print_func)(dest, str);
	}

	if (title) (*print_func)(dest, "&\n");
}

static void	GPrintJaguarZmatVar (JaguarZmatPtr j, void *dest, void (*print_func)(void *, char *))
{
	int	i, nzmat, flags;
	char	str[128];
	ZmatPtr	zmat;

	(*print_func)(dest, "&zvar\n");
	for(zmat=j->zmatlist,nzmat=0;zmat;zmat=zmat->next,nzmat++) {
		for(i=0;i<3;i++) {
			if (!zmat->varname[i][0]) continue;
			if (j->flags & JAGUAR_CARTESIAN) {
				if (IsRedundantCartZmatVar(j->zmatlist, nzmat, zmat, i)) continue;
			} else {
				switch (i) {
				case 0: flags = Z_VAR_LENGTH; break;
				case 1: flags = Z_VAR_ANGLE;  break;
				case 2: flags = Z_VAR_DIHED;  break;
				}
				if (IsRedundantIntZmatVar(j->zmatlist, nzmat, zmat, zmat->flags & flags)) continue;
			}
			sprintf(str, "%-8s = % 13.6f\n", zmat->varname[i], zmat->zval[i]);
			(*print_func)(dest, str);
		}
	}
	(*print_func)(dest, "&\n");
}

static void	GPrintJaguarZmat (JaguarZmatPtr j, void *dest, void (*print_func)(void *, char *))
{
	if (!j || !j->zmatlist) return;

	/* print header */
	if (j->header) GPrintJaguarSection(j->header, NULL, dest, print_func);

	/* print zmat */
	(*print_func)(dest, "&zmat\n");
	if (j->flags & JAGUAR_CARTESIAN) GPrintJaguarCartZmat(j, dest, print_func);
	else GPrintJaguarIntZmat(j, dest, print_func);
	(*print_func)(dest, "&\n");

	/* print zvar, etc */
	GPrintJaguarZmatVar(j, dest, print_func);

	/* other sections */
	if (j->gen)     GPrintJaguarSection(j->gen, "gen", dest, print_func);
	if (j->gvb)     GPrintJaguarSection(j->gvb, "gvb", dest, print_func);
	if (j->lmp2)    GPrintJaguarSection(j->lmp2, "lmp2", dest, print_func);
	if (j->atomic)  GPrintJaguarSection(j->atomic, "atomic", dest, print_func);
	if (j->vdw)     GPrintJaguarSection(j->vdw, "vdw", dest, print_func);
	if (j->hess)    GPrintJaguarSection(j->hess, "hess", dest, print_func);
	if (j->guess)   GPrintJaguarSection(j->guess, "guess", dest, print_func);
	if (j->pointch) GPrintJaguarSection(j->pointch, "pointch", dest, print_func);
	if (j->efields) GPrintJaguarSection(j->efields, "efields", dest, print_func);
	if (j->ham)     GPrintJaguarSection(j->ham, "ham", dest, print_func);
	if (j->orbman)  GPrintJaguarSection(j->orbman, "orbman", dest, print_func);
	if (j->echo)    GPrintJaguarSection(j->echo, "echo", dest, print_func);
	if (j->path)    GPrintJaguarSection(j->path, "path", dest, print_func);
}

void	SPrintJaguarZmat (char **str_sum, JaguarZmatPtr j)
{
	ContinuePrintString(*str_sum);
	GPrintJaguarZmat(j, NULL, SPrintString);
	*str_sum = EndPrintString();
}

void	FPrintJaguarZmat (FILE *fp, JaguarZmatPtr j)
{
	if (!fp || !j) return;
	GPrintJaguarZmat(j, (void *)fp, FPrintString);
}
