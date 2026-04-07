#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

MopacZmatPtr	FLoadMopacZmat (FILE *fp)
{
	int	n, nzmat, maxlen, ncomments, z3flag, count, done, flags, nlines;
	char	*c, *p, str[256], line[256];
	ZmatPtr	zmat;
	MopacZmatPtr	mopac;

	clearmolerror();
	if (!fp) return NULL;
	if (!(mopac = (MopacZmatPtr)calloc(1, sizeof(MopacZmat))))
		return (MopacZmatPtr)setmolerrorno(MERR_MEM);
	mopac->flags |= GAUSS_MOPAC_ZMAT;

	/* route cards */
	ncomments = 2;
	mopac->route[0] = '\0';
	maxlen = sizeof(mopac->route);
	while (1) {
		if (!READLINE(fp, str)) {
			setmolerrorno(MERR_EOF);
			break;
		}
		if (strlen(str)+strlen(mopac->route) < maxlen) strcat(mopac->route, str);
		if (strchr(str, '+')) continue;	/* addition of keyword line */
		if (strchr(str, '&')) {	/* turn next line into keyword */
			ncomments--;
			continue;
		}
		break;
	}
	if (getmolerrorno() != -1) {
		FreeGaussZmat(&mopac);
		return NULL;
	}

	/* read comments */
	mopac->title[0] = '\0';
	maxlen = sizeof(mopac->title);
	while (ncomments > 0) {
		if (!READLINE(fp, str)) {
			setmolerrorno(MERR_EOF);
			break;
		}
		if (strlen(str)+strlen(mopac->title) < maxlen) strcat(mopac->title, str);
		ncomments--;
	}
	if (getmolerrorno() != -1) {
		FreeGaussZmat(&mopac);
		return NULL;
	}

	TrimTail(mopac->route, " \t\n");
	TrimTail(mopac->title, " \t\n");

	strcpy(line, mopac->route);
	for(c=line;(p = strtok(c, " ,\t\n"));c=(char *)NULL) {
		/*
		if (strcasecmp(p, "XYZ") == 0) {
			mopac->flags |= GAUSS_CARTESIAN;
			break;
		}
		*/
	}

	mopac->charge = 0;
	mopac->multiplicity = 1;

	/* read Z-matrix (terminated by a blank line or a center number 0) */
	nzmat = 0;
	z3flag = 0;

	while (1) {
		if (!READLINE(fp, str)) {
			/*
			if (!(mopac->flags & GAUSS_CARTESIAN)) setmolerrorno(MERR_EOF);
			*/
			break;
		}
		if ((c = strrchr(str, '\n'))) *c = '\0';
		setmolerrortext(str);
		strcpy(line, str+strspn(str, " \t"));
		if (BLANKLINE(line) || line[0] == '0') {	/* End of Z-matrix */
			nlines = 0;
			while (READLINE(fp, str)) {
				nlines++;
				if (nlines > 1000) break;	/* something is wrong */
				if (strstr(str, "$EndCoord:") || strstr(str, "EndGaussZmat:")) break;
				AddtoString(&mopac->data, str);
			}
			break;
		}
		if (strncmp(line, "$EndCoord:", 10)==0 || strncmp(line, "EndGaussZmat:",14) ==0) break;

		if (!(zmat = EnterNewZmat(&mopac->zmatlist))) {
			setmolerrorno(MERR_MEM);
			break;
		}
		nzmat++;
		zmat->atom = nzmat;

		/* break up the line into tokens */
		for(c=line,count=0,done=0;(p = strtok(c, " ,\t\n"));count++,c=(char *)NULL) {
			switch (count) {
			case 0:	/* center (atomic symbol optionally follwed by number or atomic number) */
				strcpy(zmat->atomname[0], p);
				ParseZmatCenter(p, mopac->zmatlist, zmat);
				break;
			case 1:	/* bond length (floating-point number or integer) */
			case 3:	/* bond angle */
			case 5:	/* dihedral */
				n = (count==1) ? 0 : (count==3 ? 1 : 2);
				if (IsNumber(p)) sscanf(p, "%lf", &zmat->zval[n]);
				else setmolerrorno(MERR_INVALID_NUMBER);
				break;
			case 2:	/* variable(1),  const(0), rxn path (-1) */
			case 4:
			case 6:
				if (!IsInt(p) && nzmat > 4) {
					setmolerrorno(MERR_FORMAT);
					break;
				}
				flags = (count==2) ? Z_VAR_LENGTH : ((count==4) ? Z_VAR_ANGLE : Z_VAR_DIHED);
				if (*p == '1') zmat->flags |= flags;
				if (nzmat == 3) sscanf(p, "%d", &z3flag);
				if (mopac->flags & GAUSS_CARTESIAN) {
					done = 1;
					break;
				}
				break;

			case 7:	/* connectivity */
				if (done) break;
				if (nzmat <= 2) {
					zmat->zdef[0] = 1;
					done = 1;
					break;
				}
				if (!IsInt(p) && nzmat > 4) {
					setmolerrorno(MERR_FORMAT);
					break;
				}
				sscanf(p, "%d", &zmat->zdef[0]);
				break;
			case 8:
				if (done) break;
				if (!IsInt(p) && nzmat > 4) {
					setmolerrorno(MERR_FORMAT);
					break;
				}
				sscanf(p, "%d", &zmat->zdef[1]);

				if (nzmat != 3) break;

				/*
				 * MOPAC treats the 3rd Z-matrix as a special case.
				 * If a variable type for dihedral (i.e. 7th item on a line) is 2,
				 * then atom1 = 1 and atom2 = 2,
				 * else if it is 1 and abs(6th item - 2.0) < 1.0E-4,
				 * then atom1 = 2 and atom1 = 1.
				 */
				if (z3flag == 2) {
					zmat->zdef[0] = 1;
					zmat->zdef[1] = 2;
				} else if (z3flag == 1 && ABS(zmat->zval[2]-2.0) < 1.0E-4) {
					zmat->zdef[0] = 2;
					zmat->zdef[1] = 1;
				}
				done = 1;
				break;
			case 9:
				if (done) break;
				if (!IsInt(p) && nzmat > 4) {
					setmolerrorno(MERR_FORMAT);
					break;
				}
				sscanf(p, "%d", &zmat->zdef[2]);
				break;

			case 10:	/* assume charge --- but Zmat doesn't have a charge field... */
				break;
			case 11:
				strncpy(zmat->label, p, Z_LABEL_LEN);
				zmat->label[Z_LABEL_LEN] = '\0';
				break;
			case 12:
				if (IsInt(p)) {
					sscanf(p, "%d", &zmat->resseq);
					if (zmat->resseq < 0) zmat->resseq = 0;
				}
				break;
			case 13:
				strncpy(zmat->resname, p, Z_RES_LEN);
				zmat->label[Z_RES_LEN] = '\0';
				break;
			}

			if (getmolerrorno() != -1) break;
		}
		if (getmolerrorno() != -1) break;

		/* Check if this is a Cartesian coord */
		if (mopac->flags & GAUSS_CARTESIAN) continue;
		if (nzmat <= 3) continue;

		if (nzmat > 3 && zmat->zdef[0] == 0 && zmat->zdef[1] == 0 && zmat->zdef[0] == 0) {
			mopac->flags |= GAUSS_CARTESIAN;
			continue;
		}

		/* error checking */
		if (!CheckInternalCoord(nzmat, zmat->zdef[0], zmat->zdef[1], zmat->zdef[2],
			zmat->zval[0], zmat->zval[1], zmat->zval[2])) {
			setmolerrorno(MERR_INTERNAL_COORD);
			break;
		}
	}
	if (getmolerrorno() != -1) {
		FreeGaussZmat(&mopac);
		return NULL;
	}

	/* assign variable names (arbitrary) */
	for(zmat=mopac->zmatlist,nzmat=1;zmat;zmat=zmat->next,nzmat++) {
		zmat->varname[0][0] = '\0';
		zmat->varname[1][0] = '\0';
		zmat->varname[2][0] = '\0';
		if (nzmat > 1) sprintf(zmat->varname[0], "R%0d", nzmat);
		if (nzmat > 2) sprintf(zmat->varname[1], "A%0d", nzmat);
		if (nzmat > 3) sprintf(zmat->varname[2], "D%0d", nzmat);
	}

	if (mopac->flags & GAUSS_CARTESIAN) {
		ForEachZmat(mopac->zmatlist,zmat) {
			zmat->cart[0] = zmat->zval[0];
			zmat->cart[1] = zmat->zval[1];
			zmat->cart[2] = zmat->zval[2];
			zmat->zval[0] = zmat->zval[1] = zmat->zval[2] = 0.0;
			zmat->zdef[0] = zmat->zdef[1] = zmat->zdef[2] = 0;
		}
	} else {
		zmat = mopac->zmatlist;
		if (zmat) {
			zmat->zdef[0] = zmat->zdef[1] = zmat->zdef[2] = 0;
			zmat = zmat->next;
		}
		if (zmat) {
			zmat->zdef[0] = 1;
			zmat->zdef[1] = zmat->zdef[2] = 0;
			zmat = zmat->next;
		}
		if (zmat) {
			zmat->zdef[2] = 0;
			if (zmat->zdef[0] == 0 || zmat->zdef[1] == 0) {
				zmat->zdef[0] = 2;
				zmat->zdef[1] = 1;
			}
		}

		ZmatToCartesian(mopac->zmatlist);
	}

	return mopac;
}

#define DEF_MOPAC_ROUTE	"AM1 T=600 COMPFG DEBUG LARGE GEO-OK MMOK"

static char	mopac_route[512] = DEF_MOPAC_ROUTE;
char	*GetMopacRoute (void)
{
	return mopac_route;
}

void	GPrintMopacZmat (MopacZmatPtr mopac, void *dest, void (*print_func)(void *, char *))
{
	ZmatPtr	zmat;
	char	str[512], *sym;
	int	center, varlen, varang, vardih;
	double	val[3];

	if (!mopac) return;

	if (mopac->flags & GAUSS_MOPAC_ZMAT) sprintf(str, "%s", mopac->route);
	else sprintf(str, GetMopacRoute());

	if ((mopac->flags & GAUSS_CARTESIAN) && !strstr(str, "XYZ"))
		strcat(str, " XYZ");
	strcat(str, "\n");
	(*print_func)(dest, str);

	if (strcmp(mopac->title, Z_GAUSS_DEFAULT_TITLE) == 0 ||
	    strcmp(mopac->title, Z_BOSS_DEFAULT_TITLE ) == 0 ||
	    strcmp(mopac->title, Z_MCPRO_DEFAULT_TITLE) == 0) {
		strcpy(str, "\nMOPAC Z-MATRIX\n");
	} else {
		sprintf(str, "%s%s\n", strchr(mopac->title, '\n') ? "" : "\n", mopac->title);
	}
	(*print_func)(dest, str);

	for(zmat=mopac->zmatlist,center=1;zmat;zmat=zmat->next,center++) {
		if (zmat->an == 0) sym = "H";
		else if (zmat->an < 0 || zmat->an >= 98) sym = "XX";
		else sym = GetElemSymbol(zmat->an);
		varlen = (!(zmat->flags & Z_VAR_LENGTH) || center == 1) ? 0 : 1;
		varang = (!(zmat->flags & Z_VAR_ANGLE)  || center <= 2) ? 0 : 1;
		vardih = (!(zmat->flags & Z_VAR_DIHED)  || center <= 3) ? 0 : 1;

		if (zmat->zdef[0] == 0 && zmat->zdef[1] == 0 && zmat->zdef[2] == 0) {
			/* Cartesian input (Connectivity = 0 0 0) */
			val[0] = zmat->cart[0];
			val[1] = zmat->cart[1];
			val[2] = zmat->cart[2];
		} else if (center == 1) {	/* first entry */
			val[0] = val[1] = val[2] = 0;
		} else {
			val[0] = zmat->zval[0];
			val[1] = zmat->zval[1];
			val[2] = zmat->zval[2];
		}
		sprintf(str, "%-2s %12.6f %1d %11.6f %1d % 11.6f %1d %5d %3d %3d",
			sym, val[0], varlen, val[1], varang, val[2], vardih,
			zmat->zdef[0], zmat->zdef[1], zmat->zdef[2]);
		(*print_func)(dest, str);

		sprintf(str, "   0.0 %4s %4d %3s\n", zmat->label, zmat->resseq, zmat->resname);
		(*print_func)(dest, str);
	}
	sprintf(str, "%2d %12.6f %1d %11.6f %1d % 11.6f %1d %5d %3d %3d\n", 0, 0.0, 0, 0.0, 0, 0.0, 0, 0, 0, 0);
	(*print_func)(dest, str);

	if (mopac->data && (mopac->flags & GAUSS_MOPAC_ZMAT)) (*print_func)(dest, mopac->data);
}

void	SPrintMopacZmat (char **str_sum, MopacZmatPtr mopac)
{
	ContinuePrintString(*str_sum);
	GPrintMopacZmat(mopac, NULL, SPrintString);
	*str_sum = EndPrintString();
}

void	FPrintMopacZmat (FILE *fp, MopacZmatPtr mopac)
{
	if (!fp || !mopac) return;
	GPrintMopacZmat(mopac, (void *)fp, FPrintString);
}

