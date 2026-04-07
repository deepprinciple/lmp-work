#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <unistd.h>

#include "macros.h"
#include "molecule.h"

#define ISVALIDCENTER(n, max)	(n > 0 && n <= max)
#define ISVALIDBONDLEN(len)	(len > 1.0E-6)
#define ISVALIDBONDANG(ang)	(ang >= 0 && ang <= 180)

#define DEBUGLEVEL(n)	(n>0 && getmoldebuglevel() >= n)
#define PRINT_ERROR(n)	{if (DEBUGLEVEL(n) && getmolerrorno() != -1) printmolerror();}

static int	getbossversion (char *str)
{
	return BOSS_VMAX-1;
}

static ListPtr	readrangedcenter (FILE *fp, int maxcenter)
{
	int	n1, n2, i;
	char	str[256], s1[8], s2[8], *c;
	ListPtr	list, centerlist=NULL;

	clearmolerror();
	while (1) {
		if (!READLINE(fp, str)) {setmolerrorno(MERR_EOF); break;}
		if (strncmp(str, "    ", 4) == 0) break;
		if ((c = strrchr(str, '\n'))) *c = '\0';
		if (str[4] == '-') {
			if (sscanf(str, "%4c-%4c", s1, s2) != 2 ||
			   !sscanf(s1, "%d", &n1) ||
			   !sscanf(s2, "%d", &n2) ||
			   !ISVALIDCENTER(n1, maxcenter) ||
			   !ISVALIDCENTER(n2, maxcenter) ||
			    n1 > n2) {
				setmolerrortext(str);
				setmolerrorno(MERR_INVALID_CENTER);
				break;
			}
			for(i=n1;i<=n2;i++) {
				if (!(list = EnterNewList(&centerlist))) {
					setmolerrorno(MERR_MEM);
					break;
				}
				list->L_CENTER = i;
			}
		} else {
			if (!sscanf(str, "%4c", s1) || !sscanf(s1, "%d", &n1) ||
			    !ISVALIDCENTER(n1, maxcenter)) {
				setmolerrortext(str);
				setmolerrorno(MERR_INVALID_CENTER);
				break;
			}
			if (!(list = EnterNewList(&centerlist))) {
				setmolerrorno(MERR_MEM);
				break;
			}
			list->L_CENTER = n1;
		}
	}
	if (getmolerrorno() != -1) FreeList(&centerlist);
	return centerlist;
}

ListPtr	FLoadBossList (FILE *fp, char *scanfmt[], int nvar)
{
	int	i;
	char	*c, str[256], *val[8];
	ListPtr	list, bosslist=NULL;

	while (1) {
		if (!READLINE(fp, str) || strncmp(str, "    ", 4) == 0) break;
		if (!(list = EnterNewList(&bosslist))) {
			setmolerrorno(MERR_MEM);
			break;
		}
		val[0] = list->str1;
		val[1] = list->str2;
		val[2] = list->str3;
		val[3] = list->str4;
		val[4] = list->str5;
		val[5] = list->str6;
		val[6] = list->str7;
		val[7] = list->str8;
		if ((c = strrchr(str, '\n'))) *c = '\0';
		c = str;
		for(i=0;i<nvar;i++) {
			c = strscan(c, scanfmt[i], val[i]);
			if (!*c) break;
		}
	}
	if (getmolerrorno() != -1) FreeList(&bosslist);
	return bosslist;
}

static int	parse_boss_zmat_component (int component, int nzmat, ListPtr list)
{
	clearmolerror();
	switch (component) {
	case BOSS_GEOMVAR:
		if (!sscanf(list->str1, "%d", &list->L_CENTER) ||
		    !sscanf(list->str2, "%d", &list->L_PARAM) ||
		    !sscanf(list->str3, "%f", &list->L_VAL) ||
		    !ISVALIDCENTER(list->L_CENTER, nzmat) ||
		    !(list->L_PARAM >= 1 && list->L_PARAM <= 3)) {
			setmolerrorno(MERR_FORMAT);
			break;
		}

		/* error checking */
		if (list->L_PARAM == 1 && !ISVALIDBONDLEN(list->L_VAL)) {
			setmolerrorno(MERR_INVALID_BONDLEN);
			break;
		}
		if (list->L_PARAM == 2 && !ISVALIDBONDANG(list->L_VAL)) {
			setmolerrorno(MERR_INVALID_BONDANG);
			break;
		}
		break;

	case BOSS_VARBOND:
	case BOSS_VARANGLE:
		if (!sscanf(list->str1, "%d", &list->L_CENTER) ||
		    !ISVALIDCENTER(list->L_CENTER, nzmat)) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;

	case BOSS_ADDBOND:
		if (!sscanf(list->str1, "%d", &list->L_ATOM1) ||
		    !sscanf(list->str2, "%d", &list->L_ATOM2) ||
		    !ISVALIDCENTER(list->L_ATOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM2, nzmat) ||
		     list->L_ATOM1 == list->L_ATOM2) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;
	
	case BOSS_HARMONIC:
		if (!sscanf(list->str1, "%d", &list->L_ATOM1) ||
		    !sscanf(list->str2, "%d", &list->L_ATOM2) ||
		    !sscanf(list->str3, "%f", &list->L_RIJ0) ||
		    !sscanf(list->str4, "%f", &list->L_K0) ||
		    !sscanf(list->str5, "%f", &list->L_K1) ||
		    !sscanf(list->str6, "%f", &list->L_K2) ||
		    !ISVALIDCENTER(list->L_ATOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM2, nzmat) ||
		     list->L_ATOM1 == list->L_ATOM2) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;

	case BOSS_ADDANGLE:
		if (!sscanf(list->str1, "%d", &list->L_ATOM1) ||
		    !sscanf(list->str2, "%d", &list->L_ATOM2) ||
		    !sscanf(list->str3, "%d", &list->L_ATOM3) ||
		    !ISVALIDCENTER(list->L_ATOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM2, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM3, nzmat) ||
		     list->L_ATOM1 == list->L_ATOM2 ||
		     list->L_ATOM1 == list->L_ATOM3 ||
		     list->L_ATOM2 == list->L_ATOM3) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;
	
	case BOSS_VARDIHED:
		if (!sscanf(list->str1, "%d", &list->L_CENTER) ||
		    !sscanf(list->str2, "%d", &list->L_ITYPE) ||
		    !sscanf(list->str3, "%d", &list->L_FTYPE) ||
		    !sscanf(list->str4, "%f", &list->L_RANGE) ||
		    !ISVALIDCENTER(list->L_CENTER, nzmat) ||
		     (int)list->L_ITYPE < 0 ||
		     (int)list->L_FTYPE < 0) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;
	
	case BOSS_ADDDIHED:
		if (!sscanf(list->str1, "%d", &list->L_ATOM1) ||
		    !sscanf(list->str2, "%d", &list->L_ATOM2) ||
		    !sscanf(list->str3, "%d", &list->L_ATOM3) ||
		    !sscanf(list->str4, "%d", &list->L_ATOM4) ||
		    !sscanf(list->str5, "%d", &list->L_ITYPE) ||
		    !sscanf(list->str6, "%d", &list->L_FTYPE) ||
		    !ISVALIDCENTER(list->L_ATOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM2, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM3, nzmat) ||
		    !ISVALIDCENTER(list->L_ATOM4, nzmat) ||
		     (int)list->L_ITYPE < 0 ||
		     (int)list->L_FTYPE < 0) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;
	
	case BOSS_DOMAIN:
		if (!sscanf(list->str1, "%d", &list->L_FIRST_DOM1) ||
		    !sscanf(list->str2, "%d", &list->L_LAST_DOM1) ||
		    !sscanf(list->str3, "%d", &list->L_FIRST_DOM2) ||
		    !sscanf(list->str4, "%d", &list->L_LAST_DOM2) ||
		    !ISVALIDCENTER(list->L_FIRST_DOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_LAST_DOM1, nzmat) ||
		    !ISVALIDCENTER(list->L_FIRST_DOM2, nzmat) ||
		    !ISVALIDCENTER(list->L_LAST_DOM2, nzmat) ||
		     list->L_LAST_DOM1 < list->L_FIRST_DOM1 ||
		     list->L_LAST_DOM2 < list->L_FIRST_DOM2) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;

	case BOSS_CONFORM:
		if (!sscanf(list->str1, "%d", &list->L_CENTER) ||
		    !sscanf(list->str2, "%d", &list->L_PARAM) ||
		    !sscanf(list->str3, "%f", &list->L_RMIN) ||
		    !sscanf(list->str4, "%f", &list->L_RMAX) ||
		    !ISVALIDCENTER(list->L_CENTER, nzmat) ||
		    !(list->L_PARAM >= 1 && list->L_PARAM <= 3)) {
			setmolerrorno(MERR_FORMAT);
			break;
		}

		/* error checking */
		if (list->L_PARAM == 1 && (!ISVALIDBONDLEN(list->L_RMIN) || !ISVALIDBONDLEN(list->L_RMAX))) {
			setmolerrorno(MERR_INVALID_BONDLEN);
			break;
		}
		if (list->L_PARAM == 2 && (!ISVALIDBONDANG(list->L_RMIN) || !ISVALIDBONDANG(list->L_RMAX))) {
			setmolerrorno(MERR_INVALID_BONDANG);
			break;
		}
		break;

	case BOSS_LOCALHEAT:
		if (!sscanf(list->str1, "%d", &list->L_RESIDUE)) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		break;


	default:
		return 0;
	}

	if (getmolerrorno() != -1) return 0;
	else return 1;
}

static ListPtr	FLoadBossZmatComponent (FILE *fp, int component, int nzmat, int *flags_ret)
{
	int	i, nvar, n1, n2, flags;
	char	*c, s1[8], s2[8], str[256], *val[8];
	ListPtr	list, bosslist=NULL;

	static char	_sf[8][8],
		*scanfmt[8] = {_sf[0], _sf[1], _sf[2], _sf[3], _sf[4], _sf[5], _sf[6], _sf[7]};

	*flags_ret = 0;
	nvar = 0;
	flags = 0;

	switch (component) {
	case BOSS_GEOMVAR:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%12c"); nvar++;
		break;
	case BOSS_VARBOND:
	case BOSS_VARANGLE:
		strcpy(scanfmt[0], "%4c"); nvar++;
		break;
	case BOSS_ADDBOND:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		break;
	case BOSS_HARMONIC:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%10c"); nvar++;
		strcpy(scanfmt[3], "%10c"); nvar++;
		strcpy(scanfmt[4], "%10c"); nvar++;
		strcpy(scanfmt[5], "%10c"); nvar++;
		break;
	case BOSS_ADDANGLE:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%4c"); nvar++;
		break;

	case BOSS_VARDIHED:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%4c"); nvar++;
		strcpy(scanfmt[3], "%12c"); nvar++;
		break;

	case BOSS_ADDDIHED:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%4c"); nvar++;
		strcpy(scanfmt[3], "%4c"); nvar++;
		strcpy(scanfmt[4], "%4c"); nvar++;
		strcpy(scanfmt[5], "%4c"); nvar++;
		break;
	
	case BOSS_DOMAIN:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%4c"); nvar++;
		strcpy(scanfmt[3], "%4c"); nvar++;
		break;
	
	case BOSS_CONFORM:
		strcpy(scanfmt[0], "%4c"); nvar++;
		strcpy(scanfmt[1], "%4c"); nvar++;
		strcpy(scanfmt[2], "%12c"); nvar++;
		strcpy(scanfmt[3], "%12c"); nvar++;
		break;

	case BOSS_LOCALHEAT:
		strcpy(scanfmt[0], "%4c"); nvar++;
		break;

	default:
		return NULL;
	}

	while (1) {
		if (!READLINE(fp, str) || strncmp(str, "    ", 4) == 0) break;
		if ((c = strrchr(str, '\n'))) *c = '\0';
		setmolerrortext(str);
		if (strncasecmp(str, "AUTO", 4) == 0) {
			switch (component) {
			case BOSS_ADDBOND:  flags |= B_AUTO_ADDBOND;  break;
			case BOSS_ADDANGLE: flags |= B_AUTO_ADDANGLE; break;
			case BOSS_ADDDIHED: flags |= B_AUTO_ADDDIHED; break;
			}
		} else if (str[4] == '-') {
			s1[4] = s2[4] = '\0';
			if (sscanf(str, "%4c-%4c", s1, s2) != 2 ||
				!sscanf(s1, "%d", &n1) ||
				!sscanf(s2, "%d", &n2) ||
				!ISVALIDCENTER(n1, nzmat) ||
				!ISVALIDCENTER(n2, nzmat) ||
				n1 > n2) {
				setmolerrorno(MERR_INVALID_CENTER);
				break;
			}
			for(i=n1;i<=n2;i++) {
				if (!(list = EnterNewList(&bosslist))) {
					setmolerrorno(MERR_MEM);
					break;
				}
				list->L_CENTER = i;
			}
		} else {
			if (!(list = EnterNewList(&bosslist))) {
				setmolerrorno(MERR_MEM);
				break;
			}
			val[0] = list->str1;
			val[1] = list->str2;
			val[2] = list->str3;
			val[3] = list->str4;
			val[4] = list->str5;
			val[5] = list->str6;
			val[6] = list->str7;
			val[7] = list->str8;
			if ((c = strrchr(str, '\n'))) *c = '\0';
			c = str;
			for(i=0;i<nvar;i++) {
				c = strscan(c, scanfmt[i], val[i]);
				if (!*c) break;
			}

			/* parse */
			if (!parse_boss_zmat_component(component, nzmat, list)) {
				/* remove list */
				RemoveList(list, &bosslist);
				setmolerrorno(MERR_FORMAT);
			}

			if (component == BOSS_VARDIHED) flags |= B_VARDIHED_EXPLICIT;
		}
	}
	*flags_ret = flags;

	if (flags & B_AUTO_ADD) {
		if (bosslist) FreeList(&bosslist);
		return NULL;
	} else return bosslist;
}

void	FreeBossZmat (BossZmatPtr *bp)
{
	BossZmatPtr	b;
	if (!bp || !(b = *bp)) return;
	if (b->zmatlist)     FreeZmatList(&b->zmatlist);
	if (b->atomtypes)    FreeList(&b->atomtypes);
	if (b->torsiontypes) FreeList(&b->torsiontypes);

	if (b->geomvar)  FreeList(&b->geomvar);
	if (b->varbond)  FreeList(&b->varbond);
	if (b->addbond)  FreeList(&b->addbond);
	if (b->harmonic) FreeList(&b->harmonic);
	if (b->varangle) FreeList(&b->varangle);
	if (b->addangle) FreeList(&b->addangle);
	if (b->vardihed) FreeList(&b->vardihed);
	if (b->adddihed) FreeList(&b->adddihed);
	if (b->domain)   FreeList(&b->domain);
	if (b->conform)  FreeList(&b->conform);
	if (b->exclude)  FreeList(&b->exclude);

	if (b->param) FreeBossPar(&b->param);

	free((char *)b);
	*bp = NULL;
}

BossZmatPtr	CopyBossZmat (BossZmatPtr oldboss)
{
	BossZmatPtr	boss;

	if (!oldboss) return NULL;
	if (!(boss = (BossZmat *)calloc(1, sizeof(BossZmat)))) return NULL;
	BCOPY(oldboss, boss, sizeof(BossZmat));

	boss->zmatlist = CopyZmatList(oldboss->zmatlist);

	boss->atomtypes    = CopyList(oldboss->atomtypes);
	boss->torsiontypes = CopyList(oldboss->torsiontypes);

	boss->geomvar  = CopyList(oldboss->geomvar);
	boss->varbond  = CopyList(oldboss->varbond);
	boss->addbond  = CopyList(oldboss->addbond);
	boss->harmonic = CopyList(oldboss->harmonic);
	boss->varangle = CopyList(oldboss->varangle);
	boss->addangle = CopyList(oldboss->addangle);
	boss->vardihed = CopyList(oldboss->vardihed);
	boss->adddihed = CopyList(oldboss->adddihed);
	boss->domain   = CopyList(oldboss->domain);
	boss->conform  = CopyList(oldboss->conform);
	boss->exclude  = CopyList(oldboss->exclude);

	boss->param = CopyBossPar(oldboss->param);

	return boss;
}

/*
 * Kinds of BOSS input formats --- specified in the parameter file
 *    1. Z-matrix
 *    2. PDB
 *    3. Mindtool
 */

BossZmatPtr	FLoadBossZmatInput (FILE *fp, int flags)
{
	int	n, serno, itype, ftype, z1, z2, z3, nzmat, bossversion, resseq, auto_flags;
	char	str[256], label[Z_LABEL_LEN+1], resname[Z_RES_LEN+1];
	double	len, theta, phi;
	ListPtr	list;
	ZmatPtr	zmat, zmatlist;
	BossZmatPtr	boss=NULL;

	clearmolerror();
	if (!fp) return NULL;
	if (!(boss = (BossZmat *)calloc(1, sizeof(BossZmat)))) return NULL;

	if (flags & B_MCPRO) boss->flags |= B_MCPRO;
	if (flags & B_QM) boss->flags |= B_QM;
	if (flags & B_SOLVENT_ZMAT) boss->flags |= B_SOLVENT_ZMAT;

	/* first line (title) */
	if (!READLINE(fp, str)) return (BossZmatPtr)setmolerrorno(MERR_EOF);
	sscanf(str, "%[^\n]", boss->title);

	nzmat = 0;
	zmatlist = NULL;
	bossversion = 0;
	while (1) {
		if (!READLINE(fp, str)) {setmolerrorno(MERR_EOF); break;}
		setmolerrortext(str);

		if (strncasecmp(str, "$EndCoord:", 10) == 0) break;
		if (strncasecmp(str, "EndBossZmat:", 12) == 0) break;
		if (nzmat > 0 && strncmp(str, "    ", 4) == 0) break;	/* end of Z-matrix */
		if (!bossversion && !(bossversion = getbossversion(str))) {
			setmolerrorno(MERR_FORMAT);
			break;
		}
		if (strncasecmp(str, "TERZ", 3) == 0) {	/* new solute */
			zmatlist = NULL;
			continue;
		}
		label[0] = '\0';
		resname[0] = '\0';
		resseq = 0;
		strscan(str, "%4d %4c%4d %4d %4d%12lf%4d%12lf%4d%12lf% %3c%5d",
			&serno, label, &itype, &ftype, &z1, &len, &z2, &theta, &z3, &phi,
			resname, &resseq);
		TrimTail(resname, " \t\n");

		/* error checking */
if (itype == 0 || (itype < 0 && itype != -1)) {
	fprintf(stderr, "WARNING: Invalid atom type %d in Z-matrix: %s", itype, str);
}

		if (!CheckInternalCoord(nzmat+1, z1, z2, z3, len, theta, phi)) {
			setmolerrorno(MERR_INTERNAL_COORD);
			break;
		}
		if (itype != 0 && ftype == 0) ftype = itype;

		if (!(zmat = EnterNewZmat(&zmatlist))) {
			setmolerrorno(MERR_MEM);
			break;
		}
		if (!zmat->prev) {	/* first zmat */
			list = EnterNewList(&boss->zmatlist);
			list->L_ZMAT = (size_t)zmat;
		}
		nzmat++;
		zmat->atom = nzmat;

		strncpy(zmat->label, label, Z_LABEL_LEN);
		zmat->label[Z_LABEL_LEN] = '\0';
		strncpy(zmat->resname, resname, Z_RES_LEN);
		zmat->resname[Z_RES_LEN] = '\0';
		zmat->resseq = resseq;
		zmat->itype = itype;
		zmat->ftype = ftype;
		zmat->zdef[0] = z1;
		zmat->zdef[1] = z2;
		zmat->zdef[2] = z3;
		zmat->zval[0] = len;
		zmat->zval[1] = theta;
		zmat->zval[2] = phi;

		if (nzmat > 3) zmat->flags |= (phi > 0.0) ? Z_POS_DIHED : Z_NEG_DIHED;
	}

	if (getmolerrorno() != -1) {
		FreeBossZmat(&boss);
		return NULL;
	}

	/* calculate Cartesian coordinates */
	ZmatListToCartesian(boss->zmatlist);

if (!(boss->flags & B_SOLVENT_ZMAT)) {
	/* geometry variations */
	boss->geomvar = FLoadBossZmatComponent(fp, BOSS_GEOMVAR, nzmat, &auto_flags);
}

#ifdef junk
	if (getmolerrorno() != -1) {	/* set warning and just return */
		FreeList(&boss->geomvar);
		return boss;
	}
#endif

	/* variable bonds */
	boss->varbond = FLoadBossZmatComponent(fp, BOSS_VARBOND, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->varbond);
		return boss;
	}

	/* additional bonds */
	boss->addbond = FLoadBossZmatComponent(fp, BOSS_ADDBOND, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->addbond);
		return boss;
	}

if (!(boss->flags & B_SOLVENT_ZMAT)) {
	/* harmonic constraints */
	boss->harmonic = FLoadBossZmatComponent(fp, BOSS_HARMONIC, nzmat, &auto_flags);
	if (getmolerrorno() != -1) {
		FreeList(&boss->harmonic);
		return boss;
	}
}

	/* variable bond angles */
	boss->varangle = FLoadBossZmatComponent(fp, BOSS_VARANGLE, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->varangle);
		return boss;
	}

	/* additional bond angles */
	boss->addangle = FLoadBossZmatComponent(fp, BOSS_ADDANGLE, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->addangle);
		return boss;
	}

	/* variable dihedral angles */
	boss->vardihed = FLoadBossZmatComponent(fp, BOSS_VARDIHED, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->vardihed);
		return boss;
	}

	/* additional dihedral angles */
	boss->adddihed = FLoadBossZmatComponent(fp, BOSS_ADDDIHED, nzmat, &auto_flags);
	if (auto_flags) boss->flags |= auto_flags;
	if (getmolerrorno() != -1) {
		FreeList(&boss->adddihed);
		return boss;
	}

	/* domain definitions */
	boss->domain = FLoadBossZmatComponent(fp, BOSS_DOMAIN, nzmat, &auto_flags);
	if (getmolerrorno() != -1) {
		FreeList(&boss->domain);
		return boss;
	}

	/* conformational search (BOSS only) */
	if (!(boss->flags & B_MCPRO)) {
if (!(boss->flags & B_SOLVENT_ZMAT)) {
		boss->conform = FLoadBossZmatComponent(fp, BOSS_CONFORM, nzmat, &auto_flags);
		if (auto_flags) boss->flags |= auto_flags;
		if (getmolerrorno() != -1) {
			FreeList(&boss->conform);
			return boss;
		}

		/* local heating */
		boss->localheat = FLoadBossZmatComponent(fp, BOSS_LOCALHEAT, nzmat, &auto_flags);
		if (getmolerrorno() != -1) {
			FreeList(&boss->localheat);
			return boss;
		}
}

		/* check for "Final Non-Bonded Parameters for QM Atoms:" (BOSS4.1) */
		while (READLINE(fp, str)) {
			if (strncasecmp(str, "$EndCoord:", 10) == 0 || strncasecmp(str, "EndBossZmat:", 12) == 0) {
				break;
			}

			if (strstr(str, "Non-Bonded Parameters")) {
				READLINE(fp, str);
				boss->flags |= B_PARAM_IN_ZMAT;
				break;
			}
		}

	} else {	/* MCPRO excluded atoms */
		while (1) {
			int	i;
			char	*c;

			if (!READLINE(fp, str)) {setmolerrorno(MERR_EOF); break;}
			if (strncmp(str, "    ", 4) == 0) break;
			if ((c = strrchr(str, '\n'))) *c = '\0';
			setmolerrortext(str);
			c = str;
			for(i=0;i<10;i++) {
				c = strscan(c, "%4d", &n);
				if (!ISVALIDCENTER(n, nzmat)) {
					setmolerrorno(MERR_INVALID_CENTER);
					break;
				}
				if (!(list = EnterNewList(&boss->exclude))) {
					setmolerrorno(MERR_MEM);
					break;
				}
				list->L_CENTER = n;
				if (!*c) break;
			}
			if (getmolerrorno() != -1) {
				FreeList(&boss->exclude);
				return boss;
			}
		}
	}

	return boss;
}

BossZmatPtr	CreateBossZmatFromZmat (ZmatPtr zmatlist, int flags)
{
	BossZmatPtr	boss=NULL;
	ZmatPtr	newzmatlist;

	if (!zmatlist) return NULL;
	if (!(boss = (BossZmat *)calloc(1, sizeof(BossZmat)))) return NULL;
	if (flags & B_MCPRO) boss->flags |= B_MCPRO;
	sprintf(boss->title, "%s Z-Matrix", (flags & B_MCPRO) ? "MCPRO" : "BOSS");
	if (!(newzmatlist = CopyZmat(zmatlist))) {
		FreeBossZmat(&boss);
		return NULL;
	}
	boss->zmatlist = NewList();
	boss->zmatlist->L_ZMAT = (size_t)newzmatlist;

	return boss;
}

BossZmatPtr	CreateBossZmatFromGauss (GaussZmatPtr gauss)
{
	BossZmatPtr	boss;
	char	*p;

	if (!gauss || !gauss->zmatlist) return NULL;
	boss = CreateBossZmatFromZmat(gauss->zmatlist, 0);

	strcpy(boss->title, gauss->title);
	TrimStr(boss->title, "\n");
	if ((p = strchr(boss->title, '\n'))) *p = '\0';

	if (strcmp(boss->title, Z_GAUSS_DEFAULT_TITLE) == 0 ||
	    strcmp(boss->title, Z_MOPAC_DEFAULT_TITLE) == 0 ||
	    strlen(boss->title) == 0) {
		strcpy(boss->title, Z_BOSS_DEFAULT_TITLE);
	}

	return boss;
}

BossZmatPtr	CreateBossZmatFromJaguar (JaguarZmatPtr jaguar)
{
	BossZmatPtr	boss;

	if (!jaguar || !jaguar->zmatlist) return NULL;
	boss = CreateBossZmatFromZmat(jaguar->zmatlist, 0);

	return boss;
}

BossZmatPtr	CreateBossZmatFromMol (MolPtr mol, int flags, int autozmat_flags)
{
	/* autozmat_flags = 0                 --- do not add dummy
	 *                  Z_ADD_BGN_DUMMIES --- add two dummies at the center
	 */

	ZmatPtr	zmatlist;
	BossZmatPtr	boss;
	char	*title;

	if (!(zmatlist = CreateAutoZmatFromMol(mol, autozmat_flags))) return NULL;
	boss = CreateBossZmatFromZmat(zmatlist, flags);

	if (mol->header >= 0) {
		title = GetStrValue(mol->header);
		if (title && *title) {
			strcpy(boss->title, title);
			if ((title = strchr(boss->title, '\n'))) *title = '\0';
		}
	}
	return boss;
}

static void	gprint_boss_ranged_list (ListPtr list, void *dest, void (*print_func)(void *, char *))
{
	ListPtr	l;
	int	ibgn, iend;
	char	str[256];

	if (!list) return;

	for(l=list;l;) {
		ibgn = iend = l->L_CENTER;
		for(l=l->next;l;l=l->next) {
			if (l->L_CENTER != iend+1) break;
			iend++;
		}
		if (ibgn == iend) sprintf(str, "%4d\n", ibgn);
		else sprintf(str, "%04d-%04d\n", ibgn, iend);
		(*print_func)(dest, str);
	}
}

static void	GPrintBossZmat (BossZmatPtr boss, void *dest, void (*print_func)(void *, char *))
{
	ListPtr	l;
	ZmatPtr	zmat;
	int	n, nzmat, itype, ftype;
	char	*p, str[256], label[Z_LABEL_LEN+1];

	if (!boss) return;

	sprintf(str, "%s", boss->title);
	TrimStr(str, "\n");
	if ((p = strchr(str, '\n'))) *p = '\0';

	if (strcmp(str, Z_GAUSS_DEFAULT_TITLE) == 0 ||
	    strcmp(str, Z_MOPAC_DEFAULT_TITLE) == 0 ||
	    strlen(str) == 0) {
		sprintf(str, "%s", (boss->flags & B_MCPRO) ? Z_MCPRO_DEFAULT_TITLE : Z_BOSS_DEFAULT_TITLE);
	}
	strcat(str, "\n");
	(*print_func)(dest, str);

	nzmat = 0;
	ForEachList(boss->zmatlist, l) {
		if (l != boss->zmatlist) {
			sprintf(str, "TERZ\n");
			(*print_func)(dest, str);
		}
		ForEachZmat((ZmatPtr)l->L_ZMAT, zmat) {
			nzmat++;
			itype = (boss->flags & B_QM) ? zmat->an : zmat->itype;
			ftype = (boss->flags & B_QM) ?        0 : zmat->ftype;
			strncpy(label, zmat->label, Z_LABEL_LEN);

			if (zmat->an == 0) {
				if (itype == 0) {
					itype = 1;
					if (ftype != 0) ftype = 1;
				}
				strcpy(label, "H");
			}
			label[Z_LABEL_LEN] = '\0';
			TrimStr(label, " ");

			sprintf(str, "%4d %-4s%4d %4d %4d % 10.6f %4d % 10.5f %4d % 10.5f ",
				nzmat, label, itype, ftype, zmat->zdef[0], zmat->zval[0],
				zmat->zdef[1], zmat->zval[1], zmat->zdef[2], zmat->zval[2]);
			(*print_func)(dest, str);

			if (zmat->resname[0] && strcmp(zmat->resname, "UNK") != 0) {
				strncpy(label, zmat->resname, Z_RES_LEN);
				label[Z_RES_LEN] = '\0';
				TrimStr(label, " \t\n");
				sprintf(str, " %3s%5d\n", label, zmat->resseq);
			} else {
				strcpy(str, "\n");
			}

			(*print_func)(dest, str);
		}
	}

if (!(boss->flags & B_SOLVENT_ZMAT)) {
	/* geometry variation */
	sprintf(str, "                    Geometry Variations follow    (2I4,F12.6)\n");
	(*print_func)(dest, str);
	ForEachList(boss->geomvar, l) {
		sprintf(str, "%4d%4d% 11.5f0\n", l->L_CENTER, l->L_PARAM, l->L_VAL);
		(*print_func)(dest, str);
	}
}

	sprintf(str, "                    Variable Bonds follow         (I4)\n");
	(*print_func)(dest, str);
	gprint_boss_ranged_list(boss->varbond, dest, print_func);

	sprintf(str, "                    Additional Bonds follow       (2I4)\n");
	(*print_func)(dest, str);
	if (boss->flags & B_AUTO_ADDBOND) (*print_func)(dest, "AUTO\n");
	else {
		ForEachList(boss->addbond, l) {
			sprintf(str, "%4d%4d\n", l->L_ATOM1, l->L_ATOM2);
			(*print_func)(dest, str);
		}
	}

if (!(boss->flags & B_SOLVENT_ZMAT)) {
	sprintf(str, "                    Harmonic Constraints follow   (2I4,4F10.4)\n");
	(*print_func)(dest, str);
	ForEachList(boss->harmonic, l) {
		sprintf(str, "%4d%4d% 10.4f% 10.4f% 10.4f% 10.4f\n",
		l->L_ATOM1, l->L_ATOM2, l->L_RIJ0, l->L_K0, l->L_K1, l->L_K2);
		(*print_func)(dest, str);
	}
}
	
	sprintf(str, "                    Variable Bond Angles follow   (I4)\n");
	(*print_func)(dest, str);
	gprint_boss_ranged_list(boss->varangle, dest, print_func);

	sprintf(str, "                    Additional Bond Angles follow (3I4)\n");
	(*print_func)(dest, str);

	if (boss->flags & B_AUTO_ADDANGLE) (*print_func)(dest, "AUTO\n");
	else {
		ForEachList(boss->addangle, l) {
			sprintf(str, "%4d%4d%4d\n", l->L_ATOM1, l->L_ATOM2, l->L_ATOM3);
			(*print_func)(dest, str);
		}
	}

	sprintf(str, "                    Variable Dihedrals follow     (3I4,F12.6)\n");
	(*print_func)(dest, str);
	if ((boss->flags & B_QM) || !(boss->flags & B_VARDIHED_EXPLICIT)) {
		gprint_boss_ranged_list(boss->vardihed, dest, print_func);
	} else {
		ForEachList(boss->vardihed, l) {
			sprintf(str, "%4d%4d%4d% 11.5f0\n", l->L_CENTER, l->L_ITYPE, l->L_FTYPE, l->L_RANGE);
			(*print_func)(dest, str);
		}
	}

	sprintf(str, "                    Additional Dihedrals follow   (6I4)\n");
	(*print_func)(dest, str);
	if (boss->flags & B_AUTO_ADDDIHED) (*print_func)(dest, "AUTO\n");
	else {
		ForEachList(boss->adddihed, l) {
			sprintf(str, "%4d%4d%4d%4d%4d%4d\n",
			l->L_ATOM1, l->L_ATOM2, l->L_ATOM3, l->L_ATOM4, l->L_ITYPE, l->L_FTYPE);
			(*print_func)(dest, str);
		}
	}

	sprintf(str, "                    Domain Definitions follow     (4I4)\n");
	(*print_func)(dest, str);
	ForEachList(boss->domain, l) {
		sprintf(str, "%4d%4d%4d%4d\n",
			l->L_FIRST_DOM1, l->L_LAST_DOM1, l->L_FIRST_DOM2, l->L_LAST_DOM2);
		(*print_func)(dest, str);
	}

	if (boss->flags & B_MCPRO) {
		/* MCPRO only */
		sprintf(str, "                    Excluded Atoms List follows   (10I4)\n");
		(*print_func)(dest, str);

		n = 0;
		ForEachList(boss->exclude, l) {
			if (n && (n%10)==0) {
				sprintf(str, "\n");
				(*print_func)(dest, str);
			}
			sprintf(str, "%4d", l->L_CENTER);
			(*print_func)(dest, str);
			n++;
		}
		if (n) {
			sprintf(str, "\n");
			(*print_func)(dest, str);
		}
	} else {
if (!(boss->flags & B_SOLVENT_ZMAT)) {
		/* BOSS only */
		sprintf(str, "                    Conformational Search (2I4,2F12.6)\n");
		(*print_func)(dest, str);
		if (boss->flags & B_AUTO_CSEARCH) (*print_func)(dest, "AUTO");
		else {
			ForEachList(boss->conform, l) {
				sprintf(str, "%4d%4d% 11.5f0% 11.5f0\n",
					l->L_CENTER, l->L_PARAM, l->L_RMIN, l->L_RMAX);
				(*print_func)(dest, str);
			}
		}

		/* local heating */
		sprintf(str, "                    Local Heating Residues follow (I4 or I4-I4)\n");
		(*print_func)(dest, str);
		gprint_boss_ranged_list(boss->localheat, dest, print_func);
}

	}

	sprintf(str, "                    Final blank line\n");
	(*print_func)(dest, str);
}

void	SPrintBossZmat (char **str_sum, BossZmatPtr boss)
{
	BgnPrintString();
	GPrintBossZmat(boss, (void *)NULL, SPrintString);
	*str_sum = EndPrintString();
}

void	FPrintBossZmat (FILE *fp, BossZmatPtr boss)
{
	if (!fp || !boss) return;
	GPrintBossZmat(boss, (void *)fp, FPrintString);
}

#define NEED_PRINT_BOSS_ZMAT_TORSION(boss)	\
	(((boss->flags & B_VARDIHED_EXPLICIT) && boss->vardihed) || \
	(!(boss->flags & B_AUTO_ADDDIHED) && boss->adddihed))

void	GPrintBossZmatTypes (BossZmatPtr boss, void *dest, void (*print_func)(void *, char *))
{
	int	flags=0;
	char	str[256], amber[4];
	ListPtr	list;

	if (!boss || !print_func) return;

	(*print_func)(dest, " TYP AN AT   CHARGE     SIGMA    EPSILON\n");

if (!(boss->flags & B_QM)) {
	ForEachList(boss->atomtypes,list) {
		if (list->L_TYPE <= 0) continue;
		if (!((list->L_TYPE >= 800 && list->L_TYPE < 900) || (list->L_TYPE >= 9500))) continue;

		if (strlen(list->L_AMBER) < 3) sprintf(amber, "%-2s ", list->L_AMBER);
		else strncpy(amber, list->L_AMBER, 3);
		amber[3] = '\0';
		if (strncmp(amber, "  ", 2) == 0 || strncmp(amber, "??", 2) == 0) {
			sprintf(amber, "%2s ", GetElemSymbol(list->L_AN));
		}

		sprintf(str, "%4d %02d %3s% 10.6f %9.6f %9.6f     %s\n",
			list->L_TYPE, list->L_AN, amber, list->L_CHARGE,
			list->L_SIGMA, list->L_EPSILON, list->L_COMMENT);
		(*print_func)(dest, str);
		flags = 1;
	}

	if (!flags) {
		sprintf(str, "%4d %02d %3s% 10.6f %9.6f %9.6f     %s\n", 800, 0, "DUM", 0.0, 0.0, 0.0, "DUMMY");
		(*print_func)(dest, str);
	}
	ForEachList(boss->atomtypes,list) {
		if (list->L_TYPE <= 0) continue;
		if (((list->L_TYPE >= 800 && list->L_TYPE < 900) || (list->L_TYPE >= 9500))) continue;

		if (strlen(list->L_AMBER) < 3) sprintf(amber, "%-2s ", list->L_AMBER);
		else strncpy(amber, list->L_AMBER, 3);
		amber[3] = '\0';
		if (strncmp(amber, "  ", 2) == 0 || strncmp(amber, "??", 2) == 0) {
			sprintf(amber, "%2s ", GetElemSymbol(list->L_AN));
		}

		sprintf(str, "%4d %02d %3s% 10.6f %9.6f %9.6f     %s\n",
			list->L_TYPE, list->L_AN, amber, list->L_CHARGE,
			list->L_SIGMA, list->L_EPSILON, list->L_COMMENT);
		(*print_func)(dest, str);
	}
}

		(*print_func)(dest, "      Torsion Types\n");
		(*print_func)(dest, "      The following are the Fourier coefficients for dihedral angles.\n");
		(*print_func)(dest, "Type   V1         V2    f   V3    f   V4    f     Symbols   (  OPLS Types   )\n");

	if (boss->torsiontypes && NEED_PRINT_BOSS_ZMAT_TORSION(boss)) {
		ForEachList(boss->torsiontypes,list) {
			sprintf(str, "%03d  % 8.4f  % 8.4f %d% 8.4f %d% 8.4f %d   %s\n",
				list->L_TYPE, list->L_V0, list->L_V1, list->L_N1, list->L_V2, list->L_N2,
				list->L_V3, list->L_N3, list->L_COMMENT);
			(*print_func)(dest, str);
		}
	}
}

static void	GPrintBossZmatWithTypes (BossZmatPtr boss, void *dest, void (*print_func)(void *, char*))
{
	if (!boss || !print_func) return;
	if (boss->atomtypes || boss->torsiontypes) {
		GPrintBossZmat(boss, dest, print_func);
		(*print_func)(dest, "\nFinal Non-Bonded Parameters for QM Atoms:\n");
		GPrintBossZmatTypes(boss, dest, print_func);
	} else {
		GPrintBossZmat(boss, dest, print_func);
	}
}

void	SPrintBossZmatWithTypes (char **str_sum, BossZmatPtr boss)
{
	if (!str_sum || !boss) return;
	BgnPrintString();
	GPrintBossZmatWithTypes(boss, (void *)NULL, SPrintString);
	*str_sum = EndPrintString();
}

void	FPrintBossZmatWithTypes (FILE *fp, BossZmatPtr boss)
{
	if (!fp || !boss) return;
	GPrintBossZmatWithTypes(boss, (void *)fp, FPrintString);
}

static int	FReadBossParamFile (BossZmatPtr boss, FILE *fp, FILE *parfp, int flags)
{
	int	n;
	ListPtr	list, nextlist, newlist;
	ListPtr	atomtypes=NULL, torsiontypes=NULL;

	clearmolerror();
	if (parfp) {
		if (!(boss->param = FLoadBossPar(parfp, 0, flags))) {
			if (getmolerrorno() == -1) setmolerrorno(MERR_PARFORMAT);
			printmolerror();
		}
		if ((parfp != stdin) && !(flags & B_NO_REWIND_PARFILE)) rewind(parfp);
		if (!(boss->atomtypes = FLoadBossAtomTypes(parfp))) {
			if (getmolerrorno() == -1) setmolerrorno(MERR_PARFORMAT);
			printmolerror();
		}
		if (!(boss->torsiontypes = FLoadBossTorsion(parfp))) {
			if (getmolerrorno() == -1) setmolerrorno(MERR_PARFORMAT);
			printmolerror();
		}
	}

	if (fp && (flags & B_PARAM_IN_ZMAT)) {
		atomtypes = FLoadBossAtomTypesFromZmat(fp);
PRINT_ERROR(1);
		torsiontypes = FLoadBossTorsion(fp);
PRINT_ERROR(1);

		/* overrides atomtypes from the parameter file */
		if (atomtypes) {
			for(list=boss->atomtypes;list;list=nextlist) {
				nextlist = list->next;
				for(newlist=atomtypes;newlist;newlist=newlist->next) {
					if (newlist->L_TYPE == list->L_TYPE) {
						break;
					}
				}
				if (!newlist) continue;
				RemoveList(list, &boss->atomtypes);
			}

if (DEBUGLEVEL(1)) {
			fprintf(stderr, "Following atom types were read from Z-matrix.\n");
			n = 0;
			ForEachList(atomtypes,list) {
				n++;
				if (n > 1 && (n-1)%20==0) fprintf(stderr, "\n");
				fprintf(stderr, "%3d ", list->L_TYPE);
			}
			fprintf(stderr, "\n");
}

			EnterList(atomtypes, &boss->atomtypes);
		}

		/* overrides torsion types from the parameter file */
		if (torsiontypes) {
			for(list=boss->torsiontypes;list;list=nextlist) {
				nextlist = list->next;
				for(newlist=torsiontypes;newlist;newlist=newlist->next) {
					if (newlist->L_TYPE == list->L_TYPE) {
						break;
					}
				}
				if (!newlist) continue;
				RemoveList(list, &boss->torsiontypes);
			}

if (DEBUGLEVEL(1)) {
			fprintf(stderr, "Following torsion types were read from Z-matrix.\n");
			n = 0;
			ForEachList(torsiontypes,list) {
				n++;
				if (n > 1 && (n-1)%20==0) fprintf(stderr, "\n");
				fprintf(stderr, "%3d ", list->L_TYPE);
			}
			fprintf(stderr, "\n");
}

			EnterList(torsiontypes, &boss->torsiontypes);
		}

		clearmolerror();	/* we don't care if error has occured here */
	}

	return getmolerrorno() == -1;
}

void	UpdateBossZmatTypes (BossZmatPtr boss, OplsTypeTorsion *opls)
{
	ZmatPtr	zmat, zmatlist;
	ListPtr	list, newlist, dihed;
	
	ListPtr	atomtypes=NULL, torsiontypes=NULL;

	if (!boss || !opls) return;

	zmatlist = LinkZmatList(boss->zmatlist);
	ForEachList(opls->atomtypelist,list) {
		ForEachZmat(zmatlist,zmat) {
			if (list->L_TYPE == zmat->itype || list->L_TYPE == zmat->ftype) break;
		}
		if (!zmat) continue;
		newlist = NewList();
		BCOPY(list, newlist, sizeof(List));
		newlist->next = NULL;
		EnterList(newlist, &atomtypes);
	}
	ForEachList(opls->torsionlist,list) {
		ForEachList(boss->vardihed,dihed) {
			if (dihed->L_ITYPE == list->L_TYPE || dihed->L_FTYPE == list->L_TYPE) break;
		}
		if (!dihed) continue;
		newlist = NewList();
		BCOPY(list, newlist, sizeof(List));
		newlist->next = NULL;
		EnterList(newlist, &torsiontypes);
	}
	BreakZmatList(boss->zmatlist);
	
	if (atomtypes) {
		FreeList(&boss->atomtypes);
		boss->atomtypes = atomtypes;
	}
	if (torsiontypes) {
		FreeList(&boss->torsiontypes);
		boss->torsiontypes = torsiontypes;
	}
}

void	ScreenBossZmatTypes (BossZmatPtr boss)
{
	int	n, type;
	ListPtr	list, nextlist, dihed;
	ZmatPtr	zmat, zmatlist;

	if (!boss) return;
	zmatlist = LinkZmatList(boss->zmatlist);

	/* remove unnecessary atom and torsion types */
	for(list=boss->atomtypes;list;list=nextlist) {
		nextlist = list->next;
		type = list->L_TYPE;
		if (type <= 0) continue;
		ForEachZmat(zmatlist,zmat) if (zmat->itype == type || zmat->ftype == type) break;
		if (zmat) continue;

		/* this atomtype is not used in the Z-matrix */
		DeleteList(list, &boss->atomtypes);
		FreeListEntry(&list);
	}

	/* sort atom type list by type (L_TYPE) in ascending order */
	n = 0;
	ForEachList(boss->atomtypes,list) n++;
	if (n > 1) QsortList(boss->atomtypes, 1, n, M_DATA_1);

	for(list=boss->torsiontypes;list;list=nextlist) {
		nextlist = list->next;
		type = list->L_TYPE;
		if (type <= 0) continue;
		ForEachList(boss->vardihed,dihed)
			if (dihed->L_ITYPE == type || dihed->L_FTYPE == type) break;
		if (dihed) continue;
		ForEachList(boss->adddihed,dihed)
			if (dihed->L_ITYPE == type || dihed->L_FTYPE == type) break;
		if (dihed) continue;

		/* this torsion type is not used in the Z-matrix */
		DeleteList(list, &boss->torsiontypes);
		FreeListEntry(&list);
	}

	/* sort torsion type list by type (L_TYPE) in ascending order */
	n = 0;
	ForEachList(boss->torsiontypes,list) n++;
	if (n > 1) QsortList(boss->torsiontypes, 1, n, M_DATA_1);

	BreakZmatList(boss->zmatlist);
}

BossZmatPtr	FLoadBossZmat (FILE *fp, FILE *parfp, int flags)
{
	BossZmatPtr	boss;
	ZmatPtr	zmat, zmatlist;
	ListPtr	list;

	if (!fp) return NULL;
	if (!(boss = FLoadBossZmatInput(fp, flags))) return NULL;
	FReadBossParamFile(boss, fp, parfp, flags|(boss->flags & B_PARAM_IN_ZMAT));
	boss->flags &= ~B_PARAM_IN_ZMAT;
	AssignBossPar(boss->zmatlist, boss->atomtypes);

	if ((flags & B_QM) && !boss->atomtypes) {
		zmatlist = LinkZmatList(boss->zmatlist);
		ForEachZmat(zmatlist,zmat) zmat->an = zmat->itype;
		BreakZmatList(boss->zmatlist);
	}

	if ((flags & B_SOLVENT_ZMAT) && boss->atomtypes) {
		zmatlist = LinkZmatList(boss->zmatlist);
		ForEachZmat(zmatlist,zmat) {
			if (!(list = FindBossAtomType(boss->atomtypes, zmat->itype))) continue;
			zmat->charge = list->L_CHARGE;
			zmat->sigma = list->L_SIGMA;
			zmat->epsilon = list->L_EPSILON;
		}
		BreakZmatList(boss->zmatlist);
	}

	/* get rid of unused atom and torsion types from boss->atomtypes
	 * and boss->torsiontypes
	 */
	ScreenBossZmatTypes(boss);

	return boss;	/* this should be the only 'return' statement */
}

MolPtr	CreateMolFromBossZmat (BossZmatPtr boss)
{
	MolPtr	mol;
	AtomPtr	atom, conn_atom=NULL;
	BondPtr	bondlist=NULL;
	ResiduePtr	res=NULL, rinit, rterm;
	ChainPtr	chain=NULL, chainlist=NULL;
	ZmatPtr	zmat;
	ListPtr	l, atomtype_list;
	int	serno, refno, n, hetatm;
	int	chain_id = 'A'-1;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	if (!boss || !boss->zmatlist) return NULL;

	strcpy(patom->name, " ");
	chainlist = NULL;
	serno = 0;
	ForEachList(boss->zmatlist,l) {
		chain = NULL;
		res = NULL;

		chain_id++;
		strcpy(pres->name, "UNK");
		pres->chain_id = chain_id;
		pres->seq = 0;

		ForEachZmat((ZmatPtr)l->L_ZMAT,zmat) {
			serno++;
			patom->serial = serno;
			strcpy(patom->name, AN2ATOMNAME(zmat->an));

			patom->x = zmat->cart[0];
			patom->y = zmat->cart[1];
			patom->z = zmat->cart[2];
			if (zmat->resname[0]) {
				strcpy(pres->name, zmat->resname);
				TrimTail(pres->name, " \t\n");
				pres->seq = zmat->resseq;
				refno = GetResRefno(pres->name);
				hetatm = !IsAminoNucleo(refno);
			} else hetatm = 1;
/*
hetatm = 1;
*/

			atom = CreatePdbAtom(&chainlist, &chain, &res, &conn_atom, patom, &bondlist, hetatm);
zmat->atom = (int)atom;
			atom->an = zmat->an;
			atom->type = zmat->itype;
			atom->refno = GetAtomRefno(zmat->label[0] ? zmat->label : AN2ATOMNAME(atom->an));
			if ((atomtype_list = FindBossAtomType(boss->atomtypes, atom->type))) {
				strcpy(atom->ambertype, atomtype_list->L_AMBER);
				atom->charge = atomtype_list->L_CHARGE;
				atom->sigma = atomtype_list->L_SIGMA;
				atom->epsilon = atomtype_list->L_EPSILON;
			}

if (boss->flags & B_USER_CHARGE) atom->charge = zmat->charge;
if (boss->flags & B_USER_SIGMA) atom->sigma = zmat->sigma;
if (boss->flags & B_USER_EPSILON) atom->epsilon = zmat->epsilon;

		}
	}
	if (!chainlist) return NULL;

	/* update zmat->atom */
	ForEachList(boss->zmatlist,l) {
		ForEachZmat((ZmatPtr)l->L_ZMAT,zmat) {
			zmat->atom = SearchAtomInChainByPointer(chainlist, (AtomPtr)zmat->atom);
		}
	}

	/* create a molecule */
	mol = NewMol();
	mol->header = LookUpStrTable(boss->title);
	mol->chain = chainlist;
	if (bondlist) mol->bond = bondlist;
	mol->param_type = PARAM_OPLS;

/*
if (mol->bond) {
ForEachBond(mol->bond,bondlist) fprintf(stderr, "%d-%d\n", bondlist->atom1->serno, bondlist->atom2->serno);
}
*/

	GetNeighbor(mol);
	CalcBondsWithinResidues(mol);
	CalcBondsBetweenResidues(mol);
	FixConnectionTable(mol);

	/* SS list */
	ForEachChain(mol->chain,chain) {
		n = 0;
		ForEachRes(chain->residue,res) {
			n++;
			res->ss_flags = 0;
			if (!res->next) rterm = res;
		}
		if (n < 3) continue;
		rinit = chain->residue;
		if (rinit) rinit->ss_flags = R_SSINIT;
		if (rterm) rterm->ss_flags = R_SSTERM;
	}
	mol->ss = MakeSSList(mol->chain);

	return mol;
}

MolPtr	FLoadBossPdbInput (FILE *fp)
{
	clearmolerror();
	if (!fp) return NULL;
	return FLoadPdbFile(fp);
}

MolPtr	FLoadBossMindtoolInput (FILE *fp)
{
	clearmolerror();
	if (!fp) return NULL;
	return FLoadMindtool(fp, 0);
}

static void	create_geom_var (ZmatPtr ini_zmat, ZmatPtr fin_zmat, int param, ListPtr *geomvar)
{
	ListPtr	l;
	ZmatPtr	iz, fz;
	int	n;
	double	del;

	if (param < 1 || param > 3) return;
	param--;
	for(iz=ini_zmat,fz=fin_zmat,n=1;iz && fz;iz=iz->next,fz=fz->next,n++) {
		del = ABS(fz->zval[param]-iz->zval[param]);
		if (del > 1.0E-3) {
			l = EnterNewList(geomvar);
			l->L_CENTER = n;
			l->L_PARAM = param+1;
			l->L_VAL = fz->zval[param];
		}
	}
}

static int	compare_zmat_with_mol (ZmatPtr zmat, MolPtr mol)
{
	int	nzmats, natoms;
	AtomPtr	a;
	ZmatPtr	z;

	if (!zmat || !mol) return 0;
	nzmats = 0;
	ForEachZmat(zmat,z) nzmats++;

	natoms = 0;
	a = FirstAtom(mol->chain); natoms++;
	while ((a = NextAtom(a))) natoms++;

	if (nzmats != natoms) {
		return setmolerror("Number of atoms (%d) do not match the Z-Matrix (%d).", natoms, nzmats);
	}

	return 	1;
}

int	UpdateBossGeomVarWithMolCoord (BossZmatPtr boss, MolPtr mol)
{
	ZmatPtr	zmat, fin_zmat;

	if (!boss || !boss->zmatlist || !mol || !mol->chain) return 0;
	if (!(zmat = LinkZmatList(boss->zmatlist))) return 0;
	if (!compare_zmat_with_mol(zmat, mol)) {
		BreakZmatList(boss->zmatlist);
		return 0;
	}

	fin_zmat = CopyZmat(zmat);
	UpdateZmatWithMolCoord(fin_zmat, mol);

	FreeList(&boss->geomvar);
	create_geom_var(zmat, fin_zmat, 1, &boss->geomvar);
	create_geom_var(zmat, fin_zmat, 2, &boss->geomvar);
	create_geom_var(zmat, fin_zmat, 3, &boss->geomvar);

	BreakZmatList(boss->zmatlist);
	FreeZmat(&fin_zmat);
	return 1;
}

ZmatPtr	GetZmatFromBossGeomVar (BossZmatPtr boss)
{
	ZmatPtr	z, zmat;
	List	list;
	ListPtr	l;
	int	n, p;

	if (!boss) return NULL;
	if (!(zmat = LinkZmatList(boss->zmatlist))) return NULL;
	zmat = CopyZmat(zmat);
	BreakZmatList(boss->zmatlist);
	if (!zmat) return NULL;

	ForEachList(boss->geomvar,l) {
		n = l->L_CENTER;
		p = l->L_PARAM;
		if (p < 1 || p > 3) continue;
		if (!(z = SearchZmat(zmat,n))) continue;
		z->zval[p-1] = l->L_VAL;
	}

	/* swap initial and final types */
	ForEachZmat(zmat,z) {
		n = z->itype;
		z->itype = z->ftype;
		z->ftype = n;
		if (z->itype == 0) z->itype = z->ftype;
	}

	list.next = NULL;
	list.L_ZMAT = (size_t)zmat;
	AssignBossPar(&list, boss->atomtypes);

	return zmat;
}

ListPtr	GetZmatListFromBossGeomVar (BossZmatPtr boss)
{
	ZmatPtr	z, zmat;
	ListPtr	l, zmatlist;
	int	n, p;

	if (!boss || !(zmatlist = CopyZmatList(boss->zmatlist))) return NULL;
	if (!(zmat = LinkZmatList(zmatlist))) {
		FreeZmatList(&zmatlist);
		return NULL;
	}
	ForEachList(boss->geomvar,l) {
		n = l->L_CENTER;
		p = l->L_PARAM;
		if (p < 1 || p > 3) continue;
		if (!(z = SearchZmat(zmat,n))) continue;
		z->zval[p-1] = l->L_VAL;
	}
	/* swap initial and final types */
	ForEachZmat(zmat,z) {
		n = z->itype;
		z->itype = z->ftype;
		z->ftype = n;
	}

	BreakZmatList(zmatlist);
	AssignBossPar(zmatlist, boss->atomtypes);

	return zmatlist;
}

/* apply boss->geomvar (Geometry Variations) to the boss z-matrix
   (fill in zmat->zval), evaluate Cartesian coords (fill in zmat->cart),
   and update the atomic coordinates of moleclue (fill atom->x,atom->y,atom->z)
*/

int	EvalBossGeomVar (BossZmatPtr boss, MolPtr mol)
{
	ZmatPtr	zmat;
	int	ret_val;

	if (!boss || !boss->zmatlist || !mol || !mol->chain) return 0;
	if (!(zmat = LinkZmatList(boss->zmatlist))) return 0;
	ret_val = compare_zmat_with_mol(zmat, mol);
	BreakZmatList(boss->zmatlist);
	if (ret_val == 0) return 0;

	if (!(zmat = GetZmatFromBossGeomVar(boss))) return 0;

	/* update molecule coordinates */
	UpdateMolWithZmat(zmat, mol, Z_UPDATE_ATOM_ATTR);
	FreeZmat(&zmat);

	return 1;
}

GaussZmatPtr	CreateGaussZmatFromBoss (BossZmatPtr boss)
{
	int	n;
	ZmatPtr	zmat, zmatlist;
	ListPtr	list;
	GaussZmatPtr    g=NULL;

	if (!boss || !boss->zmatlist) return NULL;
	if (!(zmatlist = LinkZmatList(boss->zmatlist))) return NULL;

	if ((g = NewGaussZmat())) {
		g->zmatlist = CopyZmat(zmatlist);
		if (boss->title[0]) strcpy(g->title, boss->title);
		n = 0;
		ForEachZmat(g->zmatlist,zmat) {
			n++;
			ForEachList(boss->varbond, list) {
				if (list->L_CENTER == n) {
					zmat->flags |= Z_VAR_LENGTH;
					break;
				}
			}
			ForEachList(boss->varangle,list) {
				if (list->L_CENTER == n) {
					zmat->flags |= Z_VAR_ANGLE;
					break;
				}
			}
			ForEachList(boss->vardihed,list) {
				if (list->L_CENTER == n) {
					zmat->flags |= Z_VAR_DIHED;
					break;
				}
			}

			if (zmat->an < 0 || zmat->an >= 98) strcpy(zmat->label, "X");
			else strcpy(zmat->label, GetElemSymbol(zmat->an));
			strcpy(zmat->atomname[0], zmat->label);

			if (n > 1) sprintf(zmat->varname[0], "R%d", n);
			if (n > 2) sprintf(zmat->varname[1], "A%d", n);
			if (n > 3) sprintf(zmat->varname[2], "D%d", n);
		}
	}

	BreakZmatList(boss->zmatlist);
	return g;
}

JaguarZmatPtr	CreateJaguarZmatFromBoss (BossZmatPtr b)
{
	GaussZmatPtr	g;
	JaguarZmatPtr	j;

	if (!b || !(g = CreateGaussZmatFromBoss(b))) return NULL;
	j = CreateJaguarZmatFromGauss(g);
	FreeGaussZmat(&g);
	return j;
}

/**************************************************************************
 *******   Add All Variable/Aditional Bonds, Angles, or Dihedrals   *******
 **************************************************************************/

/* variable for addall options */
static BossZmatAddAllOpt	addall_opt = {
	1,	/* initial dihedral type */
	1,	/* final dihedral type */
	5.0,	/* variable dihedral range */
	0.0,	/* Min bond length for conformational search */
	5.0,	/* Max bond length */
	0.0,	/* Min bond angle */
	180.0,	/* Max bond angle */
	0.0,	/* Min dihedral */
	360.0	/* Max dihedral */
	};

int	AddAllBossZmatOptions (BossZmatPtr boss, MolPtr mol, unsigned long flags, BossZmatAddAllOpt *addall)
{
	int	i, first_flag=1, ndummy=0;
	ZmatPtr	zmat, z, z1, z2, z3;
	ListPtr	list, next, newlist;
	BossZmatAddAllOpt	*p;

	if (!boss) return 0;
	if (boss->flags & B_SOLVENT_ZMAT) {
		if (boss->flags & B_AUTO_ADDBOND) {
			FreeList(&boss->addbond);
			boss->flags &= ~B_AUTO_ADDBOND;
		}
		if (boss->flags & B_AUTO_ADDANGLE) {
			FreeList(&boss->addangle);
			boss->flags &= ~B_AUTO_ADDANGLE;
		}
		if (boss->flags & B_AUTO_ADDDIHED) {
			FreeList(&boss->adddihed);
			boss->flags &= ~B_AUTO_ADDDIHED;
		}
		boss->flags &= ~B_VARDIHED_EXPLICIT;
	}

	if (flags & M_BOSS_ZMAT_AUTO_ADD_BOND) {
		FreeList(&boss->addbond);
		boss->flags |= B_AUTO_ADDBOND;
		flags &= ~M_BOSS_ZMAT_ADD_BOND;
	}

	if (flags & M_BOSS_ZMAT_AUTO_ADD_ANGLE) {
		FreeList(&boss->addangle);
		boss->flags |= B_AUTO_ADDANGLE;
		flags &= ~M_BOSS_ZMAT_ADD_ANGLE;
	}

	if (flags & M_BOSS_ZMAT_AUTO_ADD_DIHED) {
		FreeList(&boss->adddihed);
		boss->flags |= B_AUTO_ADDDIHED;
		flags &= ~M_BOSS_ZMAT_ADD_DIHED;
	}

	p = addall;
	if (!p) p = &addall_opt;
	zmat = LinkZmatList(boss->zmatlist);

	if (flags & M_BOSS_ZMAT_VAR) {
		if ((flags & M_BOSS_ZMAT_VAR_DIHED) && (flags & M_BOSS_ZMAT_VAR_DIHED_EXPLICIT)) {
			boss->flags |= B_VARDIHED_EXPLICIT;
		}

		if ((flags & M_BOSS_ZMAT_VAR_BOND)  && boss->varbond)  FreeList(&boss->varbond);
		if ((flags & M_BOSS_ZMAT_VAR_ANGLE) && boss->varangle) FreeList(&boss->varangle);
		if ((flags & M_BOSS_ZMAT_VAR_DIHED) && boss->vardihed) FreeList(&boss->vardihed);
		first_flag = 1;
		ndummy = 0;
		i = 0;
		ForEachZmat(zmat,z) {
			i++;
			if (z->an < 0) {
				if (first_flag) ndummy++;
				continue;
			} else {
				first_flag = 0;
			}

			if (i <= ndummy + 1) continue;

			z1 = SearchZmat(zmat,z->zdef[0]);
			z2 = SearchZmat(zmat,z->zdef[1]);
			z3 = SearchZmat(zmat,z->zdef[2]);

			if (flags & M_BOSS_ZMAT_VAR_BOND) {
				/* add to boss->varbond */
				list = EnterNewList(&boss->varbond);
				list->L_CENTER = i;
			}
			if (i <= ndummy + 2) continue;
			if (flags & M_BOSS_ZMAT_VAR_ANGLE) {
				/* add to boss->varangle */
				list = EnterNewList(&boss->varangle);
				list->L_CENTER = i;
			}
			if (i <= ndummy + 3) continue;
			if (flags & M_BOSS_ZMAT_VAR_DIHED) {
				/* add to boss->vardihed */
				list = EnterNewList(&boss->vardihed);
				list->L_CENTER = i;
				list->L_ITYPE = p->var_dihed_itype;
				list->L_FTYPE = p->var_dihed_ftype;
				list->L_RANGE = p->var_dihed_range;
			}
		}
	}

	if (!mol) {
		BreakZmatList(boss->zmatlist);
		return 1;
	}

	if (flags & M_BOSS_ZMAT_ADD_BOND) {
		FreeList(&boss->addbond);
		list = GetAddBondList(mol->chain, zmat);
		if (list) EnterList(list, &boss->addbond);
	}
	if (flags & M_BOSS_ZMAT_ADD_ANGLE) {
		FreeList(&boss->addangle);
		list = GetAddAngleList(mol->chain, zmat);
		if (list) EnterList(list, &boss->addangle);
	}
	if (flags & M_BOSS_ZMAT_ADD_DIHED) {
		FreeList(&boss->adddihed);
		list = GetAddDihedList(mol->chain, zmat);
		if (list) EnterList(list, &boss->adddihed);
		ForEachList(boss->adddihed,list) {
			list->L_ITYPE = p->var_dihed_itype;
			list->L_FTYPE = p->var_dihed_ftype;
		}
	}
	if (flags & M_BOSS_ZMAT_CSEARCH) {
		int	tmp_flags = 0;

		if (!boss->varbond)  {
			tmp_flags |= M_BOSS_ZMAT_VAR_BOND;
			AddAllBossZmatOptions(boss, mol, M_BOSS_ZMAT_VAR_BOND,  p);
		}
		if (!boss->varangle) {
			tmp_flags |= M_BOSS_ZMAT_VAR_ANGLE;
			AddAllBossZmatOptions(boss, mol, M_BOSS_ZMAT_VAR_ANGLE, p);
		}
		if (!boss->vardihed) {
			tmp_flags |= M_BOSS_ZMAT_VAR_DIHED;
			AddAllBossZmatOptions(boss, mol, M_BOSS_ZMAT_VAR_DIHED, p);
		}

		FreeList(&boss->conform);
		ForEachList(boss->varbond,list) {
			next = list->next;
			list->next = NULL;
			newlist = CopyList(list);
			list->next = next;

			EnterList(newlist, &boss->conform);
			newlist->L_PARAM = 1;
			newlist->L_RMIN = p->var_bond_rmin;
			newlist->L_RMAX = p->var_bond_rmax;
		}
		ForEachList(boss->varangle,list) {
			next = list->next;
			list->next = NULL;
			newlist = CopyList(list);
			list->next = next;

			EnterList(newlist, &boss->conform);
			newlist->L_PARAM = 2;
			newlist->L_RMIN = p->var_angle_rmin;
			newlist->L_RMAX = p->var_angle_rmax;
		}
		ForEachList(boss->vardihed,list) {
			next = list->next;
			list->next = NULL;
			newlist = CopyList(list);
			list->next = next;

			EnterList(newlist, &boss->conform);
			newlist->L_PARAM = 3;
			newlist->L_RMIN = p->var_dihed_rmin;
			newlist->L_RMAX = p->var_dihed_rmax;
		}

		if (tmp_flags & M_BOSS_ZMAT_VAR_BOND)  FreeList(&boss->varbond);
		if (tmp_flags & M_BOSS_ZMAT_VAR_ANGLE) FreeList(&boss->varangle);
		if (tmp_flags & M_BOSS_ZMAT_VAR_DIHED) FreeList(&boss->vardihed);
	}

	BreakZmatList(boss->zmatlist);

	return 1;
}

BossZmatAddAllOpt	*GetAddAllBossZmatOptions (void)
{
	return &addall_opt;
}

void	SetImproperTorsionFlags (ListPtr dihedlist, ChainPtr chainlist, BondPtr bondlist)
{
	ListPtr	list;
	int	a1, a2, a3, a4;

	ForEachList(dihedlist,list) {
		list->L_IMPROPER = 0;
		a1 = list->L_ATOM1;
		a2 = list->L_ATOM2;
		a3 = list->L_ATOM3;
		a4 = list->L_ATOM4;

		/* A torsion is proper, if
		   a1-a2, a2-a3, and a3-a4 are bonded and
		   a1-a3 and a2-a4 are not bonded.
		*/
		if (GetSernoBond(bondlist,a1,a2) &&  GetSernoBond(bondlist,a2,a3)
			&& GetSernoBond(bondlist,a3,a4) &&
		   !GetSernoBond(bondlist,a1,a3) && !GetSernoBond(bondlist,a2,a4)) {
			list->L_IMPROPER = 0;
		} else list->L_IMPROPER = 1;
	}
}

