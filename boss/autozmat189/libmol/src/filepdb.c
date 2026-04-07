#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

/*
 read a line from file pointed by 'fp' and return PdbRecord structure.
*/
PdbRecord	readpdbrecord (FILE *fp)
{
	char	*fmt, str[82];
	int	ip;
	char	*cp;
	/*
	static PdbRecord	r, *rp = &r;
	*/
	PdbRecord	r, *rp = &r;

	BZERO((void *)rp, sizeof(PdbRecord));
	BZERO((void *)str, sizeof(str));
	if (!fgets(str, sizeof(str), fp)) {
		rp->type = PDB_END;
		return r;
	}

	if ((cp = strchr(str, '\n'))) *cp = '\0';		/* mark the end */
	else {	/* line longer than 80 chars. read in the rest and discard */
		while ((ip = getc(fp)) != '\n' && ip != EOF);
	}

	rp->type = getpdbrecordtype(str);
	fmt      = getpdbrecordscanfmt(rp->type);
	scanpdbrecord(rp, str, fmt);

/*
fprintf(stderr, "I=%s\n", str);
fprintf(stderr, "O=");
fprintpdbrecord(rp, stderr);
*/

	return r;
}

int	getpdbrecordtype (char *str)
{
	int	type, i;
	char	keywd[5];

	if (!str) return PDB_UNK;
	BZERO(keywd, sizeof(keywd));
	sscanf(str, "%4c", keywd);
	for(i=0;i<4;i++) if (!keywd[i]) keywd[i] = ' ';
	keywd[4] = '\0';

	type = PDB_UNK;
	keywd[0] = tolower(keywd[0]);
	switch (keywd[0]) {
	case 'a':
		if      (strncasecmp(keywd, "atom", 4) == 0) type = PDB_ATOM;
		else if (strncasecmp(keywd, "auth", 4) == 0) type = PDB_AUTHOR;
		else if (strncasecmp(keywd, "anis", 4) == 0) type = PDB_ANISOU;
		break;

	case 'c':
		if      (strncasecmp(keywd, "comp", 4) == 0) type = PDB_COMPND;
		else if (strncasecmp(keywd, "crys", 4) == 0) type = PDB_CRYST1;
		else if (strncasecmp(keywd, "cone", 4) == 0) type = PDB_CONECT;
		break;

	case 'e':
		if      (strncasecmp(keywd, "end ", 4) == 0) type = PDB_END;
		else if (strncasecmp(keywd, "endm", 4) == 0) type = PDB_ENDMDL;
		break;

	case 'f':
		if      (strncasecmp(keywd, "ftno", 4) == 0) type = PDB_FTNOTE;
		else if (strncasecmp(keywd, "form", 4) == 0) type = PDB_FORMUL;
		break;

	case 'h':
		if      (strncasecmp(keywd, "heta", 4) == 0) type = PDB_HETATM;
		else if (strncasecmp(keywd, "head", 4) == 0) type = PDB_HEADER;
		else if (strncasecmp(keywd, "het ", 4) == 0) type = PDB_HET;
		else if (strncasecmp(keywd, "heli", 4) == 0) type = PDB_HELIX;
		break;

	case 'j':
		if      (strncasecmp(keywd, "jrnl", 4) == 0) type = PDB_JRNL;
		break;

	case 'm':
		if      (strncasecmp(keywd, "mtri", 4) == 0) type = PDB_MTRIX;
		else if (strncasecmp(keywd, "mast", 4) == 0) type = PDB_MASTER;
		else if (strncasecmp(keywd, "mode", 4) == 0) type = PDB_MODEL;
		break;

	case 'o':
		if      (strncasecmp(keywd, "obsl", 4) == 0) type = PDB_OBSLTE;
		else if (strncasecmp(keywd, "orig", 4) == 0) type = PDB_ORIGX;
		break;

	case 'r':
		if      (strncasecmp(keywd, "rema", 4) == 0) type = PDB_REMARK;
		else if (strncasecmp(keywd, "revd", 4) == 0) type = PDB_REVDAT;
		break;

	case 's':
		if      (strncasecmp(keywd, "scal", 4) == 0) type = PDB_SCALE;
		else if (strncasecmp(keywd, "seqr", 4) == 0) type = PDB_SEQRES;
		else if (strncasecmp(keywd, "shee", 4) == 0) type = PDB_SHEET;
		else if (strncasecmp(keywd, "site", 4) == 0) type = PDB_SITE;
		else if (strncasecmp(keywd, "siga", 4) == 0) type = PDB_SIGATM;
		else if (strncasecmp(keywd, "sigu", 4) == 0) type = PDB_SIGUIJ;
		else if (strncasecmp(keywd, "sour", 4) == 0) type = PDB_SOURCE;
		else if (strncasecmp(keywd, "sprs", 4) == 0) type = PDB_SPRSDE;
		else if (strncasecmp(keywd, "ssbo", 4) == 0) type = PDB_SSBOND;
		break;

	case 't':
		if      (strncasecmp(keywd, "turn", 4) == 0) type = PDB_TURN;
		else if (strncasecmp(keywd, "tvec", 4) == 0) type = PDB_TVECT;
		else if (strncasecmp(keywd, "ter ", 4) == 0) type = PDB_TER;
		break;
	}

	return type;
}

int	scanpdbrecord (PdbRecord *record, char *str, char *fmt)
{
	PdbRecPtr	r;
	PdbSheet	*sh;

	if (!record || !fmt) return 0;
	r = &record->rec;

	switch (record->type) {
	case PDB_UNK:
		strscan(str, fmt, r->unk.text);
		break;

	case PDB_ANISOU:
	case PDB_SIGUIJ:
		strscan(str, fmt, &r->anisou.serial, r->anisou.name, &r->anisou.altloc, r->anisou.residue.name,
			&r->anisou.residue.chain_id, &r->anisou.residue.seq, &r->anisou.residue.insert,
			&r->anisou.u[0], &r->anisou.u[1], &r->anisou.u[2],
			&r->anisou.u[3], &r->anisou.u[4], &r->anisou.u[5]);
		break;

	case PDB_ATOM:
	case PDB_HETATM:
	case PDB_SIGATM:
		strscan(str, fmt, &r->atom.serial, r->atom.name, &r->atom.altloc,
			r->atom.residue.name, &r->atom.residue.chain_id,
			&r->atom.residue.seq, &r->atom.residue.insert,
			&r->atom.x, &r->atom.y, &r->atom.z,
			&r->atom.occupancy, &r->atom.tempfac, &r->atom.footnote);
		break;

	case PDB_AUTHOR:
	case PDB_COMPND:
	case PDB_JRNL:
	case PDB_SOURCE:
		strscan(str, fmt, &r->author.cont, r->author.text);
		break;

	case PDB_CONECT:
		strscan(str, fmt, &r->conect.serial,
			&r->conect.covalent[0], &r->conect.covalent[1],
			&r->conect.covalent[2], &r->conect.covalent[3],
			&r->conect.donor[0], &r->conect.donor[1], &r->conect.salt_minus,
			&r->conect.acceptor[0], &r->conect.acceptor[1], &r->conect.salt_plus);
		break;

	case PDB_CRYST1:
		strscan(str, fmt, &r->cryst1.a, &r->cryst1.b, &r->cryst1.c,
			&r->cryst1.alpha, &r->cryst1.beta, &r->cryst1.gamma,
			r->cryst1.spacegrp, &r->cryst1.z);
		break;

	case PDB_END:
		break;

	case PDB_ENDMDL:
		break;

	case PDB_FORMUL:
		strscan(str, fmt, &r->formul.component, r->formul.hetid, &r->formul.cont,
		&r->formul.exclude, r->formul.formula);
		break;

	case PDB_FTNOTE:
	case PDB_REMARK:
		strscan(str, fmt, &r->ftnote.num, r->ftnote.text);
		break;

	case PDB_HEADER:
		strscan(str, fmt, r->header.funcclass, r->header.date,
			r->header.id);
		break;

	case PDB_HELIX:
		strscan(str, fmt, &r->helix.serial, r->helix.id,
			r->helix.initial.name, &r->helix.initial.chain_id, &r->helix.initial.seq, &r->helix.initial.insert,
			r->helix.terminal.name, &r->helix.terminal.chain_id, &r->helix.terminal.seq, &r->helix.terminal.insert,
			&r->helix.helixclass, r->helix.comment);
		break;

	case PDB_HET:
		strscan(str, fmt,
			r->het.residue.name, &r->het.residue.chain_id, &r->het.residue.seq, &r->het.residue.insert,
			&r->het.natoms, r->het.text);
		break;

	case PDB_MASTER:
		strscan(str, fmt,
			&r->master.nremark,
			&r->master.nftnote,
			&r->master.nhet,
			&r->master.nhelix,
			&r->master.nsheet,
			&r->master.nturn,
			&r->master.nsite,
			&r->master.ntransform,
			&r->master.natom,
			&r->master.nter,
			&r->master.nconect,
			&r->master.nseqres);
		break;

	case PDB_MODEL:
		strscan(str, fmt, &r->model.serial);
		break;

	case PDB_MTRIX:
		strscan(str, fmt, &r->mtrix.row, &r->mtrix.serial,
			&r->mtrix.m1, &r->mtrix.m2, &r->mtrix.m3, &r->mtrix.v, &r->mtrix.given);
		break;

	case PDB_OBSLTE:
		strscan(str, fmt, &r->obslte.cont, r->obslte.date, r->obslte.old_id,
			r->obslte.new_id[0], r->obslte.new_id[1],
			r->obslte.new_id[2], r->obslte.new_id[3],
			r->obslte.new_id[4], r->obslte.new_id[2],
			r->obslte.new_id[6], r->obslte.new_id[7]);
		break;

	case PDB_ORIGX:
		strscan(str, fmt, &r->origx.row,
			&r->origx.o1, &r->origx.o2, &r->origx.o3, &r->origx.t);
		break;

	case PDB_REVDAT:
		strscan(str, fmt, &r->revdat.mod, &r->revdat.cont, r->revdat.date,
			r->revdat.id, &r->revdat.modtype, r->revdat.correct);
		break;

	case PDB_SCALE:
		strscan(str, fmt, &r->scale.row,
			&r->scale.s1, &r->scale.s2, &r->scale.s3, &r->scale.u);
		break;

	case PDB_SEQRES:
		strscan(str, fmt, &r->seqres.serial, &r->seqres.chain_id, &r->seqres.nresidues,
			r->seqres.names[0], r->seqres.names[1],  r->seqres.names[2],
			r->seqres.names[3], r->seqres.names[4],  r->seqres.names[5],
			r->seqres.names[6], r->seqres.names[7],  r->seqres.names[8],
			r->seqres.names[9], r->seqres.names[10], r->seqres.names[11],
			r->seqres.names[12]);
		break;

	case PDB_SHEET:
		sh = &r->sheet;
		strscan(str, fmt, &sh->strand, sh->id, &sh->nstrands,
			sh->initial.name, &sh->initial.chain_id, &sh->initial.seq, &sh->initial.insert,
			sh->terminal.name, &sh->terminal.chain_id, &sh->terminal.seq, &sh->terminal.insert,
			&sh->sense,
			sh->current.name,
			sh->current.residue.name,
			&sh->current.residue.chain_id,
			&sh->current.residue.seq,
			&sh->current.residue.insert,
			sh->prev.name,
			sh->prev.residue.name,
			&sh->prev.residue.chain_id,
			&sh->prev.residue.seq,
			&sh->prev.residue.insert);
		break;

	case PDB_SITE:
		strscan(str, fmt, &r->site.seq, r->site.id, &r->site.nresidues,
			r->site.residue1.name, &r->site.residue1.chain_id, &r->site.residue1.seq, &r->site.residue1.insert,
			r->site.residue2.name, &r->site.residue2.chain_id, &r->site.residue2.seq, &r->site.residue2.insert,
			r->site.residue3.name, &r->site.residue3.chain_id, &r->site.residue3.seq, &r->site.residue3.insert,
			r->site.residue4.name, &r->site.residue4.chain_id, &r->site.residue4.seq, &r->site.residue4.insert);
		break;

	case PDB_SPRSDE:
		strscan(str, fmt, &r->sprsde.cont, r->sprsde.date, r->sprsde.id,
			r->sprsde.old_id[0], r->sprsde.old_id[1],
			r->sprsde.old_id[2], r->sprsde.old_id[3],
			r->sprsde.old_id[4], r->sprsde.old_id[5],
			r->sprsde.old_id[6], r->sprsde.old_id[7]);
		break;

	case PDB_SSBOND:
		strscan(str, fmt, &r->ssbond.seq,
			r->ssbond.residue1.name, &r->ssbond.residue1.chain_id,
			&r->ssbond.residue1.seq, &r->ssbond.residue1.insert,
			r->ssbond.residue2.name, &r->ssbond.residue2.chain_id,
			&r->ssbond.residue2.seq, &r->ssbond.residue2.insert,
			r->ssbond.comment);
		break;

	case PDB_TER:
		strscan(str, fmt, &r->ter.serial,
			r->ter.residue.name, &r->ter.residue.chain_id,
			&r->ter.residue.seq, &r->ter.residue.insert);
		break;

	case PDB_TURN:
		strscan(str, fmt, &r->turn.seq, r->turn.id,
			r->turn.residue1.name, &r->turn.residue1.chain_id, &r->turn.residue1.seq, &r->turn.residue1.insert,
			r->turn.residue2.name, &r->turn.residue2.chain_id, &r->turn.residue2.seq, &r->turn.residue2.insert,
			r->turn.comment);
		break;

	case PDB_TVECT:
		strscan(str, fmt, &r->tvect.serial,
			&r->tvect.t1, &r->tvect.t2, &r->tvect.t3, r->tvect.comment);
		break;
	}

	return 1;
}

static char	*scanfmt[PDB_MAX] = {
	"%80c",							/* (j)unk */
	"%6 %5d %4c%c%3c %c%4d%c %7d%7d%7d%7d%7d%7d",		/* ANISOU */
	"%6 %5d %4c%c%3c %c%4d%c   %8lf%8lf%8lf%6f%6f %3d",	/* ATOM   */
	"%9 %c%60c",						/* AUTHOR */
	"%9 %c%60c",						/* COMPND */
	"%6 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",			/* CONECT */
	"%6 %9f%9f%9f%7f%7f%7f %11c%4d",			/* CRYST1 */
	"",							/* END    */
	"",							/* ENDMDL */
	"%9 %c%60c",						/* EXPDTA */

	"%8 %2d  %3c %2d%c%51c",				/* FORMUL */
	"%7 %3d %59c",						/* FTNOTE */
	"%10 %40c%9c%3 %4c",					/* HEADER */
	"%7 %3d %3c %3c %c %4d%c %3c %c %4d%c%2d%30c",		/* HELIX  */
	"%7 %3c  %c%4d%c  %5d%5 %40c",				/* HET    */
	"%6 %5d %4c%c%3c %c%4d%c   %8lf%8lf%8lf%6f%6f %3d",	/* HETATM */
	"%9 %c%60c",						/* JRNL   */
	"%10 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",		/* MASTER */
	"%9 %5d",						/* MODEL  */
	"%5 %d %3d%10f%10f%10f%5 %10f   %2d",			/* MTRIX  */

	"%8 %2d %9c %4c%6 %4c %4c %4c %4c %4c %4c %4c %4c",				/* OBSLTE */
	"%5 %d%4 %10f%10f%10f%5 %10f",							/* ORIGX  */
	"%7 %3d %59c",									/* REMARK */
	"%7 %3d%2d %9c %7c %c%7 %31c",							/* REVDAT */
	"%5 %d%4 %10f%10f%10f%5 %10f",							/* SCALE  */
	"%6 %4d %c %4d  %3c %3c %3c %3c %3c %3c %3c %3c %3c %3c %3c %3c %3c",		/* SEQRES */
	"%6 %4d %3c%2d %3c %c%4d%c %3c %c%4d%c%2d %4c%3c %c%4d%c %4c%3c %c%4d%c",	/* SHEET  */
	"%6 %5d %4c%c%3c %c%4d%c   %8lf%8lf%8lf%6f%6f %3d",				/* SIGATM */
	"%6 %5d %4c%c%3c %c%4d%c %7d%7d%7d%7d%7d%7d",					/* SIGUIJ */
	"%7 %3d %3c %2d %3c %c%4d%c %3c %c%4d%c %3c %c%4d%c %3c %c%4d%c",		/* SITE   */

	"%9 %c%60c",						/* SOURCE */
	"%8 %2d %9c %4c%6 %4c %4c %4c %4c %4c %4c %4c %4c",	/* SPRSDE */
	"%7 %3d %3c %c %4d%c   %3c %c %4d%c%4 %30c",		/* SSBOND */
	"%6 %5d%6 %3c %c%4d%c",					/* TER    */
	"%7 %3d %3c %3c %c%4d%c %3c %c%4d%c%4 %30c",		/* TURN   */
	"%7 %3d%10f%10f%10f%30c",				/* TVECT  */
	};

static char	*printfmt[PDB_MAX] = {
	"%s",								/* (j)unk */
	"ANISOU%5d %-4s%c%-3s %c%4d%c %7d%7d%7d%7d%7d%7d",		/* ANISOU */
	"ATOM  %5d %-4s%c%-3s %c%4d%c   %8.4f%8.4f%8.4f%6.2f%6.2f %3d",	/* ATOM   */
	"AUTHOR   %c%-60s",						/* AUTHOR */
	"COMPND   %c%-60s",						/* COMPND */
	"CONECT%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",			/* CONECT */
	"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d",		/* CRYST1 */
	"END",								/* END    */
	"ENDMDL",							/* ENDMDL */
	"EXPDTA   %c%-60s",						/* EXPDTA */

	"FORMUL  %2d  %-3s %2d%c%-51s",					/* FORMUL */
	"FTNOTE %3d %-59s",						/* FTNOTE */
	"HEADER    %-40s%-12s%-4s",					/* HEADER */
	"HELIX  %3d %3s %-3s %c %4d%c %-3s %c %4d%c%2d%-30s",		/* HELIX  */
	"HET    %-3s  %c%4d%c  %5d     %-40s",				/* HET    */
	"HETATM%5d %-4s%c%-3s %c%4d%c   %8.4f%8.4f%8.4f%6.2f%6.2f %3d",	/* HETATM */
	"JRNL     %c%-60s",						/* JRNL   */
	"MASTER    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",		/* MASTER */
	"MODEL    %5d",							/* MODEL  */
	"MTRIX%d %3d%10.5f%10.5f%10.5f     %10.5f   %2d",		/* MTRIX  */

	"OBSLTE  %2d %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s",				/* OBSLTE */
	"ORIGX%d    %10.5f%10.5f%10.5f     %10.5f",						/* ORIGX  */
	"REMARK %3d %-59s",									/* REMARK */
	"REVDAT %3d%2d %-9s %-7s %c       %-31s",						/* REVDAT */
	"SCALE%d    %10.5f%10.5f%10.5f     %10.5f",						/* SCALE  */
	"SEQRES%4d %c %4d  %-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-3s",		/* SEQRES */
	"SHEET %4d %3s%2d %-3s %c%4d%c %-3s %c%4d%c%2d %-4s%-3s %c%4d%c %-4s%-3s %c%4d%c",	/* SHEET  */
	"SIGATM%5d %-4s%c%-3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f %3d",				/* SIGATM */
	"SIGUIJ%5d %-4s%c%-3s %c%4d%c %7d%7d%7d%7d%7d%7d",					/* SIGUIJ */
	"SITE   %3d %3s %2d %-3s %c%4d%c %-3s %c%4d%c %-3s %c%4d%c %-3s %c%4d%c",		/* SITE   */

	"SOURCE   %c%-60s",									/* SOURCE */
	"SPRSDE  %2d %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s",				/* SPRSDE */
	"SSBOND %3d %3s %c %4d%c   %3s %c %4d%c    %-30s",					/* SSBOND */
	"TER   %5d      %-3s %c%4d%c",								/* TER    */
	"TURN   %3d %3s %-3s %c%4d%c %-3s %c%4d%c    %-30s",					/* TURN   */
	"TVECT  %3d%10.5f%10.5f%10.5f%-30s",							/* TVECT  */
	};

char	*getpdbrecordscanfmt (int type)
{
	if (type < 0 || type >= PDB_MAX) type = 0;
	return scanfmt[type];
}

char	*getpdbrecordprintfmt (int type)
{
	if (type < 0 || type >= PDB_MAX) type = 0;
	return printfmt[type];
}

void	fprintpdbrecord (PdbRecord *record, FILE *fp)
{
	int	i, n;
	char	*fmt;
	char	str[82];
	register PdbRec	*r;
	PdbSheet	*sh;

	if (!record) return;
	r = &record->rec;

	BZERO(str, sizeof(str));
	fmt = getpdbrecordprintfmt(record->type);

	switch (record->type) {
	case PDB_UNK:
		sprintf(str, fmt, r->unk.text);
		break;

	case PDB_ANISOU:
	case PDB_SIGUIJ:
		sprintf(str, fmt, r->anisou.serial, r->anisou.name, r->anisou.altloc, r->anisou.residue.name,
			r->anisou.residue.chain_id, r->anisou.residue.seq, r->anisou.residue.insert,
			r->anisou.u[0], r->anisou.u[1], r->anisou.u[2],
			r->anisou.u[3], r->anisou.u[4], r->anisou.u[5]);
		break;

	case PDB_ATOM:
	case PDB_HETATM:
	case PDB_SIGATM:
		sprintf(str, fmt, r->atom.serial, r->atom.name, r->atom.altloc,
			r->atom.residue.name, r->atom.residue.chain_id,
			r->atom.residue.seq, r->atom.residue.insert,
			r->atom.x, r->atom.y, r->atom.z,
			r->atom.occupancy, r->atom.tempfac, r->atom.footnote);
		break;

	case PDB_AUTHOR:
	case PDB_COMPND:
	case PDB_JRNL:
	case PDB_SOURCE:
		sprintf(str, fmt, r->author.cont, r->author.text);
		break;

	case PDB_CONECT:
		sprintf(str, fmt, r->conect.serial,
			r->conect.covalent[0], r->conect.covalent[1],
			r->conect.covalent[2], r->conect.covalent[3],
			r->conect.donor[0], r->conect.donor[1], r->conect.salt_minus,
			r->conect.acceptor[0], r->conect.acceptor[1], r->conect.salt_plus);
		break;

	case PDB_CRYST1:
		sprintf(str, fmt, r->cryst1.a, r->cryst1.b, r->cryst1.c,
			r->cryst1.alpha, r->cryst1.beta, r->cryst1.gamma,
			r->cryst1.spacegrp, r->cryst1.z);
		break;

	case PDB_END:
		sprintf(str, fmt);
		break;

	case PDB_ENDMDL:
		sprintf(str, fmt);
		break;

	case PDB_FORMUL:
		sprintf(str, fmt, r->formul.component, r->formul.hetid, r->formul.cont,
		r->formul.exclude, r->formul.formula);
		break;

	case PDB_FTNOTE:
	case PDB_REMARK:
		sprintf(str, fmt, r->ftnote.num, r->ftnote.text);
		break;

	case PDB_HEADER:
		sprintf(str, fmt, r->header.funcclass, r->header.date,
			r->header.id);
		break;

	case PDB_HELIX:
		sprintf(str, fmt, r->helix.serial, r->helix.id,
			r->helix.initial.name, r->helix.initial.chain_id, r->helix.initial.seq, r->helix.initial.insert,
			r->helix.terminal.name, r->helix.terminal.chain_id, r->helix.terminal.seq, r->helix.terminal.insert,
			r->helix.helixclass, r->helix.comment);
		break;

	case PDB_HET:
		sprintf(str, fmt,
			r->het.residue.name, r->het.residue.chain_id, r->het.residue.seq, r->het.residue.insert,
			r->het.natoms, r->het.text);
		break;

	case PDB_MASTER:
		sprintf(str, fmt,
			r->master.nremark,
			r->master.nftnote,
			r->master.nhet,
			r->master.nhelix,
			r->master.nsheet,
			r->master.nturn,
			r->master.nsite,
			r->master.ntransform,
			r->master.natom,
			r->master.nter,
			r->master.nconect,
			r->master.nseqres);
		break;

	case PDB_MODEL:
		sprintf(str, fmt, r->model.serial);
		break;

	case PDB_MTRIX:
		sprintf(str, fmt, r->mtrix.row, r->mtrix.serial,
			r->mtrix.m1, r->mtrix.m2, r->mtrix.m3, r->mtrix.v, r->mtrix.given);
		break;

	case PDB_OBSLTE:
		sprintf(str, fmt, r->obslte.cont, r->obslte.date, r->obslte.old_id,
			r->obslte.new_id[0], r->obslte.new_id[1],
			r->obslte.new_id[2], r->obslte.new_id[3],
			r->obslte.new_id[4], r->obslte.new_id[2],
			r->obslte.new_id[6], r->obslte.new_id[7]);
		break;

	case PDB_ORIGX:
		sprintf(str, fmt, r->origx.row,
			r->origx.o1, r->origx.o2, r->origx.o3, r->origx.t);
		break;

	case PDB_REVDAT:
		sprintf(str, fmt, r->revdat.mod, r->revdat.cont, r->revdat.date,
			r->revdat.id, r->revdat.modtype, r->revdat.correct);
		break;

	case PDB_SCALE:
		sprintf(str, fmt, r->scale.row,
			r->scale.s1, r->scale.s2, r->scale.s3, r->scale.u);
		break;

	case PDB_SEQRES:
		sprintf(str, fmt, r->seqres.serial, r->seqres.chain_id, r->seqres.nresidues,
			r->seqres.names[0], r->seqres.names[1],  r->seqres.names[2],
			r->seqres.names[3], r->seqres.names[4],  r->seqres.names[5],
			r->seqres.names[6], r->seqres.names[7],  r->seqres.names[8],
			r->seqres.names[9], r->seqres.names[10], r->seqres.names[11],
			r->seqres.names[12]);
		break;

	case PDB_SHEET:
		sh = &r->sheet;
		sprintf(str, fmt, sh->strand, sh->id, sh->nstrands,
			sh->initial.name, sh->initial.chain_id, sh->initial.seq, sh->initial.insert,
			sh->terminal.name, sh->terminal.chain_id, sh->terminal.seq, sh->terminal.insert,
			sh->sense,
			sh->current.name,
			sh->current.residue.name,
			sh->current.residue.chain_id,
			sh->current.residue.seq,
			sh->current.residue.insert,
			sh->prev.name,
			sh->prev.residue.name,
			sh->prev.residue.chain_id,
			sh->prev.residue.seq,
			sh->prev.residue.insert);
		break;

	case PDB_SITE:
		sprintf(str, fmt, r->site.seq, r->site.id, r->site.nresidues,
			r->site.residue1.name, r->site.residue1.chain_id, r->site.residue1.seq, r->site.residue1.insert,
			r->site.residue2.name, r->site.residue2.chain_id, r->site.residue2.seq, r->site.residue2.insert,
			r->site.residue3.name, r->site.residue3.chain_id, r->site.residue3.seq, r->site.residue3.insert,
			r->site.residue4.name, r->site.residue4.chain_id, r->site.residue4.seq, r->site.residue4.insert);
		break;

	case PDB_SPRSDE:
		sprintf(str, fmt, r->sprsde.cont, r->sprsde.date, r->sprsde.id,
			r->sprsde.old_id[0], r->sprsde.old_id[1],
			r->sprsde.old_id[2], r->sprsde.old_id[3],
			r->sprsde.old_id[4], r->sprsde.old_id[5],
			r->sprsde.old_id[6], r->sprsde.old_id[7]);
		break;

	case PDB_SSBOND:
		sprintf(str, fmt, r->ssbond.seq,
			r->ssbond.residue1.name, r->ssbond.residue1.chain_id,
			r->ssbond.residue1.seq, r->ssbond.residue1.insert,
			r->ssbond.residue2.name, r->ssbond.residue2.chain_id,
			r->ssbond.residue2.seq, r->ssbond.residue2.insert,
			r->ssbond.comment);
		break;

	case PDB_TER:
		sprintf(str, fmt, r->ter.serial,
			r->ter.residue.name, r->ter.residue.chain_id,
			r->ter.residue.seq, r->ter.residue.insert);
		break;

	case PDB_TURN:
		sprintf(str, fmt, r->turn.seq, r->turn.id,
			r->turn.residue1.name, r->turn.residue1.chain_id, r->turn.residue1.seq, r->turn.residue1.insert,
			r->turn.residue2.name, r->turn.residue2.chain_id, r->turn.residue2.seq, r->turn.residue2.insert,
			r->turn.comment);
		break;

	case PDB_TVECT:
		sprintf(str, fmt, r->tvect.serial,
			r->tvect.t1, r->tvect.t2, r->tvect.t3, r->tvect.comment);
		break;
	}

	for(i=0,n=0;i<sizeof(str);i++) {
		if (!str[i]) {
			n = i;
			str[i] = ' ';
		}
	}
	str[n] = '\0';

	fprintf(fp, "%s\n", str);

}

