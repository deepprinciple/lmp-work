#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "molecule.h"

static SS	*sslist = NULL;
PdbCryst1	cryst;	/* crystal info */

static ResiduePtr	cur_residue = NULL;
static ChainPtr	cur_chain = NULL;
static AtomPtr	conn_atom = NULL;
static char	cur_insert = '\0';

static void	ProcessSS(MolPtr m);
static void	ProcessPdbAtom(MolPtr m, PdbAtom *patom, int hetatm);
static void	ProcessPdbMol (MolPtr mlist);

#include <time.h>
static time_t	laptime = 0;
static void	bgntime (void)
{
	time(&laptime);
}

static int	endtime (void)
{
	return (int)(time(NULL) - laptime);
}

/********************************************
*** read PDB file from a file pointer "fp".
*** returns a list of molecules.
*********************************************/

MolPtr	FLoadPdbFile (FILE *fp)
{
	int	i, type, line;
	PdbRecord	_record, *record = &_record;
	PdbRec	*rec;
	PdbConect	*con;
	SS	*ss;
	AtomPtr	atom1, atom2;
	BondPtr	bond;

	Mol	*mlist = NULL, *m = NULL;

	/* clean up and initialize */
	sslist = NULL;
	cur_residue = NULL;
	cur_chain = NULL;
	conn_atom = NULL;
	cur_insert = '\0';

bgntime();
	clearmolerror();
	line = 0;
	while (1) {
		line++;
		*record = readpdbrecord(fp);
		rec = &record->rec;
		type = record->type;
		if (type == PDB_END) break;

		switch (type) {
		case PDB_UNK:	/* not a critical error. no reason to return */
			break;
		case PDB_ATOM:
		case PDB_HETATM:
			if (!m && !(m = EnterNewMol(&mlist))) {
				setmolerrortext("Error in memory allocation");
				setmolerrorno(MERR_MEM);
				break;
			}
			ProcessPdbAtom(m, &rec->atom, type == PDB_HETATM);
			break;
		case PDB_CONECT:	/* bond */
			if (!m) {
				setmolerrortext("Error in PDB file format: Line #%d", line);
				setmolerrorno(MERR_FORMAT);
				break;
			}
			con = &rec->conect;
			if (!(atom1 = GetSernoAtom(m->chain, con->serial))) break;
			if (atom1->nneighbors >= MAXNEIGHBOR) break;
			for(i=0;i<4;i++) {
				if (con->covalent[i] <= 0) continue;
				if (!(atom2 = GetSernoAtom(m->chain, con->covalent[i]))) continue;
				if (atom2->nneighbors >= MAXNEIGHBOR) continue;
				if (GetAtomBonded(m->bond, atom1, atom2)) continue;
				bond = EnterNewBond(&m->bond);
				bond->atom1 = atom1;
				bond->atom2 = atom2;
				bond->type = B_SINGLE;
			}
			break;
		case PDB_COMPND:	/* compound name */
			break;
		case PDB_CRYST1:
			BCOPY(&rec->cryst1, &cryst, sizeof(PdbCryst1));
			break;
		case PDB_ENDMDL:	/* end of the current molecule */
			break;
		case PDB_TER:
			conn_atom = NULL;
			cur_residue = NULL;
			cur_chain = NULL;
			break;
		case PDB_HEADER:
			if (!m && !(m = EnterNewMol(&mlist))) {
				setmolerrortext("Error in memory allocation");
				setmolerrorno(MERR_MEM);
				break;
			}
			if (!m->header) m->header = LookUpStrTable(rec->header.funcclass);
			break;
		case PDB_HELIX:
			ss = EnterNewSS(&sslist);
			ss->type  = R_HELIX;
			ss->chain_id = rec->helix.initial.chain_id;
			ss->init  = rec->helix.initial.seq;
			ss->term  = rec->helix.terminal.seq;
			break;
		case PDB_SHEET:
			ss = EnterNewSS(&sslist);
			ss->type  = R_SHEET;
			ss->chain_id = rec->sheet.initial.chain_id;
			ss->init  = rec->sheet.initial.seq;
			ss->term  = rec->sheet.terminal.seq;
			break;
		case PDB_TURN:
			ss = EnterNewSS(&sslist);
			ss->type  = R_TURN;
			ss->chain_id = rec->turn.residue1.chain_id;
			ss->init  = rec->turn.residue1.seq;
			ss->term  = rec->turn.residue2.seq;
			break;
		case PDB_MODEL:
			break;
		default:	/* ignore anything else */
			break;
		}

		if (getmolerrorno() != -1) break;
	}

	if (getmolerrorno() != -1) {	/* error has occurred */
		if (mlist) FreeMol(&mlist);
		if (sslist) FreeSS(&sslist);
		return NULL;
	}

	if (mlist) ProcessPdbMol(mlist);

/*
fprintf(stderr, "laptime = %d\n", endtime());
*/

	return mlist;
}

static void	ProcessPdbMol (MolPtr mlist)
{
	MolPtr		m;
	ChainPtr	c;
	ResiduePtr	r;


	/* Build bonds between backbone breaks */
	ForEachMol(mlist,m) {
		/* CalcBonds() builds only bond list. atom->{nbatom, nbbond, nneighbors}
		 * are not updated. On the other hand, CalcBondsWithinResidues
		 * needs a neighbor information. Therefore, GetNeighbor must be
		 * called before calling CalcBondsWithinResidues.
		 */
		GetNeighbor(m);
		CalcBondsWithinResidues(m);

		ForEachChain(m->chain,c) {
			for(r=c->residue;r;r=r->next) {
				if (r->backbone || !r->next) continue;
				if (r->prev && r->next->backbone) continue;
				CalcBondsBetweenTwoResidues(r, r->next, &m->bond);
			}
		}
		FixConnectionTable(m);
	}

	/* Build SS list */
	ProcessSS(mlist);
	
	/* we don't need sslist anymore */
	if (sslist) FreeSS(&sslist);
}

static void	ProcessSS (MolPtr m)
{
	Chain	*chain, *c, *cnext;
	Residue	*r;
	SS	*ss, *ssnext;
	int	prev_seqno;

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "Processing SS...\n");
}

	/* check if residue sequence numbers are in ascending order. If not,
	 * give up building an SS list.
	 */
	ForEachChainRes(m->chain,c,r) {
		if (r == c->residue) {
			prev_seqno = r->seqno;
			continue;
		}
		if (r->seqno < prev_seqno) {

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "Sequence numbers are not in ascending order in chain %c: %d -> %d\n",
		c->id, prev_seqno, r->seqno);
}

			return;
		}
		prev_seqno = r->seqno;
	}

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "Helix (=%d), Sheet (=%d), and Turn (=%d)\n", R_HELIX, R_SHEET, R_TURN);
}

	for(ss=sslist;ss;ss=ssnext) {
		ssnext = ss->next;

if (getmoldebuglevel() > 0) {
	fprintf(stderr, "type = %d chain_id = %c init = %3d term = %3d\n",
	ss->type, ss->chain_id, ss->init, ss->term);
}

		if (ss->term - ss->init + 1 < 3) {	/* delete it if too short */
			RemoveSS(ss, &sslist);
		}
	}

	for(ss=sslist;ss;ss=ss->next) {
		if ((chain = FindChainByID(m->chain, ss->chain_id))) {
			cnext = chain->next;
			chain->next = NULL;
		} else {
			chain = m->chain;
			cnext = NULL;
		}

		ForEachChain(chain,c) {
			for(r=chain->residue;r && r->seqno<ss->init;r=r->next);
			if (r) r->ss_flags |= R_SSINIT;
			for(;r && r->seqno <= ss->term;r=r->next) {
				if (r->seqno == ss->term) r->ss_flags |= R_SSTERM;
				r->ss_flags |= ss->type;
			}
		}
		if (cnext) chain->next = cnext;
	}

	m->ss = MakeSSList(m->chain);
}

static void	ProcessPdbResidue (MolPtr m, PdbResidue *pres, int hetatm)
{
	cur_insert = pres->insert;
	if (!cur_chain || cur_chain->id != pres->chain_id) cur_chain = EnterNewChain(&m->chain);
	cur_chain->id = pres->chain_id;
	cur_residue = EnterNewResidue(&cur_chain->residue);
	cur_residue->refno = GetResRefno(pres->name);
	cur_residue->seqno = pres->seq;
	cur_residue->chain = cur_chain;
	if (hetatm) cur_residue->flags |= R_HETATM;
}

#define BACKBONE_SQDIST	(7*7)

static void	ProcessPdbAtom (MolPtr m, PdbAtom *patom, int hetatm)
{
	PdbResidue	*pres;
	AtomPtr	atom, C, N;
	BondPtr	bond;
	float	x, y, z;

	pres = &patom->residue;

	x = patom->x;
	y = patom->y;
	z = patom->z;

	if (!cur_residue || cur_residue->seqno != pres->seq ||
		cur_chain->id != pres->chain_id || cur_insert != pres->insert) {
		ProcessPdbResidue(m, pres, hetatm);
	}

	atom = EnterNewAtom(&cur_residue->atom);
	if (hetatm) atom->flags |= A_HETATM;
	atom->refno   = GetPdbAtomRefno(patom->name, cur_residue->refno);
	atom->residue = cur_residue;

	atom->an      = GetPdbAtomicNumber(atom->refno, cur_residue->refno);
	atom->serno   = patom->serial;
	atom->tempfac = patom->tempfac;
	atom->altloc  = patom->altloc;
	atom->x = patom->x;
	atom->y = patom->y;
	atom->z = patom->z;

	/* create a backbone */
	/*
	if (IsProteinAlphaCarbon(cur_residue->refno, atom->refno)) {
	*/
	if (IsAlphaCarbon(atom->refno)) {
		cur_residue->backbone = atom;
		/* get the distance between the previous and current backbone atoms */
		if (conn_atom) {
			if (ATOMSQDIST(conn_atom, atom) < BACKBONE_SQDIST) {
				/* create a bond */
				/* find C in conn_atom and N in atom */
				if ((C = FindAtomInRes(conn_atom->residue, ATM_C)) &&
				    (N = FindAtomInRes(atom->residue, ATM_N))) {
					bond = EnterNewBond(&m->bond);
					bond->atom1 = C;
					bond->atom2 = N;
					bond->type = B_SINGLE;
				}
			} else {
				atom->flags |= A_BREAK;
				cur_residue->ss_flags |= R_BREAK;
			}
		}
		conn_atom = atom;
	} else if (IsNucleoSugarPhosphate(cur_residue->refno, atom->refno)) {
		cur_residue->backbone = atom;
		if (conn_atom) {
			/* create a bond */
		}
		conn_atom = atom;
	}
}

void	FPrintPdbMol (FILE *fp, MolPtr mlist)
{
	Mol	*m;
	char	*p;
	AtomPtr		atom;
	ChainPtr	chain;
	ResiduePtr	res;
	PdbAtom	*patom;
	PdbRecord	record;
	PdbResidue	*pres;
	PdbRemark	*premark;

	if (!mlist) return;
	BZERO(&record, sizeof(PdbRecord));
	patom = &record.rec.atom;
	pres = &patom->residue;

	record.type = PDB_REMARK;
	premark = &record.rec.remark;

	premark->num = 0;
	p = GetStrValue(mlist->header);
	strcpy(premark->text, "Title: ");
	if (p) {
		strncat(premark->text, p, PDB_REMARK_LEN-strlen(premark->text));
		premark->text[PDB_REMARK_LEN] = '\0';
	} else {
		strcat(premark->text, "PDB");
	}

	fprintpdbrecord(&record, fp);

	premark->num++;
	strcpy(premark->text, "Creator: XChemEdit/AUTOZMAT");
	fprintpdbrecord(&record, fp);

	premark->num++;
	strncpy(premark->text, "WL Jorgensen Group, Dept. of Chemistry, Yale University",
		PDB_REMARK_LEN);
	fprintpdbrecord(&record, fp);

	for(m=mlist;m;m=m->next) {
		/* secondary structures */
		if (m->ss) {
			SS		*ss;
			PdbSheet	*sheet;
			PdbHelix	*helix;
			PdbTurn		*turn;
			PdbResidue	*init, *term;
			Residue	*initres, *termres;
			Chain	*cnext;
			int	helix_count, sheet_count, turn_count, nstrands;

			sheet = &record.rec.sheet;
			helix = &record.rec.helix;
			turn  = &record.rec.turn;
			helix_count = 0;
			sheet_count = 0;
			turn_count  = 0;

			for(ss=m->ss,nstrands=0;ss;ss=ss->next) if (ss->type & R_SHEET) nstrands++;
			for(ss=m->ss;ss;ss=ss->next) {
				if (!(ss->type & (R_HELIX|R_SHEET|R_TURN))) continue;
				if ((chain = FindChainByID(m->chain, ss->chain_id))) {
					cnext = chain->next;
					chain->next = NULL;
				} else {
					chain = m->chain;
					cnext = NULL;
				}
				initres = FindResBySeqno(chain->residue, ss->init);
				termres = FindResBySeqno(chain->residue, ss->term);
				if (!initres || !termres) {
					if (cnext) chain->next = cnext;
					continue;
				}

				BZERO(&record, sizeof(PdbRecord));
				switch (ss->type) {
				case R_HELIX:
					helix_count++;
					init = &helix->initial;
					term = &helix->terminal;

					record.type = PDB_HELIX;
					helix->serial = helix_count;
					init->chain_id = ss->chain_id;
					init->seq = ss->init;
					term->seq = ss->term;
					break;
				case R_SHEET:
					sheet_count++;
					init = &sheet->initial;
					term = &sheet->terminal;

					record.type = PDB_SHEET;
					sheet->strand = sheet_count;
					sheet->nstrands = nstrands;
					init->chain_id = ss->chain_id;
					init->seq = ss->init;
					term->seq = ss->term;
					break;
				case R_TURN:
					turn_count++;
					init = &turn->residue1;
					term = &turn->residue2;

					record.type = PDB_TURN;
					turn->seq = turn_count;
					init->chain_id = ss->chain_id;
					init->seq = ss->init;
					term->seq = ss->term;
					break;
				}
				strcpy(init->name, GetResName(initres->refno));
				strcpy(term->name, GetResName(termres->refno));
				fprintpdbrecord(&record, fp);

				if (cnext) chain->next = cnext;
			}
		}

		/* ATOM and HETATM entries */
		for(chain=m->chain;chain;chain=chain->next) {
			pres->chain_id = chain->id;
			for(res=chain->residue;res;res=res->next) {
				pres->seq = res->seqno;
				strcpy(pres->name, GetResName(res->refno));
				pres->insert = ' ';

				if (strspn(pres->name, " \t") == strlen(pres->name)) {
					strcpy(pres->name, "UNK");
				}

				for(atom=res->atom;atom;atom=atom->next) {
					record.type = (atom->flags & A_HETATM) ? PDB_HETATM : PDB_ATOM;
					strcpy(patom->name, MakePdbAtomName(GetAtomName(atom->refno), atom->an));

					patom->serial = atom->serno;
					patom->altloc = atom->altloc;
					patom->x = atom->x;
					patom->y = atom->y;
					patom->z = atom->z;
					patom->occupancy = 0.0;
					patom->tempfac = atom->tempfac;
					patom->footnote = 0;

					fprintpdbrecord(&record, fp);
				}
			}
			fprintf(fp, "TER\n");
		}

		/* CONECT entries */
		if (m->bond) {
			PdbConect	*conect;
			int	nbonds, i;
			AtomPtr	nbatom;

			conect = &record.rec.conect;
			ForEachChainResAtom(m->chain, chain, res, atom) {
				BZERO(&record, sizeof(PdbRecord));
				record.type = PDB_CONECT;
				conect->serial = atom->serno;
				for(i=0,nbonds=0;i<atom->nneighbors;i++) {
					nbatom = atom->nbatom[i];
					if (!(atom->flags & A_HETATM) && !(nbatom->flags & A_HETATM)) continue;
					conect->covalent[nbonds] = nbatom->serno;
					nbonds++;
					if (nbonds >= 4) break;
				}
				if (nbonds > 0) fprintpdbrecord(&record, fp);
			}
		}
		fprintf(fp, "ENDMDL\n");
	}
	fprintf(fp, "END\n");
}

ResiduePtr	CreatePdbResidue (ChainPtr *top_chn, ChainPtr *cur_chn, ResiduePtr *cur_res, PdbResidue *pres, int hetatm)
{
	ResiduePtr	r;

	if (!top_chn || !cur_chn || !cur_res) return NULL;
	if (!*cur_chn || (*cur_chn)->id != pres->chain_id) {
		*cur_chn = EnterNewChain(top_chn);
		(*cur_chn)->id = pres->chain_id;
	}
	if (*cur_res && (*cur_res)->seqno != pres->seq) {
		ForEachRes((*cur_chn)->residue,r) {
			if (r->seqno == pres->seq) {
				*cur_res = r;
				return *cur_res;
			}
		}
	}
	*cur_res = EnterNewResidue(&(*cur_chn)->residue);
	(*cur_res)->refno = GetResRefno(pres->name);
	(*cur_res)->seqno = pres->seq;
	(*cur_res)->chain = *cur_chn;
	if (hetatm) (*cur_res)->flags |= R_HETATM;

	return (*cur_res);
}

/* create a PDB structure from PdbResidue, PdbAtom, and hetatm */
AtomPtr	CreatePdbAtom (ChainPtr *top_chn, ChainPtr *cur_chn, ResiduePtr *cur_res, AtomPtr *conn_atm, PdbAtom *patom, BondPtr *bondlist, int hetatm)
{
	AtomPtr	atom, C, N;
	BondPtr	bond;
	PdbResidue	*pres;

	if (!top_chn || !cur_chn || !cur_res || !conn_atm || !patom) return NULL;
	pres = &patom->residue;

	if (!*cur_chn || (*cur_chn)->id != pres->chain_id ||
	    !*cur_res || (*cur_res)->seqno != pres->seq) {
		CreatePdbResidue(top_chn, cur_chn, cur_res, pres, hetatm);
	}

	atom = EnterNewAtom(&(*cur_res)->atom);
	if (hetatm) atom->flags |= A_HETATM;
	atom->refno = GetPdbAtomRefno(patom->name, (*cur_res)->refno);
	atom->residue = *cur_res;
	atom->serno = patom->serial;
	atom->x = patom->x;
	atom->y = patom->y;
	atom->z = patom->z;

/*
fprintf(stderr, "atom = %s%d\n", GetAtomName(atom->refno), atom->serno);
*/
	/* create a backbone */
	if (IsProteinAlphaCarbon((*cur_res)->refno, atom->refno)) {
/*
fprintf(stderr, "   alpha carbon = %s%d\n", GetAtomName(atom->refno), atom->serno);
*/
		(*cur_res)->backbone = atom;
		/* get the distance between the previous and current backbone atoms */
		if (*conn_atm) {
			if (bondlist && (ATOMSQDIST((*conn_atm), atom) < BACKBONE_SQDIST)) {
				/* create a bond */
				/* find C in conn_atom and N in atom */
				if ((C = FindAtomInRes((*conn_atm)->residue, ATM_C)) &&
				    (N = FindAtomInRes(atom->residue, ATM_N))) {
					bond = EnterNewBond(bondlist);
					bond->atom1 = C;
					bond->atom2 = N;
					bond->type = B_SINGLE;
				}
			} else {
				atom->flags |= A_BREAK;
				(*cur_res)->ss_flags |= R_BREAK;
			}
		}
		(*conn_atm) = atom;
	} else if (IsNucleoSugarPhosphate((*cur_res)->refno, atom->refno)) {
		(*cur_res)->backbone = atom;
		if (*conn_atm) {
			/* create a bond */
		}
		*conn_atm = atom;
	}
	return atom;
}

