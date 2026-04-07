#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "macros.h"
#include "molecule.h"

MolPtr	FLoadZmatCartesian (FILE *fp, ZmatPtr *zmatlist_ret)
{
	int	i, an, count, natoms, first, z1, z2, z3, zmat_err=0, charge_column=0;
	char	str[256], *data[32], *p, *c;
	Chain	*chain;
	Residue	*res;
	AtomPtr	atom, conn_atom;
	MolPtr	m;
	ZmatPtr	zmat;
	double	x;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	if (!fp) return NULL;
	if (!(m = NewMol())) return NULL;
	if (zmatlist_ret) *zmatlist_ret = NULL;

	chain = NULL;
	res = NULL;
	conn_atom = NULL;

	strcpy(patom->name, " ");
	strcpy(pres->name, "UNK");
	pres->seq = 0;
	natoms = 0;
	first = 1;
	while (1) {
		if (!READLINE(fp, str)) break;
		if ((c = strchr(str, '\n'))) *c = '\0';
		if (strstr(str, "$EndCoord:")) break;
		if (BLANKLINE(str)) continue;
		if (first && str[0] == '#') {
			if (strstr(str, "#ChargeColumn")) sscanf(str, "%*[^=]=%d", &charge_column);
			continue;
		}
		setmolerrortext(str);
		c = str;
		for(i=0;i<32;i++) data[i] = NULL;
		for(count=0;(p=strtok(c, " ,\t")) && count<32;c=(char *)NULL,count++) data[count] = p;

		p = data[0];
		if (p && IsInt(p)) sscanf(p, "%d", &an);
		else {
			if (!p || ((an = GetElemNumber(p)) == 0 && strcasecmp(p, "X") != 0)) {
				/* probably illegal atomic symbol */
				if (first) {
					first = 0;
					continue;
				}
				setmolerrorno(MERR_FORMAT);
				break;
			}
		}

		/* test if the line has a valid format */
		if (count < 4 || !(IsNumber(data[1]) && IsNumber(data[2]) && IsNumber(data[3]))) {
			if (first) {
				first = 0;
				continue;
			}
			setmolerrorno(MERR_FORMAT);
			break;
		}
		first = 0;

		natoms++;
		patom->serial = natoms;
		strcpy(patom->name, AN2ATOMNAME(an));
		sscanf(data[1], "%lf", &patom->x);
		sscanf(data[2], "%lf", &patom->y);
		sscanf(data[3], "%lf", &patom->z);

		atom = CreatePdbAtom(&m->chain, &chain, &res, &conn_atom, patom, &m->bond, 1);
		atom->an = an;

		z1 = z2 = z3 = 0;
		for(i=4;i<count;i++) {
			p = data[i];
			switch (i) {
			case 4:
				if (sscanf(p, "%lf", &x) != 1) break;
				atom->charge = atom->label1 = x;
				break;
			case 5:
				if (sscanf(p, "%lf", &x) != 1) break;
				atom->label2 = x;
				if (charge_column == 6) atom->charge = x;
				break;
			case 6:
			case 7:
			case 8:
				break;
			case 9:
				sscanf(p, "%d", &z1);
				break;
			case 10:
				sscanf(p, "%d", &z2);
				break;
			case 11:
				sscanf(p, "%d", &z3);
				break;
			}
		}

		if (!zmatlist_ret || zmat_err) continue;

		/* Build a Z-matrix */
		if (!CheckInternalCoord(natoms, z1, z2, z3, 1.0, 90.0, 0.0)) {
			if (zmatlist_ret && *zmatlist_ret) FreeZmat(zmatlist_ret);
			zmat_err = 1;
			continue;
		}
		zmat = NewZmatFromAtom(atom, z1, z2, z3);
		EnterZmat(zmat, zmatlist_ret);
	}

	if (zmatlist_ret && *zmatlist_ret) {
		ForEachZmat(*zmatlist_ret,zmat) zmat->flags = 0;
		AssignZmatProp(*zmatlist_ret);
	}

	if (natoms == 0) {
		FreeMol(&m);
		if (zmatlist_ret && *zmatlist_ret) FreeZmat(zmatlist_ret);
		return NULL;
	}

	GetNeighbor(m);
	CalcBondsWithinResidues(m);
	FixConnectionTable(m);
	return m;
}

MolPtr	FLoadCartesian (FILE *fp)
{
	return FLoadZmatCartesian(fp, NULL);
}

static MolPtr	_FLoadMindtool (FILE *fp, int flags)
{
	int	i, an, count, natoms, hetatm_flag;
	char	str[256], *data[32], *p, *q, *c, *c_end, chain_id, res_name[PDB_RSIZE], atom_name[PDB_ASIZE];
	Chain	*chain;
	Residue	*res;
	AtomPtr	atom, conn_atom;
	MolPtr	m;

	static PdbAtom	_patom, *patom=(&_patom);
	static PdbResidue	*pres=(&_patom.residue);

	if (!fp) return 0;
	if (!(m = NewMol())) return 0;

	chain = NULL;
	res = NULL;
	conn_atom = NULL;

	strcpy(patom->name, " ");
	strcpy(pres->name, "UNK");
	pres->seq = 1;
	pres->chain_id = '0';
	natoms = 0;

	while (READLINE(fp, str)) {
		if ((c = strchr(str, '\n'))) *c = '\0';
		if (strstr(str, "$EndCoord:")) break;
		if (natoms == 0 && strncmp(str, "Mind", 4) == 0) continue;
		if (BLANKLINE(str)) continue;

		c = str + strspn(str, " \t");
		c_end = strrchr(str, '\0');
		if (strchr("$", *c)) continue;
		if (strncmp(c, "----", 4) == 0) break;
		if (*c == '%') {
			pres->seq++;
			continue;
		}
		if (strncmp(c, "#% Chain=", 9) == 0) {	/* chain ID */
			if (sscanf(c, "%*[^=]=%c", &chain_id) == 1 && isascii(chain_id)) pres->chain_id = chain_id;
			else pres->chain_id++;
			continue;
		}
		if (strncmp(c, "#% Residue=", 11) == 0) {	/* residue name and sequence number */
			strcpy(res_name, "   ");
			if (sscanf(c, "%*[^=]=%3c", res_name) == 1) {
				for(i=0;i<PDB_RLEN;i++) if (!isascii(res_name[i])) res_name[i] = ' ';
			}
			res_name[PDB_RLEN] = '\0';
			if (IsBlank(res_name)) strcpy(res_name, "UNK");
			strcpy(pres->name, res_name);
			continue;
		}
		if (strchr("#", *c)) continue;	/* comment */

		setmolerrortext(str);
		for(i=0;i<32;i++) data[i] = NULL;
		
		atom_name[0] = '\0';
		for(count=0;(p=strtok(c, " ,\t")) && count<7;c=(char *)NULL,count++) {
			data[count] = p;
			if (count != 6) continue;
			q = p + strlen(p) + 1;

			strcpy(atom_name, "    ");
			if ((q >= c_end && !*q) || (sscanf(q, "%4c", atom_name) != 1)) continue;
			for(i=0;i<PDB_ALEN;i++) if (!isascii(atom_name[i])) atom_name[i] = ' ';
			atom_name[PDB_ALEN] = '\0';
		}

		/* test if the line has a valid format */
		if (count < 4 || !(IsInt(data[0]) && IsNumber(data[1]) && IsNumber(data[2]) && IsNumber(data[3]))) {
			setmolerrorno(MERR_FORMAT);
			break;
		}

		natoms++;
		patom->serial = natoms;

		sscanf(data[0], "%d", &an);
		if (IsBlank(atom_name)) strcpy(patom->name, AN2ATOMNAME(an));
		else strcpy(patom->name, atom_name);

		sscanf(data[1], "%lf", &patom->x);
		sscanf(data[2], "%lf", &patom->y);
		sscanf(data[3], "%lf", &patom->z);

		hetatm_flag = !IsAminoNucleo(GetResRefno(pres->name));
		atom = CreatePdbAtom(&m->chain, &chain, &res, &conn_atom, patom, &m->bond, hetatm_flag);
		atom->an = an;
		
		/* optional */
		if (count >= 5) sscanf(data[4], "%lf", &atom->charge);
		if (count >= 6) sscanf(data[5], "%lf", &atom->sigma);
		if (count >= 7) sscanf(data[6], "%lf", &atom->epsilon);
	}

	if (natoms == 0) {
		FreeMol(&m);
		return NULL;
	}

	GetNeighbor(m);
	CalcBondsWithinResidues(m);
	if (!(flags & MINDTOOL_BREAK_PERCENT)) CalcBondsBetweenResidues(m);
	FixConnectionTable(m);

	return m;
}

MolPtr	FLoadMindtool (FILE *fp, int flags)
{
	MolPtr	m, mlist=NULL;

	if (!fp) return NULL;
	while ((m = _FLoadMindtool(fp, flags))) EnterMol(m, &mlist);
	return mlist;
}

void	FPrintCartesian (FILE *fp, ChainPtr chainlist)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr		a;
	int	natoms=0;

	ForEachChainResAtom(chainlist,c,r,a) natoms++;
	if (natoms == 0) return;
	fprintf(fp, "%d\n", natoms);
	ForEachChainResAtom(chainlist,c,r,a) {
		fprintf(fp, "%3d % 12.6f % 12.6f % 12.6f\n", (int)a->an, a->x, a->y, a->z);
	}
}

void	FPrintMindtool (FILE *fp, ChainPtr chainlist)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr		a;
	int	natoms=0;

	if ((natoms = CountAtomInChains(chainlist)) == 0) return;

	fprintf(fp, "Mind\n");
	fprintf(fp, "# The '%%' and '----' signs separate residues and molecules, resp.\n");
	fprintf(fp, "# The '#%%' sign gives chain/residue info.\n");
	fprintf(fp, "# AN       X            Y            Z       Charge     Sigma  Epsilon Name\n");

	ForEachChain(chainlist,c) {
		fprintf(fp, "#%% Chain=%c\n", isascii(c->id) ? c->id : ' ');
		ForEachRes(c->residue,r) {
			if (r != chainlist->residue) fprintf(fp, "%%\n");
			fprintf(fp, "#%% Residue=%s   Seqno=%d\n", GetResName(r->refno), r->seqno);
			ForEachAtom(r->atom,a) {
				fprintf(fp, "%3d % 12.6f % 12.6f % 12.6f % 9.5f % 8.3f % 8.3f %s\n",
					a->an, a->x, a->y, a->z, a->charge, a->sigma, a->epsilon,
					GetAtomName(a->refno));
			}
		}
	}
}

