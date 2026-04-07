#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>

#include "mol/molecule.h"
#include "get_opt.h"
#include "error.h"
#include "oplstype.h"
#include "autozmat.h"

/* global variables */
int	autozmat_flags=0;
char	default_opls_par_file[PATH_MAX] = {'\0'};

static int	debug = 0;
static int	bosszmat_flags=0;	/* use masks for AddAllBossZmatOpt */
static int	pepz_bgn_atom_type=0;

static FileFormat	file_format[FILE_MAX] = FILE_FORMAT;

static ZOption	zoption_tab[] = ZOPTION_TAB;
#define ZOPTION_SIZE	(sizeof(zoption_tab)/sizeof(ZOption))

/* general options (-g option) */
static ZOption	option_tab[] = OPTION_TAB;
#define OPTION_SIZE	(sizeof(option_tab)/sizeof(ZOption))

static int	get_file_type (char *txt)
{
	int	i, n;
	char	*p;

	n = strlen(txt);
	if (n > 3) n = 4;

	for(i=0;i<FILE_MAX;i++) {
		p = file_format[i].text;
		if (strncasecmp(p, txt, n) == 0) return i;
	}
	return -1;
}

static void	print_usage (char *prog_name)
{
	int	i;
	ZOption	*zopt;

	PRINT("Usage: %s [-hHIrvVx] [-c command_file] [-i input_format] [-o output_format]\n", prog_name);
	PRINT("          [-p parameter_file] [-z zmat_option] [-f output_file] [-e errorlog_file]\n");
#ifdef HAVE_ATOMTYPE
	PRINT("          [-a] [-P opls_dir]\n");
#endif
	PRINT("          [-t template_zmatrix_file] [-T template_zmatrix_format] [filename]\n");

	PRINT("\nDescription:\n");
	PRINT("   '%s' converts between file formats. In particular, it automates\n", prog_name);
	PRINT("   a generation of Z-matrix for BOSS, MCPRO, GAUSSIAN, MOPAC, JAGUAR, etc.\n");
	PRINT("   The algorithm used here is a breadth-first search method in the graph\n");
	PRINT("   theory with criteria in such a way that improper torsions are handled\n");
	PRINT("   efficiently. In case of writing a BOSS or MCPRO Z-matrix, two dummies\n");
	PRINT("   are added at the central atom, which will appear in the beginning of the\n");
	PRINT("   Z-matrix. Dummies are also added wherever a linearity problem occurs,\n");
	PRINT("   e.g., the nitrile carbon atom in CH3-CN. Alternatively, Z-matrix can be\n");
	PRINT("   built using template if -t option is specified.\n");
	PRINT("   Hydrogen atoms can be added to vacant valencies to C, N, O, and S if '-H'\n");
	PRINT("   option is given. In case of atomtyping-enabled 'autozmat', OPLS atom and\n");
	PRINT("   torsion types are automatically assigned for BOSS/MCPRO Z-matrix.\n");
	PRINT("   For questions and bug reports, send email to lim@doctor.chem.yale.edu.\n");
	PRINT("\nOptions:\n");
	PRINT("  -i input_format   (default = pdb)\n");
	PRINT("  -o output_format  (default = qm_boss)\n");
	PRINT("\n");
	PRINT("      These options specify input or output format.\n");
	PRINT("      Valid formats are as listed below.\n");
	for(i=0;i<FILE_MAX;i++) {
	PRINT("        %-11s (%s)%s\n",
		file_format[i].text, file_format[i].desc, i!=FILE_MAX-1 ? "," : "");
	}
	PRINT("\n");
	PRINT("  -f output_file    Saves output in output_file.\n");
	PRINT("  -e errorlog_file  Saves errors in errorlog_file.\n");
	PRINT("  -c command_file   Executes commands in command_file.\n");
	PRINT("  -p param_file     Sets the OPLS parameter file to param_file.\n");
#ifdef HAVE_ATOMTYPE
	PRINT("  -a                Instructs autozmat not to use atomtyping.\n");
	PRINT("  -P opls_dir       Sets the OPLS directory for atomtyping.\n");
#endif
	PRINT("  -t template_file  Applys template Z-matrix file in Z-matrix generation.\n");
	PRINT("  -T template_fmt   Format of the template Z-matrix specified by -t (default=gauss)\n");
	PRINT("  -I                Begins interactive session.\n");
	PRINT("  -H                Adds hydrogens to C, N, O, and S.\n");
	PRINT("  -h                Prints out this information.\n");
	PRINT("  -v                Prints out version information.\n");
	PRINT("  -V                Prints out more verbose information.\n");
	PRINT("  -x                Removes dummy atoms before printing.\n");
	PRINT("  -z zmat_option    Yields a BOSS/MCPRO Z-matrix with options which may be delimited by comma.\n");
	PRINT("                    'zmat_option' may be one or a combination of\n");
	for(i=0,zopt=zoption_tab;i<ZOPTION_SIZE;i++,zopt++) {
	PRINT("          %-14s   %s\n", zopt->text, zopt->desc);
	}

	PRINT("\nExample:\n");
	PRINT("%% %s -i mdl file.mdl > file.boss\n", prog_name);
	PRINT("   reads file.mdl as MDL MOL format and redirects the output to\n");
	PRINT("   file.boss in BOSS Z-Matirx format for QM/MM calc.\n");
	PRINT("\n");
	PRINT("%% %s -i mdl -z default file.mdl > file.boss\n", prog_name);
	PRINT("   same as above but yields a fully flexible molecule.\n");
	PRINT("\n");
	PRINT("%% %s -i mdl -z \"addbond,addangle,adddihed\" file.mdl > file.boss\n", prog_name);
	PRINT("   same as above but yields a Z-matrix with additional bonds, angles, and dihedrals.\n");
	PRINT("\n");

	PRINT("%% %s -i boss -p file.par -o pdb file.boss > file.pdb\n", prog_name);
	PRINT("   reads file.boss as BOSS Z-matrix and file.par as BOSS parameter file\n");
	PRINT("   and redirects the output to file.pdb in PDB format.\n");
	PRINT("\n");
	PRINT("%% %s -i pdb -o mdl < file.pdb > file.mdl\n", prog_name);
	PRINT("   the input is redirected from file.pdb as PDB format and\n");
	PRINT("   the output is redirected to file.mdl in MDL MOL format.\n");

	PRINT("\n");
	PRINT("%% %s -i pdb -o boss -t tpl.boss -T boss file.pdb > file.boss\n", prog_name);
	PRINT("   reads \"file.pdb\" as PDB format and creates a BOSS Z-matrix using\n");
	PRINT("   template BOSS Z-matrix \"tpl.boss\".\n");
}

static char	*base_name (char *path)
{
	char	*indx;
	if ((indx = strrchr(path, DIRSEP))) return (indx+1);
	else return path;
}

static char	_dir[PATH_MAX];
static char	*dir_name (char *path)
{
	char	*idx;
	char	*dir = _dir;

	BZERO(_dir, sizeof(_dir));
	if ((idx = strrchr(path, DIRSEP))) {
		strncpy(dir, path, idx-path);
		dir[idx-path] = '\0';
	}
	return dir;
}

static int	is_readable_file (char *path)
{
	return (access(path, F_OK) == 0 && access(path, R_OK) == 0);
}

static void	get_library_dir (char *path, char dir[])
{
	char	*p;

	dir[0] = '\0';
	p = dir_name(path);

#ifndef MSDOS
	if (*p != '/') {
#else
	if (strncmp(p+1, ":\\", 2) != 0) {
#endif

		if (*p != '.') {
			if (!(access(path, F_OK) == 0 && access(path, X_OK) == 0)) {
				/* we should check PATH. give up */
				return;
			}
			sprintf(dir, "%s%c%s", getcwd(( char *)NULL, 256), DIRSEP, p[0] ? p : ".");
		} else {
			sprintf(dir, "%s%c%s", getcwd(( char *)NULL, 256), DIRSEP, p);
		}
		p = dir;
	}

	sprintf(dir, "%s%c..%clib%c", p, DIRSEP,DIRSEP,DIRSEP);
}

static void	fprint_boss_zmat_with_flags (FILE *fp, BossZmatPtr boss, int flags)
{
	int	saved_flags;

	if (!boss) return;
	flags &= (B_MCPRO|B_QM);

	saved_flags = boss->flags;
	boss->flags &= ~(B_MCPRO|B_QM);
	boss->flags |= flags;
	FPrintBossZmatWithTypes(fp, boss);
	boss->flags = saved_flags;
}

static void	fprint_pepz (FILE *fp, ChainPtr chainlist, BossZmatPtr boss, ListPtr atomtypes, int bgn_atom_type)
{
	int	type, flags=0;
	ZmatPtr	zmat, zmatlist;

	zmatlist = LinkZmatList(boss->zmatlist);
	if (bgn_atom_type > 0) {
		flags = PEPZ_ATOMTYPE_FROM_ZMAT;
		type = bgn_atom_type;
		ForEachZmat(zmatlist,zmat) {
			zmat->itype = type;
			type++;
		}
	}
	FPrintPepz(fp, chainlist, zmatlist, boss, atomtypes, flags);
	BreakZmatList(boss->zmatlist);
}

static int	parse_bosszmat_flags (char *str, int *flags)
{
	int	i;
	char	buf[256], *p, *q;
	ZOption	*zopt;

	if (!str || !*str) return 0;
	strcpy(buf, str);
	q = buf;
	while ((p = strtok(q, ",")) != (char *)NULL) {
		q = (char *)NULL;

		if (strcasecmp(p, "none") == 0) {
			*flags = 0;
			continue;
		}
		for(i=0,zopt=zoption_tab;i<ZOPTION_SIZE;i++,zopt++) {
			if (strcasecmp(p, zopt->text) == 0) {
				*flags |= zopt->flags;
				break;
			} else if (*p == '~' && strcasecmp(p+1, zopt->text) == 0) {
				*flags &= ~zopt->flags;
				break;
			}
		}
		if (i == ZOPTION_SIZE) PRINT("Unknown option '-z %s'\n", p);
	}

	return *flags;
}

static int	parse_autozmat_flags (char *str, int *flags)
{
	int	i;
	char	buf[256], var[256], val[256], *p, *q, *idx;
	ZOption	*opt;

	if (!str || !*str) return 0;
	strcpy(buf, str);
	q = buf;
	while ((p = strtok(q, ",")) != (char *)NULL) {
		q = (char *)NULL;

		if (strcasecmp(p, "none") == 0) {
			*flags = 0;
			continue;
		}
		for(i=0,opt=option_tab;i<OPTION_SIZE;i++,opt++) {
			if ((idx = strchr(p, '='))) {	/* A=B style option */
				var[0] = val[0] = '\0';
				sscanf(idx+1, "%s", val);
				*idx = '\0';
				sscanf(p, "%s", var);
				*idx = '\n';
				if (!var[0] || !val[0]) continue;
				if (strcasecmp(var, "pepz") == 0) sscanf(val, "%d", &pepz_bgn_atom_type);
				break;
			}

			if (strcasecmp(p, opt->text) == 0) {
				*flags |= opt->flags;
				break;
			} else if (*p == '~' && strcasecmp(p+1, opt->text) == 0) {
				*flags &= ~opt->flags;
				break;
			}
		}
		if (i == OPTION_SIZE) PRINT("Unknown option '-g %s'\n", p);
	}

	return *flags;
}

void	init_mol_info (MolPtr mlist, char *in_file, int hydrogen_flags, int guess_bondtype_flags)
{
	MolPtr	m;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	BondPtr	b;

	/* initialize molecule info. Also adds hydrogens if -H option is specified. */

	ForEachMol(mlist,m) {
		GetNeighbor(m);
		SetMolBBox3(m);
		ForEachChainResAtom(m->chain,c,r,a) a->col = GetElemCPKColor(a->an);

		if (m->header < 0) m->header = LookUpStrTable(base_name(in_file));
		if (hydrogen_flags & M_ADD_H) AddHydrogenToMol(m);
		if (guess_bondtype_flags) ForEachBond(m->bond,b) b->type = GuessBondType(b->atom1, b->atom2);

if (debug) PRINT("Input molecule = %s\n", GetStrValue(m->header));
	}
}

static void	break_residue_to_molecule (MolPtr *mlist)
{
	MolPtr	newmlist=NULL, m, newm;
	ChainPtr	c;
	ResiduePtr	r, rnext;
	BondPtr	bond, nextbond;

	ForEachMol(*mlist,m) {
		/* if bonded atoms belong to different residues, delete the bond
		 */
		for(bond=m->bond;bond;bond=nextbond) {
			nextbond = bond->next;
			if (bond->atom1->residue != bond->atom2->residue) {
				DeleteBond(bond, &m->bond);
				FreeBondEntry(&bond);
			}
		}
		ForEachChain(m->chain,c) {
			for(r=c->residue;r;r=rnext) {
				rnext = r->next;

				/* create a new molecule from this residue */
				newm = EnterNewMol(&newmlist);
				newm->chain = NewChain();
				newm->chain->residue = r;
				DeleteResidue(r, &c->residue);
				r->next = NULL;
				r->prev = NULL;
				r->chain = newm->chain;
				
				/* get bonds which belong to this residue */
				for(bond=m->bond;bond;bond=nextbond) {
					nextbond = bond->next;
					if (bond->atom1->residue == r && bond->atom2->residue == r) {
						DeleteBond(bond, &m->bond);
						bond->prev = bond->next = NULL;
						EnterBond(bond, &newm->bond);
					}
				}
				
				GetNeighbor(newm);
			}
		}
	}

	FreeMol(mlist);
	*mlist = newmlist;
}

static void	get_arguments (int *ifmt, int *ofmt, int *tfmt,
	char *in_file, char *out_file, char *param_file, char *tpl_file,
	FILE **in_fp, FILE **out_fp, FILE **param_fp, int *flags)
{
	int	i, n;
	char	str[256], zopt_str[256];
	FILE	*fp;
	ZOption	*zopt;

	strcpy(str, file_format[0].text);
	for(i=1;i<FILE_MAX;i++) {
		sprintf(str, "%s, %s", str, file_format[i].text);
	}
	fprintf(stdout, "Valid input file formats: %s\n", str);
	fprintf(stdout, "Choose an input file format (%s): ", file_format[DEFAULT_IFMT].text);
	fgets(str, sizeof(str), stdin);
	TrimStr(str, " \t\n");
	if (strlen(str) == 0) strcpy(str, file_format[DEFAULT_IFMT].text);
	*ifmt = get_file_type(str);
	if (*ifmt == -1) {
		*ifmt = DEFAULT_IFMT;	/* default */

		PRINT("Unrecognized input file format: %s\n", str);
		PRINT("Using the default input file format %s (%s)...\n", 
			file_format[*ifmt].text, file_format[*ifmt].desc);
		PRINT("\n");
	}

	fprintf(stdout, "Choose an output file format (%s): ", file_format[DEFAULT_OFMT].text);
	fgets(str, sizeof(str), stdin);
	TrimStr(str, " \t\n");
	if (strlen(str) == 0) strcpy(str, file_format[DEFAULT_OFMT].text);
	*ofmt = get_file_type(str);
	if (*ofmt == -1) {
		*ofmt = DEFAULT_OFMT;	/* default */

		PRINT("Unrecognized output file format: %s\n", str);
		PRINT("Using the default output file format %s (%s)...\n", 
			file_format[*ofmt].text, file_format[*ofmt].desc);
		PRINT("\n");
	}

	do {
		in_file[0] = '\0';
		fprintf(stdout, "Enter input file name: ");
		fgets(str, sizeof(str), stdin);
		TrimStr(str, " \t\n");
		if (strlen(str) == 0) continue;

		strcpy(in_file, str);
		if (strlen(in_file) > 0) {
			if (!(fp = fopen(in_file, "r"))) {
				PRINT("Can't open file %s.\n", in_file);
			} else {
				*in_fp = fp;
				break;
			}
		}
	} while (1);

	if (NEED_PARAM_FILE(*ifmt)) {
		n = 0;
		do {
			param_file[0] = '\0';
			fprintf(stdout, "Enter parameter file name: ");
			fgets(str, sizeof(str), stdin);
			TrimStr(str, " \t\n");
			if (strlen(str) == 0) {
				n++;
				if (n == 3 && (*ifmt == FILE_BOSS_ZMAT || *ifmt == FILE_MCPRO_ZMAT)) {
					PRINT("Hmm... You don't seem to have a parameter file.\n");
					PRINT("Z-matrix will be treated as QM BOSS/MCPRO format.\n");
					if (*ifmt == FILE_BOSS_ZMAT) *ifmt = FILE_BOSS_QM_ZMAT;
					if (*ifmt == FILE_MCPRO_ZMAT) *ifmt = FILE_MCPRO_QM_ZMAT;
					break;
				}
				continue;
			}

			n = 0;
			strcpy(param_file, str);
			if (strlen(param_file) > 0) {
				if (!(fp = fopen(param_file, "r"))) {
					PRINT("Can't open file %s.\n", param_file);
				} else {
					*param_fp = fp;
					break;
				}
			}
		} while (1);
	}

	do {
		out_file[0] = '\0';
		fprintf(stdout, "Enter output file name (stdout): ");
		fgets(str, sizeof(str), stdin);
		TrimStr(str, " \t\n");
		if (strlen(str) == 0) break;
		strcpy(out_file, str);
		if (!(fp = fopen(out_file, "r"))) break;
		fclose(fp);
		if (!(fp = fopen(out_file, "w"))) {
			PRINT("Can't create %s for writing.\n", out_file);
			continue;
		}
		fclose(fp);

		PRINT("File %s already exists.\n", out_file);
		PRINT("Do you want to overwrite it (y/n)? (n): ");
		fgets(str, sizeof(str), stdin);
		TrimStr(str, " \t\n");
		if (strcasecmp(str, "y") == 0) break;
	} while (1);

	if (strlen(out_file) > 0) {
		if (!(fp = fopen(out_file, "w"))) {
			PRINT("Can't open file %s.\n", out_file);
		} else *out_fp = fp;
	}

	*flags = 0;
	zopt_str[0] = '\0';
	if (IS_BOSSZMAT_FORMAT(*ofmt)) {
		fprintf(stdout, "Z-matrix options may be one or a combination of\n");
		for(i=0,zopt=zoption_tab;i<ZOPTION_SIZE;i++,zopt++) {
			fprintf(stdout, "  %-13s  %s\n", zopt->text, zopt->desc);
		}
		fprintf(stdout, "Enter Z-matrix options if necessary: ");
		fgets(str, sizeof(str), stdin);
		TrimStr(str, " \t\n");
		if (strlen(str) > 0) {
			parse_bosszmat_flags(str, flags);
			strcpy(zopt_str, str);
		}
	}

	if (IS_ZMAT_FORMAT(*ofmt)) {
		fprintf(stdout, "Do you want to use a template Z-matrix (y/n)? (n): ");
		fgets(str, sizeof(str), stdin);
		TrimStr(str, " \t\n");
		if (strcasecmp(str, "y") == 0) {
			/* get the file name */
			n = 0;
			do {
				n++;
				tpl_file[0] = '\0';
				fprintf(stdout, "Enter template Z-matrix file name: ");
				fgets(str, sizeof(str), stdin);
				TrimStr(str, " \t\n");
				strcpy(tpl_file, str);
				if (n == 2) break;
			} while (strlen(tpl_file) == 0);

			for(i=0,str[0]='\0';i<FILE_MAX;i++) {
				if (!IS_ZMAT_FORMAT(i)) continue;
				if (!str[0]) strcpy(str, file_format[i].text);
				else sprintf(str, "%s, %s", str, file_format[i].text);
			}

			fprintf(stdout, "What is the format of Z-matrix (%s)? (%s): ",
				str, file_format[DEFAULT_TFMT].text);
			fgets(str, sizeof(str), stdin);
			TrimStr(str, " \t\n");
			if (strlen(str) == 0) strcpy(str, file_format[DEFAULT_TFMT].text);
			*tfmt = get_file_type(str);
			if (*tfmt == -1) {
				*tfmt = DEFAULT_TFMT;	/* default */

				PRINT("Unrecognized template Z-matrix format: %s\n", str);
				PRINT("Using the default output file format %s (%s)...\n", 
					file_format[*tfmt].text, file_format[*tfmt].desc);
				PRINT("\n");
			}
		}
	}

	fprintf(stdout, "\n");
	fprintf(stdout, "+--------------------------------------------------------\n");
	fprintf(stdout, "| Options you have entered are equivalent to\n");
	fprintf(stdout, "| %% autozmat -i %s -o %s", file_format[*ifmt].text, file_format[*ofmt].text);
	if (param_file[0]) {
		fprintf(stdout, " -p %s", param_file);
	}
	if (zopt_str[0]) {
		fprintf(stdout, " -z %s", zopt_str);
	}
	if (out_file[0]) {
		fprintf(stdout, " -f %s", out_file);
	}
	if (tpl_file[0]) {
		fprintf(stdout, " -t %s -T %s", tpl_file, file_format[*tfmt].text);
	}
	if (in_file[0]) {
		fprintf(stdout, " %s", in_file);
	}

	fprintf(stdout, "\n");
	fprintf(stdout, "+--------------------------------------------------------\n");
	fprintf(stdout, "\n");
}

static void	report_properties (MolPtr mlist, int flags, FILE *fp)
{
	MolPtr	m;
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a;
	float	mw;
	int	n=0;

	fprintf(fp, "\n");
	ForEachMol(mlist,m) {
		n++;
		fprintf(fp, "Molecule #%d:\n", n);
		if (flags & M_REPORT_MW) {
			mw = 0;
			ForEachChainResAtom(m->chain,c,r,a) mw += GetElemWeight(a->an);
			fprintf(fp, " MW = %.3f\n", mw);
		}
	}
}

int	main (int argc, char *argv[])
{
	int	natoms, opt, flags, hydrogen_flags=0, guess_bondtype_flags, boss_user_flags=0;
	int	ifmt=DEFAULT_IFMT, ofmt=DEFAULT_OFMT, tfmt=DEFAULT_TFMT;
	char	in_file[256], out_file[256], param_file[256], tpl_file[256], cmd_file[256];
	char	*p, lib_dir[PATH_MAX], opls_dir[PATH_MAX], str[256];
	FILE	*fp, *in_fp=NULL, *param_fp=NULL, *out_fp=stdout, *tpl_fp=NULL;

	MolPtr	m, mlist=NULL, newm, newmlist;
	GaussZmatPtr	gauss=NULL;
	JaguarZmatPtr	jaguar=NULL;
	BossZmatPtr	boss=NULL;

	param_file[0] = '\0';
	tpl_file[0] = '\0';
	cmd_file[0] = '\0';
	in_file[0] = '\0';
	out_file[0] = '\0';
	opls_dir[0] = '\0';
	lib_dir[0] = '\0';

	/* parse command-line arguments */
	while ((opt = get_opt(argc, argv, "c:e:g:z:i:o:p:P:f:t:T:hHDvVIxa")) != -1) {
		switch (opt) {
		case 'c':
			PRINT("-c option has not beem implemented yet. ignored...\n");
			strcpy(cmd_file, opt_arg);
			break;
		case 'D':
			debug = 1;
			PRINT("Debugging flag is turned on.\n");
			break;
		case 'g':
			parse_autozmat_flags(opt_arg, &autozmat_flags);
			break;
		case 'i':
			ifmt = get_file_type(opt_arg);
			if (ifmt == -1) {
				ifmt = DEFAULT_IFMT;	/* default */

				PRINT("Unrecognized input file format: %s\n", opt_arg);
				PRINT("Using the default input file format %s (%s)...\n", 
					file_format[ifmt].text, file_format[ifmt].desc);
			}

if (debug) PRINT("Input Format = %s	%s\n",
file_format[ifmt].text, file_format[ifmt].desc);

			break;
		case 'o':
			ofmt = get_file_type(opt_arg);
			if (ofmt == -1) {
				ofmt = DEFAULT_OFMT;	/* default */

				PRINT("Unrecognized output file format: %s\n", opt_arg);
				PRINT("Using the default output file format %s (%s)...\n", 
					file_format[ofmt].text, file_format[ofmt].desc);
			}

if (debug) PRINT("Output Format = %s	%s\n",
file_format[ofmt].text, file_format[ofmt].desc);

			break;

		case 'e':
			if (!seterrorfile(opt_arg)) {
				printerror();
			}
			break;

		case 'f':
			strcpy(out_file, opt_arg);
			if (!(fp = fopen(out_file, "w"))) {
				PRINT("Can't open file %s.\n", out_file);
			} else out_fp = fp;
			break;

		case 'p':
			strcpy(param_file, opt_arg);
			break;

		case 'P':
			strcpy(opls_dir, opt_arg);
			break;

		case 'z':
			parse_bosszmat_flags(opt_arg, &bosszmat_flags);
			break;
		case 'H':
			hydrogen_flags = M_ADD_H;
			break;
		case 'h':
		case '?':
			print_usage(base_name(argv[0]));
			exit(opt == 'h' ? 0 : 1);
			break;
		case 'V':	/* verbose */
			autozmat_flags |= M_VERBOSE;
			break;
		case 'v':
			PRINT("%s v%s\n", base_name(argv[0]), VERSION);
			if (autozmat_flags & M_VERBOSE) {
#ifdef HAVE_ATOMTYPE
	PRINT("Atomtyping Enabled Version.\n");
#else
	PRINT("Atomtyping Disabled Version.\n");
#endif
			}
			PRINT("(C) Copyright 1998-1999 WL Jorgensen Group, Yale University\n");
			PRINT("All Rights Reserved\n");
			exit(0);
			break;

		case 'I':
			get_arguments(&ifmt, &ofmt, &tfmt, in_file, out_file, param_file, tpl_file,
				&in_fp, &out_fp, &param_fp, &flags);
			bosszmat_flags |= flags;

			break;

		case 't':
			strcpy(tpl_file, opt_arg);
			break;

		case 'T':
			tfmt = get_file_type(opt_arg);
			if (!IS_ZMAT_FORMAT(tfmt)) {
				tfmt = DEFAULT_TFMT;

				PRINT("Unrecognized Z-matrix format: %s\n", opt_arg);
				PRINT("Using the default output Z-matrix format %s (%s)...\n", 
					file_format[tfmt].text, file_format[tfmt].desc);
			}
			break;

		case 'x':	/* remove dummies before printing */
			autozmat_flags |= M_REMOVE_DUMMY;
			break;
		
		case 'a':	/* do not use atomtype library */
			autozmat_flags |= (M_NO_ATOM_TYPING|M_NO_TORSION_TYPING);
			break;
		}
	}

	if (!in_fp) {
		if (argc-opt_ind < 1) {
			in_fp = stdin;
		} else {
			strcpy(in_file, argv[opt_ind]);
			if (!(in_fp = fopen(in_file, "r"))) {
				PRINT("ERROR: Can't open file %s.\n", in_file);
				exit(1);
			}
		}
	}

	get_library_dir(argv[0], lib_dir);
	sprintf(default_opls_par_file, "%sboss%coplsaa.par", lib_dir, DIRSEP);

	if (NEED_PARAM_FILE(ifmt)) {
		if (!param_file[0] && lib_dir[0]) {
			strcpy(param_file, default_opls_par_file);

			/* check if param_file exists */
			if (!is_readable_file(param_file)) {
				param_file[0] = '\0';
			} else {
if (autozmat_flags & M_VERBOSE) {
				PRINT("*** Parameter file not specified. Option '-p param_file' may be used.\n");
				PRINT("*** Trying with default parameter file %s ...\n", param_file);
}
			}
		}

		if (!param_file[0]) {
			PRINT("*** Input format %s may need a parameter file (use '-p param_file' option).\n",
				file_format[ifmt].text);
			PRINT("*** Proceeding without a parameter file...\n");
			PRINT("*** WARNING: Atoms may appear as dummies.\n");

		} else if (!(param_fp = fopen(param_file, "r"))) {
			PRINT("ERROR: Can't open parameter file %s.\n", param_file);
			if (in_fp != stdin) fclose(in_fp);
			exit(1);
		}
	}
	
	/* set debug levels in various libraries */
	setmolerrorfunc(seterror);
	if (autozmat_flags & M_VERBOSE) {
		setmoldebuglevel(1);
		set_opls_debug_level(2);
	}
	if (autozmat_flags & M_VERBOSE_ATOMTYPE) set_atomtype_debug_level(2);
	else set_atomtype_debug_level(-1);

	boss_user_flags = 0;
	if (autozmat_flags & M_USER_CHARGE) boss_user_flags |= B_USER_CHARGE;
	if (autozmat_flags & M_USER_SIGMA) boss_user_flags |= B_USER_SIGMA;
	if (autozmat_flags & M_USER_EPSILON) boss_user_flags |= B_USER_EPSILON;

#ifdef HAVE_ATOMTYPE
	if (opls_dir[0]) {
		sprintf(opls_dir, "OPLS_DIR=%s", opls_dir);
		putenv(opls_dir);
if (autozmat_flags & M_VERBOSE) {
		PRINT("Environment variable %s.\n", opls_dir);
}
	} else if (!(p = getenv("OPLS_DIR"))) {
		/* check if BOSS_HOME is set */
		if ((p = getenv("BOSS_HOME"))) {
			sprintf(opls_dir, "OPLS_DIR=%s%clib%copls", p, DIRSEP, DIRSEP);
			putenv(opls_dir);
if (autozmat_flags & M_VERBOSE) {
			PRINT("Environment variable BOSS_HOME=%s\n", p);
			PRINT("OPLS_DIR is set to %s.\n", opls_dir);
}
		} else {

if (autozmat_flags & M_VERBOSE) {
			PRINT("Environment variable OPLS_DIR not set.\n");
}
			if (lib_dir[0]) {
				sprintf(opls_dir, "OPLS_DIR=%s%s", lib_dir, "opls");
				putenv(opls_dir);
if (autozmat_flags & M_VERBOSE) {
				PRINT(" Set to %s ...\n", opls_dir);
}
			}
		}
	}
#endif

	/***********************************************************
	 *  Give Some Starting Info
	 ***********************************************************/

if (!(autozmat_flags & M_QUIET)) {
	PRINT("Job Info: input format=%s, output format=%s\n", file_format[ifmt].text, file_format[ofmt].text);
}

	/***********************************************************
	 *  Read Input Structure
	 ***********************************************************/

	mlist = NULL;
	boss = NULL;
	gauss = NULL;
	jaguar = NULL;

	switch (ifmt) {
	case FILE_BOSS_ZMAT:
	case FILE_BOSS_QM_ZMAT:
	case FILE_BOSS_SOLVENT_ZMAT:
	case FILE_MCPRO_ZMAT:
	case FILE_MCPRO_QM_ZMAT:
		flags = 0;
		if (ifmt == FILE_MCPRO_ZMAT || ifmt == FILE_MCPRO_QM_ZMAT) flags |= B_MCPRO;
		if (ifmt == FILE_BOSS_QM_ZMAT || ifmt == FILE_MCPRO_QM_ZMAT) flags |= B_QM;
		if (ifmt == FILE_BOSS_SOLVENT_ZMAT) flags |= B_SOLVENT_ZMAT;

		if (!(boss = FLoadBossZmat(in_fp, param_fp, flags))) break;

		if (autozmat_flags & M_USER_CHARGE) boss->flags |= B_USER_CHARGE;
		if (autozmat_flags & M_USER_SIGMA) boss->flags |= B_USER_SIGMA;
		if (autozmat_flags & M_USER_EPSILON) boss->flags |= B_USER_EPSILON;

		if (!(mlist = CreateMolFromBossZmat(boss))) {
			FreeBossZmat(&boss);
			break;
		}
		break;

	case FILE_GAUSS_ZMAT:
	case FILE_MOPAC_ZMAT:
		if (!(gauss = (ifmt == FILE_GAUSS_ZMAT) ?
			FLoadGaussZmat(in_fp) :
			FLoadMopacZmat(in_fp))) break;
		if (!(mlist = CreateMolFromZmat(gauss->zmatlist))) {
			FreeGaussZmat(&gauss);
			break;
		}
		break;

	case FILE_JAGUAR_ZMAT:
		if (!(jaguar = FLoadJaguarZmat(in_fp))) break;
		if (!(mlist = CreateMolFromZmat(jaguar->zmatlist))) {
			FreeJaguarZmat(&jaguar);
			break;
		}
		break;

	case FILE_CARTESIAN:
		mlist = FLoadCartesian(in_fp);
		break;

	case FILE_PDB:
		mlist = FLoadPdbFile(in_fp);
		break;

	case FILE_MMODEL:
		mlist = FLoadMModel(in_fp);
		break;

	case FILE_MDL:
		mlist = FLoadMdlSd(in_fp);
		break;

	case FILE_MINDTOOL:
		mlist = FLoadMindtool(in_fp, (autozmat_flags & M_BREAK_MIND) ? MINDTOOL_BREAK_PERCENT : 0);
		if (autozmat_flags & M_BREAK_MIND) break_residue_to_molecule(&mlist);
		break;

	case FILE_PEPZ:
		break;

	case FILE_TRIPOS_MOL:
		mlist = FLoadTriposMol(in_fp);
		break;

	case FILE_TRIPOS_MOL2:
		mlist = FLoadTriposMol2(in_fp);
		break;
	}

	if (!mlist) {
		if (isseterror()) printerror();
		PRINT("Failed to read %s.\n", in_file[0] ? in_file : "from stdin");
		exit(1);
	}

if (autozmat_flags & M_VERBOSE) {
	PRINT("Initial Coordinates:\n");
	ForEachMol(mlist,m) FPrintMol(stderr, m);
}

	/* Initializes molecule info. Also adds hydrogens if -H option is
	 * specified.
	 */
	if (ifmt == FILE_MMODEL || ifmt == FILE_MDL ||
	    ifmt == FILE_TRIPOS_MOL || ifmt == FILE_TRIPOS_MOL2) {
		guess_bondtype_flags = 0;
	} else {
		guess_bondtype_flags = 1;
	}
	natoms = CountAtomInChains(mlist->chain);
	init_mol_info(mlist, in_file, hydrogen_flags, guess_bondtype_flags);
	if (natoms != CountAtomInChains(mlist->chain)) {
		/* Number of atoms may have been increased by adding hydrogens.
		Z-matrix is not valid any more */
		if (gauss) FreeGaussZmat(&gauss);
		if (jaguar) FreeJaguarZmat(&jaguar);
		if (boss) FreeBossZmat(&boss);
	}

	/* apply template z-matrix */
	if (IS_ZMAT_FORMAT(ofmt) && (tpl_fp = fopen(tpl_file, "r"))) {
		BossZmatPtr	tmp_boss;
		GaussZmatPtr	tmp_gauss;
		JaguarZmatPtr	tmp_jaguar;

if (autozmat_flags & M_VERBOSE) {
	PRINT("Applying template Z-matrix %s to %s\n", tpl_file, in_file);
}

		switch (tfmt) {
		case FILE_BOSS_ZMAT:
		case FILE_BOSS_QM_ZMAT:
		case FILE_BOSS_SOLVENT_ZMAT:
		case FILE_MCPRO_ZMAT:
		case FILE_MCPRO_QM_ZMAT:
			flags = 0;
			if (tfmt == FILE_MCPRO_ZMAT || tfmt == FILE_MCPRO_QM_ZMAT) flags |= B_MCPRO;
			if (tfmt == FILE_BOSS_QM_ZMAT || tfmt == FILE_MCPRO_QM_ZMAT) flags |= B_QM;
			if (tfmt == FILE_BOSS_SOLVENT_ZMAT) flags |= B_SOLVENT_ZMAT;

			if (!(tmp_boss = FLoadBossZmat(tpl_fp, NULL, flags))) break;
			SetTemplateZmat((unsigned long)tmp_boss, Z_TYPE_BOSS);
			FreeBossZmat(&tmp_boss);
			if (!(tmp_boss = (BossZmatPtr)ApplyTemplateZmat(mlist, Z_TYPE_BOSS, boss_user_flags))) {
				printerror();
				break;
			}
			if (boss) FreeBossZmat(&boss);
			if (gauss) FreeGaussZmat(&gauss);
			if (jaguar) FreeJaguarZmat(&jaguar);
			boss = tmp_boss;
			if (autozmat_flags & M_USER_CHARGE) boss->flags |= B_USER_CHARGE;
			if (autozmat_flags & M_USER_SIGMA) boss->flags |= B_USER_SIGMA;
			if (autozmat_flags & M_USER_EPSILON) boss->flags |= B_USER_EPSILON;

			break;

		case FILE_GAUSS_ZMAT:
		case FILE_MOPAC_ZMAT:
			if (!(tmp_gauss = (tfmt == FILE_GAUSS_ZMAT) ?
				FLoadGaussZmat(tpl_fp) :
				FLoadMopacZmat(tpl_fp))) break;
			SetTemplateZmat((unsigned long)tmp_gauss, Z_TYPE_GAUSS);
			FreeGaussZmat(&tmp_gauss);
			if (!(tmp_gauss = (GaussZmatPtr)ApplyTemplateZmat(mlist, Z_TYPE_GAUSS, boss_user_flags))) {
				printerror();
				break;
			}
			if (boss) FreeBossZmat(&boss);
			if (gauss) FreeGaussZmat(&gauss);
			if (jaguar) FreeJaguarZmat(&jaguar);
			gauss = tmp_gauss;
			break;

		case FILE_JAGUAR_ZMAT:
			if (!(tmp_jaguar = FLoadJaguarZmat(tpl_fp))) break;
			SetTemplateZmat((unsigned long)tmp_jaguar, Z_TYPE_JAGUAR);
			FreeJaguarZmat(&tmp_jaguar);
			if (!(tmp_jaguar = (JaguarZmatPtr)ApplyTemplateZmat(mlist, Z_TYPE_JAGUAR, boss_user_flags))) {
				printerror();
				break;
			}
			if (boss) FreeBossZmat(&boss);
			if (gauss) FreeGaussZmat(&gauss);
			if (jaguar) FreeJaguarZmat(&jaguar);
			jaguar = tmp_jaguar;
			break;
		}
		fclose(tpl_fp);
	}

	clearerror();
	switch (ofmt) {
	case FILE_BOSS_ZMAT:
	case FILE_BOSS_QM_ZMAT:
	case FILE_BOSS_SOLVENT_ZMAT:
	case FILE_MCPRO_ZMAT:
	case FILE_MCPRO_QM_ZMAT:
	case FILE_PEPZ:

		if (!mlist) break;
		flags = 0;
		if (ofmt == FILE_MCPRO_ZMAT || ofmt == FILE_MCPRO_QM_ZMAT) flags |= B_MCPRO;
		if (ofmt == FILE_BOSS_QM_ZMAT || ofmt == FILE_MCPRO_QM_ZMAT) flags |= B_QM;
		if (ofmt == FILE_BOSS_SOLVENT_ZMAT) flags |= B_SOLVENT_ZMAT;

		if (!boss) {
			if (gauss && !(gauss->flags & GAUSS_CARTESIAN)) boss = CreateBossZmatFromGauss(gauss);
			else if (jaguar && !(jaguar->flags & JAGUAR_CARTESIAN)) boss = CreateBossZmatFromJaguar(jaguar);
			if (boss) {
				if (autozmat_flags & M_USER_CHARGE) boss->flags |= B_USER_CHARGE;
				if (autozmat_flags & M_USER_SIGMA) boss->flags |= B_USER_SIGMA;
				if (autozmat_flags & M_USER_EPSILON) boss->flags |= B_USER_EPSILON;
			}
		} else if (boss->flags & B_SOLVENT_ZMAT) {
			boss->flags &= ~B_SOLVENT_ZMAT;
		}

		if (boss) {
			if (ofmt == FILE_BOSS_SOLVENT_ZMAT) boss->flags |= B_SOLVENT_ZMAT;
			if (bosszmat_flags) AddAllBossZmatOptions(boss, mlist, bosszmat_flags, NULL);
			if (DO_TYPING(autozmat_flags)) {
				natoms = CountAtomInChains(mlist->chain);
				if (!assign_opls_parameters(boss) &&
					(autozmat_flags & M_ABORT_ON_ATOMTYPE_ERR)) {
					PRINT("Aborted.\n");
					exit(1);
				}
			}

			if (ofmt == FILE_PEPZ) {
				fprint_pepz(out_fp, mlist->chain, boss, boss->atomtypes, pepz_bgn_atom_type);
			} else {
				fprint_boss_zmat_with_flags(out_fp, boss, flags);
			}
			break;
		}
		ForEachMol(mlist,m) {
			boss = CreateBossZmatFromMol(m, flags,
				(autozmat_flags & M_NO_BGN_DUMMY) ? 0 : Z_ADD_BGN_DUMMIES);
			if (!boss) continue;
			if (autozmat_flags & M_USER_CHARGE) boss->flags |= B_USER_CHARGE;
			if (autozmat_flags & M_USER_SIGMA) boss->flags |= B_USER_SIGMA;
			if (autozmat_flags & M_USER_EPSILON) boss->flags |= B_USER_EPSILON;

			if (ofmt == FILE_BOSS_SOLVENT_ZMAT) boss->flags |= B_SOLVENT_ZMAT;
			if (bosszmat_flags) AddAllBossZmatOptions(boss, m, bosszmat_flags, NULL);

			if (DO_TYPING(autozmat_flags)) {
				natoms = CountAtomInChains(m->chain);
				if (!assign_opls_parameters(boss) &&
					(autozmat_flags & M_ABORT_ON_ATOMTYPE_ERR)) {
					PRINT("Aborted.\n");
					exit(1);
				}
			}

			if (ofmt == FILE_PEPZ) {
				fprint_pepz(out_fp, m->chain, boss, boss->atomtypes, pepz_bgn_atom_type);
			} else {
				fprint_boss_zmat_with_flags(out_fp, boss, flags);
			}
			FreeBossZmat(&boss);
		}
		break;

	case FILE_MOPAC_ZMAT:
	case FILE_GAUSS_ZMAT:
		if (!mlist) break;
		if (gauss) {
			if (ofmt == FILE_MOPAC_ZMAT) FPrintMopacZmat(out_fp, gauss);
			else FPrintGaussZmat(out_fp, gauss);
			break;
		} else if (
			(jaguar && (gauss = CreateGaussZmatFromJaguar(jaguar))) ||
			(boss && (gauss = CreateGaussZmatFromBoss(boss)))) {
			if (ofmt == FILE_MOPAC_ZMAT) FPrintMopacZmat(out_fp, gauss);
			else FPrintGaussZmat(out_fp, gauss);
			FreeGaussZmat(&gauss);
			break;
		}
		ForEachMol(mlist,m) {
			if (!(gauss = CreateGaussZmatFromMol(m))) continue;
			if (ofmt == FILE_MOPAC_ZMAT) FPrintMopacZmat(out_fp, gauss);
			else FPrintGaussZmat(out_fp, gauss);
			FreeGaussZmat(&gauss);
		}

		break;

	case FILE_JAGUAR_ZMAT:
		if (!mlist) break;
		if (jaguar) {
			FPrintJaguarZmat(out_fp, jaguar);
			break;
		} else if ((gauss && (jaguar = CreateJaguarZmatFromGauss(gauss)))) {
			FPrintJaguarZmat(out_fp, jaguar);
			FreeJaguarZmat(&jaguar);
			break;
		} else if ((boss && (jaguar = CreateJaguarZmatFromBoss(boss)))) {
			FPrintJaguarZmat(out_fp, jaguar);
			FreeJaguarZmat(&jaguar);
			break;
		}
		ForEachMol(mlist,m) {
			if (!(gauss = CreateGaussZmatFromMol(m))) continue;
			if (!(jaguar = CreateJaguarZmatFromGauss(gauss))) continue;
			FPrintJaguarZmat(out_fp, jaguar);
			FreeGaussZmat(&gauss);
			FreeJaguarZmat(&jaguar);
		}
		break;

	case FILE_CARTESIAN:
	case FILE_PDB:
	case FILE_MINDTOOL:
	case FILE_MMODEL:
	case FILE_MDL:
	case FILE_TRIPOS_MOL:
	case FILE_TRIPOS_MOL2:

		if (!mlist) break;
		ForEachMol(mlist,m) {
			FixAtomSerno(m);
		}

		/* remove dummies if -z option is set */
		if (autozmat_flags & M_REMOVE_DUMMY) {
			newmlist = NULL;
			ForEachMol(mlist,m) {
				MolPtr	mnext;

				mnext = m->next;
				m->next = NULL;
				newm = CopyMol(m);
				m->next = mnext;
				if (!newm) {
					FreeMol(&newmlist);
					newmlist = mlist;
					break;
				}

				RemoveDummiesFromMol(newm);
				EnterMol(newm, &newmlist);
			}
		} else newmlist = mlist;

		/* rebuild the connectivity information */
#ifdef HAVE_ATOMTYPE
if (DO_TYPING(autozmat_flags) && NEED_CONTAB(ofmt)) {
	ForEachMol(newmlist,m) {
		natoms = CountAtomInChains(m->chain);
		if (!build_connectivity(m)) {
			PRINT("Error occurred while building connectivity.\n");

			if (ofmt==FILE_PDB) {
				PRINT("Ignored. Use '-a' option to bypass this.\n");
				PRINT("Connectivity information may be incorrect.\n");
			} else {
				PRINT("Try '-a' option to bypass this error.\n");
				PRINT("Without atomtyping, atom types for BOSS/MCPRO will not be assigned\n");
				PRINT("and connectivity information may be inappropriate.\n");
				exit(1);
			}
		}
	}
}
#endif

if (autozmat_flags & M_VERBOSE) {
	PRINT("Final Coordinates:\n");
	ForEachMol(mlist,m) FPrintMol(stderr, m);
}
		switch (ofmt) {
		case FILE_MMODEL:
		case FILE_MDL:
		case FILE_PDB:
			switch (ofmt) {
			case FILE_MMODEL:
				FPrintMModel(out_fp, newmlist);
				break;
			case FILE_MDL:
				FPrintMdlSd(out_fp, newmlist);
				break;
			case FILE_PDB:
				FPrintPdbMol(out_fp, newmlist);
				break;
			}
			break;
		
		default:
			ForEachMol(newmlist,m) {
				switch (ofmt) {
				case FILE_TRIPOS_MOL:
					FPrintTriposMol(out_fp, m);
					break;
				case FILE_TRIPOS_MOL2:
					FPrintTriposMol2(out_fp, m);
					break;
				case FILE_CARTESIAN:
					FPrintCartesian(out_fp, m->chain);
					break;
				case FILE_MINDTOOL:
					if (m != newmlist) {
						fprintf(out_fp, "----\n");
					}
					FPrintMindtool(out_fp, m->chain);
					break;
				}
			}
			break;
		}

		if (newmlist != mlist) FreeMol(&newmlist);
		break;

	default:
		PRINT("Output in %s format not implemented yet.\n",
			file_format[ofmt].text);
		break;
	}

	if (isseterror()) printerror(); 

	/* if 'report_' options, e.g., 'report_mw', are specified
	 * do it here.
	 */
	if (autozmat_flags & M_REPORT) report_properties(mlist, (autozmat_flags & M_REPORT), out_fp);

	/* close files --- this isn't really necessary
	 * because the program is about to exit anyway...
	 */
	if (in_fp && in_fp != stdin) fclose(in_fp);
	if (out_fp && (out_fp != stdout && out_fp != stderr)) fclose(out_fp);
	if (param_fp && param_fp != stdin) fclose(param_fp);

	exit(0);
}

