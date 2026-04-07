#ifndef __BOSSPAR_H_
#define __BOSSPAR_H_
extern BossPar	*CreateBossPar (int flags);
extern BossPar	*CopyBossPar(BossPar *);
extern void	FreeBossPar (BossPar **bp);
extern BossParEntry	*FindBossParEntry (BossPar *bp, int param_id);
extern char	*GetBossParDefaultValue (BossPar *bp, int param_id);
extern char	*GetBossParValue (BossPar *bp, int param_id);
extern void	SetBossParValue (BossPar *bp, int param_id, char *val);

extern int	FGuessBossParVersion (FILE *fp);
extern BossPar	*FLoadBossPar (FILE *fp, int version, int flags);
extern void	FPrintBossPar (FILE *fp, BossPar *bp, int version);
extern ListPtr	FLoadBossAtomTypes (FILE *fp);
extern ListPtr	FLoadBossAtomTypesFromZmat (FILE *fp);
extern void	FPrintBossOneAtomType (FILE *fp, ListPtr atomtype);
extern void	FPrintBossAtomTypes (FILE *fp, ListPtr atomtypes, int version);
extern ListPtr	FLoadBossTorsion (FILE *fp);
extern void	FPrintBossOneTorsionType (FILE *fp, ListPtr torsiontype);
extern void	FPrintBossTorsionTypes (FILE *fp, ListPtr torsiontypes, int version);
extern ListPtr	FindBossAtomType (ListPtr atomtypelist, int type);
extern void	AssignBossPar (ListPtr zmatlist, ListPtr atomtypelist);

extern void	FPrintBossParFile (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes, int version);
extern void	FPrintBossParFileWithZmatTypes (FILE *fp, BossPar *bp, ListPtr atomtypes, ListPtr torsiontypes,
	ListPtr zmat_atomtypes, ListPtr zmat_torsiontypes, int version);

extern int	FLoadBossParFile (FILE *fp, int flags, BossParPtr *bp_ret, ListPtr *at_ret, ListPtr *tt_ret, int *version_ret);

extern OplsTypeTorsion	*FLoadOplsTypeTorsion (FILE *fp);
extern OplsTypeTorsion	*FLoadDefaultOplsTypeTorsion (FILE *fp);
extern OplsTypeTorsion	*LoadOplsTypeTorsion (char *filename);
extern OplsTypeTorsion	*LoadDefaultOplsTypeTorsion (char *filename);
extern OplsTypeTorsion	*GetDefaultOplsTypeTorsion (void);
extern void	FreeOplsTypeTorsion (OplsTypeTorsion **p);

extern ListPtr	FLoadOplsStretch (FILE *fp);
extern ListPtr	FLoadOplsBending (FILE *fp, int rewind_flag);
extern OplsStretchBend	*FLoadOplsStretchBend (FILE *fp);
extern OplsStretchBend	*FLoadDefaultOplsStretchBend (FILE *fp);
extern OplsStretchBend	*LoadOplsStretchBend (char *filename);
extern OplsStretchBend	*LoadDefaultOplsStretchBend (char *filename);
extern OplsStretchBend	*GetDefaultOplsStretchBend (void);
extern void	FreeOplsStretchBend (OplsStretchBend **p);

extern void	AssignDefaultBossParByCalcType (BossParPtr bp, int calc_type);
extern char	*GetBossCalcTypeText (int calc_type);

#endif
