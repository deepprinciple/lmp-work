#ifndef __BOSSINPUT_H_
#define __BOSSINPUT_H_

extern ListPtr	FLoadBossList (FILE *fp, char *scanfmt[], int nvar);
extern void	FreeBossZmat (BossZmatPtr *bp);
extern void	FPrintBossZmat (FILE *fp, BossZmatPtr boss);
extern void	SPrintBossZmat (char **str_sum, BossZmatPtr boss);
extern void	FPrintBossZmatTypes (FILE *fp, BossZmatPtr boss);
extern void	SPrintBossZmatTypes (char **str_sum, BossZmatPtr boss);
extern void	FPrintBossZmatWithTypes (FILE *fp, BossZmatPtr bos);
extern void	SPrintBossZmatWithTypes (char **str_sum, BossZmatPtr bos);
extern void	FixBossZmatAtomType (BossZmatPtr boss);
extern void	UpdateBossZmatTypes (BossZmatPtr boss, OplsTypeTorsion *opls);

extern BossZmatPtr	FLoadBossZmatInput (FILE *fp, int flags);
extern BossZmatPtr	FLoadBossZmat (FILE *fp, FILE *parfp, int flags);
extern BossZmatPtr	CopyBossZmat (BossZmatPtr oldboss);
extern BossZmatPtr	CreateBossZmatFromMol (MolPtr mol, int flags, int autozmat_flags);
extern BossZmatPtr	CreateBossZmatFromZmat (ZmatPtr zmatlist, int flags);
extern BossZmatPtr	CreateBossZmatFromGauss (GaussZmatPtr gauss);
extern BossZmatPtr	CreateBossZmatFromJaguar (JaguarZmatPtr jaguar);

extern MolPtr	CreateMolFromBossZmat (BossZmatPtr boss);
extern MolPtr	FLoadBossPdbInput (FILE *fp);
extern MolPtr	FLoadBossMindtoolInput (FILE *fp);

extern int	UpdateBossGeomVarWithMolCoord (BossZmatPtr boss, MolPtr mol);
extern ZmatPtr	GetZmatFromBossGeomVar (BossZmatPtr boss);
extern ListPtr	GetZmatListFromBossGeomVar (BossZmatPtr boss);
extern int	EvalBossGeomVar (BossZmatPtr boss, MolPtr mol);

extern GaussZmatPtr	CreateGaussZmatFromBoss (BossZmatPtr boss);
extern JaguarZmatPtr	CreateJaguarZmatFromBoss (BossZmatPtr boss);
extern void	SetImproperTorsionFlags (ListPtr dihedlist, ChainPtr chainlist, BondPtr bondlist);

extern int	AddAllBossZmatOptions (BossZmatPtr boss, MolPtr mol, unsigned long flags, BossZmatAddAllOpt *addall);
extern BossZmatAddAllOpt	*GetAddAllBossZmatOptions (void);

#endif
