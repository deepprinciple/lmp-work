#ifndef _OPLSTYPE_H_
#define _OPLSTYPE_H_

#define MAX_DOMAIN_ATOMS	1000

extern int	assign_opls_parameters (BossZmatPtr boss);
extern int	build_connectivity (MolPtr mol);
extern void	set_idle_time (int tsec);
extern void	set_opls_debug_level (int level);
extern void	set_atomtype_debug_level (int level);

#endif
