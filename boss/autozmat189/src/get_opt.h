#ifndef _GET_OPT_H_
#define _GET_OPT_H_

extern char    *opt_arg;
extern int     opt_ind;

extern int	get_opt (int argc, char *argv[], char *optstring);
extern void	clear_get_opt (void);

#endif
