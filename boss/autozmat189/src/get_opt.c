#include <stdio.h>
#include <string.h>

#include "get_opt.h"

/*********************************************************
 ********     Implementation of getopt()     *************
 *********************************************************/

char	*opt_arg;	/* Global argument pointer. */
int	opt_ind = 0;	/* Global argv index. */

static char	*scan = NULL;

void	clear_get_opt (void)
{
	scan = NULL;
	opt_ind = 0;
	opt_arg = NULL;
}

int	get_opt (int argc, char *argv[], char *optstring)
{
	register int	c;
	register char	*place;

	opt_arg = NULL;

	if (!scan || !*scan) {
		if (opt_ind == 0) opt_ind++;
		if (opt_ind >= argc) return EOF;
		place = argv[opt_ind];
		if (place[0] != '-' || place[1] == '\0') return EOF;
		opt_ind++;
		if (place[1] == '-' && place[2] == '\0') return EOF;
		scan = place+1;
	}

	c = *scan++;
	place = strchr(optstring, c);
	if (place == NULL || c == ':') {
		fprintf(stderr, "%s: unknown option %c\n", argv[0], c);
		return '?';
	}
	if (*++place == ':') {
		if (*scan != '\0') {
			opt_arg = scan, scan = NULL;
		} else {
			opt_arg = argv[opt_ind], opt_ind++;
		}
	}
	return c;
}

