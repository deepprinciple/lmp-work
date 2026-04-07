#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h>

#include "error.h"

int	mol_debug_level=0;

#define ERRORFUNC_LEN	132
#define ERRORTEXT_LEN	1024

static char	errorfunc[ERRORFUNC_LEN+1] = {'\0'};
static char	errortext[ERRORTEXT_LEN+1] = {'\0'};

int	merrno = -1;
static char	merrbuf[ERRORTEXT_LEN+1] = {'\0'};
static char	*merrtext[] = {
	"Z-MATRIX ERROR: %s",
	"File not found.",
	"Not enough memory or memory has not been allocated.",
	"End of file encountered unexpectedly.",
	"Error in the data format: %s",

	"Invalid floating-point specification: %s",
	"Invalid Number in: %s",
	"Undefined center in %s.",
	"Unknown atomic symbol in %s.",
	"Multiply defined center in %s.",
	"Illegal center name in %s.",

	"Invalid variable assignment in %s.",
	"Variable %s not found in variable or constant section.",
	"Error in the parameter format: %s",
	"Error in opening file.",
	"Error in the range: %s",

	"Error in internal coordinate specification: %s",
	"Invalid center: %s",
	"Invalid bond length: %s",
	"Invalid bond angle: %s",
	};

static int	(*molerrorfunc)(char *, ...) = NULL;
void	setmolerrorfunc (int (*func)(char *, ...))
{
	if (!func) return;
	molerrorfunc = func;
}

/* set mol library error */
int	setmolerror (char *format, ...)
{
	va_list	args;

	va_start(args, format);
	vsprintf(errortext, format, args);
	va_end(args);

	if (molerrorfunc) return (*molerrorfunc)(errortext);
	else {
		fprintf(stderr, "%s\n", errortext);
		return 0;
	}
}

int	catmolerror (char *format, ...)
{
	va_list	args;
	char	buf[1024];

	buf[0] = '\0';
	va_start(args, format);
	vsprintf(buf, format, args);
	va_end(args);

	strcat(errortext, buf);
	if (molerrorfunc) return (*molerrorfunc)(errortext);
	else {
		fprintf(stderr, "%s\n", errortext);
		return 0;
	}
}

void	clearmolerror (void)
{
	setmolerror("");
	errortext[0] = '\0';
	merrbuf[0] = '\0';
	merrno = -1;
}

void	setmolerrortext (char *format, ...)
{
	va_list	args;

	va_start(args, format);
	vsprintf(merrbuf, format, args);
	va_end(args);
}

int	setmolerrorno (int ierr)
{
	if (ierr <= 0 || ierr >= MERR_MAX) ierr = 0;
	merrno = ierr;
	return setmolerror(merrtext[ierr], merrbuf);
}

int	getmolerrorno (void)
{
	return merrno;
}

void	printmolerror (void)
{
	fprintf(stderr, "%s\n", errortext);
}

void	setmoldebuglevel (int n)
{
	mol_debug_level = n;
}

int	getmoldebuglevel (void)
{
	return mol_debug_level;
}

