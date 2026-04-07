#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h>

#define ERRORFUNC_LEN	132
#define ERRORTEXT_LEN	2048

static char	errorfunc[ERRORFUNC_LEN+1] = {'\0'};
static char	errortext[ERRORTEXT_LEN+1] = {'\0'};

static FILE	*error_fp = NULL;

/* set program error */
int	seterror (char *format, ...)
{
	va_list	args;

	va_start(args, format);
	vsprintf(errortext, format, args);
	va_end(args);

	return 0;
}

/* concatenate program error */
int	caterror (char *format, ...)
{
	va_list	args;
	char	tmptext[ERRORTEXT_LEN+1];
	int	n;
	char	*space;

	va_start(args, format);
	vsprintf(tmptext, format, args);
	va_end(args);

	n = strlen(errortext);
	if (strchr("\n ", errortext[n-1])) space = "";
	else space = " ";

	sprintf(errortext, "%s%s%s", errortext, space, tmptext);
	return 0;
}

/* copy the name of function causing error */
void	seterrorfunc (char *funcname)
{
	strncpy(errorfunc, funcname, ERRORFUNC_LEN);
}

static char	errortext_buf[ERRORTEXT_LEN+1];
char	*copyerror ()
{
	if (strlen(errortext) > 0) {
		strcpy(errortext_buf, errortext);
		return errortext_buf;
	} else return " ";	/* return a space rather than NULL, so that the program won't crash
				when vsprintf or other functions crash */
}

char	*copyerrorfunc ()
{
	static char	str[ERRORFUNC_LEN+1];

	if (strlen(errorfunc) > 0) {
		strcpy(str, errorfunc);
		return str;
	} else return NULL;
}

void	clearerror ()
{
	errortext[0] = '\0';
	errorfunc[0] = '\0';
}

int	isseterror ()
{
	if (errortext[0] && strlen(errortext) > 0) return 1;
	else return 0;
}

void	printerror ()
{
	if (!error_fp) error_fp = stderr;
	fprintf(error_fp, "\007");
	if (strlen(errorfunc) > 0) fprintf(error_fp, "%s: %s\n", errorfunc, errortext);
	else fprintf(error_fp, "%s\n", errortext);
}

void	vprint_info (char *format, ...)
{
	char	buf[ERRORTEXT_LEN];
	va_list	args;

	if (!error_fp) error_fp = stderr;
	buf[0] = '\0';
	va_start(args, format);
	vsprintf(buf, format, args);
	va_end(args);
	fprintf(error_fp, "%s", buf);
	fflush(error_fp);
}

void	print_info (char *str)
{
	if (!error_fp) error_fp = stderr;
	fprintf(error_fp, "%s", str);
	fflush(error_fp);
}

int	seterrorfile (char *path)
{
	FILE	*fp;
	if (!(fp = fopen(path, "w"))) {
		return seterror("Can't open error log file, %s.", path);
	}
	error_fp = fp;
	return 1;
}

