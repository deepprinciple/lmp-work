#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <malloc.h>

void	TrimHead (char *str, char *delim)
{
	char	*bgn, *p, *q;
	if (!str || !str[0]) return;
	bgn = str + strspn(str, delim);
	if (!*bgn) {	/* null string */
		str[0] = '\0';
		return;
	} else if (bgn == str) return;
	for(p=bgn,q=str;*p;p++,q++) *q = *p;
	*q = '\0';	/* mark the end */
}

void	TrimTail (char *str, char *delim)
{
	int	i, len;
	if (!str || !str[0]) return;
	len = strlen(str);
	for(i=len-1;i>=0;i--) {
		if (!strchr(delim, str[i])) return;
		else str[i] = '\0';
	}
}

void	TrimStr (char *str, char *delim)
{
	TrimHead(str, delim);
	TrimTail(str, delim);
}

int	IsBlank (char *str)
{
	char	*bgn;

	if (!str) return 0;
	if (!str[0]) return 1;

	bgn = str + strspn(str, " \t\n");
	return (!*bgn ? 1 : 0);
}

void	StringUpperCase (char *str)
{
	char	*p;
	if (!str || !*str) return;
	for(p=str;*p;p++) *p = toupper(*p);
}

void	StringLowerCase (char *str)
{
	char	*p;
	if (!str || !*str) return;
	for(p=str;*p;p++) *p = tolower(*p);
}

void	FPrintString (void *dest, char *str)
{
	FILE	*fp = (FILE *)dest;
	fprintf(fp, "%s", str);
}

void	AddtoString (char **str_sum, char *str)
{
	int	n;

	if (!str_sum || !str || !*str) return;
	n = strlen(str);
	if (!*str_sum) {
		*str_sum = (char *)malloc(n+1);
		strcpy(*str_sum, str);
	} else {
		n += strlen(*str_sum);
		*str_sum = (char *)realloc(*str_sum, n+1);
		strcat(*str_sum, str);
	}
}

/**************** keep printing to string ******************/
#define	STEPSIZE	1024
static char	*strtable = NULL;
static int	cur_len = 0;
static int	max_len = 0;

/* add a new string to a string table, "strtable".
 * dest is ignored.
 */
void	SPrintString (void *dest, char *str)
{
	register int	len;

	len = strlen(str);
	if (max_len - cur_len < len+1) {
		max_len += STEPSIZE;
		strtable = !strtable ? (char *)calloc(1, max_len) : (char *)realloc(strtable, max_len);
		if (!strtable) return;
	}
	strcpy(strtable + cur_len, str);
	cur_len += len;
}

void	BgnPrintString (void)
{
	if (strtable) free(strtable);
	strtable = NULL;
	cur_len = max_len = 0;
}

void	ContinuePrintString (char *str)
{
	if (!str) return;
	strtable = str;	/* must be reallocatable */
	cur_len = strlen(strtable);
	max_len = cur_len+1;
}

char	*EndPrintString (void)
{
	char	*str;
	
	if (!strtable || cur_len <= 0) return NULL;
	str = (char *)realloc(strtable, (cur_len+1)*sizeof(char));
	strtable = NULL;
	cur_len = max_len = 0;
	return str;
}

