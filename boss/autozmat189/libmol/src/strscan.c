#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "strscan.h"

/*
This is a bogus subroutine to emulate the fortran read statement.
This routine is non-destructive ("buf" remains intact).
The default field width is 1. Flag is similar to scanf.
 d	int (default = 0)
 ld	long int (default = 0)
 f	float (default = 0.0)
 lf	double (default = 0.0)
 (space) skip (similar to * in scanf)
 s	string with leading and trailing spaces stripped (default = NULL)
 c	chars as is (default = spaces)
return value: pointer to string to read next
*/
static int	is_blank (char *str, int width);

char	*strscan (char *buf, char *fmt, ...)
{
	va_list		ap;
	register int	i, width, nsuccess, blank;
	char		*p, *s, *t;
	char		c;

	if (!buf) return 0;
	nsuccess = 0;

	va_start(ap, fmt);

	for (; *fmt; fmt++) {
		if (*fmt != '%') {if (*buf) buf++; continue;}

		/* field width */
		for (++fmt,width=0;isdigit(*fmt);fmt++) width = width * 10 + *fmt - '0';
		width = width == 0 ? 1 : width;

		switch (*fmt) {
		case 'd':	/* int */
			s = buf;
			for(i=0;i<width && *buf;i++) buf++;
			c = *buf;
			*buf = '\0';
			blank = is_blank(s, width);
			*(va_arg(ap, int *)) = blank ? 0 : (int) strtol(s, &p, 10);
			if (blank || (!blank && p != s)) nsuccess++;
			*buf = c;
			break;

		case 'f':	/* float */
			s = buf;
			for(i=0;i<width && *buf;i++) buf++;
			c = *buf;
			*buf = '\0';
			blank = is_blank(s, width);
			*(va_arg(ap, float *)) = blank ? 0.0 : (float)strtod(s, &p);
			if (blank || (!blank && p != s)) nsuccess++;
			*buf = c;
			break;

		case 'l':	/* long */
			fmt++;
			if (!*fmt) break;
			switch (*fmt) {
			case 'd':	/* long integer */
				s = buf;
				for(i=0;i<width && *buf;i++) buf++;
				c = *buf;
				*buf = '\0';
				blank = is_blank(s, width);
				*(va_arg(ap, long *)) = blank ? 0 : (long) strtol(s, &p, 10);
				if (blank || (!blank && p != s)) nsuccess++;
				*buf = c;
				break;
			case 'f':	/* long double */
				s = buf;
				for(i=0;i<width && *buf;i++) buf++;
				c = *buf;
				*buf = '\0';
				blank = is_blank(s, width);
				*(va_arg(ap, double *)) = blank ? (double)0.0 : (double) strtod(s, &p);
				if (blank || (!blank && p != s)) nsuccess++;
				*buf = c;
				break;
			}
			break;

		case 's':	/* string */
			if (!*buf) width = 0;
			s = t = va_arg(ap, char *);
			if (!s) break;
			for(i=0;i<width && *buf;i++) *s++ = *buf++;
			*s = '\0';
			/* remove trailing spaces */
			while (s > t && isspace(*--s)) *s = '\0';
			if (width > 0) nsuccess++;
			break;

		case 'c':	/* chars as is */
			s = va_arg(ap, char *);
			if (!s) break;
			for(i=0;i<width;i++) s[i] = ' ';	/* default */
			for(i=0;i<width && *buf;i++) *s++ = *buf++;
			*s = '\0';
			nsuccess++;
			break;

		case ' ':	/* skip */
			for(i=0;i<width && *buf;i++) buf++;
			break;

		default:
			break;
		}
	}

	va_end(ap);

	return buf;
}

/*
This is a bogus subroutine to emulate the fortran read statement.
The default field width is 1. Flag is similar to scanf.
 d	int (default = 0)
 ld	long int (default = 0)
 f	float (default = 0.0)
 lf	double (default = 0.0)
 (space) skip (similar to * in scanf)
 s	string with leading and trailing spaces stripped (default = NULL)
 c	chars as is (default = spaces)
return value: # of successfully read variables
*/

int	filescan (FILE *fp, char *fmt, ...)
{
	va_list		ap;
	register int	i, width, nsuccess, eof, blank;
	char		*p, *s, *t;
	char		c;
	char		buf[1024];	/* maximum width = 1024 */

	if (!fp) return 0;
	nsuccess = 0;

	va_start(ap, fmt);

	for (; *fmt; fmt++) {
		eof = 0;
		if (*fmt != '%') {
			if ((int)fgetc(fp) == EOF) {eof = 1; break;} else continue;
		}

		/* field width */
		for (++fmt,width=0;isdigit(*fmt);fmt++) width = width * 10 + *fmt - '0';
		width = width == 0 ? 1 : width;

		switch (*fmt) {
		case 'd':	/* int */
			s = buf;
			for(i=0;i<width;i++) {
				c = fgetc(fp);
				if ((int)c == EOF) {eof = 1; break;}
				*s++ = c;
			}
			*s = '\0';
			if (eof) break;
			blank = is_blank(s, width);
			*(va_arg(ap, int *)) = blank ? 0 : (int) strtol(s, &p, 10);
			if (blank || (!blank && p != s)) nsuccess++;
			break;

		case 'f':	/* float */
			s = buf;
			for(i=0;i<width;i++) {
				c = fgetc(fp);
				if ((int)c == EOF) {eof = 1; break;}
				*s++ = c;
			}
			*s = '\0';
			if (eof) break;
			blank = is_blank(s, width);
			*(va_arg(ap, float *)) = blank ? 0.0 : (float)strtod(s, &p);
			if (blank || (!blank && p != s)) nsuccess++;
			break;

		case 'l':	/* long */
			fmt++;
			if (!*fmt) break;
			if (*fmt == 'd' || *fmt == 'f') {
				s = buf;
				for(i=0;i<width;i++) {
					c = fgetc(fp);
					if ((int)c == EOF) {eof = 1; break;}
					*s++ = c;
				}
				*s = '\0';
				if (eof) break;
			}
			blank = is_blank(s, width);

			switch (*fmt) {
			case 'd':	/* long integer */
				*(va_arg(ap, long *)) = blank ? 0 : (long) strtol(s, &p, 10);
				if (blank || (!blank && p != s)) nsuccess++;
				break;

			case 'f':	/* long double */
				*(va_arg(ap, double *)) = blank ? (double)0.0 : (double) strtod(s, &p);
				if (blank || (!blank && p != s)) nsuccess++;
				break;
			}
			break;

		case 's':	/* string */
			s = t = va_arg(ap, char *);
			if (!s) break;

			for(i=0;i<width;i++) {
				c = fgetc(fp);
				if ((int)c == EOF) {eof = 1; break;}
				*s++ = c;
			}
			*s = '\0';
			if (eof) break;

			/* remove trailing spaces */
			while (s > t && isspace(*--s)) *s = '\0';
			if (width > 0) nsuccess++;
			break;

		case 'c':	/* chars as is */
			s = va_arg(ap, char *);
			if (!s) break;
			for(i=0;i<width;i++) s[i] = ' ';	/* default */
			for(i=0;i<width;i++) {
				c = fgetc(fp);
				if ((int)c == EOF) {eof = 1; break;}
				*s++ = c;
			}
			*s = '\0';
			if (eof) break;
			nsuccess++;
			break;

		case ' ':	/* skip */
			for(i=0;i<width;i++) {
				c = fgetc(fp);
				if ((int)c == EOF) {eof = 1; break;}
			}
			if (eof) break;
			break;

		default:
			break;
		}
		if (eof) break;
	}

	va_end(ap);

	return nsuccess;
}

static int	is_blank (char *str, int width)
{
	register int	i;
	register char	*p;

	for(i=0,p=str;i<width && *p;i++,p++) if (*p != ' ') return 0;
	return 1;
}

/* true if a given string is a decimal number */
static int	_isinteger (char *str)
{
	register char	*p;
	register short	i;
	int	center;

	if (!str) return 0;

	for(i=0,p=str;i<strlen(str);i++,p++) {
		if (i == 0 && (*p == '+' || *p == '-')) continue;
		if (!isdigit(*p)) return 0;
	}
	if (sscanf(str, "%d", &center) != 1) return 0;

	return 1;
}

/* test if a given string is a valid numral expression */

static int	_isnumber (char *str)
{
	char	*ptr;
	double	x;
	
	if (!str) return 0;
	ptr = NULL;

	/* if ptr is not NULL, the string can not be a valid floating-point number
	expression. */

	/* strtod decomposes the input string into three parts: 1) white spaces
	2) subject string ([+-][0-9][.][0-9][Ee[+-]][0-9]) 3) final string - 
	remainder of the input string. if final string does not have zero-length, the input
	string is not a valid floating-point or integer. sscanf won't be 1 if
	str is blank. */

	strtod(str, &ptr);

	/* this doesn't work!
	if (ptr || sscanf(str, "%lf", &x) != 1) return 0;
	*/

	if (strlen(ptr) > 0 || sscanf(str, "%lf", &x) != 1) return 0;

	return 1;
}

/* floating-point number
[+-][0-9][.][0-9][Ee][+-][0-9]
+1E1, -.1E+0, +1.E-01, 1E
*/

static int	_isfloat (char *str)
{
	if (!_isnumber(str)) return 0;
	if (_isinteger(str)) return 0;
	/*
	if (!strstr(str, ".")) return 0;
	*/
	return 1;
}

/* test if a string consists of only alphabets */
int	IsAlphabet (char *str)
{
	register char	*p;
	register short	i;

	if (!str || !str[0]) return 0;
	for(i=0,p=str;*p && i<strlen(str);i++,p++) if (!isalpha(*p)) return 0;
	return 1;
}

/* test if a given string is alphanumeric [a-z,A-Z,0-9]
if sign_flag = 1, '+' or '-' sign in front of the string is ignored */

int	IsAlphaNumeric (char *str, int sign_flag)
{
	register char	*p;
	register short	i;

	if (!str) return 0;
	if (sign_flag != 1) sign_flag = 0;
	if (sign_flag) if (str[0] == '+' || str[0] == '-') str++;
	if (!isalpha(str[0])) return 0;
	for(i=1,p=str+1;i<strlen(str);i++,p++) if (!isalpha(*p) && !isdigit(*p)) return 0;
	return 1;
}

/* true if a given string is a decimal number */
int	IsInt (char *thestr)
{
	char	str[256];
	register char	*p, *q;
	register short	i;
	int	center;

	if (!thestr) return 0;
	strcpy(str, thestr);	/* work with backup */

	/* find the begining */
	q = str;
	q += strspn(q, " \t\n");
	if (!*q) return 0;

	/* remove trailers */
	if ((p = strpbrk(q, " \t\n"))) *p = '\0';

	for(i=0,p=q;*p;i++, p++) {
		if (i == 0 && (*p == '+' || *p == '-')) continue;
		if (!isdigit(*p)) return 0;
	}
	if (sscanf(str, "%d", &center) != 1) return 0;

	return 1;
}

/* test if a given string is a valid numral expression */

int	IsNumber (char *thestr)
{
	char	str[256];
	char	*ptr, *q;
	double	x;
	
	if (!thestr) return 0;
	strcpy(str, thestr);
	q = str;

	/* find the begining */
	q += strspn(q, " \t\n");
	if (!*q) return 0;

	/* remove trailers */
	if ((ptr = strpbrk(q, " \t\n"))) *ptr = '\0';

	ptr = NULL;

	/* if ptr is not NULL, the string can not be a valid floating-point number
	expression. */

	/* strtod decomposes the input string into three parts: 1) white spaces
	2) subject string ([+-][0-9][.][0-9][Ee[+-]][0-9]) 3) final string - 
	remainder of the input string. if final string does not have zero-length, the input
	string is not a valid floating-point or integer. sscanf won't be 1 if
	str is blank. */

	strtod(q, &ptr);

	if (strlen(ptr) > 0 || sscanf(str, "%lf", &x) != 1) return 0;
	return 1;
}

/* floating-point number
[+-][0-9][.][0-9][Ee][+-][0-9]
+1E1, -.1E+0, +1.E-01, 1E
*/

int	IsFloat (char *str)
{
	if (!str) return 0;

	if (!IsNumber(str)) return 0;
	if (!strstr(str, ".")) return 0;
	return 1;
}

int	GetDataType (char *thestr)
{
	char	str[256];
	char	*p, *q;

	if (!thestr || !thestr[0]) return '\0';
	strcpy(str, thestr);
	q = str;

	/* find the beginning */
	q += strspn(q, " \t\n");
	if (!*q) return 's';

	/* remove trailers */
	if ((p = strpbrk(q, " \t\n"))) *p = '\0';

	if (_isinteger(q)) return DATA_INT;
	else if (_isfloat(q)) return DATA_FLOAT;
	else return DATA_STR;
}

