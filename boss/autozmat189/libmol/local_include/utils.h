#ifndef __UTILS_H_
#define __UTILS_H_

extern void	TrimHead (char *str, char *delim);
extern void	TrimTail (char *str, char *delim);
extern void	TrimStr (char *str, char *delim);
extern void	AddtoString (char **str_sum, char *str);
extern int	IsBlank (char *str);
extern void	StringUpperCase (char *str);
extern void	StringLowerCase (char *str);

extern void	FPrintString (void *dest, char *str);
extern void	SPrintString (void *dest, char *str);
extern void	BgnPrintString (void);
extern void	ContinuePrintString (char *str);
extern char	*EndPrintString (void);

#endif
