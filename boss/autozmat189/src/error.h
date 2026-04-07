#ifndef _ERROR_H_
#define _ERROR_H_

#define PRINT	vprint_info

extern int	seterror (char *format, ...);
extern int	caterror (char *format, ...);
extern int	seterrorfunc (char *funcname);
extern char	*copyerror ();
extern char	*copyerrorfunc ();
extern int	clearerror ();
extern int	isseterror ();
extern int	printerror ();

extern void	vprint_info (char *format, ...);
extern void	print_info (char *str);
extern int	seterrorfile (char *path);

#endif
