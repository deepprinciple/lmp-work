#ifndef __MOL_ERROR_H_
#define __MOL_ERROR_H_

/* error definitions */
enum {MERR_DEFAULT,
	MERR_NOFILE,
	MERR_MEM,
	MERR_EOF,
	MERR_FORMAT,

	MERR_INVALID_FLOAT,
	MERR_INVALID_NUMBER,
	MERR_UNDEFINED_CENTER,
	MERR_UNKNOWN_ATOMSYMBOL,
	MERR_MULTI_DEF_CENTER,
	MERR_ILLEGAL_CENTERNAME,

	MERR_INVALID_VAR_ASSIGN,
	MERR_VAR_NOT_FOUND,
	MERR_PARFORMAT,
	MERR_OPENFILE,
	MERR_RANGE,

	MERR_INTERNAL_COORD,
	MERR_INVALID_CENTER,
	MERR_INVALID_BONDLEN,
	MERR_INVALID_BONDANG,
	MERR_MAX};

extern int	merrno;

extern void	setmolerrorfunc (int (*func)(char *, ...));
extern int	setmolerror (char *format, ...);
extern int	catmolerror (char *format, ...);
extern void	clearmolerror (void);
extern void	setmolerrortext (char *format, ...);
extern int	setmolerrorno (int ierr);
extern int	getmolerrorno (void);
extern void	printmolerror (void);
extern void	setmoldebuglevel (int);
extern int	getmoldebuglevel (void);

#endif
