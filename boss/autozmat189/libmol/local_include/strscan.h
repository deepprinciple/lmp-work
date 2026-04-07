#ifndef __STRSCAN_H_
#define __STRSCAN_H_

#define DATA_INT	'i'
#define DATA_FLOAT	'f'
#define DATA_STR	's'

extern char	*strscan (char *buf, char *fmt, ...);
extern int	filescan (FILE *fp, char *fmt, ...);
extern int	IsAlphabet (char *str);
extern int	IsAlphaNumeric (char *str, int sign_flag);
extern int	IsInt (char *thestr);
extern int	IsNumber (char *thestr);
extern int	IsFloat (char *str);
extern int	GetDataType (char *thestr);
#endif
