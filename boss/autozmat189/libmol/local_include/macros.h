#ifndef __MACROS_H_
#define __MACROS_H_

/* macros */
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))
#define MIN(x, y)	(((x) < (y)) ? (x) : (y))
#define ABS(x)		(((x) > 0) ? (x) : -(x))
#define ROUND(x)	(int)(((x)<0.0)?(x)-.5:(x)+.5)
#define SQUARE(x)	(x)*(x)
#define DIST(x, y)	sqrt(SQUARE(x[0]-y[0])+SQUARE(x[1]-y[1])+SQUARE(x[2]-y[2]))
#define DIST2(x, y)	sqrt(SQUARE(x[0]-y[0])+SQUARE(x[1]-y[1]))
#define SWAP(x, y)	{double tmp = (x); (x) = (y); (y) = tmp;}
#define DSWAP(x, y)	{double tmp = (x); (x) = (y); (y) = tmp;}
#define FSWAP(x, y)	{float tmp = (x); (x) = (y); (y) = tmp;}
#define ISWAP(x, y)	{int tmp = (x); (x) = (y); (y) = tmp;}
#define SIGN(x)		((x)<0.0 ? (-1.0) : (1.0))
#define SQDIST(x, y)	(SQUARE(x[0]-y[0])+SQUARE(x[1]-y[1])+SQUARE(x[2]-y[2]))
#define SQDIST2(x, y)	(SQUARE(x[0]-y[0])+SQUARE(x[1]-y[1]))
#define DEG_TO_RAD	 0.01745329251994329576
#define RAD_TO_DEG	57.29577951308232087684
#define RADIAN(x)	((x) * DEG_TO_RAD)
#define DEGREE(x)	((x) * RAD_TO_DEG)
#define READLINE(fp, str)    fgets((str), sizeof((str)), (fp))
#define VDOT(v1, v2)	(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define ANG_TO_DEBYE	4.802813198
#define DEBYE(x)	((x) * ANG_TO_DEBYE)

/* memory stuff */
#include <string.h>
#define BCOPY(s1, s2, n)	memmove((char *)s2, (char *)s1, n)
#define BZERO(s, n)		memset((char *)s, (int)'\0', n)

#define BLANKLINE(s)	(strspn(s," \t\n") == strlen(s))

/* test if (x,y) lies within a rectangle made of (x1,y1) and (x2,y2) */
#define INSIDEBOX(x,y,x1,y1,x2,y2)	(x >= x1 && y >= y1 && x <= x2 && y <= y2)
#define ATOMSQDIST(a1,a2)	(SQUARE(a1->x-a2->x)+SQUARE(a1->y-a2->y)+SQUARE(a1->z-a2->z))
#define ATOMDIST(a1,a2)		sqrt(SQUARE(a1->x-a2->x)+SQUARE(a1->y-a2->y)+SQUARE(a1->z-a2->z))

#endif
