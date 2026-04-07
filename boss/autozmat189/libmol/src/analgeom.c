#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "macros.h"
#include "molecule.h"

/*******************************************
 ***** analytical geometry with atoms ******
 *******************************************/

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define EPSILON	1.0E-6

/* vector cross product of v1 and v2 is returned in v3 */
static void	fvcross (float *v1, float *v2, float *v3)
{
	float	v1x,v1y,v1z,v2x,v2y,v2z;

	v1x = v1[0]; v1y = v1[1]; v1z = v1[2];
	v2x = v2[0]; v2y = v2[1]; v2z = v2[2];

	v3[0] = v1y*v2z - v1z*v2y;
	v3[1] = v1z*v2x - v1x*v2z;
	v3[2] = v1x*v2y - v1y*v2x;
}

double	GetAtomDistance (AtomPtr atom1, AtomPtr atom2)
{
	double	v1[3], v2[3];
	double	dist;

	if (!atom1 || !atom2) return 0.0;

	v1[0] = atom1->x; v1[1] = atom1->y; v1[2] = atom1->z;
	v2[0] = atom2->x; v2[1] = atom2->y; v2[2] = atom2->z;

	dist = (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]);
	dist = sqrt(dist);
	return dist;
}

double	GetAtomSQDist (AtomPtr atom1, AtomPtr atom2)
{
	double	v1[3], v2[3];

	if (!atom1 || !atom2) return 0.0;

	v1[0] = atom1->x; v1[1] = atom1->y; v1[2] = atom1->z;
	v2[0] = atom2->x; v2[1] = atom2->y; v2[2] = atom2->z;
	return SQDIST(v1, v2);
}

/* determine the angle made by three points */

double	GetAngle (double *v1, double *v2, double *v3)
{
	double	dist1, dist2, dist3;
	double	cosine;

	dist1 = DIST(v1, v2);
	dist2 = DIST(v2, v3);
	dist3 = DIST(v3, v1);

	/* law of cosine */
	cosine = (dist1*dist1 + dist2*dist2 - dist3*dist3)/(2*dist1*dist2);
	cosine = MIN(1.0, ABS(cosine)) * SIGN(cosine);
	return DEGREE(acos(cosine));
}

double	GetAtomAngle (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3)
{
	double	dist1, dist2, dist3;
	double	cosine;

	if (!atom1 || !atom2 || !atom3) return 0.0;

	dist1 = GetAtomDistance(atom1, atom2);
	dist2 = GetAtomDistance(atom2, atom3);
	dist3 = GetAtomDistance(atom3, atom1);

	/* law of cosine */
	cosine = (dist1*dist1 + dist2*dist2 - dist3*dist3)/(2*dist1*dist2);
	cosine = MIN(1.0, ABS(cosine)) * SIGN(cosine);
	return DEGREE(acos(cosine));
}

/* get the dihderal angle formed by four points */
/* the following code is a modified C version of the function "dihedral"
in calc.f in mindtool written by Jim Blake and Julian Tirado-Rives */

double	GetDihedral (double *v1, double *v2, double *v3, double *v4)
{
	double	xa, xc, xd, ya, yc, yd, za, zc, zd;
	double	labc, mabc, nabc, lbcd, mbcd, nbcd;
	double	mgnabc, mgnbcd, cosine, sign, dihedral;
	
	if (!v1 || !v2 || !v3 || !v4) return 0.0;

	xa = v1[0] - v2[0]; xc = v3[0] - v2[0]; xd = v4[0] - v2[0];
	ya = v1[1] - v2[1]; yc = v3[1] - v2[1]; yd = v4[1] - v2[1];
	za = v1[2] - v2[2]; zc = v3[2] - v2[2]; zd = v4[2] - v2[2];

	labc = (ya*zc) - (yc*za); mabc = (za*xc) - (zc*xa); nabc = (xa*yc) - (xc*ya);
	lbcd = (yd*zc) - (yc*zd); mbcd = (zd*xc) - (zc*xd); nbcd = (xd*yc) - (xc*yd);

	mgnabc = sqrt(labc*labc + mabc*mabc + nabc*nabc);
	mgnbcd = sqrt(lbcd*lbcd + mbcd*mbcd + nbcd*nbcd);

	cosine = (labc*lbcd + mabc*mbcd + nabc*nbcd)/(mgnabc*mgnbcd);

	if (ABS(cosine) > 1.0) cosine = SIGN(cosine);

	/* the variable "sign" is the sign, it can have values of -1, 0, +1 */
	sign = (xd-xc)*labc + (yd-yc)*mabc + (zd-zc)*nabc;
	sign = (ABS(sign) <= EPSILON) ? 1.0 : -1.0*SIGN(sign);

	dihedral = (cosine == 0.0) ? sign*M_PI/2.0 : sign*acos(cosine);
	return dihedral*RAD_TO_DEG;
}

double	GetAtomDihedral (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, AtomPtr atom4)
{
	double	v1[3], v2[3], v3[3], v4[3];

	if (!atom1 || !atom2 || !atom3 || !atom4) return 0.0;

	v1[0] = atom1->x; v1[1] = atom1->y; v1[2] = atom1->z;
	v2[0] = atom2->x; v2[1] = atom2->y; v2[2] = atom2->z;
	v3[0] = atom3->x; v3[1] = atom3->y; v3[2] = atom3->z;
	v4[0] = atom4->x; v4[1] = atom4->y; v4[2] = atom4->z;

	return GetDihedral(v1, v2, v3, v4);
}

/* make vectors out of three atoms. */
void	MakeAtomVectors (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *v1, float *v2)
{
	if (!atom1 || !atom2 || !atom3) return;

	v1[0] = atom2->x - atom1->x;
	v1[1] = atom2->y - atom1->y;
	v1[2] = atom2->z - atom1->z;

	v2[0] = atom3->x - atom1->x;
	v2[1] = atom3->y - atom1->y;
	v2[2] = atom3->z - atom1->z;
}

/* calculate the normal at point atom2 to vector va and vb which
are (atom2-atom1) and (atom3-atom2). Right-handed rule is applied.
'norm' is returned normalized.
*/

int	GetAtomNormal (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *norm)
{
	float	va[3], vb[3], len;

	if (!atom1 || !atom2 || !atom3 || !norm) return 0;
	norm[0] = norm[1] = 0; norm[2] = 1;	/* default norm = +z-axis */

	va[0] = atom2->x - atom1->x;
	va[1] = atom2->y - atom1->y;
	va[2] = atom2->z - atom1->z;
	vb[0] = atom3->x - atom2->x;
	vb[1] = atom3->y - atom2->y;
	vb[2] = atom3->z - atom2->z;

	fvcross(va, vb, norm);
	len = sqrt((double)VDOT(norm, norm));
	if (len < 1.0E-12) return 0;
	norm[0] /= len;
	norm[1] /= len;
	norm[2] /= len;

	return 1;
}

#define MAX_TOLFAC	30.0
#define DEF_TOLFAC	0.1

int	IsAtomLinear (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon)
{
	double	angle;
	if (epsilon < 0.0 || epsilon > MAX_TOLFAC) epsilon = DEF_TOLFAC;

	angle = GetAtomAngle(atom1, atom2, atom3);
	if (ABS(angle-180.0) <= epsilon || ABS(angle) <= epsilon) return 1;
	return 0;
}

int	IsAtomCis (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon)
{
	double	angle;
	if (epsilon < 0.0 || epsilon > MAX_TOLFAC) epsilon = DEF_TOLFAC;

	angle = GetAtomAngle(atom1, atom2, atom3);
	if (ABS(angle) <= epsilon) return 1;
	return 0;
}

int	IsAtomTrans (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon)
{
	double	angle;
	if (epsilon < 0.0 || epsilon > MAX_TOLFAC) epsilon = DEF_TOLFAC;

	angle = GetAtomAngle(atom1, atom2, atom3);
	if (ABS(angle-180.0) <= epsilon) return 1;
	return 0;
}

int	IsAtomPerp (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon)
{
	double	angle;

	if (epsilon < 0.0 || epsilon > MAX_TOLFAC) epsilon = DEF_TOLFAC;
	angle = GetAtomAngle(atom1, atom2, atom3);
	if (ABS(angle-90.0) <= epsilon) return 1;
	return 0;
}

#ifdef junk
int	IsLinear (double *v1, double *v2, double *v3, float epsilon)
{
	double	angle;
	if (epsilon < 0.0 || epsilon > MAX_TOLFAC) epsilon = DEF_TOLFAC;

	angle = GetAngle(v1, v2, v3);
	if (ABS(angle-180.0) <= epsilon || ABS(angle) <= epsilon) return 1;
	return 0;
}
#endif

