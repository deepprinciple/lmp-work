#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "macros.h"
#include "molecule.h"

#define VLENGTH(v)	sqrt((double)(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]))

#define	sin_7		0.12186934340514
#define	sin_15		0.25881904510252
#define sin_19		0.333642332
#define sin_30		0.5
#define sin_35		0.577216543
#define sin_36		0.5877852523

#define	cos_7		0.99258311810202
#define	cos_15		0.96669290471077
#define cos_19		0.942699737
#define cos_30		0.86602540378444
#define cos_35		0.816591186
#define cos_36		0.8090169944	

#define EPS	1.0E-6

/*************************************
 *****     global variables     ******
 *************************************/

/* Carbon */
double	MAX_SUM_ANG = 350.;	/* max sum of three single angles */
double	MAX_C1_ANG  = 114.;	/* max single bond angle */
double	MAX_C2_ANG  = 130.; 	/* max double bond angle */
double	MAX_CC2_LEN = 1.96;	/* max double C=C bond length */
double	MAX_CC3_LEN = 1.5625;	/* max triple C#C bond length */

/* Nitrogen */
double	MIN_AMIDE_ANG = 114.;	/* min angle for an amide bond */
double	MAX_CN2_LEN   = 1.30;	/* max bond length for a C=N bond */
double	MAX_AMIDE_LEN = 1.38;	/* max bond length for an amide bond */
double	MAX_CN3_LEN   = 1.18;	/* max bond length for C-N triple bond */

static double	MAX_CN2_LEN_SQ	 = 1.69;	/* MAX_CN2_LEN*MAX_CN2_LEN */
static double	MAX_AMIDE_LEN_SQ = 1.9044;
static double	MAX_CN3_LEN_SQ	 = 1.3924;

/* Oxygen */
double	MAX_CO2_LEN = 1.32;	/* max bond length for a C=O bond */
static double	MAX_CO2_LEN_SQ = 1.7424;

/* Sulfur */
double	MAX_CS2_LEN = 1.76;	/* max bond length for a C=S bond */
static double	MAX_CS2_LEN_SQ = 3.9076;

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

/* This routine is likely to be used when theatom has only one neighbor and
 * user wants to create new atoms based on these atoms. If nextatom is not
 * NULL, torsion is set to nextatom-nbatom-theatom-hydrogen1 = 0.0
 * Otherwise, a dummy neighbor atom is created on 'theatom' and the torsion
 * will be based on that dummy neighbor.
 *
 * On income, hydrogen must have already-allocated memory.
 * This routine simply fills x,y,z fields in hydrogen.
 *
 * If hydrogen->next is not NULL, the second hydrogen will be placed
 * in such a way that the torsion between between
 * nbatom-theatom-hydrogen1-hydrogen2 = 180.0.
 *
 * (or NULL)
 * nextatom                     hydrogen1
 *        \                   /
 *         \nbatom----theatom/
 *                         \_\
 *                    angle   \ hydrogen2
 *
 * theatom  --- the atom to which new atom(s) are to be added 
 * nbatom   --- neighbor atom to theatom 
 * nextatom --- neighbor to the nbatom
 * length   --- distance between theatom and new atoms (hydrogen1, hydrogen2)
 * angle    --- (in degree) angle between nbatom-theatom-hydrogen
 * hydrogen --- hydrogens to be added (one or two hydrogens) 
 */

static int	newatom_from_one_neighbor (AtomPtr theatom, AtomPtr nbatom, AtomPtr nextatom, float length, float angle, AtomPtr hydrogen)
{
	double	p[3];

	if (!theatom || !nbatom || !hydrogen) return 0;
	angle = RADIAN(angle);

	/* default values */
	hydrogen->x = theatom->x+length;
	hydrogen->y = theatom->y;
	hydrogen->z = theatom->z;
	if (hydrogen->next) {
		hydrogen->next->x = theatom->x;
		hydrogen->next->y = theatom->y+length;
		hydrogen->next->z = theatom->z;
	}
	if (GetAtomSQDist(theatom, nbatom) < 1.0E-12) return 0;	/* too close */

	if (!nextatom || GetAtomSQDist(nbatom, nextatom) < 1.0E-12 ||
		GetAtomAngle(theatom, nbatom, nextatom) > 175.0) {
		/* create a dummy point for torsion */
		p[0] = theatom->x-1.0;
		p[1] = theatom->y-1.0;
		p[2] = theatom->z-1.0;
	} else {
		p[0] = nextatom->x;
		p[1] = nextatom->y;
		p[2] = nextatom->z;
	}

	/* first hydrogen */
	InternalToCartesian(&hydrogen->x, &theatom->x, (double)length,
		&nbatom->x, (double)angle, p, (double)0.0);

	/* second hydrogen */
	if (hydrogen->next) {
		InternalToCartesian(&hydrogen->next->x, &theatom->x, (double)length,
		&nbatom->x, (double)angle, &hydrogen->x, (double)180.0);
	}

	return 1;
}

/*****************************************************
 ****  Routines for Adding Hydrogens to Carbon    ****
 *****************************************************/

/*
 * carbon --- carbon to which hydrogen is added 
 * nbatom --- neighbor atom to carbon 
 * length --- C-H bond length
 */
static AtomPtr	ca1sp (AtomPtr carbon, AtomPtr nbatom, float length)
{
	AtomPtr	hydrogen;	/* hydrogen to be added */
	double	dx, dy, dz, t;

	if (!carbon || !nbatom) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	dx = carbon->x - nbatom->x;
	dy = carbon->y - nbatom->y;
	dz = carbon->z - nbatom->z;
	t  = length / sqrt((double)(dx*dx + dy*dy + dz*dz));

	hydrogen->x = carbon->x + dx * t;
	hydrogen->y = carbon->y + dy * t;
	hydrogen->z = carbon->z + dz * t;

	return hydrogen;
}

/*
 * carbon --- carbond atom to which H is added
 * nbatom1, nbatom2 --- neighbor atoms 
 * length --- C-H bond length
 */
static AtomPtr	ca1sp2 (AtomPtr carbon, AtomPtr nbatom1, AtomPtr nbatom2, float length)
{
	AtomPtr	hydrogen;	/* hydrogen to be added */
	double	x[2], y[2], z[2];
	double	o[5], p[5], q[5];
	double	a, b, c, t;

	if (!carbon || !nbatom1 || !nbatom2) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	x[0] = carbon->x;  y[0] = carbon->y;  z[0] = carbon->z;
	o[0] = nbatom1->x; p[0] = nbatom1->y; q[0] = nbatom1->z;
	o[1] = nbatom2->x; p[1] = nbatom2->y; q[1] = nbatom2->z;

	a = o[0] - x[0];
	b = p[0] - y[0];
	c = q[0] - z[0];
	
	t = 1.0/sqrt((double)(a*a + b*b + c*c));
	
	o[2] = x[0] + a*t;
	p[2] = y[0] + b*t;
	q[2] = z[0] + c*t;
	
	a = o[1] - x[0];
	b = p[1] - y[0];
	c = q[1] - z[0];
	
	t = 1.0/sqrt((double)(a*a + b*b + c*c));
	
	o[3] = x[0] + a*t;
	p[3] = y[0] + b*t;
	q[3] = z[0] + c*t;

	x[1] = 0.5*(o[2] + o[3]);
	y[1] = 0.5*(p[2] + p[3]);
	z[1] = 0.5*(q[2] + q[3]);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = length/sqrt((double)(a*a + b*b + c*c));

	hydrogen->x = x[0] + a*t;
	hydrogen->y = y[0] + b*t;
	hydrogen->z = z[0] + c*t;

	return hydrogen;
}

static AtomPtr	ca1sp3 (AtomPtr carbon, AtomPtr nbatom1, AtomPtr nbatom2, AtomPtr nbatom3, float length)
{
	AtomPtr	hydrogen;	/* hydrogen to be added */
	double	x[2], y[2], z[2];
	double	a, b, c, t;

	if (!carbon || !nbatom1 || !nbatom2 || !nbatom3) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	x[0] = carbon->x; y[0] = carbon->y; z[0] = carbon->z;

	x[1] = (1.0/3.0)*(nbatom1->x + nbatom2->x + nbatom3->x);
	y[1] = (1.0/3.0)*(nbatom1->y + nbatom2->y + nbatom3->y);
	z[1] = (1.0/3.0)*(nbatom1->z + nbatom2->z + nbatom3->z);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = length/sqrt((double)(a*a + b*b + c*c));

	hydrogen->x = x[0] + a*t;
	hydrogen->y = y[0] + b*t;
	hydrogen->z = z[0] + c*t;

	return hydrogen;
}

/* add two hydrogens */
static AtomPtr	ca2sp2 (AtomPtr carbon, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hydrogen;	/* hydrogen to be added */

	if (!carbon || !nbatom) return NULL;

	/* create two hydrogens */
	hydrogen = NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;

	newatom_from_one_neighbor(carbon, nbatom, nextatom, length, 120.0, hydrogen);

	return hydrogen;
}

static AtomPtr	ca2sp3 (AtomPtr carbon, AtomPtr nbatom1, AtomPtr nbatom2, float length)
{
	AtomPtr	atom, hydrogen;	/* hydrogen to be added */

	float vec1[3], vec2[3], vec3[3];
	float x[5], y[5], z[5];
	float o[4], p[4], q[4];
	float a, b, c, t;

	if (!carbon || !nbatom1 || !nbatom2) return NULL;

	/* create two atoms */
	hydrogen = NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;

	MakeAtomVectors(carbon, nbatom1, nbatom2, vec1, vec2);
	fvcross(vec1, vec2, vec3);

	x[0] = carbon->x;  y[0] = carbon->y;  z[0] = carbon->z;
	o[0] = nbatom1->x; p[0] = nbatom1->y; q[0] = nbatom1->z;
	o[1] = nbatom2->x; p[1] = nbatom2->y; q[1] = nbatom2->z;
	
	a = o[0] - x[0];
	b = p[0] - y[0];
	c = q[0] - z[0];
	
	t = 1.0/sqrt((double)(a*a + b*b + c*c));
	
	o[2] = x[0] + a*t;
	p[2] = y[0] + b*t;
	q[2] = z[0] + c*t;
	
	a = o[1] - x[0];
	b = p[1] - y[0];
	c = q[1] - z[0];
	
	t = 1.0/sqrt((double)(a*a + b*b + c*c));
	
	o[3] = x[0] + a*t;
	p[3] = y[0] + b*t;
	q[3] = z[0] + c*t;

	x[1] = 0.5*(o[2] + o[3]);
	y[1] = 0.5*(p[2] + p[3]);
	z[1] = 0.5*(q[2] + q[3]);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = sin_35*length/sqrt((double)(a*a + b*b + c*c));

	x[2] = x[0] + a*t;
	y[2] = y[0] + b*t;
	z[2] = z[0] + c*t;

	t = cos_35*length/VLENGTH(vec3);

	x[3] = x[2] + vec3[0]*t;
	y[3] = y[2] + vec3[1]*t;
	z[3] = z[2] + vec3[2]*t;

	x[4] = x[2] - vec3[0]*t;
	y[4] = y[2] - vec3[1]*t;
	z[4] = z[2] - vec3[2]*t;

	atom = hydrogen;
	atom->x = x[3];
	atom->y = y[3];
	atom->z = z[3];

	atom = hydrogen->next;
	atom->x = x[4];
	atom->y = y[4];
	atom->z = z[4];

	return hydrogen;
}

static AtomPtr	ca3sp3 (AtomPtr carbon, AtomPtr nbatom, AtomPtr nextatom, float length)
{
 	AtomPtr	hydrogen, hatomlist;

	if (!carbon || !nbatom) return NULL;
	if (!(hatomlist = NewAtom())) return NULL;

	newatom_from_one_neighbor(carbon, nbatom, nextatom, length, 109.0, hatomlist);

	hydrogen = ca2sp3(carbon, nbatom, hatomlist, length);
	if (hydrogen) EnterAtom(hydrogen, &hatomlist);
	return hatomlist;
}

/*****************************************************
 ****  Routines for Adding Hydrogens to Nitrogen  ****
 *****************************************************/

static AtomPtr	a2ion (AtomPtr nitrogen, AtomPtr nbatom1, AtomPtr nbatom2, float length)
{
	AtomPtr	hydrogen;
	float	vec1[3], vec2[3], vec3[3];
	float	x[5], y[5], z[5];
	float	a, b, c, t;

	if (!nitrogen || !nbatom1 || !nbatom2) return NULL;
	hydrogen = NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;

	MakeAtomVectors(nitrogen, nbatom1, nbatom2, vec1, vec2);
	fvcross(vec1, vec2, vec3);

	x[0] = nitrogen->x;
	y[0] = nitrogen->y;
	z[0] = nitrogen->z;

	x[1] = 0.5*(nbatom1->x + nbatom2->x);
	y[1] = 0.5*(nbatom1->y + nbatom2->y);
	z[1] = 0.5*(nbatom1->z + nbatom2->z);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = sin_35*length/sqrt((double)(a*a + b*b + c*c));

	x[2] = x[0] + a*t;
	y[2] = y[0] + b*t;
	z[2] = z[0] + c*t;

	t = cos_35*length/VLENGTH(vec3);

	x[3] = x[2] + vec3[0]*t;
	y[3] = y[2] + vec3[1]*t;
	z[3] = z[2] + vec3[2]*t;

	x[4] = x[2] - vec3[0]*t;
	y[4] = y[2] - vec3[1]*t;
	z[4] = z[2] - vec3[2]*t;

	hydrogen->x = x[3];
	hydrogen->y = y[3];
	hydrogen->z = z[3];

	hydrogen->next->x = x[4];
	hydrogen->next->y = y[4];
	hydrogen->next->z = z[4];

	return hydrogen;
}

static AtomPtr	a3ion (AtomPtr nitrogen, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hydrogen, hatomlist;

	if (!nitrogen || !nbatom || !nextatom) return NULL;
	if (!(hatomlist = NewAtom())) return NULL;

	newatom_from_one_neighbor(nitrogen, nbatom, nextatom, length, 109.0, hatomlist);
	hydrogen = a2ion(nitrogen, nbatom, hatomlist, length);
	if (hydrogen) EnterAtom(hydrogen, &hatomlist);

	return hatomlist;
}

static AtomPtr	na1sp2 (AtomPtr nitrogen, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hydrogen;
	AtomPtr	hydrogen1, hydrogen2;

	if (!nitrogen || !nbatom || !nextatom) return NULL;

	/* create two hydrogens. pass them to newatom_from_one_neighbor.
	discard the first hydrogen and return the second one. */
	hydrogen = NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;

	newatom_from_one_neighbor(nitrogen, nbatom, nextatom, length, 120.0, hydrogen);

	/* free the first hydrogen and return the second one */
	hydrogen1 = hydrogen;
	hydrogen2 = hydrogen->next;
	hydrogen1->next = NULL;
	hydrogen2->prev = NULL;
	FreeAtom(&hydrogen1);

	return hydrogen2;
}

static AtomPtr	na1sp3 (AtomPtr nitrogen, AtomPtr nbatom1, AtomPtr nbatom2, float length)
{
	AtomPtr	hydrogen;
	float	vec1[3], vec2[3], vec3[3];
	float	x[6], y[6], z[6];
	float	a, b, c, t;
	float	len1, len2;

	if (!nitrogen || !nbatom1 || !nbatom2) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	MakeAtomVectors(nitrogen, nbatom1, nbatom2, vec1, vec2);
	fvcross(vec1, vec2, vec3);

	x[0] = nitrogen->x;
	y[0] = nitrogen->y;
	z[0] = nitrogen->z;

	len1 = 1.0/VLENGTH(vec1);
	x[1] = x[0] + vec1[0]*len1;
	y[1] = y[0] + vec1[1]*len1;
	z[1] = z[0] + vec1[2]*len1;

	len2 = 1.0/VLENGTH(vec2);
	x[2] = x[0] + vec2[0]*len2;
	y[2] = y[0] + vec2[1]*len2;
	z[2] = z[0] + vec2[2]*len2;

	x[3] = 0.5*(x[1] + x[2]);
	y[3] = 0.5*(y[1] + y[2]);
	z[3] = 0.5*(z[1] + z[2]);

	a = x[0] - x[3];
	b = y[0] - y[3];
	c = z[0] - z[3];

	t = sin_36*length/sqrt((double)(a*a + b*b + c*c));

	x[4] = x[0] + a*t;
	y[4] = y[0] + b*t;
	z[4] = z[0] + c*t;

	t = cos_36*length/VLENGTH(vec3);

	x[5] = x[4] + vec3[0]*t;
	y[5] = y[4] + vec3[1]*t;
	z[5] = z[4] + vec3[2]*t;

	hydrogen->x = x[5];
	hydrogen->y = y[5];
	hydrogen->z = z[5];

	return hydrogen;
}

static AtomPtr	na2sp3 (AtomPtr nitrogen, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	ATOM	tmpatom;
	AtomPtr	hydrogen;

	if (!nitrogen || !nbatom || !nextatom) return NULL;
	hydrogen = NULL;

	BZERO(&tmpatom, sizeof(tmpatom));
	newatom_from_one_neighbor(nitrogen, nbatom, nextatom, length, 109.0, &tmpatom);

	hydrogen = a2ion(nitrogen, nbatom, &tmpatom, length);

	return hydrogen;
}

static AtomPtr	a1special (AtomPtr nitrogen, AtomPtr nbatom1, AtomPtr nbatom2, float length)
{
	AtomPtr	hydrogen;
	float	x[2], y[2], z[2];
	float	a, b, c, t;

	if (!nitrogen || !nbatom1 || !nbatom2) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	x[0] = nitrogen->x;
	y[0] = nitrogen->y;
	z[0] = nitrogen->z;

	x[1] = 0.5*(nbatom1->x + nbatom2->x);
	y[1] = 0.5*(nbatom1->y + nbatom2->y);
	z[1] = 0.5*(nbatom1->z + nbatom2->z);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = length/sqrt((double)(a*a + b*b + c*c));

	hydrogen->x = x[0] + a*t;
	hydrogen->y = y[0] + b*t;
	hydrogen->z = z[0] + c*t;

	return hydrogen;
}

static AtomPtr	a2amide (AtomPtr nitrogen, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hydrogen;

	if (!nitrogen || !nbatom || !nextatom) return NULL;
	hydrogen = NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;
	if (!EnterNewAtom(&hydrogen)) return NULL;

	newatom_from_one_neighbor(nitrogen, nbatom, nextatom, length, 120.0, hydrogen);

	return hydrogen;
}

static AtomPtr	a1ion (AtomPtr nitrogen, AtomPtr nbatom1, AtomPtr nbatom2, AtomPtr nbatom3, float length)
{
	AtomPtr	hydrogen;
	float	x[2], y[2], z[2];	
	float	a, b, c, t;

	if (!nitrogen || !nbatom1 || !nbatom2 || !nbatom3) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	x[0] = nitrogen->x;
	y[0] = nitrogen->y;
	z[0] = nitrogen->z;

	x[1] = (1.0/3.0)*(nbatom1->x + nbatom2->x + nbatom3->x);
	y[1] = (1.0/3.0)*(nbatom1->y + nbatom2->y + nbatom3->y);
	z[1] = (1.0/3.0)*(nbatom1->z + nbatom2->z + nbatom3->z);

	a = x[0] - x[1];
	b = y[0] - y[1];
	c = z[0] - z[1];

	t = length/sqrt((double)(a*a + b*b + c*c));

	hydrogen->x = x[0] + a*t;
	hydrogen->y = y[0] + b*t;
	hydrogen->z = z[0] + c*t;

	return hydrogen;
}

/*****************************************************
 ****  Routines for Adding Hydrogens to Oxygen    ****
 *****************************************************/

/*
 * length --- O-H bond length 
 * angle  --- angle in degree
 */
static AtomPtr	oa1 (AtomPtr oxygen, AtomPtr nbatom, AtomPtr nextatom, float length, float angle)
{
	AtomPtr	hydrogen;

	if (!oxygen || !nbatom) return NULL;
	if (!(hydrogen = NewAtom())) return NULL;

	newatom_from_one_neighbor(oxygen, nbatom, nextatom, length, angle, hydrogen);

	return hydrogen;
}

/*********************************************************************
 ****  Executive Routines for Adding Hydrogens to C, N, O, and S  ****
 *********************************************************************/

/* ca(m)sp(n)
 * m --- number of hydrogens to be added
 * n --- atom type of carbon (sp, sp2, sp3)
 */

/* add hydrogens to carbon.
 * nextatom = the first encountered neighbor of neighbor.
 * length   = carbon-hydrogen length. returns list of hydrogen atoms.
 */
static AtomPtr	add_h2c (AtomPtr carbon, AtomPtr nbatom1, AtomPtr nbatom2, AtomPtr nbatom3, AtomPtr nextatom, float length)
{
	AtomPtr	hatomlist = NULL;
	double	angle, bondlength;

	if (!carbon) return NULL;
	hatomlist = NULL;

	switch (carbon->nneighbors) {
	case 3:
		angle = GetAtomAngle(nbatom1, carbon, nbatom2)
			+ GetAtomAngle(nbatom1, carbon, nbatom3)
			+ GetAtomAngle(nbatom3, carbon, nbatom2);

		if (angle < MAX_SUM_ANG) hatomlist = ca1sp3(carbon, nbatom1, nbatom2, nbatom3, length);
		break;

	case 2:
		angle = GetAtomAngle(nbatom1, carbon, nbatom2);

		if (angle < MAX_C1_ANG) hatomlist = ca2sp3(carbon, nbatom1, nbatom2, length);
		else if (angle < MAX_C2_ANG) hatomlist = ca1sp2(carbon, nbatom1, nbatom2, length);
		break;

	case 1:
		bondlength = GetAtomSQDist(carbon, nbatom1);	/* bondlength^2 */

		if (bondlength > MAX_CC2_LEN) hatomlist = ca3sp3(carbon, nbatom1, nextatom, length);
		else if (bondlength > MAX_CC3_LEN) hatomlist = ca2sp2(carbon, nbatom1, nextatom, length);
		else hatomlist = ca1sp(carbon, nbatom1, length);
		break;

	default:
		break;
	}

	return hatomlist;
}

static AtomPtr	add_h2n (AtomPtr nitrogen, AtomPtr nbatom1, AtomPtr nbatom2, AtomPtr nbatom3, AtomPtr nextatom, float length)
{
	AtomPtr	hatomlist;
	float	angle, bondlength;

	if (!nitrogen) return NULL;
	hatomlist = NULL;

	if (nitrogen->charge > 0.6) {
		switch (nitrogen->nneighbors) {
		case 1:
			hatomlist = a3ion(nitrogen, nbatom1, nextatom, length);
			break;
		case 2:
			hatomlist = a2ion(nitrogen, nbatom1, nbatom2, length);
			break;
		case 3:
			hatomlist = a1ion(nitrogen, nbatom1, nbatom2, nbatom3, length);
			break;
		}
		return hatomlist;
	}

	switch (nitrogen->nneighbors) {
	case 2:
		angle = GetAtomAngle(nbatom1, nitrogen, nbatom2);

		if (angle < MIN_AMIDE_ANG) hatomlist = na1sp3(nitrogen, nbatom1, nbatom2, length);
		else hatomlist = a1special(nitrogen, nbatom1, nbatom2, length);
		break;

	case 1:
		bondlength = GetAtomDistance(nitrogen, nbatom1);
		bondlength *= bondlength;

		if (bondlength > MAX_AMIDE_LEN_SQ) hatomlist = na2sp3(nitrogen, nbatom1, nextatom, length);
		else if (bondlength > MAX_CN2_LEN_SQ) hatomlist = a2amide(nitrogen, nbatom1, nextatom, length);
		else if (bondlength > MAX_CN3_LEN_SQ) hatomlist = na1sp2(nitrogen, nbatom1, nextatom, length);
		break;
	}

	return hatomlist;
}

static AtomPtr	add_h2o (AtomPtr oxygen, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hatomlist;
	float	dist; 

	if (!oxygen || oxygen->nneighbors != 1 || !nbatom) return NULL;
	hatomlist = NULL;

	dist = GetAtomDistance(oxygen, nbatom);
	dist *= dist;

	if (dist > MAX_CO2_LEN_SQ && oxygen->nneighbors == 1) hatomlist = oa1(oxygen, nbatom, nextatom, length, 105.0);
	return hatomlist;
}

static AtomPtr	add_h2s (AtomPtr sulfur, AtomPtr nbatom, AtomPtr nextatom, float length)
{
	AtomPtr	hatomlist;
	float	dist; 

	if (!sulfur || sulfur->nneighbors != 1 || !nbatom) return NULL;
	hatomlist = NULL;
	dist = GetAtomDistance(sulfur, nbatom);
	dist *= dist;

	if (dist > MAX_CS2_LEN_SQ && sulfur->nneighbors == 1) hatomlist = oa1(sulfur, nbatom, nextatom, length, 97.0);
	return hatomlist;
}

/*
 * Find an atom which can be used as a reference for torsion angle
 * to a newly created atom.
 *
 * refatom                   nbatom3
 *        \                 /
 *         nbatom1---theatom
 *                          \
 *                           nbatom2
 */

static int	find_torsion_reference_atom (AtomPtr theatom, AtomPtr *nbatom1_ret, AtomPtr *nbatom2_ret,
	AtomPtr *nbatom3_ret, AtomPtr *refatom_ret)
{
	AtomPtr	nbatom, refatom=NULL;
	int	i, j, an;

	if (!theatom || !nbatom1_ret || !nbatom2_ret || !nbatom3_ret || !refatom_ret) return 0;
	*nbatom1_ret = theatom->nbatom[0];
	*nbatom2_ret = theatom->nbatom[1];
	*nbatom3_ret = theatom->nbatom[2];
	*refatom_ret = NULL;

	an = theatom->an;

	if (an == 8) {
		for(i=0;i<theatom->nneighbors;i++) {
			nbatom = theatom->nbatom[i];
			for(j=0;j<nbatom->nneighbors;j++) {
				if (!(refatom = nbatom->nbatom[j]) || refatom == theatom) continue;
				if ((refatom->an == 8 || refatom->an == 16) &&
				    (GetAtomAngle(theatom, nbatom, refatom) < 175.0)) {
					*nbatom1_ret = nbatom;
					*refatom_ret = refatom;
					if (i == 1) {
						*nbatom2_ret = theatom->nbatom[0];
					} else if (i >= 2) {
						*nbatom2_ret = theatom->nbatom[0];
						*nbatom3_ret = theatom->nbatom[1];
					}
					return 1;
				}
			}
		}
	}

	for(i=0;i<theatom->nneighbors;i++) {
		nbatom = theatom->nbatom[i];
		for(j=0;j<nbatom->nneighbors;j++) {
			if (!(refatom = nbatom->nbatom[j]) || refatom == theatom) continue;
			if (GetAtomAngle(theatom, nbatom, refatom) < 175.0) {
				*nbatom1_ret = nbatom;
				*refatom_ret = refatom;
				if (i == 1) {
					*nbatom2_ret = theatom->nbatom[0];
				} else if (i >= 2) {
					*nbatom2_ret = theatom->nbatom[0];
					*nbatom3_ret = theatom->nbatom[1];
				}
				return 1;
			}
		}
	}

	return 0;
}

/*************************************************************
 ******** public functions ***********************************
 *************************************************************/

/* add hydrogens to all atoms in molecule. */
int	AddHydrogenToMol (MolPtr mol)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a, h, hydrogen, hatomlist;
	BondPtr	bond, hbondlist;
	int	serno;

	if (!mol || !mol->chain || !mol->chain->residue) return 0;
	hbondlist = NULL;
	ForEachChainRes(mol->chain, c, r) {
		hatomlist = NULL;
		ForEachAtom(r->atom, a) {
			if (!(hydrogen = AddHydrogenToAtom(a))) continue;
			/* create bonds between atom and newly created hydrogens */
			ForEachAtom(hydrogen, h) {
				h->residue = r;
				if (!(bond = EnterNewBond(&hbondlist))) continue;
				bond->atom1 = a;
				bond->atom2 = h;
				bond->type = B_SINGLE;
			}
			EnterAtom(hydrogen, &hatomlist);
		}
		EnterAtom(hatomlist, &r->atom);
	}
	EnterBond(hbondlist, &mol->bond);
	GetNeighbor(mol);

	/* update atom serial number */
	serno = 0;
	ForEachChainResAtom(mol->chain, c, r, a) {
		serno++;
		a->serno = serno;
	}
	return 1;
}

/*
 * add hydrogens to theatom and return the resulting atoms in hatomlist.
 * theatom = atom to which hydrogens are to be added
 */
AtomPtr	AddHydrogenToAtom (AtomPtr theatom)
{
	AtomPtr	hatomlist;	/* newly added hydrogens */
	AtomPtr	nextatom;	/* neighbor of neighbor */
	AtomPtr	atom, nbatom1, nbatom2, nbatom3;
	float	length, hrad, orad;

	if (!theatom) return NULL;
	hatomlist = NULL;

	/* find next to neighbor theatom (theatom--nbatom1--nextatom) */
	find_torsion_reference_atom(theatom, &nbatom1, &nbatom2, &nbatom3, &nextatom);

	hrad = 0.35;
	switch(theatom->an) {
		case 6:  orad = 0.70; break;
		case 7:  orad = 0.65; break;
		case 8:  orad = 0.60; break;
		case 16: orad = 1.00; break;
		default: orad = 0.70; break;
	}

	length = (hrad + orad);

	switch (theatom->an) {
		case 6:		
			hatomlist = add_h2c(theatom, nbatom1, nbatom2, nbatom3, nextatom, length);
			break;
		case 7:
			hatomlist = add_h2n(theatom, nbatom1, nbatom2, nbatom3, nextatom, length);
			break;
		case 8:
			hatomlist = add_h2o(theatom, nbatom1, nextatom, length);
			break;
		case 16:
			hatomlist = add_h2s(theatom, nbatom1, nextatom, length);
			break;
		default:
			break;
	}

	for(atom=hatomlist;atom;atom=atom->next) {
		atom->an = 1;
		atom->refno = AN2ATOMREFNO(atom->an);
		atom->col = GetElemCPKColor((int)atom->an);
	}

	return hatomlist;
}

/*
 * add a specified number of hydrogens to theatom and
 * return the resulting atoms in hatomlist.
 * theatom = atom to which hydrogens are to be added
 */
AtomPtr	AddHydrogenToAtomByNumber (AtomPtr theatom, int nhydrogens)
{
	AtomPtr	hatomlist;	/* newly added hydrogens */
	AtomPtr	nextatom;	/* neighbor of neighbor */
	AtomPtr	atom, nbatom1, nbatom2, nbatom3;
	int	nneighbors;
	float	length, hrad, orad;

	if (!theatom || nhydrogens <= 0) return NULL;
	hatomlist = NULL;
	nneighbors = theatom->nneighbors;

	/* find next to neighbor theatom (theatom--nbatom1--nextatom) */
	find_torsion_reference_atom(theatom, &nbatom1, &nbatom2, &nbatom3, &nextatom);

	hrad = 0.35;
	switch(theatom->an) {
		case 6:  orad = 0.70; break;
		case 7:  orad = 0.65; break;
		case 8:  orad = 0.60; break;
		case 16: orad = 1.00; break;
		default: orad = 0.70; break;
	}

	length = (hrad + orad);

	switch (theatom->an) {
	case 6:		
		switch (nhydrogens) {
		case 1:
			switch (nneighbors) {
			case 1:
				hatomlist = ca1sp(theatom, nbatom1, length);
				break;
			case 2:
				hatomlist = ca1sp2(theatom, nbatom1, nbatom2, length);
				break;
			case 3:
				hatomlist = ca1sp3(theatom, nbatom1, nbatom2, nbatom3, length);
				break;
			}
			break;

		case 2:
			switch (nneighbors) {
			case 1:
				hatomlist = ca2sp2(theatom, nbatom1, nextatom, length);
				break;
			case 2:
				hatomlist = ca2sp3(theatom, nbatom1, nbatom2, length);
				break;
			}
			break;

		case 3:
			switch (nneighbors) {
			case 1:
				hatomlist = ca3sp3(theatom, nbatom1, nextatom, length);
				break;
			}
			break;
		}
		break;

	case 7:
		switch (nhydrogens) {
		case 1:
			switch (nneighbors) {
			case 1:
				hatomlist = na1sp2(theatom, nbatom1, nextatom, length);
				break;
			case 2:
				hatomlist = na1sp3(theatom, nbatom1, nbatom2, length);
				break;
			case 3:
				hatomlist = a1ion(theatom, nbatom1, nbatom2, nbatom3, length);
				break;
			}
			break;

		case 2:
			switch (nneighbors) {
			case 1:
				hatomlist = na2sp3(theatom, nbatom1, nextatom, length);
				break;
			case 2:
				hatomlist = a2ion(theatom, nbatom1, nbatom2, length);
				break;
			}
			break;

		case 3:
			switch (nneighbors) {
			case 1:
				hatomlist = a3ion(theatom, nbatom1, nextatom, length);
				break;
			}
			break;
		}
		break;

	case 8:
		switch (nhydrogens) {
		case 1:
			switch (nneighbors) {
			case 1:
				hatomlist = oa1(theatom, nbatom1, nextatom, length, 105.0);
				break;
			}
			break;
		}
		break;

	case 16:
		switch (nhydrogens) {
		case 1:
			switch (nneighbors) {
			case 1:
				hatomlist = oa1(theatom, nbatom1, nextatom, length, 97.0);
				break;
			}
			break;
		}
		break;

	default:
		break;
	}

	for(atom=hatomlist;atom;atom=atom->next) {
		atom->an = 1;
		atom->refno = AN2ATOMREFNO(atom->an);
		atom->col = GetElemCPKColor((int)atom->an);
	}

	return hatomlist;
}

/* add user-specified number of hydrogens to all atoms in molecule.
   If a number is negative, add hydrogens as usual, i.e., use
   AddHydrogenToAtom().
*/
int	AddHydrogenToMolByNumber (MolPtr mol, int nC, int nN, int nO, int nS)
{
	ChainPtr	c;
	ResiduePtr	r;
	AtomPtr	a, h, hydrogen, hatomlist;
	BondPtr	bond, hbondlist;
	int	serno, n;

	if (!mol || !mol->chain || !mol->chain->residue) return 0;
	hbondlist = NULL;
	ForEachChainRes(mol->chain, c, r) {
		hatomlist = NULL;
		ForEachAtom(r->atom, a) {
			n = -1;
			switch (a->an) {
			case  6: n = nC; break;
			case  7: n = nN; break;
			case  8: n = nO; break;
			case 16: n = nS; break;
			}
			if (n < 0) {
				if (!(hydrogen = AddHydrogenToAtom(a))) continue;
			} else {
				if (!(hydrogen = AddHydrogenToAtomByNumber(a, n))) continue;
			}

			/* create bonds between atom and newly created hydrogens */
			ForEachAtom(hydrogen, h) {
				h->residue = r;
				if (!(bond = EnterNewBond(&hbondlist))) continue;
				bond->atom1 = a;
				bond->atom2 = h;
				bond->type = B_SINGLE;
			}
			EnterAtom(hydrogen, &hatomlist);
		}
		EnterAtom(hatomlist, &r->atom);
	}
	EnterBond(hbondlist, &mol->bond);
	GetNeighbor(mol);

	/* update atom serial number */
	serno = 0;
	ForEachChainResAtom(mol->chain, c, r, a) {
		serno++;
		a->serno = serno;
	}
	return 1;
}

void	SetAddHydrogenPref (int id, float val)
{
	if (id < 0 || id > H_PREF_MAX || val <= 0) return;
	switch (id) {
	case H_MAX_SUM_ANG:
		MAX_SUM_ANG = val;
		break;
	case H_MAX_C1_ANG:
		MAX_C1_ANG = val;
		break;
	case H_MAX_C2_ANG:
		MAX_C2_ANG = val;
		break;
	case H_MAX_CC2_LEN:
		MAX_CC2_LEN = val;
		break;
	case H_MAX_CC3_LEN:
		MAX_CC3_LEN = val;
		break;
	case H_MIN_AMIDE_ANG:
		MIN_AMIDE_ANG = val;
		break;
	case H_MAX_CN2_LEN:
		MAX_CN2_LEN = val;
		MAX_CN2_LEN_SQ = val*val;
		break;
	case H_MAX_AMIDE_LEN:
		MAX_AMIDE_LEN = val;
		MAX_AMIDE_LEN_SQ = val*val;
		break;
	case H_MAX_CN3_LEN:
		MAX_CN3_LEN = val;
		MAX_CN3_LEN_SQ = val*val;
		break;
	case H_MAX_CO2_LEN:
		MAX_CO2_LEN = val;
		MAX_CO2_LEN_SQ = val*val;
		break;
	case H_MAX_CS2_LEN:
		MAX_CS2_LEN = val;
		MAX_CS2_LEN_SQ = val*val;
		break;
	}
}

float	GetAddHydrogenPref (int id, float *val_ret)
{
	float	val;

	if (id < 0 || id > H_PREF_MAX) return 0.0;
	switch (id) {
	case H_MAX_SUM_ANG:
		val = MAX_SUM_ANG;
		break;
	case H_MAX_C1_ANG:
		val = MAX_C1_ANG;
		break;
	case H_MAX_C2_ANG:
		val = MAX_C2_ANG;
		break;
	case H_MAX_CC2_LEN:
		val = MAX_CC2_LEN;
		break;
	case H_MAX_CC3_LEN:
		val = MAX_CC3_LEN;
		break;
	case H_MIN_AMIDE_ANG:
		val = MIN_AMIDE_ANG;
		break;
	case H_MAX_CN2_LEN:
		val = MAX_CN2_LEN;
		break;
	case H_MAX_AMIDE_LEN:
		val = MAX_AMIDE_LEN;
		break;
	case H_MAX_CN3_LEN:
		val = MAX_CN3_LEN;
		break;
	case H_MAX_CO2_LEN:
		val = MAX_CO2_LEN;
		break;
	case H_MAX_CS2_LEN:
		val = MAX_CS2_LEN;
		break;
	}

	if (val_ret) *val_ret = val;
	return val;
}

int	GetHybridization (AtomPtr atom)
{
	int	hybrid = HYBRID_UNK;
	AtomPtr	nbatom1, nbatom2, nbatom3, nextatom;
	float	len, angle;

	if (!atom) return hybrid;
	find_torsion_reference_atom(atom, &nbatom1, &nbatom2, &nbatom3, &nextatom);

	switch (atom->an) {
	case 6:
		switch (atom->nneighbors) {
		case 3:
			angle = GetAtomAngle(nbatom1, atom, nbatom2)
				+ GetAtomAngle(nbatom1, atom, nbatom3)
				+ GetAtomAngle(nbatom3, atom, nbatom2);
			if (angle < MAX_SUM_ANG) hybrid = HYBRID_SP3;
			else hybrid = HYBRID_SP2;
			break;
		case 2:
			angle = GetAtomAngle(nbatom1, atom, nbatom2);
			if (angle < MAX_C1_ANG) hybrid = HYBRID_SP3;
			else if (angle < MAX_C2_ANG) hybrid = HYBRID_SP2;
			else hybrid = HYBRID_SP;
			break;
		case 1:
			len = GetAtomSQDist(atom, nbatom1);
			if (len > MAX_CC2_LEN) hybrid = HYBRID_SP3;
			else if (len > MAX_CC3_LEN) hybrid = HYBRID_SP2;
			else hybrid = HYBRID_SP;
			break;
		default:
			hybrid = HYBRID_SP3;
			break;
		}
		break;
	case 7:
		switch (atom->nneighbors) {
		case 4:
			hybrid = HYBRID_SP3_PLUS;
			break;
		case 3:
			hybrid = HYBRID_SP3;
			break;
		case 2:
			angle = GetAtomAngle(nbatom1, atom, nbatom2);
			if (angle < MIN_AMIDE_ANG) hybrid = HYBRID_SP3;
			else hybrid = HYBRID_SP2;
		case 1:
			len = GetAtomSQDist(atom, nbatom1);
			if (len > MAX_AMIDE_LEN_SQ) hybrid = HYBRID_SP3;
			else if (len > MAX_CN2_LEN_SQ) hybrid = HYBRID_SP2;	/* amide */
			else if (len > MAX_CN3_LEN_SQ) hybrid = HYBRID_SP2;	/* imine */
			else hybrid = HYBRID_SP;	/* nitrile or N2 */
			break;
		default:
			hybrid = HYBRID_SP3;
			break;
		}
		break;
	case 8:
		len = GetAtomSQDist(atom, nbatom1);
		if (len > MAX_CO2_LEN_SQ) hybrid = HYBRID_SP3;
		else hybrid = HYBRID_SP2;
		break;
	case 16:
		len = GetAtomSQDist(atom, nbatom1);
		if (len > MAX_CS2_LEN_SQ) hybrid = HYBRID_SP3;
		else hybrid = HYBRID_SP2;
		break;
	default:
		hybrid = HYBRID_UNK;
		break;
	}

	return hybrid;
}

