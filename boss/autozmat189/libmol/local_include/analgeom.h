#ifndef __ANALGEOM_H_
#define __ANALGEOM_H_
extern double	GetAtomDistance (AtomPtr atom1, AtomPtr atom2);
extern double	GetAtomSQDist (AtomPtr atom1, AtomPtr atom2);
extern double	GetAngle (double *v1, double *v2, double *v3);
extern double	GetAtomAngle (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3);
extern double	GetDihedral (double *v1, double *v2, double *v3, double *v4);
extern double	GetAtomDihedral (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, AtomPtr atom4);
extern void	MakeAtomVectors (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *v1, float *v2);
extern int	GetAtomNormal (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float *norm);
extern int	IsAtomLinear (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomCis (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomTrans (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsAtomPerp (AtomPtr atom1, AtomPtr atom2, AtomPtr atom3, float epsilon);
extern int	IsLinear (double *v1, double *v2, double *v3, float epsilon);
#endif
