#ifndef __LIST_H_
#define __LIST_H

/* masks for QsortList */
#define M_DATA_1	0x1
#define M_DATA_2	0x2
#define M_DATA_3	0x4
#define M_DATA_4	0x8
#define M_DATA_5	0x10
#define M_DATA_6	0x20
#define M_DATA_7	0x40
#define M_DATA_8	0x80

/**************************************
* aliases for List structure elements *
***************************************/
/* general */
#define L_DATA_1	data1
#define L_DATA_2	data2
#define L_DATA_3	data3
#define L_DATA_4	data4
#define L_DATA_5	data5
#define L_DATA_6	data6
#define L_DATA_7	data7
#define L_DATA_8	data8
#define L_VAL_1		val1
#define L_VAL_2		val2
#define L_VAL_3		val3
#define L_VAL_4		val4
#define L_VAL_5		val5
#define L_VAL_6		val6
#define L_VAL_7		val7
#define L_VAL_8		val8
#define L_STR_1		str1
#define L_STR_2		str2
#define L_STR_3		str3
#define L_STR_4		str4
#define L_STR_5		str5
#define L_STR_6		str6
#define L_STR_7		str7
#define L_STR_8		str8

/* Z-matrix */
#define L_ZMAT		data1

/* boss atom types */
#define L_TYPE		data1
#define L_AN		data2
#define L_CHARGE	val1
#define L_SIGMA		val2
#define L_EPSILON	val3
#define L_AII		val4
#define L_CII		val5
#define L_AMBER		str1
#define L_COMMENT	str8

/* boss torsional parameters (L_TYPE = data1, L_COMMENT = str8) */
#define L_N1		data2
#define L_N2		data3
#define L_N3		data4
#define L_V0		val1
#define L_V1		val2
#define L_V2		val3
#define L_V3		val4

#define L_P1		val5
#define L_P2		val6
#define L_P3		val7
#define L_ETORSION	val8

/* in comments of torsion parameters, [A1-A2-A3-A4] represents the
corresponding atom types */
#define L_A1		data5
#define L_A2		data6
#define L_A3		data7
#define L_A4		data8

/* boss bond stretching and angle bending (for strbnd file) */
#define L_K		val1
#define L_EQUIL		val2
#define L_AMBER1	str1
#define L_AMBER2	str2
#define L_AMBER3	str3
#define L_AMBER4	str4

#define L_ESTRETCH	val8
#define L_EBEND		val8

/* 1,4 and nonboneded interactions */
#define L_Q1		val1
#define L_Q2		val2
#define L_AII1		val3
#define L_AII2		val4
#define L_CII1		val5
#define L_CII2		val6
#define L_ENB_ELEC	val7
#define L_ENB_VDW	val8

/* for calculation of energy (store pointers to atoms) */
#define L_ATOMPTR1	data5
#define L_ATOMPTR2	data6
#define L_ATOMPTR3	data7
#define L_ATOMPTR4	data8

/* geometry variations */
#define L_CENTER	data1
#define L_PARAM		data2
#define L_VAL		val1

/* variable bonds  (L_CENTER = data1) */
/* variable angles (L_CENTER = data1) */

/* additional bonds, bond angles,  */
#define L_ATOM1	data1
#define L_ATOM2	data2
#define L_ATOM3	data3
#define L_ATOM4	data4

/* harmonic constraints (L_ATOM1 = data1, L_ATOM2 = data2) */
#define L_RIJ0	val1
#define L_K0	val2
#define L_K1	val3
#define L_K2	val4

/* variable dihedral angles (L_CENTER = data1) */
#define L_ITYPE	data5
#define L_FTYPE	data6
#define L_IMPROPER	data8
#define L_RANGE	val1

/* additional dihedral angles
L_ATOM1 = data1, L_ATOM2 = data2, L_ATOM3 = data3, L_ATOM4 = data4,
L_ITYPE = data5, L_FTYPE = data6, L_IMPROPER = data8 */

/* domain definitions */
#define L_FIRST_DOM1	data1
#define L_LAST_DOM1	data2
#define L_FIRST_DOM2	data3
#define L_LAST_DOM2	data4

/* conformational search (L_CENTER = data1, L_PARAM = data2) */
#define L_RMIN	val1
#define L_RMAX	val2

/* local heating */
#define L_RESIDUE	data1

/* Gaussian Z-matrix variables list (var=val) */
#define L_VARNAME	str1
#define L_VARVAL	val1
#define L_VARFLAGS	data1

/* data structure for scratch space. */
#define SL	128	/* string length */
typedef struct _List {
	size_t	data1,    data2,    data3,    data4,    data5,    data6,    data7,    data8;
	float	val1,     val2,     val3,     val4,     val5,     val6,     val7,     val8;
	char	str1[SL], str2[SL], str3[SL], str4[SL], str5[SL], str6[SL], str7[SL], str8[SL];
	struct _List	*next;
	} List, *ListPtr;


#define ForEachList(list,l)	for(l=list;l;l=l->next)

#define NEW_LIST(type) \
	(type *) calloc(1, sizeof(type));

#define APPEND_LIST(p,tmp,top) \
	if (!*top) *top = p; \
	else { \
		for(tmp=(*top);tmp->next;tmp=tmp->next); \
		tmp->next = p; \
	}

#define ENTER_LIST	APPEND_LIST

#define ENTER_NEW_LIST(p,tmp,top,type)	\
	if ((p = NEW_LIST(type))) APPEND_LIST(p,tmp,top)

#define FREE_LIST(p,tmp,top) \
	if (!top || !*top) return; \
	p = *top; \
	while (p) {tmp = p->next; free((void *)p); p = tmp;} \
	*top = NULL;

#define PRINT_LIST(p,top) \
	PRINT("List:\n"); \
	ForEachList(top,p) PRINT(" %d\n", p);

#define FOR_EACH_LIST(top,p) \
	for(p=top;p;p=p->next)

#define DELETE_LIST(p,tmp,top) \
	if (p && tmp && top && *top) { \
		if (p == *top) *top = p->next; \
		else { \
			for(tmp=(*top);tmp->next!=p;tmp=tmp->next); \
			tmp->next = p->next; \
		} \
	}


extern ListPtr	NewList ();
extern ListPtr	EnterNewList (ListPtr *top);
extern void	EnterList (ListPtr newptr, ListPtr *top);
extern void	InsertList (ListPtr *top, ListPtr node, ListPtr newlist);
extern ListPtr	CopyList (ListPtr oldlist);
extern void	PrintList (ListPtr top);
extern void	DeleteList (ListPtr del, ListPtr *top);
extern void	RemoveList (ListPtr del, ListPtr *top);
extern void	FreeList (ListPtr *top);
extern void	FreeListEntry (ListPtr *list);
extern ListPtr	SearchList (ListPtr top, int i);
extern ListPtr	FindListData (ListPtr top, int data);
extern int	FindListDataEntry (ListPtr top, int data);

extern int	CompareList (ListPtr list1, ListPtr list2, unsigned int flags);
extern void	SwapList (ListPtr listi, ListPtr listj);
extern void	QsortList (ListPtr toplist, int left, int right, int mode);
extern int	CountList (ListPtr top);

#endif


