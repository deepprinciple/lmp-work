#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "macros.h"
#include "list.h"

/*********************************************
* Functions for Handling Singly Linked Lists *
**********************************************/

ListPtr	NewList ()
{
	ListPtr	newptr;
	if (!(newptr = (ListPtr) calloc(1, sizeof(List)))) return NULL;
	return newptr;
}

ListPtr	EnterNewList (ListPtr *top)
{
	ListPtr	newptr;

	if (!(newptr = NewList())) return NULL;
	EnterList(newptr, top);
	return newptr;
}

void	EnterList (ListPtr newptr, ListPtr *top)
{
	ListPtr	slist;

	if (!newptr || !top) return;
	if (!*top) {	/* first item */
		*top = newptr;
	} else {
		/* go to the end */
		for(slist = *top; slist->next; slist = slist->next) ;

		/* append the new item */
		slist->next = newptr;
	}
}

/* newlist is inserted before node */
void	InsertList (ListPtr *top, ListPtr node, ListPtr newlist)
{
	ListPtr	newlistend, prev;

	if (!top || !node || !newlist) return;
	/* go to the end */
	for(newlistend=newlist;newlistend && newlistend->next;newlistend=newlistend->next) ;

	for(prev=(*top);prev;prev=prev->next) {
		if (prev->next == node) break;
	}

	if (node == (*top)) {
		*top = newlist;
		newlistend->next = node;
	} else {
		if (!prev) return;
		prev->next = newlist;
		newlistend->next = node;
	}
}

ListPtr	CopyList (ListPtr oldlist)
{
	ListPtr	oldptr, newptr, newlist;

	if (!oldlist) return NULL;

	newlist = NULL;
	for(oldptr=oldlist;oldptr;oldptr=oldptr->next) {
		if (!(newptr = NewList())) {
			if (newlist) FreeList(&newlist);
			return NULL;
		}
		BCOPY((void *)oldptr, (void *)newptr, (int)(sizeof(List)));
		newptr->next = NULL;
		EnterList(newptr, &newlist);
	}
	return newlist;
}

void	PrintList (ListPtr top)
{
	int	count = 1;

	while (top) {
		fprintf(stderr, "List %d = %d data = %d\n", count, (int)top, (int)top->data1);
		top = top->next; count++;
	}
	fprintf(stderr, "\n");
}

void	DeleteList (ListPtr del, ListPtr *top)
{
	register ListPtr	slist;

	if (!del || !top || !*top) return;
	/* new first item */
	if (del == *top) {
		*top = del->next;
	} else {
		/* find the item whose 'next' field points to 'del' */
		for(slist = *top; slist->next != del; slist = slist->next) ;
		slist->next = del->next;
	}
}

void	RemoveList (ListPtr del, ListPtr *top)
{
	if (!del || !top || !*top) return;
	DeleteList(del, top);
	FreeListEntry(&del);
}

void	FreeList (ListPtr *top)
{
	ListPtr	tmp, p;

	if (!top || !*top) return;
	p = *top;
	while (p) {tmp = p->next; free((void *)p); p = tmp;}
	*top = NULL;
}

void	FreeListEntry (ListPtr *list)
{
	if (!list || !*list) return;
	free((void *)*list);
	*list = NULL;
}

ListPtr	SearchList (ListPtr top, int i)
{
	register int	count = 1;

	if (i<=0) return 0;
	while (top) {
		if (i==count) return top;
		top = top->next; count++;
	}
	return NULL;
}
		
ListPtr	FindListData (ListPtr top, int data)
{
	while (top) {if (top->data1 == data) return top; top = top->next;}
	return NULL;
}

int	FindListDataEntry (ListPtr top, int data)
{
	register int	count=1;
	while (top) {
		if (top->data1 == data) return count;
		top = top->next; count++;
	}
	return 0;
}

/* compare fields of each list. if the number of entries is different,
or if fields are not the same for each list, return 0.
if (flags & M_DATA_1)	compare list->data1
if (flags & M_DATA_2)	compare list->data2
if (flags & M_DATA_3)	compare list->data3
etc.
*/
int	CompareList (ListPtr list1, ListPtr list2, unsigned int flags)
{
	if (!list1 || !list2) return 0;

	for(;list1 || list2;list1 = list1->next, list2 = list2->next) {
		if (!list1 || !list2) return 0;

		if        (flags & M_DATA_1) {
			if (list1->data1 != list2->data1) return 0;
		} else if (flags & M_DATA_2) {
			if (list1->data2 != list2->data2) return 0;
		} else if (flags & M_DATA_3) {
			if (list1->data3 != list2->data3) return 0;
		} else if (flags & M_DATA_4) {
			if (list1->data4 != list2->data4) return 0;
		} else if (flags & M_DATA_5) {
			if (list1->data5 != list2->data5) return 0;
		} else if (flags & M_DATA_6) {
			if (list1->data6 != list2->data6) return 0;
		} else if (flags & M_DATA_7) {
			if (list1->data7 != list2->data7) return 0;
		} else if (flags & M_DATA_8) {
			if (list1->data8 != list2->data8) return 0;
		}
	}
	return 1;
}

void	SwapList (ListPtr listi, ListPtr listj)
{
	List	list, *nexti, *nextj;

	nexti = listi->next;
	nextj = listj->next;
	BCOPY(listi, &list, sizeof(List));
	BCOPY(listj, listi, sizeof(List));
	BCOPY(&list, listj, sizeof(List));

	listi->next = nexti;
	listj->next = nextj;
}

static size_t	get_list_data (ListPtr list, int mode)
{
	size_t	data;

	if (!list) return (size_t)NULL;
	switch (mode) {
	case M_DATA_1: data = list->data1; break;
	case M_DATA_2: data = list->data2; break;
	case M_DATA_3: data = list->data3; break;
	case M_DATA_4: data = list->data4; break;
	case M_DATA_5: data = list->data5; break;
	case M_DATA_6: data = list->data6; break;
	case M_DATA_7: data = list->data7; break;
	case M_DATA_8: data = list->data8; break;
	default:
		data = list->data1;
		break;
	}
	return data;
}

/* sort list according to mode in ascending order */
/* mode = M_DATA_1, M_DATA_2, ... , M_DATA_8 */

void	QsortList (ListPtr toplist, int left, int right, int mode)
{
	register int	i, j;
	register size_t	data;
	register ListPtr	list, listi, listj;

	if (!toplist || left >= right) return;
	i = left; j = right;
	list = SearchList(toplist, (i+j)/2);
	data = get_list_data(list, mode);

	do {
		while ((listi = SearchList(toplist, i)) &&
			(get_list_data(listi, mode) < data) &&
			(i < right))
			i++;

		while ((listj = SearchList(toplist, j)) &&
			(get_list_data(listj, mode) > data) &&
			(j > left))
			j--;

		if (i<=j) {
			SwapList(listi, listj);
			i++;
			j--;
		}
	} while (i<=j);

	if (left < j)  QsortList(toplist, left,  j, mode);
	if (i < right) QsortList(toplist, i, right, mode);
}

int	CountList (ListPtr top)
{
	int	count=0;

	for(;top;top=top->next) count++;
	return count;
}

