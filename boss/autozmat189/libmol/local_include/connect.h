#ifndef __CONNECT_H_
#define __CONNECT_H_

extern ListPtr	LookUpBondList (ListPtr bondlist, int a1, int a2);
extern ListPtr	LookUpAngleList (ListPtr anglelist, int a1, int a2, int a3);
extern ListPtr	LookUpDihedList (ListPtr dihedlist, int a1, int a2, int a3, int a4);
extern ListPtr	LookUpImproperList (ListPtr improperlist, ListPtr dihedlist, int a1, int a2, int a3, int a4);
extern ListPtr	MakeBondList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeAngleList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeDihedList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	MakeImproperList (ChainPtr clist, ZmatPtr zmatlist, ListPtr dihedlist);
extern int	MakeConnect (ChainPtr clist, ZmatPtr zmatlist, ListPtr *listbnd, ListPtr *listang, ListPtr *listdih, ListPtr *listimp);
extern ListPtr	GetAddBondList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	GetAddAngleList (ChainPtr clist, ZmatPtr zmatlist);
extern ListPtr	GetAddDihedList (ChainPtr clist, ZmatPtr zmatlist);

#endif
