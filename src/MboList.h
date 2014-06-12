#ifndef MBO_LIST_H
#define MBO_LIST_H

typedef struct MboList *MboList;

MboList mboListCreate();
void mboListDestroy(MboList *);
void mboListCons(void *, MboList *);
void *mboListHead(MboList);
MboList mboListTail(MboList);
void mboListMap(MboList, void *(void*, void*), void*);

#endif

