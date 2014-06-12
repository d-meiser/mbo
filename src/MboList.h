#ifndef MBO_LIST_H
#define MBO_LIST_H

typedef struct MboList *MboList;

MboList mboListCreate();
void mboListDestroy(MboList);
void *mboListHead(MboList);
void *mboListTail(MboList);
void mboListMap(MboList, void *(void*, void*), void*);

#endif

