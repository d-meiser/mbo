#include <MboListImpl.h>

MboList mboListCreate()
{
	return 0;
}

void mboListDestroy(MboList *list)
{
	MboList head = *list;
	MboList tail;
	while (head) {
		tail = head->tail;
		free(head);
		head = tail;
	}
	*list = 0;
}

void mboListCons(void *node, MboList *l)
{
	MboList h = malloc(sizeof(*h));
	h->head = node;
	h->tail = *l;
	*l = h;
}

void *mboListHead(MboList l)
{
	if (l) {
		return l->head;
	} else {
		return 0;
	}
}

MboList mboListTail(MboList l)
{
	if (l) {
		return l->tail;
	} else {
		return 0;
	}
}

void mboListMap(MboList l, void *f(void *ptr, void *ctx), void *ctx)
{
}
