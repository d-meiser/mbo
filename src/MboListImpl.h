#ifndef MBO_LIST_IMPL_H
#define MBO_LIST_IMPL_H

#include <MboList.h>

struct MboList {
	void *head;
	struct MboList *tail;
};

#endif

