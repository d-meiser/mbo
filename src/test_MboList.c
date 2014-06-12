#include <MboList.h>

#include <TestUtils.h>

int testMboListCreate()
{
	int errs = 0;
	MboList l;

	l = mboListCreate();
	mboListDestroy(&l);
	CHK_EQUAL(l, 0, errs);
	return errs;
}

int testMboListCons()
{
	int errs = 0;
	MboList l, t;
	void *h;
	int data1, data2;

	l = mboListCreate();

	data1 = 10;
	mboListCons(&data1, &l);
	h = mboListHead(l);
	CHK_EQUAL(*(int*)h, 10, errs);
	t = mboListTail(l);
	CHK_EQUAL(t, 0, errs);

	data2 = 13;
	mboListCons(&data2, &l);
	h = mboListHead(l);
	CHK_EQUAL(*(int*)h, 13, errs);
	t = mboListTail(l);
	h = mboListHead(t);
	CHK_EQUAL(*(int*)h, 10, errs);
	t = mboListTail(t);
	CHK_EQUAL(t, 0, errs);

	mboListDestroy(&l);
	return errs;
}

int testMboListHead()
{
	int errs = 0;
	MboList l;
	void *h;

	l = mboListCreate();
	h = mboListHead(l);
	CHK_EQUAL(h, 0, errs);
	mboListDestroy(&l);
	return errs;
}

int testMboListTail()
{
	int errs = 0;
	MboList l, t;

	l = mboListCreate();
	t = mboListTail(l);
	CHK_EQUAL(t, 0, errs);
	mboListDestroy(&l);
	return errs;
}

int testMboListMap()
{
	int errs = 0;
	MboList l;

	l = mboListCreate();
	mboListDestroy(&l);
	return errs;
}

int main()
{
	int errs = 0;
	errs += testMboListCreate();
	errs += testMboListCons();
	errs += testMboListHead();
	errs += testMboListTail();
	errs += testMboListMap();
	return errs;
}

