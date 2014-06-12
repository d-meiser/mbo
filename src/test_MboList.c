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

void doubleItems(void *node, void *ctx)
{
	int *i = (int *)node;
	*i *= 2;
}

void accumulate(void *node, void *sum)
{
	*(int *)sum += *(int *)node;
}

int testMboListMap()
{
	int errs = 0;
	int data[3];
	int doublingResult[3];
	int i, sum, sumResult;
	MboList l, lIter;

	l = mboListCreate();

	data[0] = 10;
	mboListCons(data + 0, &l);
	data[1] = 13;
	mboListCons(data + 1, &l);
	data[2] = -3;
	mboListCons(data + 2, &l);

	mboListMap(l, &doubleItems, 0);
	doublingResult[0] = 2 * 10;
	doublingResult[1] = 2 * 13;
	doublingResult[2] = 2 * -3;
	i = 0;
	for (lIter = l; lIter != 0; lIter = mboListTail(lIter)) {
		CHK_EQUAL(*(int*)mboListHead(lIter), data[2 - i], errs);
		++i;
	}
	for (i = 0; i < 3; ++i) {
		CHK_EQUAL(doublingResult[i], data[i], errs);
	}

	sum = 0;
	mboListMap(l, &accumulate, &sum);
	sumResult = 0;
	for (i = 0; i < 3; ++i) {
		sumResult += data[i];
	}
	CHK_EQUAL(sumResult, sum, errs);

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

