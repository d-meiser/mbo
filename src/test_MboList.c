#include <MboList.h>

#include <TestUtils.h>

int testMboListCreate()
{
	int errs = 0;
	MboList l;
	l = mboListCreate();
	return errs;
}

int main()
{
	int errs = 0;
	errs += testMboListCreate();
	return errs;
}

