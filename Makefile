CC=/usr/gcc_trunk/bin/gcc
CFLAGS+=-O0 -g -Wall -fdiagnostics-color=auto -std=c99 -fsanitize=address

all: BUILD_TESTS

test_ElemOp: ElemOpImpl.o

test_ProdSpace: ProdSpaceImpl.o

test_TensProdOp: TensProdOpImpl.o

TEST_EXECUTABLES=test_ElemOp \
		 test_ProdSpace \
		 test_TensProdOp

BUILD_TESTS: ${TEST_EXECUTABLES}

.PHONY: clean
clean:
	rm -rf ${TEST_EXECUTABLES} *.o
