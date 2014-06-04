#ifndef PROD_SPACE_H
#define PROD_SPACE_H

#ifdef __cplusplus
extern "C" {
#endif

struct ProdSpace;
typedef struct ProdSpace *ProdSpace;

ProdSpace CreateProdSpace(int);
void DestroyProdSpace(ProdSpace *);
void MultToProdSpace(ProdSpace, ProdSpace *);
ProdSpace CopyProdSpace(ProdSpace);
long long DimProdSpace(ProdSpace);

#ifdef __cplusplus
}
#endif
#endif
