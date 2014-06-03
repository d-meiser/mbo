#ifndef PROD_SPACE_H
#define PROD_SPACE_H

struct ProdSpace;
typedef struct ProdSpace *ProdSpace;

ProdSpace CreateProdSpace(int);
void DestroyProdSpace(ProdSpace *);
void MultToProdSpace(ProdSpace, ProdSpace *);
ProdSpace CopyProdSpace(ProdSpace);
long long DimProdSpace(ProdSpace);

#endif
