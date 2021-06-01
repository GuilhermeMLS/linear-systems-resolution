#ifndef __SISLINEAR_H__
#define __SISLINEAR_H__

// Parâmetros para teste de convergência
#define MAXIT   50  // Número máximo de iterações em métodos iterativos
#define CONVERGE_LIMIT 50 // Maximum number of times the error can increase from one iteration to another

typedef float real_t;

typedef struct {
  unsigned int n; // tamanho do SL
  real_t erro; // critério de parada
  real_t **A; // coeficientes
  real_t *b; // termos independentes
} SistLinear_t;

// Alocaçao e desalocação de memória
SistLinear_t* alocaSistLinear (unsigned int n);
void liberaSistLinear (SistLinear_t *SL);

// Leitura e impressão de sistemas lineares
SistLinear_t *lerSistLinear (int numberOfElements);
void prnSistLinear (SistLinear_t *SL);
void prnVetor (real_t *vet, unsigned int n);

// Retorna a normaL2 do resíduo. Parâmetro 'res' deve ter o resíduo.
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res);

// Método da Eliminação de Gauss. Resultado no parâmetro 'x'
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal);

// Método de Jacobi. Valor inicial e resultado no parâmetro 'x' 
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal);

// Método de Gauss-Seidel. Valor inicial e resultado no parâmetro 'x' 
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal);

// Método de Refinamento. Valor inicial e resultado no parâmetro 'x'
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal);

void partialPivoting(SistLinear_t *linearSystem, int currentIteration);

int retrosubstitution(SistLinear_t *linearSystem, real_t *solutionArray);

real_t euclideanNorm(const real_t *vector, int size);

real_t* getResidueArray(SistLinear_t* linearSystem, real_t* solutionArray);

real_t multiplyColumns(const real_t *columnA, const real_t *columnB, int columnsSize);

real_t dotProduct(const real_t* vector, int size);

#endif // __SISLINEAR_H__

