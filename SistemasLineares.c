#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res) {

}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t* SL, real_t* x, double* tTotal) {
//    TODO: Search for NULL equations
//    if (hasNullEquation(SL)) {
//        fprintf(stderr, "Null equation detected, the linear system has infinite solutions.\n");
//        return -1;
//    }
    // improving code readability
    SistLinear_t* linearSystem = SL;
    real_t* solutionsArray = x;
    int n = linearSystem->n;
    real_t** matrix = linearSystem->A;
    real_t* b = linearSystem->b;
    real_t m = 0;
    for (int k = 0; k < n-1; k++) {
        partialPivoting(linearSystem, k);
        for (int i = 1 + k; i < n; i++) {
            m = matrix[i][k] / matrix[k][k];
            for (int j = k; j < n; j++) {
                if (j == k) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] = matrix[i][j] - matrix[k][j] * m;
                }
            }
            b[i] = b[i] - b[k] * m;
            // TODO: remove this
            printf("\ndebug\n");
            prnSistLinear(linearSystem);
            printf("\n------\n");
        }
    }
    retrosubstitution(linearSystem, solutionsArray);
    // TODO: does the system got solved?
    // Calculate time
    *tTotal = timestamp() - *tTotal;
    return 0;
}

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal) {


}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal) {


}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal) {


}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear(unsigned int n) {
    unsigned int linearSystemSize = n; //improving code readability
    SistLinear_t* linearSystem = (SistLinear_t *) malloc(sizeof(SistLinear_t));
    if (!linearSystem) {
        return NULL;
    }
    linearSystem->A = (real_t**) malloc(linearSystemSize * sizeof(real_t*));
    for (int i = 0; i < linearSystemSize; i++) {
        linearSystem->A[i] = (real_t*) malloc(linearSystemSize * sizeof(real_t));
    }
    linearSystem->b = (real_t*) malloc(linearSystemSize * sizeof(real_t));
    linearSystem->n = linearSystemSize;
    if (!(linearSystem->A) || !(linearSystem->b)) {
        liberaSistLinear(linearSystem);
        return NULL;
    }
    return linearSystem;
}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear(SistLinear_t *SL) {
    for (int i = 0; i < SL->n; i++) {
        free(SL->A[i]);
    }
    free(SL->A);
    free(SL->b);
    free(SL);
}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t* lerSistLinear() {
    unsigned int numberOfElements;
    SistLinear_t* linearSystem;
    scanf("%d", &numberOfElements);
    linearSystem = alocaSistLinear(numberOfElements);
    scanf("%f", &linearSystem->erro);
    for(int i = 0; i < numberOfElements; ++i) {
        for(int j = 0; j < numberOfElements; ++j) {
            scanf("%f", &linearSystem->A[i][j]);
        }
    }
    for(int i = 0; i < numberOfElements; ++i) {
        scanf("%f", &linearSystem->b[i]);
    }
    return linearSystem;
}


// Exibe SL na saída padrão
void prnSistLinear(SistLinear_t* SL) {
    SistLinear_t* linearSystem = SL;
    int n = linearSystem->n;
    for (int i = 0; i < n; ++i) {
        printf("\n\t");
        for (int j = 0; j < n; ++j) {
            printf(" %.6f ", linearSystem->A[i][j]);
        }
        printf("   |   %.6f", linearSystem->b[i]);
    }
    printf("\n\n");
}

// Exibe um vetor na saída padrão
void prnVetor(real_t *v, unsigned int n) {
    int i;
    printf ("\n");
    for( i = 0; i < n; i++)
        printf("%f ", v[i]);
    printf("\n\n");
}

/*!
  \brief Partial pivot of Gaussian Elimination
*/
void partialPivoting(SistLinear_t *linearSystem, int currentIteration) {
    int linearSystemSize = linearSystem->n;
    real_t** matrix = linearSystem->A;
    real_t* independentTerms = linearSystem->b;

    // 1) Find the biggest element of the column below the pivot
    int column = currentIteration;
    int greaterLineIndex = currentIteration;
    real_t biggestElement = matrix[currentIteration][currentIteration];
    for (int i = currentIteration; i < linearSystemSize; i++) {
        if (fabs(matrix[i][column]) > fabs(biggestElement)) {
            greaterLineIndex = i;
        }
    }

    // 2) Replace Line A with Line B (Line B is the Line of the biggest element)
    int lineOne = currentIteration;
    int lineTwo = greaterLineIndex;

    // 2.1) Save lineOne in an aux array
    real_t *auxLine = malloc(sizeof(real_t) * linearSystemSize);
    for (int j = 0; j < linearSystemSize; j++) {
        auxLine[j] = matrix[lineOne][j];
    }

    // 2.2) Insert lineTwo in lineOne
    for (int j = 0; j < linearSystemSize; j++) {
        matrix[lineOne][j] = matrix[lineTwo][j];
    }

    // 2.3) Insert auxLine (lineOne) in lineTwo
    for (int j = 0; j < linearSystemSize; j++) {
        matrix[lineTwo][j] = auxLine[j];
    }

    // 2.3) Replace B
    real_t auxElement = independentTerms[lineOne];
    independentTerms[lineOne] = independentTerms[lineTwo];
    independentTerms[lineTwo] = auxElement;

    free(auxLine);
}

void retrosubstitution(SistLinear_t *linearSystem, real_t *solutionArray) {
    int n = linearSystem->n;
    real_t** matrix = linearSystem->A;
    real_t* b = linearSystem->b;
    solutionArray[n - 1] =  b[n - 1] / matrix[(n - 1)][(n - 1)];
    real_t sum = 0;
    for (int i = n - 2; i >= 0; i--) {
        sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= matrix[i][j] * solutionArray[j];
        }
        solutionArray[i] = sum / matrix[i][i];
    }
}