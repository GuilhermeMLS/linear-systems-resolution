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
    res = getResidueArray(SL, x);
    return euclideanNorm(res, SL->n);
}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal) {
    //TODO: Search for NULL equations (?)
    //improving code readability
    SistLinear_t *linearSystem = SL;
    real_t *solutionsArray = x;
    int n = linearSystem->n;
    real_t **matrix = linearSystem->A;
    real_t *b = linearSystem->b;
    real_t m = 0;
    for (int k = 0; k < n - 1; k++) {
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
        }
    }
    retrosubstitution(linearSystem, solutionsArray);
    // TODO: does the system got solved?
    // Calculate time
    *tTotal = timestamp() - *tTotal;
    return 0;
}

real_t multiplyLinesForJacobiMethod(const real_t *solucao, SistLinear_t *SL, int i, unsigned int tam) {
    real_t soma = 0;
    for (int j = 0; j < tam; ++j)
        if (j != i)
            soma = soma + SL->A[i][j] * solucao[j];
    return soma;
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
int gaussJacobi(SistLinear_t *SL, real_t *x, double *tTotal) {
    SistLinear_t *linearSystem = SL;
    real_t *solution = x;
    unsigned int linearSystemSize = linearSystem->n;
    // linearSystem->x is the k + 1 vector
    // "currentSolution" is the k vector
    real_t currentSolution[linearSystemSize], diff[linearSystemSize];
    int i, numberOfIterations, errorIncreaseCounter = 0;
    // First Solution's elements are equal to 1 (as x0 is given)
    for (i = 0; i < linearSystemSize; ++i) {
        currentSolution[i] = 1;
    }
    numberOfIterations = 0;
    double currentEuclideanNorm, previousEuclideanNorm;
    do {
        for (i = 0; i < linearSystemSize; ++i) {
            if (!linearSystem->A[i][i]) {
                fprintf(stderr, "%s\n", "[Jacobi Method]: divison by zero");
                return 1;
            }
            solution[i] = (linearSystem->b[i] - multiplyLinesForJacobiMethod(currentSolution, linearSystem, i, linearSystemSize)) /
                   linearSystem->A[i][i];
        }
        for (i = 0; i < linearSystemSize; ++i) {
            diff[i] = solution[i] - currentSolution[i];
            currentSolution[i] = solution[i];
        }
        currentEuclideanNorm = euclideanNorm(diff, linearSystemSize);
        if (numberOfIterations > 0) { // Ensures that "previousEuclideanNorm" has been initialized
            errorIncreaseCounter = (previousEuclideanNorm < currentEuclideanNorm)
               ? ++errorIncreaseCounter
               : 0;
        }
        if (errorIncreaseCounter == CONVERGE_LIMIT) {
            fprintf(stderr, "[Jacobi Method] Solution dot not converge (the error has increased %d consecutive times)", CONVERGE_LIMIT);
            return -1;
        }
        if (++numberOfIterations == MAXIT) {
            fprintf(stderr, "%s\n", "[Jacobi Method] Maximum iterations number reached");
            return -2;
        }
        previousEuclideanNorm = currentEuclideanNorm;
    } while (currentEuclideanNorm > linearSystem->erro); // TODO: it was epsilon, is it correct?
    return numberOfIterations;
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
int gaussSeidel(SistLinear_t *SL, real_t *x, double *tTotal) {


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
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal) {


}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t *alocaSistLinear(unsigned int n) {
    unsigned int linearSystemSize = n; //improving code readability
    SistLinear_t *linearSystem = (SistLinear_t *) malloc(sizeof(SistLinear_t));
    if (!linearSystem) {
        return NULL;
    }
    linearSystem->A = (real_t **) malloc(linearSystemSize * sizeof(real_t *));
    for (int i = 0; i < linearSystemSize; i++) {
        linearSystem->A[i] = (real_t *) malloc(linearSystemSize * sizeof(real_t));
    }
    linearSystem->b = (real_t *) malloc(linearSystemSize * sizeof(real_t));
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
SistLinear_t *lerSistLinear() {
    unsigned int numberOfElements;
    SistLinear_t *linearSystem;
    scanf("%d", &numberOfElements);
    if (!numberOfElements || numberOfElements <= 0) {
        return NULL;
    }
    linearSystem = alocaSistLinear(numberOfElements);
    if (!linearSystem) {
        return NULL;
    }
    scanf("%f", &linearSystem->erro);
    for (int i = 0; i < numberOfElements; ++i) {
        for (int j = 0; j < numberOfElements; ++j) {
            scanf("%f", &linearSystem->A[i][j]);
        }
    }
    for (int i = 0; i < numberOfElements; ++i) {
        scanf("%f", &linearSystem->b[i]);
    }
    printf("--> n = %d, erro: %f\n", numberOfElements, linearSystem->erro);
    return linearSystem;
}


// Exibe SL na saída padrão
void prnSistLinear(SistLinear_t *SL) {
    SistLinear_t *linearSystem = SL;
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
    for (i = 0; i < n; i++)
        printf("%f ", v[i]);
}

/*!
  \brief Partial pivot of Gaussian Elimination
*/
void partialPivoting(SistLinear_t *linearSystem, int currentIteration) {
    int linearSystemSize = linearSystem->n;
    real_t **matrix = linearSystem->A;
    real_t *independentTerms = linearSystem->b;

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
    real_t **matrix = linearSystem->A;
    real_t *b = linearSystem->b;
    solutionArray[n - 1] = b[n - 1] / matrix[(n - 1)][(n - 1)];
    real_t sum = 0;
    for (int i = n - 2; i >= 0; i--) {
        sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= matrix[i][j] * solutionArray[j];
        }
        solutionArray[i] = sum / matrix[i][i];
    }
}

real_t euclideanNorm(const real_t *vector, int size) {
    return sqrt(dotProduct(vector, size));
}

real_t dotProduct(const real_t *vector, int size) {
    real_t dotProduct = 0;
    for (int i = 0; i < size; i++) {
        dotProduct += vector[i] * vector[i];
    }
    return dotProduct;
}


real_t multiplyColumns(const real_t *columnA, const real_t *columnB, int columnsSize) {
    real_t result = 0;
    for (int i = 0; i < columnsSize; i++) {
        result += columnA[i] * columnB[i];
    }
    return result;
}

real_t *getResidueArray(SistLinear_t *linearSystem, real_t *solutionArray) {
    // R = AX - b
    real_t *residueArray = malloc(linearSystem->n * sizeof(real_t));
    for (int i = 0; i < linearSystem->n; i++) {
        real_t *currentMatrixColumn = linearSystem->A[i];
        residueArray[i] = multiplyColumns(
                currentMatrixColumn,
                solutionArray,
                linearSystem->n
        ) - linearSystem->b[i];
    }
    return residueArray;
}
