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
            if (matrix[k][k] == 0) {
                fprintf(stderr, "%s\n", "[Gaussian Elimination] Error: divison by zero");
                return -1;
            }
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
    int result = retrosubstitution(linearSystem, solutionsArray);
    if (result == -1) {
        fprintf(stderr, "%s\n", "[Gaussian Elimination Retrosubstitution] Error: divison by zero");
        return -1;
    }
    // TODO: does the system got solved?
    // Calculate time
    *tTotal = timestamp() - *tTotal;
    return 0;
}

real_t multiplyLinesForJacobiMethod(const real_t *solucao, SistLinear_t *SL, int i, unsigned int tam) {
    real_t soma = 0;
    for (int j = 0; j < tam; ++j) {
        if (j != i) {
            soma = soma + SL->A[i][j] * solucao[j];
        }
    }
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
    // "currentSolution" is the k vector
    real_t currentSolution[linearSystemSize], diff[linearSystemSize];
    int i, numberOfIterations, errorIncreaseCounter = 0;
    for (i = 0; i < linearSystemSize; ++i) {
        currentSolution[i] = 1;
    }
    numberOfIterations = 0;
    double currentEuclideanNorm, previousEuclideanNorm;
    do {
        for (i = 0; i < linearSystemSize; ++i) {
            if (!linearSystem->A[i][i]) {
                fprintf(stderr, "%s\n", "[Jacobi Method] Error: divison by zero");
                return -2;
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
               ? errorIncreaseCounter + 1
               : 0;
        }
        if (errorIncreaseCounter == CONVERGE_LIMIT) {
            fprintf(stderr, "[Jacobi Method] Error: solution dot not converge (the error has increased %d consecutive times)", CONVERGE_LIMIT);
            return -1;
        }
        if (++numberOfIterations == MAXIT) {
            fprintf(stderr, "%s\n", "[Jacobi Method] Error: maximum iterations number reached");
            return -2;
        }
        previousEuclideanNorm = currentEuclideanNorm;
    } while (currentEuclideanNorm > linearSystem->erro);
    for (int i = 0; i < linearSystemSize; i++) {
        if (isnan(solution[i])) {
            fprintf(stderr, "%s\n", "[Jacobi Method] Error: system has no solution");
            return -2;
        }
    }
    // Calculate time
    *tTotal = timestamp() - *tTotal;
    return numberOfIterations;
}

real_t multilyLinesForGaussSeidelMethod(SistLinear_t *SL, const real_t *x, int i, unsigned int tam) {
    real_t soma = 0;
    for (int j = 0; j < tam; j++) {
        if (j != i) {
            soma += SL->A[i][j] * x[j];
        }
    }
    return soma;
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
    // TODO: CALCULATE TIME
    unsigned int linearSystemSize = SL->n;
    int numberOfIterations, errorIncreaseCounter = 0;
    for (int i = 0; i < linearSystemSize; i++) {
        x[i] = 1;
    }
    numberOfIterations = 0;
    double currentEuclideanNorm, previousEuclideanNorm;
    real_t diff[linearSystemSize];
    do {
        for (int i = 0; i < linearSystemSize; i++) {
            if (!SL->A[i][i]) {
                fprintf(stderr, "%s\n", "[Gauss-Seidel] Error: divison by zero");
                return -2;
            }
            diff[i] = x[i];
            x[i] = (SL->b[i] - multilyLinesForGaussSeidelMethod(SL, x, i, linearSystemSize)) / SL -> A[i][i];
        }
        for (int i = 0; i < linearSystemSize; i++) {
            diff[i] = x[i] - diff[i];
        }
        currentEuclideanNorm = euclideanNorm(diff, linearSystemSize);
        if (numberOfIterations > 0) { // Ensures that "previousEuclideanNorm" has been initialized
            errorIncreaseCounter = (previousEuclideanNorm < currentEuclideanNorm)
                ? errorIncreaseCounter + 1
                : 0;
        }
        if (errorIncreaseCounter == CONVERGE_LIMIT) {
            fprintf(stderr, "[Gauss-Seidel] Error: solution dot not converge (the error has increased %d consecutive times)", CONVERGE_LIMIT);
            return -1;
        }
        if (++numberOfIterations == MAXIT) {
            fprintf(stderr, "%s\n", "[Gauss-Seidel] Error: maximum iterations number reached");
            return -2;
        }
        previousEuclideanNorm = currentEuclideanNorm;
    } while(currentEuclideanNorm > SL->erro);
    for (int i = 0; i < linearSystemSize; i++) {
        if (isnan(x[i])) {
            fprintf(stderr, "%s\n", "[Gauss-Seidel] Error: system has no solution");
            return -2;
        }
    }
    // Calculate time
    *tTotal = timestamp() - *tTotal;
    return numberOfIterations;
}

real_t distance(const real_t *vectorA, const real_t *vectorB, const int size) {
    real_t sum = 0;
    for (int i = 0; i < size; i++) {
        sum += (vectorA[i] - vectorB[i]);
    }
    return sqrt(sum);
}

/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial paa início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal) {
    // 1) Get an inicial solution x0 by solving Ax = b; assign i = 0
    // 2) Calculate the residue r = b - Ax(i) an test the stopping criterion (a)  (is the euclidean norm smaller than SL->erro?)
    // 3) Get the "w" vector by solving Aw = r;
    // 4) Get another solution x(i+1) = x(i) + w; test the stopping criterion (b) (is the maximum distance between xi and xi+1 smaller than SL->error? If yes, stop);
    // 5) increment i and go to step (2)
    SistLinear_t *linearSystem = SL;
    SistLinear_t *auxLinearSystem = (SistLinear_t *) malloc(sizeof(SistLinear_t));
    auxLinearSystem->A = linearSystem->A;
    auxLinearSystem->b = linearSystem->b;
    auxLinearSystem->erro = linearSystem->erro;
    auxLinearSystem->n= linearSystem->n;
    int linearSystemSize = linearSystem->n;
    double_t* gaussianEliminationExecutionTime = malloc(sizeof(double_t));
    *gaussianEliminationExecutionTime = timestamp();
    real_t *currentSolution = x;
    real_t *previousSolution = malloc(linearSystemSize * sizeof(real_t));
    real_t *w = malloc(linearSystemSize * sizeof(real_t));
    real_t *residue;
    real_t residueEuclideanNorm;
    real_t currentDistance;
    // 1)
    int iterationsCounter = 1;
    do {
        // 2)
        residue = getResidueArray(linearSystem, currentSolution);
        residueEuclideanNorm = euclideanNorm(residue, linearSystemSize);
        if (residueEuclideanNorm < linearSystem->erro) {
            break;
        }
        // 3)
        *auxLinearSystem->b = *residue;
        eliminacaoGauss(auxLinearSystem, w, gaussianEliminationExecutionTime);
        // 4)
        *previousSolution = *currentSolution;
        for (int i = 0; i < linearSystemSize; i++) {
            currentSolution[i] = currentSolution[i] + w[i];
        }
        currentDistance = distance(currentSolution, previousSolution, linearSystemSize);
        if (currentDistance < linearSystem->erro) {
            break;
        }
        iterationsCounter++;
    } while (iterationsCounter < MAXIT);
    *tTotal = timestamp() - *tTotal;
    return iterationsCounter;
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
SistLinear_t *lerSistLinear(int numberOfElements) {
    SistLinear_t *linearSystem;
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

int retrosubstitution(SistLinear_t *linearSystem, real_t *solutionArray) {
    int n = linearSystem->n;
    real_t **matrix = linearSystem->A;
    real_t *b = linearSystem->b;
    if (matrix[(n - 1)][(n - 1)] == 0) {
        return -1;
    }
    solutionArray[n - 1] = b[n - 1] / matrix[(n - 1)][(n - 1)];
    real_t sum = 0;
    for (int i = n - 2; i >= 0; i--) {
        sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= matrix[i][j] * solutionArray[j];
        }
        if (matrix[i][i] == 0) {
            return -1;
        }
        solutionArray[i] = sum / matrix[i][i];
    }
    return 0;
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
        residueArray[i] = linearSystem->b[i] - multiplyColumns(
                currentMatrixColumn,
                solutionArray,
                linearSystem->n
        );
    }
    return residueArray;
}
