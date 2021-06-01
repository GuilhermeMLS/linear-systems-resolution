#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

SistLinear_t *copyLinearSystem(SistLinear_t *linearSystem) {
    SistLinear_t *copy = alocaSistLinear(linearSystem->n);
    // Copy matrix A
    for (int i = 0; i < linearSystem->n; i++) {
        for (int j = 0; j < linearSystem->n; j++) {
            copy->A[i][j] = linearSystem->A[i][j];
        }
    }
    // Copy b
    for (int i = 0; i < linearSystem->n; i++) {
        copy->b[i] = linearSystem->b[i];
    }
    // Copy erro and n
    copy->erro = linearSystem->erro;
    copy->n = linearSystem->n;
    return copy;
}


void printSolution(
    char* method,
    double executionTime,
    real_t* solutionArray,
    int linearSystemSize,
    real_t residueL2Norm,
    int numberOfIterations
) {
    if (numberOfIterations > 0) {
        printf("===> %s: %lf ms --> %d iteracoes\n", method, executionTime, numberOfIterations);
    } else {
        printf("===> %s: %lf ms\n", method, executionTime);
    }
    printf("  --> X: ");
    prnVetor(solutionArray, linearSystemSize);
    printf("\n");
    printf("  --> Norma L2 do resíduo: %f\n", residueL2Norm);
    printf("\n");
}

void refineSolution(SistLinear_t* linearSystem, real_t* solutionArray) {
    SistLinear_t *refinementLinearSystem = copyLinearSystem(linearSystem);
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* refinedSolution = malloc(refinementLinearSystem->n * sizeof(real_t));
    int numberOfIterations = refinamento(refinementLinearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(refinementLinearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(refinementLinearSystem, refinedSolution, residueArray);
    printSolution("Refinamento", *executionTime, refinedSolution, refinementLinearSystem->n, residueL2Norm, numberOfIterations);
}

void gaussianElimination(SistLinear_t* linearSystem) {
    SistLinear_t *gaussEliminationLinearSystem = copyLinearSystem(linearSystem);
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(gaussEliminationLinearSystem->n * sizeof(real_t));
    eliminacaoGauss(gaussEliminationLinearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(gaussEliminationLinearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(gaussEliminationLinearSystem, solutionArray, residueArray);
    printSolution("Eliminação Gauss", *executionTime, solutionArray, gaussEliminationLinearSystem->n, residueL2Norm, 0);
    if (residueL2Norm > 5.0) {
        refineSolution(gaussEliminationLinearSystem, solutionArray);
    }
}

void jacobiMethod(SistLinear_t* linearSystem) {
    SistLinear_t *jacobiMethodLinearSystem = copyLinearSystem(linearSystem);
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(jacobiMethodLinearSystem->n * sizeof(real_t));
    int numberOfIterations = gaussJacobi(jacobiMethodLinearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(jacobiMethodLinearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(jacobiMethodLinearSystem, solutionArray, residueArray);
    printSolution("Jacobi", *executionTime, solutionArray, jacobiMethodLinearSystem->n, residueL2Norm, numberOfIterations);
    if (residueL2Norm > 5.0) {
        refineSolution(jacobiMethodLinearSystem, solutionArray);
    }
}

void gaussSeidelMethod(SistLinear_t* linearSystem) {
    SistLinear_t *gaussSeidelMethodLinearSystem = copyLinearSystem(linearSystem);
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(gaussSeidelMethodLinearSystem->n * sizeof(real_t));
    int numberOfIterations = gaussSeidel(gaussSeidelMethodLinearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(gaussSeidelMethodLinearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(gaussSeidelMethodLinearSystem, solutionArray, residueArray);
    printSolution("Gauss-Seidel", *executionTime, solutionArray, gaussSeidelMethodLinearSystem->n, residueL2Norm, numberOfIterations);
    if (residueL2Norm > 5.0) {
        refineSolution(gaussSeidelMethodLinearSystem, solutionArray);
    }
}

void gaussianEliminationWithRefinementMethod(SistLinear_t* linearSystem) {
    SistLinear_t *gaussianEliminationWithRefinementLinearSystem = copyLinearSystem(linearSystem);
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(gaussianEliminationWithRefinementLinearSystem->n * sizeof(real_t));
    int numberOfIterations = refinamento(gaussianEliminationWithRefinementLinearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(gaussianEliminationWithRefinementLinearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(gaussianEliminationWithRefinementLinearSystem, solutionArray, residueArray);
    printSolution("Refinamento", *executionTime, solutionArray, gaussianEliminationWithRefinementLinearSystem->n, residueL2Norm, numberOfIterations);
    if (residueL2Norm > 5.0) {
        refineSolution(gaussianEliminationWithRefinementLinearSystem, solutionArray);
    }
}

int hasLinearSystemToRead(int *numberOfElements) {
    return scanf("%d", numberOfElements) != EOF;
}

int main() {
    int numberOfElements;
    int linearSystemsCounter = 1;
    SistLinear_t *linearSystem;
    while (hasLinearSystemToRead(&numberOfElements)) {
        scanf("%*c"); // ignores \n
        printf("***** Sistema %d ", linearSystemsCounter);
        linearSystem = lerSistLinear(numberOfElements);
        gaussianElimination(linearSystem);
        jacobiMethod(linearSystem);
        gaussSeidelMethod(linearSystem);
        gaussianEliminationWithRefinementMethod(linearSystem);
        linearSystemsCounter++;
    }
}