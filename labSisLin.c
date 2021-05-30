#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

#define NUMBER_OF_SYSTEMS 14 // TODO: improve this

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

void gaussianElimination(SistLinear_t* linearSystem) {
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    eliminacaoGauss(linearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(linearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(linearSystem, solutionArray, residueArray);
    printSolution("Eliminação Gauss", *executionTime, solutionArray, linearSystem->n, residueL2Norm, 0);
}

void jacobiMethod(SistLinear_t* linearSystem) {
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    int numberOfIterations = gaussJacobi(linearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(linearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(linearSystem, solutionArray, residueArray);
    printSolution("Jacobi", *executionTime, solutionArray, linearSystem->n, residueL2Norm, numberOfIterations);
}

void gaussSeidelMethod(SistLinear_t* linearSystem) {
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    int numberOfIterations = gaussSeidel(linearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(linearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(linearSystem, solutionArray, residueArray);
    printSolution("Gauss-Seidel", *executionTime, solutionArray, linearSystem->n, residueL2Norm, numberOfIterations);
}

void gaussianEliminationWithRefinementMethod(SistLinear_t* linearSystem) {
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    int numberOfIterations = refinamento(linearSystem, solutionArray, executionTime);
    real_t* residueArray = malloc(linearSystem->n * sizeof(real_t));
    real_t residueL2Norm = normaL2Residuo(linearSystem, solutionArray, residueArray);
    printSolution("Refinamento", *executionTime, solutionArray, linearSystem->n, residueL2Norm, numberOfIterations);
}

int main() {
    for (int i = 0; i < NUMBER_OF_SYSTEMS; i++) {
        printf("***** Sistema %d ", i);
        SistLinear_t* linearSystem = lerSistLinear();
        gaussianElimination(linearSystem);
        jacobiMethod(linearSystem);
        gaussSeidelMethod(linearSystem);
        gaussianEliminationWithRefinementMethod(linearSystem);
    }
}