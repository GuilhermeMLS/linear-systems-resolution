#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

void printSolution(
        char* method,
        double executionTime,
        real_t* solutionArray,
        int linearSystemSize,
        real_t residueL2Norm
) {
    printf("===> %s: %lf ms\n", method, executionTime);
    printf("  --> X:");
    prnVetor(solutionArray, linearSystemSize);
    printf("\n");
    printf("  --> Norma L2 do resíduo: %f\n", residueL2Norm);
}

void gaussianElimination(SistLinear_t* linearSystem) {
    double_t* executionTime = malloc(sizeof(double_t));
    *executionTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    eliminacaoGauss(linearSystem, solutionArray, executionTime);
    char methodName[] = "Eliminação Gauss";
    printSolution(methodName, *executionTime, solutionArray, linearSystem->n, 1);
}

int main()
{
    SistLinear_t* linearSystem = lerSistLinear();
    prnSistLinear(linearSystem);
    gaussianElimination(linearSystem);
}