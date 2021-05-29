#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main()
{
    SistLinear_t* linearSystem = lerSistLinear();
    prnSistLinear(linearSystem);
    double_t* gaussianEliminationTime = malloc(sizeof(double_t));
    *gaussianEliminationTime = timestamp();
    real_t* solutionArray = malloc(linearSystem->n * sizeof(real_t));
    eliminacaoGauss(linearSystem, solutionArray, gaussianEliminationTime);
    printf("Gaussian Elimination time: %lfms\n", *gaussianEliminationTime);
}

