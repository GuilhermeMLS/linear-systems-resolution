## Resolution of linear systems
 A program to solve linear systems through different methods: gaussian elimination, Gauss-Seidel method and Jacobi method.
 Project made to the Scientific Computation course at UFPR.
 
   > Notice: the headers and structs were given by the professor, that's why there are some parts of the project in Portuguese and with abbreviations instead of readable names.

### Compiling
```bash
$ make
```

### Running
```bash
$ ./labSisLin < sistemas.dat
```
Expected output example:

```bash
***** Sistema 0 --> n = 4, erro: 0.050000
===> Eliminação Gauss: 0.012939 ms
  --> X: 2.416667 2.583333 1.916667 4.083333 
  --> Norma L2 do resíduo: 0.000001

===> Jacobi: 0.021973 ms --> 5 iteracoes
  --> X: 2.416667 2.583333 1.916667 4.083333 
  --> Norma L2 do resíduo: 0.000001

===> Gauss-Seidel: 0.002197 ms --> 5 iteracoes
  --> X: 2.416667 2.583333 1.916667 4.083333 
  --> Norma L2 do resíduo: 0.000001

===> Refinamento: 0.001953 ms --> 1 iteracoes
  --> X: 2.416667 2.583333 1.916667 4.083333 
  --> Norma L2 do resíduo: 0.000001
```
