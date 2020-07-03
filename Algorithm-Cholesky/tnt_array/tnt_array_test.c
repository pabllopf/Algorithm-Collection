/*==========================================================================
       CREACION Y ESCRITURA EN DISCO DE LAS MATRICES QUE SE UTILIZAN EN
       LAS PRACTICA DE AN
  ==========================================================================*/

// INCLUSION DE LA LIBRERIA STANDARD PARA GESTIONAR ENTRADA-SALIDA
#include <stdio.h>
#include <stdlib.h>

// INCLUSION DE LA LIBRERIA PARA GESTIONAR LA ARITMETICA Y LOS ARRAYS
#include "../mn_aritmeticas.h"
#include "../mn_lapack.h"

main(){

  // DEFINIMOS DOS VECTORES REALES DE TAMANO 5
  Array1D< real > a(5),b(5,1.); // el segundo lo inicializamos a 1.

  // RELLENAMOS EL PRIMER VECTOR
  for(int i=0;i<a.dim();i++)  a[i]=(real) i;

  // HACEMOS LA OPERACION c[i]=a[i]-b[i]*2;
  Array1D< real > c=a-b*2;

  // IMPRIMIMOS EL RESULTADO
  c.print("c");

  // HACEMOS UNA COPIA DEL RESULTADO
  Array1D< real > d=c.copy();

  // LEEMOS DE DISCO UN ARRAY
  int dimension=mn_leer_vector("../datos/b_cholesky_1.txt",b);
  if(dimension<=0){ system("pause");  exit(1); }

  // IMPRIMIMOS EL RESULTADO
  b.print("b");


  // MANEJO MATRICES
  // LEEMOS DE DISCO UNA MATRIZ
  Array2D< real >  A;
  dimension=mn_leer_matriz("../datos/A_cholesky_1.txt",A);
  if(dimension<=0){ system("pause");  exit(1); }

  // IMPRIMIMOS A
  A.print("A");

  // CREAMOS UNA COPIA DE A
  Array2D< real >  B=A.copy();

  // CALCULAMOS LA INVERSA DE A
  Array2D< real > A_1=mn_inversa(A);

  // MULTIPLICAMOS A*A_1
  Array2D< real >  I=A_1*A;

  // IMPRIMIMOS EL RESULTADO
  I.print("I");

  // RESOLVEMOS UN SISTEMA POR EL METODO DE GAUSS
  Array1D< real > u=mn_gauss(A,b);

  // COMPROBAMOS QUE e=Au-b=0
  Array1D< real > e=A*u-b;
  e.print("e");

  // CALCULAMOS LA MULTIPLICACION DE UNA MATRIZ POR UN VECTOR
  Array1D< real > g=A*b;

  // CALCULAMOS LA TRASPUESTA DE A
   Array2D< real > A_T=A.transpose();

  // SALVAMOS EN DISCO UNA MATRIZ Y UN VECTOR
  mn_escribir_vector("vector.txt",g);
  mn_escribir_matriz("matriz.txt",A_T);

  // OPERACIONES DE SUMAR/RESTAR MATRICES y MULTPLICAR/SUMAR ESCALARES
  Array2D< real > C(2,2,(real) 1);
  Array2D< real > D(2,2,(real) 2);
  Array2D< real > E=C+(real) 1 + D*((real) 0.5);
  E.print("E");

  // HACEMOS UNA PAUSA EN EL PROGRAMA PARA PODER EXAMINAR EL RESULTADO POR CONSOLA
  system("pause");


}


