/* Funciones para la resolución de un sistema matricial por el método de Cholesky */

#include "cholesky.h"
#include "lapack.h"
#include <stdlib.h>
#include <stdio.h>

/** Función que calcula las matrices de descomposición de Cholesky. Se implementa
   en una sola matriz, es decir se devuelve B como si fuera simétrica
   para después facilitar la resolución de los sistemas triangulares.
   Devuelve una matriz vacía si termina mal.
*/
Array2D< real > mn_cholesky_factorization(const Array2D< real > &A)
{
   if(A.dim1()!=A.dim2() ) return( Array2D< real >());

   /** HACER ALUMNO */
   int q,p,k;
   int N=A.dim1();

   Array2D<real> B(N,N);

   //Aplicamos el algoritmo

   for(q=0;q<=N-1;q++){

    real sum=0.;

    //Hacemos el sumatorio a parte
    for(k=0;k<=q-1;k++){
        sum=sum+(B[q][k]*B[q][k]);
    }

    if(A[q][q]<sum){
        return Array2D<real>();
    }

    B[q][q]=sqrtl((long double) A[q][q]-sum);

    for(p=q+1;p<=N-1;p++){ //Para p = q + 1,...;N - 1

        sum=0.;

        if(B[q][q]==0){
            return Array2D<real>();
        }

        //Hacemos el sumatorio a parte
        for(k=0;k<=q-1;k++){
            sum=sum+(B[p][k]*B[q][k]);
        }

        B[p][q]=(1/(B[q][q]))*((A[p][q])-sum);

        B[q][p]=B[p][q];

    }

   }
   return (B);
}

/** FUNCION PARA RESOLVER UN SISTEMA POR EL METODO DE CHOLESKY
    DEVUELVE UN VECTOR CON LA SOLUCION. DEBE USAR LA FUNCIÓN
    mn_cholesky_factorization(const Array2D< real > &A)
    SI TERMINA MAL DEVUELVE UN VECTOR VACIO*/
Array1D< real > mn_cholesky (const Array2D< real > &A, const Array1D< real > &b)
{

 /**HACER ALUMNO */
 Array2D<real> CH = mn_cholesky_factorization(A);

 if(CH.dim1()==0){
    return Array1D<real>();
 }

 Array1D<real> z = mn_descenso(CH,b);

 if(z.dim1()==0){
    return Array1D<real>();
 }

 Array1D<real> u = mn_remonte(CH,z);

 if(u.dim1()==0){
    return Array1D<real>();
 }

 return (u);


}

// FUNCION PARA RESOLVER UN SISTEMA TRIANGULAR INFERIOR
Array1D< real > mn_descenso (const Array2D< real > &B, const Array1D< real > &b)
{
   int i,j;
   if(B.dim1()!=B.dim2() || B.dim1()!=b.dim() ) return(Array1D< real >());
   int N=B.dim1();
   Array1D< real > z(N);

   // INICIAMOS EL DESCENSO
   for (i=0;i<N;i++){
      if (B[i][i]==0) { // comprobamos que la diagonal es distinto de cero
         return(Array1D< real >());
      }
      z[i]=b[i]; // inicializamos la solucion
      for (j=0;j<i;j++){ // aplicamos la formula para calcular la solucion
         z[i]=z[i] - B[i][j]*z[j];
      }
      z[i]=z[i]/B[i][i];
  }
  return(z);
}

// FUNCION PARA RESOLVER UN SISTEMA TRIANGULAR SUPERIOR
Array1D< real > mn_remonte (const Array2D< real > &B, const Array1D< real > &z)
{
   int i,j;
   if(B.dim1()!=B.dim2()  || B.dim1()!=z.dim() ) return(Array1D< real >());
   int N=B.dim1();
   Array1D< real > DU(N);


   // INICIAMOS EL REMONTE
   for (i=N-1;i>=0;i--){
      if (B[i][i]==0){  // comprobamos que la diagonal es distinto de cero
         return(Array1D< real >());
      }
      DU[i]=z[i]; // inicializamos la solucion
      for  (j=i+1;j<N;j++){ // aplicamos la formula para calcluar la solucion
         DU[i]=DU[i] - B[i][j]*DU[j];
      }
      DU[i]=DU[i]/B[i][i];
  }
  return(DU);
}




/** FUNCION PARA CALCULAR EL DETERMINANTE DE UNA MATRIZ POR LA FACTORIZACION DE CHOLESKY
// ESTA FUNCIÓN UTILIZA QUE SI A=B*B^T ENTONCES |A|=|B|*|B| Y QUE ADEMÁS EL DETERMINANTE
// DE UNA MATRIZ TRIANGULAR ES EL PRODUCTO DE LA DIAGONAL */
real mn_determinante_factorizacion_cholesky(const Array2D< real > &A)
{
  Array2D< real > CH = mn_cholesky_factorization (A); /*Calculamos la matriz de Cholesky */
  if (CH.dim1()==0) return (-1.);
  double determinante=1.;
  for(int i=0;i<CH.dim1();i++){
    determinante*=CH[i][i];
  }
  return(determinante*determinante);
}
