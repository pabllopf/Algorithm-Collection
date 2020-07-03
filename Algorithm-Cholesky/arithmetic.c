/*========================================================================
  FUNCIONES PARA EL MANEJO DE ARITMETICAS DE PRECISION FINITA
  ======================================================================== */

// INCLUSION DE LA LIBRERIA PARA GESTIONAR LA ARITMETICA
#include "arithmetic.h"
#include <iostream>
using namespace std;

// FUNCION PARA CALCULAR LA PRECISION (UNIDAD DE REDONDEO) DE LA ARITMETICA
real mn_precision_aritmetica()
{
  // HACER ALUMNO
  real A=1.;
  real B=2.;
  int M=1;
  while(B>1.){
    A=A/2;
    B=1.+A;
    M=M+1;
  }
  return(2*A);
}

// FUNCION PARA CALCULAR EL MENOR NUMERO POSITIVO
int mn_menor_numero_positivo()
{
    // HACER ALUMNO
    real A = 1.;
    int M = 0;
    while (A > 0) {
      A = A/2;
      M =M+1;
    }
    return (M-1);
}

// FUNCION PARA CALCULAR EL MAYOR NUMERO POSITIVO
int mn_mayor_numero_positivo()
{
    // HACER ALUMNO
    real A = 1.;
    real B = 2.;
    int M = 0;
    while (B > A) {
      A = 2*A;
      B = 2*A;
      M = M+1;
    }
    return (M-1);
}

// FUNCION PARA CALCULAR LA DISTANCIA RELATIVA ENTRE 2 NUMEROS REALES
real mn_distancia(real A,real B)
{
  // HACER ALUMNO
  // PARA CALCULAR LA DISTANCIA SE DEBE UTILIZAR LA F�RMULA |B-A|/( max(|A|,|B|)+1. )
  // PARA CALCULAR EL VALOR ABSOLUTO SE PUEDE UTILIZAR LA MACRO mn_abs(x)
  // PARA CALCULAR EL M�XIMO DE 2 N�MEROS SE PUEDE UTILIZAR LA MACRO mn_max(A,B)
  return( mn_abs(B-A) / ( mn_max(mn_abs(A),mn_abs(B)) + 1. ) );
}

