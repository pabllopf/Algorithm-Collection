/* Funciones para la resolución de un sistema matricial por el método de Crout */

#include "crout.h"
#include "lapack.h"
#include <stdlib.h>
#include <stdio.h>

/** la función crout_descomposicion() calcula los vectores l[],m[] y u[] del
metodo de Crout visto en Clase a partir de una matriz tridiagonal dada por los
vectores a[] (la diagonal de la matriz) y  b[] y c[] (las codiagonales de la matriz)
la función devuelve 0 si termina correctamente y un número negativo en caso
contrario */
int crout_descomposicion(Array1D< real > &a,Array1D< real > &b,Array1D< real > &c,Array1D< real > &l,Array1D< real > &m,Array1D< real > &u)
{
  /** COMPROBAMOS DIMENSIONES */
  if(l.dim()!=a.dim() || l.dim()!=(m.dim()+1) || u.dim()!=m.dim() || u.dim()!=b.dim() || u.dim()!=c.dim()) return(-1);

  /** HACER ALUMNO */

    int N = a.dim();
    l[0] = a[0];

    if(l[0] == 0){
        return (-1);
    }
    u[0] = b[0]/l[0];

    for(int i = 1; i <= N-2; i++){
        m[i-1] = c[i-1];
        l[i] = a[i]-m[i-1]*u[i-1];
        if(l[i] == 0){
            return (-1);
        }
        u[i] = b[i]/l[i];
    }
    m[N-2] = c[N-2];
    l[N-1] = a[N-1] - m[N-2]*u[N-2];
    return (0);
}

/** la función crout_descenso() resuelve un sistema triangular inferior usando
una matriz triangular como la que sale al descomponer la matriz por el método de Crout
el vector t[] es el término independiente del sistema. La función devuelve el vector solución
devuelve un vector vacío si termina mal */
Array1D< real > crout_descenso (Array1D< real > &l,Array1D< real > &m,Array1D< real > &t)
{
   // COMPROBAMOS DIMENSIONES
  if(t.dim()!=l.dim() || t.dim()!=(m.dim()+1)) return(Array1D< real >());

  /** HACER ALUMNO */
    Array1D< real > z(l.dim());
    int N = l.dim();

    z[0] = t[0]/l[0];

    for(int i = 1; i < N; i++){

        z[i] = t[i] - (m[i-1]*z[i-1]);
        z[i] = z[i]/l[i];

    }

  return (z);


}

/** la función crout_remonte() resuelve un sistema triangular superior usando
una matriz triangular como la que sale al descomponer la matriz por el método de Crout
el vector z[] es el término independiente del sistema. La función devuelve el vector solución
devuelve un vector vacío si termina mal*/
Array1D< real > crout_remonte(Array1D< real > &u,Array1D< real > &z){
  // COMPROBAMOS DIMENSIONES
  if(z.dim()!=(u.dim()+1)) return(Array1D< real >());

  /** HACER ALUMNO */
   int N=z.dim();
   Array1D< real > v(N);

   // INICIAMOS EL ASCENSO
      v[N-1]=z[N-1]; // inicializamos la solucion

      for (int i=N-2; i>=0; i--){

         v[i]=z[i] - u[i]*v[i+1];

      }
  return(v);


}

/** la función crout_metodo_completo() resuelve un sistema tridiagonal usando el método de Crout
el vector t[] es el término independiente del sistema. La función devuelve el vector solución
utiliza las funciones anteriores para descomponer la matriz y hacer el remonte y el descenso.
Devuelve un vector vacío si termina mal */
Array1D< real > crout_metodo_completo(Array1D< real > &a,Array1D< real > &b,Array1D< real > &c,Array1D< real > &t)
{
  /** HACER ALUMNO */
  Array1D< real > l(a.dim());
  Array1D< real > m(b.dim());
  Array1D< real > u(c.dim());

  int aux = crout_descomposicion(a,b,c,l,m,u);
  if (aux == -1){return Array1D< real > ();}
  Array1D< real > z = crout_descenso(l,m,t);
  Array1D< real > v = crout_remonte(u,z);

  return (v);

}

