/*========================================================================
  FUNCIONES RELACIONADAS CON LA RESOLUCION DE SISTEMAS DE ECUACIONES
  ======================================================================== */

// INCLUSION DE LIBRERIAS NECESARIAS
#include <stdio.h>
#include "lapack.h"

// FUNCION PARA CALCULAR LA INVERSA DE UNA MATRIZ POR EL METODO DE GAUSS
Array2D< real > mn_inversa(
  const Array2D< real > &A  /* MATRIZ */)

{
  int N=A.dim1();
  if(N!=A.dim2()) return Array2D<real>();

  real max,paso,mul;
  int i,j,i_max,k;
  Array2D< real > PASO(N,2*N);
  Array2D< real > A_1(N,N);

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      PASO[i][j]=A[i][j];
      PASO[i][j+N]=0.;
    }
  }
  for(i=0;i<N;i++)
      PASO[i][i+N]=1.;

  for(i=0;i<N;i++){
    max=mn_abs(PASO[i][i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(mn_abs(PASO[j][i])>max){
         i_max=j; max=mn_abs(PASO[j][i]);
       }
    }

    if(max==0.){
      printf("Matriz No Invertible\n");
      return Array2D<real>();
    }
    if(i_max>i){
      for(k=0;k<2*N;k++){
        paso=PASO[i][k];
        PASO[i][k]=PASO[i_max][k];
        PASO[i_max][k]=paso;
      }
    }

    for(j=i+1;j<N;j++){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<2*N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  if(mn_abs(PASO[N-1][N-1])==0.){
      printf("Matriz No-Invertible\n");
      return Array2D<real>();
    }

  for(i=N-1;i>0;i--){
    for(j=i-1;j>=0;j--){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<2*N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  for(i=0;i<N;i++)
    for(j=N;j<2*N;j++)
      PASO[i][j]/=PASO[i][i];

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      A_1[i][j]=PASO[i][j+N];

  return(A_1);
}


// FUNCION PARA RESOLVER SISTEMAS POR EL METODO DE GAUSS
Array1D< real > mn_gauss(
  const Array2D< real > &A  /* MATRIZ DEL SISTEMA */,
  const Array1D< real > &b) /* VECTOR DE TERMINOS INDEPENDINENTES */
{

  real max,paso,mul,suma;
  int i,j,i_max,k;

  //COMPROBAMOS DIMENSIONES
  if( A.dim1()!=A.dim2() || A.dim1()!=b.dim() ) return(Array1D< real >());
  int N=A.dim1();

  // DEFINIMOS E INICIALIZAMOS MATRIZ Y VECTORES AUXILIARES
  Array2D< real > A2=A.copy(); // hacemos una copia de A
  Array1D< real > b2=b.copy(); // hacemos una copia de b
  Array1D< real > u(N);
  Array1D< int > piv(N); // vector para gestionar el pivoteo
  for(int k=0;k<N;k++) piv[k]=k;

  // HACER ALUMNO

  // INICIAMOS EL PROCESO DE CONVERTIR EN CERO DE LA DIAGONAL HACIA ABAJO
  for(i=0;i<N;i++){
    // calculamos el máximo de la diagonal hacia abajo
    max=mn_abs(A2[piv[i]][i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(mn_abs(A2[piv[j]][i])>max){
         i_max=j;
         max=mn_abs(A2[piv[j]][i]);
       }
    }
    if(max==0.){
      printf("Sistema no tiene Solucion 0\n");
      return(Array1D< real >());
    }
    // gestionamos el pivoteo
    if(i_max>i){
      paso=piv[i];
      piv[i]=piv[i_max];
      piv[i_max]=paso;
    }
    // combinamos las filas para hacer 0 de la diagonal hacia abajo
    for(j=i+1;j<N;j++){
      mul=-A2[piv[j]][i]/A2[piv[i]][i];
      A2[piv[j]][i]=0.;
      for(k=i+1;k<N;k++)
        A2[piv[j]][k]+=mul*A2[piv[i]][k];
      b2[piv[j]]+=mul*b2[piv[i]];
    }
  }
  if(mn_abs(A[piv[N-1]][N-1])==0.){
    printf("Sistema no tiene Solucion 1\n");
    return(Array1D< real >());
  }

  // HACEMOS  EL REMONTE PARA RESOLVER EL SISTEMA
  for(i=N-1;i>=0;i--){
    suma=b2[piv[i]];
    for(k=i+1;k<N;k++)
      suma-=A2[piv[i]][k]*u[k];
    u[i]=suma/A2[piv[i]][i];
  }

  // DEVOLVEMOS LA SOLUCIÓN
  return(u);
}


// FUNCION PARA CALCULAR EL ERROR DEL SISTEMA
real mn_error_sistema(const Array2D< real > &A, const Array1D< real > &u, const Array1D< real > &b)
{
   // HACER ALUMNO
   int i;
   if(b.dim()==0 || b.dim()!=u.dim() || b.dim()!=A.dim1() ||  b.dim()!=A.dim2()) return(0.);
   Array1D< real > e=A*u-b;

   real Sum=0.;
   for (i=0;i<b.dim();i++){
      Sum = Sum + mn_abs(e[i])/(mn_abs(b[i])+ 1.);
   }
   return (Sum/b.dim());
}


/* FUNCION PARA LEER UN VECTOR DE DISCO. RETORNA LA DIMENSION DEL VECTOR */
int mn_leer_vector(
  char *nombrefichero,
   Array1D< real > &vector)
{
  int dimension;
  float paso;
  FILE *f;

  // ABRIMOS EL FICHERO
  if(f=fopen( nombrefichero, "r"),!f){
    printf("Problema con la apertura del fichero\n");
    return -1;
  }

  // LEEMOS LA DIMENSION
  fscanf(f,"%d\n",&dimension);
  if(dimension<1) return(-2);

  // COGEMOS MEMORIA
   Array1D< real > v(dimension);

  // LEEMOS EL VECTOR
  for(int i=0;i<dimension;i++){
    fscanf(f,"%f\n",&paso);
    v[i]=paso;
  }
  fclose(f);
  vector=v.copy();
  return dimension;
}

/* FUNCION PARA ESCRIBIR UN VECTOR DE DISCO DE DIMENSION dimension Y LO ALMACENA EN vector */
int mn_escribir_vector(
  char *nombrefichero,
   Array1D< real > &vector)
{
  int i;
  FILE *f;
  int dimension=vector.dim();
  if(f=fopen( nombrefichero, "w"),!f){
    printf("Problema con la escritura del fichero\n");
    return 1;
  }
  fprintf(f,"%d\n",dimension);
  for(i=0;i<dimension;i++) fprintf(f,"%f\n",(float) vector[i]);
  fclose(f);
  return 0;
}

/* FUNCION PARA LEER UNA MATRIZ DE DISCO DE DIMENSION
dimension Y LO ALMACENA EN LA MATRIZ matriz  */
int mn_leer_matriz(
  char *nombrefichero,
  Array2D< real > &matriz)
{
  int dimension1,dimension2;
  float paso;
  FILE *f;
  if(f=fopen( nombrefichero, "r"),!f){
    printf("Problema con la apertura del fichero\n");
    return 1;
  }
  fscanf(f,"%d %d\n",&dimension1, &dimension2);
  if(dimension1<1 || dimension2<1) return(-1);

  // RESERVAMOS MEMORIA PARA LA MATRIZ
  Array2D< real > m(dimension1,dimension2);

  for(int i=0;i<dimension1;i++){
    for(int j=0;j<dimension2;j++){

       fscanf(f,"%f ",&paso);
       m[i][j]=paso;
       //printf("paso=%e\n",(double) m[i][j]);
    }
    fscanf(f,"\n");
  }
  fclose(f);
  matriz=m.copy();
  return dimension1;
}

/* FUNCION PARA ESCRIBIR UNA MATRIZ EN DISCO
dimension */
int mn_escribir_matriz(
  char *nombrefichero,
  Array2D< real > &matriz)
{
  int i,j;
  FILE *f;
  if(f=fopen( nombrefichero, "w"),!f){
    printf("Problema con la escritura del fichero\n");
    return 1;
  }
  fprintf(f,"%d %d\n",matriz.dim1(),matriz.dim2());
  for(i=0;i<matriz.dim1();i++){
    for(j=0;j<matriz.dim2();j++){
       fprintf(f,"%f ",(float) matriz[i][j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  return 0;
}






