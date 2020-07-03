/*========================================================================
  FUNCIONES RELACIONADAS CON LA RESOLUCION DE SISTEMAS DE ECUACIONES
  ======================================================================== */

// INCLUSION DE LIBRERIAS NECESARIAS
#include <stdio.h>
#include "lapack.h"

/** FUNCION PARA RESOLVER SISTEMAS POR EL METODO DE GAUSS CON PIVOTACIÓN
    DEVUELVE EL VECTOR SOLUCIÓN DEL SISTEMA SI TERMINA BIEN
    DEVUELVE UN VECTOR VACÍO EN CASO DE ERROR */
Array1D< real > mn_gauss(
  const Array2D< real > &A  /** MATRIZ DEL SISTEMA */,
  const Array1D< real > &b) /** VECTOR DE TERMINOS INDEPENDINENTES */
{
  /** COMPROBAMOS LAS DIMENSIONES DE LA MATRIZ Y EL VECTOR */
  if(A.dim1()!=A.dim2() || A.dim1()!=b.dim() || b.dim()==0) return Array1D<real>();

  /**  HACEMOS UNA COPIA DE A y b PARA MODIFICARLAS USANDO EL PROCESO DE GAUSS*/
  Array2D< real > A1=A.copy();
  Array1D< real > b1=b.copy();

  /** DECLARAMOS EL VECTOR SOLUCIÓN QUE SE DEVUELVE SI TERMINA BIEN */
  Array1D< real > u(b.dim()); // vector con la solución que se devuelve

  /** DECLARAMOS EL VECTOR PARA GESTIONAR EL PIVOTEO */
  Array1D<int> piv(b.dim());

  /** HACER ALUMNO */

    /** INICIALIZAMOS EL VECTOR DEL PIVOTE */

    int N = A.dim1();
    for(int k = 0; k < N; k++){
        piv[k] = k;
    }

 /** CONVERTIMOS EL SISTEMA EN UNO TRIANGULAR USANDO EL METODO
  DE GAUSS */

  for(int k = 0; k <= N-1; k++){

    /** DETECTAMOS EL MAXIMO DE LA DIAGONAL HACIA ABAJO */

    real max = fabs(A1[piv[k]][k]);
    int kmax = k;

    for(int j = k+1; j < N;j++){
        if(fabs(A1[piv[j]][k]) > max){
            max = fabs(A1[piv[j]][k]);
            kmax = j;
        }
    }


    /** CAMBIAMOS EL PIVOTE SI FUERA NECESARIO */

    if(kmax > k){
        int aux = piv[kmax];
        piv[kmax] = piv[k];
        piv[k] = aux;
    }

    /** PARA n=k+1 hasta N-1 hacer */
    N = A1.dim1();
    for(int n = k+1; n <= N-1;n++){
        real m = A1[piv[n]][k]  /  A1[piv[k]][k];
        A1[piv[n]][k] = 0.;

        /** Para q=k hasta N-1 hacer */
        for(int q = k+1; q <= N-1; q++){
            A1[piv[n]][q] -= m*  A1[piv[k]][q];
        }
        b1[piv[n]] -= m* b1[piv[k]];
    }


  }




  /** IMPRIMIR A1 y b1 PARA COMPROBAR QUE SON CORRECTOS */
  for(int i=0;i<b.dim();i++){
    for(int j=0;j<b.dim();j++){
      printf("A'%d%d=%1.0lf ",i,j,A1[piv[i]][j]);
    }
    printf("\n");
  }
  printf("\n");
  for(int j=0;j<b.dim();j++){
    printf("b'[%d]=%1.0lf ",j,b1[piv[j]]);
  }
  printf("\n\n");
  //system("pause");

  return(u);
}


// FUNCION PARA CALCULAR EL ERROR DEL SISTEMA
real mn_error_sistema(const Array2D< real > &A, const Array1D< real > &u, const Array1D< real > &b)
{
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












