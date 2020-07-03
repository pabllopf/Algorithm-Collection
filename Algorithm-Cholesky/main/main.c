/* Pablo Perdomo Falcón 01-12-2019*/
#include <stdio.h>
#include "../cholesky.h"
#include "../lapack.h"
#include <stdlib.h>

int main()
{
  Array2D< real > A;
  Array1D< real > b, b2;

  int i,N;
  printf("For the first example matrix, the decomposition results\n");
  printf("from Cholesky are on the class slides\n\n");

  for(int cont=2;cont<5;cont++){
    char nameA[200];
    char nameb[200];
    if(cont==1){
      sprintf(nameA,"../data/A_gauss_%d.txt",1);
      sprintf(nameb,"../data/b_gauss_%d.txt",1);
    }
    else if(cont<5){
      sprintf(nameA,"../data/A_cholesky_%d.txt",cont-1);
      sprintf(nameb,"../data/b_cholesky_%d.txt",cont-1);
    }
    else if(cont==5){
      sprintf(nameA,"../data/A_%d.txt",10);
      sprintf(nameb,"../data/b_%d.txt",10);
    }
    else if(cont==6){
      sprintf(nameA,"../data/A_%d.txt",100);
      sprintf(nameb,"../data/b_%d.txt",100);
    }
    else if(cont==7){
      sprintf(nameA,"../data/A_%d.txt",1000);
      sprintf(nameb,"../data/b_%d.txt",1000);
    }

    N=mn_leer_matriz(nameA,A);
    N=mn_leer_vector(nameb,b);

    if(N<10){
      A.print("A");
      b.print("b");
    }

    Array1D< real > u=mn_gauss(A,b);

    if (u.dim()>0){
      real error=mn_error_sistema(A,u,b);
      printf ("Gauss: The system error is: %e\n\n",(double) error);
      if(N<10){
        printf("The solution by the GAUSS method is\n");
        u.print("  u_gauss");
      }
    }
    else{
      printf("The system cannot be resolved by the GAUSS method\n");
    }

    Array1D< real > v=mn_cholesky(A,b);

    if (v.dim()>0)
    {
      real error=mn_error_sistema(A,v,b);
      printf ("Cholesky: The system error is: %e\n\n",(double) error);
      if(N<10){
        printf("The solution by the method of CHOLESKY is\n");
        v.print("  u_chole");
      }
    }

    else{
      printf("The system cannot be resolved by the CHOLESKY method\n");
    }

    printf ("----------------------------------------------\n");
    system("pause");
  }
    return 0;
}


