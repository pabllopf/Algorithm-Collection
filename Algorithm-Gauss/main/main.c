/* Pablo Perdomo Falcon 04-12-2019 */
#include <stdio.h>
#include "../lapack.h"
#include <stdlib.h>

int main()
{
  Array2D< real > A, A2;
  Array1D< real > b, b2,u,u2;
  int N;

  N=mn_leer_matriz("../data/A_gauss_1.txt",A);
  N=mn_leer_vector("../data/b_gauss_1.txt",b);
  N=mn_leer_matriz("../data/A_gauss_s.txt",A2);
  N=mn_leer_vector("../data/b_gauss_s.txt",b2);
  N=mn_leer_vector("../data/u_gauss_s.txt",u2);

  printf("SYSTEM MATRIX A\n");
  A.print("A");
  printf("VECTOR OF INDEPENDENT TERMS b\n");
  b.print("b");

  printf("RESULTS TO BE OBTAINED \n\n");
  A2.print("A'");
  b2.print("b'");
  printf("SOLUTION:\n");
  u2.print("u_correcta");

  printf("GAUSS METHOD RESULTS \n\n");
  u=mn_gauss(A,b);
  printf("SOLUTION:\n");
  u.print("u_gauss");


  system("pause");
  return 0;
}



