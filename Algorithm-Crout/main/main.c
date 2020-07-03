/* Pablo Perdomo Falcón 03-12-2019 */
#include <stdio.h>
#include "../crout.h"
#include "../lapack.h"
#include <stdlib.h>

int main()
{
  int N;
  {
    printf("EXAMPLE 1 CROUT METHOD: \n\n");
    Array1D< real > a,b,c,t,l(3,0.),m(2,0.),u(2,0.);
    N=an_leer_vector("../data/a1.txt",a);
    N=an_leer_vector("../data/b1.txt",b);
    N=an_leer_vector("../data/c1.txt",c);
    N=an_leer_vector("../data/t1.txt",t);

    printf("We check the crout decomposition function()\n");
    int error=crout_descomposicion(a,b,c,l,m,u);
    printf("vector l[] calculated : (%lf,%lf,%lf) \n",l[0],l[1],l[2]);
    printf("vector l[] correct : (%lf,%lf,%lf) \n\n",2.,1.5,4./3.);
    printf("vector m[] calculated : (%lf,%lf) \n",m[0],m[1]);
    printf("vector m[] correct : (%lf,%lf) \n\n",-1.,-1.);
    printf("vector u[] calculated : (%lf,%lf) \n",u[0],u[1]);
    printf("vector u[] correct : (%lf,%lf) \n\n\n",-0.5,-2./3.);

    printf("We check the crout descent function()\n");
    Array1D< real > z=crout_descenso (l,m,t);
    if(z.dim()==3){
      printf("vector z[] calculated: (%lf,%lf,%lf) \n",z[0],z[1],z[2]);
    }
    else { printf("error dimension vector\n"); }
    printf("vector z[] correct : (%lf,%lf,%lf) \n\n\n",0.5,1./3.,1.);

    printf("We check the full method crout function()\n");
    Array1D< real > v=crout_metodo_completo(a,b,c,t);
    if(v.dim()==3){
      printf("solution v[] calculated : (%lf,%lf,%lf) \n",v[0],v[1],v[2]);
    }
    else { printf("error dimension vector\n"); }
    printf("solution v[] correct : (%lf,%lf,%lf) \n\n\n",1.,1.,1.);
    system("pause");
  }

  {
    printf("EXAMPLE 2 CROUT METHOD:\n\n");
    Array1D< real > a,b,c,t,l(3,0.),m(2,0.),u(2,0.);
    N=an_leer_vector("../data/a2.txt",a);
    N=an_leer_vector("../data/b2.txt",b);
    N=an_leer_vector("../data/c2.txt",c);
    N=an_leer_vector("../data/t2.txt",t);

    printf("We check the crout descent function()\n");
    int error=crout_descomposicion(a,b,c,l,m,u);
    printf("vector l[] calculated : (%lf,%lf,%lf) \n",l[0],l[1],l[2]);
    printf("vector l[] correct : (%lf,%lf,%lf) \n\n",1.,2.,3.);
    printf("vector m[] calculated : (%lf,%lf) \n",m[0],m[1]);
    printf("vector m[] correct : (%lf,%lf) \n\n",4.,5.);
    printf("vector u[] calculated : (%lf,%lf) \n",u[0],u[1]);
    printf("vector u[] correct : (%lf,%lf) \n\n\n",6.,7.);

    printf("We check the crout descent function()\n");
    Array1D< real > z=crout_descenso (l,m,t);
    if(z.dim()==3){
      printf("vector z[] calculated : (%lf,%lf,%lf) \n",z[0],z[1],z[2]);
    }
    else { printf("error dimension vector\n"); }
    printf("vector z[] correct : (%lf,%lf,%lf) \n\n\n",15.,9.,1.);

    printf("We check the full method crout function()\n");
    Array1D< real > v=crout_metodo_completo(a,b,c,t);
    if(v.dim()==3){
      printf("solution v[] calculated : (%lf,%lf,%lf) \n",v[0],v[1],v[2]);
    }
    else { printf("error dimension vector\n"); }
    printf("solution v[] correct : (%lf,%lf,%lf) \n\n\n",3.,2.,1.);
    system("pause");
  }

  return 0;
}
