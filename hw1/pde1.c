#include <stdio.h>



//dx/dt=-aY+f(t)


//void *f1(float *f, float *t){return t;}


void *Y_np1(float *o, float *Y_n, float *h_n, float *t){
  *o=*Y_n+*h_n*(*t);
}


int main(){

  float t=4,xnp1=83,xn=33,hn=3333;

 float o;
 Y_np1(&o,&xn,&hn,&t);
printf("%f",o);


}
