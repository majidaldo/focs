#include <stdio.h>
#include <math.h>


//dx/dt=-aY+f(t)


void *Y_np1(float *o, float *Y_n, float *h_n, float *f){
  *o=*Y_n+*h_n*(*f);
}
void *ftx(float *o, float *t){*o=*t;}
void *gty(float *o, float *t){*o=*t;}
void *F(float *o, float *ab, float *Gr,float *XY, float *fg){
  *o=-(*ab)*(*Gr)*(*XY)+*fg;
    }

int main(){

 unsigned int i=0;
 float Y=100.0,X=100.0,t=0.0,f,g; //vars
 float h_n=.01,a=.2,b=.2; //consts a is effectiveness of Y, b is X
 float Cv=1.0; //conventional
 for(i=0;i<200;i++){
   if(X<=0){break;}//or could have it continue but have a floor of zero
   if(Y<=0){break;}

   printf("%f %f %f\n",t,X, Y);  

   t=t+h_n;

   ftx(&f,&t);
   F(&f,&a,&Cv,&Y,&f);
   gty(&g,&t);
   F(&g,&b,&Cv,&X,&g);

   Y_np1(&X,&X,&h_n,&f);//soln of 
   Y_np1(&Y,&Y,&h_n,&g);

  
 }

}
