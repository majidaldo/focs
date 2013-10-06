#include <stdio.h>

/*
this code is probably longer (more complicated?) but it's
flexible and I think more readable than having +/-1 i,j
all over the place
*/

typedef  float nt;// number type
typedef unsigned int it; //index type


//THE GRID

/*assuming grid is made from bottom left
2 234567
1 678901
0 012345
 
  012345
*/
typedef struct cartgrid{
  nt* u;
  it xn;  it yn;
  nt Dx;  nt Dy;} cg;
#define xiyi2i(xi,yi,xn) (yi)*(xn)+(xi)
//navigating the grid
nt *gt(cg *g, it *xi, it *yi
       ,int *Dxi, int *Dyi){return &(g->u)[xiyi2i(*xi+*Dxi,*yi+*Dyi,g->xn)];}
nt *up(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi,  *yi+1,g->xn)];}
nt *dn(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi,  *yi-1,g->xn)];}
nt *lf(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi-1,*yi  ,g->xn)];}
nt *rt(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi+1,*yi  ,g->xn)];}
nt *at(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi,  *yi  ,g->xn)];}
#undef xiyi2i



//DERIVATIVES

typedef void (*fp)(nt *, cg *, it *, it *);
void fwddx (nt *o, cg *g, it *xi, it *yi){*o=(*rt(g,xi,yi)-*at(g,xi,yi))/(  g->Dx);}
void bwddx (nt *o, cg *g, it *xi, it *yi){*o=(*at(g,xi,yi)-*lf(g,xi,yi))/(  g->Dx);}
void ctrdy (nt *o, cg *g, it *xi, it *yi){*o=(*up(g,xi,yi)-*dn(g,xi,yi))/(2*g->Dy);}
void ctrdx2(nt *o, cg *g, it *xi, it *yi){*o=(*rt(g,xi,yi)
					  -2*(*at(g,xi,yi))
					  +   *lf(g,xi,yi))/      ( (g->Dx)*(g->Dx));}
//                                                             C has no builtin pwr :/
void ctrdy2(nt *o, cg *g, it *xi, it *yi){*o=(*up(g,xi,yi)
					  -2*(*at(g,xi,yi))
					  +   *dn(g,xi,yi))/      ( (g->Dy)*(g->Dy));}

void del(nt *O, cg *g
	 , it *xi, it *yi
	 , fp dfx, fp dfy
	 ){
  (*dfx)(&O[0],g,xi,yi);
  (*dfy)(&O[1],g,xi,yi);} //C is complicated!
void del2(nt *o, cg *g
	  , it *xi, it *yi
	 , fp dfx2, fp dfy2
	 ){
  (*dfx2)(  o,g,xi,yi); nt oo=*o;
  (*dfy2)(&oo,g,xi,yi);    *o=*o+oo;}
void dot(nt *o, nt *X1, nt *X2){
  *o=(X1[0]*X2[0])
    +(X1[1]*X2[1]);}


//EQUATIONS

void v(nt *O, cg *g, it *xi, it *yi){ 
  O[0]=1;O[1]=0;}			 

void evalat(nt *o, cg *g
	    ,nt *Dt
	    ,it *xi, it *yi
	    ,nt *V
	    ,nt *k
	    ,fp dfx,  fp dfy
	    ,fp dfx2, fp dfy2){
  // dT/dt+v.delT=kdel2T   terms tr1+tr2=tr3
  nt       tr3; del2(&tr3   ,g,xi,yi,dfx2,dfy2); tr3=tr3*(*k);
  nt tr2tr2[2]; del ( tr2tr2,g,xi,yi,dfx ,dfy );
  nt       tr2; dot(&tr2,tr2tr2,V);
  *o=(tr3-tr2)*(*Dt)+*at(g,xi,yi);
  // printf("\n(%d,%d):%f=%f(%f-%f)+%f",*xi,*yi,*o,*Dt,tr3,tr2,*at(g,xi,yi));
       }


int main(){

  nt us[100]={0};
  cg T;
  T.u=us;
  T.xn=10; T.yn=10;
  T.Dx=10; T.Dy=10;
  it txi=5,tyi=5;
  *at(&T,&txi,&tyi)=40;
  printf("%f",*at(&T,&txi,&tyi));
  nt Dt=1;
  nt k=1;


  //nt d; fwddx(&d,&g,&txi,&tyi);
  nt V[2]={1,0};

  cg T2=T; nt u2[100]={0}; T2.u=u2;
  //int myarray[100]={0}; //sets all elements to 0
  //int myarray[100]={1}; //only sets first element to one. the rest are 0s(!!!)

  unsigned int xi;unsigned int yi;
  unsigned int i;
  for(i=0;i<3;i++){
    for(yi=1;yi<9;yi++){
      for(xi=1;xi<9;xi++){
	fp fpdx=&bwddx; //fp fpdy=&ctrdy;
	evalat(at(&T2,&xi,&yi),&T
	       ,&Dt
	       ,&xi,&yi
	       ,V
	       ,&k
	       ,fpdx, &ctrdy
	       ,&ctrdx2 ,&ctrdy2);
	printf("\n(%d,%d):%f->%f",xi,yi,*at(&T,&xi,&yi),*at(&T2,&xi,&yi));
      }
    }
    //printf("55o%f 55n%f",*at(&T,&txi,&tyi),*at(&T2,&txi,&tyi));
    nt *newT=T2.u, *oldT=T.u;
    T2.u=oldT; T.u=newT;//make newT old and viceversa
    printf("\n");
  }

}



//dT/dt+v.delT= k del2T
//      upwind  ctr diff



void sigma(float *o, float *v, float *dt, float *dx){
  *o=*v*(*dt)/(*dx);//dx 
}

  /*
void uwu_np1p(float *o		
	     , float *u_n_ij
	     , vec2 *v, float *dt, float *dx //sigma
	     ,float *u_n_im1j){
float sigx; sigma(&sigx, &(v->x), dt, dx);
float sigy; sigma(&sigy, &v->y, dt, dx);
*o=*u_n_ij-sigx*(*u_n_ij-*u_n_im1j)-sigy*(*u_n_ij-*u_n_ij
}
  */
