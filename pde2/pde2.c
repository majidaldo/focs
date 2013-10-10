#include <stdio.h>
#include <stdlib.h> //for malloc()
#include <math.h> //for round()

/*
this code is probably longer (more complicated?) and not optimized.
 but it's flexible and I think more readable than having +/-1 i,j
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
  it xn;  it yn; //number of points
  nt Dx;  nt Dy;} cg;
#define xiyi2i(xi,yi,xn) (yi)*(xn)+(xi)
//navigating the grid
nt *gt(cg *g, it *xi, it *yi
       ,int *Dxi
       ,int *Dyi)            {return &(g->u)[xiyi2i(*xi+*Dxi,*yi+*Dyi,g->xn)];}
nt *up(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi     ,*yi+1   ,g->xn)];}
nt *dn(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi     ,*yi-1   ,g->xn)];}
nt *lf(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi-1   ,*yi     ,g->xn)];}
nt *rt(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi+1   ,*yi     ,g->xn)];}
nt *at(cg *g, it *xi, it *yi){return &(g->u)[xiyi2i(*xi     ,*yi     ,g->xn)];}
#undef xiyi2i



//DERIVATIVES

typedef void (*fp)(nt *, cg *, it *, it *);//eval something on the grid at a pt
void fwddx (nt *o, cg *g, it *xi, it *yi){*o=(*rt(g,xi,yi)-*at(g,xi,yi))/(  g->Dx);}
void bwddx (nt *o, cg *g, it *xi, it *yi){*o=(*at(g,xi,yi)-*lf(g,xi,yi))/(  g->Dx);}
void ctrdy (nt *o, cg *g, it *xi, it *yi){*o=(*up(g,xi,yi)-*dn(g,xi,yi))/(2*g->Dy);}
void ctrdx2(nt *o, cg *g, it *xi, it *yi){*o=(*rt(g,xi,yi)
					  -2*(*at(g,xi,yi))
					  +   *lf(g,xi,yi))/      ( (g->Dx)*(g->Dx));
//if(*xi==13 && *yi==2){printf("\ncdx2 %f %f %f %f %f",*o,*rt(g,xi,yi),*at(g,xi,yi),*lf(g,xi,yi),(g->Dx));}
}
//                                                             C has no builtin pwr :/
void ctrdy2(nt *o, cg *g, it *xi, it *yi){*o=(*up(g,xi,yi)
					  -2*(*at(g,xi,yi))
					  +   *dn(g,xi,yi))/      ( (g->Dy)*(g->Dy));
//if(*xi==13 && *yi==2){printf("\ncdy2 %f %f %f %f %f",*o,*up(g,xi,yi),*at(g,xi,yi),*dn(g,xi,yi),(g->Dy));}				  
					  }

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




//UTILS

void xiyi2xy(nt* O, cg *g, it *xi, it *yi){
  O[0]=g->Dx*(*xi);
  O[1]=g->Dy*(*yi);
}
void xy2xiyi(it* O, cg *g, nt *x, nt *y){
O[0]=round(*x/g->Dx);
O[1]=round(*y/g->Dy);
}




//EQUATIONS

//todo

//make v func ptr
typedef void(*vfp)(nt*, cg*, it*, it*, nt*);
void vxpara(nt *O, cg *g
	 , it *xi, it *yi
	 , nt *params){
/*
y being the dependent var
In[2]:= p=a x^2 + b x +c ==y
Out[2]= c+b x+a x^2==y
In[12]:= FullSimplify[Solve[{
p/.{y->0,x->xv+.5aa}
,p/.{y->0,x->xv-.5aa}
,p/.{y->ym,x->xv}}
,{a,b,c}]]
Out[12]{{
a->0. -(4. ym)/aa^2
,b->(8. xv ym)/aa^2
,c->(1. -(4. xv^2)/aa^2) ym}}
*/
#define vmax params[0]
#define a    params[1]
#define yv   params[2] 
  nt X[2]; xiyi2xy(X,g,xi,yi);
 
  nt aa2=(a)*(a);
  O[0]=      (-4/(aa2))     *(X[1]*X[1]) \
            +(8*yv/aa2)     * X[1]       \
            +(1-(4*yv*yv)/aa2);
  O[0]=O[0]*vmax;//todo elsewhere zero?
  O[1]=0;
#undef yv
#undef a
#undef vmax
}

void vxconst(nt *O, cg *g
	 , it *xi, it *yi
	 , nt *params){
  O[0]=params[0];O[1]=0;
}

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
  //printf("\n(%d,%d):%f",*xi,*yi,*o);
       }

//useless just put logic in the main loops
void evalover(it x1, it x2
	      ,it y1, it y2
	      ,cg *T, cg *T2
	      ,vfp v, nt *vparams
	      ,nt k, nt Dt
	      ,fp dfx, fp dfy
	      ,fp dfx2,fp dfy2){
  it xi,yi;
  nt V[2];
  for(yi=y1;yi<=y2;yi++){
    for(xi=x1;xi<=x2;xi++){
      (*v)(V,T,&xi,&yi,vparams);
	evalat(at(T2,&xi,&yi),T
	       ,&Dt
	       ,&xi,&yi
	       ,V
	       ,&k
	       ,dfx, dfy
	       ,dfx2 ,dfy2);      
    }
  }
}

cg makegrid(it xn, it yn
	    ,nt Lx, nt Ly){
  cg g;
  g.xn=xn;
  g.yn=yn;
  g.Dx=Lx/(xn-1);
  g.Dy=Ly/(yn-1);
  g.u=malloc(xn*yn*sizeof(nt));
  return g;
}






int main(){


  nt Lx=1.000,Ly=.500;
  it nx=30,ny=15;
  cg T=makegrid(nx,ny,Lx,Ly);

  nt T0=0;
  it Tit;
  //set to init temp
  for(Tit=0;Tit<(T.yn*T.xn);Tit++){T.u[Tit]=T0;}



  nt Dt=.01;
  nt k=.02;

  //needle coords.
#define sx .1*Lx
#define sy .1*Ly
#define xn .2*Lx
#define yn .4*Ly
  nt LL[2]={ xn-.5*sx //lower left
	     ,yn-.5*sy};
  nt UR[2]={ xn+.5*sx
	     ,yn+.5*sy};
  it LLnXi[2]; xy2xiyi(LLnXi,&T,&LL[0],&LL[1]);
  it URnXi[2]; xy2xiyi(URnXi,&T,&UR[0],&UR[1]);
#undef sx
#undef sy
#undef xn
#undef yn

  nt Tn=1;
  //set needle temp
  for(Tit =LLnXi[0];Tit <=URnXi[0]; Tit++){it Tit2;
    for(Tit2=LLnXi[1];Tit2<=URnXi[1]; Tit2++){
      *at(&T,&Tit,&Tit2)=Tn;}}

  printf("needle at (%d,%d) (%d,%d)\n",LLnXi[0],LLnXi[1],URnXi[0],URnXi[1]);


  //vessel coords
#define a  .2*Ly
#define yv .7*Ly
  LL[0]=0; //lower left
  LL[1]=yv-.5*a;
  UR[0]=Lx;
  UR[1]=yv+.5*a;
  it LLvXi[2]; xy2xiyi(LLvXi,&T,&LL[0],&LL[1]);
  it URvXi[2]; xy2xiyi(URvXi,&T,&UR[0],&UR[1]);

  nt vmax=1.5;
  nt vpp[3]={vmax,a,yv};
  nt vcp[1]={0};

#undef a
#undef yv

  printf("blood vessel at (%d,%d) (%d,%d)\n",LLvXi[0],LLvXi[1],URvXi[0],URvXi[1]);


  cg T2=makegrid(nx,ny,Lx,Ly);
  //int myarray[100]={0}; //sets all elements to 0
  //int myarray[100]={1}; //only sets first element to one. the rest are 0s(!!!)
 
  FILE *fo=fopen("t.txt","w");



  it xi,yi;
  unsigned int nf=3000;


  for(Tit=0;Tit<nf;Tit++){
    fprintf(fo,"FRAME\n");
    //printf(fo,"FRAME\n");

    fp afp;
    for(yi=0;yi<(ny);yi++){
      for(xi=0;xi<(nx);xi++){

	//apply boundary condition at vessel out
	if((xi==URvXi[0]) && (LLvXi[1]<=yi && yi<=URvXi[1])){
	  *at(&T2,&xi,&yi)=*lf(&T,&xi,&yi);//dT/dn=0
	  //fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	}
	//all other edges leave at init temp
	else if(xi==0 || xi==(nx-1) || yi==0 || yi==(ny-1)){
	  //fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	}
	//in vessel
	else if( (LLvXi[0]<=xi && xi<=(URvXi[0]-1)) && //the -1 shouldn't matter here
		 LLvXi[1]<=yi && yi<=URvXi[1]){
	  nt V[2];
	  vxpara(V,&T,&xi,&yi,vpp);
	  if(V[0]>0){afp=&bwddx;}
	  else{      afp=&fwddx;}
	  evalat(at(&T2,&xi,&yi),&T
		 ,&Dt
		 ,&xi,&yi
		 ,V
		 ,&k
		 ,afp, &ctrdy//vy=0 so this shouldn't matter
		 ,&ctrdx2 ,&ctrdy2);
	  //fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	}
	//in needle
	else if( (LLnXi[0]<=xi && xi<=(URnXi[0])) && 
		 LLnXi[1]<=yi && yi<=URnXi[1]){
	  //a waste of computing but it't not my priority
	  //set needle temp
	  *at(&T2,&xi,&yi)=Tn;
	  //fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	}
	//elsewhere (not needle or boundary or vessel
	else{
	  nt V[2]={0,0};
	  evalat(at(&T2,&xi,&yi),&T
		 ,&Dt
		 ,&xi,&yi
		 ,V
		 ,&k
		 ,&bwddx, &ctrdy //both shouldn't matter
		 ,&ctrdx2 ,&ctrdy2);
	  //fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	}
	//lesson learned: fp didn't help much. could have just passed a number

	if(*at(&T2,&xi,&yi)>abs((Tn)*1.5) || *at(&T2,&xi,&yi)<abs(T0*.5)){
	  fprintf(stderr,"ABORTED: ridiculous T calculated  %f",*at(&T2,&xi,&yi));
	  exit(1);}
	{
	  fprintf(fo,"%d %d %f\n",xi,yi,*at(&T2,&xi,&yi));
	} 


     }
    }
 
    nt *newT=T2.u, *oldT=T.u;
    T2.u=oldT; T.u=newT;//make newT old and viceversa
  
  }

  //cleanup
  free(T2.u);
  free(T.u);
  fclose(fo);

}






//the following is unreadable!
//dT/dt+v.delT= k del2T
//      upwind  ctr diff

//void sigma(float *o, float *v, float *dt, float *dx){
// *o=*v*(*dt)/(*dx);//dx 
//}

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
