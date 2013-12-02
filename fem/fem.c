#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
?:
-dot prod rhs?
-0s for missing phis?
-4,4,4
*/


typedef float mf;//myfloat
typedef unsigned int mui;

typedef enum      { floatt=0,   uintt=1} typeenum;
mui typesizes[2]={  sizeof(mf), sizeof(mui)};


//MATRIX DEFS

typedef struct { 
  void *data;
  mui nrows,ncols; 
  typeenum typenum;
} structmatrix;

structmatrix makematrix(mui nrows, mui ncols, typeenum tn){
  structmatrix mat;
  mat.nrows=nrows;
  mat.ncols=ncols;
  mat.typenum=tn;
  mat.data=malloc(nrows*ncols*typesizes[tn]);//don't forget to call free()
  return mat; 
}
void *idx(mui *ri, mui *ci, structmatrix *amat){
  if((*ri>=amat->nrows) || (*ci>=amat->ncols)){ printf("outofbounds %u %u",*ri,*ci);exit(1);}
  //col mjr
#define poffset typesizes[amat->typenum]*((*ri)+(amat->nrows)*(*ci))
  return amat->data+poffset;
#undef poffset
}

void initfltmatrix(mf iv, structmatrix *m){//could just do memset
  mui ri,ci;
  mf *v;
  for(ci=0;ci<m->ncols;ci++){
  for(ri=0;ri<m->nrows;ri++){
    v=(mf*) idx(&ri,&ci,m); *v=iv; //how do i assign directly?
  }}
}

typedef struct {
  structmatrix coords;
  structmatrix conn;
} gridstruct ;

mui MAXLINE=50;
FILE *readto(FILE *fp, char *line){
  //returns pts to line after line
  char readline[MAXLINE];
  while(fscanf(fp,"%s",readline)!=EOF){
    if(strcmp(readline,line)==0){return fp;}
  }
  return fp;
}

structmatrix readcoords(FILE *fp){
  mui npts;
  char ptsstr[]="POINTS";
  fp=readto(fp,ptsstr);
  fscanf(fp,"%u",&npts);

  structmatrix m=makematrix(npts,2,floatt);
  mui ri,xc=0,yc=1;mf dontcare;//not storing 0s (3rd dim)
  for(ri=0;ri<npts;ri++){
    fscanf(fp, "%f %f %f",(mf*)idx(&ri,&xc,&m)
	                 ,(mf*)idx(&ri,&yc,&m),&dontcare);
  }

  return m;
}
structmatrix readconn(FILE *fp){
  mui ntri;
  char meshstr[]="mesh";
  fp=readto(fp,meshstr);
  fscanf(fp,"%u",&ntri);

  structmatrix m=makematrix(ntri,3,uintt);
  mui ri,ac=0,bc=1,cc=2;
  mui a,b,c,*pa,*pb,*pc;
  for(ri=0;ri<ntri;ri++){
    fscanf(fp,"%u %u %u",&a,&b,&c);
    a--;b--;c--;//to make node index starts at 0
    pa=( idx(&ri,&ac,&m)); *pa=a;
    pb=( idx(&ri,&bc,&m)); *pb=b;
    pc=( idx(&ri,&cc,&m)); *pc=c;
  }
  
  return m; 
}


structmatrix readphi(FILE *fp, mui nphi){
  structmatrix m=makematrix(nphi,1,floatt);

  mui ri=0,c=0;
  for (ri=0;ri<nphi;ri++){
    fscanf(fp,"%f",(mf*)idx(&ri,&c,&m));
  }

  return m;
}





typedef struct { structmatrix D; mf A;} dercalcout;
dercalcout calcdernA(mui *ti,structmatrix *conn, structmatrix *coords){
  //should not have written it as a loop
  //gives the derivative of element. div by 2A to get correct deivative
  //mui ti//,ci;//ti for triangle index
  mui ac=0,bc=1,cc=2,xc=0,yc=1;
  mui a,b,c;
  mf abx,aby,acx,acy;
  mf *n1x,*n2x,*n3x,*n1y,*n2y,*n3y;
#define abcxy(bc,xcyc) \
                      +(*(mf*)idx(&bc,&xcyc,coords))	\
                      -(*(mf*)idx(&a ,&xcyc,coords))

  structmatrix D=makematrix(conn->ncols,2,floatt);//acol for x,y
  //for(ti=0;ti<A.nrows;ti++){
    a=(*(mui*)idx(ti,&ac,conn));
    b=(*(mui*)idx(ti,&bc,conn));
    c=(*(mui*)idx(ti,&cc,conn));
    abx=abcxy(b,xc);
    aby=abcxy(b,yc);
    acx=abcxy(c,xc);
    acy=abcxy(c,yc);
    n1x=(mf*)idx(&ac,&xc,&D); *n1x=-acy+aby;
    n2x=(mf*)idx(&bc,&xc,&D); *n2x= acy;
    n3x=(mf*)idx(&cc,&xc,&D); *n3x=-aby;
    n1y=(mf*)idx(&ac,&yc,&D); *n1y= acx-abx;
    n2y=(mf*)idx(&bc,&yc,&D); *n2y=-acx;
    n3y=(mf*)idx(&cc,&yc,&D); *n3y= abx;

    // }

    dercalcout o;
    o.D=D;
    o.A=fabs(abx*acy-acx*aby)/2.0;//getting negative A
    return o;

#undef abcxy
}




structmatrix matrixmul(structmatrix *A,structmatrix *B){
  if(A->ncols!=B->nrows)
    {printf("incompatible matrices for multiplication");exit(1);}

  structmatrix C=makematrix(A->nrows,B->ncols,floatt);
#define  Bv(ri,ci) *(mf*) idx(&ri,&ci, B)
#define  Av(ri,ci) *(mf*) idx(&ri,&ci, A)
#define  Cv(ri,ci) *(mf*) idx(&ri,&ci,&C)

  mui ri,ci,ici;
  for(ri=0;ri<C.nrows;ri++){//for a row,col in product, it's the summation of
    for(ci=0;ci<C.ncols;ci++){//...row A times col B

      mf sum=0;
      for(ici=0;ici<A->ncols;ici++){
	//double a=Av(ri,ici),b=Bv(ci,ici);
	;sum+=(Av(ri,ici)*Bv(ici,ci));}

      Cv(ri,ci)=sum;
    }
  }

  return C;

#undef Av
#undef Bv
#undef Cv
}


structmatrix calcd( structmatrix *Mc
		    ,structmatrix *Ml
		    ,structmatrix *u){

  structmatrix Mccpy=makematrix(Mc->nrows,Mc->ncols,floatt);
  memcpy(Mccpy.data,(Mc->data),sizeof(floatt)*(Mc->nrows*Mc->ncols)); 
  structmatrix d=makematrix(Mc->nrows,1,floatt);

  //Mc - Ml
  mui ri,zero=0;
  mf *mcv;
  for(ri=0;ri<Mc->nrows;ri++){
             mcv=idx(&ri,&ri  ,&Mccpy);
    *mcv=  *(mf*)idx(&ri,&ri  ,Mc)
         -(*(mf*)idx(&ri,&zero,Ml));
  }

  d=matrixmul(&Mccpy,u);
  free(Mccpy.data);

  return d;	   
}
structmatrix calcup1( structmatrix *r
		      ,structmatrix *d
		      ,structmatrix *Mlm1){

  structmatrix up1=makematrix(r->nrows,1,floatt);


  mui ri,zero=0;//(r+d)*Ml^-1
  mf *pup1;
  for(ri=0;ri<r->nrows;ri++){
    pup1=  idx(&ri,&zero,&up1);
    *pup1=(  *(mf*)idx(&ri,&zero,r)
	    +*(mf*)idx(&ri,&zero,d) )
           *(*(mf*)idx(&ri,&zero,Mlm1));
  }

  return up1;
}


void printmat(structmatrix *mat){
  mui ri,ci;
  for(ri=0;ri<mat->nrows;ri++){
    for(ci=0;ci<mat->ncols;ci++){
      if(mat->typenum==floatt){printf("\n%u %u %f",ri,ci,*(mf*)idx(&ri,&ci,mat));}
      else {printf("\n%u %u %u",ri,ci,*(mui*)idx(&ri,&ci,mat));}
    }
  }
}
void writemat(structmatrix *mat, char *fname){
  FILE *f=fopen(fname,"w");
  mui ri,ci;
  for(ri=0;ri<mat->nrows;ri++){
    for(ci=0;ci<mat->ncols;ci++){
      if(mat->typenum==floatt){fprintf(f,"\n%u %u %f",ri,ci,*(mf*)idx(&ri,&ci,mat));}
      else {fprintf(f,"\n%u %u %u",ri,ci,*(mui*)idx(&ri,&ci,mat));}
    }
  }
}



int main(){
	
  printf("asdfgrgrgrhrhrdhd\n");

  FILE *gfp=   fopen("grid0"  ,"r");
  FILE *phifp=fopen( "phi0.0","r");

  //READING
  structmatrix coords=readcoords(gfp);
  structmatrix conn=    readconn(gfp);
  structmatrix phi=     readphi(phifp,coords.nrows);

  //ALLOC GLOBAL MATRICES
  structmatrix Mc=makematrix(coords.nrows,coords.nrows,floatt);
  structmatrix Ml=makematrix(coords.nrows,1           ,floatt);
  structmatrix rx=makematrix(coords.nrows,1           ,floatt);
  structmatrix ry=makematrix(coords.nrows,1           ,floatt);
  initfltmatrix(0,&Mc);//set 0s
  initfltmatrix(0,&Ml);
  initfltmatrix(0,&rx);
  initfltmatrix(0,&ry);

  //LOCAL MATRICES
  mf aNi[3][1]={ {1}
		,{1}
		,{1}};//  *A/3
  structmatrix INi   =makematrix(3,1,floatt);
  mf aNiNj[3][3]={ {2,1,1}
		  ,{1,2,1}
		  ,{1,1,2}}; // *A/12
  structmatrix INiNj =makematrix(3,3,floatt);
  /* mf aNiNjl[3][3]={ {4,0,0}// not used */
  /* 		   ,{0,4,0} */
  /* 		    ,{0,0,4}};//   *A/12 */
  /* structmatrix INiNjl=makematrix(3,3,floatt); */
  mui ri,ci;mf *v;
  for(ri=0;ri<3;ri++){
  for(ci=0;ci<1;ci++){ v=(mf*)idx(&ri,&ci,&INi)   ; *v=  aNi[ri][ci];}}
  for(ri=0;ri<3;ri++){
  for(ci=0;ci<3;ci++){ v=(mf*)idx(&ri,&ci,&INiNj) ; *v=aNiNj[ri][ci];}}
  /* for(ri=0;ri<3;ri++){ */
  /* for(ci=0;ci<3;ci++){ v=(mf*)idx(&ri,&ci,&INiNjl) ;*v=aNiNjl[ri][ci];}} */


  //FILL GLOBAL MATRICES
  mui cri,xc=0,yc=1,zero=0;
  mui sri,sci;
  mui abcr,abcc;
  mf *src,*Mcdst,*Mldst,*rxdst,*rydst;
  mf *pphi;
  mf *pDx,*pDy,*pA;
  mf delNjphijx;//dot prod delNj*phij
  mf delNjphijy;
  mui idp;
  dercalcout AnD;
  //connectivity loop
  for(cri=0;cri<conn.nrows;cri++){//0,1,2,3,4,5,6,7,8,9,10....

  AnD=calcdernA(&cri,&conn,&coords);// 1/2A
  pA=&(AnD.A);

  
    //lhs stuff: loop over the M submat
    for(sri=0;sri<INiNj.nrows;sri++){//0,1,2
    for(sci=0;sci<INiNj.ncols;sci++){//0,1,2
      src=          idx(&sri, &sci ,&INiNj);//a value inside the submat
      abcc=*(mui*)  idx(&cri, &sci ,&conn);//the col index that gets a,b,or c
      abcr=*(mui*)  idx(&cri, &sri ,&conn);
      Mcdst=        idx(&abcr,&abcc,&Mc);// ??? is this ok?
      *Mcdst+=(*src)*(*pA)/12.;
	}
    abcr=*(mui*)  idx(&cri  ,&sri  ,&conn);//the col index that get a,b, or c
    Mldst=        idx(&abcr ,&zero ,&Ml);
    *Mldst+=4.0*(*pA/12);
    //rhs stuff
    delNjphijx=0;
    delNjphijy=0;
    for(idp=0;idp<INi.nrows;idp++){//0,1,2 //could be integrated in previous loop
      //but i wanted to isolate lhs and rhs avoid confusion
      abcr=*(mui*)  idx(&cri  ,&idp  ,&conn);//node num
      pphi=         idx(&abcr ,&zero ,&phi);//ptr to phi
      pDx=          idx(&idp  ,&xc   ,&AnD.D);//ptr to derivative
      pDy=          idx(&idp  ,&yc   ,&AnD.D);
      delNjphijx+=(*pphi*(*pDx));
      delNjphijy+=(*pphi*(*pDy));
    }
    
    src=          idx(&sri  ,&zero ,&INi); //value inside the submat
    rxdst=        idx(&abcr ,&zero ,&rx);  //ptr to destination
    rydst=        idx(&abcr ,&zero ,&ry);  //
    *rxdst+=(*src)*(*pA/3.) * (delNjphijx/(2*(*pA)));//could elim A here
    *rydst+=(*src)*(*pA/3.) * (delNjphijy/(2*(*pA)));
  
}

  free(AnD.D.data);

  }//conn row iter


  //SOLVING PROCEDURE
  //make inverse matrix
  structmatrix Mlm1=makematrix(Ml.nrows,Ml.ncols,floatt);
  mf *inv;
  for(ri=0;ri<Mlm1.nrows;ri++){
    inv=          idx(&ri,&zero,&Mlm1);
   *inv=pow(*(mf*)idx(&ri,&zero,&Ml),-1);
  }

  structmatrix vx =makematrix(coords.nrows,1,floatt);
  structmatrix vy =makematrix(coords.nrows,1,floatt); 
  structmatrix up1=makematrix(coords.nrows,1,floatt); 
  structmatrix d;
  initfltmatrix(0,&vx);
  initfltmatrix(0,&vy);

  mui i,equals,ii,SOLNNITER=100;
  mf EQUALTOL=.01;

 
  //iterloop
  for(i=0;i<SOLNNITER;i++){
    d=calcd(&Mc,&Ml,&vx);
    up1=calcup1(&ry,&d,&Mlm1);
    free(d.data);
    equals=0;
    for(ii=0;ii<up1.nrows;ii++){
      if( fabs(   (    *(mf*)idx(&ii,&zero,&up1)
  		      -*(mf*)idx(&ii,&zero,&vx)
  		      )				\
  		  /   (*(mf*)idx(&ii,&zero,&vx))
  		  ) < EQUALTOL)
  	{equals++;}
    }
    free(vx.data);
    vx=up1;
    if(equals==up1.nrows){printf("\nx soln converged");break;}
  }
  if(i==SOLNNITER){printf("\nx soln did not converge");}

  //too lazy to make a function so i'm repeating code. bad majid.
  //iterloop
  for(i=0;i<SOLNNITER;i++){
    d=calcd(&Mc,&Ml,&vy);
    up1=calcup1(&ry,&d,&Mlm1);
    free(d.data);
    equals=0;
    for(ii=0;ii<up1.nrows;ii++){
      if( fabs(   (    *(mf*)idx(&ii,&zero,&up1)
  		      -*(mf*)idx(&ii,&zero,&vy)
  		      )				\
  		  /   (*(mf*)idx(&ii,&zero,&vy))
  		  ) < EQUALTOL)
  	{equals++;}
    }
    free(vy.data);
    vy=up1;
    if(equals==up1.nrows){printf("\ny soln converged");break;}
  }
  if(i==SOLNNITER){printf("\ny soln did not converge");}



  //OUTPUT
  writemat(&coords,"coords");
  writemat(&vx,"vx");
  writemat(&vy,"vy");
  writemat(&phi,"phi");

   //mui tr=911,tc=912;printf("%f",*(float*) idx(&tr,&tc,&Mc));

  free(Mlm1.data);
  free(Ml.data);
  free(Mc.data);
  free(rx.data);
  free(ry.data);
  free(coords.data);
  free(conn.data);
  free(phi.data);
  free(vx.data);
  free(vy.data);

  return 0;
};


