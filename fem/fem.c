#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef float mf;//myfloat
typedef unsigned int mui;

typedef enum      { floatt=0,   uintt=1} typeenum;
mui typesizes[2]={  sizeof(mf), sizeof(mui)};


//MATRIX DEFS

typedef struct { 
  void *data;
  mui nrows,ncols; 
  typeenum typenum;
} matrixstruct;

matrixstruct makematrix(mui nrows, mui ncols, typeenum tn){
  matrixstruct mat;
  mat.nrows=nrows;
  mat.ncols=ncols;
  mat.typenum=tn;
  mat.data=malloc(nrows*ncols*typesizes[tn]);//don't forget to call free()
  return mat; 
}
void *idx(mui *ri, mui *ci, matrixstruct *amat){
  if((*ri>=amat->nrows) || (*ci>=amat->ncols)){ printf("outofbounds %u %u",*ri,*ci);exit(1);}
  //col mjr
#define poffset typesizes[amat->typenum]*((*ri)+(amat->nrows)*(*ci))
  return amat->data+poffset;
#undef poffset
}

void initfltmatrix(mf iv, matrixstruct *m){
  mui ri,ci;
  mf *v;
  for(ci=0;ci<m->ncols;ci++){
  for(ri=0;ri<m->nrows;ri++){
    v=(mf*) idx(&ri,&ci,m); *v=iv; //how do i assign directly?
  }}
}

typedef struct {
  matrixstruct coords;
  matrixstruct conn;
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

matrixstruct readcoords(FILE *fp){
  mui npts;
  char ptsstr[]="POINTS";
  fp=readto(fp,ptsstr);
  fscanf(fp,"%u",&npts);

  matrixstruct m=makematrix(npts,2,floatt);
  mui ri,xc=0,yc=1;mf dontcare;//not storing 0 (3rd num)
  for(ri=0;ri<npts;ri++){
    fscanf(fp, "%f %f %f",idx(&ri,&xc,&m)
                         ,idx(&ri,&yc,&m),&dontcare);
  }

  return m;
}
matrixstruct readconn(FILE *fp){
  mui ntri;
  char meshstr[]="mesh";
  fp=readto(fp,meshstr);
  fscanf(fp,"%u",&ntri);

  matrixstruct m=makematrix(ntri,3,uintt);
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


matrixstruct readphi(FILE *fp, mui nphi){
  matrixstruct m=makematrix(nphi,1,floatt);

  mui ri=0,c=0;
  for (ri=0;ri<nphi;ri++){
    fscanf(fp,"%f",idx(&ri,&c,&m));
  }

  return m;
}


void putsubmat(matrixstruct *mat, matrixstruct *submat
	       , matrixstruct *conn ,matrixstruct *A,mf mtiple){
  mui cri,zero=0;
  mui sri,sci;
  mui abcr,abcc;
  mf *src,*dst;
  for(cri=0;cri<conn->nrows;cri++){//0,1,2,3,4,5,6,7,8,9,10....

    //loop over submat
    for(sri=0;sri<submat->nrows;sri++){//0,1,2
    for(sci=0;sci<submat->ncols;sci++){//0,1,2 or 0
      src=         idx(&sri, &sci ,submat);//a value inside the submat
      abcc=*(mui*) idx(&cri, &sci ,conn);//the col index that gets a,b,or c
      abcr=*(mui*) idx(&cri, &sri ,conn);
      dst=         idx(&abcr,&abcc,mat);// ??? is this ok?
      *dst=(*src)*(*(mf*)idx(&cri,&zero,A))*mtiple;
      //printf("\n %f %u %u",*dst,abcr,abcc);
	}}

  }

}


matrixstruct calcareas(matrixstruct *conn, matrixstruct *coords){
  mui ti,ci;
  mui ac=0,bc=1,cc=2,xc=0,yc=1,zero=0;
  mui a,b,c;
  mf abx,aby,acx,acy;
  mf *area;
#define abcxy(bc,xcyc) \
                      (*(mf*)idx(&bc,&xcyc,coords))	\
                     -(*(mf*)idx(&a ,&xcyc,coords))

  matrixstruct A=makematrix(conn->nrows,1,floatt);
  for(ti=0;ti<A.nrows;ti++){
    a=(*(mui*)idx(&ti,&ac,conn));
    b=(*(mui*)idx(&ti,&bc,conn));
    c=(*(mui*)idx(&ti,&cc,conn));
    abx=abcxy(b,xc);
    aby=abcxy(b,yc);
    acx=abcxy(c,xc);
    acy=abcxy(c,yc);
    area=(mf*)idx(&ti,&zero,&A);
    *area=fabs(abx*acy-aby*acx);
    /* abx=(*(mf*)idx(&b,&xc,coords)) */
    /*    -(*(mf*)idx(&a,&xc,coords)) */
    /* aby=(*(mf*)idx(&b,&yc,coords)) */
    /*    -(*(mf*)idx(&a,&yc,coords)) */
 
    }

#undef abcxy
  return A;
}




void printmat(matrixstruct *mat){
  mui ri,ci;
  for(ri=0;ri<mat->nrows;ri++){
    for(ci=0;ci<mat->ncols;ci++){
      if(mat->typenum==floatt){printf("\n%u %u %f",ri,ci,*(mf*)idx(&ri,&ci,mat));}
      else {printf("\n%u %u %u",ri,ci,*(mui*)idx(&ri,&ci,mat));}
    }
  }
}
//todo writemat



int main(){
	
  printf("asdfgrgrgrhrhrdhd\n");

  FILE *gfp=   fopen("grid0"  ,"r");
  FILE *phixfp=fopen( "phi0.0","r");
  FILE *phiyfp=fopen( "phi0.1","r");

  //reading
  matrixstruct coords=readcoords(gfp);
  matrixstruct conn=    readconn(gfp);
  matrixstruct phix=     readphi(phixfp,coords.nrows);
  matrixstruct phiy=     readphi(phiyfp,coords.nrows);

  //global matrix
  matrixstruct M=makematrix(coords.nrows,coords.nrows,floatt);
  initfltmatrix(0,&M);//set 0s
  
  //areas
  matrixstruct A=calcareas(&conn,&coords);


  //fill
  mf aNi[3][1]={ {1}
		,{1}
		,{1}};//  *A/3
  matrixstruct INi=makematrix(3,1,floatt);
  mf aNiNj[3][3]={ {2,1,1}
		  ,{1,2,1}
		  ,{1,1,2}}; // *A/12
  matrixstruct INiNj=makematrix(3,3,floatt);
  mui ri,ci;mf *v;
  for(ri=0;ri<3;ri++){
  for(ci=0;ci<1;ci++){ v=(mf*)idx(&ri,&ci,&INi)   ; *v=  aNi[ri][ci];}}
  for(ri=0;ri<3;ri++){
  for(ci=0;ci<3;ci++){ v=(mf*)idx(&ri,&ci,&INiNj) ; *v=aNiNj[ri][ci];}}

  putsubmat(&M,&INiNj,&conn,&A,(1.0/12));

  

  //printmat(&M);

  //freealloc
  
  mui tr=911,tc=912;printf("%f",*(float*) idx(&tr,&tc,&M));
};


