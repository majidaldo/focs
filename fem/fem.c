#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

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
  mui ri,a=0,b=1,c=2;
  for(ri=0;ri<ntri;ri++){
    fscanf(fp,"%u %u %u",idx(&ri,&a,&m)
  	                ,idx(&ri,&b,&m)
                        ,idx(&ri,&b,&m));
  }

  return m;
}


int main(){
	
  printf("asdfgrgrgrhrhrdhd\n");

  FILE *gfp=fopen("grid0","r");

  matrixstruct coords=readcoords(gfp);
  matrixstruct conn=    readconn(gfp);

  //mui tr=2,tc=1;printf("%f",*(float*) idx(&tr,&tc,&coords));

  //char tc[]="POINTS";  
  /* gf=readto(gf,tc); */
  /* fscanf(gf,"%s",tc); */
  /* printf(tc); */

};

