#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>




typedef float mf;//myfloat
typedef unsigned int mui;


//CELL DEFS

typedef enum { normal, cancer, complx, necrotic } typecellenum;
typedef struct structcell structcell;//why like this??
struct structcell  {
  structcell *next; //ptrs to neigbors
  //llcells llnbrs;
  typecellenum typenum;
  mf A;
} ;
/* typedef struct { */
/*   structcell *next; */
/* } llcells; */


//MATRIX DEFS

typedef enum      { floatt=0,   uintt=1} mattypeenum;
mui typesizes[2]={  sizeof(mf), sizeof(mui)};

typedef struct { 
  void *data;
  mui nrows,ncols; 
  mattypeenum typenum;
} structmatrix;

structmatrix makematrix(mui nrows, mui ncols, mattypeenum tn){
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

/* void initfltmatrix(mf iv, structmatrix *m){//could just do memset */
/*   mui ri,ci; */
/*   mf *v; */
/*   for(ci=0;ci<m->ncols;ci++){ */
/*   for(ri=0;ri<m->nrows;ri++){ */
/*     v=(mf*) idx(&ri,&ci,m); *v=iv; //how do i assign directly? */
/*   }} */
/* } */

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

  structmatrix m=makematrix(npts,3,floatt);
  mui ri,xc=0,yc=1,zc=2;
  for(ri=0;ri<npts;ri++){
    fscanf(fp, "%f %f %f",(mf*)idx(&ri,&xc,&m)
	                 ,(mf*)idx(&ri,&yc,&m)
	                 ,(mf*)idx(&ri,&zc,&m));
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





mf calcA(mui *ti,structmatrix *conn, structmatrix *coords){
  //mui ti//,ci;//ti for triangle index
  mui ac=0,bc=1,cc=2,xc=0,yc=1,zc=2;
  mui a,b,c;
  mf abx,aby,abz,acx,acy,acz;
#define abcxyz(bc,xcyczc) \
                      +(*(mf*)idx(&bc,&xcyczc,coords))	\
                      -(*(mf*)idx(&a ,&xcyczc,coords))

    a=(*(mui*)idx(ti,&ac,conn));
    b=(*(mui*)idx(ti,&bc,conn));
    c=(*(mui*)idx(ti,&cc,conn));
    abx=abcxyz(b,xc);
    aby=abcxyz(b,yc);
    abz=abcxyz(b,zc);
    acx=abcxyz(c,xc);
    acy=abcxyz(c,yc);
    acz=abcxyz(c,zc);
   
    mf z=  +abx*acy
           -aby*acx ; 
    mf y=  +abz*acx
      -abx*acz;
    mf x=  +aby*acz
	   -abz*acy;
    mf A=pow(x*x+y*y+z*z,.5)/2;
    return A;

#undef abcxyz
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


void printmat(structmatrix *mat){
  mui ri,ci;
  for(ri=0;ri<mat->nrows;ri++){
    for(ci=0;ci<mat->ncols;ci++){
      if(mat->typenum==floatt){printf("\n%u %u %f",ri,ci,*(mf*)idx(&ri,&ci,mat));}
      else {printf("\n%u %u %u",ri,ci,*(mui*)idx(&ri,&ci,mat));}
    }
  }
}




void *popcells(structmatrix *conn,structmatrix *coords){
  mui zero=0;

  //array of pointers to structcells
  structcell *pcells=malloc(conn->nrows*sizeof(structcell));
  mui icell,cellii,iabc,abcii,lookfor,at;
  structcell *pcc;
  //llcells cll;

  //for acell: initit
  //for anode in acell:
  //look for it in all other cell
  //append it to neighbors list

  //iterate over cells
  for(icell=0
	;icell<3//icell<conn->nrows
	;icell++){
    //cell initializations
    pcc=&pcells[icell];
    pcc->A=calcA(&icell,conn,coords); 
    pcc->typenum=normal;
    // cll.next=malloc(self)//kkeep adding to list thru mallo(ptr to cell)

    //for a node in this cell
    for(iabc=0;iabc<3;iabc++){
      lookfor=*(mui*)idx(&iabc,&zero,conn);

      //look for it in other cells
      for(cellii=0;cellii<conn->nrows;cellii++){
	if(cellii==icell){continue;}
	for(abcii=0;abcii<3;abcii++){
	  //innermost loop
	  at=*(mui*)idx(&abcii,&zero,conn);
	  if(at==lookfor){
	    //pcll=exit this loop when found
	  }

	}
      }
    

    }
    

  }    
    

}



int main(){

  FILE *gfp=   fopen("face"  ,"r");

  //READING
  structmatrix coords=readcoords(gfp);
  structmatrix conn=    readconn(gfp);

  //ALLOC GLOBAL MATRICES

  mui cri;
  mf A;

  popcells(&conn,&coords);

  /* //connectivity loop */
  /* for(cri=0 */
  /* 	;cri<3//;cri<conn.nrows */
  /* 	;cri++){//0,1,2,3,4,5,6,7,8,9,10.... */

  /*   A=calcA(&cri,&conn,&coords);// 1/2A */
  /*   printf("\n%f",A); */

  /* }//conn row iter */



  //OUTPUT
  //writemat(&coords,"coords");

  //free link ptrs

  return 0;
};


