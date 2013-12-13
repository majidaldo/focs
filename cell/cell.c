#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>



typedef float mf;//myfloat
typedef unsigned int mui;


//CELL DEFS

typedef enum { normal=0, cancer=1, complx=2, necrotic=3 } typecellenum;
//numbered to make sure it stays the same for post processing
typedef struct structcell structcell;//why like this??
typedef struct structlink structlink;
struct structcell  {
  mui nnbrs;
  structlink *nbrs;
  typecellenum typenum;
  mf A;
} ;
struct structlink {///begins with self as first member
  structlink *next;
  structcell *curr;
} ;


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



structcell *popcells(structmatrix *conn,structmatrix *coords){

  //array of pointers to structcells
  structcell *pcells=malloc(conn->nrows*sizeof(structcell));
  mui icell,cellii,iabc,abcii			\
    ,lookfor,at;
  structcell *pcc;
  structlink *plink,*pnewlink;
  //for acell: initit
  //for anode in acell:
  //look for it in all other cell
  //append it to neighbors list

  //iterate over cells
  for(icell=0
	;icell<conn->nrows
	;icell++){
    //cell initializations
    pcc=&pcells[icell];
    pcc->A=calcA(&icell,conn,coords); 
    pcc->typenum=normal;
    plink=malloc(sizeof(structlink));
    pcc->nbrs=plink;
    plink->curr=pcc;//ptr to self as a start
    pcc->nnbrs=0;

    
    for(cellii=0;cellii<conn->nrows;cellii++){//look at other nodes
      if(cellii==icell){continue;}//except self

      for(iabc=0;iabc<3;iabc++){
	lookfor=*(mui*)idx(&icell,&iabc,conn);//look for this node */
    	for(abcii=0;abcii<3;abcii++){ 
 	  at=*(mui*)idx(&cellii,&abcii,conn); 
   	  if(at==lookfor){
    	    pnewlink=malloc(sizeof(structlink));
    	    pnewlink->curr=&pcells[cellii];//aah i love C!
    	    plink->next=pnewlink;
    	    plink=pnewlink;
    	    (pcc->nnbrs)++;
    	    //membership est. exit loop 
	    goto memest;//i had too! needed to break out of two loops
	  }
	}
      }
    memest: ;


    }//look at other nodes


  }

    
  return pcells;
}


bool prob_true(double p){ 
    return rand() < p * (RAND_MAX+1.0);
}


void normal2cancer(structlink *pll, mf *p2cancer){
  mf allarea,nbrcancerarea=0;
  structcell *pcell=pll->curr;//self
  structlink *pcl;
  structcell *pcc;
 
  pcl=pcell->nbrs;//
  pcc=pcl->curr;//self
  allarea=pcell->A;
  mui cellni; 
  for(cellni=0;cellni<pcell->nnbrs;cellni++){//neighbor loop
    pcl=pcl->next;
    pcc=pcl->curr;
	
    if(pcc->typenum==cancer){nbrcancerarea+=pcc->A;}
    allarea+=pcc->A;

  }
  
  bool tf;
  if(nbrcancerarea==0){tf=false;}
  else{tf= prob_true(*p2cancer*(nbrcancerarea/allarea));}
  //return tf;
  //tf=true;//code check
  if(tf==true){pcell->typenum=cancer;}
}


void cancer2complex(structlink *pll, mf *p2complex){
  mui nbrnotnormal=0;
  structcell *pcell=pll->curr;//self
  structlink *pcl;
  structcell *pcc;

  pcl=pcell->nbrs;//
  pcc=pcl->curr;//self
  mui cellni;
  for(cellni=0;cellni<pcell->nnbrs;cellni++){//neighbor loop
    pcl=pcl->next;
    pcc=pcl->curr;
	
    if(pcc->typenum!=normal){nbrnotnormal++;}

  }

  bool tf;
  if(nbrnotnormal==pcell->nnbrs){tf=prob_true(*p2complex);}
  else{tf=false;}
 //return tf;
  //tf=true;//code check
  if(tf==true){pcell->typenum=complx;}
}


void complex2necrotic(structlink *pll, mf *p2necrotic){
  structcell *pcell=pll->curr;//self

  bool tf;
  if(complx==pcell->typenum){tf=prob_true(*p2necrotic);}
  else{tf=false;}
 //return tf;
  //tf=true;//code check
  if(tf==true){pcell->typenum=necrotic;}

}



void writecells(FILE *f, structcell *pcells, mui ncells){
  mui wi;
  for(wi=0;wi<ncells;wi++){
    fprintf(f,"%d ",pcells[wi].typenum);
  }

}


int main(int argc, char *argv[]){
  if(argv[1]==NULL || argv[2]==NULL || argv[3]==NULL)
    {printf("provide 3 args for probilities (factors): \
normal->cancer, cancer->complex, complex->necrotic. ");
      exit(1);}
  mf p2cancer;   sscanf(argv[1],"%f",&p2cancer);
  mf p2complex;  sscanf(argv[2],"%f",&p2complex);
  mf p2necrotic; sscanf(argv[3],"%f",&p2necrotic);

  mui randn=13;
  if(argv[4]==NULL){
    printf("4th optional arg for random number seed.\
defaults to %d",randn);}
  else{sscanf(argv[4],"%d",&randn);}

  FILE *gfp=   fopen("face"  ,"r");

  //READING
  structmatrix coords=readcoords(gfp);
  structmatrix conn=    readconn(gfp);
  //POPULATE WITH INITIALS
  structcell *pcells=popcells(&conn,&coords);
  //cancerize one
  pcells[0].typenum=cancer;
  
  srand(randn);

  //ITERATE
  mui icell,iter;
  mui necrotics=0;
  structcell *pcc;

  const mui MAXITER=10000;
  FILE *fo=fopen("cells","w");//char fname[5];
  printf("\nnum of necrotic cells\n");
  //for an iteration
  for(iter=0;iter<MAXITER;iter++){

    //sprintf(fname, "%d.cell", iter);
    writecells(fo,pcells,conn.nrows);fprintf(fo,"\n");
    if(necrotics==conn.nrows){break;}


    //---------
    for(icell=0;icell<conn.nrows;icell++){

      pcc=&pcells[icell];
      //rules
      if     (pcc->typenum==normal){   normal2cancer(pcc->nbrs,&p2cancer);}
      else if(pcc->typenum==cancer){  cancer2complex(pcc->nbrs,&p2complex);}
      else if(pcc->typenum==complx){complex2necrotic(pcc->nbrs,&p2necrotic);
	if(pcc->typenum==necrotic){necrotics++;}//add to necrotics if it changed
      }
      //else it's necrotic no need to do anything
    
    }//cell iterator brace
    //---------

    printf("%d ",necrotics);

  }
  if(iter<MAXITER){
    printf("\nfinished in %d iterations",iter);
    FILE *io=fopen("results","a+");
    fprintf(io,"%f %f %f %d %d\n"
	    ,p2cancer,p2complex,p2necrotic,randn,iter);
  }
  else{printf("\nmax iteration reached =%d",MAXITER);}
 
  //other OUTPUT
  writemat(&coords,"coords");
  writemat(&conn  ,"conn");

  //FREE mallocs
  //free links
  structlink *ol,*cl;
  mui ni;
  for(icell=0;icell<conn.nrows;icell++){
    ol=pcells[icell].nbrs;//lnk to self
    for(ni=0;ni<(pcells[icell].nnbrs)+1;ni++){//+1 bc strting w slf
      //advance
      cl=ol->next;
      //burn the bridge
      free(ol);
      ol=cl;
    }//finally burn yourself
    //free(ol);//this segfaults.not sure if i wiped all the list
  }
  free(coords.data);
  free(conn.data);
  free(pcells);

  return 0; 
}

