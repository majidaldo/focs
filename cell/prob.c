#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

bool prob_true(double p){
    return rand() < p * (RAND_MAX+1.0);
}
unsigned int lcgrand(){
//UINT_MAX
#define a 123
#define b 321
return ; 
#undef a
#undef b
}

int main(){

srand(13);

unsigned int i, st=0;
for (i=0;i<1000;i++){
if(true==prob_true(.2)){printf("|");st++;}
else {printf(".");}
}


printf("\n %d", st);

}