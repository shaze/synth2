

/* (c) University of the Witwatersrand, Johannesburg on behalf of the Pan-African 
   Bioinformatics Network for H3Africa   */

#include "hash.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int  table[hash_size];

const static int  mask = hash_size-1;


int rs2hash ( char * name) {
  int result=0;
  for(int i=0; name[i]!=0; i++) {
    result = (result<<1);
    result = (result+name[i]) ^ (name[i]<<8)+1027*(name[i]);
  }
  result=result&mask;
  return result;
}


void hash_add(char * snp, int p) {
  int ind, i=1;
  ind = rs2hash(snp);
  while(table[ind] != -1) {
    ind=(ind+i)&mask;
    i=i+2;
  }
  if (strcmp("snp-known61733550",snp)==0)
    printf("Hash index %d\n",ind);
  table[ind]=p;
}


int  hash_find(char ** bim_table, char * snp) {
  int ind, orig,i=1;
  orig=ind = rs2hash(snp);
  while((table[ind]!=-1) && (strcmp(bim_table[table[ind]],snp) != 0)) {
    ind=(ind+i)&mask;
    i=i+2;
  }
  if (table[ind]==-1) {
    printf("No entry <%s> %d in hash table\n", snp, orig);
    exit(-17);
  }
  return table[ind];
}
  

void hash_init () {
  memset(table,-1,hash_size);
}
