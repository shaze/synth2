

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>


#define max_cluster 3



static int  num_s, num_v, num_new, num_clusters, bsize, new_bsize;
static unsigned char * old_bed;

static int mask[4] = {3,12,48,192};

static int chrom_start[27];


int getBim(char * bim, int *chrom_start) {
  // get the start position of each chromosome in the bim file
  // if there aren't any in that region, set to -1
  FILE *f;
  int chrom, cm, pos, num_read; // data in bim file
  int next_chrom, line_num;
  char rs[8192], maj[32], min[32];

  
  f = fopen(bim,"r");
  num_read = fscanf(f,"%d%s%d%d%s%s\n",&chrom, rs, &cm, &pos, maj, min);
  next_chrom = line_num = 0;
  while (num_read == 6) {
      while(next_chrom<chrom) { // missing chroms in bim file
	  chrom_start[next_chrom]=-1;
	  next_chrom++;
      }
      if (chrom==next_chrom)  {
	 chrom_start[next_chrom] = line_num;
	 next_chrom++;
     }
      line_num++;
      num_read = fscanf(f,"%d%s%d%d%s%s\n",&chrom, rs, &cm, &pos, maj, min);
  }
  while (next_chrom<27) {
      chrom_start[next_chrom]=-1;
      next_chrom ++;
  }
  /*
  for(next_chrom = 0; next_chrom<27; next_chrom++) {
    printf("chrom[%d] = %d\n",next_chrom, chrom_start[next_chrom]);
    }*/
  fclose(f);
  return line_num;
}



unsigned char * readBed(char *bed) {
    FILE *bf;
    int num_read;
    unsigned char magic[3], *bdata;

    bf = fopen(bed, "r");
    num_read =  fread(magic, 1, 3, bf);
    assert ( (magic[0]==108) && (magic[1]==27) && (magic[2]==1));
    bdata = calloc(num_v,bsize);
    num_read = fread(bdata, bsize, num_v, bf);
    return bdata;
}

unsigned char * allocBed(int num_new) {
    unsigned char  *bdata;
    bdata = malloc(3+num_v*new_bsize);
    bzero(bdata,3+num_v*new_bsize);
    bdata[0]=108;
    bdata[1]=27;
    bdata[2]=1;
    return bdata;
}

int getFam(char * fam) {
  FILE *f;
  char dummy[4096];
  int num=0, res;
  f = fopen(fam,"r");
  while (fscanf(f,"%[^\n]\n",dummy)  != -1) num++;
  fclose(f);
  return num;
}  

int * getClusters(char * cluster_fname, int *num) {
  FILE *f;
  int *cluster, d, nrd,j,m;
  char dummy[4096],c;
  *num=0;
  f = fopen(cluster_fname,"r");
  while (fscanf(f,"%[^\n]\n",dummy)  != -1) (*num)++;
  fclose(f);
  cluster = (int *) malloc(sizeof(int)*(*num)*max_cluster);
  f = fopen(cluster_fname,"r");  
  for (int i=0; i<*num; i++) {
     j=0;
      do {
	 m=fscanf(f,"%d%[\n]",&d,&c);
	 cluster[i*max_cluster+j]=d;
	 j++;
      } while (m<2);
      while(j<max_cluster) {
	 cluster[i*max_cluster+j]=-1;
	 j++;
      }
  }
  return cluster;
}

void updateChoices(int * choice, int v, int * clusters, int *rnds, int gen) {
   // Update which individual will chosen as source for the new person
   // based upon genome position. At the moment, this is done at
   // the start of a new chromosome -- to keep consistency we don't change
   // after chromosome 23
   static int chrom = 0 ;  // NB: static
   int next_choice;
   if ((v<=chrom_start[chrom]) || (chrom>23)) return; // no update
   chrom++;
   for (int i=0; i<gen; i++) {
      next_choice = choice[i]+1;
      if ((next_choice==max_cluster) || (clusters[rnds[i]*max_cluster+next_choice]==-1))
	 choice[i]=0;
      else {
	 choice[i]=next_choice;
      }
   }
}
	     
  

void generateNewVariant(unsigned char * new_bed,int *clusters, int v, int *rnds, int gen,
			             int * choice) {
   int    src_idx, src_pr, trg_idx, trg_pr, src;
   int   shift;
   unsigned char *old, *trg, data;

   updateChoices(choice, v, clusters, rnds, gen);
   // index into the BED file for this variant's block
   old = 3+ old_bed + v*bsize; // (3 is for magic numbers)
   trg = 3+ new_bed+ v*new_bsize;
   for(int i=0; i<num_new; i++) {
        // we're producing a smaller set, so the new "person" is at a different
      // position, both byte and then postion within the byte
      trg_idx = i>>2;            // byte for the new person 
      trg_pr  = i & 0x03;       // where in the byte
      src       = clusters[rnds[i]*max_cluster+choice[i]];
      src_idx = src>>2;        // byte of old person
      src_pr  = src & 0x03;    //  where in the byte 
      // get data we need
      data     = old[src_idx] & mask[src_pr];
      // put in right place -- probably need to shift it to left or right
      // need to shift by multiples of two
      shift = (src_pr-trg_pr)<<1; // *2 
      if (shift>0)
	 trg[trg_idx] |= (data << shift);
      else
	 trg[trg_idx] |= (data >> -shift);
   }
}


int * shuffled(int gen) {
   // generate a shuffled data set
   int * rnds, i, swap, next;
   rnds = malloc(sizeof(int)*gen);
   for (i=0; i<gen; i++) rnds[i]=i;
   for (i=0; i<gen-1; i++) {
      next = arc4random_uniform(gen);
      swap = rnds[gen-i-1];
      rnds[gen-i-1]= rnds[next];
      rnds[next] = swap;
   }
   return rnds;
}

int * generateSamples(unsigned char * new_bed, int * clusters) {
   int gen, *rnds, i, next,swap, v;
   int choice[num_new]; // for each person who the current source is 

   bzero(choice, num_new*sizeof(int));

   gen = num_new < num_clusters ? num_new : num_clusters;
   rnds  = shuffled(num_s);
   for (v=0; v<num_v; v++) 
      generateNewVariant(new_bed, clusters, v, rnds, gen, choice);
   return rnds;
}

void outputData(char * base, unsigned char * new_bed, int * is_case) {
  FILE *g;
  char outfname[4096];
  int   phe;

  strcpy(outfname,base);
  strcat(outfname,".bed");
  g = fopen(outfname,"w");
  fwrite(new_bed,1,3+new_bsize*num_v,g);
  fclose(g);
  strcpy(outfname,base);
  strcat(outfname,".fam");
  g = fopen(outfname,"w");
  for(int i=0; i<num_new; i++) {
     phe = 1+is_case[i];
     fprintf(g,"D%04d\tD%04d\t0\t0\t0\t%d\n",i,i,phe);
  }
  fclose(g);
}


unsigned char get_prob(float prob_cut) {
   float coin;
   unsigned int   num=0, res;
   coin =  (float) arc4random_uniform(100000)/100000;
   if (coin < prob_cut) num++;
   coin =  (float) arc4random_uniform(100000)/100000;
   if (coin < prob_cut) num++;
    switch (num) {
      case 0: res=3;
	 break;
      case 1: res=2;
	 break;
      case 2: res=0;
   }
   return res;
}

int * mutateData(char * mute_file, unsigned char * new_bed, int num_new) {
   /* mute_file contains a list of mutations -- each line a triple
    * index in bim file, prob of mutation in case, prob of mutation in control */
   FILE *f;
   int * rnd, v, *is_case;
   int    p_index; // which byte in the block
   int    p_pr=0;     // which two bit of  the byte
   float  case_p, cntl_p, prob, prob_cut;
   unsigned char  data;

   is_case = malloc(num_new*sizeof(int));
   bzero(is_case, num_new*sizeof(int));
   rnd = shuffled(num_new);
   // make the first half cases
   for(int i=0;i<num_new/2;i++) {
      is_case[rnd[i]]=1;
   }
   f = fopen(mute_file, "r");
   while (fscanf(f,"%d %f %f\n", &v, &case_p, &cntl_p) != -1) {
      data =  0;
      p_pr = 0;
      for(int p=0; p<num_new; p++) {
	 prob_cut = is_case[p] ? case_p : cntl_p;
	 p_index = p>>2;
	 // Add to data -- the lowest two-bits in the byte contain the
	 // lowest numbered individual
	 data =  data  | (get_prob(prob_cut) << (2*p_pr));
	 p_pr++;
         if (p_pr == 4) {
	    int bptr = 3+v*new_bsize+p_index;
  	    new_bed[bptr] = data;
	    p_pr=0;
	    data = 0;
	 }
      }
   }
   return is_case;
}


int main(int argc, char *argv[]) {

   int  *clusters, *is_case;
  unsigned char * new_bed;
  char  bed[4096], bim[4096], fam[4096];
  strcpy(bed, argv[1]);
  strcpy(bim, argv[1]);
  strcpy(fam, argv[1]);  
  strcat(bed, ".bed");
  strcat(bim, ".bim");
  strcat(fam, ".fam");  

  num_new         = atoi(argv[2]);
  new_bsize        = (num_new+3)/4;
  num_s             = getFam(fam);
  bsize               = (num_s+3)/4;
  num_v            = getBim(bim, chrom_start);
  old_bed           = readBed(bed);
  new_bed         = allocBed(num_new);
  clusters           = getClusters(argv[3],&num_clusters);

  generateSamples(new_bed, clusters);

  is_case = mutateData(argv[4],new_bed, num_new);
  
  outputData(argv[5], new_bed, is_case);
 
}
