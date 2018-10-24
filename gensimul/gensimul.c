

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
typedef __u_char u_char;

#include <sys/types.h>
#include <bsd/stdlib.h>

#define max_cluster 3



typedef struct {
  long  num_v; // num variants in file
  long  num_s; // num sample;
  unsigned char magic[3];
  unsigned char *bed;
} bed_file;
  

static bed_file old_bed, new_bed;
static int      num_clusters;

static int mask[4] = {3,12,48,192};

static int chrom_start[27];
static char ** orig_fam;

int getBlockSize(bed_file *b) {
  return (b->num_s+3)>>2;
}

int getBim(char * bim, int *chrom_start) {
  // get the start position of each chromosome in the bim file
  // if there aren't any in that region, set to -1
  FILE *f;
  float cm;
  int chrom, pos, num_read; // data in bim file
  int next_chrom, line_num;
  char rs[8192], maj[32], min[32];

  
  f = fopen(bim,"r");
  if (f==NULL) {
    printf("Can't open file %s\n",bim);
    exit(10);
  }
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
      num_read = fscanf(f,"%d%s%f%d%s%s\n",&chrom, rs, &cm, &pos, maj, min);
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



unsigned char * readBed(bed_file * bed, char *bname) {
    FILE *bf;
    int num_read;
    long  bsize;

    bsize = getBlockSize(bed);
    bf = fopen(bname, "r");
    if (bf==NULL) {
      printf("Can't open file %s\n",bed);
      exit(11);
    }
    num_read =  fread(bed->magic, 1, 3, bf);
    assert ( (bed->magic[0]==108) && (bed->magic[1]==27) && (bed->magic[2]==1));
    bed->bed = malloc(bed->num_v*bsize);
    num_read = fread(bed->bed, bsize, bed->num_v, bf);
}

unsigned char * allocBed(bed_file * bed) {
    bed->bed = malloc(bed->num_v*getBlockSize(bed));
    bzero(bed->bed,bed->num_v*getBlockSize(bed));
    bed->magic[0]=108;
    bed->magic[1]=27;
    bed->magic[2]=1;
}

int getFam(char * fam) {
  FILE *f;
  char dummy[4096], fid[64], iid[64];
  int num=0, res;
  f = fopen(fam,"r");
  if (f==NULL) {
      printf("Can't open file %s\n",fam);
      exit(12);
  }
  while (fscanf(f,"%[^\n]\n",dummy)  != -1) num++;
  fclose(f);
  orig_fam = malloc(sizeof(char *)*num);
  f = fopen(fam,"r");
  if (f==NULL) {
      printf("Can't open file %s\n",fam);
      exit(12);
  }
  num =0; 
  while ((res=fscanf(f,"%s %s%[^\n]\n",fid,iid,dummy))  != -1) {
    orig_fam[num]=malloc(128);
    sprintf(orig_fam[num],"%s\t%s",fid,iid);
    num++;
  }
  fclose(f);
  return num;
}  



int * getClusters(char * cluster_fname, int *num) {
  FILE *f;
  int *cluster, d, nrd,j,m, dummy[8192];
  char c;
  *num=0;
  f = fopen(cluster_fname,"r");
  if (f==NULL) {
      printf("Can't open file %s\n",cluster_fname);
      exit(13);
  }
  while (fscanf(f,"%[^\n]\n",dummy)  != -1) (*num)++;
  fclose(f);
  cluster = (int *) malloc(sizeof(int)*(*num)*max_cluster);
  f = fopen(cluster_fname,"r");  
  for (int i=0; i<*num; i++) {
     j=0;
      do {
	 m=fscanf(f,"%d%[\n]",&d,dummy);
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
	     
  

void generateNewVariant(bed_file * new_bed,int *clusters, int v, int *rnds, int gen,
			             int * choice) {
   int    src_idx, src_pr, trg_idx, trg_pr, src;
   int   shift;
   unsigned char *old, *trg, data;

   updateChoices(choice, v, clusters, rnds, gen);
   // index into the BED file for this variant's block
   old = old_bed.bed  + (long) v*getBlockSize(&old_bed); 
   trg = new_bed->bed + (long) v*getBlockSize(new_bed); 
   for(int i=0; i<new_bed->num_s; i++) {
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
      shift = (trg_pr-src_pr)<<1; // *2 
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

int * generateSamples(bed_file * new_bed, int * clusters) {
   int gen, *rnds, i, next,swap, v;
   int choice[new_bed->num_s]; // for each person who the current source is 

   bzero(choice, new_bed->num_s*sizeof(int));

   gen = new_bed->num_s < num_clusters ? new_bed->num_s : num_clusters;
   rnds  = shuffled(num_clusters);
   for (v=0; v<new_bed->num_v; v++)  {
      generateNewVariant(new_bed, clusters, v, rnds, gen, choice);

   }
   return rnds;
}

void outputData(char * base, bed_file * new_bed, int * is_case) {
  FILE *g;
  char outfname[4096];
  int   phe;

  strcpy(outfname,base);
  strcat(outfname,".bed");
  g = fopen(outfname,"w");
  fwrite(new_bed->magic,1,3,g);
  fwrite(new_bed->bed,1,getBlockSize(new_bed)*new_bed->num_v,g);
  fclose(g);
  strcpy(outfname,base);
  strcat(outfname,".fam");
  g = fopen(outfname,"w");
  for(int i=0; i<new_bed->num_s; i++) {
     phe = 1+is_case[i];
     fprintf(g,"%s\t0\t0\t0\t%d\n",orig_fam[i],phe);
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

int * mutateData(char * mute_file, bed_file * new_bed) {
   /* mute_file contains a list of mutations -- each line a triple
    * index in bim file, prob of mutation in case, prob of mutation in control */
   FILE *f;
   int * rnd, v, *is_case;
   int    p_index; // which byte in the block
   int    p_pr=0;     // which two bit of  the byte
   float  case_p, cntl_p, prob, prob_cut;
   unsigned char  data;

   is_case = malloc(new_bed->num_s*sizeof(int));
   bzero(is_case, new_bed->num_s*sizeof(int));
   rnd = shuffled(new_bed->num_s);
   // make the first half cases
   for(int i=0;i<new_bed->num_s/2;i++) {
      is_case[rnd[i]]=1;
   }
   f = fopen(mute_file, "r");
   while (fscanf(f,"%d %f %f\n", &v, &case_p, &cntl_p) == 3) {
      data =  0;
      p_pr = 0;
      for(int p=0; p<new_bed->num_s; p++) {
	 prob_cut = is_case[p] ? case_p : cntl_p;
	 p_index = p>>2;
	 // Add to data -- the lowest two-bits in the byte contain the
	 // lowest numbered individual
	 data =  data  | (get_prob(prob_cut) << (2*p_pr));
	 p_pr++;
         if (p_pr == 4) {
	   long bptr = (long) 3+((long) v)*((long) getBlockSize(new_bed))+p_index;
  	    new_bed->bed[bptr] = data;
	    p_pr=0;
	    data = 0;
	 }
      }
   }
   return is_case;
}


int main(int argc, char *argv[]) {


  int  *clusters, *is_case;
  char  bed[4096], bim[4096], fam[4096];
  strcpy(bed, argv[1]);
  strcpy(bim, argv[1]);
  strcpy(fam, argv[1]);  
  strcat(bed, ".bed");
  strcat(bim, ".bim");
  strcat(fam, ".fam");  



  new_bed.num_s   = atoi(argv[2]);;
  old_bed.num_s   = getFam(fam);


  old_bed.num_v   = new_bed.num_v = getBim(bim, chrom_start);


  readBed(&old_bed,bed);

  allocBed(&new_bed);

  clusters           = getClusters(argv[3],&num_clusters);
  generateSamples(&new_bed, clusters);
  is_case = mutateData(argv[4],&new_bed);
  outputData(argv[5], &new_bed, is_case);
 
}
