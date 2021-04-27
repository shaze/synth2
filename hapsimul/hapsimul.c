

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
typedef  unsigned char u_char;

#include <sys/types.h>
#include <stdlib.h>

#define max_cluster 3



static char * allowed_options="e:n:h:";

static struct option long_options[] = {
  {"error-rate",required_argument,0,'e'},
  {"num-gen",required_argument,0,'n'},
  {"haplo-len",required_argument,0,'H'},
  {"gen-pheno",required_argument,0,'g'},  
  {0,0,0,0}
};

/* Global variables -- set by default and perhaps program options and never modified */

static double opt_err_rate=0.001;     // New errors introduced
static int   opt_num_gen=1000;    // Number of individuals generated
static int   opt_hap_len  = 50;       // Size of haplo window
static char   *open_gen_pheno = 0;

static char * input_fname, * output_fname;



void process_options(int argc, char * argv[]) {
  int opt;
  char outfname[255];

  while ((opt = getopt_long(argc,argv,allowed_options, long_options,0)) != -1) {
    switch (opt) {
    case 'e':
      opt_err_rate = atof(optarg);
      break;
    case 'n':
      opt_num_gen=atoi(optarg);
      break;
    case 'h' :
      opt_hap_len=atoi(optarg);
      break;
    }
  }
  if (argc-optind != 2) {
    printf("We need an input and output file -- check arguments carefully "
	   " (either you didn't specify this, or forgot to give a value to one "
	   " of the other arguments");
    exit(-10);
  }
  input_fname=argv[optind];
  output_fname=argv[optind+1];
  if (strcmp(input_fname,output_fname)==0) {
    printf("Input and output file names cannot be the same!");
    exit(-11);
  }
}

static int  num_s, num_v, num_clusters, bsize, new_bsize;
static unsigned char * old_bed;

// for extracting out the bits in a byte (it's in binary to help visualisation)
static int mask[4] = {0b000000011,0b00001100,0b00110000, 0b11000000};

static int chrom_start[27];
  

int getBim(char * bim, int *chrom_start) {
  // get the start position of each chromosome in the bim file
  // if there aren't any in that region, set to -1
  FILE *f;
  float cm;
  int chrom, pos, num_read; // data in bim file
  int next_chrom, line_num;
  char rs[512], maj[32100], min[32100];

  
  f = fopen(bim,"r");
  num_read = fscanf(f,"%d%s%f%d%s%s\n",&chrom, rs, &cm, &pos, maj, min);
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
      num_read = fscanf(f,"%d%511s%f%d%32000s%32000s\n",
			&chrom, rs, &cm, &pos, maj, min);
      if ((strlen(min)>=31998) || (strlen(maj)>=31998)){
	printf("There is a allele of greater than length 32000 in the bim file (line %d)\n",
	       line_num);
	exit(-13);
      }
      if ((strlen(rs)>510)) {
	printf("There is a SNP ID  of greater than length 500 in the bim file");
	exit(-14);
      }
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
    unsigned char magic[3], *bdata;

    bf = fopen(bed, "r");
    fread(magic, 1, 3, bf);
    assert ( (magic[0]==108) && (magic[1]==27) && (magic[2]==1));
    bdata = calloc(num_v,bsize);
    fread(bdata, bsize, num_v, bf);
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


  

void generateNewPerson(unsigned char * new_bed, int pers) {
  int    src_idx, src_pr, trg_idx, trg_pr, src, v;
   int   shift;
   unsigned char *old, *trg, data;

   v=0;
   while (v<num_v) {
     src=arc4random_uniform(num_s); // pick source for current haplo
     src_idx = src>>2;        // byte of old person
     src_pr  = src & 0x03;    //  where in the byte 
     for(int h=0; h<opt_hap_len ; h++) {
	 if (v==num_v) return;
	 old = old_bed + (long) v*bsize;
	 trg = 3+ new_bed+ (long) v*new_bsize; // (3 is for magic numbers)
	 trg_idx = pers>>2;            // byte for the new person 
	 trg_pr  = pers & 0x03;       // where in the byte
	 data     = old[src_idx] & mask[src_pr];
	 // put in right place -- probably need to shift it to left or right
	 // as unlikely to be in same position in the byte
	 // need to shift by multiples of two
	 shift = (trg_pr-src_pr)<<1; // *2 
	 if (shift>0)
	     trg[trg_idx] |= (data << shift);
	 else
	     trg[trg_idx] |= (data >> -shift);
	 v++;
       }
   }
}


void generateSamples(unsigned char * new_bed) {
   int  pers;
   for (pers=0; pers<opt_num_gen; pers++) 
     generateNewPerson(new_bed, pers);
}

void outputData(char * base, unsigned char * new_bed, char * is_case) {
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
  for(int i=0; i<opt_num_gen; i++) {
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

int * shuffled(int n) {
  return malloc(n*sizeof(int));
}

char  * mutateData(char * mute_file, unsigned char * new_bed, int num_new) {
   /* mute_file contains a list of mutations -- each line a triple
    * index in bim file, prob of mutation in case, prob of mutation in control */
   FILE *f;
   int * rnd, v;
   char *is_case;
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
   while (fscanf(f,"%d %f %f\n", &v, &case_p, &cntl_p) == 3) {
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
	    long bptr = (long) 3+((long) v)*((long) new_bsize)+p_index;
  	    new_bed[bptr] = data;
	    p_pr=0;
	    data = 0;
	 }
      }
   }
   return is_case;
}


int main(int argc, char *argv[]) {

  int  *clusters;
  char *is_case;
  unsigned char * new_bed;
  char  bed[4096], bim[4096], fam[4096];

  process_options(argc, argv);

  // max fname is 4096 including suffix!
  snprintf(bed,4096,"%s.bed%c",input_fname,0);
  snprintf(bim,4096,"%s.bim%c",input_fname,0);
  snprintf(fam,4096,"%s.fam%c",input_fname,0);  

  new_bsize        = (opt_num_gen+3)/4;
  num_s             = getFam(fam);
  bsize               = (num_s+3)/4;
  num_v            = getBim(bim, chrom_start);
  old_bed           = readBed(bed);
  new_bed         = allocBed(opt_num_gen);

  generateSamples(new_bed);

  is_case = malloc(opt_num_gen);
  bzero(is_case,opt_num_gen);
    //mutateData(argv[4],new_bed, num_new);
  
  outputData(output_fname, new_bed, is_case);
 
}
