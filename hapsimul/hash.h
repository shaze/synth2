

#define hash_size (1<<30)


int  hash_find(char ** bim_table, char * snp);
void hash_add(char * snp, int p);
void hash_init();