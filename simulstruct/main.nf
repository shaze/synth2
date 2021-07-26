


outname=params.out
source = params.source
out_dir=params.out_dir
expansion = params.expansion

coreCh = Channel
    .fromFilePairs("${source}.{bed,bim,fam}",size:3) { f -> f.simpleName }

coreCh.into { core1Ch; core2Ch; core3Ch }

process prune {
   input:
      set val(base), file(plinks) from core1Ch
   output:
      set val(base),
          file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into (pCh1, pCh2)
   script:
      out = "${base}-prune"  
   """
     plink --autosome --allow-no-sex --bfile $base --indep-pairwise 50 10 0.2 --out tmp
     plink --allow-no-sex --bfile $base --extract tmp.prune.in \
           --make-bed --out $out
   """
}


process makePCA {
   cpus 4
   input:
      set val(base), file(bed), file(bim), file(fam) from pCh1
   output:
      file("${base}.eigenvec.var") into evec_ch
   script:
     simple=bed.simpleName
    """
      plink --threads 4 --bfile $simple  --pca tabs var-wts --out $base
    """
}


process makeGenome {
   cpus 4
   input:
      set val(base), file(bed), file(bim), file(fam) from pCh2
   output:
      file ("${base}.genome.gz") into genomeCh
   script:
      simple=bed.simpleName
   """
      plink --threads 4 --bfile $simple --genome gz --out $base
   """

}


process makeClusters {
    cpus 4
    input:
       set val(base),  file(plinks) from core2Ch
       file(genome) from genomeCh
    output:
       file ("cluster.cluster2") into clusterCh
    script:
    """
       plink --bfile $base --read-genome $genome --cluster --K ${params.K} --out cluster
    """
}

process massageClusters {
   cpus 1
   input:
      file(cluster) from clusterCh
   output:
      file("cluster_*") into refined_cluster_ch mode 'flatten'
   script:
   """
      refine_cluster.py $cluster 60
   """

}

ready_data_ch = core3Ch.combine(refined_cluster_ch).map { data -> [data[0],data[1],data[2]] }



process extractCluster{
  input:
     set val(base), file(plinks), file(cluster) from ready_data_ch
  output:
     set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into cluster_data_ch
  script:
     out = cluster.simpleName
     """
     plink --allow-no-sex --bfile $base --keep $cluster --make-bed --out $out
     """
}

process simulateCluster{
   input:
      set file(bed), file(bim), file(fam), file(evec) from cluster_data_ch.combine(evec_ch)
   output:
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into simulated_data_ch
   script:
      input=bed.simpleName
      out="${input}_sim"
      """
      N=`getsize.py $fam $expansion`
      hapsimul --fam-prefix $input --num-gen \$N --pca $evec 1600 -H 128 $input $out
      """
}


process combinedPlink {
     input :
        file (all_plinks) from simulated_data_ch.flatten().toList()
     output:
       set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  \
           pca_out_ch
       set file("${out}.bed"), file("${out}.bim"), file("tmp.fam") into  \
           fam_munge_ch	   
     publishDir params.out_dir, pattern: "${out}.{bed,bim}"
     script:
       out=outname
       """
        ls *bed > bed.lst
	ls *bim > bim.lst
        ls *fam > fam.lst
        paste bed.lst bim.lst fam.lst > merge.lst
        plink --allow-no-sex --merge-list merge.lst --make-bed  --out $out
        cp ${out}.fam tmp.fam
      """
 }


process mungeFam {
  input:
    set file(bed), file(bim), file(fam) from fam_munge_ch
  output:
     file("${out}.fam") into res_ch
  publishDir params.out_dir
  script:
   base = bed.simpleName
   out=outname
   """
     mkdir res
     plink --bed $bed --bim $bim --fam $fam --make-bed --impute-sex --out res/tmp
     rename_fam.py --sex-error 0.001 res/tmp.fam ${outname}.fam 
   """
}



process pruneOut {
   input:
      file(plinks) from  pca_out_ch
   output:
      set  file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into final_prune_ch
   script:
      base=plinks[0].simpleName
      out = "${base}-prune"  
   """
     plink --autosome --allow-no-sex --bfile $base --indep-pairwise 50 10 0.2 --out tmp
     plink --allow-no-sex --bfile $base --extract tmp.prune.in \
           --make-bed --out $out
   """
}


process makePCAOut {
   cpus 4
   input:
      set file(bed), file(bim), file(fam) from final_prune_ch
   output:
      file("*eigen*")
   publishDir params.out_dir, pattern: "${simple}*.eig*"     
   script:
     simple=bed.simpleName
    """
      plink --threads 4 --bfile $simple  --pca tabs var-wts --out $simple
    """

}

