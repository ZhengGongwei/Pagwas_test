library(scPagwas)
library(Seurat)

Args <- commandArgs(T)
gwas = print(Args[1])

tcell <- readRDS("/share2/pub/zhenggw/zhenggw/tcell_scPagwas/tcell.rds")
Idents(tcell) <- "leiden_name_0418"

indexs <- strsplit(gwas, "/", fixed= T)
indexs <- strsplit(indexs[[1]][8], ".", fixed= T)
index= indexs[[1]][1]

Pagwas<-scPagwas_main(Pagwas =NULL,
                   gwas_data = gwas,
                   Single_data = tcell,
                   output.prefix = index,
                   output.dirs="result1",
                   seruat_return=T,
                   Pathway_list=Genes_by_pathway_kegg,
                   ncores=5,
                   assay="RNA",
                   singlecell=T, 
                   celltype=T,
                   block_annotation = block_annotation,
                   chrom_ld = chrom_ld)

filename <- paste(index,".RData",sep='')
save(Pagwas,file = filename)