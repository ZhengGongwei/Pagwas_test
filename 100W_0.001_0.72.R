library("dplyr")
library("data.table")
library("Seurat")
library("Matrix")
library("stringr") 
library("parallel")
library("irlba")
library("glmnet")
library(GenomicRanges)
library(utils)
library(ggplot2) 
library(ggthemes) 
library(ggpubr)
library(forestplot)
require(bigstatsr)

source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/GWAS_summary_input.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Link_pathway_blocks_gwas.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Tissue_eqtls_Input4.0.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/link_pwpca_block_nes2.0.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/scPagwas_Umap.r")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/sub_functions.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Pathway_annotation_input.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Pathway_pcascore_run.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Single_data_input.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/scPagwas_perform_inference.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Pagwas_perform_inference.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/Pagwas_main.R")
source("/share/pub/dengcy/GWAS_Multiomics/pagwas/R/bootstrap_Pvalue_Bar.R")

#数据准备
#1.导入处理好的基因注释文件
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/gtf_df.RData")
#2.导入处理好的通路文件
load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.pathway_kegg.RData")
ld_folder="/share/pub/dengcy/Singlecell/COVID19/data/LD"
#single cell
load("/share/pub/dengcy/GWAS_Multiomics/test/brain/mouse_brain_seu2.RData")
#GWAS
load("/share2/pub/zhenggw/zhenggw/Pagwas_test/100W_0.001_0.72.RData")
#"data_SCZ_New_100W_p0.001_a0.72"     "mouse_brain_seu2"

CSPF_timeEffect<- lapply(c(50,100,150,200,250,300),function(num1){#pathway
							lapply(c(2000,4000,6000,8000,10000),function(num2){#nfeatures	

								 		  message("cells:5000")
								 		  message(paste0("Pathways:",num1))
								 		  message(paste0("nfeatures:",num2))

										  seed_result<-lapply(1001:1021,function(i){
										  message(paste0("seed:",i))
										  set.seed(i)

										  pathway<-genes.by.pathway_kegg[sample(1:length(genes.by.pathway_kegg),num1)]

	  									  cc<-sample(1:ncol(mouse_brain_seu2),5000)
										  mouse_brain_test<-mouse_brain_seu2[,cc]

										  timestart<-Sys.time()

										  Pagwas <- Pagwas_main(Pagwas = NULL,
										          gwas_data = data_SCZ_New_100W_p0.001_a0.72,
										          add_eqtls = "FALSE",
										          block_annotation = gtf_df,
										          Single_data = mouse_brain_test,
										          nfeatures = num2,
										          Pathway_list= pathway,
										          ld_folder = ld_folder,
										          simp_results = T,
										          iters = 200
										          )

										  bootstrap_results<- Pagwas$bootstrap_results
										  bootstrap_results<-bootstrap_results[-1,]
										  timeend<-Sys.time()
										  runningtime<-timeend-timestart
										  return(list(bootstrap_results,runningtime))
										  })
										 names(seed_result)<- paste0("Random_",1:20)
										 return(seed_result)
	})
  })

name <-	lapply(c(50,100,150,200,250,300),function(k){#pathway
  				lapply(c("2k","4K","6K","8K","10K"),function(m){#nfeatures	
					    a <- paste0("Cell:5k_SNP:10W","_Path:",k,"_features:",m)
					    print(a)
  })
})

names(CSPF_timeEffect)<-name
save(CSPF_timeEffect,file="/share2/pub/zhenggw/zhenggw/Pagwas_test/timeEffect_100W_0.001_0.72.RData")

