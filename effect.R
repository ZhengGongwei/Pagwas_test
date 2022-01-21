library("dplyr")
library("foreach")
library("data.table")
library("Seurat")
library("Matrix")
library(stringr) 
#library("parallel")
# library('doParallel')
library("irlba")
library("glmnet")
#library(tidyfst)
library(GenomicRanges)
library(utils)

library(ggplot2) 
library(ggthemes) 
library(ggpubr)
load("timeEffect_100W_p0.01_a_0.RData")
Path<-c(1,2,3,4,5,6)
features<-c(1,2,3,4,5)

re_df_effect<-lapply(Path,function(i){
	re_list<- lapply(features,function(j){
		effectl<-unlist(lapply(1:20,function(n){
				c1 <- CSPF_timeEffect[[i]][[j]][[21]][[1]]$annotation
				c2 <- CSPF_timeEffect[[i]][[j]][[21]][[1]]$bt_value
				df1 <- data.frame(c1,c2)
				colnames(df1) <- c("annotation","bt_value")
				c3 <- CSPF_timeEffect[[i]][[j]][[n]][[1]]$annotation
				c4 <- CSPF_timeEffect[[i]][[j]][[n]][[1]]$bt_value
				df2 <- data.frame(c3,c4)
				colnames(df2) <- c("annotation","bt_value")
				df <- merge(df1,df2,by="annotation")
			cor(df$bt_value.x,df$bt_value.y,method="spearman")
			}))
		return(effectl)
		})
	re_df<-as.data.frame(re_list)
	colnames(re_df)<- paste0("Features",2*features,"K")
	re_df$Path<- as.character(rep(50*i,nrow(re_df)))
	re_df$seed<-1:20
	return(re_df)
})

re_gg_effect<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),re_df_effect)
re_gg_effect<- reshape2::melt(re_gg_effect,id.vars=c("Path","seed"))
re_gg_effect$Path<-as.character(re_gg_effect$Path)
re_gg_effect

re_gg_effect$Path<-factor(re_gg_effect$Path,levels = c("50","100" ,"150","200","250","300"))
p4<- ggplot(re_gg_effect,aes(x=Path,y=value,fill=variable))+
#geom_violin(trim=FALSE,color="white")+
geom_boxplot(outlier.size=0.1,width= 0.3,position=position_dodge(0.7),alpha=0.6)+ #箱式图异常值大小调整
geom_jitter(color="black" ,size=0.5,alpha=0.3)+theme(legend.position="none")+ #不需要图例
theme_bw()+labs( x = "Random Number of Pathways", y="effect",title = "effect_100W_p0.01_0")
pdf("effect_100W_p0.01_0.pdf")
print(p4)
dev.off()