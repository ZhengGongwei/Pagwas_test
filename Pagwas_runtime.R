#time
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

load("timeEffect_100W_0.001_0.72.RData")
Path<-c(1,2,3,4,5,6)
features<-c(1,2,3,4,5)
name <- names(CSPF_timeEffect)
re_df_time<-lapply(Path,function(i){
  re_list<- lapply(features,function(j){
    timel<-unlist(lapply(1:20,function(n) CSPF_timeEffect[name[i]][[1]][[j]][[n]][[2]]))
    return(timel)
    })
re_df<-as.data.frame(re_list)
colnames(re_df)<- paste0("Features",2*features,"K")
re_df$Path<- as.character(rep(50*i,nrow(re_df)))
re_df$seed<-1:20

return(re_df)
   })
#re_df_time
re_gg_time<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),re_df_time)
re_gg_time<- reshape2::melt(re_gg_time,id.vars=c("Path","seed"))
re_gg_time$Path<-as.character(re_gg_time$Path)
re_gg_time

for(i in 1:600)
{
	if(re_gg_time$value[i] > 30)
	{
		re_gg_time$value[i] <- re_gg_time$value[i]/60
	}
}

re_gg_time$Path<-factor(re_gg_time$Path,levels = c("50","100" ,"150","200","250","300"))
p4<- ggplot(re_gg_time,aes(x=Path,y=value,fill=variable))+
#geom_violin(trim=FALSE,color="white")+
geom_boxplot(outlier.size=0.1,width= 0.3,position=position_dodge(0.7),alpha=0.6)+ #箱式图异常值大小调整
geom_jitter(color="black" ,size=0.5,alpha=0.3)+theme(legend.position="none")+ #不需要图例
theme_bw()+labs( x = "Random Number of Pathways", y="time.minutes",title = "100W_p0.001_0.72")
pdf("100W_p0.001_0.72.pdf")
print(p4)
dev.off()
