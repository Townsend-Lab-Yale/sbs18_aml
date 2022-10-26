library(data.table)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

#define some of the whole-figure aesthetic choices
fusColors<-c("#eec900","#bebada","#87ceeb") # rr, other, none
sigColors<-c("#000000","#999999","#cccccc","#ff6a6a")
repsss<-c("RUNX1::RUNX1T1", "Other fusions", "No fusions") 
longSigNames<-c("Other signatures","Spontaneous deamination of 5-methylcytosine (SBS1)",
                "Unknown, clock-like (SBS5)","Damage by reactive oxygen species (SBS18)")
sigsWeCareAbout<-c("SBS1","SBS18","SBS5")


##Panel A
#read in data
gunMP<-read.csv(sep="\t",file="results/gunnarsson_MP_attribution_counts.txt")
gunCEW<-read.csv(sep="\t",file="results/population_scaled_effects.txt")
allPatients<-gunMP$Unique_Patient_Identifier #get list of all patients
allFig2ADat<-data.frame(matrix(ncol=(length(sigsWeCareAbout)+3))) #create data structure to hold everything
colnames(allFig2ADat)<-c(sigsWeCareAbout,"Other","WeightType","ID") #name columns

#for every patient, get the proportion SW and CEW and store the data in form for ggplot to use.  
for(i in allPatients){
  thisSig<-gunMP[gunMP$Unique_Patient_Identifier==i,]
  thisCEW<-gunCEW[gunCEW$Unique_Patient_Identifier==i,]
  thisCEW<-thisCEW[,grepl(pattern = "nrsi",colnames(thisCEW))]
  
  CEWProc<-data.frame(matrix(ncol = (length(sigsWeCareAbout)+1),nrow=1))
  CEWProc[1,]<-rep(0,(length(sigsWeCareAbout)+1))
  colnames(CEWProc)<-c(sigsWeCareAbout,"Other")
  for (j in 1:nrow(thisCEW)) {
    thisVar<-thisCEW[j,]
    colnames(thisVar)<-gsub("_nrsi","",colnames(thisVar))
    for(k in colnames(thisVar)){
      if(k %in%sigsWeCareAbout){
        CEWProc[,k]<-CEWProc[,k]+thisVar[1,k]
      }else{
        CEWProc[,"Other"]<-CEWProc[,"Other"]+thisVar[1,k]
      }
      
    }
  }
  CEWProp<-CEWProc[1,]/sum(CEWProc[1,])
  
  SWProc<-CEWProc
  SWProc[1,]<-0
  for(j in 3:ncol(thisSig)){
    thisCol<-colnames(thisSig)[j]
    if(thisCol%in%sigsWeCareAbout){
      SWProc[1,thisCol]<-(SWProc[1,thisCol]+thisSig[1,j])
    }else{
      SWProc[1,"Other"]<-(SWProc[1,"Other"]+thisSig[1,j])
    }
  }
  SWProp<-SWProc[1,]/sum(SWProc[1,])
  
  fig2ADat<-rbind(SWProp,CEWProp)
  fig2ADat$WeightType<-c("SW","CEW")
  fig2ADat$ID<-i
  allFig2ADat<-rbind(allFig2ADat,fig2ADat)
}

# remove empty placeholder row
allFig2ADat<-allFig2ADat[-1,]

#order the plots by CEW for Signature 18, then name all the fusions so they're consistent across all data/plots
toGetId<-allFig2ADat[which(allFig2ADat$WeightType=="CEW"),]
ID_order<-toGetId[order(-toGetId$SBS18),"ID"]
fill_order<-NULL
for(i in ID_order){fill_order<-c(fill_order,gunMP[gunMP$Unique_Patient_Identifier==i,"fusion"])}
fill_order[fill_order=="rr"]<-fusColors[1]
fill_order[fill_order=="other"]<-fusColors[2]
fill_order[fill_order=="negative"]<-fusColors[3]

##I am making the data for Figure2A 'ploty', putting it in a form for ggplot2. 
plotyFig2ADat<-data.frame(matrix(ncol=4))
colnames(plotyFig2ADat)<-c("ID","WeightType","Weight","Signature")
curInd<-1
for(i in 1:(length(sigsWeCareAbout)+1)){
  sigName<-colnames(allFig2ADat)[i]
  toAdd<-allFig2ADat[,c("ID","WeightType",sigName)]
  toAdd$Signature<-sigName
  colnames(toAdd)[3]<-"Weight"
  plotyFig2ADat<-rbind(plotyFig2ADat,toAdd)
}
#again, remove the dumb column because I did data structure bad.
plotyFig2ADat<-plotyFig2ADat[-1,]

#order our data for plotting by factors w/levels
plotyFig2ADat$Signature<-factor(plotyFig2ADat$Signature,levels=c("Other","SBS1","SBS5","SBS18"))
plotyFig2ADat$ID<-factor(plotyFig2ADat$ID,levels=ID_order)
plotyFig2ADat$WeightType<-factor(plotyFig2ADat$WeightType,levels=c("SW","CEW"))

#make the base plot for panel a
pa<-ggplot(plotyFig2ADat, 
       aes(x=WeightType,y=Weight,fill = Signature)) + 
  geom_bar(stat = "identity",color="black",show.legend = T) + 
  theme_classic() + 
  facet_wrap(~ID,nrow=1) + 
  scale_fill_manual(values = sigColors,limits=force, labels=longSigNames) + 
  labs(y="Weight proportion", x="Signature weights and cancer effect weights for Gunnarsson WGS samples") +
  theme(axis.text = element_text(size = 7)) +
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10))#+   

#get the legend for use in the bottom part of the whole figure
legend1<-get_legend(pa)
#remove the legend now that we have copied it.
pa<-pa+ theme(legend.position = "none")

#iterate over all of the facet panels, coloring them by fusion status
ga<-ggplot_gtable(ggplot_build(pa))
stripr<-which(grepl('strip-t',ga$layout$name))
fills<-fill_order
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', ga$grobs[[i]]$grobs[[1]]$childrenOrder))
  ga$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#plot(ga)

#Panel B
#mostly code from Mandell, with some edits here and there
effects_by_sig = fread('results/effect_share_by_signature.txt')
effects_by_sig = effects_by_sig[cohort == 'gunnarsson'] # Just Gunnarsson WGS samples (not TCGA)
#get names to be consistent with other bits of the plot
effects_by_sig$fusion<-gsub("rr","RUNX1::RUNX1T1",effects_by_sig$fusion)
effects_by_sig$fusion<-gsub("other","Other fusions",effects_by_sig$fusion)
effects_by_sig$fusion<-gsub("negative","No fusions",effects_by_sig$fusion)

#order the fusions, so that coloring remains consistent from panel to panel
effects_by_sig$fusion<-factor(effects_by_sig$fusion,levels=repsss)

#plot panel b
pb<-ggplot(effects_by_sig[sigsWeCareAbout, on = 'signature'],  aes(sig_weight, effect_share, color = fusion)) + 
  facet_wrap(~signature) + geom_point(show.legend = F) + 
  xlab('Mutational source weight') + ylab('Cancer effect weight') +
  scale_color_hue(labels=c("RUNX1::RUNX1T1", "Other fusions", "No fusions"), 
                  breaks = c('RUNX1::RUNX1T1', 'Other fusions', 'No fusions')) +
  scale_color_manual(values = fusColors,limits=force) + 
  labs(color = 'gene fusion phenotype') +
  theme_grey() + 
  geom_abline(color = 'darkgrey', linetype = 'dashed')

#color the facet panels by signature
gb<-ggplot_gtable(ggplot_build(pb))
stripr<-which(grepl('strip-t',gb$layout$name))
fills<-sigColors[c(2,4,3)] #order, removing black for "other"
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gb$grobs[[i]]$grobs[[1]]$childrenOrder))
  gb$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Panel C
#get average groupwise signature contributions for Gunnarsson
gun_fusion<-data.frame(matrix(ncol=4))
colnames(gun_fusion)<-c("fusion","signature","effect_share","sig_weight")
tmp = fread('results/effect_share_by_signature.txt')
for(i in unique(tmp$signature)){
  toAddnGun<-tmp[signature == i & cohort == 'gunnarsson', lapply(.SD, mean), by = 'fusion', .SDcols = c('effect_share', 'sig_weight')]
  toAddnGun$signature<-i
  toAddnGun[c(1,4,2,3)]
  gun_fusion<-rbind(gun_fusion,toAddnGun)
}
# Remove empty row
gun_fusion<-gun_fusion[-1,]


# Repeat with the TCGA data
tcga_fusion<-read.csv("results/tcga_fusion_group_effect_shares.txt", sep="\t")

panCa<-data.frame(matrix(ncol=4))
colnames(panCa)<-c("Fusion_ID","WeightType","Weight","Signature")

for(i in 1:nrow(gun_fusion)){
  thisRow<-gun_fusion[i,]
  thisSW<-thisRow$sig_weight
  thisCEW<-thisRow$effect_share
  thisSig<-thisRow$signature
  thisFus<-thisRow$fusion
  panCa<-rbind(panCa,c(thisFus,"SW",thisSW,thisSig))
  panCa<-rbind(panCa,c(thisFus,"CEW",thisCEW,thisSig))
}
panCa<-panCa[-1,]
panCa$Weight<-as.numeric(panCa$Weight)
#get proportion, rather than absolute amounts.
#Ca is the first bit of panel C, Cb is the second (now known as panel D)
#there is a lot of redundant code between the two that should be functionized, but didn't do so.
panCaa<-data.frame(matrix(ncol=4))
colnames(panCaa)<-c("Fusion_ID","WeightType","Weight","Signature")
panCaa<-panCa[which(panCa$Signature%in%sigsWeCareAbout),]
to_proc<-panCa[which(!panCa$Signature%in%sigsWeCareAbout),]
negSW<-sum(to_proc[(to_proc$Fusion_ID=="negative"&to_proc$WeightType=="SW"),"Weight"])
negCEW<-sum(to_proc[to_proc$Fusion_ID=="negative"&to_proc$WeightType=="CEW","Weight"])
rrSW<-sum(to_proc[to_proc$Fusion_ID=="rr"&to_proc$WeightType=="SW","Weight"])
rrCEW<-sum(to_proc[to_proc$Fusion_ID=="rr"&to_proc$WeightType=="CEW","Weight"])
otherSW<-sum(to_proc[to_proc$Fusion_ID=="other"&to_proc$WeightType=="SW","Weight"])
otherCEW<-sum(to_proc[to_proc$Fusion_ID=="other"&to_proc$WeightType=="CEW","Weight"])

panCaa<-rbind(panCaa,c("negative","SW",negSW,"Other Signatures"))
panCaa<-rbind(panCaa,c("negative","CEW",negCEW,"Other Signatures"))
panCaa<-rbind(panCaa,c("rr","SW",rrSW,"Other Signatures"))
panCaa<-rbind(panCaa,c("rr","CEW",rrCEW,"Other Signatures"))
panCaa<-rbind(panCaa,c("other","SW",otherSW,"Other Signatures"))
panCaa<-rbind(panCaa,c("other","CEW",otherCEW,"Other Signatures"))

panCaa$Weight<-as.numeric(panCaa$Weight)

for(i in 1:nrow(panCaa)){
  relRows<-panCaa[(panCaa$Fusion_ID==panCaa[1,"Fusion_ID"]&panCaa$WeightType==panCaa[1,"WeightType"]),]
  panCaa[i,"Weight"]<-(panCaa[i,"Weight"]/sum(relRows$Weight))
}
#name the fusion statuses to be consistent with other plots
panCaa$Fusion_ID<-gsub("rr",repsss[1],panCaa$Fusion_ID)
panCaa$Fusion_ID<-gsub("other",repsss[2],panCaa$Fusion_ID)
panCaa$Fusion_ID<-gsub("negative",repsss[3],panCaa$Fusion_ID)
#use factors to order the data for plotting.
panCaa$WeightType<-factor(panCaa$WeightType,levels=c("SW","CEW"))
panCaa$Fusion_ID<-factor(panCaa$Fusion_ID,levels=repsss)
panCaa$Signature<-factor(panCaa$Signature,levels=c("Other Signatures","SBS1","SBS5","SBS18"))

#make panel c
pc1<-ggplot(data = panCaa,
       aes(x=WeightType,
           y=Weight,
           fill=Signature)) + 
  geom_bar(stat="identity",color="black",show.legend = F) + 
  theme_classic() + 
  facet_wrap(~ Fusion_ID, nrow=1) +
  scale_fill_manual(values = sigColors,limits=force) + 
  guides(fill = "none") + 
  theme_classic() + 
  labs(y="Weight proportion", x= "Signature weights and cancer effect weights for Gunnarsson\ntumors agglomerated by fusion status") 

#color facet panels by fusion type.
gc1<-ggplot_gtable(ggplot_build(pc1))
stripr<-which(grepl('strip-t',gc1$layout$name))
fills<-fusColors
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gc1$grobs[[i]]$grobs[[1]]$childrenOrder))
  gc1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#panel D, same as C but for TCGA data
panCb<-data.frame(matrix(ncol=4))
colnames(panCb)<-c("Fusion_ID","WeightType","Weight","Signature")

for(i in 1:nrow(tcga_fusion)){
  thisRow<-tcga_fusion[i,]
  thisSW<-thisRow$sig_weight
  thisCEW<-thisRow$effect_share
  thisSig<-thisRow$signature
  thisFus<-thisRow$fusion
  panCb<-rbind(panCb,c(thisFus,"SW",thisSW,thisSig))
  panCb<-rbind(panCb,c(thisFus,"CEW",thisCEW,thisSig))
}
panCb<-panCb[-1,]
panCb$Weight<-as.numeric(panCb$Weight)

panCbb<-data.frame(matrix(ncol=4))
colnames(panCbb)<-c("Fusion_ID","WeightType","Weight","Signature")
panCbb<-panCb[which(panCb$Signature%in%sigsWeCareAbout),]
to_proc<-panCb[which(!panCb$Signature%in%sigsWeCareAbout),]
negSW<-sum(to_proc[(to_proc$Fusion_ID=="negative"&to_proc$WeightType=="SW"),"Weight"])
negCEW<-sum(to_proc[to_proc$Fusion_ID=="negative"&to_proc$WeightType=="CEW","Weight"])
rrSW<-sum(to_proc[to_proc$Fusion_ID=="rr"&to_proc$WeightType=="SW","Weight"])
rrCEW<-sum(to_proc[to_proc$Fusion_ID=="rr"&to_proc$WeightType=="CEW","Weight"])
otherSW<-sum(to_proc[to_proc$Fusion_ID=="other"&to_proc$WeightType=="SW","Weight"])
otherCEW<-sum(to_proc[to_proc$Fusion_ID=="other"&to_proc$WeightType=="CEW","Weight"])

panCbb<-rbind(panCbb,c("negative","SW",negSW,"Other Signatures"))
panCbb<-rbind(panCbb,c("negative","CEW",negCEW,"Other Signatures"))
panCbb<-rbind(panCbb,c("rr","SW",rrSW,"Other Signatures"))
panCbb<-rbind(panCbb,c("rr","CEW",rrCEW,"Other Signatures"))
panCbb<-rbind(panCbb,c("other","SW",otherSW,"Other Signatures"))
panCbb<-rbind(panCbb,c("other","CEW",otherCEW,"Other Signatures"))

panCbb$Weight<-as.numeric(panCbb$Weight)

for(i in 1:nrow(panCbb)){
  relRows<-panCbb[(panCbb$Fusion_ID==panCbb[1,"Fusion_ID"]&panCbb$WeightType==panCbb[1,"WeightType"]),]
  panCbb[i,"Weight"]<-(panCbb[i,"Weight"]/sum(relRows$Weight))
}

panCbb$Fusion_ID<-gsub("rr",repsss[1],panCbb$Fusion_ID)
panCbb$Fusion_ID<-gsub("other",repsss[2],panCbb$Fusion_ID)
panCbb$Fusion_ID<-gsub("negative",repsss[3],panCbb$Fusion_ID)

panCbb$WeightType<-factor(panCbb$WeightType,levels=c("SW","CEW"))
panCbb$Fusion_ID<-factor(panCbb$Fusion_ID,levels=repsss)
panCbb$Signature<-factor(panCbb$Signature,levels=c("Other Signatures","SBS1","SBS5","SBS18"))

pc2<-ggplot(data = panCbb,
       aes(x=WeightType,
           y=Weight,
           fill=Signature)) + 
  geom_bar(stat="identity",color="black",show.legend = F) + 
  theme_classic() + 
  facet_wrap(~ Fusion_ID, nrow=1) +
  scale_fill_manual(values = sigColors,limits=force) + 
  labs(y="Weight proportion", x= "Signature Weights and cancer effect weights for TCGA tumors\nagglomerated by fusion status") 

gc2<-ggplot_gtable(ggplot_build(pc2))
stripr<-which(grepl('strip-t',gc2$layout$name))
fills<-fusColors
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gc2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gc2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# Need legend for gene fusion phenotypes.
# Using panel B gives dots, not boxes.
# Instead, save a dummy graph and never render,
# storing its legend. 
dummyDat<-data.frame(matrix(ncol=3,nrow = 3,data = (1/12)))
colnames(dummyDat)<-c("dumX","dumY","Gene fusion phenotype")
dummyDat$`Gene fusion phenotype`<-factor(repsss,levels = repsss)

lplot2p<-ggplot(data = dummyDat,
                  aes(x=as.character(dumX),
                      y=dumY,
                      fill=`Gene fusion phenotype`))+
  geom_bar(stat="identity",color="black",show.legend = T)+
  scale_fill_manual(values = fusColors,limits=force)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10)) 

legend2<-get_legend(lplot2p)
lplot2p<-lplot2p+theme(legend.title=element_text(size=11),legend.text=element_text(size=9)) 
lplot2p<-lplot2p+theme(legend.margin=margin(l=20))
legend2a<-get_legend(lplot2p)

grid = plot_grid(ga,gb,plot_grid(gc1,gc2,ncol=2,labels=c("C","D")),
                 plot_grid(legend1,NULL,legend2,ncol=3,rel_widths = c(1,-.6,1)),ncol = 1, labels = c("A","B"),
                 rel_heights = c(1,1,1,0.5))
cowplot::save_plot('Figure2.png', grid, base_height = 13, base_width = 12.5, bg = 'white')


