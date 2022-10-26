library(MutationalPatterns)
library(ggplot2)
library(cowplot)
library(grid)
library(dplyr)
library(plyr)
library(gridExtra)
library(ggtext)

set.seed("20221331")#set seed from year + time of day


#gets the order of the sample ID we want to graph (by SBS18 content)
getID_order<-function(){
  #read in data
  gunMP<-read.csv(sep="\t",file="results/gunnarsson_MP_attribution_counts.txt") #in principle don't need this file, but easier
  gunCEW<-read.csv(sep="\t",file="results/population_scaled_effects.txt")#main data from Mandell
  allPatients<-gunMP$Unique_Patient_Identifier #get list of all patients
  sigsWeCareAbout<-c("SBS1","SBS18","SBS5")
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
  
  #because the way I set up the initial structure is dumb, must remove the NA row (first)
  allFig2ADat<-allFig2ADat[-1,]
  
  #order the plots by CEW for Signature 18, then name all the fusions so they're consistent across all data/plots
  toGetId<-allFig2ADat[which(allFig2ADat$WeightType=="CEW"),]
  ID_order<-toGetId[order(-toGetId$SBS18),"ID"]
  return(ID_order)
  
}


sigColors<-c("#000000","#999999","#cccccc","#ff6a6a")

getJustStatus<-read.csv("results/gunnarsson_MP_attribution_counts.txt",sep="\t")[,c(1,2)]
observedSigs<-unique(read.csv(file="results/effect_share_by_signature.txt",sep="\t")$signature)

fit_res_boot = as.matrix(read.table('results/raw_mp_out.txt'))

fit_rel<-fit_res_boot/rowSums(fit_res_boot)
contri_tb <- fit_rel %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
  tidyr::gather(key = "sig", value = "contri", -exp) %>% 
  dplyr::mutate(sample = gsub("_[^_]+$", "", exp), sample = factor(sample, 
                                                                   levels = unique(sample)), sig = factor(sig, levels = unique(sig)))

ID_order<- getID_order()
getJustStatus[order(match(getJustStatus$Unique_Patient_Identifier,ID_order)),]
fusColors<-c("#eec900","#bebada","#87ceeb") #rr, other, none (yellow,purple, blue)
fill_order<-getJustStatus[order(match(getJustStatus$Unique_Patient_Identifier,ID_order)),"fusion"]
fill_order[fill_order=="rr"]<-fusColors[1]
fill_order[fill_order=="other"]<-fusColors[2]
fill_order[fill_order=="negative"]<-fusColors[3]

contri_boots <- fit_res_boot/rowSums(fit_res_boot)

nr_sigs <- length(unique(contri_tb$sig))
contri_tb3 <- contri_tb %>% dplyr::group_by(sample, sig) %>% 
  dplyr::summarise(mean = mean(contri[contri != 0]), 
                   percentage = sum(contri != 0)/dplyr::n()) %>% 
  dplyr::ungroup() %>% dplyr::filter(!is.na(mean)) %>% 
  dplyr::mutate(sample = factor(sample, levels = rev(levels(sample))))

max_dot_size <- dplyr::case_when(nr_sigs >= 40 ~ 5, nr_sigs >= 
                                   30 ~ 7, nr_sigs >= 20 ~ 8, nr_sigs >= 10 ~ 10, TRUE ~ 15)

contri_tb4<-NULL
for(i in rev(ID_order)){
  if(is.null(contri_tb4)){
    contri_tb4<-contri_tb3[which(contri_tb3$sample==i),]
  } else{
    contri_tb4<-bind_rows(contri_tb4,contri_tb3[which(contri_tb3$sample==i),])
  }
}
contri_tb4$sample<-factor(contri_tb4$sample,levels=rev(ID_order))
tb4_sigColor<-NULL
for(i in levels(contri_tb4$sig)){
  if(i=="SBS1"){
    tb4_sigColor<-c(tb4_sigColor,sigColors[2])
  }else if(i=="SBS5"){
    tb4_sigColor<-c(tb4_sigColor,sigColors[3])
  }else if(i=="SBS18"){
    tb4_sigColor<-c(tb4_sigColor,sigColors[4])
  }else{
    tb4_sigColor<-c(tb4_sigColor,sigColors[1])
  }
  
}

#Remove misc signatures for clarity
toKeepSigs<-unique(levels(contri_tb4$sig))
toKeepSigs<-c("SBS1","SBS5","SBS18","SBS40")
contri_tb4<-contri_tb4[which(contri_tb4$sig%in%toKeepSigs),]


#make incomplete versions of the figure to save their legend to use later
fig <- ggplot(contri_tb4, aes(x = sig, y = sample)) + 
  geom_point(aes(size = mean), shape = 21, fill = 'black') +
  scale_size_continuous(range = c(1, max_dot_size)) + 
  labs(size = "Mean contribution", fill = "Detection\nfrequency", y = "Gunnarsson sample") +
  theme(legend.title=element_text(size=19),legend.text=element_text(size=15),
        legend.justification = 'left')

meanLegend<-get_legend(fig) 

fig <- ggplot(contri_tb4, aes(x = sig, y = sample)) + 
  geom_point(aes(fill = percentage), shape = 21) +
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,1)) + 
  scale_size_continuous(range = c(1, max_dot_size)) + 
  labs(size = "Mean contribution", fill = "Detection\nfrequency", y = "Gunnarsson sample") +
  theme(legend.title=element_text(size=19),legend.text=element_text(size=15),
        legend.justification = 'left')

bootLegend<-get_legend(fig)

#make complete version of the figure for later
fig <- ggplot(contri_tb4, aes(x = sig, y = sample)) + 
  geom_point(aes(fill = percentage, size = mean), shape = 21) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,1)) + 
  scale_size_continuous(range = c(1, max_dot_size)) + 
  labs(size = "Mean contribution", fill = "Detection\nfrequency", y = "Gunnarsson sample") +
  theme(legend.title=element_text(size=19),legend.text=element_text(size=15),
        legend.justification = 'left')


fig <- fig + labs(x = "Signature") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,  size = 10, hjust = 1, vjust = 0.5), 
        text = element_text(size = 12), 
        strip.text.y = element_text(size = 8))

fig <- fig + theme(panel.grid.major = element_line(colour = "gray92"))


element_custom <- function(...) {
  structure(list(...), class = c("element_custom", "element_blank"))
}

element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, y=y, gp=gpar(col=element$colour))
  padding <- unit(1,"line")
  rg <- rectGrob(y=y,width=grobWidth(tg)+padding, height=unit(1,"line")+padding, 
                 gp=gpar(fill = element$fill, col=NA, alpha=0.5))
  gTree(children=gList(rg, tg), width=grobWidth(tg) + padding, cl="custom_axis")
}
widthDetails.custom_axis <- function(x) x$width + unit(2,"mm") # fudge

fig<-fig+theme(axis.text.y = element_custom(colour = 1, fill=rev(fill_order)))
fig<-fig+theme(legend.margin=margin(t=0, r=0, b=0, l=0,unit = "cm"))
fig <- fig+theme(legend.position = "none") +
  theme(plot.margin = margin(0.2,0,1,0.2, "cm"))


# need temp directory to store custom label images
label_dir = paste0(tempdir(), '/sig_labs')

if(!dir.exists(label_dir)){dir.create(label_dir)}
for(i in 1:length(levels(contri_tb4$sig))){
  rects <- data.frame(x = 1,
                      colors = tb4_sigColor[i],
                      text = gsub("SBS","",levels(contri_tb4$sig)[i]))
  text_col<-"white"
  if(tb4_sigColor[i]==sigColors[4]){text_col<-"black"}
  if(tb4_sigColor[i]==sigColors[3]){text_col<-"black"}
  re<-ggplot(rects, aes(x, y = 0, fill = colors, label = text)) +
    geom_tile(width = 1.1, height = .8) + # make square tiles
    geom_text(color = text_col, size=12) + 
    scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
    coord_fixed() + # make sure tiles are square
    theme_void() + # remove any axis markings+
    theme(text = element_text(size = 20))
  ggsave(paste0(label_dir, "/",i,"_rec.png"), re, width = 1.1, height = .85)
}


labelss<-paste0("<img src='", label_dir, "/", 1:length(tb4_sigColor),"_rec.png' width='50'  /><br>")
labelss<-setNames(labelss,levels(contri_tb4$sig))

fig <- fig + scale_x_discrete(name = NULL, labels = labelss)
fig<-fig+theme(axis.title.y = element_text(size = 28))
fig<-fig+theme(axis.text.x  = element_markdown(angle=0,hjust=0.5))
fig<-add_sub(fig, "Single Base Signature (SBS)",y  = 0.5, vjust = 0,size = 28)


