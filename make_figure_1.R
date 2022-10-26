library(data.table)
library(ggplot2)

# Run from project directory.
# Dotplot script and signature attribution results expected to be present as shown.
dotplot_script = "scripts/signature_stability_dotplot.R"
gn_counts = fread("results/gunnarsson_MP_attribution_counts.txt")

fusion_labels = c('rr' = 'RUNX1-RUNX1T1', 'other' = 'Other fusions', 'negative' = 'No fusions')
fusion_colors = c("#eec900","#bebada","#87ceeb") # rr, other, negative

gn_counts[, total := rowSums(.SD), .SDcols = patterns("SBS")]
gn_counts[, prop_SBS18 := SBS18 / total]

set.seed(912) # for reproducible jitter
gn_boxplot = ggplot(gn_counts, aes(x = fusion, y = prop_SBS18)) + 
  geom_boxplot(aes(fill = fusion), outlier.shape = NA) +
  geom_jitter(width = .1, size = 2, color = 'grey39') +  
  scale_fill_manual(breaks = c('rr', 'other', 'negative'), values = fusion_colors) +
  scale_x_discrete(limits = c('rr', 'other', 'negative'), 
                   labels = fusion_labels) + 
  scale_y_continuous(limits = c(NA, 0.4)) +
  xlab("Gene-fusion group") +
  ylab("Proportion of SNVs attributed to SBS18") +
  theme_classic() + theme(legend.position = "none")+
  theme(plot.margin = margin(0.2,0.2,1,0.8, "cm"))

gn_boxplot<-gn_boxplot+theme(axis.text=element_text(size=24),
                             axis.title=element_text(size=28))

# Significance testing (Mann-Whitney U)
wilcox.test(gn_counts[fusion == 'rr', prop_SBS18],
            gn_counts[fusion == 'other', prop_SBS18],
            alternative = "two.sided")$p.value
# P = .007


wilcox.test(gn_counts[fusion == 'rr', prop_SBS18],
            gn_counts[fusion == 'negative', prop_SBS18],
            alternative = "two.sided")$p.value
# P = .006



# Run code to generate dotplot (saved to "fig" variable)
source(dotplot_script)
dotplot = fig

# Put some stuff together or legend (some legend material that is shared with figure 2)
repsss<-c("RUNX1::RUNX1T1", "Other fusions", "No fusions")
dummyDat<-data.frame(matrix(ncol=3,nrow = 3,data = (1/12)))
colnames(dummyDat)<-c("dumX","dumY","Gene fusion phenotype")
dummyDat$`Gene fusion phenotype`<-factor(repsss,levels = repsss)

fusion_legend<-ggplot(data = dummyDat,
                aes(x=as.character(dumX),
                    y=dumY,
                    fill=`Gene fusion phenotype`))+
  geom_bar(stat="identity",color="black",show.legend = T)+
  scale_fill_manual(values = fusColors,limits=force, 
                    guide = guide_legend(title = 'Gene fusion\nphenotype'))+
  theme(legend.title=element_text(size=19),legend.text=element_text(size=15),
        legend.justification = 'left')

fusion_legend = get_legend(fusion_legend)


plots = plot_grid(gn_boxplot, dotplot, rel_widths = c(1.74, 1), align = 'hv', axis = 'b', 
                  nrow = 1,  hjust = c(-.5, 0.35), label_size = 32, labels = c("A","B"))
legends = plot_grid(meanLegend, bootLegend, fusion_legend, nrow=3,
                    rel_heights = c(1,-.4,1,-.4,1))

fig1 = plot_grid(plots, legends, nrow = 1, rel_widths = c(6.85, 1))
                 
                 
save_plot("Figure1.png", fig1, base_height = 10.2, base_width = 28, bg = 'white')


# For reporting cumulative attribution of SBS1, SBS5, SBS18, SBS40,
# read in gn_counts again and calculate.
gn_counts = fread("results/gunnarsson_MP_attribution_counts.txt")
gn_counts[, total := rowSums(.SD), .SDcols = patterns("SBS")]
sbs_cols = names(gn_counts)[names(gn_counts) %like% 'SBS']

# Demonstrating that it's fine to show just SBS1, SBS5, SBS18, SBS40 in signature plots.
# (They cover 87.5% of attribution, and next highest accounts for 3%.)
round(sort(decreasing = T, 
           colSums(gn_counts[, .SD, .SDcols = patterns('SBS')]))/ sum(gn_counts$total), 3)[1:6]
# SBS5   SBS1  SBS40  SBS18  SBS19 SBS10c 
# 0.362  0.280  0.118  0.115  0.026  0.013 


