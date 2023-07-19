#loading libraries
library(phyloseq)
library(microbiome)
library(ggthemes) # additional themes fro ggplot2
library(ggpubr)
library(vegan)
library(repr)
library(gridExtra) # gridding plots
library(viridis)
library(reshape2)
library(knitr)
library(microbiome)
library(sjPlot)
set.seed(10000)
#library(dplyr)

options(repr.plot.width=12, repr.plot.height=8)


theme_glab <- function(base_size = 30,
                    base_family = "",
                    base_line_size = base_size / 180,
                    base_rect_size = base_size / 180) {
   
    font <- "Helvetica" #assign font family up front
   
    theme_bw(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
        legend.background =  element_blank(),
        legend.title =       element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65),
                                         hjust = 0),
        legend.text =        element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65)),
        legend.key.size =    unit(0.8, "lines"),
     
      plot.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        hjust = 0),
       
      axis.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
      axis.text = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
       
      plot.caption = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.7),
        hjust = 1),
       
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)),

     
      complete = TRUE
    )
}

data<-readRDS("shark_data.Rds")
rdata<-readRDS("shark_data_relative_abundance.Rds")

shark1<-subset_samples(data, species==c("Prionace_glauca"))
shark2<-subset_samples(data, species==c("Somniosus_rostratus"))
sh<-merge_phyloseq(shark1,shark2)


plot_richness(sh, x="species", measures =c("Shannon")) +
   geom_boxplot(aes(fill= species)) +
   theme_glab()+
     scale_fill_viridis_d(alpha=0.8)+
  ylim(1.5, 5)+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))

plot_richness(sh, x="zone", measures =c("Shannon")) +
   geom_boxplot(aes(fill= species),alpha=0.8) +
   theme_glab()+
     scale_fill_viridis_d(alpha=0.6)+
  ylim(1.5, 3.1)+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  scale_fill_manual(values = c("Prionace_glauca" = "#482677FF",
                                "Somniosus_rostratus"="#95D840FF"))

## Testing for differences in alpha diversity between the Locations
zone_alpha_div <- data.frame(estimate_richness(sh, measures=c("Shannon")), 
                                 data.frame(sample_data(sh)$zone))

kruskal.test(zone_alpha_div$Shannon~zone_alpha_div$sample_data.sh..zone)

## Testing for differences in alpha diversity between the Locations
species_alpha_div <- data.frame(estimate_richness(sh, measures=c("Shannon")), 
                                 data.frame(sample_data(sh)$species))

kruskal.test(species_alpha_div$Shannon~species_alpha_div$sample_data.sh..species)

species_order<-list("Prionace_glauca","Somniosus_rostratus","Water")

plot_richness(data, x="type", measures =c("Shannon")) +
   geom_boxplot(aes(fill= type)) +
   theme_glab()+
     scale_fill_viridis_d(alpha=0.8)+
  ylim(1.5, 5)+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  scale_fill_manual(values = c("Environment" = "#39568CFF",
                                "Shark"="#FDE725FF"),labels = c("Sea water", "Sharks"))

## Testing for differences in alpha diversity between the Locations
species_alpha_div <- data.frame(estimate_richness(data, measures=c("Shannon")), 
                                 data.frame(sample_data(data)$type))

kruskal.test(species_alpha_div$Shannon~species_alpha_div$sample_data.data..type)

plot_richness(data, x="species", measures =c("Shannon")) +
   geom_boxplot(aes(fill= species)) +
   theme_glab()+
     scale_fill_viridis_d(alpha=0.8)+
  ylim(1.5, 5)+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))

alpha <- data.frame("Shannon" = phyloseq::estimate_richness(data, measures = "Shannon"),
"Species" = phyloseq::sample_data(data)$species)

alpha %>%
  group_by(Species) %>%
  dplyr::summarise(mean_shannon = mean(Shannon))

plot_richness(data, x="zone", measures =c("Shannon")) +
   geom_boxplot(aes(fill= zone)) +
   theme_glab()+
     scale_fill_viridis_d(alpha=0.8)+
  ylim(1.5, 5)+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))

## Testing for differences in alpha diversity between the Locations
species_alpha_div <- data.frame(estimate_richness(data, measures=c("Shannon")), 
                                 data.frame(sample_data(data)$species))

kruskal.test(species_alpha_div$Shannon~species_alpha_div$sample_data.data..species)

## Testing for differences in alpha diversity between the Locations
species_alpha_div <- data.frame(estimate_richness(data, measures=c("Shannon")), 
                                 data.frame(sample_data(data)$zone))

kruskal.test(species_alpha_div$Shannon~species_alpha_div$sample_data.data..zone)

data_ra_family = tax_glom(rdata, "Family", NArm = TRUE)
data_ra_order = tax_glom(rdata, "Order", NArm = TRUE)
data_ra_class = tax_glom(rdata, "Class", NArm = TRUE)
data_ra_phylum = tax_glom(rdata, "Phylum", NArm = TRUE)

order<-list("Water","Prionace_glauca","Somniosus_rostratus")

prova<-plot_bar(data_ra_phylum, fill="Phylum") + 
#scale_fill_manual(values=col_vector)+
facet_grid(~species, drop = T, scales = "free_x", space="free_x")+
#theme_glab() + 
#theme_hc(base_size=20) +
labs(y= "Abundance")+
theme(axis.text=element_text(size=18),
      axis.text.x=element_text(angle=20,hjust=1,size=0),
      #axis.ticks.x=element_line(vjust=5),
      strip.text.x=element_text(size=20),
      legend.position="bottom", 
      legend.text=element_text(size=14),
      legend.title=element_text(size=16),
      panel.grid.minor=element_blank(), 
      panel.grid.major=element_blank(),
      panel.background=element_blank(),
      axis.title=element_text(size=25, color = rgb(100, 100, 100, maxColorValue = 255)),
      axis.title.x=element_text(vjust=-0.5)
)

prova$data$species <- factor(prova$data$species , levels = order)
prova

save_plot("prova.svg", fig = prova, width=40, height=34)


plot_bar(data_ra_phylum, fill="Phylum" , x="sample" , title = "Diversity at Phylum level") +
theme_glab()+
theme(axis.text.x=element_text(angle=90))

w<-subset_samples(rdata, species=="Water")
w <- prune_taxa(taxa_sums(w) > 0, w)
table(tax_table(w)[, "Phylum"])


b<-subset_samples(rdata, species=="Prionace_glauca")
b <- prune_taxa(taxa_sums(b) > 0, b)
table(tax_table(b)[, "Phylum"])


k<-subset_samples(rdata, species=="Somniosus_rostratus")
k <- prune_taxa(taxa_sums(k) > 0, k)
table(tax_table(k)[, "Phylum"])

prova<-plot_bar(data_ra_class, fill="Class") + 
#scale_fill_manual(values=col_vector)+
facet_grid(~species, drop = T, scales = "free_x", space="free_x")+
#theme_glab() + 
#theme_hc(base_size=20) +
labs(y= "Abundance")+
theme(axis.text=element_text(size=18),
      axis.text.x=element_text(angle=20,hjust=1,size=0),
      #axis.ticks.x=element_line(vjust=5),
      strip.text.x=element_text(size=20),
      legend.position="bottom", 
      legend.text=element_text(size=14),
      legend.title=element_text(size=16),
      panel.grid.minor=element_blank(), 
      panel.grid.major=element_blank(),
      panel.background=element_blank(),
      axis.title=element_text(size=25, color = rgb(100, 100, 100, maxColorValue = 255)),
      axis.title.x=element_text(vjust=-0.5)
)

prova$data$species <- factor(prova$data$species , levels = order)
prova

save_plot("prova2.svg", fig = prova, width=40, height=34)

table(tax_table(w)[, "Family"])

table(tax_table(b)[, "Family"])

table(tax_table(k)[, "Family"])

# Extract the Top 10 Phyla by relative abundance
top10_family <- sort(taxa_sums(data_ra_family), TRUE)[1:10]

family_top10 <- prune_taxa(names(top10_family), data_ra_family)

prova<-plot_bar(family_top10, fill="Family") + 
#scale_fill_manual(values=col_vector)+
facet_grid(~species, drop = T, scales = "free_x", space="free_x")+
#theme_glab() + 
#theme_hc(base_size=20) +
labs(y= "Abundance")+
theme(axis.text=element_text(size=18),
      axis.text.x=element_text(angle=20,hjust=1,size=0),
      #axis.ticks.x=element_line(vjust=5),
      strip.text.x=element_text(size=20),
      legend.position="bottom", 
      legend.text=element_text(size=14),
      legend.title=element_text(size=16),
      panel.grid.minor=element_blank(), 
      panel.grid.major=element_blank(),
      panel.background=element_blank(),
      axis.title=element_text(size=25, color = rgb(100, 100, 100, maxColorValue = 255)),
      axis.title.x=element_text(vjust=-0.5)
)

prova$data$species <- factor(prova$data$species , levels = order)
prova

save_plot("prova2.svg", fig = prova, width=40, height=34)


# Plot the top 10 families
plot_bar(family_top10, fill="Family", x="sample" , title = "Diversity of the top 10 family") +
    theme_glab()+
theme(axis.text.x=element_text(angle=90))

plot_bar(family_top10, fill="Family", x="zone" , title = "Diversity of the top 10 family") +
    theme_glab()+
theme(axis.text.x=element_text(angle=90))

plot_bar(family_top10, fill="Family", x="species" , title = "Diversity of the top 10 family") +
    theme_glab()+
theme(axis.text.x=element_text(angle=90))





set.seed(1000) # to make the analysis reproducible

prok_dist_wjac <- distance(data, method = "jaccard")
prok_dist_unjac <- distance(data, method = "jaccard", binary = TRUE)

prok_nmds_jw <- ordinate(data,prok_dist_wjac, method = "NMDS",trymax=1000)

prok_nmds_juw <- ordinate(data,prok_dist_unjac, method = "NMDS",trymax=999)

# nmds W jaccard, removing other contaminants
c<-plot_ordination(data, prok_nmds_jw, type="samples") +
#geom_text(aes(label= sample), size=4, hjust=0.4,vjust=2.5) +
geom_point(aes(fill=type),shape=21,size=9,color="black",stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()+
  scale_fill_manual(values = c("Environment" = "#39568CFF",
                                "Shark"="#FDE725FF"),labels = c("Sea water", "Sharks"))
c
#save_plot("c.svg", fig = c, width=40, height=34)


# nmds W jaccard, removing other contaminants
plot_ordination(data, prok_nmds_juw, type="samples",title="nMDS unweighted Jaccard similarity") +
geom_text(aes(label= sample), size=4, hjust=0.2,vjust=2) +
geom_point(aes(fill=species),shape=21,size=7,color="black",stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(data, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(data)$species)

shark1<-subset_samples(data, species==c("Prionace_glauca"))
shark2<-subset_samples(data, species==c("Somniosus_rostratus"))
sh<-merge_phyloseq(shark1,shark2)

shark_dist_wjac <- distance(sh, method = "jaccard")
shark_dist_unjac <- distance(sh, method = "jaccard", binary = TRUE)

shark_nmds_jw <- ordinate(sh,shark_dist_wjac, method = "NMDS",trymax=1000)

shark_nmds_juw <- ordinate(sh,shark_dist_unjac, method = "NMDS",trymax=999)

# nmds W jaccard, removing other contaminants
plot_ordination(sh, shark_nmds_jw, type="samples",title="nMDS weighted Jaccard similarity") +
#geom_text(aes(label= sample), size=4, hjust=0.4,vjust=2.5) +
geom_point(aes(fill=species),shape=21,size=8,color="black",stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()+
  scale_fill_manual(values = c("Prionace_glauca" = "#482677FF",
                                "Somniosus_rostratus"="#95D840FF"))


# nmds W jaccard, removing other contaminants
plot_ordination(sh, shark_nmds_juw, type="samples",title="nMDS unweighted Jaccard similarity") +
#geom_text(aes(label= sample), size=4, hjust=0.2,vjust=2) +
geom_point(aes(fill=species),shape=21,size=7,color="black",stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(sh, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(sh)$species)

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(sh, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(sh)$zone)

plot_ordination(sh, shark_nmds_jw, type="samples") +
geom_point(aes(color= species, shape= zone),size=8,stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
scale_color_viridis_d()+
  scale_fill_manual(values = c("Prionace_glauca" = "#482677FF",
                                "Somniosus_rostratus"="#95D840FF"))+ 
  scale_shape_manual(values = c(15, 16, 17, 18)) 

# nmds W jaccard, removing other contaminants
d<-plot_ordination(sh, shark_nmds_jw, type="samples") +
geom_point(aes(color= species, shape= zone),size=8,stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
scale_fill_viridis_d()+
  scale_color_manual(values = c("Prionace_glauca" = "#482677FF",
                                "Somniosus_rostratus"="#95D840FF"))+ 
  scale_shape_manual(values = c(15, 16, 17, 18)) 
d
#save_plot("d.svg", fig = d, width=40, height=34)


# nmds W jaccard, removing other contaminants
plot_ordination(sh, shark_nmds_jw, type="samples",title="nMDS weighted Jaccard similarity") +
#geom_text(aes(label= sample), size=4, hjust=0.4,vjust=2.5) +
geom_point(aes(color= species, shape= zone),size=8,stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
     scale_color_viridis_d()

# nmds W jaccard, removing other contaminants
plot_ordination(sh, shark_nmds_juw, type="samples",title="nMDS unweighted Jaccard similarity") +
#geom_text(aes(label= sample), size=4, hjust=0.2,vjust=2) +
geom_point(aes(fill=zone),shape=21,size=7,color="black",stroke=0.25) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()





prok_pcoa_jw <- ordinate(data,prok_dist_wjac, method = "PCoA",trymax=100)

prok_pcoa_juw <- ordinate(data,prok_dist_unjac, method = "PCoA",trymax=100)

#PCoA W jaccard
plot_ordination(data, prok_pcoa_jw, type="samples",title="PCoA weighted Jaccard similarity") +
geom_point(aes(fill=species),shape=21,size=8,color="black",stroke=0.25) + 
#scale_fill_manual(values=c("#fb9b06","#fde725","#6ece58","#1f9e89","#31688e","#481668"),
 #                 breaks=c("SPO","SAF","PF","SACCF","SBDY","RSCS")) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()


#PCoA Uw jaccard
plot_ordination(data, prok_pcoa_juw, type="samples",title="PCoA unweighted Jaccard similarity") +
geom_point(aes(fill=species),shape=21,size=7,color="black",stroke=0.25) +
#scale_fill_manual(values=c("#fb9b06","#fde725","#6ece58","#1f9e89","#31688e","#481668"),
 #                 breaks=c("SPO","SAF","PF","SACCF","SBDY","RSCS")) +
theme_glab() + theme(legend.position = "right") +
guides(fill = guide_legend(override.aes = list(shape = 21) ),
            shape = guide_legend(override.aes = list(fill = "black")))+
     scale_fill_viridis_d()

save.image()
