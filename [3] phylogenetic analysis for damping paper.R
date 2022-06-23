## =====================================
## 3, PGLS analyses in damping paper ####
## =====================================
rm(list=ls())
library(ape)
library(phytools)
library(reshape)
library(caper)
library(tidyverse)
library(picante)
library(ggtext)
library(ggpmisc)

# install.packages("/Users/shajiang/Downloads/phytools_0.7-80.tar.gz", repos = NULL, type="source")
setwd("/Users/shajiang/Documents/AAA @Stanford/research/Notes_AgeStructure")
getwd()

full_uni <- read.csv2("./phylogenetic tree/match_list_uncorrected_newversion v2.csv",
                      sep = ",")
all_data <- read.csv2("./output data/full animal and plant data v2.csv")
length(unique(all_data$SpeciesAccepted))
final_tree_read<-read.tree(file="./phylogenetic tree/final_tree_uncorrected_newversion v2.tre")
final_tree_read$tip.label <-  gsub("_", " ",final_tree_read$tip.label)
species_tree = as.vector(unique(final_tree_read$tip.label))

species_name_correct <- as.vector(all_data$SpeciesAccepted)

## check the variance of each variable within species
check <- all_data %>%
  group_by(SpeciesAccepted, db_sep)%>%
  summarise(number = n(),
            Tc.var = var(Tc),
            sigma.var = var(sigma),
            damping.time.var = var(damping.time),
            Tc.mean = mean(Tc),
            sigma.mean = mean(sigma),
            damping.time.mean = mean(damping.time),
            Tc.cv = sqrt(Tc.var)/Tc.mean,
            sigma.cv = sqrt(sigma.var)/sigma.mean,
            damping.time.cv = sqrt(damping.time.var)/damping.time.mean)
summary(check$Tc.mean)
summary(check$sigma.mean)
summary(check$Tc.var)
summary(check$sigma.var)
summary(check$Tc.cv)
summary(check$sigma.cv)
summary(check$damping.time.cv)

hist(log(check$number))
summary(check$number)

unique(all_data$Class)
length(unique(filter(all_data, Class %in% c("Actinopterygii", "Aves", "Mammalia", "Reptilia"))$SpeciesAccepted))

matrix_select<-all_data %>%
  group_by(SpeciesAccepted, db_sep)%>%
  mutate(number = n())

## use the median value of damping time to choose one matrix for each species
median_index = function(x) {
  lx = length(x)
  if (lx %% 2 == 1) {
    return(median(x))
  }
  return(median(c(x,0))) # for even number, add 0 in the end
}

median_data_sep <- all_data %>%
  dplyr::select(SpeciesAccepted, Tc, sigma, damping.time, db_sep, db_taxa,db_source,
                Class, omega, alpha, damping.cal, damping.approx, Order, ProjectionInterval)%>%
  group_by(SpeciesAccepted)%>%
  mutate(sigma.median = median_index(sigma))%>%
  distinct(SpeciesAccepted,.keep_all=TRUE)%>%
  mutate(logTc = log(Tc),
         logsigma = log(sigma),
         logdamping.time = log(damping.time),
         logreptime = log(omega - alpha))%>%
  filter(SpeciesAccepted %in% species_tree)%>%
  as.data.frame()

length(unique(median_data_sep$SpeciesAccepted))
summary(median_data_sep$ProjectionInterval) # NAs are for the GMO dataset,for which the number should be 1
table(median_data_sep$ProjectionInterval)
table(median_data_sep$db_source,
      median_data_sep$ProjectionInterval)

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(ProjectionInterval))

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(SpeciesAccepted))

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(SpeciesAccepted))

# Check if the tree is rooted
tree <- final_tree_read

is.rooted(tree)

# If it is binary
is.binary(tree)

# Make the node labels unique
tree<- makeNodeLabel(tree)

# PGLS analyses -------------------------------------------------------------------------------------------------
#### PGLS function ####
# need to have a grp variable in input data
PGLS_fun = function(data, variable, data_group, outputlist){
  output_combine=data.frame()
  db_combine=data.frame()

  for (d in 1:(length(data_group))) {
    data_filter <- filter(data, grp %in% data_group[d])
    comp_data <- comparative.data(phy = final_tree_read,
                                  data = data_filter,
                                  names.col = SpeciesAccepted,
                                  vcv = TRUE,
                                  vcv.dim = 3)
    rowname_data <- rownames(comp_data$data)
    count=0
    for (expl in 1:(length(variable)-1)) {
      for (resp in (expl+1):(length(variable))) {
        count=count+1
        # take log
        explanatory=as.numeric(as.matrix(comp_data$data[variable[expl]]))
        response=as.numeric(as.matrix(comp_data$data[variable[resp]]))
        mod=pgls(explanatory ~ response, comp_data, lambda='ML', param.CI = 0.95)

        output=c(data_group[d],
                 paste(variable[expl],"~",variable[resp]),
                 summary(mod)$coeff[1,1],
                 summary(mod)$coeff[2,],
                 mod$param[1],
                 mod$param[2],
                 mod$param.CI$lambda$ci.val[1],
                 mod$param.CI$lambda$ci.val[2],
                 mod$mlVals[2],
                 summary(mod)$r.squared,
                 unique(data_filter$db_taxa))

        db = comp_data[["data"]] %>%
          mutate(residual = mod$residual,
                 variable = paste(variable[expl],"~",variable[resp]))

        print(paste("PGLS model: ",variable[expl]," ~ ",variable[resp],sep=""))

        db_combine = rbind(db_combine,db)
        output_combine = rbind(output_combine,output)
      }
    }
  }
  colnames(output_combine) <- c(outputlist[-length(outputlist)])
  output_combine[,"Padjusted"]=p.adjust(output_combine[,"P"],"BH")
  output_combine <- as.data.frame(output_combine)%>%
    mutate_at(c(3:13,15), as.numeric)%>%
    arrange(grp)%>%
    as.data.frame()

  return(list(output_combine, db_combine))
}
variable <- c("logdamping.time", "logsigma", "logTc")

median_data_sep = median_data_sep %>%
  mutate(grp = db_sep)
data_group <- unique(median_data_sep$grp)
# We create an empty lists where we will store the results
outputlist=c("grp","variable","intercept",
             "slope","SE","t","P",
             "Kappa",
             "Pagel lambda","lambda CI lower","lambda CI upper",
             "LambdaOptimization",
             "R2","db_taxa", "Padjusted")
result = PGLS_fun(median_data_sep, variable, data_group, outputlist)
saveTable = result[[1]]%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)

residual = result[[2]]
# Save the results
write.csv(saveTable, file = "./phylogenetic tree/PGLsList_uncorrected by seperate groups (median).csv") # PGLS raw correlations

########## New figures in paper ##########
pos = median_data_sep %>% group_by(db_sep) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.9 * max(log(sigma)),
            damp.pos = 0.9 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.15,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.15)
#### (1) Tc and sigma #####
summ <- saveTable %>%
  group_by(grp) %>%
  filter(variable %in% "logsigma ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

# CI for slope
0.9661649-0.03977229
0.9661649+0.03977229

1.0624987-0.02397057
1.0624987+0.02397057

1.1221496-0.01721676
1.1221496+0.01721676

f_labels <- data.frame(
  db_sep = summ$grp,
  label1 = c(paste0("y = ",summ$slope, "x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_Tc_sigma <- ggplot(median_data_sep,
                     aes(log(Tc), log(sigma)))+
  geom_point(shape = 21, size = 2.5, alpha = 0.7)+
  # scale_shape_manual(values=myshape)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(S)))+
  facet_wrap(. ~db_sep, nrow = 3, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -4, y = pos$sigma.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -4, y = pos$sigma.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)

p_Tc_sigma
ggsave("./plot/Tc and sigma PGLS.png", p_Tc_sigma , width = 6, height = 8)
ggsave("./plot/Figure_1.pdf", p_Tc_sigma , width = 6, height = 8)

 #### (2) Tc and damping #####
summ <- saveTable %>%
  group_by(grp) %>%
  filter(variable %in% "logdamping.time ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

# CI for slope
0.7800803-0.04312197
0.7800803+0.04312197

0.8216637-0.06731833
0.8216637+0.06731833

0.6526575-0.04745337
0.6526575+0.04745337

f_labels <- data.frame(
  db_sep = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_Tc_damping <- ggplot(median_data_sep,
                       aes(log(Tc), log(damping.time)))+
  geom_point(shape = 21, size = 2.5, alpha = 0.7)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(tau)))+
  facet_wrap(. ~db_sep, nrow = 3, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -4, y = pos$damp.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -4, y = pos$damp.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)

p_Tc_damping
ggsave("./plot/Tc and damping PGLS.png", p_Tc_damping , width = 6, height = 8)
ggsave("./plot/Figure_2.pdf", p_Tc_damping , width = 6, height = 8)

########## New figures in appendix ##########
#### (1) Figure by Class #####
# only keep Classes with enough data points. Use PGLS
median_data_sep2 <- median_data_sep%>%
  group_by(Class)%>%
  mutate(num = n(),
         grp = Class)%>%
  dplyr::filter(num >=20)%>%
  ungroup()%>%
  as.data.frame()

variable <- c("logsigma", "logdamping.time", "logTc")
data_group = unique(median_data_sep2$grp)
saveTable_Class = PGLS_fun(median_data_sep2, variable, data_group, outputlist)[[1]]%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)


#### (1.1) Tc and sigma by Class#####
# for animal
db_taxa_select = c("Animal")
pos = median_data_sep2 %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.97 * max(log(sigma)),
            damp.pos = 0.97 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.1,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.1)

summ <- saveTable_Class %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  filter(variable %in% "logsigma ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  Class = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)


p_Tc_sigma_class1 <- ggplot(median_data_sep2 %>%
                              filter(db_taxa %in% db_taxa_select),
                            aes(log(Tc), log(sigma)))+
  geom_point(shape = 21, size = 2.5, alpha = 1)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(S)))+
  facet_wrap( ~Class, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -4, y = pos$sigma.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -4, y = pos$sigma.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)

print(p_Tc_sigma_class1)
ggsave(paste0("./plot/Tc and sigma ",db_taxa_select," facet by Class PGLS.png"), p_Tc_sigma_class1, width = 8, height = 7)
ggsave("./plot/Figure_C1.pdf", p_Tc_sigma_class1 , width = 8, height = 7)

# for plant
db_taxa_select = c("Plant")
pos = median_data_sep2 %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.97 * max(log(sigma)),
            damp.pos = 0.97 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.08,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.08)


summ <- saveTable_Class %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  filter(variable %in% "logsigma ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  Class = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_Tc_sigma_class2 <- ggplot(median_data_sep2 %>%
                              filter(db_taxa %in% db_taxa_select),
                            aes(log(Tc), log(sigma)))+
  geom_point(alpha = 1, shape = 21, size = 2.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(S)))+
  facet_wrap( ~Class, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = 0.4, y = pos$sigma.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = 0.4, y = pos$sigma.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)
print(p_Tc_sigma_class2)
ggsave(paste0("./plot/Tc and sigma ",db_taxa_select," facet by Class PGLS.png"), p_Tc_sigma_class2, width = 8, height = 4.5)
ggsave("./plot/Figure_C2.pdf", p_Tc_sigma_class2 , width = 8, height = 7)


#### (1.2) Tc and damping by Class #####
# for animal
db_taxa_select = c("Animal")
pos = median_data_sep2 %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.97 * max(log(sigma)),
            damp.pos = 0.97 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.1,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.1)


summ <- saveTable_Class %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  filter(variable %in% "logdamping.time ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  Class = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_Tc_damping_class1 <- ggplot(median_data_sep2 %>%
                              filter(db_taxa %in% db_taxa_select),
                            aes(log(Tc), log(damping.time)))+
  geom_point(alpha = 1, shape = 21, size = 2.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(tau)))+
  facet_wrap( ~Class, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -4, y = pos$damp.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -4, y = pos$damp.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)
print(p_Tc_damping_class1)
ggsave(paste0("./plot/Tc and damping ",db_taxa_select," facet by Class PGLS.png"), p_Tc_damping_class1, width = 8, height = 7)
ggsave("./plot/Figure_C4.pdf", p_Tc_damping_class1 , width = 8, height = 7)


# for plant
db_taxa_select = c("Plant")
pos = median_data_sep2 %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.97 * max(log(sigma)),
            damp.pos = 0.97 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.1,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.1)

summ <- saveTable_Class %>%
  filter(db_taxa %in% db_taxa_select)%>%
  group_by(grp) %>%
  filter(variable %in% "logdamping.time ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  Class = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
)

p_Tc_damping_class2 <- ggplot(median_data_sep2 %>%
                              filter(db_taxa %in% db_taxa_select),
                            aes(log(Tc), log(damping.time)))+
  geom_point(alpha = 1, shape = 21, size = 2.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(tau)))+
  facet_wrap( ~Class, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = 0.4, y = pos$damp.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = 0.4, y = pos$damp.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)
print(p_Tc_damping_class2)
ggsave(paste0("./plot/Tc and damping ",db_taxa_select," facet by Class PGLS.png"), p_Tc_damping_class2, width = 8, height = 4.5)
ggsave("./plot/Figure_C5.pdf", p_Tc_damping_class2 , width = 8, height = 4.5)

#### (2) sigma and damping #####
pos = median_data_sep %>% group_by(db_sep) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.9 * max(log(sigma)),
            damp.pos = 0.9 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.15,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.15)

summ <- saveTable %>%
  group_by(grp) %>%
  filter(variable %in% "logdamping.time ~ logsigma")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

# CI for slope
0.6026386-0.04642599
0.6026386+0.04642599

0.5245011-0.07077411
0.5245011+0.07077411

0.3795539-0.04529756
0.3795539+0.04529756

f_labels <- data.frame(
  db_sep = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_sigma_damping <- ggplot(median_data_sep,
                       aes(log(sigma), log(damping.time)))+
  geom_point(shape = 21, size = 2.5, alpha = 0.7)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(S)))+ylab(expression(log(tau)))+
  facet_wrap(. ~db_sep, nrow = 3, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -5.5, y = pos$damp.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -5.5, y = pos$damp.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)

p_sigma_damping
ggsave("./plot/sigma and damping PGLS.png", p_sigma_damping , width = 6, height = 8)
ggsave("./plot/Figure_C7.pdf", p_sigma_damping , width = 6, height = 8)

#### (3) Reproductive span versus generation time ####
# Use PGLS for age-structured animals
variable<- c("logreptime", "logTc")
data_group = c("Animal by age")
saveTable_ageanimal = PGLS_fun(filter(median_data_sep,grp == data_group), variable, data_group, outputlist)[[1]]%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)


pos = median_data_sep %>%
  group_by(grp) %>%
  filter(grp == data_group)%>%
  mutate(reptime = case_when(
    (omega - alpha) >0 ~ omega - alpha,
    TRUE ~ 0.000000001
  ))%>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            reptime.pos = 0.97 * max(log(reptime)),
            reptime.pos2 = reptime.pos - (max(log(reptime) - min(log(reptime))))*0.06)


summ <- saveTable_ageanimal %>%
  group_by(grp) %>%
  filter(variable %in% "logreptime ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  db_sep = summ$grp,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

p_Tc_rep <- ggplot(filter(median_data_sep,grp == data_group),
                       aes(logTc, logreptime))+
  geom_point(alpha = 1, shape = 21, size = 2.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = c(0.80, 0.25),
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(T["c"])))+ylab(expression(log(omega - alpha)))+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -3.2, y = pos$reptime.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -3.2, y = pos$reptime.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)
p_Tc_rep
ggsave("./plot/Tc and reproduction period age structured animal PGLS.png", p_Tc_rep , width = 6, height = 5)
ggsave("./plot/Figure_C3.pdf", p_Tc_rep , width = 6, height = 5)

#### (4) Comparison of tau ####
# No need to use PGLS
db_sep_list= c("Animal by age", "Animal by stage", "Plant by stage")
lm1 = lm(log(1/damping.approx) ~ log(1/damping.cal) , data = filter(all_data, db_sep %in%db_sep_list[1]))
lm2 = lm(log(1/damping.approx) ~ log(1/damping.cal) , data = filter(all_data, db_sep %in%db_sep_list[2]))
lm3 = lm(log(1/damping.approx) ~ log(1/damping.cal) , data = filter(all_data, db_sep %in%db_sep_list[3]))

pos = all_data %>%
  group_by(db_sep) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            dampapprox.pos = 0.93 * max(log(1/damping.approx)),
            dampapprox.pos2 = dampapprox.pos - (max(log(1/damping.approx) - min(log(1/damping.approx))))*0.13)


summ <- data.frame(
  db_sep = db_sep_list,
  slope = c(summary(lm1)$coefficients[2,1], summary(lm2)$coefficients[2,1], summary(lm3)$coefficients[2,1]),
  intercept = c(summary(lm1)$coefficients[1,1], summary(lm2)$coefficients[1,1], summary(lm3)$coefficients[1,1]),
  SE = c(summary(lm1)$coefficients[2,2], summary(lm2)$coefficients[2,2], summary(lm3)$coefficients[2,2]),
  R2 = c(summary(lm1)$r.squared, summary(lm2)$r.squared, summary(lm3)$r.squared),
  Padjusted = c(summary(lm1)$coefficients[2,4], summary(lm2)$coefficients[2,4], summary(lm3)$coefficients[2,4])
)%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)


f_labels <- data.frame(
  db_sep = db_sep_list,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)


p_damping <- ggplot(all_data,
                    aes(log(1/damping.cal), log(1/damping.approx))) +
  geom_point(shape = 21, alpha = 0.7, size = 2.5)+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
  #              parse = TRUE, color = "black", rr.digits = 2, coef.digits = 3,size = 5) +
  # scale_shape_manual(values=myshape)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15))+
  facet_wrap(. ~db_sep, nrow = 3, scales = "free_y") +
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(paste(log(tau), " calculated from PPM")))+
  ylab(expression(paste(log(tau), " from analytical approximation")))+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -6, y = pos$dampapprox.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -6, y = pos$dampapprox.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)

p_damping
ggsave("./plot/damping calcuated and approximated.png", p_damping , width = 6, height = 8)
ggsave("./plot/Figure_C10.pdf", p_damping , width = 6, height = 8)


#### (5) residual of PGLS (Tc~damping) and sigma #####
variable<- c("logdamping.time", "logTc")
median_data_sep = median_data_sep%>%
  mutate(grp = db_sep)
data_group = unique(median_data_sep$grp)
residual = PGLS_fun(median_data_sep, variable, data_group, outputlist)[[2]]%>%
  filter(variable %in% "logdamping.time ~ logTc")

db_sep_list= c("Animal by age", "Animal by stage", "Plant by stage")
lm1 = lm(residual ~ log(sigma) , data = filter(residual, db_sep %in%db_sep_list[1]))
lm2 = lm(residual ~ log(sigma) , data = filter(residual, db_sep %in%db_sep_list[2]))
lm3 = lm(residual ~ log(sigma) , data = filter(residual, db_sep %in%db_sep_list[3]))

pos = residual %>%
  group_by(db_sep) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            residual.pos = 0.93 * max(residual),
            residual.pos2 = residual.pos - (max(residual - min(residual)))*0.13)

summ <- data.frame(
  db_sep = db_sep_list,
  slope = c(summary(lm1)$coefficients[2,1], summary(lm2)$coefficients[2,1], summary(lm3)$coefficients[2,1]),
  intercept = c(summary(lm1)$coefficients[1,1], summary(lm2)$coefficients[1,1], summary(lm3)$coefficients[1,1]),
  SE = c(summary(lm1)$coefficients[2,2], summary(lm2)$coefficients[2,2], summary(lm3)$coefficients[2,2]),
  R2 = c(summary(lm1)$r.squared, summary(lm2)$r.squared, summary(lm3)$r.squared),
  Padjusted = c(summary(lm1)$coefficients[2,4], summary(lm2)$coefficients[2,4], summary(lm3)$coefficients[2,4])
)%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)

f_labels <- data.frame(
  db_sep = db_sep_list,
  label1 = c(paste0("y = ",summ$slope,"x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
)

p_residual <- ggplot(residual,
                     aes(log(sigma), residual))+
  geom_point(shape = 21, alpha = 1, size = 2.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(S)))+ylab("residuals")+
  facet_wrap(. ~db_sep, nrow = 3, scales = "free_y")+
  geom_richtext(data = f_labels, aes(label = label1), size = 5,
                x = -5.5, y = pos$residual.pos,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)+
  geom_richtext(data = f_labels, aes(label = label2), size = 5,
                x = -5.5, y = pos$residual.pos2,
                hjust = 0,
                color = "black",
                fill = NA,
                label.colour = NA)
p_residual
ggsave("./plot/residuals of Tc and damping and sigma PGLS.png", p_residual , width = 6, height = 8)
ggsave("./plot/Figure_C6.pdf", p_residual , width = 6, height = 8)

#### (7) phylo PCA #####
db_sep_list= c("Animal by age", "Animal by stage", "Plant by stage")
# change i to choose different dataset
i=1
i=2
i=3
pca.data <- median_data_sep%>%
  filter(db_sep %in%db_sep_list[i])%>%
  dplyr::select(logTc, logsigma, logdamping.time, SpeciesAccepted)
rownames(pca.data) = pca.data$SpeciesAccepted
pca.data = pca.data%>%
  dplyr::select(-"SpeciesAccepted")

median.data.pca<-phyl.pca(final_tree_read,pca.data)#method = "lambda"
# median.data.pca2 <- prcomp(pca.data, center =  TRUE, scale. = TRUE) #simple PCA

summary(median.data.pca)
print(median.data.pca)

######################################
## plot for pPCA
scale = 2
s <- data.frame(summary(median.data.pca)$importance)
## a function to convert number to percentage
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

name_position_x =c(-1.1, 1.4, -0.5)
name_position_y = matrix(c(-0, 0, 0.3, #Tc
                    -0.2, -0.2, -0.3, #sigma
                    0.1, 0.1, 0.2),#tau
                    byrow = T, nrow = 3)

data.plot = median_data_sep%>%
  filter(db_sep %in%db_sep_list[i])%>%
  mutate(PC1 = median.data.pca$S[,1],
         PC2 = median.data.pca$S[,2])
p <- ggplot(data.plot,
       aes(PC1, PC2, color = Class))+
  geom_point(shape = 20, size = 4, alpha= 0.6)+
  ## for Tc
  geom_segment(aes(x = 0, y = 0,
                  xend = median.data.pca$L[1,1]*scale, yend = median.data.pca$L[1,2]*scale),
              arrow = arrow(length = unit(0.3, "cm")), color = "red")+
  geom_text(
    label= expression(log(T["c"])),
      x=median.data.pca$L[1,1]*scale+name_position_x[i],
    y=median.data.pca$L[1,2]*scale+name_position_y[1,i],
    color = "red",
    nudge_y = 0.5,
    size = 7
  )+
  ## for sigma
  geom_segment(aes(x = 0, y = 0,
                   xend = median.data.pca$L[2,1]*scale, yend = median.data.pca$L[2,2]*scale),
               arrow = arrow(length = unit(0.3, "cm")), color = "red")+
  geom_text(
    label= expression(log(S)),
    x=median.data.pca$L[2,1]*scale+name_position_x[i],
    y=median.data.pca$L[2,2]*scale+name_position_y[2,i],
    color = "red",
    nudge_y = 0.5,
    size = 7
  )+
  ## for tau
  geom_segment(aes(x = 0, y = 0,
                   xend = median.data.pca$L[3,1]*scale, yend = median.data.pca$L[3,2]*scale),
               arrow = arrow(length = unit(0.3, "cm")), color = "red")+
  geom_text(
    label=expression(log(tau)),
    x=median.data.pca$L[3,1]*scale+name_position_x[i],
    y=median.data.pca$L[3,2]*scale+name_position_y[3,i],
    color = "red",
    nudge_y = 0.5,
    size = 7
  )+
  geom_hline(yintercept=0, linetype=2, col = 'black')+
  geom_vline(xintercept=0, linetype=2, col = 'black')+
  xlab(paste0("PPC1 (",percent(s$PC1[2]),")"))+
  ylab(paste0("PPC2 (",percent(s$PC2[2]),")"))+
  ## theme
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20),
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

p
ggsave(paste0("./plot/PPCA for ",db_sep_list[i],".png"), p , width = 8, height = 6)
figure_list <- c("Figure_3", "Figure_C8", "Figure_C9")
ggsave(paste0("./plot/",figure_list[i],".pdf"), p , width = 8, height = 6)

