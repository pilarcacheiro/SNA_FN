###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(plyr);library(dplyr);library(tidyr);library(stringr);library(readr);library(ggplot2);
library(igraph);library(intergraph)
library(sna)

## set dir ####################################################################################################################

fn.dir <-  ".../SNA_FN/Data/"

###############################################################################################################################

## import network #############################################################################################################

fn <- read.delim(paste0(fn.dir,"FN_SNA.txt"),stringsAsFactors = F) %>%
  tibble::column_to_rownames(var="X") %>%
  replace(.,is.na(.),0)

fnm <- as.matrix(fn)

fnnet <- graph.adjacency(fnm,mode="undirected",weighted=TRUE,diag=FALSE)

## import attributes ##########################################################################################################


## network with attributes
## layout_with_dh

V(fnnet)$name

attributes <- read.delim(paste0(fn.dir,"FN_SNA_attributes.txt"),stringsAsFactors = F) 

V(fnnet)$gender <- as.character(attributes$gender[match(V(fnnet)$name,attributes$name)])
V(fnnet)$shape <- ifelse(V(fnnet)$gender=="F","square", "circle")



set.seed(1000)
coords = layout_with_dh(fnnet)



## community detection ####################################################################################################

fnnet_com_lou <- cluster_louvain(fnnet)

fnnet_com_lou$membership
fnnet_com_lou$names

V(fnnet)$membership <- fnnet_com_lou$membership
V(fnnet)$color <- ifelse(V(fnnet)$membership=="1","#E69F00",
                         ifelse(V(fnnet)$membership=="2","#CC79A7",
                                ifelse(V(fnnet)$membership=="3","#56B4E9", "#009E73")))


## network plot ##########################################################################################################
  
par(mar=c(0,0,0,0)+.05)
plot.igraph(fnnet,vertex.label=V(fnnet)$Name,layout=coords, 
            vertex.size =5,vertex.color = V(fnnet)$color,vertex.shape = V(fnnet)$shape,
            vertex.label.family ="Helvetica",vertex.label.cex = 0.6,
            vertex.label.color="black",edge.color="gray72",edge.width=E(fnnet)$weight/4,
            edge.curved=FALSE,frame=FALSE,
            ylim=c(-1.1,1.1),
            xlim=c(-1.1,1.1), asp = 0)
legend("bottomleft", legend=c("Female", "Male"), 
       pch=c(0,1),  pt.cex=2,cex=1, bty="n")


## brokerage rules ##########################################################################################################

detach("package:igraph") 
fnnet2 <- asNetwork(fnnet)

# raw and normalized scores rounded to 2 digits

brokerage.raw <- as.data.frame(round(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "membership"))$raw.nli, 2)) %>%
  rownames_to_column("Alter") %>%
  rename(Coordinator = w_I, Itinerant = w_O, Representative = b_IO, Gatekeeper = b_OI, Liaison = b_O) %>%
  select(-t) %>%
  gather("Brokerage","Value",-Alter) %>%
  left_join(fnnet_membership.df,by=c("Alter" ="Name")) %>%
  mutate(Membership = as.factor(Membership))

  
brokerage.norm <- as.data.frame(round(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "membership"))$z.nli,2)) %>%
  rownames_to_column("Alter") %>%
  rename(Coordinator = w_I, Itinerant = w_O, Representative = b_IO, Gatekeeper = b_OI, Liaison = b_O) %>%
  select(-t) %>%
  gather("Brokerage","Value",-Alter) %>%
  left_join(fnnet_membership.df,by=c("Alter" ="Name")) %>%
  mutate(Membership = as.factor(Membership))

## brokerage plots ##########################################################################################################

values.raw <- brokerage.raw %>% group_by(Brokerage) %>% 
  filter(Alter %in% c("Ada.Lovelace","Parthenope.Nightingale","Harriet.Martineau","Lord.Palmerston")) %>%
  filter(Brokerage!="Itinerant")


ggplot(brokerage.raw %>% group_by(Brokerage) %>% 
         filter(!Value %in% boxplot(brokerage.raw$Value)$out) %>%
         filter(Brokerage!="Itinerant"),aes(x=Membership, y = Value,fill=Membership)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Brokerage,scales="free_y",ncol =2) +
  scale_fill_manual(values=c("#E69F00", "#CC79A7", "#56B4E9", "#009E73")) +
  ylab("raw brokerage scores") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.position="none") +
  theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=10,margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  geom_text(data = values.raw, aes(x = Membership, y = Value+5, label = Alter), size = 3) +
  geom_point(data = values.raw, aes(x = Membership, y = Value))


values.norm <- brokerage.norm %>% group_by(Brokerage) %>% 
  filter(Alter %in% c("Ada.Lovelace","Parthenope.Nightingale","Harriet.Martineau","Lord.Palmerston")) %>%
  filter(Brokerage!="Itinerant")

ggplot(brokerage.norm %>% group_by(Brokerage) %>% 
         filter(!Value %in% boxplot(brokerage.norm$Value)$out) %>%
         filter(Brokerage!="Itinerant"),aes(x=Membership, y = Value,fill=Membership)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Brokerage,scales="free_y",ncol =2) +
  scale_fill_manual(values=c("#E69F00", "#CC79A7", "#56B4E9", "#009E73")) +
  ylab("normalized brokerage scores") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.position="none") +
  theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_text(data = values.norm, aes(x = Membership, y = Value +1, label = Alter), size = 3) +
  geom_point(data = values.norm, aes(x = Membership, y = Value))

###############################################################################################################################
###############################################################################################################################

