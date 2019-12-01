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

V(fnnet)$color <- fnnet_com_lou$membership
V(fnnet)$color <- ifelse(V(fnnet)$color=="1","#E69F00",
                         ifelse(V(fnnet)$color=="2","#CC79A7",
                                ifelse(V(fnnet)$color=="3","#56B4E9", "#009E73")))


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
brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "gender"))$raw.nli

# normalized scores
brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "gender"))$z.nli   
# Normalized, rounded to 2 digits
round(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "gender"))$z.nli, 2) 

brokerage.df <- as.data.frame(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "gender"))$raw.nli)
  
brokerage.df.barplot  <- brokerage.df %>%
  rownames_to_column("Alter") %>%
  arrange(-t) %>%
  filter(t<1804 & t>50) %>%
  rename(Coordinator = w_I, Itinerant = w_O, Representative = b_IO, Gatekeeper = b_OI, Liaison = b_O) %>%
  select(-t) %>%
  select(-Liaison) %>%
  gather("Brokerage","Value",-Alter) %>%
  filter(Alter %in% c("Richard.Milnes","Parthenope.Nightingale","Arthur.Hugh.Clough","Aunt.Mai.Smith","Selina.Bracebridge"))

## brokerage plot ##########################################################################################################

ggplot(brokerage.df.barplot,aes(x=Brokerage, y = Value,fill=Brokerage)) + 
  geom_bar(stat="identity") +
  facet_wrap(~Alter,nrow=1) +
  scale_fill_manual(values=c("#999999", "#F0E442", "#0072B2", "#009E73")) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.position="none") +
  theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
