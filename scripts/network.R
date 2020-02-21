###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(plyr);library(dplyr);library(tidyr)
library(tibble);library(stringr);library(readr)
library(ggplot2)
library(igraph)
library(intergraph)
library(sna)
library(visNetwork)

## import network #############################################################################################################

fn <- read.delim("./data/FN_SNA.txt",stringsAsFactors = F) %>%
  tibble::column_to_rownames(var="X") %>%
  replace(.,is.na(.),0)

fnm <- as.matrix(fn)

fnnet <- graph.adjacency(fnm,mode="undirected",weighted=TRUE,diag=FALSE)

## import attributes ##########################################################################################################

## network with attributes

attributes <- read.delim("./data/FN_SNA_attributes.txt",stringsAsFactors = F) 

V(fnnet)$gender <- as.character(attributes$gender[match(V(fnnet)$name,attributes$name)])
V(fnnet)$shape <- ifelse(V(fnnet)$gender=="F","square", "circle")
V(fnnet)$label <- str_replace_all(V(fnnet)$name, "[.]", " ")



## community detection ####################################################################################################

fnnet.com.lou <- cluster_louvain(fnnet)

fnnet.com.lou$names <- str_replace_all(fnnet.com.lou$names, "[.]", " ")

V(fnnet)$membership <- fnnet.com.lou$membership
V(fnnet)$color <- ifelse(V(fnnet)$membership=="1","#E69F00",
                         ifelse(V(fnnet)$membership=="2","#CC79A7",
                                ifelse(V(fnnet)$membership=="3","#56B4E9", "#009E73")))

fnnet.membership.df <- data.frame(Name = fnnet.com.lou$names,
                                  Membership = fnnet.com.lou$membership)


## modularity ############################################################################################################

modularity(fnnet,membership(fnnet.com.lou))

## network plot ##########################################################################################################

# set node size to highlight FN and relevant alters

V(fnnet)$size <- c(9,5,5,5,5,5,7,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,5,5,5,5,5,5,
                   5,7,5,5,5,5,5,5,5,5,5,7,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)


# select network layout
# place vertices on the plane using the force-directed layout algorithm by Fruchterman and Reingold

set.seed (5250)
fn.layout <- layout.fruchterman.reingold(fnnet)


# plot
 
par(mar=c(0,0,0,0)+.05)
plot.igraph(fnnet,vertex.label=V(fnnet)$label,layout=fn.layout, 
            vertex.size =V(fnnet)$size,vertex.color = V(fnnet)$color,vertex.shape = V(fnnet)$shape,
            vertex.label.family ="Helvetica",vertex.label.cex = 0.75,
            vertex.label.color="black",edge.color="gray",
            edge.width=E(fnnet)$weight/6,
            edge.curved=FALSE,frame=FALSE,
            ylim=c(-1.1,1.1),
            xlim=c(-1.1,1.1), asp = 0)
legend("bottomleft", legend=c("government figures (female)",
                              "government figures (male)",
                              "family and close friends (female)",
                              "family and close friends (male)",
                              "social reformers (female)",
                              "social reformers (male)",
                              "statisticians, scientists and academics (female)",
                              "statisticians, scientists and academics (male)"), 
       pch=c(22,21,22,21,22,21,22,21), 
       pt.bg =c("#56B4E9","#56B4E9","#E69F00","#E69F00",
                "#009E73","#009E73","#CC79A7","#CC79A7"),
       pt.cex=1.1,cex=0.7,bty="n")


## network interactive plot ############################################################################################


V(fnnet)$shape <- ifelse(V(fnnet)$gender=="F","square", "dot")
V(fnnet)$name <- V(fnnet)$label
V(fnnet)$size <- c(30,20,20,20,20,20,25,20,20,20,20,20,20,20,20,20,20,20,20,20,20,25,20,20,20,20,20,20,20,25,20,20,
                   20,20,20,20,20,20,20,25,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20)

legend.nodes <- data.frame (label = c("government figures\n(female)",
                                      "government figures\n(male)",
                                      "family and close friends\n(female)",
                                      "family and close friends\n(male)",
                                      "social reformers\n(female)",
                                      "social reformers\n(male)",
                                      "statisticians, scientists and academics\n(female)",
                                      "statisticians, scientists and academics\n(male)"),
                            shape = c("square","dot","square","dot","square","dot",
                                      "square","dot"),size = rep(8,8),
                            color.background = c("#56B4E9","#56B4E9","#E69F00","#E69F00",
                                                 "#009E73","#009E73","#CC79A7","#CC79A7"),
                            color.border = c("#56B4E9","#56B4E9","#E69F00","#E69F00",
                                             "#009E73","#009E73","#CC79A7","#CC79A7"))


fn.egonet <- visIgraph(fnnet) %>%
  visEdges(color="grey",width =0.3) %>%
  visLegend(addNodes = legend.nodes, useGroups = FALSE,width = 0.13, position = "right") %>%
  visIgraphLayout(type = "full",layout = "layout_with_fr",randomSeed = 55055) 


visSave(fn.egonet, file = "fn.egonet.html")


## brokerage rules ##########################################################################################################

detach("package:igraph") 
fnnet2 <- asNetwork(fnnet)

# raw and normalized scores rounded to 2 digits

brokerage.raw <- as.data.frame(round(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "membership"))$raw.nli, 2)) %>%
  rownames_to_column("Alter") %>%
  rename(Coordinator = w_I, Itinerant = w_O, Representative = b_IO, Gatekeeper = b_OI, Liaison = b_O) %>%
  select(-t) %>%
  gather("Brokerage","Value",-Alter) %>%
  left_join(fnnet.membership.df,by=c("Alter" ="Name")) %>%
  mutate(Membership = as.factor(Membership)) %>%
  mutate(Alter = str_replace_all(Alter, "[.]", " "))

  
brokerage.norm <- as.data.frame(round(brokerage(fnnet2, cl=get.vertex.attribute(fnnet2, "membership"))$z.nli,2)) %>%
  rownames_to_column("Alter") %>%
  rename(Coordinator = w_I, Itinerant = w_O, Representative = b_IO, Gatekeeper = b_OI, Liaison = b_O) %>%
  select(-t) %>%
  gather("Brokerage","Value",-Alter) %>%
  left_join(fnnet.membership.df,by=c("Alter" ="Name")) %>%
  mutate(Membership = as.factor(Membership))

## brokerage plot ç##########################################################################################################

values.raw.1 <- brokerage.raw %>% group_by(Brokerage) %>% 
  filter(Alter %in% c("Parthenope Nightingale","Harriet Martineau","Lord Palmerston")) %>%
  filter(Brokerage!="Itinerant") 

values.raw.2<- brokerage.raw %>% group_by(Brokerage) %>% 
  filter(Alter=="Dr William Farr" & Brokerage=="Liaison" |
           Alter=="Ada Lovelace" & Brokerage %in% c("Representative","Gatekeeper","Coordinator"))  
           
values.raw.3<- brokerage.raw %>% group_by(Brokerage) %>% 
  filter(Alter=="Dr William Farr" &  Brokerage %in% c("Representative","Gatekeeper","Coordinator") |
           Alter=="Ada Lovelace" &  Brokerage=="Liaison")  


## fixed scale, FN value not captured in plot

ggplot(brokerage.raw %>% group_by(Brokerage) %>% 
         filter(Alter!="Florence Nightingale") %>%
         filter(Brokerage!="Itinerant"),aes(x=Membership, y = Value,fill=Membership)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  facet_wrap(~Brokerage,scales="fixed",ncol =2) +
  scale_fill_manual(labels = c("family and close friends", 
                               "statisticians, scientists and academics",
                               "government figures",
                               "social reformers"),
                    values=c("#E69F00", "#CC79A7", "#56B4E9", "#009E73")) +
  ylab("raw brokerage scores") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.position="bottom") +
  theme(strip.background = element_blank())+
  theme(panel.border = element_rect(colour = "black")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_text(data = values.raw.1, aes(x = Membership, y = Value + 50, label = Alter), size = 3) +
  geom_text(data = values.raw.2, aes(x = Membership, y = Value + 50, label = Alter), size = 3) +
  geom_text(data = values.raw.3, aes(x = Membership, y = Value + 40, label = Alter), size = 3) +
  geom_point(data = values.raw.1, aes(x = Membership, y = Value), shape=15,size = 2) +
  geom_point(data = values.raw.2, aes(x = Membership, y = Value), shape=15,size = 2) +
  geom_point(data = values.raw.3, aes(x = Membership, y = Value), shape=15,size = 2) 

##############################################################################################################################
###############################################################################################################################

