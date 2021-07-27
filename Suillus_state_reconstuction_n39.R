#state reconstuction for Nhu's Paper
#started 26.Jul.2021
#version #Running R version 4.0.1

#set directory
#setwd("")

#make reproducible
set.seed(666)

#load libs
library("phytools")
library("ape")
library("quadprog")
library("devtools")
library("phangorn")
library("ggtree")

state_df_iq<- read.csv("Suillus_trait_states_IQ.csv")

#read in tree
phy_iq<- as.character(state_df_iq[1,4])
phy_iq<- read.tree(text = phy_iq)

#read in the state data
host_state_df <- data.frame(taxa=state_df_iq$X, host_state=(state_df_iq[,3]))

x2<- host_state_df[,2] 
x3<-setNames(host_state_df[,2],host_state_df[,1])


##re-root
#OG tree
basic_tree <- 
  ggtree(tr = phy_iq, 
         layout  = 'rectangular')
plot(basic_tree + geom_tiplab(size = 1, align = TRUE, linesize = .25, linetype = 3))

#root w/ outgroup
#phy_iq$tip.label #8
phy_iq_rooted<- phytools::reroot(phy_iq, 8, interactive=FALSE, position =.01)

basic_tree2 <- 
  ggtree(tr = phy_iq_rooted, 
         layout  = 'rectangular')
basic_tree2 + geom_tiplab(size = 1, align = TRUE, linesize = .25, linetype = 3)

#make a simmap to connect chr states to tips
phy_iq<-make.simmap(phy_iq_rooted, x3)

#link the host states
states_iq<-getStates(phy_iq,"tips")

#set color pallet: 
#(#012623 = dark green = )
#(#035941 = light green = )
#(#AFBF36 = lime green = )
#(#F2622E = orange = )
#(#F20505 = red =)
#palette = c("#012623", "#035941", "#AFBF36", "#F2622E", "#F20505")

#plot IQ tree
#fitER<-ace(states_iq,phy_iq,model="ER",type="discrete")
#add neglagble length to root ro avoid error
dst<-phy_iq
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(phy_iq))*1e-6

##choose model
# equal rates model 
fitER<-fitMk(dst,states_iq,model="ER")
fitER
plot(fitER)

# SYM
fitSYM<-fitMk(dst,states_iq,model="SYM")
fitSYM
plot(fitSYM, show.zeros=FALSE)

# ARD 
fitARD<-fitMk(dst,states_iq,model="ARD")
fitARD
plot(fitARD,show.zeros=FALSE)

#get AIC
AIC_iq<-setNames(sapply(list(fitER,fitSYM,fitARD),AIC),c("ER","SYM","ARD"))
AIC_iq
aic.w(AIC_iq)
#ER wins
#ER        SYM        ARD 
#0.8396497 0.1598033 0.0005470 

#fit with ER model
fitER<-ace(states_iq,dst,type="discrete",model="ER")
fitER

#plot tree
#cols for pies
cols<-setNames(palette(c("#F2E6D8", "#AD8286", "#64273B", "#3E4C63", "#82ABBA"))[1:length(unique(states_iq))],sort(unique(states_iq)))
#cols for branches
branch_col<- setNames(palette(c("#231F20", "#231F20"))[1:length(unique(states_iq))],sort(unique(states_iq)))


#note - there's a bug in phytools, where it doesnt like to set colors- 
#if the colors show up in the pies, run the color settign function again until it looks right. 

#ladderize
phy_iq_l<- ladderize(phy_iq, right=FALSE)

#plot basic tree
#plotTree(phy_iq_l,type="phylogram",fsize=0.7, ftype="i", colors = branch_col)
#nodelabels(node=1:phy_iq$Nnode+Ntip(phy_iq),
#           pie=fitER$lik.anc,piecol=cols,cex=0.5)
#tiplabels(pie=to.matrix(states_iq,sort(unique(states_iq))),piecol=cols,cex=0.2)
#add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
#                  y=-max(nodeHeights(phy_iq)),fsize=0.8)


##replace JGI names with full species names w/ project codes
#first, sort the input replacement file to match the order of the tip labels 
state_df_iq_sorted<- state_df_iq[match(phy_iq_l$tip.label, state_df_iq$X),]

#now replace tip labels 
phy_iq_l$tip.label <- state_df_iq_sorted$X.1

#replot that shit
###plot tree with tips aligned
par(fg="transparent")
cw<-reorder(phy_iq_l,"cladewise")
plotTree(cw, fsize=0.4, lwd = .9)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
par(fg="black")
text(rep(max(obj$xx[1:Ntip(cw)]),Ntip(cw)),obj$yy[1:Ntip(cw)],
     labels=cw$tip.label,font=3,pos=4, cex=0.5)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])),
                           rep(obj$yy[i],2),lty="dotted")
#add pies
nodelabels(node=1:phy_iq$Nnode+Ntip(phy_iq),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
#add ovserved states
tiplabels(pie=to.matrix(states_iq,sort(unique(states_iq))),piecol=cols,cex=0.2)
#add branch support values
nodelabels(phy_iq$node.label,adj=c(1.5,-0.8),frame="none", cex=0.3)
#add trait (host)legend
add.simmap.legend(colors=cols,fsize=0.8, prompt=TRUE)#place legend with mouse before saving plot
#save output
dev.copy(pdf,"Suillus_tree_w_reconstruction.pdf")
dev.off()

