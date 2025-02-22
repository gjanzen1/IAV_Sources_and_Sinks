
################
### This script generates statistics and plots related to predictions of interstate spread of IAV.
### Run script1 and script3 before running this script.
################

setwd("C:/Users/garrett.janzen/OneDrive - USDA/Projects/IAV_Env_Eco")
load("IAV_Sources_and_Sinks.RData") # required to run this script

##########################

library("ggplot2")
library("gplots")
library("reshape2")
library("markovchain")
library("diagram")
library("lattice")
library("ade4")
library("vegan")
library("igraph")

##########################

#### Choose your preferred version of constellation: #### 
# data_sub <- data_ss[,c("Date","State","Constellation"),]
# data_sub <- data_ss[,c("Date","State","UID_simple"),]
data_sub <- data_ss[,c("Date","State","UID_complex"),]

#######################################################
############# Considering all states ###################
#### !!! You must also modify the 3rd line of the for-loop below, choosing `constellations`, `UIDs_simple`, or `UIDs_complex`,
#### !!! and must modify the line below, choosing `$Constellation`, `$UID_simple`, or `$UID_complex`

statestringlist <- list()
statestringlist_repeat <- list()
statestring_repeat <- i <- r <- k <- NULL

for(i in 1:length(UIDs_complex)){
  const <- UIDs_complex[i]
  dat <- data_sub[which(data_sub$UID_complex == const),];dim(dat)
  statestring_repeat <- dat$State
  statestringlist_repeat <- append(statestringlist_repeat, list(statestring_repeat))
  names(statestringlist_repeat)[i] <- const
  statestring_repeat <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring <- dat$State
  statestringlist <- append(statestringlist, list(statestring))
  names(statestringlist)[i] <- const
  statestring <- NULL
}

statestringlist_xx <- statestringlist
statestringlist_repeat_xx <- statestringlist_repeat
data_dmax <- as.Date(max(data$Date))

for(i in 1:length(statestringlist_xx)){
  const <- names(statestringlist_xx)[i]
  
  dat <- data_sub[which(data_sub$UID_complex == const),];dim(dat)
  # Below, we append "XX" to the state string if the constellation's final detection is before the cut-off date, with uncertainty window
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_repeat_xx[[i]], "XX")), list(c(statestringlist_repeat_xx[[i]]))) ### GGG < or > ?
  statestringlist_repeat_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  xxstring <- ifelse(max(dat$Date) < data_dmax-window, list(c(statestringlist_xx[[i]], "XX")), list(c(statestringlist_xx[[i]]))) ### GGG < or > ?
  statestringlist_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
}

m <- do.call("rbind", lapply(statestringlist, function(x) cbind(head(x, -1), tail(x, -1))))
mc <- markovchainFit(m)
est <- mc$estimate
tm <- est@transitionMatrix

m_repeat <- do.call("rbind", lapply(statestringlist_repeat, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat <- markovchainFit(m_repeat)
est_repeat <- mc_repeat$estimate
tm_repeat <- est_repeat@transitionMatrix

m_xx <- do.call("rbind", lapply(statestringlist_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_xx <- markovchainFit(m_xx)
est_xx <- mc_xx$estimate
tm_xx <- est_xx@transitionMatrix

m_repeat_xx <- do.call("rbind", lapply(statestringlist_repeat_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_xx <- markovchainFit(m_repeat_xx)
est_repeat_xx <- mc_repeat_xx$estimate
tm_repeat_xx <- est_repeat_xx@transitionMatrix

tm_melt <- melt(tm)
colnames(tm_melt) <- c("Origin", "Destination", "Probability")
tm_melt_repeat <- melt(tm_repeat)
colnames(tm_melt_repeat) <- c("Origin", "Destination", "Probability")

hmplot <- ggplot(tm_melt, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot
hmplot_repeat <- ggplot(tm_melt_repeat, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat
# 
# pdf("Plots/script4_network.pdf")
# plot(est)
# dev.off()
# pdf("Plots/script4_tm.pdf", width=9.5)
# hmplot
# dev.off()
# png("Plots/script4_network.png")
# plot(est)
# dev.off()
# png("Plots/script4_tm.png", width=750)
# hmplot
# dev.off()
# 
# pdf("Plots/script4_network_repeat.pdf")
# plot(est_repeat)
# dev.off()
# pdf("Plots/script4_tm_repeat.pdf", width=9.5)
# hmplot_repeat
# dev.off()
# png("Plots/script4_network_repeat.png")
# plot(est_repeat)
# dev.off()
# png("Plots/script4_tm_repeat.png", width=750)
# hmplot_repeat
# dev.off()

tm_thresh <- ifelse(tm < 0.25, 0, tm)
tm_thresh_repeat <- ifelse(tm_repeat < 0.25, 0, tm_repeat)

# dev.off()
# pdf("Plots/script4_visnetwork.pdf")
# vis <- plotmat(round(t(tm), digits=2), relsize=1.1,
#                box.type="circle", box.size=0.02,
#                shadow.col="grey", shadow.size=0.005);vis
# dev.off()
# pdf("Plots/script4_visnetwork_thresh.pdf")
# vis_thresh <- plotmat(round(t(tm_thresh), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_thresh
# dev.off()
# dev.off()
# png("Plots/script4_visnetwork.png")
# vis <- plotmat(round(t(tm), digits=2), relsize=1.1,
#                box.type="circle", box.size=0.02,
#                shadow.col="grey", shadow.size=0.005);vis
# dev.off()
# png("Plots/script4_visnetwork_thresh.png")
# vis_thresh <- plotmat(round(t(tm_thresh), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_thresh
# dev.off()
# 
# dev.off()
# pdf("Plots/script4_visnetwork_repeat.pdf")
# vis_repeat <- plotmat(round(t(tm_repeat), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_repeat
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_repeat.pdf")
# vis_thresh_repeat <- plotmat(round(t(tm_thresh_repeat), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_thresh_repeat
# dev.off()
# 
# dev.off()
# png("Plots/script4_visnetwork_repeat.png")
# vis_repeat <- plotmat(round(t(tm_repeat), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_repeat
# dev.off()
# png("Plots/script4_visnetwork_thresh_repeat.png")
# vis_thresh_repeat <- plotmat(round(t(tm_thresh_repeat), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_thresh_repeat
# dev.off()

#########
#Using tm, but dropping rows/colummns of states we don't want to see.
#This is NOT the same as the _subset runs, which drop states before transitions are counted.
#This version ONLY effects the graph, not the data itself.
mcstatereduced <- c("IL","NC","MN","MO","IN","NE","IA","OH","OK","SD","XX") #top 10 states

tm_reduced <- tm
tm_reduced <- tm_reduced[which(colnames(tm_reduced) %in% mcstatereduced), which(rownames(tm_reduced) %in% mcstatereduced)]

tm_reduced_thresh <- tm_reduced
tm_reduced_thresh <- ifelse(tm_reduced <= 0.14, 0, tm_reduced) 

# pdf("Plots/script4_visnetwork_thresh_reduced.pdf")
# vis_thresh <- plotmat(round(t(tm_reduced_thresh), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005, cex.txt = 1.3, dtext = 0.6,
#                       cex.main = 1.2, txt.col="black",box.col="beige", lcol="black");vis_thresh
# dev.off()
# png("Plots/script4_visnetwork_thresh_reduced.png", width=1100, height=750)
# vis_thresh <- plotmat(round(t(tm_reduced_thresh), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005, cex.txt = 1.3, dtext = 0.6,
#                       cex.main = 1.2, txt.col="black",box.col="beige", lcol="black");vis_thresh
# dev.off()
#########

tm_melt_xx <- melt(tm_xx)
colnames(tm_melt_xx) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_xx <- melt(tm_repeat_xx)
colnames(tm_melt_repeat_xx) <- c("Origin", "Destination", "Probability")
#We want to remove any rows that have an origin of XX, for the purpose of plots
tm_melt_xx_noxxorigin <- tm_melt_xx[which(tm_melt_xx$Origin != "XX"),];dim(tm_melt_xx);dim(tm_melt_xx_noxxorigin)
tm_melt_repeat_xx_noxxorigin <- tm_melt_repeat_xx[which(tm_melt_repeat_xx$Origin != "XX"),];dim(tm_melt_repeat_xx);dim(tm_melt_repeat_xx_noxxorigin)

hmplot_xx <- ggplot(tm_melt_xx_noxxorigin, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_xx
hmplot_repeat_xx <- ggplot(tm_melt_repeat_xx_noxxorigin, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix");hmplot_repeat_xx

# pdf("Plots/script4_network_xx.pdf")
# plot(est_xx)
# dev.off()
# pdf("Plots/script4_tm_xx.pdf", width=9.5)
# hmplot_xx
# dev.off()
# png("Plots/script4_network_xx.png")
# plot(est_xx)
# dev.off()
# png("Plots/script4_tm_xx.png", width=750)
# hmplot_xx
# dev.off()

# pdf("Plots/script4_network_repeat_xx.pdf")
# plot(est_repeat_xx)
# dev.off()
pdf("Plots/script4_tm_repeat_xx.pdf", width=9.5)
hmplot_repeat_xx
dev.off()
# png("Plots/script4_network_repeat_xx.png")
# plot(est_repeat_xx)
# dev.off()
png("Plots/script4_tm_repeat_xx.png", width=750)
hmplot_repeat_xx
dev.off()

tm_thresh_xx <- ifelse(tm_xx < 0.25, 0, tm_xx)
tm_thresh_repeat_xx <- ifelse(tm_repeat_xx < 0.25, 0, tm_repeat_xx)

# dev.off()
# pdf("Plots/script4_visnetwork_xx.pdf")
# vis_xx <- plotmat(round(t(tm_xx), digits=2), relsize=1.1,
#                   box.type="circle", box.size=0.02,
#                   shadow.col="grey", shadow.size=0.005);vis_xx
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_xx.pdf")
# vis_thresh_xx <- plotmat(round(t(tm_thresh_xx), digits=2), relsize=1.1,
#                          box.type="circle", box.size=0.02,
#                          shadow.col="grey", shadow.size=0.005);vis_thresh_xx
# dev.off()
# 
# dev.off()
# png("Plots/script4_visnetwork_xx.png")
# vis_xx <- plotmat(round(t(tm_xx), digits=2), relsize=1.1,
#                   box.type="circle", box.size=0.02,
#                   shadow.col="grey", shadow.size=0.005);vis_xx
# dev.off()
# png("Plots/script4_visnetwork_thresh_xx.png")
# vis_thresh_xx <- plotmat(round(t(tm_thresh_xx), digits=2), relsize=1.1,
#                          box.type="circle", box.size=0.02,
#                          shadow.col="grey", shadow.size=0.005);vis_thresh_xx
# dev.off()

##################################################################
############# Considering only specific states ###################
#### !!! You must also modify the 3rd line of the for-loop below, choosing `constellations`, `UIDs_simple`, or `UIDs_complex`,
#### !!! and must modify the line below, choosing `$Constellation`, `$UID_simple`, or `$UID_complex`

statestringlist_subset <- list()
statestringlist_repeat_subset <- list()
statestringlist_subset_xx <- list()
statestringlist_repeat_subset_xx <- list()

statestring_repeat_subset <- i <- r <- k <- NULL

#select your desired states by modifying this string:
# mcstatesubset <- c("(IL|NC|MN|MO|IN|NE|IA|XX)")
mcstatesubset <- c("(IA|MN|NC|IL|IN|NE|MO|OH|OK|SD|XX)")

for(i in 1:length(UIDs_complex)){
  # i <- 14
  const <- UIDs_complex[i]
  dat <- data_sub[which(data_sub$UID_complex == const),];dim(dat)
  statestring_repeat_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_repeat_subset <- statestring_repeat_subset[grepl(paste0(mcstatesubset), statestring_repeat_subset)]
  statestringlist_repeat_subset <- append(statestringlist_repeat_subset, list(statestring_repeat_subset))
  names(statestringlist_repeat_subset)[i] <- const
  statestring_repeat_subset <- NULL
  
  dat <- dat[which(!duplicated(dat$State)),];dim(dat)
  statestring_subset <- dat$State
  ### To focus on a subset of states, filter here:
  statestring_subset <- statestring_subset[grepl(paste0(mcstatesubset), statestring_subset)]
  statestringlist_subset <- append(statestringlist_subset, list(statestring_subset))
  names(statestringlist_subset)[i] <- const
  statestring_subset <- NULL
  
  statestringlist_repeat_subset_xx[[i]] <- statestringlist_repeat_subset[[i]]
  xxstring <- ifelse("XX" %in% statestringlist_repeat_xx[[i]], list(c(statestringlist_repeat_subset_xx[[i]], "XX")), list(c(statestringlist_repeat_subset_xx[[i]])))
  statestringlist_repeat_subset_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
  
  statestringlist_subset_xx[[i]] <- statestringlist_subset[[i]]
  xxstring <- ifelse("XX" %in% statestringlist_xx[[i]], list(c(statestringlist_subset_xx[[i]], "XX")), list(c(statestringlist_subset_xx[[i]])))
  statestringlist_subset_xx[[i]] <- unlist(xxstring)
  xxstring <- NULL
}

m_subset <- do.call("rbind", lapply(statestringlist_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_subset <- markovchainFit(m_subset)
est_subset <- mc_subset$estimate
tm_subset <- est_subset@transitionMatrix

m_repeat_subset <- do.call("rbind", lapply(statestringlist_repeat_subset, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_subset <- markovchainFit(m_repeat_subset)
est_repeat_subset <- mc_repeat_subset$estimate
tm_repeat_subset <- est_repeat_subset@transitionMatrix

m_subset_xx <- do.call("rbind", lapply(statestringlist_subset_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_subset_xx <- markovchainFit(m_subset_xx)
est_subset_xx <- mc_subset_xx$estimate
tm_subset_xx <- est_subset_xx@transitionMatrix

m_repeat_subset_xx <- do.call("rbind", lapply(statestringlist_repeat_subset_xx, function(x) cbind(head(x, -1), tail(x, -1))))
mc_repeat_subset_xx <- markovchainFit(m_repeat_subset_xx)
est_repeat_subset_xx <- mc_repeat_subset_xx$estimate
tm_repeat_subset_xx <- est_repeat_subset_xx@transitionMatrix

tm_melt_subset <- melt(tm_subset)
colnames(tm_melt_subset) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_subset <- melt(tm_repeat_subset)
colnames(tm_melt_repeat_subset) <- c("Origin", "Destination", "Probability")

hmplot_subset <- ggplot(tm_melt_subset, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_tile() +
  labs(title="Markov chain transition matrix") +
  theme(axis.text.y = element_text(size = 22),
        text = element_text(size=26));hmplot_subset
hmplot_repeat_subset <- ggplot(tm_melt_repeat_subset, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix") +
  theme(axis.text.y = element_text(size = 22),
        text = element_text(size=26));hmplot_repeat_subset

# pdf("Plots/script4_network_subset.pdf")
# plot(est_subset)
# dev.off()
# pdf("Plots/script4_tm_subset.pdf", width=9.5)
# hmplot_subset
# dev.off()

# pdf("Plots/script4_network_repeat_subset.pdf")
# plot(est_repeat_subset)
# dev.off()
# pdf("Plots/script4_tm_repeat_subset.pdf", width=9.5)
# hmplot_repeat_subset
# dev.off()

# png("Plots/script4_network_subset.png", width=750, height=750)
# plot(est_subset)
# dev.off()
# png("Plots/script4_tm_subset.png", width=750, height=650)
# hmplot_subset
# dev.off()

# png("Plots/script4_network_repeat_subset.png", width=750, height=750)
# plot(est_repeat_subset)
# dev.off()
# png("Plots/script4_tm_repeat_subset.png", width=750, height=650)
# hmplot_repeat_subset
# dev.off()

###

tm_melt_subset_xx <- melt(tm_subset_xx)
colnames(tm_melt_subset_xx) <- c("Origin", "Destination", "Probability")
tm_melt_repeat_subset_xx <- melt(tm_repeat_subset_xx)
colnames(tm_melt_repeat_subset_xx) <- c("Origin", "Destination", "Probability")

hmplot_subset_xx <- ggplot(tm_melt_subset_xx, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_tile() +
  labs(title="Markov chain transition matrix") +
  theme(axis.text.y = element_text(size = 22),
        text = element_text(size=26));hmplot_subset_xx
hmplot_repeat_subset_xx <- ggplot(tm_melt_repeat_subset_xx, aes(Destination, Origin, fill = Probability)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low="#CCFFFF", high="dodgerblue4") +
  geom_tile() +
  labs(title="Markov chain transition matrix") +
  theme(axis.text.y = element_text(size = 22),
        text = element_text(size=26));hmplot_repeat_subset_xx

# pdf("Plots/script4_network_subset_xx.pdf")
# plot(est_subset_xx)
# dev.off()
# pdf("Plots/script4_tm_subset_xx.pdf", width=9.5)
# hmplot_subset_xx
# dev.off()

# pdf("Plots/script4_network_repeat_subset_xx.pdf")
# plot(est_repeat_subset_xx)
# dev.off()
# pdf("Plots/script4_tm_repeat_subset_xx.pdf", width=9.5)
# hmplot_repeat_subset_xx
# dev.off()

# png("Plots/script4_network_subset_xx.png", width=750, height=750)
# plot(est_subset_xx)
# dev.off()
# png("Plots/script4_tm_subset_xx.png", width=750, height=650)
# hmplot_subset_xx
# dev.off()

# png("Plots/script4_network_repeat_subset_xx.png", width=750, height=750)
# plot(est_repeat_subset_xx)
# dev.off()
# png("Plots/script4_tm_repeat_subset_xx.png", width=750, height=650)
# hmplot_repeat_subset_xx
# dev.off()

###You would not use this reduced model for these sorts of predictions.
# #Predict a new constellation
# newconst <- c("MS")
# prior_state <- tail(newconst, 1)
# predict(est_subset, prior_state) #transmission to a new state
# predict(est_repeat_subset, prior_state) #next detection location

tm_thresh_subset <- ifelse(tm_subset <= 0.18, 0, tm_subset) #mean(tm_subset) + sd(tm_subset)
tm_thresh_repeat_subset <- ifelse(tm_repeat_subset <= 0.18, 0, tm_repeat_subset) #mean(tm_subset) + sd(tm_subset)

# dev.off()
# pdf("Plots/script4_visnetwork_subset.pdf")
# vis_subset <- plotmat(round(t(tm_subset), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_subset
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_subset.pdf")
# vis_thresh_subset <- plotmat(round(t(tm_thresh_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_thresh_subset
# dev.off()
# dev.off()
# png("Plots/script4_visnetwork_subset.png")
# vis_subset <- plotmat(round(t(tm_subset), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_subset
# dev.off()
# png("Plots/script4_visnetwork_thresh_subset.png")
# vis_thresh_subset <- plotmat(round(t(tm_thresh_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_thresh_subset
# dev.off()
# 
# dev.off()
# pdf("Plots/script4_visnetwork_repeat_subset.pdf")
# vis_repeat_subset <- plotmat(round(t(tm_repeat_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_repeat_subset
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_repeat_subset.pdf")
# vis_thresh_repeat_subset <- plotmat(round(t(tm_thresh_repeat_subset), digits=2), relsize=1.1,
#                                     box.type="circle", box.size=0.02,
#                                     shadow.col="grey", shadow.size=0.005);vis_thresh_repeat_subset
# dev.off()
# dev.off()
# png("Plots/script4_visnetwork_repeat_subset.png")
# vis_repeat_subset <- plotmat(round(t(tm_repeat_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_repeat_subset
# dev.off()
# png("Plots/script4_visnetwork_thresh_repeat_subset.png")
# vis_thresh_repeat_subset <- plotmat(round(t(tm_thresh_repeat_subset), digits=2), relsize=1.1,
#                                     box.type="circle", box.size=0.02,
#                                     shadow.col="grey", shadow.size=0.005);vis_thresh_repeat_subset
# dev.off()
# 
# dev.off()
# png("Plots/script4_visnetwork_subset.png", width=750, height=750)
# vis_subset <- plotmat(round(t(tm_subset), digits=2), relsize=1.1,
#                       box.type="circle", box.size=0.02,
#                       shadow.col="grey", shadow.size=0.005);vis_subset
# dev.off()
# png("Plots/script4_visnetwork_thresh_subset.png", width=750, height=750)
# vis_thresh_subset <- plotmat(round(t(tm_thresh_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_thresh_subset
# dev.off()
# 
# dev.off()
# png("Plots/script4_visnetwork_repeat_subset.png", width=750, height=750)
# vis_repeat_subset <- plotmat(round(t(tm_repeat_subset), digits=2), relsize=1.1,
#                              box.type="circle", box.size=0.02,
#                              shadow.col="grey", shadow.size=0.005);vis_repeat_subset
# dev.off()
# png("Plots/script4_visnetwork_thresh_repeat_subset.png", width=950, height=950)
# vis_thresh_repeat_subset <- plotmat(round(t(tm_thresh_repeat_subset), digits=2), relsize=1.1,
#                                     box.type="circle", box.size=0.02,
#                                     shadow.col="grey", shadow.size=0.005,
#                                     cex = 1.4, cex.txt = 1.4, box.cex=1.6,
#                                     arr.length=0.5, arr.lcol="darkseagreen", arr.col="darkseagreen");vis_thresh_repeat_subset
# dev.off()

#####

tm_thresh_subset_xx <- ifelse(tm_subset_xx <= 0.12, 0, tm_subset_xx) #mean(tm_subset) + sd(tm_subset)
tm_thresh_repeat_subset_xx <- ifelse(tm_repeat_subset_xx <= 0.12, 0, tm_repeat_subset_xx) #mean(tm_subset) + sd(tm_subset)

# dev.off()
# pdf("Plots/script4_visnetwork_subset_xx.pdf")
# vis_subset_xx <- plotmat(round(t(tm_subset_xx), digits=2), relsize=1.1,
#                          box.type="circle", box.size=0.02,
#                          shadow.col="grey", shadow.size=0.005);vis_subset_xx
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_subset_xx.pdf")
# vis_thresh_subset_xx <- plotmat(round(t(tm_thresh_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_thresh_subset_xx
# dev.off()
# dev.off()
# png("Plots/script4_visnetwork_subset_xx.png")
# vis_subset_xx <- plotmat(round(t(tm_subset_xx), digits=2), relsize=1.1,
#                          box.type="circle", box.size=0.02,
#                          shadow.col="grey", shadow.size=0.005);vis_subset_xx
# dev.off()
# png("Plots/script4_visnetwork_thresh_subset_xx.png")
# vis_thresh_subset_xx <- plotmat(round(t(tm_thresh_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_thresh_subset_xx
# dev.off()
# 
# dev.off()
# pdf("Plots/script4_visnetwork_repeat_subset_xx.pdf")
# vis_repeat_subset_xx <- plotmat(round(t(tm_repeat_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_repeat_subset_xx
# dev.off()
# pdf("Plots/script4_visnetwork_thresh_repeat_subset_xx.pdf")
# vis_thresh_repeat_subset_xx <- plotmat(round(t(tm_thresh_repeat_subset_xx), digits=2), relsize=1.1,
#                                        box.type="circle", box.size=0.02,
#                                        shadow.col="grey", shadow.size=0.005);vis_thresh_repeat_subset_xx
# dev.off()
# dev.off()
# png("Plots/script4_visnetwork_repeat_subset_xx.png")
# vis_repeat_subset_xx <- plotmat(round(t(tm_repeat_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_repeat_subset_xx
# dev.off()
# png("Plots/script4_visnetwork_thresh_repeat_subset_xx.png")
# vis_thresh_repeat_subset_xx <- plotmat(round(t(tm_thresh_repeat_subset_xx), digits=2), relsize=1.1,
#                                        box.type="circle", box.size=0.02,
#                                        shadow.col="grey", shadow.size=0.005);vis_thresh_repeat_subset_xx
# dev.off()
# 
# 
# dev.off()
# png("Plots/script4_visnetwork_subset_xx.png", width=750, height=750)
# vis_subset_xx <- plotmat(round(t(tm_subset_xx), digits=2), relsize=1.1,
#                          box.type="circle", box.size=0.02,
#                          shadow.col="grey", shadow.size=0.005);vis_subset_xx
# dev.off()
# png("Plots/script4_visnetwork_thresh_subset_xx.png", width=750, height=750)
# vis_thresh_subset_xx <- plotmat(round(t(tm_thresh_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_thresh_subset_xx
# dev.off()
# 
# dev.off()
# png("Plots/script4_visnetwork_repeat_subset_xx.png", width=750, height=750)
# vis_repeat_subset_xx <- plotmat(round(t(tm_repeat_subset_xx), digits=2), relsize=1.1,
#                                 box.type="circle", box.size=0.02,
#                                 shadow.col="grey", shadow.size=0.005);vis_repeat_subset_xx
# dev.off()
# png("Plots/script4_visnetwork_thresh_repeat_subset_xx.png", width=950, height=950)
# vis_thresh_repeat_subset_xx <- plotmat(round(t(tm_thresh_repeat_subset_xx), digits=2), relsize=1.1,
#                                        box.type="circle", box.size=0.02,
#                                        shadow.col="grey", shadow.size=0.005,
#                                        cex = 1.4, cex.txt = 1.4, box.cex=1.6,
#                                        arr.length=0.5, arr.lcol="darkseagreen", arr.col="darkseagreen");vis_thresh_repeat_subset_xx
# dev.off()

#####################################
#### Predict a new constellation ####

#Select one of the state strings from the lists, or make one up:
names(statestringlist_repeat)
newconst <- statestringlist_repeat[sample(1:length(statestringlist_repeat), 1)];newconst
newconst <- statestringlist_repeat_xx$H1N2TVVPPT;newconst
newconst <- c("NE", "OK", "MO", "KS");newconst

prior_state <- tail(unlist(newconst), 2)[1];prior_state
tm_predict_repeat <- tm_repeat[which(rownames(tm_repeat)==prior_state),]
tm_predict <- tm[which(rownames(tm)==prior_state),which(colnames(tm) %!in% newconst)]
tm_predict_repeat_new <- tm_predict_repeat[which(names(tm_predict_repeat) %!in% newconst)] # novel states from the _repeat model

# terra::predict(est, prior_state) #transmission to a new state # This method is flawed because it only returns one state when there are ties.
# terra::predict(est_repeat, prior_state) #next detection location

print("Most likely state of next detection:");head(sort(round(tm_predict_repeat, digits=3), decreasing=TRUE));names(which(tm_predict_repeat == max(tm_predict_repeat))) #next detection location
print("Most likely next novel state:");head(sort(round(tm_predict_repeat_new, digits=3), decreasing=TRUE));names(which(tm_predict_repeat_new == max(tm_predict_repeat_new))) #next novel detection location

#####################################
### Summary statistics of tm_repeat

rowSums(tm_repeat)
sort(colSums(tm_repeat))
mean(tm_repeat)
tail(sort(colMeans(tm_repeat)), n=9)
tail(sort(colMeans(tm_repeat_xx)), n=9)

#probability of next detection in same state vs. new state
tm_repeat_nodiag <- tm_repeat
diag(tm_repeat_nodiag) <- NA

mean(tm_repeat_nodiag, na.rm=T) # mean probability of being detected in a specific different state from before
mean(rowSums(tm_repeat_nodiag, na.rm=T)) # mean probability of not being detected next in same state as before
mean(diag(tm_repeat)) # mean probability of being detected in the same state as before

#probability of next detection in Iowa vs. anywhere else
mean(tm_repeat[,which(colnames(tm_repeat) == "IA")]) # mean probability of being detected next in Iowa
mean(tm_repeat_xx[,which(colnames(tm_repeat_xx) == "IA")]) # mean probability of being detected next in Iowa, allowing XX
mean(rowSums(tm_repeat[,which(colnames(tm_repeat) != "IA")])) # mean probability of not being detected next in Iowa
mean(tm_repeat[,which(colnames(tm_repeat) != "IA")]) # mean probability of being detected next in any other state

#probability of next "detection" being an XX (i.e., probability that a given detection is the last of its kind)
mean(tm_repeat_xx[,which(colnames(tm_repeat_xx) == "XX")])

#probability that you end up at XX, based on your FIRST/STARTING state
statestringlist_repeat_xx

i <- j <- NULL
state_counter_list_master <- xx_counter_list_master <- list()
state_counter <- 0
xx_counter <- 0
for(i in 1:length(unique(data$State))){ #for each state,
  statei <- unique(data$State)[i]
  state_counter_list <- xx_counter_list <- list()
  for(k in 1:length(statestringlist_repeat_xx)){ #for each virus,
    statestringlist_repeat_xx_k <- statestringlist_repeat_xx[[k]]
    if(statestringlist_repeat_xx_k[1] == statei){
      state_counter <- state_counter+1
    }
    if(statestringlist_repeat_xx_k[1] == statei & statestringlist_repeat_xx_k[length(statestringlist_repeat_xx_k)] == "XX"){
      xx_counter <- xx_counter+1
    }
  }
  state_counter_list <- append(state_counter_list,state_counter)
  xx_counter_list <- append(xx_counter_list,xx_counter)
  state_counter_list_master[[i]] <- unlist(state_counter_list)
  xx_counter_list_master[[i]] <- unlist(xx_counter_list)
  state_counter <- 0
  xx_counter <- 0
}

names(xx_counter_list_master) <- unique(data$State)
names(state_counter_list_master) <- unique(data$State)
dudd <- cbind(unlist(state_counter_list_master),unlist(xx_counter_list_master))
dudd <-  unlist(xx_counter_list_master)/unlist(state_counter_list_master)
names(dudd) <- unique(data$State)
sort(dudd)
plot(dudd)

#####################################
### Evaluating transition matrices for accuracy

statestringlist_no1 <- statestringlist[sapply(statestringlist, length) > 1];length(statestringlist_no1)
statestringlist_repeat_no1 <- statestringlist_repeat[sapply(statestringlist_repeat, length) > 1];length(statestringlist_repeat_no1)
statestringlist_xx_no1 <- statestringlist_xx[sapply(statestringlist_xx, length) > 1];length(statestringlist_xx_no1)                             #exclude, of no value: almost all constellations end in XX
statestringlist_repeat_xx_no1 <- statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, length) > 1];length(statestringlist_repeat_xx_no1)
statestringlist_subset_no1 <- statestringlist_subset[sapply(statestringlist_subset, length) > 1];length(statestringlist_subset_no1)
statestringlist_repeat_subset_no1 <- statestringlist_repeat_subset[sapply(statestringlist_repeat_subset, length) > 1];length(statestringlist_repeat_subset_no1)

statestring_masterlist <- list(statestringlist_no1,                #no repeated states, minus strings with length of 1, XX excluded
                               statestringlist_repeat_no1,         #repeated states, minus strings with length of 1, XX excluded
                               statestringlist_xx_no1,             #no repeated states, minus strings with length of 1, XX included
                               statestringlist_repeat_xx_no1,      #repeated states, minus strings with length of 1, XX included
                               statestringlist_subset_no1,         #top pork producers, no repeats, minus strings with length of 1, XX excluded
                               statestringlist_repeat_subset_no1)  #top pork producers, repeated states, minus strings with length of 1, XX excluded
descriptor_repeat <- c("no",
                       "repeat",
                       "no",
                       "repeat",
                       "no",
                       "repeat")
descriptor_xx <- c("no",
                   "no",
                   "xx",
                   "xx",
                   "no",
                   "no")
mantel_mean_list <- list()
counter_list_master <- list()
trials_list_master <- list()
counter <- 0 #needs to be defined here, in case the first constellation tries to add 1 to counter. It is reset to 0 appropriately later in the script.

#Logic:
#For each of these four lists,
#make the 4 training and 4 testing lists
#make the 4 mc eval models
#rectify the 4 testing lists with setdiff and the for loop i made before
#then do the actual accuracy test

i <- j <- k <- NULL
for(i in 1:length(statestring_masterlist)){
  statestringlisti <- statestring_masterlist[[i]]
  tm_eval_list <- list()
  statestring_test_list <- list()
  
  #make the four evaluatory tms for statestring i:
  for(j in 1:4){
    statestringlist_train <- statestringlisti[-seq(j, length(statestringlisti), 4)]
    statestringlist_test <- statestringlisti[seq(j, length(statestringlisti), 4)]
    
    m_eval <- do.call("rbind", lapply(statestringlist_train, function(x) cbind(head(x, -1), tail(x, -1))))
    mc_eval <- markovchainFit(m_eval)
    est_eval <- mc_eval$estimate
    tm_eval <- est_eval@transitionMatrix
    tm_melt_eval <- melt(tm_eval)
    colnames(tm_melt_eval) <- c("Origin", "Destination", "Probability")
    hmplot_eval <- ggplot(tm_melt_eval, aes(Destination, Origin, fill = Probability)) +
      geom_tile(color = "black") +
      scale_fill_gradient(low = "white", high = "orange") +
      geom_tile() +
      labs(title="Markov chain transition matrix");hmplot_eval
    
    tm_eval_list[[j]] <- tm_eval
    statestring_test_list[[j]] <- statestringlist_test
  }
  
  #fix the evaluatory tms by making sure they contain the same states:
  states_eval1 <- setdiff(union(colnames(as.data.frame(tm_eval_list[2])), union(colnames(as.data.frame(tm_eval_list[3])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[1]))) 
  states_eval2 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[3])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[2])))
  states_eval3 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[2])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[3])))
  states_eval4 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[2])),
                                                                                colnames(as.data.frame(tm_eval_list[3])))), colnames(as.data.frame(tm_eval_list[4])))
  # print(states_eval1);  print(states_eval2);  print(states_eval3);  print(states_eval4)
  tm_eval1_df <- as.data.frame(tm_eval_list[[1]])
  tm_eval2_df <- as.data.frame(tm_eval_list[[2]])
  tm_eval3_df <- as.data.frame(tm_eval_list[[3]])
  tm_eval4_df <- as.data.frame(tm_eval_list[[4]])
  
  states_eval_list <- list(states_eval1, states_eval2, states_eval3, states_eval4)
  tm_eval_list_df <- list(tm_eval1_df, tm_eval2_df, tm_eval3_df, tm_eval4_df)
  
  j <- k <- NULL
  for(j in 1:length(tm_eval_list_df)){                     #for each evaluatory tm,
    tmj <- tm_eval_list_df[[j]]                            #name the evaluatory tm j "tmj"
    if(length(states_eval_list[[j]]) > 0){                 #if tmj is missing more than 0 states,
      for(k in 1:length(states_eval_list[[j]])){           #for each missing state k,
        # k <- 1
        st <- states_eval_list[[j]][k]                     #define the state
        tmj$newcol <- 0                                    #add new column
        tmj <- rbind(tmj, 0)                               #add new row
        colnames(tmj)[ncol(tmj)] <- st
        rownames(tmj)[nrow(tmj)] <- st
        st <- NULL
      }
      tm_eval_list_df[[j]] <- tmj                          #overwrite the element in the list with the fixed tmj object
    }
  }
  
  new_order <- sort(colnames(tm_eval_list_df[[1]]))
  for(j in 1:4){
    tm_eval_list_df[[j]] <- tm_eval_list_df[[j]][new_order, new_order]
  }

  #begin the testing:
  j <- k <- NULL
  counter_list <- list()                                                                            #start a list to store the counters
  trials_list <- list()                                                                              #start a list to store the trials
  
  for(j in 1:length(statestring_test_list)){                                                        #for each 1 of 4 eval/testing group pairs,
    tm_eval_colmeans <- colMeans(tm_eval_list_df[[j]])
    
    for(k in 1:length(statestring_test_list[[j]])){                                                 #for each constellation in the test group,
      constellation <- names(statestring_test_list[[j]][k])                                            #pick the constellation,
      constellation_statestring <- as.data.frame(unlist(statestring_test_list[[j]][k]))[,1]         #pull that constellation's state sequence,
      state_final_actual <- tail(constellation_statestring, 1);state_final_actual                   #take the final state,
      prior_state <- tail(constellation_statestring, 2)[1];prior_state                              #and the second to final state.
      
      tm_eval_sub <- tm_eval_list_df[[j]][which(rownames(tm_eval_list_df[[j]])==prior_state),];tm_eval_sub    #tm_eval row of the prior state
      
      if(descriptor_repeat[i] == "no"){                                                                                             #if we are disallowing repeats,
        tm_eval_sub <- tm_eval_sub[names(tm_eval_sub) %!in% head(constellation_statestring, n=length(constellation_statestring)-1)] #the pool of possible next states cannot be in the constellation state string,
                                                                                                                                    #minus the final state
      }
      
      state_final_predict <- colnames(tm_eval_sub)[which(tm_eval_sub == max(tm_eval_sub))];state_final_predict
      
      if(length(state_final_predict) > 1){                                                                    #if it predicts 2+ tied states,
        state_final_predict <- tm_eval_colmeans[state_final_predict]                                          #pull colmeans of tied states,
        state_final_predict <- names(state_final_predict[state_final_predict == max(state_final_predict)])[1] #and pick the tied state with higher colmean. If still tied, pick the first one.
      }
      
      if(state_final_predict == state_final_actual){
        counter <- counter+1                                                                        #add to counter iff prediction matches actual
      }
      
    }
    
    counter_list <- append(counter_list,counter)
    trials_list <- append(trials_list, k)
    counter <- 0
  }
  counter_list_master[[i]] <- unlist(counter_list)
  trials_list_master[[i]] <- unlist(trials_list)
}

accuracy <- sum(counter_list_master[[1]])/sum(trials_list_master[[1]]);accuracy                              #predict final state in string, barring repeats, barring XX
# accuracy_xx <- sum(counter_list_master[[3]])/sum(trials_list_master[[3]]);accuracy_xx                        #predict final state in string, barring repeats, permitting XX
accuracy_subset <- sum(counter_list_master[[5]])/sum(trials_list_master[[5]]);accuracy_subset                #predict final state in string, given top 10 swine states, barring repeats, barring XX

accuracy_repeat <- sum(counter_list_master[[2]])/sum(trials_list_master[[2]]);accuracy_repeat                #predict final state in string, permitting repeats, barring XX
# accuracy_xx_repeat <- sum(counter_list_master[[4]])/sum(trials_list_master[[4]]);accuracy_xx_repeat          #predict final state in string, permitting repeats, permitting XX
accuracy_subset_repeat <- sum(counter_list_master[[6]])/sum(trials_list_master[[6]]);accuracy_subset_repeat  #predict final state in string, given top 10 swine states, permitting repeats, barring XX

#####################################################
### Final Tests for predicting extinction events (XX)
### Test 1: Given the first detection of a new virus, predict whether it will ever be detected again (ie, whether next state is XX or notXX)
###### Requires state lists with 2+ length, counting XX
### Test 2: Given the second detection of a virus, predict whether it will ever be detected again (ie, whether next state is XX or notXX)
###### Requires state lists with 3+ length, counting XX
### Test 3: Given a random detection somewhere in a string of detections, predict whether it will ever be detected again (ie, whether next state is XX or notXX)
###### Requires state lists with 2+ length, counting XX
### These three tests work whether we are looking at all detections or only novel states.
### There will be 3x2=6 accuracy metrics. Each metric is contrasted with the accuracy of always guessing XX and the accuracy of always guessing notXX.
### Tests 1 and 3 can use the same lists, whereas Test 2 has to be more selective.
### Tests 1 and 3 must include viruses with string lengths of 2 or more. XX is allowed, not required.
### Test 2 must include viruses with string lengths of 3 or more. XX is allowed, not required. 

statestringlist_test1 <- statestringlist_test3 <- statestringlist_xx_no1;table(lengths(statestringlist_test1))
statestringlist_repeat_test1 <- statestringlist_repeat_test3 <- statestringlist_repeat_xx_no1;table(lengths(statestringlist_repeat_test1))
statestringlist_test2 <- statestringlist_xx[sapply(statestringlist_xx, length) > 2];table(lengths(statestringlist_test2))
statestringlist_repeat_test2 <- statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, length) > 2];table(lengths(statestringlist_repeat_test2))

statestring_masterlist_xx <- list(statestringlist_test1,
                                  statestringlist_repeat_test1,
                                  statestringlist_test2,
                                  statestringlist_repeat_test2,
                                  statestringlist_test3,
                                  statestringlist_repeat_test3)
descriptor_repeat <- c("no",
                       "repeat",
                       "no",
                       "repeat",
                       "no",
                       "repeat")
descriptor_test <- c("test1",
                     "test1",
                     "test2",
                     "test2",
                     "test3",
                     "test3")
counter <- counter_guessXX <- counter_guessnotXX <- 0
counter_list_master <- counter_list_guessXX_master <- counter_list_guessnotXX_master <- trials_list_master <- list()
i <- j <- k <- NULL

set.seed(10);random_numbers <- sample(0:100, size = 1000, replace = TRUE)

for(i in 1:length(statestring_masterlist_xx)){
  statestringlisti <- statestring_masterlist_xx[[i]]
  tm_eval_list <- list()
  statestring_test_list <- list()
  #make the four evaluatory tms for statestring i:
  for(j in 1:4){
    statestringlist_train <- statestringlisti[-seq(j, length(statestringlisti), 4)]
    statestringlist_test <- statestringlisti[seq(j, length(statestringlisti), 4)]
    m_eval <- do.call("rbind", lapply(statestringlist_train, function(x) cbind(head(x, -1), tail(x, -1))))
    mc_eval <- markovchainFit(m_eval)
    est_eval <- mc_eval$estimate
    tm_eval <- est_eval@transitionMatrix
    tm_melt_eval <- melt(tm_eval)
    colnames(tm_melt_eval) <- c("Origin", "Destination", "Probability")
    hmplot_eval <- ggplot(tm_melt_eval, aes(Destination, Origin, fill = Probability)) +
      geom_tile(color = "black") +
      scale_fill_gradient(low = "white", high = "orange") +
      geom_tile() +
      labs(title="Markov chain transition matrix");hmplot_eval
    
    tm_eval_list[[j]] <- tm_eval
    statestring_test_list[[j]] <- statestringlist_test
  }
  #fix the evaluatory tms by making sure they contain the same states:
  states_eval1 <- setdiff(union(colnames(as.data.frame(tm_eval_list[2])), union(colnames(as.data.frame(tm_eval_list[3])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[1]))) 
  states_eval2 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[3])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[2])))
  states_eval3 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[2])),
                                                                                colnames(as.data.frame(tm_eval_list[4])))), colnames(as.data.frame(tm_eval_list[3])))
  states_eval4 <- setdiff(union(colnames(as.data.frame(tm_eval_list[1])), union(colnames(as.data.frame(tm_eval_list[2])),
                                                                                colnames(as.data.frame(tm_eval_list[3])))), colnames(as.data.frame(tm_eval_list[4])))
  # print(states_eval1);  print(states_eval2);  print(states_eval3);  print(states_eval4)
  tm_eval1_df <- as.data.frame(tm_eval_list[[1]])
  tm_eval2_df <- as.data.frame(tm_eval_list[[2]])
  tm_eval3_df <- as.data.frame(tm_eval_list[[3]])
  tm_eval4_df <- as.data.frame(tm_eval_list[[4]])
  states_eval_list <- list(states_eval1, states_eval2, states_eval3, states_eval4)
  tm_eval_list_df <- list(tm_eval1_df, tm_eval2_df, tm_eval3_df, tm_eval4_df)
  j <- k <- NULL
  for(j in 1:length(tm_eval_list_df)){                     #for each evaluatory tm,
    tmj <- tm_eval_list_df[[j]]                            #name the evaluatory tm j "tmj"
    if(length(states_eval_list[[j]]) > 0){                 #if tmj is missing more than 0 states,
      for(k in 1:length(states_eval_list[[j]])){           #for each missing state k,
        st <- states_eval_list[[j]][k]                     #define the state
        tmj$newcol <- 0                                    #add new column
        tmj <- rbind(tmj, 0)                               #add new row
        colnames(tmj)[ncol(tmj)] <- st
        rownames(tmj)[nrow(tmj)] <- st
        st <- NULL
      }
      tm_eval_list_df[[j]] <- tmj                          #overwrite the element in the list with the fixed tmj object
    }
  }
  new_order <- sort(colnames(tm_eval_list_df[[1]]))
  for(j in 1:4){
    tm_eval_list_df[[j]] <- tm_eval_list_df[[j]][new_order, new_order]
  }
  
  #begin the testing:
  j <- k <- NULL
  counter_list <- list()                                                                            #start a list to store the counters
  counter_list_guessXX <- list()                                                                    #start a list to store the counters
  counter_list_guessnotXX <- list()                                                                 #start a list to store the counters
  trials_list <- list()                                                                             #start a list to store the trials
  for(j in 1:length(statestring_test_list)){                                                        #for each 1 of 4 eval/testing group pairs,
    tm_eval_colmeans <- colMeans(tm_eval_list_df[[j]])
    for(k in 1:length(statestring_test_list[[j]])){                                                 #for each constellation in the test group,
      constellation <- names(statestring_test_list[[j]][k])                                         #pick the constellation,
      constellation_statestring <- as.data.frame(unlist(statestring_test_list[[j]][k]))[,1]         #pull that constellation's state sequence,
      set.seed(random_numbers[k]);rand <- sample(1:(length(constellation_statestring)-1), 1)        #pick a random number 1:(length of string-1)

      switch(descriptor_test[i],
             "test1" = {
               state_current <- constellation_statestring[1]
               state_next <- constellation_statestring[2]
             },
             "test2" = {
               state_current <- constellation_statestring[2]
               state_next <- constellation_statestring[3]
             },
             "test3" = {
               state_current <- constellation_statestring[rand]
               state_next <- constellation_statestring[rand+1]
             }
      )
      
      tm_eval_sub <- tm_eval_list_df[[j]][which(rownames(tm_eval_list_df[[j]])==state_current),];tm_eval_sub                      #tm_eval row of the prior state
      
      #if we are disallowing repeats, remove certain elements from consideration:
      if(descriptor_repeat[i] == "no"){
        if(descriptor_test[i] == "test1"){
          tm_eval_sub <- tm_eval_sub[names(tm_eval_sub) %!in% head(constellation_statestring, n=1)]
        }
        if(descriptor_test[i] == "test2"){
          tm_eval_sub <- tm_eval_sub[names(tm_eval_sub) %!in% head(constellation_statestring, n=2)]
        }
        if(descriptor_test[i] == "test3"){
          tm_eval_sub <- tm_eval_sub[names(tm_eval_sub) %!in% head(constellation_statestring, n=rand)]
        }
      }
      
      state_next_predict <- colnames(tm_eval_sub)[which(tm_eval_sub == max(tm_eval_sub))];state_next_predict
      
      if(length(state_next_predict) > 1){                                                                   #if it predicts 2+ tied states,
        state_next_predict <- tm_eval_colmeans[state_next_predict]                                          #pull colmeans of tied states,
        state_next_predict <- names(state_next_predict[state_next_predict == max(state_next_predict)])[1] #and pick the tied state with higher colmean. If still tied, pick the first one.
      }
      
      mcguess <- if(state_next_predict == "XX") {
        "XX"
      } else {
        "notXX"
      }
      
      if ((mcguess == "XX" && state_next == "XX") || 
          (mcguess == "notXX" && state_next != "XX")) {
        counter <- counter + 1
      }
      
      counter_guessXX <- counter_guessXX + (state_next == "XX")
      counter_guessnotXX <- counter_guessnotXX + (state_next != "XX")
    }
    counter_list <- append(counter_list,counter)
    counter_list_guessXX <- append(counter_list_guessXX,counter_guessXX)
    counter_list_guessnotXX <- append(counter_list_guessnotXX,counter_guessnotXX)
    trials_list <- append(trials_list, k)
    counter <- counter_guessXX <- counter_guessnotXX <- 0    
  }
  counter_list_master[[i]] <- unlist(counter_list)
  counter_list_guessXX_master[[i]] <- unlist(counter_list_guessXX)
  counter_list_guessnotXX_master[[i]] <- unlist(counter_list_guessnotXX)
  trials_list_master[[i]] <- unlist(trials_list)
}

###########################################
###########################################
### Test 1: Given the first detection of a new virus, predict whether it will ever be detected again (ie, whether next state is XX or notXX)
### Test 2: Given the second detection of a virus, predict whether it will ever be detected again (ie, whether next state is XX or notXX)
### Test 3: Given a random detection somewhere in a string of detections, predict whether it will ever be detected again (ie, whether next state is XX or notXX)

accuracy_test1 <- sum(counter_list_master[[1]])/sum(trials_list_master[[1]]);accuracy_test1                                    #predict whether new virus will spread to a new state
accuracy_guessxx_test1 <- sum(counter_list_guessXX_master[[1]])/sum(trials_list_master[[1]]);accuracy_guessxx_test1            #guessing it will not spread to a new state
accuracy_guessnotxx_test1 <- sum(counter_list_guessnotXX_master[[1]])/sum(trials_list_master[[1]]);accuracy_guessnotxx_test1   #guessing it will spread to a new state

accuracy_test1_repeat <- sum(counter_list_master[[2]])/sum(trials_list_master[[2]]);accuracy_test1_repeat                                    #predict whether new virus will go extinct
accuracy_guessxx_test1_repeat <- sum(counter_list_guessXX_master[[2]])/sum(trials_list_master[[2]]);accuracy_guessxx_test1_repeat            #guessing it will go extinct
accuracy_guessnotxx_test1_repeat <- sum(counter_list_guessnotXX_master[[2]])/sum(trials_list_master[[2]]);accuracy_guessnotxx_test1_repeat   #guessing it will not go extinct

###############

accuracy_test2 <- sum(counter_list_master[[3]])/sum(trials_list_master[[3]]);accuracy_test2                                    #predict whether virus with 2 detections will spread to a new state
accuracy_guessxx_test2 <- sum(counter_list_guessXX_master[[3]])/sum(trials_list_master[[3]]);accuracy_guessxx_test2            #guessing it will not spread to a new state
accuracy_guessnotxx_test2 <- sum(counter_list_guessnotXX_master[[3]])/sum(trials_list_master[[3]]);accuracy_guessnotxx_test2   #guessing it will spread to a new state

accuracy_test2_repeat <- sum(counter_list_master[[4]])/sum(trials_list_master[[4]]);accuracy_test2_repeat                                    #predict whether virus with 2 detections will go extinct
accuracy_guessxx_test2_repeat <- sum(counter_list_guessXX_master[[4]])/sum(trials_list_master[[4]]);accuracy_guessxx_test2_repeat            #guessing it will go extinct
accuracy_guessnotxx_test2_repeat <- sum(counter_list_guessnotXX_master[[4]])/sum(trials_list_master[[4]]);accuracy_guessnotxx_test2_repeat   #guessing it will not go extinct

################

# accuracy_test3  <- sum(counter_list_master[[5]])/sum(trials_list_master[[5]]);accuracy_test3                                    #predict whether virus with unknown history will spread to a new state
# accuracy_guessxx_test3  <- sum(counter_list_guessXX_master[[5]] )/sum(trials_list_master[[5]]);accuracy_guessxx_test3           #guessing it will not spread to a new state
# accuracy_guessnotxx_test3  <- sum(counter_list_guessnotXX_master[[5]] )/sum(trials_list_master[[5]]);accuracy_guessnotxx_test3  #guessing it will spread to a new state

accuracy_test3_repeat <- sum(counter_list_master[[6]])/sum(trials_list_master[[6]] );accuracy_test3_repeat                                   #predict whether virus with unknown history will go extinct
accuracy_guessxx_test3_repeat <- sum(counter_list_guessXX_master[[6]])/sum(trials_list_master[[6]]);accuracy_guessxx_test3_repeat            #guessing it will go extinct
accuracy_guessnotxx_test3_repeat <- sum(counter_list_guessnotXX_master[[6]])/sum(trials_list_master[[6]]);accuracy_guessnotxx_test3_repeat   #guessing it will not go extinct

###########################################
###########################################
### Extra insights:

#what fraction of viruses have only been detected once?
length(statestringlist_repeat[sapply(statestringlist_repeat, function(lst)length(lst) == 1)])/length(statestringlist_repeat) 
#what fraction of viruses are only ever detected once, then go extinct?
length(statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, function(lst)length(lst) == 2 && tail(lst, 1) == "XX")])/length(statestringlist_repeat_xx)
#what fraction of single-detection viruses are we confident are extinct?
length(statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, function(lst)length(lst) == 2 && tail(lst, 1) == "XX")])/
  length(statestringlist_repeat[sapply(statestringlist_repeat, function(lst)length(lst) == 1)])

#what fraction of viruses have only been detected twice?
length(statestringlist_repeat[sapply(statestringlist_repeat, function(lst)length(lst) == 2)])/length(statestringlist_repeat) 
#what fraction of viruses are only ever detected twice, then go extinct?
length(statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, function(lst)length(lst) == 3 && tail(lst, 1) == "XX")])/length(statestringlist_repeat_xx)

#what fraction of viruses went extinct at some point?
length(statestringlist_repeat_xx[sapply(statestringlist_repeat_xx, function(lst) tail(lst, 1) == "XX")])/length(statestringlist_repeat_xx)                       

#what fraction of viruses are only ever detected in a single state?
length(statestringlist[sapply(statestringlist, function(lst)length(lst) == 1)])/length(statestringlist) 
#what fraction of viruses are only ever detected in a single state, then go extinct?
length(statestringlist_xx[sapply(statestringlist_xx, function(lst)length(lst) == 2 && tail(lst, 1) == "XX")])/length(statestringlist_xx)
#what fraction of viruses only detected in a single state are we confident went extinct?
length(statestringlist_xx[sapply(statestringlist_xx, function(lst)length(lst) == 2 && tail(lst, 1) == "XX")])/
  length(statestringlist[sapply(statestringlist, function(lst)length(lst) == 1)]) 

###################

#######
# save.image("IAV_Sources_and_Sinks.RData")
#######

# https://datascience.blog.wzb.eu/2018/05/31/three-ways-of-visualizing-a-graph-on-a-map/

###########################
#Modify the following line to select the tm you wish to plot, labeling it tm_map:
tm_map <- tm_reduced
#Set the value of the lower threshold for edge inclusion. 0 works for smaller tms.
crit_thresh <- 0.14
###########################

colnames(centroids)[2] <- "lon"
centroids <- merge(centroids, data_usda_agg, all.x=TRUE, all.y=FALSE, by.x="postal_code", by.y="State")
centroids$state <- tolower(centroids$state)
centroids_formerge <- centroids[,4:5]
colnames(centroids_formerge) <- c("state","Swine Inventory (million)")
centroids_formerge$`Swine Inventory (million)` <- centroids_formerge$`Swine Inventory (million)`/1000000

tm_try <- as.data.frame(c(tm_map))
colnames(tm_try) <- "probability"
tm_try$to <- rep(rownames(tm_map),each=ncol(tm_map))
tm_try$from <- rep(colnames(tm_map),nrow(tm_map))
tm_try2 <- merge(tm_try, centroids, all.x=TRUE, all.y=FALSE, by.x="from", by.y="postal_code");head(tm_try2)
tm_try2$Production <- NULL
colnames(tm_try2)[4:5] <- c("y","x")
tm_try3 <- merge(tm_try2, centroids, all.x=TRUE, all.y=FALSE, by.x="to", by.y="postal_code");head(tm_try3)
tm_try3$Production <- NULL
colnames(tm_try3)[7:8] <- c("yend","xend")
tm_try4 <- tm_try3[tm_try3$y!=tm_try3$yend,]
tm_try4 <- tm_try4[tm_try4$x!=tm_try4$xend,]
tm_plot <- plot(tm_try4$probability, size=2);tm_plot
tm_try5 <- tm_try4[tm_try4$probability >= crit_thresh,] #if the probability is below the critical threshold defined above, drop it.
tm_try5$Production.y <- NULL
rm(tm_try, tm_try2, tm_try3, tm_try4)

tm_try5$Category <- "other"
tm_try5$Category <- ifelse(tm_try5$to == "IA", "to Iowa", tm_try5$Category)
tm_try5$Category <- ifelse(tm_try5$from == "IA", "from Iowa", tm_try5$Category)
tm_try5$Category <- factor(tm_try5$Category, levels=c("to Iowa","from Iowa","other"))
colnames(tm_try5)[3] <- "Probability"

centroids_try <- centroids[which(centroids$postal_code %in% tm_try5$to | centroids$postal_code %in% tm_try5$from),]
g <- graph_from_data_frame(tm_try5, directed = FALSE, vertices = centroids_try)

all_states  <- map_data('state')
all_states$id  <- 1:nrow(all_states)
stateData <- merge(all_states,centroids_formerge,by.x="region",by.y="state",all.x=TRUE, all.y=FALSE)
stateData <- stateData[order(stateData$id), ]
table(stateData$`Swine Inventory (million)`);plot(stateData$`Swine Inventory (million)`)

usa_shapes <- geom_polygon(data = stateData, aes(x = long, y = lat, group = group, fill = `Swine Inventory (million)`), size = 0.15)
centroids_try$Production = degree(g)

mapcoords <- coord_fixed(xlim = c(-103, -76.5), ylim = c(33, 48))

maptheme <- theme(panel.grid = element_blank()) + 
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm')) +
  theme(legend.title = element_text(size = 17), 
        legend.text = element_text(size = 10))

map1 <- ggplot(centroids_try) + usa_shapes + maptheme + mapcoords + 
  scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3", breaks=c(0,3,6,9,12)) +
  geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
             aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Probability),
             curvature = 0.33,
             alpha = 0.85, show.legend=TRUE) +
  scale_color_gradient(name="Transition probability",low="#CCFFFF", high="dodgerblue4", breaks=c(0.1,0.2,0.3,0.4)) +
  scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
  geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
            hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
            size = 7, color = "black", fontface = "bold");map1

map2 <- ggplot(centroids_try) + usa_shapes + maptheme + mapcoords +
  scale_fill_gradient(name="Swine inventory (million)",low="grey30",high="pink3") +
  geom_curve(data=tm_try5, arrow = arrow(length = unit(0.03, "npc"), type="closed"),
             aes(x = x, y = y, xend = xend, yend = yend, size=Probability, color=Category),
             curvature = 0.33,
             alpha = 0.85, show.legend=TRUE) +
  scale_color_manual(values=c("#C8102E", "gold", "gray75")) +
  guides(color = guide_legend(override.aes = list(size = 20))) +
  scale_size_continuous(guide = FALSE, range = c(1, 4)) + 
  geom_text(data = centroids_try, aes(x = lon, y = lat, label = postal_code), 
            hjust = 0, nudge_x = 0.2, nudge_y = -0.2,
            size = 7, color = "black", fontface = "bold");map2

pdf("Plots/script4_map1_probability.pdf", height=7.2, width=12)
map1
dev.off()
png("Plots/script4_map1_probability.png", height=720, width=1200)
map1
dev.off()
# pdf("Plots/script4_map2_direction.pdf", height=7.2, width=12)
# map2
# dev.off()
# png("Plots/script4_map2_direction.png", height=720, width=1200)
# map2
# dev.off()

##############################
##############################
##############################
