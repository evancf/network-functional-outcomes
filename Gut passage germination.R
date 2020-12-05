# Packages ---------------------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Load libraries
packages <- c("rjags", "bipartite", "plyr", "plotrix")
ipak(packages)



# Load data ---------------------------------------------------------------------

load("gut passage.RData")


# The data object contains:

# 1 - Several vectors that together make a dataframe. Let's look at it this way:
gut.data <- as.data.frame(data[1:6]) # These are the columns necessary for
# the analysis, but let's also show the treatments
gut.data$treatment <- rep(c("whole", "mech", "gut"), 
                          times = c(length(data$wf.indices),
                                    length(data$mech.indices),
                                    length(data$bird.indices)))
head(gut.data)
# The columns are:
# obs.surv = number of germinants found in the cell throughout the study

# portion.of.fruit = portion of the fruit that was in the cell (some were broken up)

# bird.sp = numeric species of bird. Note that bird names are saved in bird.sp.list

# tree.sp = equivalent for plants

# nested.bird.id = numeric value for the bird individual - note that this
# is nested in the sense that there is fruit dove 1,2,3... and starling 1,2,3...

# nested.tree.id = equivalent for plants

# 2 - A couple vectors describing how many levels of the bird and plant "random
# effects" there were to help with indexing while fitting random effects

# 3 - Several vectors that help deal with multi-seeded fruits:

# lambda.added = this is the mean number of what we describe as "added seeds". This
# accounts for multi-seeded fruits, and we use a Poisson distribution to
# characterize distribution of seed number. But because the Poisson (and other relevant
# distributions) include 0, we need to add this distribution to 1. In other words,
# the average number of seeds is this distribution plus 1. If lambda.added = 5, then
# the mean number of seeds in fruits of that species is 6.

# variable.wf.seed.spp = ones indicate the species that have a variable number of
# seeds per fruit, zeros indicate species with only one seed per cell

# n.min.potential.germ = the lowest number of potential germinants in a cell
# for each species. When there are fused/not-separable seeds, these get planted 
# together, and a single gut/mech seed can produce multiple germinants. Here, only Premna


# 4 - bt.combos. This helps with indexing so that the interaction terms
# get fit for bird-tree combinations we actually observed. Can imagine this
# like a bird ID versus tree ID matrix. So [1,1] is Pipturus & Starling

# 5 - number of bird (n.bird.sp) and plant (n.tree.sp) speces and number of
# bird-plant combinations (n.combos)

# 6 - a few vectors with the indices that relate to the whole fruit, mechanically
# depulped and bird gut-passed treatments.




# Write model ------------------------------------------------------------------

# Note that the model handles indexing with three loops - one for seeds in each
# treatment. This makes the appropriate effects correspond only to the
# appropriate treatments (don't want a bird gut passage effect estimated for
# whole fruit seeds).

sink("gut germ model.txt")
cat("
    model {
    
    ### First loop over the whole fruits (no bird, yes tree random effects)
    for(i in wf.indices){
    
    # Process model
    obs.surv[i] ~ dbinom(p[i], round(n.min.seeds.in.fruit[tree.sp[i]] + n.added.seed.true[i] * variable.wf.seed.spp[tree.sp[i]] * portion.of.fruit[i]))
    logit(p[i]) <- b.whole[tree.sp[i]] + b.tree.id[tree.sp[i],nested.tree.id[i]] 
    
    # Data model
    n.added.seed.true[i] ~ dpois(lambda.added[tree.sp[i]])
    
    }
    
    
    
    ### Second loop over the mech seeds (no bird, yes tree random effects)
    for(i in mech.indices){
    
    obs.surv[i] ~ dbinom(p[i], n.min.potential.germ[tree.sp[i]])
    logit(p[i]) <- b.whole[tree.sp[i]] + b.tree.id[tree.sp[i],nested.tree.id[i]] + b.deinhib[tree.sp[i]]
    
    }
    
    
    
    
    ### Third loop over the bird seeds (yes bird, yes tree random effects)
    for(i in bird.indices){
    
    obs.surv[i] ~ dbinom(p[i], n.min.potential.germ[tree.sp[i]])
    logit(p[i]) <- b.whole[tree.sp[i]] + b.tree.id[tree.sp[i],nested.tree.id[i]] + b.deinhib[tree.sp[i]] + b.scar.int + b.scar.bird.sp[bird.sp[i]] + b.scar.tree.sp[tree.sp[i]] + b.scar.bt.inter[bird.sp[i], tree.sp[i]] + b.bird.id[bird.sp[i],nested.bird.id[i]]
    
    }
    
    
    
    ### Priors
    
    # 1: tree id
    for(j in 1:n.tree.sp){
    for(k in 1:tree.id.nested.length[j]){
    b.tree.id[j,k] ~ dnorm(0, ran.tree.sp.tau[j])
    }
    ran.tree.sp.tau[j] <- pow(ran.tree.sp.sig[j], -2)
    ran.tree.sp.sig[j] ~ dunif(0,10)
    }
    
    
    
    
    # 2: bird id
    for(j in 1:n.bird.sp){
    for(k in 1:bird.id.nested.length[j]){
    b.bird.id[j,k] ~ dnorm(0, ran.bird.sp.tau[j])
    }
    ran.bird.sp.tau[j] <- pow(ran.bird.sp.sig[j], -2)
    ran.bird.sp.sig[j] ~ dunif(0,10)
    }
    
    
    
    # 3: whole fruit effect
    for(j in 1:n.tree.sp){
    b.whole[j] ~ dnorm(mu.whole,tau.whole)
    }
    mu.whole ~ dt(0, 1/2.5^2, 1)
    tau.whole <- pow(sig.whole, -2)
    sig.whole ~ dunif(0,10)
    
    
    # 4: deinhibition effect
    for(j in 1:n.tree.sp){
    b.deinhib[j] ~ dnorm(mu.deinhib,tau.deinhib)
    }
    mu.deinhib ~ dt(0, 1/2.5^2, 1)
    tau.deinhib <- pow(sig.deinhib, -2)
    sig.deinhib ~ dunif(0,10)
    
    
    # 5: scarification effects
    
    # Overall scar effect
    b.scar.int ~ dt(0, 1/2.5^2, 1)
    
    # Birds (bird species-level effect of bird sp id on germination)
    b.scar.bird.sp[1] <- 0
    for(j in 2:n.bird.sp){
    b.scar.bird.sp[j] ~ dnorm(0, tau.scar.bird.sp)
    }
    tau.scar.bird.sp <- pow(sig.scar.bird.sp, -2)
    sig.scar.bird.sp ~ dunif(0,10)
    
    # Trees (tree species-level effect of tree sp id on germination)
    b.scar.tree.sp[1] <- 0
    for(j in 2:n.tree.sp){
    b.scar.tree.sp[j] ~ dnorm(0, tau.scar.tree.sp)
    }
    tau.scar.tree.sp <- pow(sig.scar.tree.sp, -2)
    sig.scar.tree.sp ~ dunif(0,10)
    
    # Interactions (combination specific effect on germination)
    b.scar.bt.inter[bt.combos[1,1],bt.combos[1,2]] <- 0
    for(j in 2:n.combos){
    b.scar.bt.inter[bt.combos[j,1],bt.combos[j,2]] ~ dnorm(0, tau.scar.bt.inter)
    }
    tau.scar.bt.inter <- pow(sig.scar.bt.inter, -2)
    sig.scar.bt.inter ~ dunif(0,10)
    
    
    
    
    ### Derived quantities
    
    # 1: Fruit
    for(j in 1:n.tree.sp){
    out.fruit[j] <- ilogit(b.whole[j])
    }
    
    # 2: Mech
    for(j in 1:n.tree.sp){
    out.mech[j] <- ilogit(b.whole[j] + b.deinhib[j])
    }
    
    # 3: Gut passed
    for(j in 1:n.combos){
    out.gut[bt.combos[j,1],bt.combos[j,2]] <- ilogit(b.whole[bt.combos[j,2]] + b.deinhib[bt.combos[j,2]] + b.scar.int + b.scar.bird.sp[bt.combos[j,1]] + b.scar.tree.sp[bt.combos[j,2]] + b.scar.bt.inter[bt.combos[j,1], bt.combos[j,2]])
    }
    
    
    
    
    # Derived Ratios
    
    # Mech v. fruit
    for(j in 1:n.tree.sp){
    mech.v.fruit[j] <- out.mech[j] / out.fruit[j]
    }
    
    # Gut v. fruit
    for(j in 1:n.combos){
    gut.v.fruit[j] <- out.gut[bt.combos[j,1],bt.combos[j,2]] / out.fruit[bt.combos[j,2]]
    }
    
    # Gut v. mech
    for(j in 1:n.combos){
    gut.v.mech[j] <- out.gut[bt.combos[j,1],bt.combos[j,2]] / out.mech[bt.combos[j,2]]
    }
    
    
    } # End of model
    
    ",fill=TRUE)
sink()





# Run model and sample ---------------------------------------------------------

gp.model <- jags.model("gut germ model.txt", 
                       data = data, 
                       inits = inits,
                       n.chains = 3, n.adapt = 5000)

update(gp.model, n.iter = 10000)





# Save model output ------------------------------------------------------------

zj <- jags.samples(gp.model,data=data, inits=inits,
                   variable.names=c("mech.v.fruit",
                                    "gut.v.fruit",
                                    "gut.v.mech",
                                    "out.fruit",
                                    "out.mech",
                                    "out.gut"), 
                   n.iter=40000,
                   thin=100)





# Log response ratio approach: All treatments versus whole fruit ---------------

q975 <- function(x) quantile(x, 0.975, na.rm = T)
q025 <- function(x) quantile(x, 0.025, na.rm = T)

bt.combos <- data$bt.combos

# Set up for matrices that give quantiles for these ratios
gut.v.fruit.q025 <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
gut.v.fruit.med <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
gut.v.fruit.q975 <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
# Now fill in the real data from the jags output
gut.v.fruit.q025[bt.combos] <- summary(zj$gut.v.fruit, FUN=q025)$stat
gut.v.fruit.med[bt.combos] <- summary(zj$gut.v.fruit, FUN=median)$stat
gut.v.fruit.q975[bt.combos] <- summary(zj$gut.v.fruit, FUN=q975)$stat

# Fill in some vectors for the mech v fruit ratios
mech.v.fruit.q025 <- summary(zj$mech.v.fruit, FUN=q025)$stat
mech.v.fruit.med <- summary(zj$mech.v.fruit, FUN=median)$stat
mech.v.fruit.q975 <- summary(zj$mech.v.fruit, FUN=q975)$stat

# Set up for matrices that give quantiles for these ratios
gut.v.mech.q025 <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
gut.v.mech.med <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
gut.v.mech.q975 <- matrix(NA, nrow=max(bt.combos[,1]), ncol=max(bt.combos[,2]))
# Now fill in the real data from the jags output
gut.v.mech.q025[bt.combos] <- summary(zj$gut.v.mech, FUN=q025)$stat
gut.v.mech.med[bt.combos] <- summary(zj$gut.v.mech, FUN=median)$stat
gut.v.mech.q975[bt.combos] <- summary(zj$gut.v.mech, FUN=q975)$stat



# Single panel depulped versus whole fruit -------------------------------------

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# up to now have been using EBL legacy names - need to change this.
mf.order <- order(-mech.v.fruit.med)
mf.names <- firstup(tree.sp.list[mf.order])
mf.names <- gsub("Chili", "Capsicum", mf.names)
mf.names <- gsub("Pouteria", "Planchonella", mf.names)
mf.names <- gsub("Guamia", "Meiogyne", mf.names)

par(mfrow=c(1,1), pin = c(7,4))

plot(1:length(tree.sp.list), log10(mech.v.fruit.med[mf.order]),
     xaxt = "n",
     yaxt = "n",
     xlim = c(1,length(tree.sp.list)),
     pch = 16,
     ylim = log10(c(0.0625,8)),
     frame = F,
     xlab = "",
     ylab = "")
text(1:length(tree.sp.list)+0.5, log10(0.04), 
     mf.names,
     font = 3,
     xpd=TRUE, srt=45, pos=2)
mtext(side = 2, line = 4.5,
      "Deinhibition effect",
      cex = 1.2)

mtext(side = 2, line = 3.5,
      "(depulped : whole fruit germination)")

axis(1, at = 1:length(tree.sp.list), labels = rep("",length(tree.sp.list)))
axis(2, at = log10(c(0.0625,0.125, 0.25,0.5,1,2,4,8)), 
     labels = c("0.0625","0.125","0.25","0.5","1","2","4","8"),
     las = 2)
segments(x0 = 1:length(tree.sp.list),
         y0 = log10(mech.v.fruit.q025[mf.order]),
         y1 = log10(mech.v.fruit.q975[mf.order]),
         lend = "butt",
         lwd = 2)
abline(h=log10(1), lty = 2)



# Log response ratio approach: gut passed versus depulped ----------------------

xvals.spaced <- 1:4
counter <- 2
spacer <- 4

while(counter < 19){
  xvals.spaced <- c(xvals.spaced, max(xvals.spaced + spacer):max(xvals.spaced + spacer + 3))
  counter <- counter + 1
}


gm.order <- order(colMeans(gut.v.mech.med, na.rm = T), decreasing = T)
bars <- log10(c(as.vector(gut.v.mech.med[,gm.order])))
ci.lo <- log10(c(as.vector(gut.v.mech.q025[,gm.order])))
ci.hi <- log10(c(as.vector(gut.v.mech.q975[,gm.order])))

gm.names <- firstup(tree.sp.list[gm.order])
gm.names <- gsub("Chili", "Capsicum", gm.names)
gm.names <- gsub("Pouteria", "Planchonella", gm.names)
gm.names <- gsub("Guamia", "Meiogyne", gm.names)


bird.pch <- 21:24

par(pin=c(7,4), mfrow=c(1,1))
plot(xvals.spaced, bars, frame.plot = F,
     ylim = log10(c(0.16,4)),#c(min(c(0,ci.lo), na.rm = T), max(c(0,ci.hi), na.rm = T)),# 
     #log="y",
     pch = NA,
     bg = 1,
     cex = 0.8,
     #type = "h",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     las = 1)

mtext(side = 2, line = 4.5,
      "Scarification effect",
      cex = 1.2)

mtext(side = 2, line = 3.5,
      "(gut-passed : depulped germination)")

xvals.spaced.species <- seq(2.5,mean(xvals.spaced[(length(xvals.spaced)-3):(length(xvals.spaced)-2)]),length.out = length(tree.sp.list))

text(xvals.spaced.species - 2.5,
     log10(5), 
     gm.names, 
     cex = 0.8, font = 3, srt = 45,
     pos = 4,
     xpd = T)

segments(x0 = c(1-spacer/2, xvals.spaced.species + spacer),
         y0 = log10(0.25)-0.08,
         y1 = log10(4)+0.08,
         lwd = 0.5,
         xpd = T,
         col = "grey70")
abline(h = 0, lty = 2)

legend(0.5,log10(0.5), 
       pch=bird.pch,
       legend=c("starling","fruit dove", "g. white-eye", "b. white-eye"))

segments(x0 = xvals.spaced, y0 = ci.lo, y1 = ci.hi,
         lwd = 1.5,
         lend = "butt")
points(xvals.spaced, bars,
       pch = bird.pch,
       bg = "white",
       cex = 0.9)
#axis(1, at=c(xvals.spaced.species), labels=rep("", length(xvals.spaced.species)))
axis(2, at=log10(c(0.25, 0.5, 1, 2, 4)), labels = c(0.25, 0.5, 1, 2, 4), las = 1)



# Log response ratio approach: gut passed versus whole -------------------------

xvals.spaced <- 1:4
counter <- 2
spacer <- 4

while(counter < 19){
  xvals.spaced <- c(xvals.spaced, max(xvals.spaced + spacer):max(xvals.spaced + spacer + 3))
  counter <- counter + 1
}


gf.order <- order(colMeans(gut.v.fruit.med, na.rm = T), decreasing = T)
bars <- log10(c(as.vector(gut.v.fruit.med[,gf.order])))
ci.lo <- log10(c(as.vector(gut.v.fruit.q025[,gf.order])))
ci.hi <- log10(c(as.vector(gut.v.fruit.q975[,gf.order])))


bird.pch <- 21:24

par(pin=c(7,4), mfrow=c(1,1))
plot(xvals.spaced, bars, frame.plot = F,
     ylim = log10(c(0.0625,16)),
     #log="y",
     pch = NA,
     bg = 1,
     cex = 0.8,
     #type = "h",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     las = 1)

mtext(side = 2, line = 4.5,
      "Gut passage effect",
      cex = 1.2)

mtext(side = 2, line = 3.5,
      "(gut-passed : whole germination)")

xvals.spaced.species <- seq(2.5,mean(xvals.spaced[(length(xvals.spaced)-3):(length(xvals.spaced)-2)]),length.out = length(tree.sp.list))

gf.names <- firstup(tree.sp.list[gf.order])
gf.names <- gsub("Chili", "Capsicum", gf.names)
gf.names <- gsub("Pouteria", "Planchonella", gf.names)
gf.names <- gsub("Guamia", "Meiogyne", gf.names)


text(xvals.spaced.species - 2.5,
     log10(12), 
     gf.names, 
     cex = 0.8, font = 3, srt = 45,
     pos = 4,
     xpd = T)

segments(x0 = c(1-spacer/2, xvals.spaced.species + spacer),
         y0 = log10(0.05),
         y1 = log10(12),
         lwd = 0.5,
         xpd = T,
         col = "grey70")
abline(h = 0, lty = 2)

legend(0.5,log10(0.5), 
       pch=bird.pch,
       legend=c("starling","fruit dove", "g. white-eye", "b. white-eye"))

segments(x0 = xvals.spaced, y0 = ci.lo, y1 = ci.hi,
         lwd = 1.5,
         lend = "butt")
points(xvals.spaced, bars,
       pch = bird.pch,
       bg = "white",
       cex = 0.9)
#axis(1, at=c(xvals.spaced.species), labels=rep("", length(xvals.spaced.species)))
axis(2, at = log10(c(0.0625,0.125, 0.25,0.5,1,2,4,8,16)), 
     labels = c("0.0625","0.125","0.25","0.5","1","2","4","8","16"),
     las = 2)





par(mfrow=c(1,1))
gf.mat <- t(gut.v.fruit.med)
gf.mat[which(is.na(gf.mat))] <- 0
gf.mat <- as.data.frame(gf.mat)
colnames(gf.mat) <- bird.sp.list
rownames(gf.mat) <- tree.sp.list


# This is ugly code - there must be a better way...

rotate <- function(x) t(apply(x, 2, rev))

wtgd.destroyed.spp <- c("pipturus","aidia","ficus t.",
                        "ficus p.",  "psychotria", "premna",
                        "passiflora", "carica", "pouteria",
                        "momordica", "guamia", "coccinia")
wtgd.row <- ifelse(tree.sp.list %in% wtgd.destroyed.spp, 0.01, NA)

gut.v.fruit.med.with.wtgd <- matrix(rbind(gut.v.fruit.med, wtgd.row),nrow=5)

myvec <- round(log10(gut.v.fruit.med.with.wtgd)*100)
myvec
myvec <- sortweb(myvec)
plant.sp.order <- order(-colSums(ifelse(is.na(myvec),0,1)))
myvec <- myvec[c(1,2,5,3,4),plant.sp.order]
myvec <- rotate(rotate(myvec))
myvec <- apply(myvec,2,rev)

par(pin=c(5/3.2,length(tree.sp.list)/3.2))
image(myvec,
      col = smoothColors("red",
                         200,#-min(myvec, na.rm = T)-1,
                         "lightgrey",
                         max(myvec, na.rm = T)-1,
                         "blue"),
      frame = F,
      xaxt = "n",
      yaxt = "n",
      breaks=-201:(max(myvec, na.rm=T)+1))#min(myvec, na.rm = T):(max(myvec, na.rm=T)+1))
text(seq(from = 0,to = 1, length = 5)-0.125, 1.05, c("starling", "fruit dove", "ground dove", "g. white-eye", "b. white-eye"), xpd=TRUE, srt=35, pos=4)

net.plants <- firstup(tree.sp.list[rev(plant.sp.order)])
net.plants <- gsub("Chili", "Capsicum", net.plants)
net.plants <- gsub("Pouteria", "Planchonella", net.plants)
net.plants <- gsub("Guamia", "Meiogyne", net.plants)
text(-0.125,seq(from = 0,to = 1, length = length(tree.sp.list)), net.plants, xpd=TRUE, pos=2, font = 3)


legend(x=0.5,
       y=0.22,
       bty="n",
       pch = c(15,15,15,15),
       pt.cex = 2,
       legend = rev(c("<0.01",0.5,1,2,5)),
       col = rev(c("red",rgb(1,0,0,0.5),"lightgrey",rgb(0,0,1,0.5),"blue")),
       xpd=T)
text(0.4, 0.24, xpd = T, pos = 4,
     "gut-passed : fruit\ngermination ratio")