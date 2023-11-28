# require(ggmap)
require(ggplot2)
require(wesanderson)
require(LaplacesDemon)
require(fields)
require(this.path)

RNGkind(sample.kind = "Rejection")
rm(list = ls())
setwd(this.path::here())
pic.format = ".jpeg"

# Some hyperparameters to specify the dataset and the fitted model
dat = "NepalExtended" # "Simulated" # 
pen.type = "L2"
if (dat == "Simulated") {cv = F} else if (dat == "NepalExtended") {cv = T}
seed = 1
model = "MSP"


# Load dataset and the fitted model
load(paste0(dat, ".Rdata"))
load(paste0("Modeling/", dat, "_", model, "_", pen.type, 
            "_seed", seed, "_", cv, "_Results.Rdata"))
# if (dat == "NepalExtended") load(paste0("NepalMap.Rdata"))


# Load functions
source('Functions/Utils.R')
source('Functions/Fit.R')
source('Functions/Objectives.R')
# source('Functions/Lambda_tuning.R')
# source('Functions/Merge_subr.R')
# source('Functions/Algorithm1.R')


# Prepare some variables
attach(models.save)

# Modeling time
run.time.base = model.base$time
run.time.merge = model.merge$RT
cat(" Modeling time for the base model (mainly for the lambda-tuning):", run.time.base/60, "mins\n",
    "Modeling time for the merged model:", run.time.merge/60, "mins\n")

attach(hyper)
attach(model.merge)

set.seed(seed)
# dependence parameters estimates
est.base = est.ns.traj[[1]]
est.merge = est.ns.traj[[length(est.ns.traj)]]
# partitions
reg.idx.base = reg.idx.traj[[1]]
reg.idx.merge = reg.idx.traj[[length(reg.idx.traj)]]
# number of regions
R.base = R.traj[[1]]
R.merge = R.traj[[length(R.traj)]]
# lambda
lambda.base = lambda.traj[[1]]
lambda.merge = lambda.traj[[length(lambda.traj)]]
# adjacent
adj.set.base = adj.set.traj[[1]]
adj.set.merge = adj.set.traj[[length(adj.set.traj)]]

# Prepare dataframe for plot
est.base.mat = vec2mat(est.base, R.base, lambda.base)
est.merge.mat = vec2mat(est.merge, R.merge, lambda.merge)
if (dat == "NepalExtended") {
  set.seed(4)
  clu <- kmeans(logit(elev/(max(elev)+1)), 3, nstart = 25, iter.max = 100000) #, max(elev)+1
  elev.clu = clu$cluster
  sapply(1:3, function(i) mean(elev[which(clu$cluster == i)]))
  df = data.frame(s1 = S[,1], s2 = S[,2], 
                  elev = elev, elev.clu = elev.clu,
                  parti.base = reg.idx.base,
                  parti.merge = reg.idx.merge,
                  sig.base = est.base.mat[1,][reg.idx.base], 
                  phi.base = est.base.mat[2,][reg.idx.base],
                  sig.merge = est.merge.mat[1,][reg.idx.merge],
                  phi.merge = est.merge.mat[2,][reg.idx.merge])
  width1 = (range(S[,1])[2] - range(S[,1])[1])*2-7; height1 = range(S[,2])[2] - range(S[,2])[1]
  width2 = 14; height2 = 14
  wid1 = 7
  sill.limits = c(2, 8)
  range.limits = c(3, 9)
} else if (dat == "Simulated") {
  df = data.frame(s1 = S[,1], s2 = S[,2], 
                  parti.base = reg.idx.base,
                  parti.merge = reg.idx.merge,
                  sig.base = est.base.mat[1,][reg.idx.base], 
                  phi.base = est.base.mat[2,][reg.idx.base],
                  sig.merge = est.merge.mat[1,][reg.idx.merge],
                  phi.merge = est.merge.mat[2,][reg.idx.merge])
  width1 = 11; height1 = 5
  width2 = 14; height2 = 14
  wid1 = 5
  sill.limits = c(min(c(df$sig.base, df$sig.merge)), max(c(df$sig.base, df$sig.merge)))
  range.limits = c(min(c(df$phi.base, df$phi.merge)), max(c(df$phi.base, df$phi.merge)))
}


# colors
mycolor0 = sample(colors(), R.base)
mycolor1 = c('yellow', 'goldenrod1', 'darkorange', 'orangered', 'red2', 
             'darkolivegreen', 'darkgreen', 'green4', 'green', 'cornflowerblue', 
             'chocolate4', 'darkred', 'seagreen1', 'deepskyblue', 'deepskyblue4', 
             'blue', 'navy', 'purple4', 'grey32', 'grey74', 'plum1', 'hotpink', 
             'deeppink1', 'darkorchid1', 'black')
mycolor2 = c('#9fa887', '#e5dbbc', '#ab742c')[3:1]

####################################################################################
# Plot: General Plots for Data Application #########################################
####################################################################################
# This part will be skipped if the dataset is simulated, i.e., dat == "Simulated" 

if (dat == "NepalExtended") {
  p.elev <- ggplot(data=df, aes(x = s1, y = s2, fill = elev), width = 0.25, height = 0.25) +
    geom_tile() +
    scale_fill_gradientn("Elevation (m)", colours = terrain.colors(10)) +
    theme(plot.title = element_text(hjust = 0.5, size=25),
          legend.key.size = unit(0.1, "in"), axis.title=element_text(size=18),
          axis.text=element_text(size=14), legend.text = element_text(size=12),
          legend.title = element_text(size=12), legend.position = c(0.85, 0.75)) +
    xlab("Longitude") + ylab("Latitude") 
    
  print(p.elev)
  
  p.clu = ggplot(data=df, aes(x = s1, y = s2, colour = as.factor(elev.clu))) +
    geom_point(shape = 15, size = 2) +
    theme(legend.key.size = unit(0.1, "in"), axis.title=element_text(size=18),
          axis.text=element_text(size=14), legend.text = element_text(size=12),
          legend.title = element_text(size=12), legend.position = c(0.85, 0.75)) +
    xlab("Longitude") + ylab("Latitude") +
    scale_colour_manual("Subdomain", labels=c("high","medium", "low"),
                        values=mycolor2) 
  print(p.clu)
  
  p.general = ggpubr::ggarrange(p.elev, p.clu, nrow = 1)
  p.general
  ggsave(paste0("pic/", dat, "_elev", pic.format),
         plot = p.general, width = width1, height = height1, units = "cm")
}


####################################################################################
# Plot: Partitions #################################################################
####################################################################################

if (dat == "Simulated") {
  true.parti = sapply(1:nrow(partition), function(i) {
    parti.i = as.numeric(partition[i,]-1)
    parti.i[1] + parti.i[2]*2 + 1
  })
  p.true = ggplot(data=df, aes(x = s1, y = s2, colour = as.factor(true.parti))) +
    geom_point(shape = 15, size = 2) +
    scale_colour_manual(values=mycolor1) +
    theme(legend.position="none", axis.title=element_text(size=18),
          axis.text=element_text(size=10)) +
    xlab("Longitude") + ylab("Latitude")
  ggsave(paste0("pic/", dat, "_trueparti_", model, pic.format),
         plot = p.true, width = width1/2, height = height1, units = "cm")
}

p.basepar = ggplot(data=df, aes(x = s1, y = s2, colour = as.factor(parti.base))) +
  geom_point(shape = 15, size = 2) +
  scale_colour_manual(values=mycolor0) +
  theme(legend.position="none", axis.title=element_text(size=18),
        axis.text=element_text(size=10)) +
  xlab("Longitude") + ylab("Latitude")

p.mergepar = ggplot(data=df, aes(x = s1, y = s2, colour = as.factor(parti.merge))) +
  geom_point(shape = 15, size = 2) +
  scale_colour_manual(values=mycolor1) +
  theme(legend.position="none", axis.title=element_text(size=18),
        axis.text=element_text(size=10)) +
  xlab("Longitude") + ylab("Latitude")

p.parti = ggpubr::ggarrange(p.basepar, p.mergepar, nrow = 1)
p.parti
ggsave(paste0("pic/", dat, "_parti_", model, "_", pen.type, pic.format),
       plot = p.parti, width = width1, height = height1, units = "cm")


####################################################################################
# Plot: Dependence Parameter Estimates #############################################
####################################################################################

axis.text.size = 8
axis.title.size = 8
legend.size = 8
legend.key = 0.2

p.sig.base = ggplot(df, mapping = aes(x = s1, y = s2, fill = sig.base), 
                    width = 0.25, height = 0.25) + geom_tile() +
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis(option = "viridis", limits = sill.limits, oob = scales::squish,
                     name="Sill Est") +
  theme(plot.title = element_text(hjust = 0.5, size=25),
        legend.key.size = unit(legend.key, "in"), legend.text = element_text(size=legend.size),
        legend.title = element_text(size=legend.size), legend.position = "none",
        axis.title=element_text(size=axis.title.size),
        axis.text=element_text(size=axis.text.size))

p.phi.base = ggplot(df, mapping = aes(x = s1, y = s2, fill = phi.base), 
                    width = 0.25, height = 0.25) + geom_tile() +
  xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5, size=25),
        axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  scale_fill_viridis(option = "viridis", limits = range.limits, oob = scales::squish,
                     name="Range Est") +
  theme(legend.key.size = unit(legend.key, "in"), legend.text = element_text(size=legend.size),
        legend.title = element_text(size=legend.size), legend.position = "none",
        axis.title=element_text(size=axis.title.size),
        axis.text=element_text(size=axis.text.size))

p.sig.merge = ggplot(df, mapping = aes(x = s1, y = s2, fill = sig.merge), 
                     width = 0.25, height = 0.25) + geom_tile() +
  xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5, size=25),
        axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  scale_fill_viridis(option = "viridis", limits = sill.limits, oob = scales::squish,
                     name="Sill Est") +
  theme(legend.key.size = unit(legend.key, "in"), legend.text = element_text(size=legend.size),
        legend.title = element_text(size=legend.size), legend.position = "none",
        axis.title=element_text(size=axis.title.size),
        axis.text=element_text(size=axis.text.size))

p.phi.merge = ggplot(df, mapping = aes(x = s1, y = s2, fill = phi.merge), 
                     width = 0.25, height = 0.25) + geom_tile() +
  xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5, size=25),
        axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  scale_fill_viridis(option = "viridis", limits = range.limits, oob = scales::squish,
                     name="Range Est") +
  theme(legend.key.size = unit(legend.key, "in"), legend.text = element_text(size=legend.size),
        legend.title = element_text(size=legend.size), legend.position = "none",
        axis.title=element_text(size=axis.title.size),
        axis.text=element_text(size=axis.text.size))

# Joint plot ###################################################################
library(gridExtra)
library(grid)

plots = list(p.sig.base, p.phi.base, p.sig.merge, p.phi.merge)
gSig <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
legend1 <- grobTree(gSig[[which(sapply(gSig, function(x) x$name) == "guide-box")]])
titles1 = list(paste0("Model B ", pen.type), paste0("Model M ", pen.type))
gPhi <- ggplotGrob(plots[[2]] + theme(legend.position="bottom"))$grobs
legend2 <- grobTree(gPhi[[which(sapply(gPhi, function(x) x$name) == "guide-box")]])
titles2 = list(paste0("Model B ", pen.type), paste0("Model M ", pen.type))
titles0 = list("Sill", "Range")
titles1 <- lapply(1:2, function(ii) textGrob(paste0(titles1[[ii]])))
titles2 <- lapply(1:2, function(ii) textGrob(paste0(titles2[[ii]])))
titles0 <- lapply(1:2, function(ii) textGrob(paste0(titles0[[ii]])))

hei1 <- unit(5, "cm")
hei2 <- unit(2, "cm")
hei3 <- unit(0.1, "cm")
hei4 <- unit(0.8, "cm")
wid = unit(wid1, "cm")

layout <- matrix(c(11, 12,
                   1, 2,
                   7, 9, 
                   3, 4,
                   8, 10,
                   5, 6), ncol=2, byrow=TRUE)
p.est = grid.arrange(grobs = c(plots, list(legend1), list(legend2), titles1, titles2, titles0),
                     layout_matrix = layout,
                     heights = unit.c(hei4, hei1, hei3, hei1, hei3, hei2),
                     widths = unit.c(wid, wid))
ggsave(paste0("pic/", dat, "_est_", model, "_", pen.type, pic.format),
       plot = p.est, width = width2, height = height2, units = "cm")

####################################################################################
# PL, PPL ##########################################################################
####################################################################################
if (model == "MSP") {
  pl.s = model.s$pl.s
  pl.base = -PL.ns(log(est.base), sets$D.valid, sets$Z.valid, 
                   reg.idx.base[sets$validIdx], "None", 
                   adj.set.base, lambda.base, model = model)
  pl.merge = -PL.ns(log(est.merge), sets$D.valid, sets$Z.valid, 
                    reg.idx.merge[sets$validIdx], "None", 
                    adj.set.merge, lambda.merge, model = model)
  ppl.base = ppl.ns.traj[[1]]
  ppl.merged = ppl.ns.traj[[length(ppl.ns.traj)]]
} else if (model == "IMSP") {
  dat.valid = -1/log(1-exp(-sets$Z.valid))
  pl.s = -PL.s(log(model.s$est.s), sets$D.valid, dat.valid, model = "IMSP_trans")
  pl.base = -PL.ns(log(est.base), sets$D.valid, dat.valid, 
                   reg.idx.base[sets$validIdx], "None", 
                   adj.set.base, lambda.base, model = "IMSP_trans")
  pl.merge = -PL.ns(log(est.merge), sets$D.valid, dat.valid, 
                    reg.idx.merge[sets$validIdx], "None", 
                    adj.set.merge, lambda.merge, model = "IMSP_trans")
  fold_n = hyper$fold_n
  ppl.base = mean(sapply(1:fold_n, function(i){
    dat.hold = -1/log(1-exp(-sets$Z.hold[[i]]))
    -PL.ns(log(est.base), sets$D.hold[[i]],
           dat.hold, reg.idx.base[sets$holdIdx[[i]]], pen.type, adj.set.base,
           lambda.base, model = "IMSP_trans")
  }))
  ppl.merged = mean(sapply(1:fold_n, function(i){
    dat.hold = -1/log(1-exp(-sets$Z.hold[[i]]))
    -PL.ns(log(est.merge), sets$D.hold[[i]],
           dat.hold, reg.idx.merge[sets$holdIdx[[i]]], pen.type, adj.set.merge,
           lambda.merge, model = "IMSP_trans")
  }))
}

df.summary1 = data.frame(Model_S = c(round(pl.s), "--"),
                        Model_B = round(c(pl.base, ppl.base)),
                        Model_M = round(c(pl.merge, ppl.merged)))
rownames(df.summary1) = c("PL", "PPL")
df.summary1


if (dat == "Simulated") {
  RMSE_Sill.sta = sqrt(mean((model.s$est.s[1] - sig.all)^2))
  RMSE_Sill.base = sqrt(mean((df$sig.base - sig.all)^2))
  RMSE_Sill.merge = sqrt(mean((df$sig.merge - sig.all)^2))
  RMSE_Range.sta = sqrt(mean((model.s$est.s[2] - sig.all)^2))
  RMSE_Range.base = sqrt(mean((df$phi.base - phi.all)^2))
  RMSE_Range.merge = sqrt(mean((df$phi.merge - phi.all)^2))
  
  df.summary2 = data.frame(Model_S = round(c(RMSE_Sill.sta, RMSE_Range.sta), 2),
                           Model_B = round(c(RMSE_Sill.base, RMSE_Range.base), 2),
                           Model_M = round(c(RMSE_Sill.merge, RMSE_Range.merge), 2))
  rownames(df.summary2) = c("RMSE_Sill", "RMSE_Range")
  df.summary2
}

