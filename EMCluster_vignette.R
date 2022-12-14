library(EMCluster)
library(PPtree)
library(MASS)

##########################################################
# EM Cluster
##########################################################
rm(list = ls(all.names = TRUE))

x <- as.matrix(crabs[, 4:8])
p <- ncol(x)
min.n <- p * (p + 1) / 2
min.n.iter <- 10
.EMC$short.iter <- 5000

ret <- list()
ret.save <- list()
k <- 1
i.rep <- 0
llhdval.curr <- -Inf
K.max <- 6
repeat{
  # Initial.
  i.rep <- i.rep + 1
  seed <- 1234 + k + i.rep
  set.seed(seed)
  
  ret[[k]] <- list()
  ret[[k]]$seed <- seed
  
  # Run RndEM.
  method <- "Rnd.EM"
  repeat{
    tmp.init <- init.EM(x, nclass = k, min.n = min.n, min.n.iter = min.n.iter,
                        method = method)
    if(any(tmp.init$nc <= min.n)){
      next
    }
    var <- LTSigma2variance(tmp.init$LTSigma)
    flag <- 0
    for(k.var in 1:dim(var)[3]){
      # Check degenerate.
      tmp <- try(solve(var[,, k.var]), silent = TRUE)
      if(class(tmp)[1] == "try-error"){
        flag <- 1
        break
      }
    }
    if(flag == 0){
      ret[[k]]$Rnd <- tmp.init
      break
    }
  }
  
  # Check llhdval.
  flag <- 1
  if(ret[[k]]$Rnd$llhdval > llhdval.curr){
    llhdval.curr <- ret[[k]]$Rnd$llhdval
    ret.save[[k]] <- ret[[k]]$Rnd
    ret.save[[k]]$case.save <- "Rnd"
    flag <- 0
  } 
  if(flag == 0){
    cat("k = ", k, ", llhdval = ", llhdval.curr, ", ", sep = "")
    cat("nc =", ret.save[[k]]$nc, "\n")
    k <- k + 1
  }
  
  if(k > K.max){
    break
  }
}
save(list = c("x", "ret", "ret.save"), file = "C:/Users/user/Desktop/UNIST_internship/data/crabs.RndEM.rda")
####################################################
# adjR
##########################################################
rm(list = ls())

load("C:/Users/user/Desktop/UNIST_internship/data/crabs.RndEM.rda")

crabs.id <- rep(1L, nrow(crabs))
crabs.id[crabs$sp == "B" & crabs$sex == "F"] <- 2L
crabs.id[crabs$sp == "O" & crabs$sex == "M"] <- 3L
crabs.id[crabs$sp == "O" & crabs$sex == "F"] <- 4L

for(k0 in 1:6){
  emobj0 <- ret.save[[k0]]
  
  cat("K = ", k0, "\n", sep = "")
  print(RRand(emobj0$class, crabs.id))
}

####################################################
# test
##########################################################
rm(list = ls())

load("C:/Users/user/Desktop/UNIST_internship/data/crabs.RndEM.rda")

tau <- 0.5
n.mc.base <- 1000

ret.test <- NULL
for(k0 in 1:5){
  for(ka in (k0+1):6){
    seed <- 1234 + k0 + ka
    set.seed(seed)
    n.mc <- n.mc.base * max(c(k0, ka))
    
    emobj0 <- ret.save[[k0]]
    emobja <- ret.save[[ka]]
    
    ll.0 <- emobj0$llhdval
    ll.a <- emobja$llhdval
    delta.hat <- emobja$llhdval - emobj0$llhdval
    E.delta <- get.E.delta(x, emobj0, emobja, tau = tau, n.mc = n.mc)
    E.chi2.0 <- get.E.chi2(x, emobj0, emobja, "0", tau = tau, n.mc = n.mc,
                           verbose = FALSE)
    E.chi2.a <- get.E.chi2(x, emobj0, emobja, "a", tau = tau, n.mc = n.mc,
                           verbose = FALSE)
    T <- 2 * (delta.hat - E.delta)
    
    pv0 <- pchisq.my(T, E.chi2.0[1], E.chi2.0[2], lower.tail = FALSE)
    pva <- pchisq.my(T, E.chi2.a[1], E.chi2.a[2], lower.tail = FALSE)
    pv <- pv0 * tau + pva * (1 - tau)
    
    cat("  H0: K = ", k0, " vs Ha: K = ", ka, "\n",
        "    ll0 = ", ll.0, ", lla = ", ll.a, "\n",
        "    df0 = ", E.chi2.0[1], ", nc0 = ", E.chi2.0[2],
        ", dfa = ", E.chi2.a[1], ", nca = ", E.chi2.a[2], "\n",
        "    delta.hat = ", delta.hat, ", E.delta = ", E.delta,
        ", T = ", T, "\n",
        "    pv0 = ", pv0, ", pva = ", pva, " pv = ", pv,
        "\n", sep = "")
    
    ret.test <- rbind(ret.test,
                      c(k0, ka, ll.0, ll.a, E.chi2.0, E.chi2.a, delta.hat,
                        E.delta, T, pv0, pva, pv, seed))
  }
}

colnames(ret.test) <- c("K0", "Ka", "ll0", "lla",
                        "df0", "nc0", "dfa", "nca",
                        "delta.hat", "E.delta", "T",
                        "pv0", "pva", "pv", "seed")

save(list = c("ret.test"), file = "C:/Users/user/Desktop/UNIST_internship/data/crabs.test.rda")

print(ret.test)
##########################################################
# plot
##########################################################
rm(list = ls())

load("C:/Users/user/Desktop/UNIST_internship/data/crabs.RndEM.rda")

x <- as.matrix(crabs[, 4:8])

da.s.all <- list()
for(k0 in 2:5){
  emobj <- ret.save[[k0]]
  var <- LTSigma2variance(emobj$LTSigma)
  x.pp <- x
  Sigma.pp <- var
  for(k.var in 1:dim(var)[3]){
    tmp <- eigen(var[,, k.var])
    Sigma.k.inv <- tmp$vector %*% diag(sqrt(1/tmp$values)) %*% t(tmp$vector)
    
    tmp <- x[emobj$class == k.var,]
    tmp.mu <- colMeans(tmp)
    tmp <- t(t(tmp) - tmp.mu) %*% Sigma.k.inv
    tmp <- t(t(tmp) + as.vector(tmp.mu))
    x.pp[emobj$class == k.var,] <- tmp
    
    # Sigma.pp[,, k.var] <- Sigma.k.inv
  }
  set.seed(1234)
  tmp.pp <- PP.optimize.random("LDA", 2, x.pp, emobj$class, std = FALSE)
  # plot(y %*% tmp.pp$proj.best, xlab = "PP1", ylab = "PP2",
  #      col = ret.save[[k0]]$class, pch = ret.save[[k0]]$class)
  x.new <- x %*% tmp.pp$proj.best
  mu.new <- emobj$Mu %*% tmp.pp$proj.best
  var.new <- array(0, dim = c(2, 2, dim(var)[3]))
  for(k.var in 1:dim(var)[3]){
    tmp <- t(tmp.pp$proj.best) %*% Sigma.pp[,, k.var] %*% tmp.pp$proj.best
    var.new[,, k.var] <- tmp
  }
  da.s <- list(pi = emobj$pi,
               Mu = mu.new,
               LTSigma = variance2LTSigma(var.new),
               class = emobj$class,
               nclass = emobj$nclass)
  da.s.all[[k0]] <- list(model = da.s, x = x.new)
  
  # color.class <- 1:11
  # postscript(file = paste("./plot/iris_pp_k=", k0, ".eps", sep = ""),
  #            height = 6, width = 6, horizontal = FALSE)
  # plotem(da.s, x.new, xlab = "PP1", ylab = "PP2",
  #        main = paste("K = ", k0, sep = ""))
  # dev.off()
}

save(da.s.all, x, file = "C:/Users/user/Desktop/UNIST_internship/data/crabs.pp.rda")

##########################################################
# plot contour
##########################################################
rm(list = ls())

load("C:/Users/user/Desktop/UNIST_internship/data/crabs.pp.rda")

for(k0 in 2:5){
  Pi <- da.s.all[[k0]]$model$pi
  Mu <- da.s.all[[k0]]$model$Mu
  S <- LTSigma2variance(da.s.all[[k0]]$model$LTSigma)
  class <- da.s.all[[k0]]$model$class
  da <- da.s.all[[k0]]$x
  
  # crabs has weird projection and needs a rotation to get a better view.
  if(k0 == 2) angle <- 40 / 180 * pi
  if(k0 == 3) angle <- 10 / 180 * pi
  if(k0 == 4) angle <- 45 / 180 * pi
  if(k0 == 5) angle <- -20 / 180 * pi
  if(k0 == 6) angle <- 0 / 180 * pi
  
  pdf(paste("C:/Users/user/Desktop/UNIST_internship/plot/crabs_ppcont_k=", k0, ".pdf", sep = ""),
      width = 6, height = 6)
  plotppcontour(da, Pi, Mu, S, class)
  dev.off()
}

##########################################################
# qmap
##########################################################
rm(list = ls())

load("C:/Users/user/Desktop/UNIST_internship/data/crabs.test.rda")

### quant.av
K.min <- min(ret.test[, 1]) 
K.max <- max(ret.test[, 2]) 
pv.un <- ret.test[, 14]
pv.un <- pv.un[ret.test[, 1] <= K.max & ret.test[, 2] <= K.max]

#postscript("./plot/qmap.iris.ps", width = 5, height = 5, horizontal = FALSE)
pdf("./plot/qmap.crabs.pdf", width = 5, height = 5)
plotq(pv.un, dim = K.max - 1, start = K.min)
dev.off()
