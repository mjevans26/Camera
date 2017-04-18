EN.test.lam <- spNmix.lam(EN14$y, EN14$X, M = 1000, niters = 50000, xlims = EN14$xlims, ylims = EN14$ylims, tune=c(10, 1, 2),cov=EN14$forest)
EN.test.sig <- spNmix.sig(EN14$y, EN14$X, M = 1000, niters = 50000, xlims = EN14$xlims, ylims = EN14$ylims, tune=c(0.2,10,0.15,5), cov=EN14$forest)
EN.10.15.2 <- spNmix(EN14$y, EN14$X, M = 1000, niters = 50000, xlims = EN14$xlims, ylims = EN14$ylims, tune = c(10, 15, 2))

