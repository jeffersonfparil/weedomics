library(testthat)
library(glmnetUtils)
library(raster)
library(parallel)
library(ExtDist)

dirname_src = utils::getSrcDirectory(function(){})[1]
if (is.na(dirname_src)) {
    dirname_src = "."
}
setwd(dirname_src)
source("./data_loading_and_merging.r")

PERF = function(y_test, y_hat) {
    n = length(y_test)
    if ((var(y_hat)==0) | (var(y_test)==0)) {
        cor = 0.0
    } else {
        cor = cor(y_test, y_hat)
    }
    e = y_test - y_hat
    mbe = ( sum(e)/n )
    mae = ( sum(abs(e))/n )
    mse = ( sum(e^2)/n )
    rmse = ( sqrt(sum(e^2)/n) )
    return(list(cor=cor, mbe=mbe, mae=mae, mse=mse, rmse=rmse, y_hat=y_hat, y_test=y_test))
}

OLS_train_test = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    b_hat = tryCatch(t(x_train) %*% solve(x_train %*% t(x_train)) %*% y_train, error=function(e){
                    SVD = svd(x_train %*% t(x_train))
                    vec_sv = rep(0, length(SVD$d))
                    for (i in 1:length(vec_sv)) {
                        if (SVD$d[i] <= .Machine$double.eps)  {
                            vec_sv[i] = 0.0
                        } else {
                            vec_sv[i] = 1/SVD$d[i]
                        }
                    }
                    return(t(x_train) %*% SVD$v %*% diag(vec_sv) %*% t(SVD$u) %*% y_train)
                })
    y_hat = (x_test %*% b_hat)[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

LASSO_train_test = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    nfolds = 10
    mod_lasso = cva.glmnet(
        x_train,
        y_train,
        alpha = 1.0,
        nfolds = nfolds,
        foldid = sample(rep(seq_len(nfolds), length = nrow(x_train))),
        outerParallel = NULL,
        checkInnerParallel = TRUE
        )
    y_hat = predict(mod_lasso$modlist[[1]], newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

RIDGE_train_test = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    nfolds = 10
    mod_lasso = cva.glmnet(
        x_train,
        y_train,
        alpha = 0.0,
        nfolds = nfolds,
        foldid = sample(rep(seq_len(nfolds), length = nrow(x_train))),
        outerParallel = NULL,
        checkInnerParallel = TRUE
        )
    y_hat = predict(mod_lasso$modlist[[1]], newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

ELASTIC_train_test = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    cores = parallel::makeCluster(detectCores())
    nfolds = 10
    vec_alphas = sort(unique(c(seq(0, 1, len=11)^3, 1-seq(0, 1, len=11)^3)))
    mod_elastic = cva.glmnet(
        x_train,
        y_train,
        alpha = vec_alphas,
        nfolds = nfolds,
        foldid = sample(rep(seq_len(nfolds), length = nrow(x_train))),
        outerParallel = cores,
        checkInnerParallel = TRUE
        )
    vec_mse = c()
    for (i in (1:length(vec_alphas))) {
        # i = 1
        mod = mod_elastic$modlist[[i]]
        y_hat = predict(mod, newx=x_test, s="lambda.min")[,1]
        vec_mse = c(vec_mse, PERF(y_test, y_hat)$mse)
    }
    cbind(vec_alphas, vec_mse)
    idx = which(vec_mse==min(vec_mse, na.rm=TRUE))
    y_hat = predict(mod_elastic$modlist[[idx]], newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

KFOLD_CV = function(x, y, r=5, k=10) {
    x = as.matrix(x)
    y = as.matrix(y)
    n = length(y)
    s = floor(n/k)
    if (s < 5) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("The number of cross-validation folds is too big to have a fold size of at least 5.")
        print("We expect to have at least 5 expected reponses to get good estimates of prediction performance per cross-validation fold.")
        print("Please reduce the number of folds.")
        return("ERROR!")
    }
    list_rep_fold = list()
    for (rep in 1:r) {
        for (fold in 1:k) {
            eval(parse(text=paste0("list_rep_fold$`", rep, "-", fold, "` = c(rep, fold)")))
        }
    }
    perf = parallel::mclapply(list_rep_fold, function(rep_fold) {
                                # rep_fold = c(1, 2)
                                rep = rep_fold[1]
                                fold = rep_fold[2]
                                set.seed(rep)
                                idx = sample(c(1:n), n, replace=FALSE)
                                i = (fold-1)*s + 1
                                j = fold*s
                                if (fold == k) {
                                    j = n
                                }
                                bool_test = c(1:n) %in% c(i:j)
                                bool_train = !bool_test
                                idx_train = idx[bool_train]
                                idx_test = idx[bool_test]
                                ols = tryCatch(OLS_train_test(x, y, idx_train, idx_test), error=function(e){NA})
                                lasso = tryCatch(LASSO_train_test(x, y, idx_train, idx_test), error=function(e){NA})
                                ridge = tryCatch(RIDGE_train_test(x, y, idx_train, idx_test), error=function(e){NA})
                                elastic = tryCatch(ELASTIC_train_test(x, y, idx_train, idx_test), error=function(e){NA})
                                return(list(ols=ols, lasso=lasso, ridge=ridge, elastic=elastic))
                            }, mc.cores=parallel::detectCores())
    ols = list()
    lasso = list()
    ridge = list()
    elastic = list()
    for (rep in 1:r) {
        eval(parse(text=paste0("ols$rep_", rep, " = list()")))
        eval(parse(text=paste0("lasso$rep_", rep, " = list()")))
        eval(parse(text=paste0("ridge$rep_", rep, " = list()")))
        eval(parse(text=paste0("elastic$rep_", rep, " = list()")))
        for (fold in 1:k) {
            eval(parse(text=paste0("ols$rep_", rep, "$fold_", fold, " = perf$`", rep, "-", fold, "`$ols")))
            eval(parse(text=paste0("lasso$rep_", rep, "$fold_", fold, " = perf$`", rep, "-", fold, "`$lasso")))
            eval(parse(text=paste0("ridge$rep_", rep, "$fold_", fold, " = perf$`", rep, "-", fold, "`$ridge")))
            eval(parse(text=paste0("elastic$rep_", rep, "$fold_", fold, " = perf$`", rep, "-", fold, "`$elastic")))
        }
    }
    out = list(ols=ols,
                lasso=lasso,
                ridge=ridge,
                elastic=elastic)
    return(out)
}

GWAS_OLS_MULTIPLEREG = function(G, y, K, n_chr=7, fname_manhattan_plot_png="../res/Glyphosate_resistance_GWAS_OLS_multiplereg.png") {
    n = nrow(G)
    l = ncol(G)
    k = ncol(K)
    if (n != length(y)) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Genotype and phenotype data have different number of observations (populations).")
        return("ERROR!")
    }
        if (n != nrow(K)) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Genotype and kinship/covariate data have different number of observations (populations).")
        return("ERROR!")
    }
    allele_names = colnames(G)
    chr = as.factor(unlist(lapply(strsplit(allele_names, "-"), FUN=function(x){x[1]})))
    pos = as.numeric(unlist(lapply(strsplit(allele_names, "-"), FUN=function(x){x[2]})))
    allele = unlist(lapply(strsplit(allele_names, "-"), FUN=function(x){x[3]}))
    nloci = length(unique(paste0(chr, "-", pos)))
    freq = colMeans(G)
    ### Prepare the explanatory variables
    X = cbind(rep(1, n), K, G); p = ncol(X)
    b_hat = t(X) %*% solve(X %*% t(X)) %*% y
    allele_effects = b_hat[(p-l+1):p]
    LL_gauss = sum(-log10(dnorm(allele_effects, mean=0, sd=sd(allele_effects))))
    LL_laplace = sum(-log10(dLaplace(allele_effects, mean=0, b=sd(abs(allele_effects)))))
    if (LL_gauss < LL_laplace) {
        pval = 2 * pnorm(abs(allele_effects), mean=0, sd=sd(allele_effects), lower.tail=FALSE)
        statistic_dist = "Gaussian"
    } else {
        pval = 2*(1.0 - ExtDist::pLaplace(abs(allele_effects), mu=0, b=sd(abs(allele_effects))))
        statistic_dist = "Laplace"
    }
    lod = -log10(pval + 1e-12)
    # ### Define the output GWAS dataframe
    df_gwas = data.frame(X.chr, pos, alleles, freq, pheno="pheno_0", statistic, pvalue, lod)
    ### Preparing output
    df_other_effects = data.frame(id=c("Intercept", paste0("K", 1:k)), effects=b_hat[1:(k+1)])
    return(list(df_gwas=df_gwas,
                fname_manhattan_plot_png=fname_manhattan_plot_png,
                df_other_effects=df_other_effects))
}

GWAS_PLOT = function(df_gwas, n, fname_manhattan_plot_png, n_chr=7, statistic_dist="Gaussian", alpha=0.05, fname_genome_annotation="", flank_size_kb=0, plot_png=TRUE) {
    ### Define the output GWAS dataframe
    # df_gwas = data.frame(X.chr, pos, alleles, freq, pheno="pheno_0", statistic, pvalue, lod)
    ### Sort by loci
    df_gwas = df_gwas[order(df_gwas$alleles), ]
    df_gwas = df_gwas[order(df_gwas$pos), ]
    df_gwas = df_gwas[order(df_gwas$X.chr), ]
    ### Add consecutive locus position ids and colours corresponding to each chromosome
    df_gwas$consecutive_pos = NA
    df_gwas$colours_chr = NA
    vec_chr = sort(unique(df_gwas$X.chr))
    # vec_col = rep(rainbow(n_chr), times=ceiling(length(vec_chr)/n_chr))
    vec_col = rep(c("#1f78b4", "#b2df8a"), times=ceiling(length(vec_chr)/2))
    minimum = min(df_gwas$pos[df_gwas$X.chr == vec_chr[1]], na.rm=TRUE)
    df_gwas$consecutive_pos[df_gwas$X.chr == vec_chr[1]] = df_gwas$pos[df_gwas$X.chr == vec_chr[1]] - minimum
    df_gwas$colours_chr[df_gwas$X.chr == vec_chr[1]] = vec_col[1]
    maximum = max(df_gwas$consecutive_pos[df_gwas$X.chr == vec_chr[1]], na.rm=TRUE)
    for (i in 2:length(vec_chr)) {
        # i = 2
        minimum = min(df_gwas$pos[df_gwas$X.chr == vec_chr[i]], na.rm=TRUE)
        df_gwas$consecutive_pos[df_gwas$X.chr == vec_chr[i]] = (df_gwas$pos[df_gwas$X.chr == vec_chr[i]] - minimum) + maximum + 1
        df_gwas$colours_chr[df_gwas$X.chr == vec_chr[i]] = vec_col[i]
        maximum = max(df_gwas$consecutive_pos[df_gwas$X.chr == vec_chr[i]], na.rm=TRUE)
    }
    ### Output
    GWAS = list(df_gwas=df_gwas, fname_manhattan_plot_png=fname_manhattan_plot_png)
    ### Plot
    if (plot_png==TRUE) {
        png(fname_manhattan_plot_png, width=2600, height=1000)
        layout(matrix(c(1,1,1,2,1,1,1,3), nrow=2, byrow=TRUE))
        par(cex=2)
        # Manhattan plot
        plot(x=df_gwas$consecutive_pos, y=df_gwas$lod, type='n', xlab="Genome", ylab="-log10(p)", xaxt="n", las=1,
            main=gsub(".png", "", gsub("_", " ", basename(fname_manhattan_plot_png))))
        grid()
        chr_lab_pos = c()
        for (chr in vec_chr) {
            df_gwas_sub = df_gwas[df_gwas$X.chr == chr, ]
            chr_lab_pos = c(chr_lab_pos, median(df_gwas_sub$consecutive_pos))
            points(x=df_gwas_sub$consecutive_pos, y=df_gwas_sub$lod, pch=19, col=df_gwas_sub$colours_chr)
        }
        l = nrow(df_gwas)
        threshold = -log10(alpha/l)
        abline(h=threshold, col="red", lty=2)

        if (fname_genome_annotation != "") {
            GWAS$df_gwas_peak_genes = LOAD_GENOME_ANNOTATION_AND_MAP_WITH_GWAS(GWAS=GWAS,
                fname_genome_annotation=fname_genome_annotation,
                alpha=alpha,
                flank_size_kb=flank_size_kb)
        } else {
            GWAS$df_gwas_peak_genes = data.frame()
        }

        axis(side=1, at=chr_lab_pos, lab=vec_chr)
        # QQ-plot
        par(mar=c(5, 5, 1, 1))
        plot(x=sort(-log10(seq(0, 1, length=nrow(df_gwas)))), y=sort(df_gwas$lod), xlim=range(df_gwas$lod), ylim=range(df_gwas$lod),
            xlab="Expected -log10(p)", ylab="Observed -log10(p)")
        grid()
        lines(range(df_gwas$lod), range(df_gwas$lod), col="red", lty=2)
        legend("bottomright", statistic_dist, bty="n")
        # Allele effects distribution
        par(mar=c(5, 5, 1, 1))
        hist(df_gwas$statistic, main="", xlab="Estimated additive genetic effects")
        nloci = length(unique(paste0(df_gwas$X.chr, "-", df_gwas$pos)))
        legend("topright", legend=c(paste0("n=", n), paste0("p=", l), paste0("nloci=", nloci)), bty="n")
        grid()
        dev.off()
    }
    ### Output
    return(GWAS)
}

LAPLACE_DIST = function(x, mu=0, sigma=1, log.p=FALSE, fun=c("d", "p", "q", "r")[4]) {
    if (fun == "d") {
        ### Probability density function
        out = -log(2) - log(sigma) - (abs(x-mu)/sigma)
    } else if (fun == "p") {
        ### Cumulative distribution function
        idx_upper_tail = x >= mu
        out = -log(2) + (x-mu)/sigma
        out[idx_upper_tail] = -log(2) + (-x[idx_upper_tail]-mu)/sigma
    } else if (fun == "q") {
        ### Quantile function
        log.p = TRUE
        if(x<0.5) {
            out = mu + sigma*log(2*x)
        } else {
            out = mu - sigma*log(2*(1-x))
        } 
    } else {
        ### Random sampling from the Laplace distribution
        log.p = TRUE
        u = runif(x)
        idx_upper_tail = u >= 0.5
        out = mu + sigma*log(2*u)
        out[idx_upper_tail] = mu - sigma*log(2*(1-u[idx_upper_tail]))
    }
    if (log.p == FALSE){
        out = exp(out)
    }
    return(out)
}

ESTIMATE_GWALPHA_PVALUES = function(statistic) {
    # statistic = df_gwas$statistic
    ### Probably repeat r times with 10% of the statistic
    x = sample(statistic, size=100, replace=FALSE)
    vec_names = c("Gaussian", "Laplace", "Exponential_with_squared_effects", "Exponential_with_abs_effects")
    opt_gauss = optim(par=c(1, 1), method="L-BFGS-B", lower=c(-Inf, 1e-9), upper=c(+Inf, +Inf), fn=function(params){-sum(dnorm(x, mean=params[1], sd=params[2], log=TRUE))})
    opt_lapla = optim(par=c(1, 1), method="L-BFGS-B", lower=c(-Inf, 1e-12), upper=c(+Inf, +Inf), fn=function(params){-sum(LAPLACE_DIST(x=x, mu=params[1], sigma=params[2], log=TRUE, fun="d"))})
    opt_expon = optim(par=c(1), method="L-BFGS-B", lower=1e-20, upper=+Inf, fn=function(params){-sum(dexp(x^2, rate=params[1], log=TRUE))})
    opt_expab = optim(par=c(1), method="L-BFGS-B", lower=1e-20, upper=+Inf, fn=function(params){-sum(dexp(abs(x), rate=params[1], log=TRUE))})
    vec_negloglik = c(opt_gauss$value, opt_lapla$value, opt_expon$value, opt_expab$value)
    vec_negloglik[c(2)] = 10000
    vec_negloglik[c(3)] = 10000
    vec_negloglik[c(4)] = 10000
    idx = which(vec_negloglik == min(vec_negloglik))
    if (vec_names[idx] == "Gaussian") {
        distribution = "Gaussian"
        pval = 2*pnorm(abs(statistic), mean=opt_gauss$par[1], sd=opt_gauss$par[2], lower.tail=FALSE)
    } else if (vec_names[idx] == "Laplace") {
        distribution = "Laplace"
        pval = 2*LAPLACE_DIST(x=abs(statistic), mu=opt_lapla$par[1], sigma=opt_lapla$par[2], fun="p")
    } else if (vec_names[idx] == "Exponential_with_squared_effects") {
        distribution = "Exponential_with_squared_effects"
        pval = pexp(statistic^2, rate=opt_expon$par[1], lower.tail=FALSE)
    } else if (vec_names[idx] == "Exponential_with_abs_effects") {
        distribution = "Exponential_with_abs_effects"
        pval = pexp(abs(statistic), rate=opt_expon$par[1], lower.tail=FALSE)
    }
    return(list(pval=pval, distribution=distribution))
}
##################
### UNIT TESTS ###
##################
model_fitting = function() {
    seed = 42069
    set.seed(seed)
    n = 100
    p = 1000
    q = 10
    maf = 1e-4
    h2 = 0.75
    X_sim = matrix(runif(n*p, min=maf, max=1-maf), nrow=n)
    b = rep(0, p)
    idx_b = sort(sample(c(1:p), q))
    b[idx_b] = 1.0
    xb = X_sim %*% b
    v_xb = var(xb)
    v_e = (v_xb/h2) - v_xb
    e = rnorm(n, mean=0, sd=sqrt(v_e))
    y = xb + e
    y_sim = scale(y, center=T, scale=T)[,1]
    idx_train = 1:90
    idx_test = 91:100

    test_that(
        "PERF", {
            print("PERF:")
            set.seed(seed)
            expect_equal(PERF(y_test=y_sim, y_hat=y_sim*+2)$cor, +1, tolerance=1e-7)
            expect_equal(PERF(y_test=y_sim, y_hat=y_sim*-1)$cor, -1, tolerance=1e-7)
            expect_equal(PERF(y_test=y_sim, y_hat=xb)$mae, +5, tolerance=0.5)
            expect_equal(PERF(y_test=y_sim, y_hat=xb)$mbe, -5, tolerance=0.5)
            expect_equal(PERF(y_test=y_sim, y_hat=xb)$rmse, +5, tolerance=0.5)
        }
    )

    test_that(
        "OLS_train_test", {
            print("OLS_train_test:")
            set.seed(seed)
            cor1 = OLS_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = OLS_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 1.0, tolerance=1e-7)
            expect_equal(cor2, 0.191562494269293, tolerance=1e-7)
        }
    )

    test_that(
        "LASSO_train_test", {
            print("LASSO_train_test:")
            set.seed(seed)
            cor1 = LASSO_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = LASSO_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.97135207990405, tolerance=1e-7)
            expect_equal(cor2, 0.62379540969574, tolerance=1e-7)
        }
    )

    test_that(
        "RIDGE_train_test", {
            print("RIDGE_train_test:")
            set.seed(seed)
            cor1 = RIDGE_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = RIDGE_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.996749570837874, tolerance=1e-7)
            expect_equal(cor2, 0.139118883736874, tolerance=1e-7)
        }
    )

    test_that(
        "ELASTIC_train_test", {
            print("ELASTIC_train_test:")
            set.seed(seed)
            cor1 = ELASTIC_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = ELASTIC_train_test(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.9998965, tolerance=1e-7)
            expect_equal(cor2, 0.6237954, tolerance=1e-7)
        }
    )

    test_that(
        "KFOLD_CV", {
            print("KFOLD_CV:")
            kfold = KFOLD_CV(x=X_sim, y=y_sim, r=2, k=10)
            ols_corr = mean(unlist(lapply(kfold$ols, FUN=function(x){lapply(x, FUN=function(y){y$cor})})))
            elastic_corr = mean(unlist(lapply(kfold$elastic, FUN=function(x){lapply(x, FUN=function(y){y$cor})})))
            ols_rmse = mean(unlist(lapply(kfold$ols, FUN=function(x){lapply(x, FUN=function(y){y$rmse})})))
            elastic_rmse = mean(unlist(lapply(kfold$elastic, FUN=function(x){lapply(x, FUN=function(y){y$rmse})})))
            expect_equal(ols_corr,     0.1555020, tolerance=1e-7)
            expect_equal(elastic_corr, 0.5628153, tolerance=1e-7)
            expect_equal(ols_rmse,     0.9575458, tolerance=1e-7)
            expect_equal(elastic_rmse, 0.8135351, tolerance=1e-7)
        }
    )
}

# options(digits.secs=7)
# start_time = Sys.time()
# model_fitting()
# end_time = Sys.time()
# print(end_time - start_time)