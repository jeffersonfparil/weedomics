#################################################################################################################
#################################################################################################################
#################################################################################################################

### Trying to make a better catwalk, i.e. modelling platform

library(testthat)
library(glmnet)
library(raster)
library(parallel)

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

fn_ols = function(x, y, idx_train, idx_test) {
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

fn_ridge = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    mod_ridge = cv.glmnet(x=x_train, y=y_train, alpha=0.0)
    y_hat = predict(mod_ridge, newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

fn_lasso = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    mod_lasso = cv.glmnet(x=x_train, y=y_train, alpha=1.0)
    y_hat = predict(mod_lasso, newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

fn_elastic = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]
    y_test = y[idx_test]
    mod_elastic = cv.glmnet(x=x_train, y=y_train)
    y_hat = predict(mod_elastic, newx=x_test, s="lambda.min")[,1]
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    return(out)
}

fn_bayesa = function(x, y, idx_train, idx_test) {
    prefix = paste0("BayesA_", gsub(" ", "", date()), round(runif(1)*1e9), "-")
    y_new = y
    y_new[idx_test] = NA
    fit = BGLR::BGLR(y=y_new, ETA=list(list(X=x, model='BayesA')), nIter=6000,burnIn=1000, saveAt=prefix, verbose=FALSE)
    out = PERF(y[idx_test], fit$yHat[idx_test])
    out$idx_test = idx_test
    return(out)
}

fn_bayesb = function(x, y, idx_train, idx_test) {
    prefix = paste0("BayesB_", gsub(" ", "", date()), round(runif(1)*1e9), "-")
    y_new = y
    y_new[idx_test] = NA
    fit = BGLR::BGLR(y=y_new, ETA=list(list(X=x, model='BayesB')), nIter=6000,burnIn=1000, saveAt=prefix, verbose=FALSE)
    out = PERF(y[idx_test], fit$yHat[idx_test])
    out$idx_test = idx_test
    return(out)
}

fn_bayesc = function(x, y, idx_train, idx_test) {
    prefix = paste0("BayesC_", gsub(" ", "", date()), round(runif(1)*1e9), "-")
    y_new = y
    y_new[idx_test] = NA
    fit = BGLR::BGLR(y=y_new, ETA=list(list(X=x, model='BayesC')), nIter=6000,burnIn=1000, saveAt=prefix, verbose=FALSE)
    out = PERF(y[idx_test], fit$yHat[idx_test])
    out$idx_test = idx_test
    return(out)
}

fn_logistic_ridge = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]  ### y is a binary vector
    y_test = y[idx_test]    ### y is a binary vector
    mod_elastic = cv.glmnet(x=x_train, y=y_train, family="binomial", alpha=0.0)
    r_hat = predict(mod_elastic, type="response", newx=x_test, s="lambda.min")[,1]
    y_hat = ifelse(r_hat >= 0.5, 1, 0)
    l_hat = predict(mod_elastic, type="link", newx=x_test, s="lambda.min")[,1]
    model = glm(y_hat ~ l_hat, family=binomial(link="logit"))
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    out$y_hat = y_hat ### predicted binary response
    out$l_hat = l_hat ### linear predictors
    return(out)
}

fn_logistic_lasso = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]  ### y is a binary vector
    y_test = y[idx_test]    ### y is a binary vector
    mod_elastic = cv.glmnet(x=x_train, y=y_train, family="binomial", alpha=1.0)
    r_hat = predict(mod_elastic, type="response", newx=x_test, s="lambda.min")[,1]
    y_hat = ifelse(r_hat >= 0.5, 1, 0)
    l_hat = predict(mod_elastic, type="link", newx=x_test, s="lambda.min")[,1]
    model = glm(y_hat ~ l_hat, family=binomial(link="logit"))
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    out$y_hat = y_hat ### predicted binary response
    out$l_hat = l_hat ### linear predictors
    return(out)
}

fn_logistic_elastic = function(x, y, idx_train, idx_test) {
    x_train = x[idx_train, ]
    x_test = x[idx_test, ]
    y_train = y[idx_train]  ### y is a binary vector
    y_test = y[idx_test]    ### y is a binary vector
    mod_elastic = cv.glmnet(x=x_train, y=y_train, family="binomial")
    r_hat = predict(mod_elastic, type="response", newx=x_test, s="lambda.min")[,1]
    y_hat = ifelse(r_hat >= 0.5, 1, 0)
    l_hat = predict(mod_elastic, type="link", newx=x_test, s="lambda.min")[,1]
    model = glm(y_hat ~ l_hat, family=binomial(link="logit"))
    out = PERF(y_test, y_hat)
    out$idx_test = idx_test
    out$y_hat = y_hat ### predicted binary response
    out$l_hat = l_hat ### linear predictors
    return(out)
}

KFOLD_CV = function(x, y, r=5, k=10, vec_models=c("ols", "ridge", "lasso", "elastic", "bayesa", "bayesb", "bayesc", "logistic_ridge", "logistic_lasso", "logistic_elastic"), logistic_y_bin_threshold=0.5) {
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
                                for (model in vec_models) {
                                    if (grepl("logistic", model)) {
                                        y_bin = ifelse(y >= logistic_y_bin_threshold, 1.0, 0.0)
                                        eval(parse(text=paste0(model, " = tryCatch(fn_", model, "(x, y_bin, idx_train, idx_test), error=function(e){NA})")))
                                    } else {
                                        eval(parse(text=paste0(model, " = tryCatch(fn_", model, "(x, y, idx_train, idx_test), error=function(e){NA})")))
                                    }
                                }
                                out = eval(parse(text=paste0("list(", paste(paste0(vec_models, "=", vec_models), collapse=","), ")")))
                                return(out)
                            }, mc.cores=parallel::detectCores())
    for (model in vec_models) {
        eval(parse(text=paste0(model, " = list()")))
        for (rep in 1:r) {
            eval(parse(text=paste0(model, "$rep_", rep, " = list()")))
            for (fold in 1:k) {
                eval(parse(text=paste0(model, "$rep_", rep, "$fold_", fold, " = perf$`", rep, "-", fold, "`$", model)))
            }
        }
    }
    out = eval(parse(text=paste0("list(", paste(paste0(vec_models, "=", vec_models), collapse=","), ")")))
    return(out)
}

LOAD_PHENOTYPES = function(fname_phenotype="res/phenotype_data.csv", batch="all", phenotype_names="all") {
    # fname_phenotype = "res/phenotype_data.csv"
    # batch = c("all", 2019)[1]
    # batch = c("all", 2019)[2]
    # phenotype_names = list(c("all"), c("Glyphosate"), c("Glyphosate", "Sulfometuron"))[[1]]
    # phenotype_names = list(c("all"), c("Glyphosate"), c("Glyphosate", "Sulfometuron"))[[2]]
    # phenotype_names = list(c("all"), c("Glyphosate"), c("Glyphosate", "Sulfometuron"))[[3]]
    df_phenotype = read.csv(fname_phenotype)
    if (batch!="all") {
        idx_rows = df_phenotype$Batch == batch
        if (sum(idx_rows) == 0) {
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print("ERROR!")
            print("Batch ID to include is incorrect. Please select from 2018 and 2019.")
            return("ERROR!")
        }
        df_phenotype = df_phenotype[idx_rows, ]
    }
    variable_names = c("X.Population", "Pool_size", "Latitude", "Longitude", phenotype_names)
    if (phenotype_names[1]=="all") {
        phenotype_names = eval(parse(text=paste0("colnames(df_phenotype)[!(", paste(paste0("(colnames(df_phenotype)=='",  c("Batch", variable_names[1:4]), "')"), collapse=" | "), ")]")))
        variable_names = c("X.Population", "Pool_size", "Latitude", "Longitude", phenotype_names)
    }
    idx_cols = eval(parse(text=paste0("which(", paste(paste0("(colnames(df_phenotype) == '", variable_names, "')"), collapse=" | "), ")")))
    if (length(idx_cols) < 5) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Please select a valid phenotype name. Select from:")
        print(eval(parse(text=paste0("colnames(df_phenotype)[!(", paste(paste0("(colnames(df_phenotype)=='",  c("Batch", variable_names[1:4]), "')"), collapse=" | "), ")]"))))
        return("ERROR!")
    }
    idx_rows = eval(parse(text=paste(paste0("!is.na(df_phenotype[ ,which(colnames(df_phenotype) == '", phenotype_names, "')])"), collapse=" & ")))
    df_phenotype = df_phenotype[idx_rows, idx_cols]
    return(df_phenotype)
}

LOAD_GENOTYPES = function(fname_genotype="res/genotype_data-allele_frequencies-p_minus_one.csv", maf=0.01) {
    # fname_genotype = "res/genotype_data-allele_frequencies-p_minus_one.csv"
    # maf = 0.01
    df_must_transpose = read.csv(fname_genotype)
    vec_locus_id = c("X.chr", "pos", "allele")
    idx_cols = eval(parse(text=paste0("which(", paste(paste0("(colnames(df_must_transpose)!='", vec_locus_id, "')"), collapse=" & "), ")")))
    df_genotype = data.frame(X.Population=colnames(df_must_transpose)[idx_cols], 
                             t(df_must_transpose[, idx_cols, ]))
    colnames(df_genotype) = c("X.Population", 
                    eval(parse(text=paste0("paste0(", 
                        paste(paste0("df_must_transpose$", vec_locus_id), collapse=", '-', "), 
                    ")")))
                    )
    rownames(df_genotype) = df_genotype$X.Population
    idx_cols = c(2:ncol(df_genotype))
    idx_cols = c(1, which(!((colMeans(df_genotype[ ,idx_cols]) < maf) |
                           (colMeans(df_genotype[ ,idx_cols]) > (1-maf))))
                )
    if (length(idx_cols)==1) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print(paste0("No locus retained after filtering for a minimum allele frequency (MAF) of ", maf, "."))
        print("Please reduce MAF.")
        return("ERROR!")
    }
    df_genotype = df_genotype[, idx_cols]
    return(df_genotype)
}

LOAD_ENVIRONMENTAL_DATA = function(dname_climate_grid_data="res/bom_grid_data_1991_to_2020", dname_soil_grid_data="res/fao_soil_data_version1.2", extension_name="asc") {
    # dname_climate_grid_data="res/bom_grid_data_1991_to_2020"
    # dname_soil_grid_data="res/fao_soil_data_version1.2"
    # extension_name="asc"
    fnames_envi = list.files(dname_climate_grid_data)[grep(paste0(extension_name, "$"), list.files(dname_climate_grid_data))]
    fnames_soil = list.files(dname_soil_grid_data)[grep(paste0(extension_name, "$"), list.files(dname_soil_grid_data))]
    if (length(fnames_envi)==0) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print(paste0("Climatic grid data with *.", extension_name, " extension name not found in directory: ", dname_climate_grid_data, "."))
        return("ERROR!")
    }
    if (length(fnames_soil)==0) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print(paste0("Soil grid data with *.", extension_name, " extension name not found in directory: ", dname_soil_grid_data, "."))
        return("ERROR!")
    }
    fnames_envi = paste0(dname_climate_grid_data, "/", fnames_envi)
    fnames_soil = paste0(dname_soil_grid_data, "/", fnames_soil)
    raster_layers = list()
    for (f in c(fnames_envi, fnames_soil)) {
        # f = c(fnames_envi, fnames_soil)[1]
        name = gsub(".asc", "", basename(f))
        eval(parse(text=paste0("raster_layers$`", name, "` = raster::raster(f)")))
    }
    return(raster_layers)
}

MERGE_PHENOTYPE_WITH_GENOTYPE_DATA = function(df_phenotype, df_genotype) {
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
    # df_genotype = LOAD_GENOTYPES(maf=0.2)
    if (sum(colnames(df_phenotype)=="X.Population") < 1){
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Population name column, i.e. `X.Population` not found in the phenotype dataframe.")
        return("ERROR!")
    }
    if (sum(colnames(df_genotype)=="X.Population") < 1){
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Population name column, i.e. `X.Population` not found in the genotype dataframe.")
        return("ERROR!")
    }
    df_phe_gen = merge(df_phenotype, df_genotype, by="X.Population")
    if (nrow(df_phe_gen) < 1) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("No data retained after merging. Please check the phenotype and genotype data.")
        return("ERROR!")
    }
    return(df_phe_gen)
}

MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA = function(df_phenotype, raster_layers) {
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
    # raster_layers = LOAD_ENVIRONMENTAL_DATA()
    df_phe_env = df_phenotype
    if (sum(colnames(df_phe_env)=="Latitude") < 1){
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Latitude column, i.e. `Latitude` not found in the phenotype dataframe.")
        return("ERROR!")
    }
    if (sum(colnames(df_phe_env)=="Longitude") < 1){
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print("Longitude column, i.e. `Longitude` not found in the phenotype dataframe.")
        return("ERROR!")
    }
    idx_rows = !is.na(df_phe_env$Latitude) & !is.na(df_phe_env$Longitude)
    df_phe_env = df_phe_env[idx_rows, ]
    list_extent = list()
    epsilon = 1e-8
    for (i in 1:nrow(df_phe_env)) {
        # i = 1
        x = df_phe_env$Longitude[i]
        y = df_phe_env$Latitude[i]
        ext = as(extent(x-epsilon, x+epsilon, y-epsilon, y+epsilon), 'SpatialPolygons')
        crs(ext) = "+proj=longlat +datum=WGS84 +no_defs"
        eval(parse(text=paste0("list_extent$`", i, "` = ext")))
    }
    list_envi = parallel::mclapply(1:length(raster_layers), function(k) {
                    z = c()
                    for (i in 1:nrow(df_phe_env)) {
                        # i = 1
                        cropped = raster::crop(raster_layers[[k]], list_extent[[i]])
                        z = c(z, cropped@data@values)
                    }
                    return(z)
                }, mc.cores=parallel::detectCores())
    names_envi = names(raster_layers)
    for (k in 1:length(names_envi)) {
        # k = 1
        eval(parse(text=paste0("df_phe_env$`", names_envi[k], "` = list_envi[[k]]")))
    }
    return(df_phe_env)
}

MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA = function(df_phenotype, df_genotype, raster_layers) {
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
    # df_genotype = LOAD_GENOTYPES(maf=0.2)
    # raster_layers = LOAD_ENVIRONMENTAL_DATA()
    df_phe_gen = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype, df_genotype)
    df_phe_env = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype, raster_layers)
    names_envi = names(raster_layers)
    idx_cols = which(colnames(df_phe_env) %in% c("X.Population", names_envi))
    df_phe_env = df_phe_env[, idx_cols]
    df_phe_gen_env = merge(df_phe_gen, df_phe_env, by="X.Population")
    return(df_phe_gen_env)
}

##################
### UNIT TESTS ###
##################
tests = function () {
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
        "fn_ols", {
            print("fn_ols:")
            set.seed(seed)
            cor1 = fn_ols(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = fn_ols(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 1.0, tolerance=1e-7)
            expect_equal(cor2, 0.191562494269293, tolerance=1e-7)
        }
    )


    test_that(
        "fn_lasso", {
            print("fn_lasso:")
            set.seed(seed)
            cor1 = fn_lasso(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = fn_lasso(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.97135207990405, tolerance=1e-7)
            expect_equal(cor2, 0.62379540969574, tolerance=1e-7)
        }
    )

    test_that(
        "fn_ridge", {
            print("fn_ridge:")
            set.seed(seed)
            cor1 = fn_ridge(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = fn_ridge(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.996749570837874, tolerance=1e-7)
            expect_equal(cor2, 0.139118883736874, tolerance=1e-7)
        }
    )

    test_that(
        "fn_elastic", {
            print("fn_elastic:")
            set.seed(seed)
            cor1 = fn_elastic(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = fn_elastic(x=X_sim, y=y_sim, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 0.97135207990405, tolerance=1e-7)
            expect_equal(cor2, 0.62379540969574, tolerance=1e-7)
        }
    )

    ### Add unit tests for Bayesian models and logistic models

    test_that(
        "fn_logistic_elastic", {
            print("fn_logistic_elastic:")
            set.seed(seed)
            y_bin = ifelse(y_sim >= 0.5, 1.0, 0.0)
            cor1 = fn_logistic_elastic(x=X_sim, y=y_bin, idx_train=idx_train, idx_test=idx_train)$cor
            cor2 = fn_logistic_elastic(x=X_sim, y=y_bin, idx_train=idx_train, idx_test=idx_test)$cor
            expect_equal(cor1, 1.0000000, tolerance=1e-7)
            expect_equal(cor2, 0.2721655, tolerance=1e-7)
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

    test_that(
        "LOAD_PHENOTYPES", {
            print("LOAD_PHENOTYPES:")
            df1 = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Sulfometuron"))
            df2 = LOAD_PHENOTYPES(batch="SHOULD_ERROR")
            df3 = LOAD_PHENOTYPES(phenotype_names="SHOULD_ERROR_TOO")
            expect_equal(dim(df1), c(47, 6), tolerance=1e-7)
            expect_equal(df2, "ERROR!")
            expect_equal(df3, "ERROR!")
        }
    )

    test_that(
        "LOAD_GENOTYPES", {
            print("LOAD_GENOTYPES:")
            df1 = LOAD_GENOTYPES(maf=0.1)
            df2 = LOAD_GENOTYPES(maf=0.5)
            expect_equal(dim(df1), c(119, 855), tolerance=1e-7)
            expect_equal(df2, "ERROR!")
        }
    )

    test_that(
        "LOAD_ENVIRONMENTAL_DATA", {
            print("LOAD_ENVIRONMENTAL_DATA:")
            df1 = LOAD_ENVIRONMENTAL_DATA()
            df2 = LOAD_ENVIRONMENTAL_DATA(dname_climate_grid_data="SHOULD_ERROR")
            df3 = LOAD_ENVIRONMENTAL_DATA(dname_soil_grid_data="SHOULD_ERROR_TOO")
            expect_equal(names(df1)[c(1, 4, 12, 13, 18)], c("evapora", "rh1500h", "nutavail", "nutreten", "workable"), tolerance=1e-7)
            expect_equal(df3, "ERROR!")
            expect_equal(df2, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_GENOTYPE_DATA", {
            print("MERGE_PHENOTYPE_WITH_GENOTYPE_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
            df_genotype = LOAD_GENOTYPES(maf=0.2)
            df1 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=df_phenotype, df_genotype=df_genotype)
            df2 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=data.frame(), df_genotype=df_genotype)
            df3 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=data.frame(X.Population=NA), df_genotype=df_genotype)
            expect_equal(dim(df1), c(113, 322), tolerance=1e-7)
            expect_equal(df3, "ERROR!")
            expect_equal(df2, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA", {
            print("MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
            raster_layers = LOAD_ENVIRONMENTAL_DATA()
            df1 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=df_phenotype, raster_layers=raster_layers)
            df2 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=data.frame(Longitude=NA), raster_layers=raster_layers)
            df3 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=data.frame(Latitude=NA), raster_layers=raster_layers)
            expect_equal(dim(df1), c(113, 24), tolerance=1e-7)
            expect_equal(df2, "ERROR!")
            expect_equal(df3, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA", {
            print("MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate", "Clethodim"))
            df_genotype = LOAD_GENOTYPES(maf=0.2)
            raster_layers = LOAD_ENVIRONMENTAL_DATA()
            df1 = MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA(df_phenotype=df_phenotype, df_genotype=df_genotype, raster_layers=raster_layers)
            expect_equal(dim(df1), c(113, 340), tolerance=1e-7)
        }
    )
}
# tests()
