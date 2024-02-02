setwd("/data-weedomics-1/weedomics/exploratory")
source("src/modelling_and_crossfold_validation.r")

fname_phenotype = "res/phenotype_data.csv"
fname_genotype = "res/genotype_data-allele_frequencies-p_minus_one.csv"
dname_climate_grid_data = "res/bom_grid_data_1991_to_2020"
dname_soil_grid_data = "res/fao_soil_data_version1.2"
fname_csiro_soil_data = "res/CSIRO_AusSoil_complete.csv"

df_genotype = LOAD_GENOTYPES(fname_genotype=fname_genotype,
                                maf=0.01)
raster_layers = LOAD_ENVIRONMENTAL_DATA(dname_climate_grid_data=dname_climate_grid_data,
                                dname_soil_grid_data=dname_soil_grid_data,
                                extension_name="asc")
df_covariates_csiro_soil_data = read.csv(file=fname_csiro_soil_data)
colnames(df_covariates_csiro_soil_data)[grepl("Accession", colnames(df_covariates_csiro_soil_data))] = "X.Population"

vec_names_ids = c("X.Population", "Batch", "Pool_size", "Latitude", "Longitude", "long", "lat")
vec_names_phenotypes = colnames(read.csv(fname_phenotype))
vec_names_phenotypes = vec_names_phenotypes[!(vec_names_phenotypes %in% vec_names_ids)]

for (herbi in vec_names_phenotypes) {
    # herbi = vec_names_phenotypes[2]
    df_phenotype = LOAD_PHENOTYPES(fname_phenotype=fname_phenotype, batch="all", phenotype_names=herbi)
    df_phenotype = merge(df_phenotype, df_covariates_csiro_soil_data, by="X.Population")
    df = MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA(df_phenotype=df_phenotype, df_genotype=df_genotype, raster_layers=raster_layers)
    vec_names_csiro_soil_data = colnames(df_covariates_csiro_soil_data)
    vec_names_csiro_soil_data = vec_names_csiro_soil_data[!(vec_names_csiro_soil_data %in% vec_names_ids)]
    vec_names_phenotypes = colnames(df_phenotype)
    vec_names_phenotypes = vec_names_phenotypes[!(vec_names_phenotypes %in% vec_names_ids)]
    vec_names_phenotypes = vec_names_phenotypes[!(vec_names_phenotypes %in% vec_names_csiro_soil_data)]
    vec_names_loci = colnames(df_genotype)
    vec_names_loci = vec_names_loci[!(vec_names_loci %in% vec_names_ids)]
    # vec_names_fao_bom_data = colnames(df)
    # vec_names_fao_bom_data = vec_names_fao_bom_data[!(vec_names_fao_bom_data %in% c(vec_names_ids, vec_names_csiro_soil_data, vec_names_phenotypes, vec_names_loci))]
    y = df[, colnames(df) %in% herbi]
    G = df[, colnames(df) %in% vec_names_loci]
    S = df[, colnames(df) %in% vec_names_csiro_soil_data]

    if (var(y, na.rm=TRUE) < 0) {
        print("No variance in the phenotype data.")
        next
    }
    vec_var_loci = apply(G, MARGIN=2, FUN=var, na.rm=TRUE)
    G = G[, vec_var_loci != 0]
    vec_var_csiro_soil_data = apply(S, MARGIN=2, FUN=var)
    S = S[, !((vec_var_csiro_soil_data == 0) | (is.na(vec_var_csiro_soil_data)))]

    y = scale(y, scale=TRUE, center=TRUE)
    G = scale(G, scale=TRUE, center=TRUE)
    S = scale(S, scale=TRUE, center=TRUE)

    for (explanatory_variables in c("G", "E_and_G")) {
        # explanatory_variables = "G"
        options(digits.secs=7)
        start_time = Sys.time()
        if (explanatory_variables == "G") {
            eval(parse(text=paste0("kfold_", explanatory_variables, " = KFOLD_CV(x=G, y=y, r=10, k=10)")))
        } else if (explanatory_variables == "E_and_G") {
            eval(parse(text=paste0("kfold_", explanatory_variables, " = KFOLD_CV(x=cbind(S, G), y=y, r=10, k=10)")))
        }
        end_time = Sys.time()
        print(end_time - start_time)

        vec_models = c("ols", "lasso", "ridge", "elastic", "bayesa", "bayesb", "bayesc")
        col_bars = rgb(0.1, 0.1, 0.1, alpha=0.5)
        col_points = rgb(0.9, 0.2, 0.1, alpha=0.5)
        col_line = rgb(0.9, 0.2, 0.1, alpha=0.5)
        for (model in vec_models) {
            # model = "ols"
            y_idx = suppressWarnings(eval(parse(text=paste0("unlist(lapply(kfold_", explanatory_variables, "$", model, ", FUN=function(x){lapply(x, FUN=function(y){if(!is.na(y)){y$idx_test}})}))"))))
            y_true = suppressWarnings(eval(parse(text=paste0("unlist(lapply(kfold_", explanatory_variables, "$", model, ", FUN=function(x){lapply(x, FUN=function(y){if(!is.na(y)){y$y_test}})}))"))))
            y_pred = suppressWarnings(eval(parse(text=paste0("unlist(lapply(kfold_", explanatory_variables, "$", model, ", FUN=function(x){lapply(x, FUN=function(y){if(!is.na(y)){y$y_hat}})}))"))))

            y_true = (y_true * attr(y, "scaled:scale")) + attr(y, "scaled:center")
            y_pred = (y_pred * attr(y, "scaled:scale")) + attr(y, "scaled:center")

            df_y = data.frame(y_idx, y_true, y_pred)

            y_true = aggregate(y_true ~ y_idx, data=df_y, FUN=mean)
            y_pred = aggregate(y_pred ~ y_idx, data=df_y, FUN=mean)
            y_pred_sd = aggregate(y_pred ~ y_idx, data=df_y, FUN=sd)
            df_y = merge(y_true, y_pred, by="y_idx"); colnames(df_y) = c("y_idx", "y_true", "y_pred_mean")
            df_y = merge(df_y, y_pred_sd, by="y_idx"); colnames(df_y) = c("y_idx", "y_true", "y_pred_mean", "y_pred_sd")

            corr = cor(df_y$y_true, df_y$y_pred_mean)
            mae = mean(abs(df_y$y_true - df_y$y_pred_mean))
            rmse = sqrt(mean((df_y$y_true - df_y$y_pred_mean)^2))
            r2 = 1.00 - (sum((df_y$y_true - df_y$y_pred_mean)^2) / sum((df_y$y_true - mean(df_y$y_true))^2))

            svg(paste0("res/genomic_prediction-", model, "-", herbi, ".svg"))
            plot(x=0, y=0, xlim=range(df_y$y_true), ylim=range(c(df_y$y_pred_mean-df_y$y_pred_sd, df_y$y_pred_mean+df_y$y_pred_sd)), type="n", 
                xlab="Observed (survival rate, %)", 
                ylab="Predicted (survival rate, %)", 
                main=paste0(herbi, " Resistance\n(", model, "; y~", explanatory_variables, ")"))
            grid()
            arrows(x0=df_y$y_true, x1=df_y$y_true, y0=df_y$y_pred_mean-df_y$y_pred_sd, y1=df_y$y_pred_mean+df_y$y_pred_sd, code=3, col=col_bars, angle=90, length=0.0, lty=1)
            points(x=df_y$y_true, y=df_y$y_pred_mean, pch=19, col=col_points)
            abline(lm(y_pred_mean ~ y_true, data=df_y), lty=2, col=col_line)
            legend("bottomright", legend=c(
                paste0("Pearson's correlation = ", round(100*corr, 2), "%"),
                paste0("Mean absolute error = ", round(mae, 2)),
                paste0("Root mean squared error = ", round(rmse, 2)),
                paste0("Coefficient of determination = ", round(100*r2, 2), "%")
            ), bty="n")
            dev.off()
        }
    }
}
