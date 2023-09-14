library(testthat)
library(raster) # sudo apt install libgdal-dev
library(parallel)

dirname_src = utils::getSrcDirectory(function(){})[1]
if (is.na(dirname_src)) {
    dirname_src = "."
}
setwd(dirname_src)


LOAD_PHENOTYPES = function(fname_phenotype="../res/phenotype_data.csv", batch="all", phenotype_names="all") {
    # fname_phenotype = "../res/phenotype_data.csv"
    # batch = c("all", 2019)[1]
    # batch = c("all", 2019)[2]
    # phenotype_names = list(c("all"), c("Glyphosate"), c("Glyphosate", "Sulfometuron"))[[1]]
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
        print("Please validate the expected column names. Expected:")
        print(variable_names)
        print("From the input phenotype data:")
        print(eval(parse(text=paste0("colnames(df_phenotype)[!(", paste(paste0("(colnames(df_phenotype)=='",  c("Batch", variable_names[1:4]), "')"), collapse=" | "), ")]"))))
        return("ERROR!")
    }
    idx_rows = eval(parse(text=paste(paste0("!is.na(df_phenotype[ ,which(colnames(df_phenotype) == '", phenotype_names, "')])"), collapse=" & ")))
    df_phenotype = df_phenotype[idx_rows, idx_cols]
    return(df_phenotype)
}

LOAD_GENOTYPES = function(fname_genotype="../res/genotype_data-allele_frequencies-p_minus_one.csv", maf=0.01) {
    # fname_genotype = "../res/genotype_data-allele_frequencies-p_minus_one.csv"
    # maf = 0.001
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
    freqs = colMeans(df_genotype[ ,idx_cols])
    idx_cols = c(1, 1 + which(!(freqs < maf) | (freqs > (1-maf))))
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

LOAD_ENVIRONMENTAL_DATA = function(dname_climate_grid_data="../res/bom_grid_data_1991_to_2020", dname_soil_grid_data="../res/fao_soil_data_version1.2", extension_name="asc") {
    # dname_climate_grid_data="../res/bom_grid_data_1991_to_2020"
    # dname_soil_grid_data="../res/fao_soil_data_version1.2"
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

LOAD_KINSHIP_AND_SORT = function(population_names_sorted, fname_fst="../res/fst.csv", col_id_population_names=1) {
    FST = read.csv(fname_fst)
    idx = c(1:ncol(FST))[c(1:ncol(FST)) != col_id_population_names]
    if (length(idx) == ncol(FST)) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print(paste0("Column containing the population names at '", col_id_population_names, "' not found."))
        return("ERROR!")
    }
    K = as.matrix(FST[, idx])
    if (nrow(K) != ncol(K)) {
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("ERROR!")
        print(paste0("Fst matrix not square. Please make sure the file '", fname_fst, "' has n rows and n+1 columns, where the extra column correspond to the population names."))
        return("ERROR!")
    }
    idx = FST[, col_id_population_names] %in% population_names_sorted
    K = K[idx, idx]
    idx = c(1:nrow(K))[order(FST[,1][idx])]
    K = K[idx, idx]
}

LOAD_GENOME_ANNOTATION_AND_MAP_WITH_GWAS = function(GWAS, fname_genome_annotation="../res/genome_annotation.gff", alpha=0.05, flank_size_kb=0) {
    gff = read.delim(fname_genome_annotation, sep="\t", header=FALSE, comment.char="#")
    p = nrow(GWAS$df_gwas)
    bonferroni_threshold = -log10(alpha/p)
    df = GWAS$df_gwas[GWAS$df_gwas$lod >= bonferroni_threshold, ]
    if (nrow(df)==0) {
        df_gwas_peak_genes = data.frame(chr=c(), pos=c(), freq=c(), eff=c(), lod=c(), gene=c(), gene_start_pos=c(), gene_end_pos=c())
        return(df_gwas_peak_genes)
    }
    vec_chr = c()
    vec_pos = c()
    vec_freq = c()
    vec_eff = c()
    vec_lod = c()
    vec_gene = c()
    vec_gene_start_pos = c()
    vec_gene_end_pos = c()
    if (length(df$chr)==0) {
        df$chr = df$X.chr
        GWAS$df_gwas$chr = GWAS$df_gwas$X.chr
    }
    if (length(df$allele_effects)==0) {
        df$allele_effects = df$statistic
        GWAS$df_gwas$allele_effects = GWAS$df_gwas$statistic
    }
    for (i in 1:nrow(df)) {
        # i = 1
        chr = as.character(df$chr[i])
        pos = df$pos[i]
        pos = df$pos[i]
        freq = df$freq[i]
        eff = df$allele_effects[i]
        lod = df$lod[i]
        idx = (gff$V1==chr) & (((gff$V4-(flank_size_kb*1e3)) <= pos) & ((gff$V5+(flank_size_kb*1e3)) >= pos) | 
                               ((gff$V4-(flank_size_kb*1e3)) >= pos) & ((gff$V5+(flank_size_kb*1e3)) <= pos)) & (gff$V3 == "gene")
        if (sum(idx) > 0) {
            labels = unlist(strsplit(gff$V9[idx][1], ";"))
            gene = gsub("Name=", "", labels[grepl("Name=", labels)])
            vec_chr = c(vec_chr, chr)
            vec_pos = c(vec_pos, pos)
            vec_freq = c(vec_freq, freq)
            vec_eff = c(vec_eff, eff)
            vec_lod = c(vec_lod, lod)
            vec_gene = c(vec_gene, gene)
            vec_gene_start_pos = c(vec_gene_start_pos, gff$V4[idx][1])
            vec_gene_end_pos = c(vec_gene_end_pos, gff$V5[idx][1])
        }
    }
    df_gwas_peak_genes = data.frame(chr=vec_chr, pos=vec_pos, freq=vec_freq, eff=vec_eff, lod=vec_lod, gene=vec_gene, gene_start_pos=vec_gene_start_pos, gene_end_pos=vec_gene_end_pos)
    # ### Manhattan plot with peak labels
    # png(fname_manhattan_plot_with_peak_labels_png, width=2000, height=1000)
    # par(cex=2)
    # plot(x=GWAS$df_gwas$consecutive_pos, y=GWAS$df_gwas$lod, type='n', xlab="Genome", ylab="-log10(p)", xaxt="n", las=1,
    #     main=gsub(".png", "", gsub("_", " ", basename(fname_manhattan_plot_with_peak_labels_png))))
    # grid()
    # vec_chr = unique(GWAS$df_gwas$chr)
    # chr_lab_pos = c()
    # for (chr in vec_chr) {
    #     df_gwas_sub = GWAS$df_gwas[GWAS$df_gwas$chr == chr, ]
    #     chr_lab_pos = c(chr_lab_pos, median(df_gwas_sub$consecutive_pos))
    #     points(x=df_gwas_sub$consecutive_pos, y=df_gwas_sub$lod, pch=19, col=df_gwas_sub$colours_chr)
    # }
    ### Append peak labels
    if (nrow(df_gwas_peak_genes) > 0) {
        for (i in 1:nrow(df_gwas_peak_genes)) {
            # i = 1
            idx = (GWAS$df_gwas$chr %in% df_gwas_peak_genes$chr[i]) & (GWAS$df_gwas$pos %in% df_gwas_peak_genes$pos[i])
            idx_maxlod = GWAS$df_gwas$lod[idx]==max(GWAS$df_gwas$lod[idx]) ### use the allele with the highest LOD
            lod = GWAS$df_gwas$lod[idx][idx_maxlod]
            consecutive_pos = GWAS$df_gwas$consecutive_pos[idx][idx_maxlod]
            points(x=consecutive_pos, y=lod, pch=1, cex=4)
            text(x=consecutive_pos, y=lod, pos=2, label=df_gwas_peak_genes$gene[i])
        }
    }
    # ### Append Bonferroni threshold
    # abline(h=bonferroni_threshold, col="red", lty=2)
    # axis(side=1, at=chr_lab_pos, lab=vec_chr)
    # dev.off()
    # return(list(df_gwas_peak_genes=df_gwas_peak_genes,
    #             fname_manhattan_plot_with_peak_labels_png=fname_manhattan_plot_with_peak_labels_png))
    return(df_gwas_peak_genes)
}

MERGE_PHENOTYPE_WITH_GENOTYPE_DATA = function(df_phenotype, df_genotype) {
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
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
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
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
    # df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
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

MERGE_PHENOTYPE_AND_FST = function(df_phenotype, fname_fst="../res/fst.csv") {
    # df_phenotype = LOAD_PHENOTYPES(fname_phenotype="phenotype_data.csv", batch="all", phenotype_names="Glyphosate")
    df_fst = read.csv(fname_fst)
    colnames(df_fst)[1] = "X.Population"
    ### Find the populations with both Fst and phenotype data and sort
    df_tmp = merge(df_phenotype, df_fst, by="X.Population")
    vec_pop = colnames(df_tmp)[(ncol(df_phenotype)+1):ncol(df_tmp)]
    vec_pop = vec_pop[vec_pop %in% df_tmp$X.Population]
    idx_sort = order(vec_pop)
    ### Rebuild df_fst so that the columns are sorted
    mat_fst_tmp = as.matrix(df_fst[,2:ncol(df_fst)])
    df_fst = data.frame(X.Population=df_fst$X.Population[idx_sort], mat_fst_tmp[idx_sort, idx_sort])
    ### Merge and sort so the rows and columns are both sorted and the resulting Fst matrix is symmetric
    df_phe_fst = merge(df_phenotype, df_fst, by="X.Population")
    df_phe_fst = df_phe_fst[order(df_phe_fst$X.Population), ]
    mat_fst = as.matrix(df_phe_fst[, (ncol(df_phenotype)+1):ncol(df_phe_fst)])
    rownames(mat_fst) = df_phe_fst$X.Population
    return(list(df_phe_fst=df_phe_fst, mat_fst=mat_fst))
}


##################
### UNIT TESTS ###
##################
data_loading_and_merging = function() {
    test_that(
        "LOAD_PHENOTYPES", {
            print("LOAD_PHENOTYPES:")
            df1 = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            df2 = LOAD_PHENOTYPES(batch="SHOULD_ERROR")
            df3 = LOAD_PHENOTYPES(phenotype_names="SHOULD_ERROR_TOO")
            expect_equal(dim(df1), c(113, 5), tolerance=1e-7)
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
        "LOAD_KINSHIP_AND_SORT", {
            print("LOAD_KINSHIP_AND_SORT:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            population_names_sorted = df_phenotype$X.Population
            df1 = LOAD_KINSHIP_AND_SORT(population_names_sorted)
            df2 = LOAD_KINSHIP_AND_SORT(population_names_sorted, col_id_population_names=100000000)
            df3 = LOAD_KINSHIP_AND_SORT(population_names_sorted, fname_fst="../res/phenotype_data.csv")
            expect_equal(dim(df1), c(113, 113))
            expect_equal(df3, "ERROR!")
            expect_equal(df2, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_GENOTYPE_DATA", {
            print("MERGE_PHENOTYPE_WITH_GENOTYPE_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            df_genotype = LOAD_GENOTYPES(maf=0.2)
            df1 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=df_phenotype, df_genotype=df_genotype)
            df2 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=data.frame(), df_genotype=df_genotype)
            df3 = MERGE_PHENOTYPE_WITH_GENOTYPE_DATA(df_phenotype=data.frame(X.Population=NA), df_genotype=df_genotype)
            expect_equal(dim(df1), c(113, 321), tolerance=1e-7)
            expect_equal(df3, "ERROR!")
            expect_equal(df2, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA", {
            print("MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            raster_layers = LOAD_ENVIRONMENTAL_DATA()
            df1 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=df_phenotype, raster_layers=raster_layers)
            df2 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=data.frame(Longitude=NA), raster_layers=raster_layers)
            df3 = MERGE_PHENOTYPE_WITH_ENVIRONMENTAL_DATA(df_phenotype=data.frame(Latitude=NA), raster_layers=raster_layers)
            expect_equal(dim(df1), c(113, 23), tolerance=1e-7)
            expect_equal(df2, "ERROR!")
            expect_equal(df3, "ERROR!")
        }
    )

    test_that(
        "MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA", {
            print("MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            df_genotype = LOAD_GENOTYPES(maf=0.2)
            raster_layers = LOAD_ENVIRONMENTAL_DATA()
            df1 = MERGE_PHENOTYPE_WITH_GENOTYPE_AND_ENVIRONMENTAL_DATA(df_phenotype=df_phenotype, df_genotype=df_genotype, raster_layers=raster_layers)
            expect_equal(dim(df1), c(113, 339), tolerance=1e-7)
        }
    )

    test_that(
        "MERGE_PHENOTYPE_AND_FST", {
            print("MERGE_PHENOTYPE_AND_FST:")
            df_phenotype = LOAD_PHENOTYPES(phenotype_names=c("Glyphosate"))
            list_df1_mat1 = MERGE_PHENOTYPE_AND_FST(df_phenotype=df_phenotype)
            df1 = list_df1_mat1$df_phe_fst
            mat1 = list_df1_mat1$mat_fst
            expect_equal(dim(df1), c(113, 118), tolerance=1e-7)
            expect_equal(dim(mat1), c(113, 113), tolerance=1e-7)
        }
    )
}

# options(digits.secs=7)
# start_time = Sys.time()
# data_loading_and_merging()
# end_time = Sys.time()
# print(end_time - start_time)