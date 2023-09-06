dirname_src = utils::getSrcDirectory(function(){})[1]
if (is.na(dirname_src)) {
    dirname_src = "."
}
setwd(dirname_src)
source("./data_loading_and_merging.r")

PLOT_TAJIMA = function(df_tajima_pheno, pop0, pop1, 
                       chromosomome_identifier="NC_",
                       pop0_lab=paste0(pop0," (Glyphosate-susceptible)"),
                       pop1_lab=paste0(pop1, " (Glyphosate-resistant)"),
                       prefix="tajima_d-maf0.0_cov10_win10kb",
                       vec_col=c("#1f78b4", "#b2df8a"),
                       plot_svg=TRUE) {
    ### Extract the Tajima's D of each population
    idx_col = grepl("Window", colnames(df_tajima_pheno))
    y0 = t(df_tajima_pheno[df_tajima_pheno$X.Population==pop0, idx_col])[,1]
    y1 = t(df_tajima_pheno[df_tajima_pheno$X.Population==pop1, idx_col])[,1]
    ### Normalise Tajima's D
    y0 = scale(y0, center=TRUE, scale=TRUE)
    y1 = scale(y1, center=TRUE, scale=TRUE)
    ### Identify and keep only the windows within chromosomes
    vec_windows = colnames(df_tajima_pheno)[idx_col]
    vec_windows = gsub(chromosomome_identifier, gsub("_", "", chromosomome_identifier), vec_windows) ### Remove underscores because we will be parsing window names using "_" as delimiter
    chromosomome_identifier_no_ = gsub("_", "", chromosomome_identifier)
    idx_chr = grepl(chromosomome_identifier_no_, vec_windows)
    vec_windows = vec_windows[idx_chr]
    y0 = y0[idx_chr]
    y1 = y1[idx_chr]
    ### Merge into a data frame
    W = matrix(unlist(strsplit(vec_windows, "_")), byrow=TRUE, ncol=3)
    Tajimas_D = data.frame(X.chr=gsub("Window.", "", W[,1]),
                           pos_ini=as.numeric(W[,2]), pos_fin=as.numeric(W[,3]),
                           pos=round(as.numeric(W[,2]) + ((as.numeric(W[,3])-as.numeric(W[,2]))/2)),
                           y0=y0, y1=y1)
    ### Transform chromosome names wuth the original identifier
    Tajimas_D$X.chr = gsub(chromosomome_identifier_no_, chromosomome_identifier, Tajimas_D$X.chr)
    ### Add consecutive locus position ids and colours corresponding to each chromosome
    Tajimas_D$consecutive_pos = NA
    Tajimas_D$colours_chr = NA
    vec_chr = sort(unique(Tajimas_D$X.chr))
    vec_col = rep(vec_col, times=ceiling(length(vec_chr)/length(vec_col)))
    minimum = min(Tajimas_D$pos[Tajimas_D$X.chr == vec_chr[1]], na.rm=TRUE)
    Tajimas_D$consecutive_pos[Tajimas_D$X.chr == vec_chr[1]] = Tajimas_D$pos[Tajimas_D$X.chr == vec_chr[1]] - minimum
    Tajimas_D$colours_chr[Tajimas_D$X.chr == vec_chr[1]] = vec_col[1]
    maximum = max(Tajimas_D$consecutive_pos[Tajimas_D$X.chr == vec_chr[1]], na.rm=TRUE)
    for (i in 2:length(vec_chr)) {
        # i = 2
        minimum = min(Tajimas_D$pos[Tajimas_D$X.chr == vec_chr[i]], na.rm=TRUE)
        Tajimas_D$consecutive_pos[Tajimas_D$X.chr == vec_chr[i]] = (Tajimas_D$pos[Tajimas_D$X.chr == vec_chr[i]] - minimum) + maximum + 1
        Tajimas_D$colours_chr[Tajimas_D$X.chr == vec_chr[i]] = vec_col[i]
        maximum = max(Tajimas_D$consecutive_pos[Tajimas_D$X.chr == vec_chr[i]], na.rm=TRUE)
    }
    ### Find the chromosome name label positions
    vec_chr_xlab_pos = aggregate(consecutive_pos ~ X.chr, data=Tajimas_D, FUN=median)
    ### Calculte the Tajima's D difference between the pair of populations but keep only the Tajima's D on the pop1 that is at least -2.00 (i.e. significant-ish selective sweep)
    Tajimas_D$diff = Tajimas_D$y0 - Tajimas_D$y1
    Tajimas_D$diff[Tajimas_D$y1 > -2] = 0
    ### Plot
    if (plot_svg==TRUE) {
        fname_svg = paste0(prefix, "-", pop0, "-", pop1, ".svg")
        print(fname_svg)
        svg(fname_svg, width=11, height=7)
        layout(matrix(seq(1,3), ncol=1))
        plot(x=Tajimas_D$consecutive_pos, Tajimas_D$y0, type="n", col=Tajimas_D$colours_chr, main=pop0_lab, xaxt="n", xlab="Genome", ylab="Tajima's D"); grid()
        axis(1, at=vec_chr_xlab_pos$consecutive_pos, lab=vec_chr_xlab_pos$X.chr)
        for (i in 1:nrow(vec_chr_xlab_pos)) {
            # i = 1
            idx = Tajimas_D$X.chr == vec_chr_xlab_pos$X.chr[i]
            lines(x=Tajimas_D$consecutive_pos[idx], y=Tajimas_D$y0[idx], col=Tajimas_D$colours_chr[idx], lwd=2)
        }
        plot(x=Tajimas_D$consecutive_pos, Tajimas_D$y1, type="n", col=Tajimas_D$colours_chr, main=pop1_lab, xaxt="n", xlab="Genome", ylab="Tajima's D"); grid()
        axis(1, at=vec_chr_xlab_pos$consecutive_pos, lab=vec_chr_xlab_pos$X.chr)
        for (i in 1:nrow(vec_chr_xlab_pos)) {
            # i = 1
            idx = Tajimas_D$X.chr == vec_chr_xlab_pos$X.chr[i]
            lines(x=Tajimas_D$consecutive_pos[idx], y=Tajimas_D$y1[idx], col=Tajimas_D$colours_chr[idx], lwd=2)
        }
        plot(x=Tajimas_D$consecutive_pos, Tajimas_D$diff, type="n", col=Tajimas_D$colours_chr, main=paste0(pop0, " - ", pop1), xaxt="n", xlab="Genome", ylab="F( Tajima's D )"); grid()
        axis(1, at=vec_chr_xlab_pos$consecutive_pos, lab=vec_chr_xlab_pos$X.chr)
        for (i in 1:nrow(vec_chr_xlab_pos)) {
            # i = 1
            idx = Tajimas_D$X.chr == vec_chr_xlab_pos$X.chr[i]
            lines(x=Tajimas_D$consecutive_pos[idx], y=Tajimas_D$diff[idx], col=Tajimas_D$colours_chr[idx], lwd=2)
        }
        dev.off()
    }
    return(Tajimas_D)
}


### DOES NOT HAVE A TEST SUITE YET! 20230905
IDENTIFY_GENES_WITHIN_PEAKS = function(Tajimas_D, fname_genome_annotation="../res/genome_annotation.gff",
                                       margin_bp=1e5) {
    gff = read.delim(fname_genome_annotation, sep="\t", header=FALSE, comment.char="#")
    peaks = Tajimas_D[Tajimas_D$diff > 0, ]
    peaks$gene = NA
    for (i in 1:nrow(peaks)) {
        # i = 2
        pos_chr = peaks$X.chr[i]
        pos_ini = peaks$pos_ini[i]
        pos_fin = peaks$pos_fin[i]
        idx = (gff$V1 == pos_chr) & 
            (gff$V3 == "gene") &
            (gff$V4-margin_bp <= pos_ini) &
            (gff$V5+margin_bp >= pos_fin)
        if (sum(idx) > 1) {
            print(i)
            print(W_peaks[i, ])
            print(grepl("product", gff$V9[idx]))
            peaks$gene[i] = gff$V9[idx][1]
        }
    }
    return(peaks)
}

##################
### UNIT TESTS ###
##################
popgen = function() {
    seed = 42069
    set.seed(seed)
    n_windows = 50
    window_size = 10000
    vec_pop = c("pop0", "pop1")
    vec_chr_scaff = c("chr_A", "chr_B", "chr_C", "scaff_x", "scaff_y")
    vec_chr = paste0("Window.", rep(vec_chr_scaff, each=n_windows/(length(vec_chr_scaff))))
    pos_ini = (0:(n_windows-1) * window_size) + 1
    pos_fin = 1:n_windows * window_size
    vec_windows = paste(vec_chr, pos_ini, pos_fin, sep="_")
    y0 = rnorm(n_windows)
    y1 = rnorm(n_windows)
    df_tajima_pheno = data.frame(vec_pop, rbind(y0, y1))
    colnames(df_tajima_pheno) = c("X.Population", vec_windows)
    str(df_tajima_pheno)

    test_that(
        "PLOT_TAJIMA", {
            print("PLOT_TAJIMA:")
            set.seed(seed)
            out = PLOT_TAJIMA(df_tajima_pheno, "pop0", "pop1", 
                              chromosomome_identifier="chr_",
                              prefix="test",
                              plot_svg=FALSE)
            expect_equal(nrow(out), 30, tolerance=1e-7)
            expect_equal(ncol(out),  8, tolerance=1e-7)
            expect_equal(out$diff[1],  2.242691, tolerance=1e-6)
            expect_equal(out$diff[2],  0.000000, tolerance=1e-6)
            expect_equal(out$diff[3],  0.000000, tolerance=1e-6)
            expect_equal(out$diff[4],  0.000000, tolerance=1e-6)
            expect_equal(out$diff[5],  2.367536, tolerance=1e-6)
        }
    )
}

# options(digits.secs=7)
# start_time = Sys.time()
# popgen()
# end_time = Sys.time()
# print(end_time - start_time)