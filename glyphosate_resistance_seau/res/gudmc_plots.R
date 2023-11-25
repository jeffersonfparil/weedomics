dir = getwd()
source("../src/data_loading_and_merging.r")
setwd(dir)
dat = read.csv("gudmc-maf0.0_cov10_win10kb_slide10kb_minlocwin10.csv")
phen = LOAD_PHENOTYPES(fname_phenotype="phenotype_data.csv", batch="all", phenotype_names="Glyphosate")
threshold_resistant_population = 75
threshold_susceptible_population = 0
alpha = 0.05

vec_pop_a = sort(unique(dat$pop_a))
vec_pop_b = sort(unique(dat$pop_b))

pb = txtProgressBar(min=0, max=length(vec_pop_a)*length(vec_pop_b), style=3)
counter = 0
for (pop_a in vec_pop_a) {
    for (pop_b in vec_pop_b) {
        # pop_a="ACC062"; pop_b="ACC041"
        counter = counter + 1
        setTxtProgressBar(pb, counter)
        # pop_a = unique(dat$pop_a)[1]; pop_b = unique(dat$pop_b)[93]
        if (pop_a == pop_b) {
            next
        }
        idx_phen_pop_a = phen$X.Population == pop_a
        idx_phen_pop_b = phen$X.Population == pop_b
        if ((sum(idx_phen_pop_a)==0) | (sum(idx_phen_pop_b)==0)) {
            next
        }
        df = dat[((dat$pop_a==pop_a) & (dat$pop_b==pop_b)), ]
        y_pop_a = phen$Glyphosate[idx_phen_pop_a]
        y_pop_b = phen$Glyphosate[idx_phen_pop_b]
        ### Use only population pairs where pop_b is resistant and the difference between resistance is at least 50% (divergent evolution)
        ### ... or if both populations are resistant (convergent evolution).
        if ((y_pop_b < threshold_resistant_population) & (y_pop_a < threshold_resistant_population)) {
            next
        }
        df = df[grepl("^NC", df$chr), ]
        vec_chr = unique(df$chr)
        vec_chr_pos = c()
        vec_consecutive_pos = c(0)
        for (chr in vec_chr) {
            # chr = vec_chr[1]
            idx_chr = which(df$chr == chr)
            vec_pos = df$pos_ini[idx_chr] + (df$pos_fin[idx_chr] - df$pos_ini[idx_chr])/2
            vec_pos = vec_pos - min(vec_pos) + 1 + tail(vec_consecutive_pos, n=1)
            vec_chr_pos = c(vec_chr_pos, median(vec_pos))
            vec_consecutive_pos = c(vec_consecutive_pos, vec_pos)
        }
        vec_consecutive_pos = vec_consecutive_pos[-1]
        df$vec_consecutive_pos = vec_consecutive_pos
        x = df$vec_consecutive_pos
        lod_threshold = -log10(0.01)
        tajima = df$tajima_d_pop_b
        dfst = df$fst_delta
        tajima_pval = df$tajima_width_one_tail_pval_pop_b
        dfst_pval = df$fst_delta_one_tail_pval
        tajima_idx = which((-log10(tajima_pval) >= lod_threshold) & (tajima != 0))
        dfst_idx = which((-log10(dfst_pval) >= lod_threshold) & (dfst != 0))
        overlap_idx = which(((-log10(tajima_pval) >= lod_threshold) & (tajima != 0)) & ((-log10(dfst_pval) >= lod_threshold) & (dfst != 0)))
        labels = paste0(df$chr, "_", df$pos_ini, "-", df$pos_fin)
        svg(paste0("gudmc_plots-", pop_b, "_x_", pop_a, ".svg"), width=14, height=7)
        layout(matrix(1:2, nrow=2, ncol=1))
        plot(x=x, y=tajima, xaxt="n", main=paste0("Tajima's D (", pop_b, "[R])"), xlab="", ylab="", type="l", las=2); grid()
        axis(1, at=vec_chr_pos, label=vec_chr)
        ### significantly deviating Tajima's D width compared with the genomewide mean
        tajima_idx_width_sigsmall = df$tajima_width_pop_b[tajima_idx] < mean(df$tajima_width_pop_b)
        tajima_idx_width_siglarge = df$tajima_width_pop_b[tajima_idx] > mean(df$tajima_width_pop_b)
        ### de novo mutation (green points)
        points(x=x[tajima_idx][tajima_idx_width_sigsmall], y=tajima[tajima_idx][tajima_idx_width_sigsmall], cex=2, col="green")
        ### standing genetic variation (red points)
        points(x=x[tajima_idx][tajima_idx_width_siglarge], y=tajima[tajima_idx][tajima_idx_width_siglarge], cex=2, col="red")
        ### troughs and peaks with significantly deviated Fst in the plot below
        pop_a_res = "[S]"
        if (y_pop_a >= threshold_resistant_population) {
            pop_a_res = "[R]"
        }
        plot(x=x, y=dfst, xaxt="n", main=paste0("Fst deviation from genomewide mean\n(", pop_a, pop_a_res, " vs ", pop_b, "[R])"), xlab="", ylab="", type="l", las=2); grid()
        axis(1, at=vec_chr_pos, label=vec_chr)
        if (length(overlap_idx) > 0) {
            text(x=x[overlap_idx], y=dfst[overlap_idx], cex=0.5, pos=2, label=labels[overlap_idx])
            if (dfst[overlap_idx] > 0) {
                points(x=x[overlap_idx], y=dfst[overlap_idx], col="blue")
            } else if (dfst[overlap_idx] < 0) {
                points(x=x[overlap_idx], y=dfst[overlap_idx], col="orange")
            }
        }
        legend("topright", legend=c("de novo mutation (narrow Tajima'D)",
                                    "standing genetic variabtion (wide Tajima's D)",
                                    "independent emergence (higher than the mean Fst)",
                                    "migration (lower than the mean Fst)",
                                    "shared ancestry (~mean Fst, no colour)"),
                           fill=c("green",
                                  "red",
                                  "blue",
                                  "orange",
                                  "white"), bty="n", cex=0.75)
        dev.off()
    }
}
close(pb)
