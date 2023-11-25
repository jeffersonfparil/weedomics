dir = getwd()
source("../src/data_loading_and_merging.r")
setwd(dir)
dat = read.csv("gudmc-maf0.0_cov10_win10kb_slide10kb_minlocwin10.csv")
phen = LOAD_PHENOTYPES(fname_phenotype="phenotype_data.csv", batch="all", phenotype_names="Glyphosate")
threshold_resistant_population = 75
threshold_susceptible_population = 0
alpha = 0.05

vec_pop_a = c()
vec_pop_b = c()
vec_pop_a_resistant = c()
vec_pop_b_resistant = c()
vec_resistance_absolute_difference = c()
vec_n_tajima_troughs_in_pop_b = c()
vec_n_significantly_narrower_tajima_troughs_in_pop_b = c()
vec_n_significantly_wider_tajima_troughs_in_pop_b = c()
vec_tajima_troughs_with_significantly_lower_pairwise_fst = c()
vec_tajima_troughs_with_significantly_higher_pairwise_fst = c()
pb = txtProgressBar(min=0, max=length(unique(dat$pop_a))*length(unique(dat$pop_b)), style=3)
counter = 0
for (pop_a in unique(dat$pop_a)) {
    for (pop_b in unique(dat$pop_b)) {
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
        idx_sig_tajima_troughs = (df$tajima_d_pop_b < mean(df$tajima_d_pop_b)) & (df$tajima_width_one_tail_pval_pop_b <= alpha)
        tajima_width_minus_mean = df$tajima_width_pop_b - mean(df$tajima_width_pop_b, na.rm=TRUE)
        vec_pop_a = c(vec_pop_a, pop_a)
        vec_pop_b = c(vec_pop_b, pop_b)
        vec_pop_a_resistant = c(vec_pop_a_resistant, y_pop_a >= threshold_resistant_population)
        vec_pop_b_resistant = c(vec_pop_b_resistant, y_pop_b >= threshold_resistant_population)
        vec_resistance_absolute_difference = c(vec_resistance_absolute_difference, abs(y_pop_a - y_pop_b))
        vec_n_tajima_troughs_in_pop_b = c(vec_n_tajima_troughs_in_pop_b, sum(idx_sig_tajima_troughs))
        vec_n_significantly_narrower_tajima_troughs_in_pop_b = c(vec_n_significantly_narrower_tajima_troughs_in_pop_b, sum(idx_sig_tajima_troughs & (tajima_width_minus_mean <= 0)))
        vec_n_significantly_wider_tajima_troughs_in_pop_b = c(vec_n_significantly_wider_tajima_troughs_in_pop_b, sum(idx_sig_tajima_troughs & (tajima_width_minus_mean > 0)))
        vec_tajima_troughs_with_significantly_lower_pairwise_fst = c(vec_tajima_troughs_with_significantly_lower_pairwise_fst, sum(idx_sig_tajima_troughs & (df$fst_delta_one_tail_pval <= alpha) & (df$fst_delta < 0)))
        vec_tajima_troughs_with_significantly_higher_pairwise_fst = c(vec_tajima_troughs_with_significantly_higher_pairwise_fst, sum(idx_sig_tajima_troughs & (df$fst_delta_one_tail_pval <= alpha) & (df$fst_delta > 0)))
    }
}
close(pb)
df_counts = data.frame(pop_a=vec_pop_a, pop_b=vec_pop_b, pop_a_resistant=vec_pop_a_resistant, pop_b_resistant=vec_pop_b_resistant, resistance_absolute_difference=vec_resistance_absolute_difference, n_tajima_troughs_in_pop_b=vec_n_tajima_troughs_in_pop_b, n_significantly_narrower_tajima_troughs_in_pop_b=vec_n_significantly_narrower_tajima_troughs_in_pop_b, n_significantly_wider_tajima_troughs_in_pop_b=vec_n_significantly_wider_tajima_troughs_in_pop_b, tajima_troughs_with_significantly_lower_pairwise_fst=vec_tajima_troughs_with_significantly_lower_pairwise_fst, tajima_troughs_with_significantly_higher_pairwise_fst=vec_tajima_troughs_with_significantly_higher_pairwise_fst)

### Convergent evolution
idx = (df_counts$pop_a_resistant==TRUE) & (df_counts$pop_b_resistant==TRUE)
df_counts_sub = df_counts[idx, ]
mean_prob_de_novo = mean(df_counts_sub$n_significantly_narrower_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_de_novo = sd(df_counts_sub$n_significantly_narrower_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_standing_genetic_variation = mean(df_counts_sub$n_significantly_wider_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_standing_genetic_variation = sd(df_counts_sub$n_significantly_wider_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_independent_emergence = mean(df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_independent_emergence = sd(df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_migration_from_one_pop = mean(df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_migration_from_one_pop = sd(df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_shared_ancestry = mean(1-((df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst + df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst) / df_counts_sub$n_tajima_troughs_in_pop_b))
sd_prob_shared_ancestry = sd(1-((df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst + df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst) / df_counts_sub$n_tajima_troughs_in_pop_b))
convergent_evolution = list(mean_prob_de_novo=mean_prob_de_novo, sd_prob_de_novo=sd_prob_de_novo, mean_prob_standing_genetic_variation=mean_prob_standing_genetic_variation, sd_prob_standing_genetic_variation=sd_prob_standing_genetic_variation, mean_prob_independent_emergence=mean_prob_independent_emergence, sd_prob_independent_emergence=sd_prob_independent_emergence, mean_prob_migration_from_one_pop=mean_prob_migration_from_one_pop, sd_prob_migration_from_one_pop=sd_prob_migration_from_one_pop, mean_prob_shared_ancestry=mean_prob_shared_ancestry, sd_prob_shared_ancestry=sd_prob_shared_ancestry)

### Divergent evolution
idx = (df_counts$pop_a_resistant==FALSE) & (df_counts$pop_b_resistant==TRUE)
df_counts_sub = df_counts[idx, ]
mean_prob_de_novo = mean(df_counts_sub$n_significantly_narrower_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_de_novo = sd(df_counts_sub$n_significantly_narrower_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_standing_genetic_variation = mean(df_counts_sub$n_significantly_wider_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_standing_genetic_variation = sd(df_counts_sub$n_significantly_wider_tajima_troughs_in_pop_b / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_independent_emergence = mean(df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_independent_emergence = sd(df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_migration_from_one_pop = mean(df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
sd_prob_migration_from_one_pop = sd(df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst / df_counts_sub$n_tajima_troughs_in_pop_b)
mean_prob_shared_ancestry = mean(1-((df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst + df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst) / df_counts_sub$n_tajima_troughs_in_pop_b))
sd_prob_shared_ancestry = sd(1-((df_counts_sub$tajima_troughs_with_significantly_lower_pairwise_fst + df_counts_sub$tajima_troughs_with_significantly_higher_pairwise_fst) / df_counts_sub$n_tajima_troughs_in_pop_b))
divergent_evolution = list(mean_prob_de_novo=mean_prob_de_novo, sd_prob_de_novo=sd_prob_de_novo, mean_prob_standing_genetic_variation=mean_prob_standing_genetic_variation, sd_prob_standing_genetic_variation=sd_prob_standing_genetic_variation, mean_prob_independent_emergence=mean_prob_independent_emergence, sd_prob_independent_emergence=sd_prob_independent_emergence, mean_prob_migration_from_one_pop=mean_prob_migration_from_one_pop, sd_prob_migration_from_one_pop=sd_prob_migration_from_one_pop, mean_prob_shared_ancestry=mean_prob_shared_ancestry, sd_prob_shared_ancestry=sd_prob_shared_ancestry)

out = data.frame(evolution=c("convergent", "divergent"), 
                 rbind(unlist(convergent_evolution), unlist(divergent_evolution)))

write.table(out, file="gudmc_summary_results.csv", row.names=F, col.names=T, sep=",", quote=F)
