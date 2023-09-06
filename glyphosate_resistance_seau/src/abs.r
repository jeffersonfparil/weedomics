### Approximate Bayesian Computation
library(EasyABC)
prior <- list(c("unif", 1e-9, 1e-6))
observed <- 262
ABC_SLiM <- ABC_sequential(method="Lenormand",
                           use_seed=TRUE,
                           model=runSLiM,
                           prior=prior,
                           summary_stat_target=observed,
                           nb_simul=1000)



### Working directory and input files
dir="/data-weedomics-1/weedomics/glyphosate_resistance_seau/res"
setwd(dir)
fname_heterozygosity_LOC124669605_5pop_csv = "heterozygosity_LOC124669605_5pop.csv"

### Fxed population genetic paratmeters
### Mutation rate: 9.01 x 10^(-9) - estimates from Oryza sativa (not using Arabipsis thaliana estimates of 4.35 x 10(-9)) from Wang et a, 2019
u = 9.01e-9

df = read.csv(fname_heterozygosity_LOC124669605_5pop_csv)
df = df[c(2,3,1,5,4), ] ### sort manually by resistance levels
vec_pi = df$Mean_across_windows
vec_Ne = vec_pi / (4*u)

df$Ne = vec_Ne

df_sim = read.csv("slim_log.tmp")
svg("test.svg")
plot(x=df_sim$cycle, y=df_sim$FST_p1xp2, typ="l")
grid()
dev.off()
