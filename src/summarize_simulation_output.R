library(tidyverse)
library(data.table)

results <- read.table("mnt/mpox_seqs/simulations/pure_samples/aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

# read lineages file created using nextclade
lineages <- fread("mnt/mpox_seqs/lineages_B1.tsv",header = FALSE)
# df operations
results_comb_names <- results %>% separate(lineages, into = c("obs_lin1","obs_lin2"), sep = " ")
results_comb_names <- results_comb_names %>% separate(abundances, into = c("obs_abun1","obs_abun2"), sep = " ")
results_comb_names <-results_comb_names[,c(1,3,4,5,6)]
results_comb_names <- results_comb_names %>% 
  separate(X, into = as.character(1:4),sep="_")
results_comb_names$exp_abun1 <- as.numeric(str_remove(results_comb_names$X, ".vcf"))
results_comb_names<-results_comb_names[,-4]
colnames(results_comb_names) <- c("isolate1","exp_abun1","isolate2","obs_lin1","obs_lin2","obs_abun1","obs_abun2","exp_abun2")
# Fix order for abundances )
results_comb_names <- results_comb_names %>% mutate(real_obs_abun1 = case_when(obs_abun1 > obs_abun2~ obs_abun1,
                                                                               obs_abun2 > obs_abun1 ~ obs_abun2, obs_abun1 == obs_abun2~ obs_abun1))
results_comb_names <- results_comb_names %>% mutate(real_obs_abun2 = case_when(obs_abun1 < obs_abun2~ obs_abun1,
                                                                               obs_abun2 < obs_abun1 ~ obs_abun2, obs_abun1 == obs_abun2~ obs_abun1))
results_comb_names <- results_comb_names %>% mutate(exp_abun1 = as.numeric(exp_abun1))
results_comb_names <- results_comb_names %>% mutate(exp_abun2 = as.numeric(exp_abun2))
results_comb_names <- results_comb_names %>% mutate(real_exp_abun1 = case_when(exp_abun1 > exp_abun2~ exp_abun1,
                                                                               exp_abun2 > exp_abun1 ~ exp_abun2, exp_abun1 == exp_abun2~ exp_abun1))
results_comb_names <- results_comb_names %>% mutate(real_exp_abun2 = case_when(exp_abun1 < exp_abun2~ exp_abun1,
                                                                               exp_abun2 < exp_abun1 ~ exp_abun2, exp_abun1 == exp_abun2~ exp_abun1))
results_comb_names <- results_comb_names %>% mutate(res1 = abs(as.numeric(real_obs_abun1) - as.numeric(real_exp_abun1)))
results_comb_names <- results_comb_names %>% mutate(res2 = abs(as.numeric(real_obs_abun2) - as.numeric(real_exp_abun2)))

results$isolate <-str_remove(results$X, ".vcf") 
 lineages  %>%
  inner_join(results, by =c("V1" ="isolate")) %>% select(V1,V2,lineages,abundances)%>%
   write_csv("Downloads/B1_pure_simulation.csv")
 
colnames(results_comb_names)[1:2] <- c("isolate1","exp_lin1")

results_comb_names <- lineages %>%
  inner_join(results_comb_names, by =c("V1" ="isolate2")) 
colnames(results_comb_names)[1:2] <- c("isolate2","exp_lin2")
# reorder, change column names
results_comb_names <- results_comb_names %>%
  select(isolate1,isolate2,exp_lin1,exp_lin2,obs_lin1,obs_lin2,
         real_exp_abun1,real_exp_abun2,real_obs_abun1,real_obs_abun2,res1,res2)
colnames(results_comb_names)[7:10] <-c("exp_abun1","exp_abun2","obs_abun1","obs_abun2")
results_comb_names %>% write_csv("Downloads/B1_rep_simulation_output.csv")
