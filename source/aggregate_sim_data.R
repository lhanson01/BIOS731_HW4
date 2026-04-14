load("./results/n_100_sim_1-50.RDA")

n_chain <- 3
results_paths <- here::here("results", list.files(here::here("results")))

full_results <- vector(mode = "list", length = 0)
for(path in results_paths){
  load(path)  
  n <- length(results[[1]]$Gibbs$results_tib$labels[[1]])
  job_tib <- vector(mode="list", length = length(results))
  for(i in 1:length(results)){
    sim_tib <- results[[i]]$Gibbs$results_tib
    temp_list <- vector(mode = "list", length = n_chain)
    for(j in 1:length(temp_list)){
      temp_list[[j]] <- sim_tib[1,] %>% mutate(chain = j)
    }
    job_tib[[i]] <- bind_rows(temp_list) %>% mutate(n = n)
  }
full_results[[path]] <- bind_rows(job_tib)
}

full_results_tib <- bind_rows(full_results)

mu_col <- names(full_results_tib)[grep("mu", names(full_results_tib))]
mu_true <- c(0, 5, 10, 20)
for(k in 1:4){
  full_results_tib <- full_results_tib %>% 
    mutate(!!paste(mu_col[k], "Covers") := unlist(map(CIs, 
      function(CI_mat){
        covers <- unlist(mu_true[k] >= CI_mat[1,k] & mu_true[k] <= CI_mat[2,k])
      }
      )) 
    )
}

full_results_ag <- full_results_tib %>%
  select(-labels, -CIs) %>%
  group_by(n) %>%
  summarise(mu_1_mean = mean(mu_1),
            mu_2_mean = mean(mu_2),
            mu_3_mean = mean(mu_3),
            mu_4_mean = mean(mu_4),
            mu_1_coverage = mean(`mu_1 Covers`),
            mu_2_coverage = mean(`mu_2 Covers`),
            mu_3_coverage = mean(`mu_3 Covers`),
            mu_4_coverage = mean(`mu_4 Covers`),
            mean_comp_time = mean(comp_time)) %>%
  mutate(bias_mu_1 = mu_1_mean - mu_true[1],
         bias_mu_2 = mu_2_mean - mu_true[2],
         bias_mu_3 = mu_3_mean - mu_true[3],
         bias_mu_4 = mu_4_mean - mu_true[4])
  


