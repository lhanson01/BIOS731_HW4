library(purrr)

n_chain <- 3
results_paths <- here::here("results", list.files(here::here("results")))

full_results <- vector(mode = "list", length = 0)
for(path in results_paths){
  load(path)  
  n <- length(results[[1]]$Gibbs$results_tib$labels[[1]])
  job_tib <- vector(mode="list", length = length(results))
  for(i in 1:length(results)){
    sim_tib_gibbs <- results[[i]]$Gibbs$results_tib %>%
      rename_with(~ paste0(.x, "_gibbs"), starts_with("mu"))
    vi_res <- as_tibble(t(c(results[[i]]$VI$means, results[[i]]$VI$comp_time) %>% 
      `names<-`(c("mu_1_VI","mu_2_VI","mu_3_VI","mu_4_VI","comp_time_VI" ))))
    temp_list <- vector(mode = "list", length = n_chain)
    for(j in 1:length(temp_list)){
      print(j)
      temp_list[[j]] <- bind_cols(sim_tib_gibbs[1,], vi_res) %>% mutate(chain = j)
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

# full_results_longer <- full_results_tib %>%
#   pivot_longer(
#     cols = matches("mu_(\\d+)_(.*)"),
#     names_to = c("mu_k","Method"),
#     names_pattern = "mu_(\\d+)_(.*)"
#   )

full_results_ag <- full_results_tib %>%
  select(-labels, -CIs) %>%
  group_by(n) %>%
  summarise(mu_1_mean_gibbs = mean(mu_1_gibbs),
            mu_2_mean_gibbs = mean(mu_2_gibbs),
            mu_3_mean_gibbs = mean(mu_3_gibbs),
            mu_4_mean_gibbs = mean(mu_4_gibbs),
            mu_1_mean_VI = mean(mu_1_VI),
            mu_2_mean_VI = mean(mu_2_VI),
            mu_3_mean_VI = mean(mu_3_VI),
            mu_4_mean_VI = mean(mu_4_VI),
            mu_1_SE_gibbs = sd(mu_1_gibbs),
            mu_2_SE_gibbs = sd(mu_2_gibbs),
            mu_3_SE_gibbs = sd(mu_3_gibbs),
            mu_4_SE_gibbs = sd(mu_4_gibbs),
            mu_1_mean_VI = mean(mu_1_VI),
            mu_2_mean_VI = mean(mu_2_VI),
            mu_3_mean_VI = mean(mu_3_VI),
            mu_4_mean_VI = mean(mu_4_VI),
            mu_1_SE_VI = 3*sd(mu_1_VI),
            mu_2_SE_VI = 3*sd(mu_2_VI),
            mu_3_SE_VI = 3*sd(mu_3_VI),
            mu_4_SE_VI = 3*sd(mu_4_VI),
            mu_1_coverage = mean(`mu_1_gibbs Covers`),
            mu_2_coverage = mean(`mu_2_gibbs Covers`),
            mu_3_coverage = mean(`mu_3_gibbs Covers`),
            mu_4_coverage = mean(`mu_4_gibbs Covers`),
            mean_comp_time_gibbs = mean(comp_time),
            mean_comp_time_VI = mean(comp_time_VI)) %>%
  mutate(bias_mu_1_gibbs = mu_1_mean_gibbs - mu_true[1],
         bias_mu_2_gibbs = mu_2_mean_gibbs - mu_true[2],
         bias_mu_3_gibbs = mu_3_mean_gibbs - mu_true[3],
         bias_mu_4_gibbs = mu_4_mean_gibbs - mu_true[4],
         bias_mu_1_VI = mu_1_mean_VI - mu_true[1],
         bias_mu_2_VI = mu_2_mean_VI - mu_true[2],
         bias_mu_3_VI = mu_3_mean_VI - mu_true[3],
         bias_mu_4_VI = mu_4_mean_VI - mu_true[4],
         mu_1_coverage_SE = sqrt(mu_1_coverage*(1-mu_1_coverage)/n),
         mu_2_coverage_SE = sqrt(mu_2_coverage*(1-mu_2_coverage)/n),
         mu_3_coverage_SE = sqrt(mu_1_coverage*(1-mu_3_coverage)/n),
         mu_4_coverage_SE = sqrt(mu_1_coverage*(1-mu_4_coverage)/n),
         )


bias_frame <- full_results_ag %>%
  select(n, matches("bias"), matches("mu_\\d+_SE")) %>%
  pivot_longer(cols = matches("bias_mu_\\d+"),
               names_to = c("mu_k","Method"),
               names_pattern = "bias_mu_(\\d+)_(.*)",
               values_to = "Bias") %>%
  mutate(n = factor(n))

bias_plot <- ggplot(bias_frame, aes(fill = mu_k, y = Bias, x = n)) +
         geom_col(position = "dodge", stat = "identity") + 
  facet_wrap(vars(Method)) +
  labs(title = "Bias by n and mu_k")

ggsave("./output/bias_plot.png",
       bias_plot)

coverage_frame <- full_results_ag %>%
  select(n, matches("coverage")) %>%
  pivot_longer(cols = matches("coverage_SE$"), names_to = "k", 
               names_pattern = "(.*)_coverage_SE", values_to = "coverage_SE") %>%
  pivot_longer(cols = matches("coverage$"), names_to = "k2", 
               names_pattern = "(.*)_coverage", values_to = "coverage") %>%
  filter(k == k2) %>%
  select(-k2) %>%
  mutate(LL = coverage - 1.96*coverage_SE,
         UL = coverage + 1.96*coverage_SE,
         n = factor(n))

cov_plot <- ggplot(coverage_frame, aes(fill = k, y = coverage, x = n)) +
  geom_col(position = "dodge", stat = "identity") +
  coord_cartesian(ylim=c(0.7,1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_errorbar(aes(ymin = LL, ymax = UL), 
                  position = position_dodge(0.9), width = 0.5, color = "grey25") +
  labs(title = "Coverage by mu and n")

ggsave("./output/cov_plot.png",
       cov_plot)




