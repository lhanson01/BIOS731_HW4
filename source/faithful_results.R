mu_chains_full <- data.frame(mu_1_chain_1 = gibbs_results$mu_hist[[1]][,1],
                             mu_1_chain_2 = gibbs_results$mu_hist[[2]][,1],
                             mu_1_chain_3 = gibbs_results$mu_hist[[3]][,1],
                             mu_2_chain_1 = gibbs_results$mu_hist[[1]][,2],
                             mu_2_chain_2 = gibbs_results$mu_hist[[2]][,2],
                             mu_2_chain_3 = gibbs_results$mu_hist[[3]][,2]) %>%
  pivot_longer(
    cols = matches("mu_(\\d+)_chain_(\\d+)*"),
    names_to = c("mu","Chain"),
    names_pattern = "mu_(\\d+)_chain_(\\d+)") 


mu_chains_average <- mu_chains_full %>%
  mutate(average_var = ceiling((seq(0.1,10000, length.out = 10000)/10))) %>%
  group_by(average_var) %>%
  summarise(mu_1_ch1_short = mean(mu_1_chain_1),
            mu_1_ch2_short = mean(mu_1_chain_2),
            mu_1_ch3_short = mean(mu_1_chain_3),
            mu_2_ch1_short = mean(mu_2_chain_1),
            mu_2_ch2_short = mean(mu_2_chain_2),
            mu_2_ch3_short = mean(mu_2_chain_3)) %>%
  mutate(x = seq_along(mu_1_ch1_short)) %>%
  pivot_longer(
    cols = matches("mu_(\\d+)_ch(\\d+)*"),
    names_to = c("mu","Chain"),
    names_pattern = "mu_(\\d+)_ch(\\d+)") 

mu_1_chains <- mu_chains %>% filter(mu == 1)
mu_2_chains <- mu_chains %>% filter(mu == 2)

mu_1_trace_plot <- ggplot(data = mu_1_chains, aes(x = x, y = value, color = Chain) ) +
  geom_line(linewidth = 0.1) +
  facet_wrap(vars(Chain)) +
  labs(title = "Mu_1 Trace Plot Across 3 Chains") + 
  xlab("Sample (averaged in groups of 10)")

mu_2_trace_plot <- ggplot(data = mu_2_chains, aes(x = x, y = value, color = Chain)) +
  geom_line(linewidth = 0.1) +
  facet_wrap(vars(Chain)) +
  labs(title = "Mu_2 Trace Plot Across 3 Chains") + 
  xlab("Sample (averaged in groups of 10)")

mu_1_lag_plot <- forecast::ggAcf(
  mu_chains_full %>% filter(mu == 1, Chain == 1) %>% select(value),
  type = "correlation") + 
  labs(title = "Mu 1 Chain 1 lag plot") + ylab("Correlation")

mu_2_lag_plot <- forecast::ggAcf(
  mu_chains_full %>% filter(mu == 2, Chain == 1) %>% select(value),
  type = "correlation") + 
  labs(title = "Mu 2 Chain 1 lag plot") + ylab("Correlation")


chain1_ess <- effectiveSize(gibbs_results$mu_hist[[1]])
chain2_ess <- effectiveSize(gibbs_results$mu_hist[[2]])
chain3_ess <- effectiveSize(gibbs_results$mu_hist[[1]])

ggsave(
  mu_1_trace_plot,
  file = "./output/mu_1_trace_plot.png",
  width = 12,
  height = 5
)

ggsave(
  mu_2_trace_plot,
  file = "./output/mu_2_trace_plot.png",
  width = 12,
  height = 5
)

ggsave(
  mu_1_lag_plot,
  file = "./output/mu_1_lag_plot.png"
)

ggsave(
  mu_2_lag_plot,
  file = "./output/mu_2_lag_plot.png"
)
