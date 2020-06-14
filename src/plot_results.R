experiment_name = "main_simulation"
metrics_filename = sprintf("../results/%s/metrics_%s.tsv", experiment_name, experiment_name)
df_metrics = read_tsv(metrics_filename, col_types = "cccddidi", comment = "#")

df_metrics = df_metrics %>% 
  mutate(method = factor(method, labels = c("LOCO CRT", "Knockoffs", "Marginal CRT", "HRT")),
         rho = factor(rho, levels = c(0.2, 0.5, 0.8), labels = c("rho = 0.2", "rho = 0.5", "rho = 0.8")),
         k = factor(k, levels = c(25, 50, 100), labels = c("k = 25", "k = 50", "k = 100")))

# Power of FDR-controlling methods
p = df_metrics %>% 
  group_by(method, error_rate, metric, A, k, rho) %>%
  summarise(value_mean = mean(value), value_sd = sd(value), value_se = sd(value)/sqrt(n())) %>%
  ungroup() %>%
  filter(error_rate == "FDR", metric == "power") %>% 
  ggplot(aes(x = A, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("blue", "forestgreen", "dodgerblue", "darkgoldenrod")) +
  facet_grid(k ~ rho, scales = "free") + theme_bw() + 
  xlab("Non-null coefficient magnitude") + ylab("Power") + 
  ggtitle("Controlling the false discovery rate") + 
  theme(legend.title = element_blank(), legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
plot(p)

# Power of FWER-controlling methods
p = df_metrics%>% 
  group_by(method, error_rate, metric, A, k, rho) %>%
  summarise(value_mean = mean(value), value_sd = sd(value), value_se = sd(value)/sqrt(n())) %>%
  ungroup() %>%
  filter(error_rate == "FWER", metric == "power", method != "Knockoffs") %>% 
  ggplot(aes(x = A, y = value_mean, group = method, colour = method)) + 
  geom_line() + geom_point() + 
  scale_colour_manual(values = c("blue", "dodgerblue", "darkgoldenrod")) +
  facet_grid(k ~ rho, scales = "free") + theme_bw() + xlab("Non-null coefficient magnitude") + ylab("Power") + 
  ggtitle("Controlling the familywise error rate") + 
  theme(legend.title = element_blank(), legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
plot(p)

# Knockoff variability
experiment_name = "knockoffs_variability"
rejections_filename = sprintf("../results/%s/rejections_%s.tsv", experiment_name, experiment_name)
df_rejections = read_tsv(metrics_filename, col_types = "cccddidi", comment = "#")

p = 500
buffer = ceiling(p/200)
df_rejections = df_rejections %>%
  mutate(nonnull = 
           k == 25 & variable %in% round(seq(buffer, p-buffer, length.out = 25)) |
           k == 50 & variable %in% round(seq(buffer, p-buffer, length.out = 50)) | 
           k == 100 & variable %in% round(seq(buffer, p-buffer, length.out = 100)))

get_Jaccard_and_power = function(df,...){
  num_rejections = integer(num_realizations)
  names(num_rejections) = realization_names
  total_rejections = df %>% pull(method) %>% table()
  num_rejections[names(total_rejections)] = total_rejections
  num_intersections = matrix(0, num_realizations, num_realizations, 
                             dimnames = list(realization_names, realization_names))
  power = numeric(num_realizations)
  names(power) = realization_names
  df_power = df %>% group_by(method) %>% summarise(power = sum(nonnull)/unique(k))
  power[df_power$method] = df_power$power
  index_start = 1
  while(index_start < nrow(df)){
    index_end = index_start + sum(df$variable == df$variable[index_start])-1
    num_intersections[df$method[index_start:index_end], df$method[index_start:index_end]] = 
      num_intersections[df$method[index_start:index_end], df$method[index_start:index_end]] + 1
    index_start = index_end + 1
  }
  num_union = outer(num_rejections, num_rejections, "+") - num_intersections
  Jaccard = num_intersections/num_union
  Jaccard[is.na(Jaccard)] = 1
  return(tibble(Jaccard = mean(Jaccard[upper.tri(Jaccard)]), power = mean(power)))
}

df_Jaccard_power = df_statistics %>% 
  filter(rejection) %>% 
  group_by(A, k, rho, rep) %>% 
  arrange(variable, .by_group = TRUE) %>% 
  group_modify(get_Jaccard_and_power, keep = TRUE) %>%
  ungroup() %>%
  mutate(k = factor(k, levels = c(25, 50, 100), labels = c("k = 25", "k = 50", "k = 100")),
         rho = factor(rho, levels = c(0.2, 0.5, 0.8), labels = c("rho = 0.2", "rho = 0.5", "rho = 0.8")))

p = df_Jaccard_power %>% 
  ggplot(aes(x = power, y = Jaccard)) + geom_point() + 
  geom_hline(yintercept = 9/11, linetype = "dashed") + 
  facet_grid(k ~ rho, scales = "free_x") + theme_bw() + ylim(0,1) + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) + 
  xlab("Power") + ylab("Expected Jaccard Index")
plot(p)

p = df_Jaccard_power %>% 
  filter(k == "k = 50", rho == "rho = 0.5") %>%
  ggplot(aes(x = power, y = Jaccard)) + geom_point() + 
  geom_hline(yintercept = 9/11, linetype = "dashed") + 
  geom_smooth(span = 0.25) + 
  theme_bw() + ylim(0,1) + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "0.25", "0.5", "0.75", "1")) + 
  ggtitle("Stability of knockoff rejections") +
  xlab("Power") + ylab("Expected Jaccard Index") +
  theme(plot.title = element_text(hjust = 0.5))
plot(p)