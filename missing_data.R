# Assessing and handeling missing data

# Describing missing data
library(VIM)

png(file.path(results_dir, "missing.png"), width = 9, height = 6.6, units = "in", res = 150)
missing <- aggr(panel, numbers = TRUE,
     sortVars = TRUE, labels = names(panel),
     cex.axis = .5, gap = 3)
dev.off()

# Describing response time (should have been a exclusion criteria(?))
## minutes
summary(panel$answer_time_ms)/60000
