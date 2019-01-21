library(dplyr)

beta_t <- 30 * c(1.24, 1.14, 1.16, 1.31, 1.24, 1.12, 1.06, 1.02, 0.94, 0.98, 1.06, 1.08, 0.96,
				 0.92, 0.92, 0.86, 0.76, 0.63, 0.62, 0.83, 1.13, 1.20, 1.11, 1.02, 1.04, 1.08)

alpha <- 0.97

measles_data <- read.csv("../data/measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_list <- measles_UK %>%
	split(as.character(.$loc))

ld <- measles_list$LONDON

beta <- beta_t[ld$biweek]

dd <- data.frame(
	linear=beta * (ld$cases/0.52)^alpha/ld$pop,
	exponential=1-exp(-beta * (ld$cases/0.52)^alpha/ld$pop)
)


