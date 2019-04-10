library(dplyr)
library(ggplot2); theme_set(theme_bw())

measles_data <- read.csv("../data/measlesUKUS.csv")

bs <- measles_data %>%
	filter(loc=="BOSTON", year >= 1930, year < 1940)

ld <- measles_data %>%
	filter(loc=="LONDON", year > 1955)

alldata <- rbind(bs, ld) %>%
	mutate(loc=factor(loc, levels=c("LONDON", "BOSTON"), labels=c("London", "Boston")))

g1 <- ggplot(alldata) +
	geom_point(aes(decimalYear, cases), size=1, shape=1) +
	scale_y_continuous("Cases") +
	scale_x_continuous("Date", breaks=seq(1930, 1970, by=1), expand=c(0,0)) +
	facet_wrap(~loc, scale="free", ncol=1) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)

ggsave("fig1.pdf", g1, width=6, height=4)
