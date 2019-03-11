library(dplyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

source("../tsir_track/tsir_track_data.R")

ttd <- tsir_track_list %>%
	bind_rows(.id="citation")

ttdg <- ttd %>%
	group_by(citation) %>%
	summarize(
		year=mean(year),
		grenfell=any(author=="Grenfell")
	) %>%
	mutate(
		grenfell=ifelse(grenfell, "Includes Bryan Grenfell", "Does not include Bryan Grenfell")
	)

ttds <- ttd %>%
	group_by(citation, disease) %>%
	summarize(
		year=mean(year)
	) %>%
	mutate(
		weight=1/length(unique(disease))
	)

ttds_disease <- ttds %>%
	group_by(disease) %>%
	mutate(n=sum(weight)) %>%
	ungroup %>%
	mutate(disease=ifelse(n <= 2, "Others", disease),
		   disease=factor(disease, levels=c("Others", "Varicella", "Rubella", "HFMD", "Measles"))) 

g1 <- ggplot(ttdg) +
	geom_bar(aes(year, fill=grenfell)) +
	scale_x_continuous("Year") +
	scale_y_continuous("Number of studies") +
	theme(
		legend.title = element_blank(),
		legend.position = c(0.16, 0.86)
	)

g2 <- ggplot(ttds_disease) + 
	geom_bar(aes(year, fill=disease, weight=weight, group=interaction(year, disease))) +
	scale_x_continuous("Year") +
	scale_y_continuous("Number of studies") +
	theme(
		legend.title = element_blank(),
		legend.position = c(0.261, 0.9),
		legend.direction = "horizontal",
		legend.text = element_text(margin = margin(r = 3, unit = "pt"))
	)

gtot <- arrangeGrob(g1, g2, ncol=1)

ggsave("tsir_track_fig.pdf", gtot, width=8, height=6)
