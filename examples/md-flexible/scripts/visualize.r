#!/usr/bin/env Rscript

## Dependencies
library(tidyverse)


## Parse arguments
argv <- commandArgs(TRUE)
argc <- length(argv)

datalist = vector("list", length=argc)
for (i in 1:argc) {
	file <- argv[i]
	freq <- as.numeric(gsub("\\D", "", file))

	datalist[[i]] <- read_csv(file) |>
		filter(Iteration %% freq == 0) |> # Only get full respa iterations
		mutate(type = basename(file))
}


data <- dplyr::bind_rows(datalist)
plot <- ggplot(data, aes(x=Iteration, y=`Potential Energy`, colour=factor(type))) + geom_line(linewidth=0.5)

ggsave("potential_energy.png", width=20, height=10)
