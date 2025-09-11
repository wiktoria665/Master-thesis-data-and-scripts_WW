require(flowCore)
require(ggcyto)
require(tidyverse)
require(knitr)
require(gridExtra)
require(openCyto)
require(flowClust)
require(readr)
require(flowStats)
require(data.table)
require(openxlsx)


### analyse flow cytometry results to obtain green cell fraction
plates <- c("01", "02", "05", "06", "07", "08")
path <- "~/cyto/experiment2/"
excel_file <- "exp2-results.xlsx"
merged_counts <- data.frame(matrix(nrow = 0, ncol = 8))
merged_fractions <- data.frame(matrix(nrow = 0, ncol = 5))

for (i in plates) {
	fs <-
		read.flowSet(
			path = path,
			pattern = paste(i, "-Well.*.fcs", sep = ""),
			alter.names = TRUE,
			truncate_max_range = FALSE
		)
	pData(fs)$well <-
		paste(
			gsub("(.*-).*-.*.fcs", "\\1", sampleNames(fs)),
			gsub(".*-.*-(.*).fcs", "\\1", sampleNames(fs)),
			sep = ""
		) # extract well from name and add new 'well' column
	colnames(fs)
	
	
	##transformation
	truncate <-
		truncateTransform(transformationId = "defaultTruncateTransform", a = 1200) # to get rid of negative and low fluorescence
	l10 <-
		logTransform(
			transformationId = "Log10transform",
			logbase = 10,
			r = 1,
			d = 1
		) # to separate yellow/green
	my_trans <-
		transformList(c("FL1.A", "FL1.H", "FL2.A", "FL2.H"),
									list(truncate, truncate, truncate, truncate))
	fs.t <- transform(fs, my_trans)
	my_trans <-
		transformList(c("FL1.A", "FL1.H", "FL2.A", "FL2.H"), list(l10, l10, l10, l10))
	fs.t <- transform(fs.t, my_trans)
	
	## create "big" data in one experiment - use it for gating only
	# data.1frame <- fs.t[[1]]
	# exprs(data.1frame) <- fsApply(fs.t, function(x) {
	# 	x <- exprs(x)
	# 	return(x)
	# })
	
	## final gating set
	gs <- GatingSet(fs.t)
	
	## automatic gating for singletons (from previous analyses)
	singlet.gate <- readRDS("singlet-gate.rds")
	singlet.gate@filterId <- "singlet"
	
	## check gate on merged data
	# ggcyto(data.1frame, aes(x = FSC.A, y = FSC.H), subset = "root") +
	# 	geom_hex(bins = 500) +
	# 	geom_gate(singlet.gate) +
	# 	ggcyto_par_set(limits = list(x = c(0, 1.5e6), y = c(0, 1.5e6))) +
	# 	geom_stats(adjust = 0.8, location = "plot") # check gate
	
	gs_pop_add(gs, singlet.gate, name = "singlet")
	recompute(gs) # recompute GatingSet
	
	## plot first data with first gate for separately for all wells
	a = 1
	for (j in seq(8)) {
		b = min(c(a + 11, 96))
		ggcyto(gs[[a:b]], aes(x = FSC.A, y = FSC.H), subset = "root") +
			facet_wrap(~ well, ncol = 4) +
			ggcyto_par_set(limits = list(x = c(0, 1.5e6), y = c(0, 1.5e6))) +
			geom_hex(bins = 200) +
			geom_gate("singlet") +
			geom_stats(adjust = 0.8, location = "plot") +
			theme_minimal()
		
		a <- a + 11
		
		ggsave(paste(
			path,
			"gate-singlet-plate",
			i,
			"-",
			as.character(j),
			".pdf",
			sep = ""
		))
	}
	
	## add green and yellow gates (manually set)
	g.green <-
		polygonGate(
			filterId = "green",
			"FL1.A" = c(3, 6, 6, 3),
			"FL2.A" = c(3.05, 6.05, 6.35, 3.35)
		) # define gate
	g.yellow <-
		polygonGate(
			filterId = "yellow",
			"FL1.A" = c(3, 6, 6, 3),
			"FL2.A" = c(3.35, 6.35, 6.65, 3.65)
		) # define gate
	
	## check gate on merged data
	# ggcyto(data.1frame, aes(x = FL1.A, y = FL2.A), subset = "singlet") +
	# 	geom_hex(bins = 200) +
	# 	geom_gate(g.yellow, colour = "yellow") +
	# 	geom_gate(g.green, colour = "green") +
	# 	ggcyto_par_set(limits = "data") +
	# 	geom_stats(adjust = 0.8, location = "data") # check gate
	
	
	gs_pop_add(gs, g.yellow, parent = "singlet") # add gate to GatingSet
	gs_pop_add(gs, g.green, parent = "singlet") # add gate to GatingSet
	recompute(gs) # recompute GatingSet
	
	a = 1
	for (j in seq(8)) {
		b = min(c(a + 11, 96))
		ggcyto(gs[[a:b]], aes(x = FL1.A, y = FL2.A), subset = "singlet") +
			facet_wrap(~ well, ncol = 4) +
			geom_hex(bins = 200) +
			geom_gate(g.yellow, colour = "yellow") +
			geom_gate(g.green, colour = "green") +
			ggcyto_par_set(limits = list(x = c(3, 6.7), y = c(3, 6.7))) +
			theme_minimal()
		
		a <- a + 8
		ggsave(
			paste(
			  path,
				"gate-green-yellow-plate",
				i,
				"-",
				as.character(j),
				".pdf",
				sep = ""
			)
		)
	}
	
	
	ps <- gs_pop_get_count_with_meta(gs)
	ps <- ps %>%
		mutate("percent_of_parent" = Count / ParentCount)
	
	
	ps.wide <- ps %>%
		select(name, well, Population, Count, ParentCount) %>%
		pivot_wider(names_from = c(Population),
								values_from = c(Count, ParentCount)) %>%
		as.data.frame()
	
	names(ps.wide) <-
		str_replace(str_replace_all(names(ps.wide), "/", "."), "_.", "_")
	ps.wide <- ps.wide %>%
		mutate("green_fraction" = Count_singlet.green / (Count_singlet.yellow + Count_singlet.green)) %>%
		select(name,
					 well,
					 Count_singlet.green,
					 Count_singlet.yellow,
					 green_fraction) %>%
		mutate("plate" = gsub("-[A-H][0-9]*", "", well)) %>%
		mutate(well = gsub("[0-9][0-9]-", "", well)) %>%
		rename("count_singlet.green" = Count_singlet.green,
					 "count_singlet.yellow" = Count_singlet.yellow) %>%
		select(plate,
					 well,
					 count_singlet.green,
					 count_singlet.yellow,
					 green_fraction)
	
	merged_counts <- rbind(merged_counts, ps)
	merged_fractions <- rbind(merged_fractions, ps.wide)
}
colnames(merged_counts)	<- colnames(ps)
colnames(merged_fractions)	<- colnames(ps.wide)

## write results to excel
excel <-
	buildWorkbook(
		merged_counts,
		sheetName = "counts",
		colNames = TRUE,
		rowNames = FALSE
	)

addWorksheet(
	excel,
	sheetName = "fractions"
)

writeData(
	excel, 
	"fractions", 
	merged_fractions, 
	rowNames = FALSE, 
	colNames = TRUE
)


saveWorkbook(excel, paste(path, excel_file, sep=""), overwrite = TRUE)



