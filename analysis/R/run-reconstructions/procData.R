processData <- function()
{
	# Steps
	# 1. Read in and process datasets
	# 2. Create R data files with arrays for fitRR.r
	# 3. Plots

	# --- 1. Read in and process datasets ----

	# -- aggregate run-size counts by day, N_gyd --
	counts <- read.csv(here("analysis/data/raw/borderCounts.csv"))
	names(counts)[names(counts)=='count_type'] <- 'gear'
	counts$gear <- as.character(counts$gear)

	# -- GSI sub-stock sampling by day, n_sgyd --
	gsi <- read.csv(here("analysis/data/raw/border-gsi-table-2024-update-full.csv")) %>%
	  rename( sample_num=fish )
	stockID <- read.csv(here("analysis/data/raw/stockIDs.csv")) %>% arrange(plotOrder)
	stockID$stockNum <- stockID$plotOrder
	gsi <- dplyr::left_join(gsi, stockID, by='CU_no')

	# change gear name to gillnet or fishWheel
	gsi$gear <- as.character(gsi$gear)
	gsi$gear[gsi$gear=='Fish Wheel'] 		<- 'fishWheel'
	gsi$gear[gsi$gear=='Test Fishery'] 		<- 'eagle'

	# remove rows with no julian_gsi or errors
	gsi <- subset(gsi, !is.na(julian))
	gsi <- subset(gsi, julian <300 & julian >100)

	# assign zeros for NA probabilities
	# gsi$prob[is.na(gsi$prob)] <- 0

	# External run size indices for fish wheel years only (<2005)
	borderPass <- read.csv(here("analysis/data/raw/border-passage.csv"))

	# Add julian day adjustment for samples/counts from fishwheel site to scale everything relative to Eagle sonar site, which is about 48 km downstream from fishwheel locations (i.e. approx 1-day travel for Chinook)
	fwDayAdj <- 1
	counts$julian[counts$count_type=='fishWheel'] <- counts$julian[counts$count_type=='fishWheel'] +1

	gsi$julian_date[gsi$data_label=='YukonRetro'] <- gsi$julian_date[gsi$data_label=='YukonRetro'] +1

	# --- 2. Create R data files with arrays for fitRR.r ---

	# Generate array with proportions by stock, day, year, gear
	stockNames <- stockID$CU_name
	stockRegion <- stockID$CU_no
	gsiGear <- c('eagle','fishWheel')
	fDay <- 160
	lDay <- 285
	days <- fDay:lDay
	fYear <- min(gsi$year)
	lYear <- max(gsi$year)
	yrs <- fYear:lYear

	# Generate array for counts by stock, day, year, and gear
	n_sdtg <- array(dim=c(length(stockNames),
						  length(fDay:lDay),
						  length(fYear:lYear),
						  length(gsiGear)))

	dimnames(n_sdtg) <- list(stockNames=stockNames,
							 julianDay=days,
							 year = yrs,
							 gears=gsiGear)

	# Fill n_sdtg array
	# loop over yrs
	for(t in 1:length(yrs))
	{
		# loop over gears
		for (g in 1:length(gsiGear))
		{
			tmp <- subset(gsi, year==yrs[t] &
								gear == gsiGear[g] &
								!is.na(julian) &
								!is.na(CU_no) &
								prob>0)

			# if no data for given gear, skip
			if(nrow(tmp)==0)
				n_sdtg[,,t,g] <-NA

			# Check if any days duplicated across samples
			smpls <- unique(tmp[,c('julian','sample_num')])

			if(any(table(smpls$sample_num)>1))
			{

			 errors <- smpls$sample_num[duplicated(smpls$sample_num)]
			  for(err in errors)
			  {
			  	errDays <- table(tmp$julian[tmp$sample_num==err]) %>%
			  			sort(decreasing=TRUE)

			  	tmp$julian[tmp$sample_num==err] <- as.integer(names(errDays)[1])
			  }
			}

			gsiDat_gt <- tmp

			#loop over days
			for(d in 1:length(days))
			{
				tmp <- subset(gsiDat_gt, year== yrs[t] &
									gear == gsiGear[g] &
									julian == days[d] &
									!is.na(julian) &
									!is.na(CU_no) &
									prob>0)


				# if no data for given day, skip
				if(nrow(tmp)==0)
				{
					n_sdtg[,d,t,g] <- NA
				}
				else
				{

					# Renormalize across samples that sum of probs !=1
					nSmpls <- length(unique(tmp$sample_num))
					if( sum(tmp$prob) != nSmpls)
						for(smp in unique(tmp$sample_num))
						{
							nProbs <- tmp[tmp$sample_num==smp,]
							if(sum(nProbs$prob) != 1)
							{
								normProbs <- nProbs$prob/sum(nProbs$prob)
								tmp$prob[tmp$sample_num==smp] <- normProbs
							}
						}

					# calculate proportions for each stocks
					n_s <- dplyr::summarize(group_by(tmp,stockNum),
									 expCounts = sum(prob))
					n_s <- dplyr::left_join(stockID,n_s, by='stockNum')
					n_s$expCounts[is.na(n_s$expCounts)] <-0


					# Check if sum of probs adds to sample nums
					if(round(sum(n_s$expCounts),10) != nSmpls )
						browser(cat('ERROR: sum of normalized GSI probs != sample size'))

					n_sdtg[,d,t,g] <- n_s$expCounts

				}

			}
		}
	}


	# Generate array for index counts
	idxGear <- c('eagle','fishWheel')
	E_dtg <- array(dim=c( length(fDay:lDay),
						  length(fYear:lYear),
						  length(idxGear)))

	dimnames(E_dtg) <- list( julianDay=days,
							 year = yrs,
							 gears=idxGear)


	# loop over yrs
	for(t in 1:length(yrs))
	{
		# loop over gears
		for (g in 1:length(idxGear))
		{
			tmp <- subset(counts, year==yrs[t] &
						  gear == idxGear[g] )



			# if no data for given gear, skip
			if(nrow(tmp)==0)
				next()

			#loop over days
			for(d in 1:length(days))
			{
				tmp <- subset(counts, year== yrs[t] &
									gear == idxGear[g] &
									julian == days[d])


				# if no data for given day, skip
				if(nrow(tmp)==0)
					E_dtg [d,t,g] <- NA
				else
				{
					E_dtg [d,t,g] <- sum(tmp$count, na.rm=T)
				}


			}
		}
	}

	# save list
	chinookYkData <- list(	E_dtg = E_dtg,
													borderPass = borderPass,
						   						n_sdtg = n_sdtg,
						   						stockNames = stockNames,
													stockRegion = stockRegion,
													gears = idxGear,
													fDay = fDay,
													lDay = lDay,
													days = days,
													fYear = fYear,
													lYear = lYear,
													yrs = yrs )

	save(chinookYkData, file=here('analysis/data/generated/chinookYkData.Rdata'))
}


