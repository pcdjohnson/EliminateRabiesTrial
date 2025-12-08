# clear the memory
rm(list = ls())

##########################################################################################
# Project: Eliminating human rabies: impact of enhanced vaccination coverage
# Purpose: Program to randomise allocation of interventions to wards
# Written by: Paul Johnson, IBAHCM, University of Glasgow
# Version: DRAFT_V02
# Date: 2019-11-07
# Changes from DRAFT_V01: testing on new list of wards
##########################################################################################

##########################################################################################
# INSTRUCTIONS:
# 1. Generate a new random seed in the box below, via the random.org URL provided.
# 2. Replace the existing random seed of 0 at the line "rand.seed <- 0" with the new seed.
# 3. Set the working directory to the location of this program, e.g. :
#    Session -> Set Working Directory -> To Source File Location.
# 4. Run the script in R, e.g. in RStudio click "Source".
# 5. Email the output file "Data/ElimRabiesRandomisation.csv" to felix.lankester@wsu.edu.
#    Do not copy the email to the blinded trial statistician, Paul Johnson.
##########################################################################################

##########################################################################################
# Random seed:
rand.seed <- 0
set.seed(rand.seed)
# this number was generated from the URL:
# https://www.random.org/integers/?num=1&min=1&max=1000000000&col=1&base=10&format=html&rnd=new
# at:
# Timestamp: 
##########################################################################################

##########################################################################################
# Background and Methods #
# ########################
# This program reads in a data file ("Mara_wards CDP Wards.csv") listing wards
# in Mara region, Tanzania, that are to be randomized in a 1:1 ratio to 
# two methods of mass dog vaccination delivery: Pulsed and Continuous.
# Before randomisation, three exclusion criteria are applied. Wards with
#   1. Wards already being vaccinated by the existing ring vaccination program, identified 
#      by "CDP" in the column "CDP WARD".
#   2. Wards that do not have an extension officer, identified by "None" or missing value in
#      the column "Personell".
#   3. Wards that are part of the pilot study, identified by any non-NA value in
#      the column "PILOT.WARD".
#
# Wards are then randomised to pulsed or continuous vaccination by permuted block randomisation,
# stratified by district.
##########################################################################################

##########################################################################################
# start programme

# choose block size
block.size <- 6

# read in the wards data
wards.all <- read.csv("randomisation/Data/Mara_wards_FINAL.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))
nrow(wards.all)

# drop rows where there is no ward name, these contain additional data for wards named in 
# row above
wards.all <- wards.all[!is.na(wards.all$Ward), ]
nrow(wards.all)

# extract the district name
wards.all$district <- (sapply(strsplit(wards.all$LGA, " "), "[", 1))
table(wards.all$district, exclude = NULL)

# extract the ward name
wards.all$ward <- (sapply(strsplit(wards.all$Ward, " "), function(x) paste(x[-1], collapse = "")))

table(wards.all$ward, exclude = NULL)
table(table(wards.all$ward, exclude = NULL)) # ward names are not unique - need to add district to make unique ID

# clean up district and ward names (restrict to lower case letters)
for(n in c("ward", "district")) {
  wards.all[, n] <- gsub(" ", "", gsub("[^[:alpha:] ]", "", tolower(wards.all[, n])))
  wards.all[!grepl('^[A-Za-z]+$', wards.all[, n]), n] <- NA
} ; rm(n)

# create unique ID by joining district and ward name
wards.all$id <- paste(wards.all$district, wards.all$ward, sep = "-")

# are there any duplicated ward names?
length(wards.all$id)
length(unique(wards.all$id))
wards.all$id[duplicated(wards.all$id)]
# no

# apply the exclusion criteria 
#   1. exclude "CDP" in the column "CDP WARD".
#   2. exclude "None" and NA in the column "Personell".
wards <- 
  wards.all[!wards.all$CDP.WARD %in% "CDP" & 
              !wards.all$Personell %in% c("None", NA) & 
              wards.all$PILOT.WARD %in% NA, ]
nrow(wards.all)
nrow(wards)
rm(wards.all)

table(wards$Personell, exclude = NULL)

# sort wards data by district then ward
wards <- wards[order(wards$id), ]

# drop columns that are not required
wards$Region <- wards$Name <- wards$Personell <- wards$CDP.WARD <- 
  wards$LGA <- wards$Ward <- wards$PILOT.WARD <- NULL
head(wards)

# in order to stratify by district, we need the number of wards in each district
strat.size <- table(wards$district)
strat.size

# add randomisation order
for(s in names(strat.size)) {
  wards$rand.order[wards$district == s] <- sample(1:strat.size[s])
}; rm(s)

# re-order wards by district then randomisation order
wards <- wards[order(wards$district, wards$rand.order), ]

# randomisation

# the aim is to randomise the n = 
nrow(wards)
sum(strat.size)
# ...wards in equal numbers to two groups, stratified by district (= stratum), 
# using permuted block randomization. this method ensures that equal 
# numbers are allocated to each intervention group within each stratum, 
# with the wards randomly allocated to groups within each stratum, within 
# blocks of a specified size

allocation.list <-
  lapply(names(strat.size), function(b) {
    # n wards in stratum
    n <- strat.size[b]
    # n short of perfect blocks
    n.short <- block.size - (n %% block.size)
    # make list of allocations, for n + n.short wards
    # giving complete blocks
    n.extra <- n + n.short
    # divide wards into 2 groups labelled 1 and 2
    allocation <- 1:n.extra %% 2 + 1
    # permute the allocations randomly, within blocks of size = block.size 
    # and trim the extra wards off the end
    allocation.matrix <- matrix(allocation, nrow = block.size)
    rm(allocation)
    allocation <- c(apply(allocation.matrix, 2, sample))[1:n]
    # randomly select which intervention to assign to 1 and which to 2
    groups <- sample(c("Pulsed", "Continuous"))
    # create randomised sequence
    out <- groups[allocation]
    # check numbers are as close to equal as possible
    print(diff(table(groups[allocation])))
    # output randomised sequence
    out
  })
names(allocation.list) <- names(strat.size)

# add the intervention allocations to the wards data set
for(d in names(strat.size)) {
  wards$intervention[wards$district == d] <- allocation.list[[d]]
}

# check that the numbers in each group are reasonably even within each district
tapply(wards$intervention, wards$district, function(x) diff(table(x)))

# check the number of wards allocated to each treatment overall
table(wards$intervention)

# add random seed and date-time stamp to the output 
wards$rand.seed <- rand.seed
wards$date.time.stamp <- Sys.time()

# program name
wards$program.name <- 
  list.files()[substr(list.files(), nchar(list.files()) - 1, nchar(list.files())) == ".R"]

# add user name
wards$program.user <- Sys.info()["user"]

# R version 
wards$r.version <- version$version.string

# loaded packages
wards$packages.loaded <-
  if(is.null(sessionInfo()$otherPkgs)) "none" else
    paste0(names(sessionInfo()$otherPkgs), unlist(sapply(sessionInfo()$otherPkgs, "[", "Version")),
           collapse = ";")

head(wards)
  
# export data table with randomised intervention allocations as csv
write.csv(wards, file = "randomisation/Data/ElimRabiesRandomisation.csv", row.names = FALSE)

# end programme
##########################################################################################
