## Program start ----
#
#    Eliminate Rabies Trial analysis
#    Started: 24 Jan 2024
#


## Settings and global options ----

# Remove objects, close graphics devices
rm(list=ls())
graphics.off()

# Set timer running
start.time <- Sys.time()

# Global settings

# Is this the final run, with larger and slower bootstrap and MCMC samples,
# and with RNG seeds set?
final.run <- TRUE

# Define vaccination coverage
# Use vaccinated_M3 as the primary outcome - more reliable, although also lots of missing data
# but include analysis of an alternative definition of coverage, vaccinated_M4,
vaccination.definition <- c(main = "Vaccinated_M3", alternative = "Vaccinated_M4")[1]

# Make simulations reproducible by setting RNG seed?
set.rng.seed <- final.run

# Iteration multiplier (higher values give slower and more accurate results)
iter.mult <- ifelse(final.run, 10, 0.1)

# Set minimum acceptable coverage threshold
cov.thresh <- 0.4

# Need more time to load large data files
options(timeout = 1000) 

# mclapply needs this for set.seed to work
RNGkind("L'Ecuyer-CMRG") 

# Turn off all open graphics devices
graphics.off()

# Choose colours to represent the two trial arms on plots
# (taken from the coat of arms of Tanzania 
# (https://en.wikipedia.org/wiki/Coat_of_arms_of_Tanzania#/media/File:Coat_of_arms_of_Tanzania.svg) 
# and checked for colour-blind visibility here: https://boydorr.gla.ac.uk/lucanelli/colorquest/ )
arm.colours <- c(`Team-based` = "#233E9B", `Community-based` = "#DA1728")


## Set report properties, working directory, filenames, word macro, etc ----

# Text for title page of tables output
study.title <- "Eliminate Rabies Trial"             
subtitle <- 
  paste0("Eliminating human rabies: impact of enhanced vaccination coverage (",
         names(vaccination.definition), " coverage definition)")
report.purpose<-"Analysis"
brief.purpose<-"analysis"
draft.or.final<-ifelse(final.run, "Final", "Draft")
report.version<-"1.0"
report.status<-paste(draft.or.final,"Report Version",report.version)    
author<-"Paul Johnson"

# Set program directory and filename
root.dir<-""
prog.directory<-root.dir
prog.name<-
  paste(
    tolower(gsub(" ","_",study.title)),"_",
    brief.purpose,"_report_v",
    report.version,
    sep="")
prog.file<-paste0(prog.directory,prog.name,".R")


# Set output file for tables
file.out.directory <- 
  ifelse(vaccination.definition == "Vaccinated_M3", paste0(prog.directory,brief.purpose,"_report"), 
         paste0(prog.directory,brief.purpose,"_report_alt_outcome_definition"))
file.out.name <-
  paste0(prog.name, 
         ifelse(vaccination.definition == "Vaccinated_M3", "", "_alt_outcome_definition"),
         ".doc")
file.out<-paste(file.out.directory,file.out.name,sep="/")

# Create directories for figures and header/footer. 
# Delete old folder first to prevent old figures being used in latest version of report.
figures.directory<-
  paste0(
    file.out.directory,"/",substr(file.out.name,1,nchar(file.out.name)-4),
    "_figures")
header.directory<-
  paste0(file.out.directory,"/",substr(file.out.name,1,nchar(file.out.name)-4),"_files")
if(file.exists(figures.directory)) unlink(figures.directory,recursive=TRUE)
if(file.exists(header.directory))  unlink(header.directory,recursive=TRUE)
lapply(list(figures.directory,header.directory),dir.create, recursive = TRUE)


## Additional libraries ----

# Load functions 
source(paste0(prog.directory,"functions/elim_rabies_functions.R"))

# Load packages
need.packages <-
  c("lubridate", "foreign", "Hmisc", "gdata", "lme4", "glmmTMB", "DHARMa", "ggplot2", "sjPlot",
    "parallel", "forcats", "brms", "bayesplot", "coda", "see", "vioplot", "sf", "spdep")
sapply(need.packages,install.load)


## Import and organise data ----

### Import, process and inspect data sets ----

# Load survey data
dat <- read.csv("data/T3_Trial_Data.csv")
(data.names <- names(dat))
nrow(dat)

# Sort survey data
dat <- dat[order(dat$Round, dat$ward, dat$sub_village, dat$Household.file.id), ]

# Load vaccination data
vax <- read.csv("data/T3_Vac_Trial_Data.csv")
dim(vax)
names(vax)
head(vax)
table(dat$pet.id %in% vax$pet.id)
table(vax$pet.id %in% dat$pet.id)
table(vax$Household.file.id %in% dat$Household.file.id)
table(vax$registration.number %in% dat$registration.number[dat$Vaccinated_M3 %in% 1])
table(dat$registration.number[dat$Vaccinated_M3 %in% 1] %in% vax$registration.number)

# Make trial arm a factor
vax$Trial_Arm <- factor(vax$Trial_Arm, c("Pulse", "Continuous"), names(arm.colours))

# Create ward ID:
vax$ward <- factor(paste(vax$ward, substr(vax$district, 1, 1), sep = "-"))

# Format vax dates
vax$date.vaccinated <- as.Date(vax$date.vaccinated, format = "%d/%m/%Y")
vax$day.vaccinated <- as.numeric(vax$date.vaccinated - min(vax$date.vaccinated)) + 1


# Load estimated number of dogs per ward
dog <- read.csv("data/T3_Ward_Pop_Data.csv")
dim(dog)
names(dog)
head(dog)
# Create ward factor matching other data sets
dog$ward <- factor(paste(dog$ward_name, substr(dog$dist_name, 1, 1), sep = "-"))
# Check all trial wards are present
setdiff(levels(dog$ward), levels(vax$ward))
setdiff(levels(vax$ward), levels(dog$ward))
# Drop non-trial wards
dog <- droplevels(dog[as.character(dog$ward) %in% levels(vax$ward), ])
dim(dog)

# Check for missing pet IDs
miss.petid <- is.na(dat$pet.id)
sum(miss.petid)
dat[miss.petid, ]
nrow(dat)
dat <- dat[!miss.petid, ]
nrow(dat)

# How much missing data? (% correct on 25 Aug 2023)
cbind(data.names, paste0(round(100 * colSums(is.na(dat[, data.names]))/nrow(dat), 1), "%"))
# Trial arm will be 100% missing, intentionally due to blinding,
# until the trial closes.
# Reg number manual is 99% missing, because it is used only if the 
#   registration number was missing/wrong.
# Date vaccinated is 67% missing, only available if the dog has a fully
#   complete and correct vaccination card at the time of the survey.
# Vaccinated_M1Y, Vaccinated_M1, and Vaccinated_M3 are all approx 1/3 
#   missing. Vaccinated_M2 has no missing data. These are 4 methods of 
#   classifying dogs as vaccinated or unvaccinated. See data dictionary at
#   end of this script for definitions.
# All other columns have no missing data

# Create factors

# Use only the first word of district, in order to use the original six
# districts used in the design and randomisation

dat$district6 <- sapply(strsplit(dat$district, " "), "[", 1)
dat$district <- factor(dat$district)

# Age of dogs
dat$age.years <- dat$age.years + dat$age.months/12
label(dat$age.years) <- "Age (years)"

# Round - change language from round (which risks confusion with vaccination round)
# to annual survey visit (V1, V2)
# Change round from 1-6 to Year1-V1, Year1-V2, ..., Year3-V2
table(dat$Round, exclude = NULL)
rounds <- 
  c(R1 = "Y1-V1", R2 = "Y1-V2", 
    R3 = "Y2-V1", R4 = "Y2-V2", 
    R5 = "Y3-V1", R6 = "Y3-V2")
dat$YearVisit <- factor(rounds[dat$Round])
table(dat$YearVisit, dat$Round, exclude = NULL)
label(dat$YearVisit) <- "Year-visit"

# Create year factor
dat$Year <- factor(substr(dat$YearVisit, 1, 2))
label(dat$Year) <- "Year"

# Create new visit factor
dat$Visit <- factor(substr(dat$YearVisit, 4, 5))
table(dat$Year, dat$Visit)
label(dat$Visit) <- "Visit"

# Month of vaccination
dat$month.vaccinated <- 
  factor(month(dat$date.vaccinated, label = TRUE), 
         c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),
         ordered = FALSE)
label(dat$month.vaccinated) <- "Month vaccinated"

# Precise dates of assessments
dat$date.added <- as.Date(dat$date.added, format = "%d/%m/%Y")

# Create a variable for the day of the trial year on which each dog was surveyed
for(y in levels(dat$Year)) {
  y.num <- as.numeric(gsub("Y", "", y))
  dat$YearDay[dat$Year == y] <- (dat$date.added[dat$Year == y] - as.Date("2020-11-01") - 365 * (y.num - 1))
}

# Box plot of day of year vaccinated by trial year (1, 2 or 3) and visit (1 or 2).
# The dashed line at end of year shows the extent to which vaccinations leaked into 
# the following year.
boxplot(YearDay ~ YearVisit, data = dat, ylim = c(0, max(dat$YearDay)))
abline(h = 365, lty = 2)


# Calculate the date range for each of the 6 surveys (2 visits X 3 years)
visit.durations <- tapply(dat$date.added, dat$Visit, function(x) diff(range(x)))
round(12 *  visit.durations/ 365, 1)
round(12 *  mean(visit.durations)/ 365, 1)
mean.visit.dates <- as.Date(round(tapply(dat$date.added, dat$YearVisit, mean)))
visit.date.ranges <- 
  cbind(FirstSurvey = as.character(as.Date(round(tapply(dat$date.added, dat$YearVisit, min)))),
        LastSurvey = as.character(as.Date(round(tapply(dat$date.added, dat$YearVisit, max)))))
rownames(visit.date.ranges) <- levels(dat$YearVisit)
visit.date.ranges

# Duration of each survey in months
round(12 * diff(mean.visit.dates)/365, 1)
round(12 * mean(diff(mean.visit.dates)[c(1, 3, 5)])/365, 1)

# Make duplicated ward names unique by adding the district name
length(unique(dat$ward))
dat$ward <- factor(paste(dat$ward, substr(dat$district, 1, 1), sep = "-"))
table(dat$ward)
nlevels(dat$ward)

# Check that ward names match between survey data (dat) and vaccinations data (vax)
setdiff(levels(dat$ward), levels(vax$ward))
setdiff(levels(vax$ward), levels(dat$ward))

# Create unique sub-village names:
length(unique(dat$sub_village))
dat$sub_village <- factor(paste(dat$sub_village, dat$ward, sep = "-"))
nlevels(dat$sub_village)

# Create a factor for sub-village-visit, for fitting a random effect to.
dat$sub_village.rnd <- factor(paste(dat$sub_village, dat$YearVisit, sep = "-"))

# Create a factor for ward-visit
dat$ward.rnd <- 
  factor(paste(dat$ward, dat$YearVisit, sep = "-"),
         apply(expand.grid(levels(dat$ward), levels(dat$YearVisit)), 1, paste, collapse = "-"))

# Create a factor for district-visit, for fitting a random effect to.
dat$district.rnd <- factor(paste(dat$district, dat$YearVisit, sep = "-"))

# Land use factor
dat$LandUse <- factor(dat$LandUse)
label(dat$LandUse) <- "Land use"

# Sex factor
dat$Sex <- factor(dat$Sex)
label(dat$Sex) <- "Sex"

# Sex and pregnancy factor
dat$Sex.pregnant <- 
  factor(paste0(dat$Sex, dat$pregnant), 
         c("FemaleFALSE", "FemaleTRUE", "MaleFALSE"),
         c("Non-pregnant female", "Pregnant female", "Male"))
label(dat$Sex.pregnant) <- "Sex and pregnancy"

# Inspect data
names(dat)
head(dat)
sample.rows(dat)
dim(dat)

# Get list of ward-visit with fewer than 30 dogs, and those with more than 60
ward.rnd.tab <- table(dat$ward.rnd[dat$YearVisit != "Y3-V2"], exclude = NULL)
ward.rnd.tab <- ward.rnd.tab[-grep("-Y3-V2", names(ward.rnd.tab))]
ward.rnd.eq0 <- sort(names(ward.rnd.tab)[ward.rnd.tab == 0])
ward.rnd.lt30 <- sort(names(ward.rnd.tab)[ward.rnd.tab < 30 & ward.rnd.tab > 0])
ward.rnd.gt60 <- sort(names(ward.rnd.tab)[ward.rnd.tab > 60])
#unique(dat[dat$ward %in% substr(ward.rnd.eq0, 1, nchar(ward.rnd.eq0) - 6), c("Village", "ward", "latitude")])
sort(table(dat$ward[dat$ward.rnd %in% ward.rnd.lt30]))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.lt30 & dat$LandUse == "Rural", c("ward", "YearVisit")])))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.lt30 & dat$LandUse == "Semi-Urban", c("ward", "YearVisit")])))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.lt30 & dat$LandUse == "Urban", c("ward", "YearVisit")])))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.gt60 & dat$LandUse == "Rural", c("ward", "YearVisit")])))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.gt60 & dat$LandUse == "Semi-Urban", c("ward", "YearVisit")])))
table(unique(droplevels(dat[dat$ward.rnd %in% ward.rnd.gt60 & dat$LandUse == "Urban", c("ward", "YearVisit")])))

# Ward-visits that didn't survey any dogs:
ward.rnd.eq0

# Calculate the proportion of entries that are unique for each column
round(sapply(dat, function(x) length(unique(x)))/nrow(dat), 5)

# How many dogs have been recorded more than once?
# Pet ID is not a true pet ID, because dogs are assigned a new pet ID at every
# survey visit, even if they've been surveyed before
table(table(dat$pet.id))

# How many households have been sampled?
length(unique(dat$Household.file.id))

# Distribution on no of dogs across HH?
table(table(dat$Household.file.id))

# Distribution on no of dogs by visit?
barplot(table(dat$YearVisit), ylab = "No of dogs", xlab = "Visit")

# Number of dogs by ward
barplot(rev(sort(table(dat$ward))), ylab = "No of dogs", xlab = "Ward", names.arg = FALSE)

# Create a data set of wards
wards <- data.frame(name = levels(dat$ward))
rownames(wards) <- wards$name

# Add district factors to the wards data
wards$district <- factor(tapply(as.character(dat$district), dat$ward, unique)[rownames(wards)])
wards$district6 <- factor(tapply(dat$district6, dat$ward, unique)[rownames(wards)])
label(wards$district) <- label(wards$district6) <- "District"

# Add land use to the wards data frame
wards$LandUse <- 
  factor(tapply(as.character(dat$LandUse), dat$ward, unique)[rownames(wards)],
         levels(dat$LandUse))
label(wards$LandUse) <- "Land use"

# Re-order wards by district order
wards <- wards[order(wards$district), ]

# Add the total number of dogs vaccinated in each ward
wards$n.vaccinated <- table(vax$ward)[rownames(wards)]
label(wards$n.vaccinated) <- "N vaccinations"

# Add the total number of dogs in each ward
wards$n.dogs <- dog$DogPop[match(rownames(wards), dog$ward)]
label(wards$n.dogs) <- "N dogs (estimated)"

# Add the intervention allocations to the wards data set
wards$Trial_Arm <- 
  factor(tapply(dat$Trial_Arm, dat$ward, unique)[rownames(wards)], 
         c("Pulse", "Continuous" ), names(arm.colours))

# Check that the numbers in each group are reasonably even within each district
tapply(wards$Trial_Arm, wards$district6, function(x) diff(table(x)))

# Check the number of wards allocated to each treatment overall
table(wards$Trial_Arm)
table(wards$Trial_Arm, wards$district6)

# Total number of vaccinations delivered in each ward
par(mfrow = c(2, 1))
hist(wards$n.vaccinated)
boxplot(n.vaccinated ~ Trial_Arm, data = wards)
title("Total number of vaccinations delivered\nin each ward, by trial arm")
par(mfrow = c(1, 1))

# Make Trial_Arm a factor in the survey data frame
dat$Trial_Arm <- factor(dat$Trial_Arm, c("Pulse", "Continuous"), names(arm.colours))
table(dat$Trial_Arm, exclude = NULL)

# Create a visit-arm factor for later predictions
dat$Trial_Arm.rnd <- 
  factor(paste(dat$YearVisit, dat$Trial_Arm, sep = "-"),
         paste0(
           c("Y1-V1-", "Y1-V1-", 
             "Y1-V2-", "Y1-V2-", 
             "Y2-V1-", "Y2-V1-", 
             "Y2-V2-", "Y2-V2-", 
             "Y3-V1-", "Y3-V1-", 
             "Y3-V2-", "Y3-V2-"), 
           names(arm.colours)))


### Choose the outcome measure (define coverage) ----

# Define coverage following the choice of outcome measure in the global settings
dat$Vaccinated <- dat[, vaccination.definition]

# Create categorical variables for vaccinated / not / missing
dat$Vaccinated.cat <- dat$Vaccinated
dat$Vaccinated.cat[is.na(dat$Vaccinated.cat)] <- -1
dat$Vaccinated.cat <- factor(dat$Vaccinated.cat, c(0, 1, -1), c("Not vaccinated", "Vaccinated", "Missing"))
dat$Vaccinated.cat.na <- factor(dat$Vaccinated.cat, c("Not vaccinated", "Vaccinated"))
table(dat$Vaccinated.cat.na, dat$Vaccinated, exclude = NULL)
table(dat$Vaccinated.cat, dat$Vaccinated, exclude = NULL)
label(dat$Vaccinated.cat.na) <- label(dat$Vaccinated.cat) <- 
  label(dat$Vaccinated) <- "Primary coverage outcome"

# Secondary outcome measure: vaccinated_M4
# Vaccinated dogs are counted as in the primary outcome, M3, 
# but dogs who are claimed to be vaccinated but with no certificate are counted as not
# vaccinated and therefore added to the denominator
dat$Vaccinated2 <- dat$Vaccinated_M4
dat$Vaccinated2.cat <- dat$Vaccinated2
dat$Vaccinated2.cat <- factor(dat$Vaccinated2.cat, c(0, 1), c("Not vaccinated", "Vaccinated"))
label(dat$Vaccinated2.cat) <- label(dat$Vaccinated2) <- "Secondary coverage outcome"
table(Primary = dat$Vaccinated.cat, Secondary = dat$Vaccinated2.cat, exclude = NULL)


# Add coverage data to the wards data set
wards$Coverage <- tapply(dat$Vaccinated, dat$ward, mean, na.rm = TRUE)[rownames(wards)]
label(wards$Coverage) <- "Mean coverage (primary outcome)"
wards$Coverage2 <- tapply(dat$Vaccinated2, dat$ward, mean)[rownames(wards)]
label(wards$Coverage2) <- "Mean coverage (secondary outcome)"

# Add n dogs vaxxed and unvaxxed under both primary and secondary outcomes.
# These are actually the same but calculate them separately for completeness.
wards$Vaccinated <- tapply(dat$Vaccinated, dat$ward, sum, na.rm = TRUE)[rownames(wards)]
label(wards$Vaccinated) <- "N dogs vaccinated"
wards$Vaccinated2 <- tapply(dat$Vaccinated2, dat$ward, sum)[rownames(wards)]
label(wards$Vaccinated2) <- "N dogs vaccinated (secondary outcome)"
wards$Total <- tapply(dat$Vaccinated, dat$ward, function(x) sum(!is.na(x)))[rownames(wards)]
label(wards$Total) <- "N dogs surveyed"
wards$Total2 <- tapply(dat$Vaccinated2, dat$ward, function(x) sum(!is.na(x)))[rownames(wards)]
label(wards$Total2) <- "N dogs surveyed (secondary outcome)"

# Show the number of dogs vaccinated by subvillage and visit-year for a random ward
tab <- 
  with(dat[dat$ward == sample(levels(dat$ward), 1), ], {
    table(sub_village, YearVisit)
  })
tab[rowSums(tab) > 0, ]

# Compare coverage between the two outcome definitions
plot(wards$Coverage, wards$Coverage2, cex = 4 * wards$Vaccinated/max(wards$Vaccinated))
abline(0, 1)
cor.test(wards$Coverage, wards$Coverage2)

# Compare survey-based vaccination coverage 
# with coverage based on the number of vaccines delivered
plot(wards$Coverage, wards$n.vaccinated/(3*wards$n.dogs), xlim = 0:1, ylim = 0:1)
axis(2)
abline(0, 1)
cor.test(wards$Coverage, wards$n.vaccinated/(3*wards$n.dogs))

# Add other characteristics to the wards data frame
head(dat)
wards$Proportion_female <- tapply(dat$Sex == "Female", dat$ward, mean, na.rm = TRUE)[rownames(wards)]
wards$age.years <- tapply(dat$age.years, dat$ward, mean, na.rm = TRUE)[rownames(wards)]
label(wards$age.years) <- "Mean age (years)"
wards$SV_V_Distance <- tapply(dat$SV_V_Distance, dat$ward, mean, na.rm = TRUE)[rownames(wards)]
label(wards$SV_V_Distance) <- "Mean subvillage-village distance (km)"
wards$Completeness <- 100 * rowSums(table(dat$ward, dat$YearVisit) > 0)[rownames(wards)] / nlevels(dat$YearVisit)
label(wards$Completeness) <- "Survey completeness (%)"
wards$lat <- tapply(dat$SV_latitude, dat$ward, mean, na.rm = TRUE)[rownames(wards)]
wards$lon <- tapply(dat$SV_longitude, dat$ward, mean, na.rm = TRUE)[rownames(wards)]


# Compare raw coverage between arms, by land use category
tapply(wards$Coverage, list(wards$LandUse, wards$Trial_Arm), mean)
table(wards$LandUse, wards$Trial_Arm)

# Calculate the number of households surveyed by ward and visit-year
no.of.hh.surveyed <- 
  tapply(dat, dat[, c("ward", "YearVisit")], 
         function(x) sum(colSums(table(droplevels(x)$Household.file.id, droplevels(x)$sub_village) > 0.5)),
         simplify = FALSE)
summary(unlist(no.of.hh.surveyed))
sort(unlist(no.of.hh.surveyed))
table(unlist(no.of.hh.surveyed) >= 25)
mean(unlist(no.of.hh.surveyed) >= 25)


# Create a data set of wards x visits
wards.rnd <- 
  data.frame(ward.rnd = apply(expand.grid(levels(dat$ward), levels(dat$YearVisit)), 1, paste, collapse = "-"))
nrow(wards.rnd)/6
rownames(wards.rnd) <- wards.rnd$ward.rnd
wards.rnd$YearVisit <- 
  factor(substr(wards.rnd$ward.rnd, nchar(wards.rnd$ward.rnd) - 4, nchar(wards.rnd$ward.rnd)))
wards.rnd$Year <- factor(substr(wards.rnd$YearVisit, 1, 2))
wards.rnd$Visit <- factor(substr(wards.rnd$YearVisit, 4, 5))
wards.rnd$ward <- 
  factor(substr(wards.rnd$ward.rnd, 1, nchar(wards.rnd$ward.rnd) - 6))
wards.rnd$ward.rnd <- factor(wards.rnd$ward.rnd)
wards.rnd$ward.Year <- factor(paste(wards.rnd$ward, wards.rnd$Year, sep = "-"))
wards.rnd$district <- factor(wards[as.character(wards.rnd$ward), "district"])
wards.rnd$district.rnd <- factor(paste(wards.rnd$district, wards.rnd$YearVisit, sep = "-"))
head(wards.rnd)
dim(wards.rnd)

# Add trial arm
wards.rnd$Trial_Arm <- factor(wards[as.character(wards.rnd$ward), "Trial_Arm"], levels(dat$Trial_Arm))

# Make a factor for trial arm x year x visit
wards.rnd$YearVisitArm <- 
  factor(paste(wards.rnd$YearVisit, wards.rnd$Trial_Arm, sep = "-"),
         levels(dat$Trial_Arm.rnd))

# Add mean survey date (only for dogs with primary outcome data)
wards.rnd$mean.date.added <-
  tapply(dat$date.added[!is.na(dat$Vaccinated)], dat$ward.rnd[!is.na(dat$Vaccinated)], 
         function(x) {
           as.character(as.Date(round(mean(x))))
         })[rownames(wards.rnd)]

# Add mean survey date as days since start of trial year
wards.rnd$YearDay <- round(tapply(dat$YearDay, dat$ward.rnd, mean, na.rm = TRUE))

# Impute missing YearDay as mean of wards in that year and visit
wards.rnd.NA.YearDay <- round(tapply(wards.rnd$YearDay, wards.rnd$YearVisit, mean, na.rm = TRUE))
wards.rnd$YearDay[is.na(wards.rnd$YearDay)] <- 
  wards.rnd.NA.YearDay[match(wards.rnd$YearVisit[is.na(wards.rnd$YearDay)], names(wards.rnd.NA.YearDay))]

# Create a survey date that is averaged over wards but not over year or visit. This will
# be used for extrapolating across the year for an average ward.
wards.rnd$YearDay.mean <-
  round(tapply(wards.rnd$YearDay, wards.rnd$YearVisit, mean))[as.character(wards.rnd$YearVisit)]

# Convert day of trial year to month of trial year, rounded to nearest month
wards.rnd$YearMonth.mean <- round(wards.rnd$YearDay.mean * 12 / 365 + 0.5) 

# What was the mean month of each visit
(mean.visit.month <- tapply(wards.rnd$YearMonth.mean, wards.rnd$Visit, mean))

# Add coverage (primary outcome)
wards.rnd[, c("NotVaccinated", "Vaccinated")] <-
  table(dat$ward.rnd, dat$Vaccinated)[rownames(wards.rnd), c("0", "1")]
wards.rnd$N <- wards.rnd$Vaccinated + wards.rnd$NotVaccinated
wards.rnd$Coverage <- wards.rnd$Vaccinated / wards.rnd$N

# Add coverage (secondary outcome)
wards.rnd[, c("NotVaccinated2", "Vaccinated2")] <-
  table(dat$ward.rnd, dat$Vaccinated2)[rownames(wards.rnd), c("0", "1")]
wards.rnd$N2 <- wards.rnd$Vaccinated2 + wards.rnd$NotVaccinated2
wards.rnd$Coverage2 <- wards.rnd$Vaccinated2 / wards.rnd$N2

# Remove rows with no coverage data
wards.rnd.reduced <- wards.rnd[wards.rnd$N > 0, ]
wards.rnd.reduced$ward.short <- 
  factor(sapply(strsplit(as.character(wards.rnd.reduced$ward), "\\-"), "[", 1))


# What proportion of wards were surveyed at each visit, by arm:
cbind(table(wards.rnd.reduced$YearVisit, wards.rnd.reduced$Trial_Arm), 
      Total = table(wards.rnd.reduced$YearVisit),
      `% wards with data` = 
        paste0(round(100 * table(wards.rnd.reduced$YearVisit)/nlevels(wards.rnd.reduced$ward)), "%"))


# Fill in missing subvillage distances where the subvillage distance is available for other
# HH in that SV, otherwise use the mean distance across
sv.without.coords <- levels(droplevels(dat$sub_village[is.na(dat$SV_V_Distance)]))
dat[dat$sub_village %in% sv.without.coords, ]

# Take a look at the subvillage-village distances
sv.v.distances <- tapply(dat$SV_V_Distance, dat$sub_village, function(x) unique(na.omit(x)))
hist(unlist(sv.v.distances), breaks = 50)
summary(unlist(sv.v.distances))
hist(scale(log10(dat$SV_V_Distance)))
dat$SV_V_Distance.log10.scaled <- scale(log10(dat$SV_V_Distance))

# Coverage proportion by arm and visit
covtab <- prop.table(table(dat$Trial_Arm, dat$YearVisit, dat$Vaccinated), 1:2)[, , "1"]
barplot(covtab, beside = TRUE, legend.text = TRUE, 
        args.legend = list(x = "topleft", bty = "n", horiz = TRUE))

# Years vaccinated (only available for visits 5 & 6 = year 3)
table(YearsVaccinated = dat$years.vaccinated, 
      Visit = dat$YearVisit, 
      Vaccinated = dat$Vaccinated, 
      exclude = NULL)

# Ever vaccinated: this round / earlier round / never
year.vax.start.dates <- as.Date(c(Y1 = "2020-11-01", Y2 = "2021-11-01", Y3 = "2022-11-01"))
dat$ever.vacc.calc <- factor(NA, c("Never", "Earlier year", "This year"))

# For year 1, there was no vaccination in previous years
dat$ever.vacc.calc[dat$Year %in% "Y1"] <- 
  levels(dat$ever.vacc.calc)[1 + 2 * dat$Vaccinated[dat$Year %in% "Y1"]]

# For year 2
dat$ever.vacc.calc[dat$Year %in% "Y2"] <-
  levels(dat$ever.vacc.calc)[2 + (as.Date(dat$date.vaccinated)[dat$Year %in% "Y2"] >= year.vax.start.dates["Y2"])]
dat$ever.vacc.calc[dat$Year %in% "Y2" & dat$Vaccinated %in% 0] <- "Never"

# For year 3
dat$ever.vacc.calc[dat$Year %in% "Y3"] <-
  levels(dat$ever.vacc.calc)[2 + (as.Date(dat$date.vaccinated)[dat$Year %in% "Y3"] >= year.vax.start.dates["Y3"])]
dat$ever.vacc.calc[dat$Year %in% "Y3" & dat$Vaccinated %in% 0] <- "Never"

# Check year 3 calculated ever.vacc against survey ever.vacc
dat$ever.vacc.y3 <- factor(NA, levels(dat$ever.vacc.calc))
dat$ever.vacc.y3[dat$Year %in% "Y3"] <- 
  levels(dat$ever.vacc.y3)[2 + (as.integer(substr(dat$years.vaccinated, nchar(dat$years.vaccinated), nchar(dat$years.vaccinated)))[dat$Year %in% "Y3"] > 2.5)]
dat$ever.vacc.y3[dat$Year %in% "Y3" & dat$Vaccinated %in% 0] <- "Never"
table(dat$ever.vacc.y3, exclude = NULL)
table(dat$years.vaccinated[dat$Year %in% "Y3"], 
      dat$ever.vacc.calc[dat$Year %in% "Y3"], exclude = NULL)
table(SurveyQ = dat$ever.vacc.y3, CalcFromDate = dat$ever.vacc.calc)
table(EverVacc = dat$ever.vacc.calc, Vaccinated = dat$Vaccinated, dat$Year, exclude = TRUE)

# Define labels and definitions for ward.rnd variables
label(wards.rnd$ward.rnd) <- "Ward-year-visit"
label(wards.rnd$YearVisit) <- "Year-visit"
label(wards.rnd$Year) <- "Year"
label(wards.rnd$Visit) <- "Visit"
label(wards.rnd$ward) <- "Ward"
attr(wards.rnd$ward, "definition") <- "The initial of the district is appended, to make ward names unique"
label(wards.rnd$district) <- "District"
label(wards.rnd$district.rnd) <- "District-year-visit"
label(wards.rnd$Trial_Arm) <- "Trial arm"
label(wards.rnd$YearVisitArm) <- "Year-visit-arm"
label(wards.rnd$mean.date.added) <- "Mean survey date"
attr(wards.rnd$mean.date.added, "definition") <- "Mean survey date among dogs with primary outcome data"
label(wards.rnd$NotVaccinated) <- "No of dogs not vaccinated (primary outcome)"
label(wards.rnd$Vaccinated) <- "No of dogs vaccinated (primary outcome)"
label(wards.rnd$N) <- "No of dogs surveyed (primary outcome)"
label(wards.rnd$Coverage) <- "Coverage (primary outcome)"
attr(wards.rnd$Coverage, "definition") <- "Vaccinated/N"
label(wards.rnd$NotVaccinated2) <- "No of dogs not vaccinated (secondary outcome)"
label(wards.rnd$Vaccinated2) <- "No of dogs vaccinated (secondary outcome)"
label(wards.rnd$N2) <- "No of dogs surveyed (secondary outcome)"
label(wards.rnd$Coverage2) <- "Coverage (secondary outcome)"
attr(wards.rnd$Coverage2, "definition") <- "Vaccinated2/N2"

# Apply variable names as labels, as a default where a label hasn't already been applied
for(nm in names(dat)) {
  if(label(dat[, nm]) == "") {
    label(dat[, nm]) <- nm
    label(dat[, nm]) <- gsub("_", " ", nm)
  } 
}; rm(nm)
for(nm in names(wards)) {
  if(label(wards[, nm]) == "") {
    label(wards[, nm]) <- nm
    label(wards[, nm]) <- gsub("_", " ", nm)
  }
}; rm(nm) 
for(nm in names(wards.rnd)) {
  if(label(wards.rnd[, nm]) == "") {
    label(wards.rnd[, nm]) <- nm
    label(wards.rnd[, nm]) <- gsub("_", " ", nm)
  } 
}; rm(nm) 

# Check labels
label(dat)
label(wards)
label(wards.rnd)

## Perform analyses ----

### Primary efficacy analysis ----


# Mean vaccination coverage will be compared between the two arms using a
# generalized linear mixed-effects regression model (GLMM) fitted using 
# maximum likelihood. Logit-normal random effects will be fitted to account 
# for variation in coverage between districts, wards, subvillages and households. 
# For each of these levels, except household, additional visit-specific random 
# effects will be included to allow for unexplained variation in coverage that 
# is inconsistent between visits. Intervention arm, assessment visit (2 or 11 months)
# and year will be modelled as fixed effects. All three two-way interactions and a
# three-way interaction between year, visit, and arm will be fitted. The effect 
# of fitting these interactions will be to allow coverage to differ between 
# the six visits (year × visit), and to allow the intervention effect to 
# differ between the six visits (arm × year, arm × visit, arm × year × visit). 

# The primary analysis will be a test of the null hypothesis of equal coverage
# (Arm One (team-based) : Arm Two (community-based) odds ratio = 1) 
# between community-based and team-based delivery against the two-sided alternative hypothesis of 
# unequal coverage. 
# The full model, including the main effect of 
# trial arm and the arm × year, arm × visit, and arm × year × visit interactions, 
# will be compared with a null model not including trial arm or its interactions with visit or year. 
# In effect, the null model allows coverage to vary over the six time points 
# (the six household surveys), while the full model includes an additional 
# six parameters that allow coverage to vary independently between trial arms at 
# each of the six time points. Under the null hypothesis, therefore, there is no
# intervention effect at any of the six time points, while the alternative hypothesis 
# is that there is an intervention effect at one or more time points. 


# Create list for storing analysis data sets and model outputs
primary.analysis <- 
  list(analysis.data = droplevels(dat[!is.na(dat$Vaccinated), ]))

# Fit primary analysis full model
primary.analysis$full.model.formula <- 
  Vaccinated ~ 
    Visit + Year + Trial_Arm + 
    Visit:Year + Visit:Trial_Arm + Year:Trial_Arm + 
    Visit:Year:Trial_Arm +
    (1 | Household.file.id) + 
    (1 | sub_village.rnd) + (1 | sub_village) + 
    (1 | ward.rnd) + (1 | ward) + 
    (1 | district.rnd) + (1 | district)
# If using the alternative definition of vaccination coverage, drop Household.file.id 
# random effect due to lack of information in data on its variance
if(vaccination.definition == "Vaccinated_M4") {
  primary.analysis$full.model.formula <- 
    update(primary.analysis$full.model.formula, ~ . - (1 | Household.file.id))
}
primary.analysis$full.model <- 
  glmmTMB(primary.analysis$full.model.formula, 
          family = binomial, data = primary.analysis$analysis.data)
summary(primary.analysis$full.model)

# Check that the whole analysis data set was used
stopifnot(nrow(primary.analysis$analysis.data) == nrow(model.matrix(primary.analysis$full.model)))

# Fit primary analysis null model
primary.analysis$null.model <- 
  update(primary.analysis$full.model, 
         ~ . - Trial_Arm - Visit:Trial_Arm - Year:Trial_Arm - Visit:Year:Trial_Arm)
summary(primary.analysis$null.model)

# Test intervention effect
primary.analysis$likelihood.ratio.test <-
  anova(primary.analysis$full.model, primary.analysis$null.model)
primary.analysis$likelihood.ratio.test

# Primary analysis p-value
primary.analysis$intervention.p.value <- 
  primary.analysis$likelihood.ratio.test[2, "Pr(>Chisq)"]

# Make table storing treatment effect at each of the six year-visits,
# with estimates and 95% CI
effect.data <- data.frame(YearVisit = factor(levels(dat$YearVisit), levels(dat$YearVisit)))
rownames(effect.data) <- levels(dat$YearVisit)

# Create template data set for storing predicted coverage at each time point by arm
pred.tab <- 
  expand.grid(Trial_Arm = levels(dat$Trial_Arm), 
              Visit = levels(dat$Visit), 
              Year = levels(dat$Year))
pred.tab$YearVisit <- factor(paste(pred.tab$Year, pred.tab$Visit, sep = "-"), levels(dat$YearVisit))
rownames(pred.tab) <- paste(pred.tab$Year, pred.tab$Visit, pred.tab$Trial_Arm, sep = "-")
pred.tab$x <- as.numeric(pred.tab$YearVisit) + (as.numeric(pred.tab$Trial_Arm) - 1.5)/20

# Create matrix for extracting log OR at each time point 
X <- model.matrix(formula(primary.analysis$full.model, fixed.only = TRUE), 
                  data = cbind(pred.tab, Vaccinated = 1))
A <- X[pred.tab$Trial_Arm == "Community-based", ] - X[pred.tab$Trial_Arm == "Team-based", ]

# Calculate ORs, their SEs, 95% CIs, and p-values
# Adjust ORs using Zeger's approximate bias correction, which is OR^z, or equivalently exp(logOR*z), where:
primary.analysis$re.var <- unlist(VarCorr(primary.analysis$full.model)$cond)
V <- sum(primary.analysis$re.var)
z <- 1/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V)
effect.data$InterventionOR.log <- (A %*% fixef(primary.analysis$full.model)$cond)[, 1]
effect.data$InterventionOR.log.se <- sqrt(diag(A %*% vcov(primary.analysis$full.model)$cond %*% t(A)))
effect.data$InterventionOR.log.ci.lo <- effect.data$InterventionOR.log - qnorm(0.975) * effect.data$InterventionOR.log.se
effect.data$InterventionOR.log.ci.hi <- effect.data$InterventionOR.log + qnorm(0.975) * effect.data$InterventionOR.log.se
effect.data[, c("InterventionOR", "InterventionOR.ci.lo", "InterventionOR.ci.hi")]  <- 
  exp(effect.data[, c("InterventionOR.log", "InterventionOR.log.ci.lo", "InterventionOR.log.ci.hi")] * z)
effect.data$Intervention.p.value <- 2 * (1 - pnorm(abs(effect.data$InterventionOR.log)/effect.data$InterventionOR.log.se))
effect.data$x <-  as.numeric(effect.data$YearVisit)

# Add effect estimate table to primary analysis list
primary.analysis$effect.estimates <- effect.data
rm(X, A)

# Calculate predicted coverage at each of the six time points in each arm, with 95% CIs,
# adjusted for Jensen's inequality
full.model.prediction <- predict(primary.analysis$full.model, newdata = pred.tab, 
                                 re.form = ~ 0, type = "link", se.fit = TRUE)
pred.tab$logodds <- full.model.prediction$fit
pred.tab$logodds.se <- full.model.prediction$se.fit
pred.tab$logodds.ci.lo <- pred.tab$logodds - qnorm(0.975) * pred.tab$logodds.se
pred.tab$logodds.ci.hi <- pred.tab$logodds + qnorm(0.975) * pred.tab$logodds.se
pred.tab$Coverage <- jensen.logit.adjust(plogis(pred.tab$logodds), V = V, method = "zeger")
pred.tab$Coverage.ci.lo <- jensen.logit.adjust(plogis(pred.tab$logodds.ci.lo), V = V, method = "zeger")
pred.tab$Coverage.ci.hi <- jensen.logit.adjust(plogis(pred.tab$logodds.ci.hi), V = V, method = "zeger")

# Add simple mean coverage, for comparison
pred.tab$Coverage.simple[pred.tab$Trial_Arm == "Team-based"] <- 
  covtab["Team-based", ]
pred.tab$Coverage.simple[pred.tab$Trial_Arm == "Community-based"] <- 
  covtab["Community-based", ]


### Secondary efficacy analyses ----

# 1.
# Vaccination coverage at the eleven-month time point within each annual vaccination
#cycle will be compared between arms in order to test the hypothesis that eleven-month
# vaccination coverage in the continous delivery arm is higher than the coverage at 
# the same time-point in the Team-based delivery arm. A model including the fixed effects 
# of year, trial arm, and their interaction will be compared with a model with only year 
# as a fixed effect. Both models will be fitted to data from the second (11-month) visit 
# only, excluding data from the first (2-month) visit.

# Create list for storing model outputs
secondary.analyses <- 
  list(no1 = list(),
       no2 = list(),
       no3 = list(),
       no4 = list())

# Fit secondary analysis 1 null model.
# This is the full model from the primary analysis, reparameterised so that 
# the single fixed effect is a 9-category factor combining Year, Visit and Arm,
# but lacking terms that allow coverage to differ between arms at the second visit 
# of each year
secondary.analyses$no1$null.model <- 
  update(primary.analysis$null.model, 
         ~ . + 
           I(Trial_Arm.rnd == "Y1-V1-Community-based") + 
           I(Trial_Arm.rnd == "Y2-V1-Community-based") + 
           I(Trial_Arm.rnd == "Y3-V1-Community-based"))


# Check that the two analysis data sets match in size
stopifnot(nrow(model.matrix(secondary.analyses$no1$null.model)) == 
            nrow(model.matrix(primary.analysis$full.model)))

# Test secondary analysis 1 intervention effect
secondary.analyses$no1$likelihood.ratio.test <-
  anova(primary.analysis$full.model, secondary.analyses$no1$null.model)

# Secondary analysis 1 p-value
secondary.analyses$no1$intervention.p.value <- 
  secondary.analyses$no1$likelihood.ratio.test[2, "Pr(>Chisq)"]

# 2.	
# The full model will be used to compare the likelihood of coverage dipping below 
# the coverage threshold in each arm by implementing the following plan:
#  a.	Use parametric bootstrapping from the fitted full model to generate nsim samples
#      from the sampling distribution of the model parameters, capturing uncertainty in 
#      the model parameter estimates (fixed effects and random effect variances).
#  b.	For each of the nsim sets of parameters, simulate a new set of wards with 
#      ward-level coverage, including random effect variation, giving the distribution
#      of predicted coverage in each arm for a random ward, taking account of parameter 
#      estimation uncertainty and random effect variation.
#  c.	Calculate the proportion of wards with coverage below threshold for each arm.


# Because we want to predict coverage across the whole year, we need to fit a new model
# where visits 1 and 2 are replaced by the exact day of the year when the dog was surveyed.
# The vaccination year starts on November 1st and ends on October 31st.
secondary.analyses$no2$full.model.formula <- 
  Vaccinated ~ 
    YearDay + Year + Trial_Arm + 
    YearDay:Year + YearDay:Trial_Arm + Year:Trial_Arm + 
    YearDay:Year:Trial_Arm +
    (1 | Household.file.id) + 
    (1 | sub_village) + (1 | sub_village.rnd) + 
    (1 | ward) + (1 | ward.rnd) + 
    (1 | district) + (1 | district.rnd)
# If using the alternative definition of vaccination coverage, drop Household.file.id 
# random effect due to lack of information in data on its variance
if(vaccination.definition == "Vaccinated_M4") {
  secondary.analyses$no2$full.model.formula <- 
    update(secondary.analyses$no2$full.model.formula, ~ . - (1 | Household.file.id))
}

secondary.analyses$no2$full.model <- 
  glmmTMB(secondary.analyses$no2$full.model.formula, 
          family = binomial, data = primary.analysis$analysis.data)
logLik(secondary.analyses$no2$full.model)
summary(secondary.analyses$no2$full.model)

# Test the 3-way interaction
secondary.analyses$no2$no.3way.ixn.model <- 
  update(secondary.analyses$no2$full.model, ~ . - YearDay:Year:Trial_Arm)
secondary.analyses$no2$no.3way.ixn.p.value <- 
  anova(secondary.analyses$no2$full.model, secondary.analyses$no2$no.3way.ixn.model)[2, "Pr(>Chisq)"]
secondary.analyses$no2$drop.2way.ixns <-
  drop1(secondary.analyses$no2$no.3way.ixn.model, test = "Chisq")
secondary.analyses$no2$no.2way.ixn.p.values <- 
  secondary.analyses$no2$drop.2way.ixns[-1, "Pr(>Chi)"]
names(secondary.analyses$no2$no.2way.ixn.p.values) <- rownames(secondary.analyses$no2$drop.2way.ixns)[-1]

# Drop non-significant 2-way interaction
secondary.analyses$no2$final.model <- 
  update(secondary.analyses$no2$no.3way.ixn.model, ~ . - YearDay:Year)

# Run this line to show that all remaining higher order terms are significant at P < 0.05
drop1(secondary.analyses$no2$final.model, test = "Chisq")

# Create data set for predicting monthly coverage by arm and year
newdata.by.month.and.arm <-
  expand.grid(YearMonth = 1:12,
              Year = levels(dat$Year),
              ward = levels(dat$ward))
newdata.by.month.and.arm$Trial_Arm <- 
  wards$Trial_Arm[match(newdata.by.month.and.arm$ward, wards$name)]
dim(newdata.by.month.and.arm)

# Calculate day from month
newdata.by.month.and.arm$YearDay <- 
  round(365 * (newdata.by.month.and.arm$YearMonth - 0.5) / 12)
newdata.by.month.and.arm$YearMonth <- 
  factor(newdata.by.month.and.arm$YearMonth, 1:12)

# Add visit and district
newdata.by.month.and.arm$Visit <- 
  factor(newdata.by.month.and.arm$YearDay > 365/2, c(FALSE, TRUE), levels(dat$Visit))
newdata.by.month.and.arm$ward.rnd <- 
  factor(paste(newdata.by.month.and.arm$ward, 
               newdata.by.month.and.arm$Year, 
               newdata.by.month.and.arm$Visit, sep = "-"),
         levels(dat$ward.rnd))
newdata.by.month.and.arm$district <- 
  wards[as.character(newdata.by.month.and.arm$ward), "district"]
newdata.by.month.and.arm$district.rnd <- 
  wards.rnd[as.character(newdata.by.month.and.arm$ward.rnd), "district.rnd"]

# Add ward-Year 
newdata.by.month.and.arm$ward.Year <- 
  paste(newdata.by.month.and.arm$ward, newdata.by.month.and.arm$Year, sep = "-")

# Choose the number of simulations
secondary.analyses$no2$nsim <- max(10, 1000 * iter.mult)
# 1000 takes ~14 min on laptop (Apple M1 Max with 10 cores)

# Choose which random effects to include in the simulations.
# We are predicting coverage at the ward-6-visit level, 
# so need to include all random effects at that level and above
which.re.to.include <- c("ward.rnd", "ward", "district.rnd", "district")

# Immediately before simulating from the model, set a random number seed using
# https://www.random.org/integers/?num=1&min=0&max=1000000000&col=1&base=10&format=html&rnd=new
# Here are your random numbers: 606505365
# Timestamp: 2024-04-22 12:58:09 UTC
if(set.rng.seed) set.seed(606505365)

# Simulate nsim responses from the primary analysis model, including random effect variation
sim.response1 <- simulate(primary.analysis$full.model, nsim = secondary.analyses$no2$nsim)
dim(sim.response1)

# Simulate nsim responses from the secondary analysis #2 model, including random effect variation
sim.response2 <- simulate(secondary.analyses$no2$final.model, nsim = secondary.analyses$no2$nsim)
dim(sim.response2)

# Get no of vaccinated dogs from first simulation, to check consistency with set seed
sum(sim.response1[[1]][, 1])

# First set RNG seed
# Here are your random numbers: 841687013
# Timestamp: 2024-04-22 13:12:03 UTC
if(set.rng.seed) set.seed(841687013)

# refit the primary analysis full model to each simulated response vector,
# looping nsim times 
# within each loop, get predicted coverage at each time point in each arm, using pred.tab as a template
sim.fits.list <-
  mclapply(1:secondary.analyses$no2$nsim, function(i) {
    
    # print progress
    if(i/secondary.analyses$no2$nsim == round(i/secondary.analyses$no2$nsim, 2)) 
      message_parallel(paste0(100 * (i/secondary.analyses$no2$nsim), "%"))
    
    # 2a. Parametric bootstrapping
    
    # fit full model to simulated data
    fit1 <- refit(primary.analysis$full.model, newresp = sim.response1[[i]][, 1])
    
    # get random effect variance estimates
    re.est1 <- unlist(VarCorr(fit1)$cond)
    
    # get predicted coverage at each time point in each arm
    pred.cov.by.survey.and.arm <- 
      jensen.logit.adjust(predict(fit1, newdata = pred.tab, re.form = ~ 0, type = "response"),
                          V = sum(re.est1), method = "zeger")
    
    # Calculate treatment effect at each year-visit
    stopifnot(all(pred.tab$YearVisit[pred.tab$Trial_Arm == "Community-based"] == 
                    pred.tab$YearVisit[pred.tab$Trial_Arm == "Team-based"]))
    intervention.effect.abs <- 
      pred.cov.by.survey.and.arm[pred.tab$Trial_Arm == "Community-based"] - 
      pred.cov.by.survey.and.arm[pred.tab$Trial_Arm == "Team-based"]
    intervention.effect.oddsratio <- 
      exp(qlogis(pred.cov.by.survey.and.arm[pred.tab$Trial_Arm == "Community-based"]) - 
            qlogis(pred.cov.by.survey.and.arm[pred.tab$Trial_Arm == "Team-based"]))
    
    #  2b. Simulate a new set of wards with ward-level coverage, including random effect variation
    
    # fit full model to simulated data
    fit2 <- refit(secondary.analyses$no2$final.model, newresp = sim.response2[[i]][, 1])
    
    # get random effect variance estimates
    re.est2 <- unlist(VarCorr(fit2)$cond)
    
    # Get fixed effect predictions for this parametric bootstrap sample for each
    # ward
    predicted.logit.coverage <- 
      predict(fit2, re.form = ~ 0, newdata = wards.rnd[, c("Year", "YearDay", "Trial_Arm")])
    # Add variation at every random effect level down to ward.
    # Note that these are predictions for new unsampled wards, not sampled wards.
    for(re in which.re.to.include) {
      print(re)
      re.residuals <- rnorm(nlevels(wards.rnd[, re]), sd = sqrt(re.est2[re]))
      predicted.logit.coverage <- predicted.logit.coverage + re.residuals[as.numeric(wards.rnd[, re])]
    }
    
    # The vector predicted.logit.coverage contains log odds predicted coverage on the mean
    # survey date for each ward at each of the six visits (two time points per year). We want to 
    # extrapolate and interpolate from those time points to the midpoint of each month, 
    # for each ward-year.
    newdata.by.month.and.arm$predicted.logit.coverage <- NA
    for(wy in levels(wards.rnd$ward.Year)) {
      Y <- predicted.logit.coverage[wards.rnd$ward.Year %in% wy]
      X <- wards.rnd$YearDay.mean[wards.rnd$ward.Year %in% wy]
      new.X <- newdata.by.month.and.arm$YearDay[newdata.by.month.and.arm$ward.Year %in% wy]
      newdata.by.month.and.arm$predicted.logit.coverage[newdata.by.month.and.arm$ward.Year %in% wy] <- 
        predict(lm(Y ~ X), newdata = data.frame(X = new.X))
    }
    pred.cov.by.month.and.arm <- newdata.by.month.and.arm$predicted.logit.coverage
    
    # adjust the predicted coverages to correct for bias due to Jensen's inequality
    pred.cov.by.month.and.arm <-  
      jensen.logit.adjust(plogis(pred.cov.by.month.and.arm),
                          V = sum(re.est2[!names(re.est2) %in% which.re.to.include]),
                          method = "zeger")
    list(pred.cov.by.survey.and.arm = pred.cov.by.survey.and.arm,
         pred.cov.by.month.and.arm = pred.cov.by.month.and.arm, 
         intervention.effect.abs = intervention.effect.abs,
         intervention.effect.oddsratio = intervention.effect.oddsratio)
  }, mc.cores = detectCores() - 1)


# Extract tables of predicted coverage
pred.cov.by.month.and.arm <- 
  do.call("cbind", lapply(sim.fits.list, "[[", "pred.cov.by.month.and.arm"))
dim(pred.cov.by.month.and.arm)
pred.cov.by.survey.and.arm <- 
  do.call("cbind", lapply(sim.fits.list, "[[", "pred.cov.by.survey.and.arm"))
dim(pred.cov.by.survey.and.arm)

# sum elements of both tables to check reproducibility when set.rng.seed is TRUE
round(sum(pred.cov.by.month.and.arm) + sum(pred.cov.by.survey.and.arm))

# pred.cov.by.month.and.arm is a table with nsim columns and
# n_year x n_month x n_arm rows
nlevels(newdata.by.month.and.arm$Year) * 
  nlevels(newdata.by.month.and.arm$ward) * 
  nlevels(newdata.by.month.and.arm$YearMonth)
secondary.analyses$no2$nsim
dim(pred.cov.by.month.and.arm)

# The rows align with the year-month-arm rows of the newdata.by.month.and.arm data frame.
# Each column is a simulated coverage for a *new* ward at a given year and month,
# in a given trial arm. 

# pred.cov.by.survey.and.arm is a table with nsim columns and
# n_arm x n_visit x n_year rows
nlevels(wards.rnd$Trial_Arm) * nlevels(wards.rnd$Visit) * nlevels(wards.rnd$Year) 
secondary.analyses$no2$nsim
dim(pred.cov.by.survey.and.arm)

# For comparison with the max likelihood effect estimates and Wald 95% CIs, 
# calculate mean coverage and 95% CI at each time point in each arm, 
# adjusted for Jensen's inequality
pred.tab$CoveragePB <- apply(pred.cov.by.survey.and.arm, 1, mean)
pred.tab$CoveragePB.ci.lo <- apply(pred.cov.by.survey.and.arm, 1, quantile, 0.025)
pred.tab$CoveragePB.ci.hi <- apply(pred.cov.by.survey.and.arm, 1, quantile, 0.975)


# plot coverage estimates
estimate.ci.plot(x.term = "x", y.term = "CoveragePB", 
                 y.ci.lo = "CoveragePB.ci.lo", y.ci.hi = "CoveragePB.ci.hi", 
                 plotdata = pred.tab, xlab = "Visit", ylab = "Coverage", 
                 pch = 20 + as.numeric(pred.tab$Trial_Arm), 
                 ylim = 0:1,
                 bg = as.numeric(pred.tab$Trial_Arm) + 4, 
                 main = "Parametric bootstrap estimates",
                 x.at = unique(round(pred.tab$x)), x.labels = levels(pred.tab$YearVisit))
legend("topleft", legend = levels(dat$Trial_Arm), pch = 20 + 1:nlevels(dat$Trial_Arm),
       pt.bg = 1:nlevels(dat$Trial_Arm) + 4)

estimate.ci.plot(x.term = "x", y.term = "Coverage", y.ci.lo = "Coverage.ci.lo", 
                 y.ci.hi = "Coverage.ci.hi", 
                 plotdata = pred.tab, xlab = "Visit", ylab = "Coverage", 
                 pch = 20 + as.numeric(pred.tab$Trial_Arm), 
                 ylim = 0:1,
                 bg = as.numeric(pred.tab$Trial_Arm) + 4, 
                 main = "Mean and Wald CI",
                 x.at = unique(round(pred.tab$x)), x.labels = levels(pred.tab$YearVisit))
legend("topleft", legend = levels(dat$Trial_Arm), pch = 20 + 1:nlevels(dat$Trial_Arm),
       pt.bg = 1:nlevels(dat$Trial_Arm) + 4)

estimate.ci.plot(x.term = "x", y.term = "Coverage.simple", 
                 plotdata = pred.tab, xlab = "Visit", ylab = "Coverage", 
                 pch = 20 + as.numeric(pred.tab$Trial_Arm), 
                 ylim = 0:1,
                 bg = as.numeric(pred.tab$Trial_Arm) + 4, 
                 main = "Simple mean coverage",
                 x.at = unique(round(pred.tab$x)), x.labels = levels(pred.tab$YearVisit))
legend("topleft", legend = levels(dat$Trial_Arm), pch = 20 + 1:nlevels(dat$Trial_Arm),
       pt.bg = 1:nlevels(dat$Trial_Arm) + 4)

# The estimates and CI are almost identical between the standard mean and Wald CIs and the bootstrapped CIs,
# and the Wald CIs are much faster to produce, so stick with these.
pred.tab[, c("Coverage.simple", "Coverage", "CoveragePB")]
par(mfrow = c(2, 2))
plot(CoveragePB ~ Coverage, data = pred.tab, xlim = 0:1, ylim = 0:1)
abline(0, 1)
plot(Coverage.simple ~ Coverage, data = pred.tab, xlim = 0:1, ylim = 0:1)
abline(0, 1)
plot(Coverage.simple ~ CoveragePB, data = pred.tab, xlim = 0:1, ylim = 0:1)
abline(0, 1)
par(mfrow = c(1, 1))

# Delete all the redundant coverage estimates
pred.tab$logodds <- pred.tab$logodds.se <- pred.tab$logodds.ci.lo <- 
  pred.tab$logodds.ci.hi  <- pred.tab$CoveragePB <- pred.tab$CoveragePB.ci.lo <- 
  pred.tab$CoveragePB.ci.hi <- pred.tab$Coverage.simple <- NULL

# Add the predicted coverages to the effect estimates table, and format the numbers for output
primary.analysis$effect.estimates$Team <-
  paste0(my.format(pred.tab$Coverage[pred.tab$Trial_Arm == "Team-based"], ndp = 2), " (",
         my.format(pred.tab$Coverage.ci.lo[pred.tab$Trial_Arm == "Team-based"], ndp = 2), ", ",
         my.format(pred.tab$Coverage.ci.hi[pred.tab$Trial_Arm == "Team-based"], ndp = 2), ")") 
primary.analysis$effect.estimates$Community <-
  paste0(my.format(pred.tab$Coverage[pred.tab$Trial_Arm == "Community-based"], ndp = 2), " (",
         my.format(pred.tab$Coverage.ci.lo[pred.tab$Trial_Arm == "Community-based"], ndp = 2), ", ",
         my.format(pred.tab$Coverage.ci.hi[pred.tab$Trial_Arm == "Community-based"], ndp = 2), ")") 
primary.analysis$effect.estimates$`Odds ratio (95% CI)` <- 
  paste0(my.format(primary.analysis$effect.estimates$InterventionOR, ndp = 2), " (",
         my.format(primary.analysis$effect.estimates$InterventionOR.ci.lo, ndp = 2), ", ",
         my.format(primary.analysis$effect.estimates$InterventionOR.ci.hi, ndp = 2), ")") 
primary.analysis$effect.estimates$`<i>p</i>-value` <- 
  p.format(primary.analysis$effect.estimates$Intervention.p.value)
primary.analysis$effect.estimates$`Year-visit` <- primary.analysis$effect.estimates$YearVisit
primary.analysis$effect.estimates <- 
  primary.analysis$effect.estimates[, c("Year-visit", "Team", "Community", "Odds ratio (95% CI)", "<i>p</i>-value")]

# Add the coverage estimates to the primary analysis list
primary.analysis$predicted.coverage <- pred.tab

# Calculate mean and 95% CI coverage across the three years of the trial from
# the parametric bootstrap samples
boot.mean.coverage.by.arm <-
  apply(pred.cov.by.month.and.arm, 2, 
        function(x) tapply(x, newdata.by.month.and.arm$Trial_Arm, mean))
mean.coverage <- my.format(apply(boot.mean.coverage.by.arm, 1, mean), ndp = 2)
mean.coverage.ci95 <- 
  paste0(mean.coverage, " (",
         apply(my.format(apply(boot.mean.coverage.by.arm, 1, quantile, c(0.025, 0.975)), ndp = 2), 
               2, paste, collapse = ", "), ")")
names(mean.coverage.ci95) <- sapply(strsplit(names(mean.coverage), "-"), "[", 1)

# Calculate mean and 95% odds ratios across the three years of the trial from
# the parametric bootstrap samples
boot.log.odds.ratio <-
  apply(boot.mean.coverage.by.arm, 2, function(x) diff(qlogis(x)))
boot.odds.ratio.ci95 <- 
  paste0(my.format(exp(mean(boot.log.odds.ratio)), ndp = 2), " (",
         paste0(my.format(exp(quantile(boot.log.odds.ratio, c(0.025, 0.975))), ndp = 2), collapse = ", "),
         ")")
boot.log.odds.ratio.p.value <- 
  p.format(2 * min(c(mean(boot.log.odds.ratio > 0), mean(boot.log.odds.ratio < 0))))

primary.analysis$effect.estimates.allyears <- 
  rbind(c(`Year-visit` = "All years", 
          mean.coverage.ci95,
          `Odds ratio (95% CI)` = boot.odds.ratio.ci95,
          `<i>p</i>-value` = boot.log.odds.ratio.p.value))

# Secondary 2c. Calculate the proportion of wards at each time point with coverage below threshold for each arm.

# Extract from the simulation results the distribution of coverage in each ward at each visit 
coverage.dist <-
  data.frame(Coverage = c(pred.cov.by.month.and.arm), 
             YearMonth = newdata.by.month.and.arm$YearMonth,
             Year = newdata.by.month.and.arm$Year,
             Trial_Arm = newdata.by.month.and.arm$Trial_Arm,
             Simulation = rep(1:ncol(pred.cov.by.month.and.arm), each = nrow(pred.cov.by.month.and.arm)))
head(coverage.dist)
dim(coverage.dist)
table(coverage.dist$Simulation, coverage.dist$YearMonth, coverage.dist$Trial_Arm)

# Get mean coverage by month, with 95% CI from parametric bootstrap
newdata.by.month.and.arm$Coverage.mean <- apply(pred.cov.by.month.and.arm, 1, mean)
newdata.by.month.and.arm$Coverage.ci.lo <- apply(pred.cov.by.month.and.arm, 1, quantile, 0.025)
newdata.by.month.and.arm$Coverage.ci.hi <- apply(pred.cov.by.month.and.arm, 1, quantile, 0.975)
secondary.analyses$no2$cov.dist.by.month.and.arm.mean <- 
  tapply(newdata.by.month.and.arm$Coverage.mean, 
         newdata.by.month.and.arm[, c("YearMonth", "Trial_Arm")], 
         mean)
secondary.analyses$no2$cov.dist.by.month.and.arm.ci.lo <- 
  tapply(newdata.by.month.and.arm$Coverage.ci.lo, 
         newdata.by.month.and.arm[, c("YearMonth", "Trial_Arm")], 
         mean)
secondary.analyses$no2$cov.dist.by.month.and.arm.ci.hi <- 
  tapply(newdata.by.month.and.arm$Coverage.ci.hi, 
         newdata.by.month.and.arm[, c("YearMonth", "Trial_Arm")], 
         mean)

secondary.analyses$no2$cov.dist.table <- 
  secondary.analyses$no2$cov.dist.by.month.and.arm.mean
for(arm in levels(dat$Trial_Arm)) {
  secondary.analyses$no2$cov.dist.table[, arm] <-
    my.format(secondary.analyses$no2$cov.dist.by.month.and.arm.mean[, arm], ndp = 2)
  secondary.analyses$no2$cov.dist.table[, arm] <- 
    paste0(secondary.analyses$no2$cov.dist.table[, arm], 
           " (",
           my.format(secondary.analyses$no2$cov.dist.by.month.and.arm.ci.lo[, arm], ndp = 2),
           ", ",
           my.format(secondary.analyses$no2$cov.dist.by.month.and.arm.ci.hi[, arm], ndp = 2),
           ")")
}

# Get mean coverage < threshold by month, with 95% CI from parametric bootstrap
Coverage.lt40.boot <- 
  sapply(1:ncol(pred.cov.by.month.and.arm), function(j) {
    tapply(pred.cov.by.month.and.arm[, j] < cov.thresh, newdata.by.month.and.arm[, c("YearMonth", "Trial_Arm")], mean)
  })
secondary.analyses$no2$cov.lt40.by.month.and.arm.mean <- 
  matrix(apply(Coverage.lt40.boot, 1, mean), ncol = nlevels(dat$Trial_Arm),
         dimnames = dimnames(secondary.analyses$no2$cov.dist.by.month.and.arm.mean))
secondary.analyses$no2$cov.lt40.by.month.and.arm.ci.lo <- 
  matrix(apply(Coverage.lt40.boot, 1, quantile, 0.025), ncol = nlevels(dat$Trial_Arm),
         dimnames = dimnames(secondary.analyses$no2$cov.dist.by.month.and.arm.mean))
secondary.analyses$no2$cov.lt40.by.month.and.arm.ci.hi <- 
  matrix(apply(Coverage.lt40.boot, 1, quantile, 0.975), ncol = nlevels(dat$Trial_Arm),
         dimnames = dimnames(secondary.analyses$no2$cov.dist.by.month.and.arm.mean))
secondary.analyses$no2$cov.lt40.table <- 
  secondary.analyses$no2$cov.lt40.by.month.and.arm.mean
for(arm in levels(dat$Trial_Arm)) {
  secondary.analyses$no2$cov.lt40.table[, arm] <-
    my.format(secondary.analyses$no2$cov.lt40.by.month.and.arm.mean[, arm], ndp = 2)
  secondary.analyses$no2$cov.lt40.table[, arm] <- 
    paste0(secondary.analyses$no2$cov.lt40.table[, arm], 
           " (",
           my.format(secondary.analyses$no2$cov.lt40.by.month.and.arm.ci.lo[, arm], ndp = 2),
           ", ",
           my.format(secondary.analyses$no2$cov.lt40.by.month.and.arm.ci.hi[, arm], ndp = 2),
           ")")
}

# Estimate the probability of a random (new) ward being below the coverage threshold
lt.thresh.Trial_Arm <- 
  data.frame(
    Estimate = 
      tapply(coverage.dist$Coverage < cov.thresh, coverage.dist[, c("Trial_Arm")], mean),
    Estimate.ci.lo = 
      apply(tapply(coverage.dist$Coverage < cov.thresh, coverage.dist[, c("Trial_Arm", "Simulation")], mean), 1, 
            quantile, 0.025),
    Estimate.ci.hi = 
      apply(tapply(coverage.dist$Coverage < cov.thresh, coverage.dist[, c("Trial_Arm", "Simulation")], mean), 1, 
            quantile, 0.975))

# This is a single estimate of P(coverage < threshold) per arm, averaged over the fixed effects of 
# Month and Year and the random effects of ward, ward.rnd, district, district.rnd
lt.thresh.Trial_Arm

# Add these proportion estimates to the secondary analyses no 2 list
secondary.analyses$no2$prop.estimates.overall <-
  data.frame(rownames(lt.thresh.Trial_Arm),
             paste0(my.format(lt.thresh.Trial_Arm$Estimate, ndp = 2),  " (",
                    my.format(lt.thresh.Trial_Arm$Estimate.ci.lo, ndp = 2), ", ",
                    my.format(lt.thresh.Trial_Arm$Estimate.ci.hi, ndp = 2), ")"))
names(secondary.analyses$no2$prop.estimates.overall) <- c("Trial Arm", paste0("P(Coverage&lt;", cov.thresh, ") (95% CI)"))
secondary.analyses$no2$prop.estimates.overall

# Secondary 3

# Test hypothesis that variation in vaccination coverage over the three-year time
# frame of the study was lower in the community-based arm than in the team-based arm.
# Do this by comparing the variance of the five differences in predicted coverages 
# over the six visits. A high variance would be expected if the amount of change 
# from one visit to the next fluctuated a lot. A steady trend, whether increasing, decreasing, 
# or flat, would give a low variance (i.e. the data is effectively de-trended). 
# The test statistic is the ratio of these variances, community-based:team-based, where the null
# value is 1. 

# Function to calculate the ratio of mean absolute differences in log odds of coverage
# along a vector between two groups
var.diff.ratio <-
  function(x, arm) {
    or <- exp(tapply(x, arm, function(y) mean(abs(diff(qlogis(y))))))
    c(or, OR.Difference = as.vector(exp(diff(log(or)))))
  }

# Calculate the point estimate of the community-based:team-based ratio
secondary.analyses$no3$var.diff.ratio.community.to.team <- 
  var.diff.ratio(x = pred.tab$Coverage, arm = pred.tab$Trial_Arm)

# Calculate the 95% CI from the bootstrapped predicted coverages
secondary.analyses$no3$var.diff.ratio.community.to.team.ci <- 
  t(apply(apply(pred.cov.by.survey.and.arm, 2, var.diff.ratio, pred.tab$Trial_Arm), 1, quantile, c(0.025, 0.975)))
round(cbind(Estimate = secondary.analyses$no3$var.diff.ratio.community.to.team,
            secondary.analyses$no3$var.diff.ratio.community.to.team.ci), 2)


# Secondary 4

# Test hypothesis that variation in vaccination coverage among villages/wards 
# was lower in the continous arm than in the team arm.
# Do this by allowing the inter-ward random effect variance 
# to vary by arm - i.e. estimate arm-specific residual variances at the 
# ward level. Use MCMC in brms, which is a front end for fitting
# GLMMs using Stan (https://mc-stan.org/)


# Define the model formula
secondary.analyses$no4$formula <- 
  brmsformula(Vaccinated ~ 
                Intercept +
                Visit + Year + Trial_Arm + 
                Visit:Year + Visit:Trial_Arm + Year:Trial_Arm + 
                Visit:Year:Trial_Arm +
                (1 | Household.file.id) + 
                (1 | sub_village.rnd) + 
                (1 | sub_village) + 
                (1 | ward.rnd) + 
                (1 | gr(ward, by = Trial_Arm)) + 
                (1 | district.rnd) + 
                (1 | district) - 1,
              family = bernoulli())
# If using the alternative definition of vaccination coverage, drop Household.file.id 
# random effect due to lack of information in data on its variance
if(vaccination.definition == "Vaccinated_M4") {
  secondary.analyses$no4$formula <- 
    update(secondary.analyses$no4$formula, ~ . - (1 | Household.file.id))
}


# Fit the model using brm

# Set random seed using this URL:
# https://www.random.org/integers/?num=1&min=0&max=1000000000&col=1&base=10&format=html&rnd=new
# Here are your random numbers: 517169365
# Timestamp: 2024-05-01 10:10:13 UTC
secondary.analyses$no4$rng.seed <-
  ifelse(set.rng.seed, 517169365, NA) 

secondary.analyses$no4$full.model <- 
  brm(secondary.analyses$no4$formula, 
      data = primary.analysis$analysis.data, 
      prior = prior(normal(0, 10), class = b) + 
        prior(cauchy(0, 2), class = sd),
      iter = 2000 * iter.mult,
      cores = detectCores(),
      seed = secondary.analyses$no4$rng.seed)


# Model diagnostics, optionally
if(FALSE) {
  
  # Check the fit of the rpimary analysis model
  simulateResiduals(primary.analysis$full.model, plot = TRUE, n = 1000, seed = 606505365)
  
  # Check visually that the chains look healthy, an informal check
  plot(secondary.analyses$no4$full.model)
  
  # PPC using DHARMa
  model.check <- 
    createDHARMa(
      simulatedResponse = t(posterior_predict(secondary.analyses$no4$full.model)),
      observedResponse = primary.analysis$analysis.data$Vaccinated,
      fittedPredictedResponse = apply(t(posterior_epred(secondary.analyses$no4$full.model)), 1, mean),
      integerResponse = TRUE)
  plot(model.check)
  
  # PPC using pp_check, ppc_bars_grouped
  pp_check(secondary.analyses$no4$full.model, ndraws = min(100, 10 * iter.mult))
  ppc_bars_grouped(y = as.vector(primary.analysis$analysis.data$Vaccinated),
                   yrep = posterior_predict(secondary.analyses$no4$full.model, ndraws = min(100, 10 * iter.mult)), 
                   group = primary.analysis$analysis.data$Trial_Arm.rnd,
                   prob = 0.95, freq = FALSE)
}

# Extract posterior MCMC sample 
secondary.analyses$no4$mcmc <- 
  as.mcmc(secondary.analyses$no4$full.model, combine_chains = TRUE)
hist(secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmCommunity-based"]^2 / 
       secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmTeam-based"]^2, xlim = c(0, 20),
     breaks = 400000)

# Get median and 95% CrI from the posterior variance ratio
secondary.analyses$no4$var.ratio.median <-
  median(secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmCommunity-based"]^2 / 
           secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmTeam-based"]^2)
secondary.analyses$no4$var.ratio.ci95 <-
  HPDinterval(secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmCommunity-based"]^2 / 
                secondary.analyses$no4$mcmc[, "sd_ward__Intercept:Trial_ArmTeam-based"]^2)

# Collect secondary analyses 3 and 4 in the secondary analyses list
secondary.analyses$no3no4$Estimates <-
  structure(
    c(paste0(my.format(secondary.analyses$no3$var.diff.ratio.community.to.team["OR.Difference"], ndp = 2), " (",
             my.format(secondary.analyses$no3$var.diff.ratio.community.to.team.ci["OR.Difference", "2.5%"], ndp = 2), ", ",
             my.format(secondary.analyses$no3$var.diff.ratio.community.to.team.ci["OR.Difference", "97.5%"], ndp = 2), ")"),
      paste0(my.format(secondary.analyses$no4$var.ratio.median, ndp = 2), " (",
             my.format(secondary.analyses$no4$var.ratio.ci95[, "lower"], ndp = 2), ", ",
             my.format(secondary.analyses$no4$var.ratio.ci95[, "upper"], ndp = 2), ")")),
    dim = 2:1,
    dimnames = 
      list(c("Secondary 3 (variation over time)", "Secondary 4 (variation over space)"),
           "Community-based:team-based ratio (95% CI)"))

# Additional secondary analyses not specified in SAP:

# Secondary 5

# Compare change in coverage from V1 to V2 between arms 
# (does coverage diverge between the arms from V1 to V2, 
# e.g. fall in the team-based arm and not in the Community-based arm?) 
# This is an interaction between arm and visit.

# First need to ask: is this divergence in coverage consistent over the three years? 
# This is an interaction between arm, visit and year. 
secondary.analyses$no5$no.3way.ixn.model <- 
  update(primary.analysis$full.model, ~ . - Visit:Year:Trial_Arm)
secondary.analyses$no5$likelihood.ratio.test.3way.ixn <-
  anova(primary.analysis$full.model, secondary.analyses$no5$no.3way.ixn.model)
secondary.analyses$no5$p.value.3way.ixn <- secondary.analyses$no5$likelihood.ratio.test.3way.ixn[2, "Pr(>Chisq)"]
secondary.analyses$no5$p.value.3way.ixn

# The arm x visit x year interaction is not significant, so there's no evidence that the divergence
# between the arms from V1 to V2 differs between years. Therefore drop the 3-way interaction and test 
# the three 2-way interactions:
drop1.three2way.ixn <- drop1(secondary.analyses$no5$no.3way.ixn.model, test = "Chisq")
secondary.analyses$no5$p.value.three2way.ixn <- drop1.three2way.ixn[-1, "Pr(>Chi)"]
names(secondary.analyses$no5$p.value.three2way.ixn) <- rownames(drop1.three2way.ixn)[-1]
secondary.analyses$no5$p.value.three2way.ixn

# The year:visit interaction is non-significant marginally, so drop it
secondary.analyses$no5$full.model.noVisitXYearIxn <- 
  update(secondary.analyses$no5$no.3way.ixn.model, ~ . - Visit:Year)
# Check that the remaining terms are significant, and store the p-values
drop1.two2way.ixn <- drop1(secondary.analyses$no5$full.model.noVisitXYearIxn, test = "Chisq")
drop1.two2way.ixn
secondary.analyses$no5$p.value.two2way.ixn <- drop1.two2way.ixn[-1, "Pr(>Chi)"]
names(secondary.analyses$no5$p.value.two2way.ixn) <- rownames(drop1.two2way.ixn)[-1]
secondary.analyses$no5$p.value.two2way.ixn

# ...they are, so this is the final model for the secondary analyses
secondary.analyses$no5$final.model <- secondary.analyses$no5$full.model.noVisitXYearIxn 
secondary.analyses$no5$full.model.noVisitXYearIxn <- NULL

# Therefore we can now estimate, and test for, a single 
# V1-V2 divergence between the arms, which is the arm x visit interaction:
secondary.analyses$no5$estimate.ci.pval <- 
  c(confint(secondary.analyses$no5$final.model, method = "wald")["VisitV2:Trial_ArmCommunity-based", c(3, 1, 2)],
    p.value = as.vector(secondary.analyses$no5$p.value.two2way.ixn["Visit:Trial_Arm"]))
secondary.analyses$no5$effect.estimates <-
  cbind(InteractionOR = paste0(my.format(exp(secondary.analyses$no5$estimate.ci.pval[1:3] * z), ndp = 2), 
                               c(" (", ", ", ")"), collapse = ""),
        `<i>p</i>-value` = p.format(secondary.analyses$no5$estimate.ci.pval["p.value"]))

# Secondary 6

# New table:
# Estimated V1 coverage required to have V2 coverage > threshold, separately for each arm, with 95% CIs.
# Do we need to do this for each year? No, because there were no significant for differences between years
# in arm-specific V1-V2 odds ratios:
secondary.analyses$no5$p.value.3way.ixn

# Having dropped the three way interaction, we asked if the visit x year interaction 
# was significant. We know from secondary analysis 5 above that it isn't:
secondary.analyses$no5$p.value.three2way.ixn
# but having dropped it the remaining 2-way interactions are significant:
secondary.analyses$no5$p.value.two2way.ixn

# Both interactions are significant, so the final secondary #6 model has two interactions
secondary.analyses$no6$final.model <- secondary.analyses$no5$final.model

# Use the reduced model to estimate a single decline in coverage (visit effect) in each arm:

# Create contrasts matrix for extracting V1-V2 log OR by arm
X <- model.matrix(formula(secondary.analyses$no6$final.model, fixed.only = TRUE), 
                  data = cbind(pred.tab[pred.tab$Year == "Y1", ], Vaccinated = 1))
A <- X[X[, "VisitV2"] == 1, ] - X[X[, "VisitV2"] == 0, ]

# There are significant differences between arms (two V1-V2 effects, not one),
# as we know from the interaction test done above:
secondary.analyses$no5$p.value.two2way.ixn["Visit:Trial_Arm"]

# Therefore calculate the two log ORs, ORs, and p-values
secondary.analyses$no6$effect.data <- data.frame(Trial_Arm = gsub("Y1-V2-", "", rownames(A)))
rownames(secondary.analyses$no6$effect.data) <- rownames(A)
secondary.analyses$no6$effect.data$V2V1.logodds.ratio <- 
  (A %*% fixef(secondary.analyses$no6$final.model)$cond)[, 1]
secondary.analyses$no6$effect.data$V2V1.logodds.ratio.se <- 
  sqrt(diag(A %*% vcov(secondary.analyses$no6$final.model)$cond %*% t(A)))
secondary.analyses$no6$effect.data$V2V1.logodds.ratio.p.value <- 
  2 * (1 - pnorm(abs(secondary.analyses$no6$effect.data$V2V1.logodds.ratio)/
                   secondary.analyses$no6$effect.data$V2V1.logodds.ratio.se))
secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.lo <- 
  secondary.analyses$no6$effect.data$V2V1.logodds.ratio - qnorm(0.975) * 
  secondary.analyses$no6$effect.data$V2V1.logodds.ratio.se
secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.hi <- 
  secondary.analyses$no6$effect.data$V2V1.logodds.ratio + qnorm(0.975) * 
  secondary.analyses$no6$effect.data$V2V1.logodds.ratio.se
secondary.analyses$no6$effect.data$V2V1.odds.ratio <-
  exp(secondary.analyses$no6$effect.data$V2V1.logodds.ratio * z)
secondary.analyses$no6$effect.data$V2V1.odds.ratio.ci.lo <-
  exp(secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.lo * z)
secondary.analyses$no6$effect.data$V2V1.odds.ratio.ci.hi <-
  exp(secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.hi * z)


# What would coverage have to be at the start of the year (month 0) in order to have 
# coverage > threshold at month 12 (end of year)
# given these V1-V2 effects?
# If we were to calculate the required coverage at V1 to get a target coverage at V2:
# V1.logodds + V1V2.logOR = V2.logodds, therefore
# V1.logodds = V2.logodds - V1V2.logOR

# However, the time difference between 
# Because the mean visit month for V2 was 12...
mean.visit.month
# ...we can treat V2 and 12-month coverage as equal. 
# However, the average V1 month was 5 months, so to extrapolate back to 
# month 0, we need to "boost" the log odds ratio by multiplying it 
# by 12/diff(mean.visit.month):
secondary.analyses$no6$effect.data$V1.cov.for.V2.eq.thresh <-
  plogis(qlogis(cov.thresh) - 
           secondary.analyses$no6$effect.data$V2V1.logodds.ratio * 12/diff(mean.visit.month) * z)
secondary.analyses$no6$effect.data$V1.cov.for.V2.eq.thresh.ci.lo <-
  plogis(qlogis(cov.thresh) - 
           secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.hi * 12/diff(mean.visit.month) * z)
secondary.analyses$no6$effect.data$V1.cov.for.V2.eq.thresh.ci.hi <-
  plogis(qlogis(cov.thresh) - 
           secondary.analyses$no6$effect.data$V2V1.logodds.ratio.ci.lo * 12/diff(mean.visit.month) * z)
secondary.analyses$no6$effect.data$Trial_Arm <- 
  factor(secondary.analyses$no6$effect.data$Trial_Arm, levels(dat$Trial_Arm))
secondary.analyses$no6$effect.data$x <-
  (as.numeric(secondary.analyses$no6$effect.data$Trial_Arm) - 1.5)/20

# Format the results of secondary analysis 6 into a table
secondary.analyses$no6$Table <- secondary.analyses$no6$effect.data
secondary.analyses$no6$Table$`Visit effect OR (95% CI)` <- 
  paste0(my.format(secondary.analyses$no6$Table$V2V1.odds.ratio, ndp = 2), " (",
         my.format(secondary.analyses$no6$Table$V2V1.odds.ratio.ci.lo, ndp = 2), ", ",
         my.format(secondary.analyses$no6$Table$V2V1.odds.ratio.ci.hi, ndp = 2), ")") 
secondary.analyses$no6$Table$`<i>p</i>-value` <- 
  p.format(secondary.analyses$no6$Table$V2V1.logodds.ratio.p.value)
secondary.analyses$no6$Table$`Start-of-year coverage for end-of-year coverage &ge; threshold` <-
  paste0(my.format(secondary.analyses$no6$Table$V1.cov.for.V2.eq.thresh, ndp = 2), " (",
         my.format(secondary.analyses$no6$Table$V1.cov.for.V2.eq.thresh.ci.lo, ndp = 2), ", ",
         my.format(secondary.analyses$no6$Table$V1.cov.for.V2.eq.thresh.ci.hi, ndp = 2), ")") 
secondary.analyses$no6$Table <-
  secondary.analyses$no6$Table[, c("Trial_Arm", 
                                   "Visit effect OR (95% CI)", 
                                   "<i>p</i>-value", 
                                   "Start-of-year coverage for end-of-year coverage &ge; threshold")]
secondary.analyses$no6$Table <- secondary.analyses$no6$Table[order(secondary.analyses$no6$Table$Trial_Arm), ]

rm(A, X)


# Additional analyses following peer review

# Check the sensitivity of the results to adjusting for land use, noting that 
# the stratified randomisation resulted in there being a higher proportion of 
# urban wards in the team-based arm than in the community-based arm.
table(wards$LandUse, wards$Trial_Arm)
secondary.analyses$no7$final.model <- 
  update(primary.analysis$full.model, ~ . + LandUse)

# Test land use effect
secondary.analyses$no7$likelihood.ratio.test <-
  anova(primary.analysis$full.model, secondary.analyses$no7$final.model)
secondary.analyses$no7$likelihood.ratio.test

# Land use effect analysis p-value
secondary.analyses$no7$landuse.p.value <- 
  secondary.analyses$no7$likelihood.ratio.test[2, "Pr(>Chisq)"]

# Test for spatial autocorrelation between wards in the residuals of the primary
# analysis model by creating a spatial version of the wards data set with ward 
# coordinates caluclated as mean sub-village coordinates, and mean residuals,
# then calculating and testing Moran's I across k = 1 to 10 nearest neighbours.
# The wide range of nearest neighbours tests across a range of scales, 
# from very local (k = 1) to larger scale (k = 10)
wards.spatial <- st_as_sf(wards, coords = c("lon","lat"), crs = 4326)
wards.spatial$mean.resid <- 
  tapply(resid(primary.analysis$full.model),
         primary.analysis$analysis.data$ward, mean, na.rm = TRUE)[rownames(wards)]
moran.stats <-
  t(sapply(1:10, function(k) {
    knn <- knearneigh(st_coordinates(wards.spatial), k = k, longlat = TRUE)
    nb <- knn2nb(knn)
    out <- moran.test(wards.spatial$mean.resid, nb2listw(nb, style = "W"))
    c(k = k, i = out$statistic, p = out$p.value)
  }))
colnames(moran.stats) <- c("kNN", "Moran's <i>I</i>", "<i>p</i>-value")
moran.stats[, 2] <-  round(moran.stats[, 2], 2)
moran.stats[, "<i>p</i>-value"] <-  p.format(moran.stats[, "<i>p</i>-value"])


## Output tables and figures ----

print("Outputting report")

# Initialise table numbering

tabnum <- 1

# Initialise figure numbering

fignum <- 2 # (Fig 1 is not produced by this script)

# Make header and footer file

cat(
  "<html>",
  "<head>",
  "</head>",
  "<body>",
  "<font  face='Times New Roman' size=2 color='#808080'>",
  "<div style='mso-element:header' id=h1>",
  "<p align=right style='margin-top:0pt; margin-bottom:0em;'>",paste0(file.out.name,"x"),"</p>",
  "</div>",
  "<div style='mso-element:footer' id=f1>",
  "<p align=right style='margin-top:0pt; margin-bottom:0em;'>",
  "Page <span style='mso-field-code:\" PAGE\"'></span> of <span style='mso-field-code:\" NUMPAGES \"'></span>",
  "</p>",
  "</div>",
  "</body>",
  "</html>",
  sep="\n",
  file=paste0(header.directory,"/header.html"))

# do title page

cat(
  "<html>",
  "<head>",
  paste0("<link rel=File-List href=\"",strsplit(header.directory, "/")[[1]][2],"/\">"),
  "<style>",
  "@page Section1",
  "     {size:612.0pt 792.0pt;",
  "     mso-page-orientation:landscape;",  # tying to force the word doc to have landscape orientation
  "     margin:2.5cm 2.5cm 2.5cm 2.5cm ;",
  "     mso-header-margin:35.4pt;",
  "     mso-footer-margin:35.4pt;",
  paste0("     mso-header:url(\"",substr(file.out.name,1,nchar(file.out.name)-4),"_files/header.html\") h1;"),
  paste0("     mso-footer:url(\"",substr(file.out.name,1,nchar(file.out.name)-4),"_files/header.html\") f1;"),
  "     }", 
  "div.Section1",
  "     {page:Section1;}",
  "@page Section2",
  "     {size:612.0pt 792.0pt;",
  "     mso-page-orientation:landscape;",
  "     margin: 2.5cm 2.5cm 2.5cm 2.5cm ;",
  "     mso-header-margin:35.4pt;",
  "     mso-footer-margin:35.4pt;",
  paste0("     mso-header:url(\"",substr(file.out.name,1,nchar(file.out.name)-4),"_files/header.html\") h1;"),
  paste0("     mso-footer:url(\"",substr(file.out.name,1,nchar(file.out.name)-4),"_files/header.html\") f1;"),
  "     }", 
  "div.Section2",
  "  {page:Section2;}",
  "</style>",
  "</head>",
  "<body>",
  "<font  face='Times New Roman' size=2>",
  "<div class=Section2>",
  "<h2>",paste(study.title,report.purpose,"Report"),"</h2>\n",
  "<h3>",subtitle,"</h3>\n",
  "<p><em>Report description:</em> <br/>",report.status,"</p>\n",
  "<p><em>Prepared by:</em> <br/>",author,"</p>\n",
  "<p><em>Last run on:</em> <br/>",date(),"</p>\n",
  "<p><em>Created by program:</em> <br/>",prog.file,"</p>\n",
  "<p><em>Created using software:</em> <br/>",R.Version()$version.string,
  " for ",capwords(.Platform$OS.type),
  "with packages",
  paste(need.packages[-length(need.packages)],collapse=", ")," and ",
  need.packages[length(need.packages)],
  "</p>\n",
  "<br style='page-break-after:always'><br/>",      
  sep="\n",file=file.out,append=FALSE)

### Tables ----

#### Tables for main text ----

# Relevant characteristics of wards, by trial arm
Table<-
  demog.tab(
    x.data = wards[, c("district6", "LandUse", 
                       "Total",
                       "n.dogs",
                       "n.vaccinated")],
    cohort = wards$Trial_Arm, total=TRUE,
    max.char=80, n.cohort=TRUE, n.cohort.na = FALSE,
    test.list=FALSE,
    sum.index = 2:4)
#ViewTable.HTML(Table)
output.html.table(Table,
                  paste("Characteristics of wards, overall and by trial arm."),
                  table.number=tabnum,file=file.out,
                  fullWidth=TRUE,new.page.after=FALSE,font.size=9)
tabnum<-tabnum+1; rm(Table)


# Coverage estimates and OR estimates from the full model
tab.heads <-
  cbindTable(colnames(primary.analysis$effect.estimates)[1],
             rbindTable("Coverage (95% CI)", 
                        cbindTable(colnames(primary.analysis$effect.estimates)[2], colnames(primary.analysis$effect.estimates)[3])),
             colnames(primary.analysis$effect.estimates)[4],
             colnames(primary.analysis$effect.estimates)[5],
             "<i>p</i>-value<br/>(primary analysis)", 
             "<i>p</i>-value<br/>(secondary analysis 1)")
Table <- 
  list(tab.heads = tab.heads,
       tab.rows =
         rbindTable(
           cbindTable(cbindMatrix(as.matrix(primary.analysis$effect.estimates)), 
                      p.format(primary.analysis$intervention.p.value), 
                      p.format(secondary.analyses$no1$intervention.p.value)),
           cbindTable(cbindMatrix(primary.analysis$effect.estimates.allyears), "", "")))

#ViewTable.HTML(Table)

output.html.table(Table,
                  paste("Primary analysis and secondary analysis 1. Coverage and intervention odds ratio estimates (95% CI) at each survey",
                        "estimated from the GLMM reported in the previous table. <i>p</i>-values at each survey are from likelihood ratio tests of the",
                        "primary analysis null hypothesis of no intervention effect at any of the six surveys and secondary",
                        "analysis 1 null hypothesis of no intervention effect at any of the three visit 2 surveys.",
                        "Mean coverage (95% CI), odds ratio (95% CI) and <i>p</i>-value across all years were estimated from",
                        secondary.analyses$no2$nsim, "parametric bootstrap samples.",
                        "Population: all dogs with primary outcome data."),
                  table.number=tabnum,file=file.out,
                  fullWidth=TRUE,new.page=TRUE,
                  new.page.after=FALSE,font.size=9)
tabnum<-tabnum+1; rm(Table, tab.heads)

# Secondary analysis 2

# Table of P(coverage < threshold) by month and arm
Table <- 
  list(tab.heads = 
         cbindTable("Trial month", 
                    rbindTable(colnames(secondary.analyses$no2$prop.estimates.overall)[2], 
                               cbindMatrix(colnames(secondary.analyses$no2$cov.lt40.table)))),
       tab.rows = cbindMatrix(secondary.analyses$no2$cov.lt40.table, use.rownames = TRUE))

Table$tab.rows <- 
  rbindTable(Table$tab.rows, 
             cbindMatrix(c("All months", secondary.analyses$no2$prop.estimates.overall[, 2])))

#ViewTable.HTML(Table)

output.html.table(Table,
                  paste("Secondary analysis 2. Estimated probability (95% CI) of village-level coverage below",
                        paste0(100 * cov.thresh, "% by month and averaged over the trial year, for each trial arm."),
                        "The distribution of coverage among villages was estimated from",
                        secondary.analyses$no2$nsim, "parametric bootstrap samples. ",
                        "Population: all dogs with primary outcome data."),
                  table.number=tabnum,file=file.out,
                  fullWidth=TRUE,new.page.after=FALSE,font.size=9)
tabnum<-tabnum+1; rm(Table)

# Secondary analyses 3 and 4
Table <- 
  list(tab.heads = cbindMatrix(rbind(c("Analysis", colnames(secondary.analyses$no3no4$Estimates)))),
       tab.rows = cbindMatrix(secondary.analyses$no3no4$Estimates, use.rownames = TRUE))
#ViewTable.HTML(Table)
output.html.table(Table,
                  paste("Secondary analyses 3 and 4. Estimated median (95% CI) Community-based:Team-based temporal and spatial inter-village variance ratio.",
                        "Population: all dogs with primary outcome data."),
                  table.number="for main text",file=file.out,
                  fullWidth=TRUE,new.page.after=FALSE,font.size=9)
#tabnum<-tabnum+1; 
rm(Table)



# Secondary analysis 5
Table <- 
  list(tab.heads = cbindMatrix(rbind(colnames(secondary.analyses$no5$effect.estimates))),
       tab.rows = cbindMatrix(secondary.analyses$no5$effect.estimates))
#ViewTable.HTML(Table)
output.html.table(Table,
                  paste0("Secondary analysis 5. Estimate (95% CI) V2:V1 ratio of Community-based:Team-based odds ratios. ", 
                         "This ratio represents the interaction between trial arm and visit; the <i>p</i>-value tests ",
                         "the null hypothesis of equal intervention effect between visit 1 and visit 2. ",
                         "This interaction effect did not differ between the three years ",
                         "(trial arm &times; visit &times; year interaction <i>p</i>-value = ",
                         p.format(secondary.analyses$no5$p.value.3way.ixn), "). ",
                         "Population: all dogs with primary outcome data."),
                  table.number="for main text",file=file.out,
                  fullWidth=TRUE,new.page.after=FALSE,font.size=9)
#tabnum<-tabnum+1; 
rm(Table)


# Secondary analysis 6

# Table showing V2:V1 odds ratios and coverage required at V1 to keep coverage over threshold at V2
Table <- 
  list(tab.heads = cbindMatrix(rbind(colnames(secondary.analyses$no6$Table))),
       tab.rows = cbindMatrix(as.matrix(secondary.analyses$no6$Table)))
#ViewTable.HTML(Table)
output.html.table(Table,
                  paste0("Secondary analysis 6. Estimate (95% CI) of V2:V1 odds ratios and of minimum coverage ",
                         "at visit 1 required for coverage at visit 2 to exceed ", 
                         100 * cov.thresh, "%, by trial arm. ", 
                         "V2:V1 odds ratios differed significantly between trial arms (interaction P ",
                         p.format(secondary.analyses$no5$p.value.two2way.ixn["Visit:Trial_Arm"]),
                         "). Population: all dogs with primary outcome data."),
                  table.number=tabnum,file=file.out,
                  fullWidth=TRUE,new.page.after=TRUE,font.size=9)
tabnum<-tabnum+1; rm(Table)

#### Supplementary tables ----

# Re-initialise table numbering
tabnum <- 1

# Estimates from full model
Table <- 
  tab_model(primary.analysis$null.model, primary.analysis$full.model, 
            transform = NULL, emph.p = FALSE, show.se = TRUE, show.icc = FALSE,
            ci.hyphen = ", ",
            dv.labels = c("Primary analysis null hypothesis model", "Primary analysis alternative hypothesis model"),
            title = paste0("<b>Table S", tabnum, 
                           "</b>. Estimates of fixed effects (log odds and log odds ratios) and ",
                           "random effects (variances) from the GLMMs fitted for the primary analysis. ",
                           "Number of observations, number of each random effect level, and marginal and conditional ",
                           "R<sup>2</sup> are also presented. ", 
                           "The null hypothesis of no intervention effect was rejected (",
                           "&chi;<sup>2</sup>(", diff(primary.analysis$likelihood.ratio.test[, "Df"]),
                           "), <i>p</i> ", gsub("<", "&lt; ", p.format(primary.analysis$intervention.p.value)), "). ",
                           "Population: all dogs with primary outcome data."),
            CSS = list(css.table = "font-size:9", 
                       css.caption = "+font-weight: normal"), 
            show.intercept = TRUE)
# reduce table cell padding
for(page.element in c("page.style", "page.content", "page.complete", "knitr")) {
  Table[[page.element]] <- gsub("padding:0.2cm", "padding:0.05cm", Table[[page.element]])
  Table[[page.element]] <- gsub("padding-bottom:0.1cm", "padding-bottom:0.01cm", Table[[page.element]])
  Table[[page.element]] <- gsub("padding-top:0.1cm", "padding-top:0.01cm", Table[[page.element]])
}
cat(Table$knitr, file=file.out, sep="\n", append=T)
cat('<p>&nbsp;</p>', file = file.out, sep="\n", append = TRUE)
cat('<p>&nbsp;</p>', file = file.out, sep="\n", append = TRUE)
tabnum<-tabnum+1; rm(Table)


# Compare full primary analysis model with same model adjusted for land use category
Table <- 
  tab_model(primary.analysis$full.model, secondary.analyses$no7$final.model, 
            transform = NULL, emph.p = FALSE, show.se = TRUE, show.icc = FALSE,
            ci.hyphen = ", ",
            dv.labels = c("Primary analysis full model", "Primary analysis full model adjusted for land use"),
            title = paste0("<b>Table S", tabnum, 
                           "</b>. Estimates of fixed effects (log odds and log odds ratios) and ",
                           "random effects (variances) from the primary analysis GLMM and the same model adjusted for land use category. ",
                           "Number of observations, number of each random effect level, and marginal and conditional ",
                           "R<sup>2</sup> are also presented. ", 
                           "Adding land use category to the primary analysis model did not significantly improve fit (",
                           "&chi;<sup>2</sup>(", diff(secondary.analyses$no7$likelihood.ratio.test[, "Df"]),
                           "), <i>p</i> = ", p.format(secondary.analyses$no7$landuse.p.value), "). ",
                           "Population: all dogs with primary outcome data."),
            CSS = list(css.table = "font-size:9", 
                       css.caption = "+font-weight: normal"), 
            show.intercept = TRUE)
# reduce table cell padding
for(page.element in c("page.style", "page.content", "page.complete", "knitr")) {
  Table[[page.element]] <- gsub("padding:0.2cm", "padding:0.05cm", Table[[page.element]])
  Table[[page.element]] <- gsub("padding-bottom:0.1cm", "padding-bottom:0.01cm", Table[[page.element]])
  Table[[page.element]] <- gsub("padding-top:0.1cm", "padding-top:0.01cm", Table[[page.element]])
}
cat(Table$knitr, file=file.out, sep="\n", append=T)
cat('<p>&nbsp;</p>', file = file.out, sep="\n", append = TRUE)
tabnum<-tabnum+1; rm(Table)

# Relevant characteristics of dogs, by trial arm
t1.vars <- c("age.years", "Sex", "YearVisit",
             "Vaccinated.cat.na", "Vaccinated2.cat")
names(dat)[!names(dat) %in% t1.vars]
Table<-
  demog.tab(
    x.data = dat[, t1.vars],
    cohort = dat$Trial_Arm, total=TRUE,
    max.char=80, n.cohort=TRUE, test.list=FALSE)
#ViewTable.HTML(Table)   
output.html.table(Table,
                  paste("Characteristics of surveyed dogs, overall and by trial arm. Population: all dogs."),
                  table.number=paste0("S", tabnum),file=file.out,
                  fullWidth=TRUE, new.page.after=TRUE, font.size=9)
tabnum<-tabnum+1; rm(Table)

# Moran's I test for spatial autocorrelation
tab.heads <-
  cbindMatrix(colnames(moran.stats))
Table <- 
  list(tab.heads = tab.heads,
       tab.rows = cbindMatrix(moran.stats))
#ViewTable.HTML(Table)
output.html.table(Table,
                  paste("Moran's test for spatial autocorrelation in the residuals of the primary analysis full model,",
                        "applied at a range of spatial scales defined by the number of nearest neighbours (kNN)."),
                  table.number=paste0("S", tabnum),file=file.out,
                  fullWidth=TRUE, new.page.after=TRUE, font.size=9)
tabnum<-tabnum+1; rm(Table)


### Figures ----

# Plot of the number of vaccinations and the number of dogs surveyed per day, by arm
height <- 1300
width <- 1600
tiff(fig.name(fignum, figures.directory = figures.directory, fmt = ".tiff"), 
    height = height, width = width, units = "px",
    res = 300)
old.par <- par(mar = c(5.1, 4.1, 2.1, 4.1), mfrow = c(1, 1), mgp = c(1.8, 0.5, 0))
vax$day.vaccinated.fac <- factor(vax$day.vaccinated, min(vax$day.vaccinated):max(vax$day.vaccinated))
vax.tab <- table(vax$Trial_Arm, vax$day.vaccinated.fac)
vax.ymax <- 3500
vax.plot.alpha <- 0.8
vax.plot.lwd <- 1.5
vax.cex <- 0.7
plot(colnames(vax.tab), vax.tab["Community-based", ], type = "l", 
     col = alpha(arm.colours["Community-based"], vax.plot.alpha), lwd = vax.plot.lwd,
     ylab = "N vaccinations (lower series)", xlab = "Month",
     axes = FALSE,
     ylim = c(0, vax.ymax),
     xlim = range(c(dat$day.added, vax$day.vaccinated)), cex.lab = vax.cex)
lines(colnames(vax.tab), vax.tab["Team-based", ], col = alpha(arm.colours["Team-based"], vax.plot.alpha), 
      lwd = vax.plot.lwd)
legend("right", legend = names(arm.colours), lwd = vax.plot.lwd, 
       col = alpha(arm.colours, vax.plot.alpha), bty = "n", cex = vax.cex)
vax.plot.x.lab <- pretty(as.numeric(colnames(vax.tab)) * 12 /365, n = 10, bty = "n")
vax.plot.x.at <- vax.plot.x.lab * 365 / 12
axis(1, at = vax.plot.x.at, labels = vax.plot.x.lab, cex.axis = vax.cex)
axis(2, cex.axis = vax.cex)
axis(3, at = mean.visit.dates - min(vax$date.vaccinated) + 1,
     labels = rep(c(expression("V"[1]), expression("V"[2])), 3), cex.axis = vax.cex)
axis(4, at = pretty(c(0, vax.ymax)), labels = (vax.ymax - pretty(c(0, vax.ymax)))/10,
     cex.axis = vax.cex)
box()
vax$month.vaccinated <- substr(vax$date.vaccinated, 1, 7)
table(vax$month.vaccinated, vax$Trial_Arm)

dat$day.added <- as.numeric(dat$date.added - min(vax$date.vaccinated)) + 1
dat$day.added.fac <- factor(dat$day.added, min(dat$day.added):max(dat$day.added))

survey.tab <- table(dat$Trial_Arm, dat$day.added.fac)
lines(colnames(survey.tab), vax.ymax - 10 * survey.tab["Community-based", ] - 8, 
      col = alpha(arm.colours["Community-based"], vax.plot.alpha/1.5), lwd = vax.plot.lwd)
lines(colnames(survey.tab), vax.ymax - 10 * survey.tab["Team-based", ], 
      col = alpha(arm.colours["Team-based"], vax.plot.alpha/1.5), lwd = vax.plot.lwd)

mtext("N dogs surveyed (upper series)", side = 4, line = 2, cex = vax.cex)
par(old.par)
dev.off()

cat(
  "<p align='justify'>",     
  paste0("<img src='",paste(fig.name(fignum, figures.directory = strsplit(figures.directory, "/")[[1]][2], fmt = ".tiff")),
         "' height=", height/5, " width=", width/5, " >\n"),
  paste0(
    "<br><br><b>Figure ",fignum,
    ".</b> Number of rabies vaccinations performed (lower series; total = ", sum(vax.tab),
    ") and number of dogs surveyed (upper series; total = ", sum(survey.tab), 
    ") per day, by trial arm.</p><br>"),
  file=file.out,sep="\n",append=T)
fignum<-fignum+1


# Figure of coverage estimates by arm and survey time point
height <- 1350
width <- 1575
tiff(fig.name(fignum, figures.directory = figures.directory, fmt = ".tiff"), 
    height = height, width = width, units = "px",
    res = 300)
old.par<-par(mfrow = c(1, 1), mgp = c(1.8, 0.6, 0))
estimate.ci.plot(x.term = "x", y.term = "Coverage", 
                 y.ci.lo = "Coverage.ci.lo", y.ci.hi = "Coverage.ci.hi", 
                 plotdata = pred.tab, xlab = "Visit", ylab = "Coverage \u00B1 95% CI", 
                 pch = 20 + as.numeric(pred.tab$Trial_Arm), 
                 ylim = 0:1,
                 bg = arm.colours, main = "",
                 x.at = unique(round(pred.tab$x)), 
                 x.labels = sapply(sapply(add_brackets(levels(pred.tab$YearVisit)), 
                                          function(x) substitute(parse(text = a), list(a = x))), eval),
                 cex = 1.2)
lines(Coverage ~ x, data = pred.tab[pred.tab$Trial_Arm == "Team-based", ], type = "b", cex = 0)
lines(Coverage ~ x, data = pred.tab[pred.tab$Trial_Arm == "Community-based", ], type = "b", cex = 0)
legend("topleft", legend = levels(dat$Trial_Arm), pch = 20 + 1:nlevels(dat$Trial_Arm),
       pt.bg = arm.colours, pt.cex = 1.5, bty = "n")
abline(h = cov.thresh, lty = 2)
par(old.par)
dev.off()
cat(
  "<p align='justify'>",     
  paste0("<img src='",paste(fig.name(fignum, figures.directory = strsplit(figures.directory, "/")[[1]][2], fmt = ".tiff")),
         "' height=", height/5, " width=", width/5, " >\n"),
  paste0(
    "<br><br><b>Figure ",fignum,
    ".</b> Estimated coverage &plusmn; 95% confidence limits at each survey time point, by trial arm. ",
    "Population: all dogs with primary outcome data.</p><br>"),
  file=file.out,sep="\n",append=T)
fignum<-fignum+1


# Violin plot showing probability distribution of coverage by trial month for each arm
coverage.dist$YearMonthArm <- 
  factor(paste(coverage.dist$YearMonth, coverage.dist$Trial_Arm),
         paste(rep(levels(coverage.dist$YearMonth), nlevels(coverage.dist$Trial_Arm)), 
               rep(levels(coverage.dist$Trial_Arm), each = nlevels(coverage.dist$YearMonth))))
height <- 1000
width <- 2000
tiff(fig.name(fignum, figures.directory = figures.directory, fmt = ".tiff"), 
    height = height, width = width, units = "px",
    res = 300)
old.par<-par(mfrow = c(1, 1), mar = c(3.1, 3.1, 1.1, 1.1), mgp = c(1.8, 0.6, 0))
vioplot(Coverage ~ YearMonthArm, data = coverage.dist, axes = FALSE,
        ylim = c(0, 1.05), 
        xlim = c(1.3, nlevels(coverage.dist$YearMonthArm) - 0.3),
        col = rep(arm.colours, each = nlevels(coverage.dist$YearMonth)),
        rectCol = "grey",
        lineCol = "grey",
        cex.axis = 0.001, xlab = "Trial Month") 
axis(2, cex.axis = 0.7)
axis(1, at = 1:nlevels(coverage.dist$YearMonthArm), 
     labels = rep(levels(coverage.dist$YearMonth), nlevels(coverage.dist$Trial_Arm)),
     cex.axis = 0.7)
abline(h = cov.thresh, lty = 2)
abline(v = 12.5)
box()
text(x = c(6, 18), y = rep(1.05, 2), labels = names(arm.colours))
par(old.par)
dev.off()
cat(
  "<p align='justify'>",     
  paste0("<img src='",paste(fig.name(fignum, figures.directory = strsplit(figures.directory, "/")[[1]][2], fmt = ".tiff")),
         "' height=", height/5, " width=", width/5, " >\n"),
  paste0(
    "<br><br><b>Figure ",fignum,
    ".</b> Estimated probability distribution of coverage by trial month for each arm. ",
    "The coverage target of ", 100 * cov.thresh, "% is shown by a dashed line. ",
    "White circles indicate medians and quartiles and range are shown by box-and-whisker plots.</p><br>"),
  file=file.out,sep="\n",append=T)
fignum<-fignum+1


# Figure of odds ratios estimates by survey time point
height <- 1000
width <- 1700
tiff(fig.name(fignum, figures.directory = figures.directory, fmt = ".tiff"), 
    height = height, width = width, units = "px",
    res = 300)
old.par <- par(mfrow = c(1, 1), mar = c(3.1, 3.1, 1.1, 1.1), mgp = c(1.8, 0.6, 0))
estimate.ci.plot(x.term = "x", y.term = "InterventionOR", 
                 y.ci.lo = "InterventionOR.ci.lo", y.ci.hi = "InterventionOR.ci.hi", 
                 plotdata = effect.data, xlab = "", 
                 ylab = "Odds ratio \u00B1 95% CI",
                 main = "",
                 ylim = c(range(c(0.9, effect.data$InterventionOR.ci.lo, effect.data$InterventionOR.ci.hi))),
                 pch = 21, 
                 bg = "white",
                 x.at = effect.data$x,
                 x.labels = sapply(sapply(add_brackets(levels(effect.data$YearVisit)), 
                                          function(x) substitute(parse(text = a), list(a = x))), eval),
                 cex = 1.5)
lines(InterventionOR ~ x, data = effect.data, type = "b", cex = 0)
abline(h = 1, lty = 2)
if(names(vaccination.definition) == "main") axis(side = 2, at = 1)
par(old.par)
dev.off()
cat(
  "<p align='justify'>",     
  paste0("<img src='",paste(fig.name(fignum, figures.directory = strsplit(figures.directory, "/")[[1]][2], fmt = ".tiff")),
         "' height=", height/5, " width=", width/5, " >\n"),
  paste0(
    "<br><br><b>Figure ",fignum,
    ".</b> Estimated intervention effect (Community-based:Team-based) odds ratios &plusmn; 95% confidence limits at each survey time point. ",
    "Population: all dogs with primary outcome data.</p><br>"),
  file=file.out,sep="\n",append=T)
fignum<-fignum+1

#### Supplementary figures ----

# Re-initialise figure numbering
fignum <- 1

# Figures of raw coverage estimates by survey time point and ward, one per arm
height <- 1000
width <- 1600
for(arm in names(arm.colours)) {
  print(arm.colours[arm])
  tiff(fig.name(paste0("S", fignum), figures.directory = figures.directory, fmt = ".tiff"), 
      height = height, width = width, units = "px",
      res = 300)
  raw.cov.plot <-
    ggplot(data = wards.rnd.reduced[wards.rnd.reduced$Trial_Arm == arm, ], 
           aes(x = YearVisit, y = Coverage, group = ward, shape = Trial_Arm)) +
    geom_point(aes(size = sqrt(I(N)/20)), color = arm.colours[arm], show.legend = FALSE) +
    geom_line(linewidth = 0.2) + 
    facet_wrap(~ district + ward.short) +
    xlab("") +
    guides(col = "none") + 
    theme_bw(base_size = 4) +
    theme(axis.text.x = element_text(angle = 90)) 
  print(raw.cov.plot)
  dev.off()
  cat(
    "<p align='justify'>",     
    paste0("<img src='",paste(fig.name(paste0("S", fignum), figures.directory = strsplit(figures.directory, "/")[[1]][2], fmt = ".tiff")),
           "' height=", height/6, " width=", width/6, " >\n"),
    paste0(
      "<br><br><b>Figure S",fignum,
      ".</b> Coverage by ward and survey time point in the ", arm,
      " arm.</p><br>"),
    file=file.out,sep="\n",append=T)
  fignum<-fignum+1
  rm(raw.cov.plot, arm)
}


# Close html file
cat("</div>","</body>","</html>",file=file.out,sep="\n",append=T)

file.rename(file.out, paste0(file.out, ".html"))
system(paste("open", paste0(file.out, ".html")))
Sys.sleep(10)
file.rename(paste0(file.out, ".html"), file.out)
system(paste("open", file.out))

# Report run time
finish.time <- Sys.time() 
finish.time - start.time
