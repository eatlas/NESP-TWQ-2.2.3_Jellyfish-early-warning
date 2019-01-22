# Irukandji analysis
# 28/11/2017
# Anthony J. Richardson
# This work is licensed under a Creative Commons Attribution 4.0 International License. 
# (https://creativecommons.org/licenses/by/4.0/)
#
# Treat each row as a sting - ignore NUMBER_STINGS as a variable because it is the # stings on a person
# Approach:
# 1. Where do we have good data through time? From 1985 (currently missing 150 obs)
# Rectify wind data N:S, E:W and NS=EW. Do cumulative sums for each day 1-14 days before
# Make a dataframe for analysis
# The unit of analysis is # stings for each Y, M, D, Region (so different Lat and Lon in same Region treated the same)
# 2. Calculate number of stings each day
# 2. Put in 0 stings for every day (without data) from 1985
# 3. Do model for each region to start with: FNQ, NQ, CQ (and separate Mainland and Islands/Reefs)
# Choose best model on minimising AIC
# 4. Model #Stings(0, 1, 2, ...) ~ Month + Wind (1-14 days: NS-EW) + Tide (Incoming, Outgoing) + SST + Curator (Gerswhin, Fenner, Southcott)
# Or model as 0, 1 response for each day
# Use negative binomial because lots of zeros

library(MASS) # If use negative binomial
library(splines)
library(lubridate)
library(ggplot2)
library(gridExtra) # For multiple plots on a page
library(maps)
library(mapdata) # For hi-res Australia map
library(tidyverse)
library(effects)
library(visreg)
library(Amelia) # Visually plotting MD
library(GGally) # Good pairs plot
library(circular)
library(clifro) # Plotting windroses
rm(list = ls())
source("~ric325/Dropbox/Documents/IMOS/ClaireDavies/Harm.r")

#################################
### Create function for mode ###
getmode <- function(v) {
  if(max(v) == 0)
    v <- 15 # 15h00
  v <- v[v != 0] # remove zeros
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#################################
SwitchDirn <- function(theta) {
for(i in 1:length(theta)){
  if(theta[i] < 90){
    theta[i] <- 90 - theta[i] + 0}
  else if(theta[i] < 180 & theta[i] >= 90){
    theta[i] <- 180 - theta[i] + 270}
  else if(theta[i] < 270 & theta[i] >= 180){
    theta[i] <- 270 - theta[i] + 180}
  else if(theta[i] < 360 & theta[i] >= 270){
    theta[i] <- 360 - theta[i] + 90}
  }
SwitchDirn <- theta
}

################################################
PlotRose <- function(Speed, Direction, Title){
  SpeedCuts <- c(3, 6, 9, 12, 19)
  clifro::windrose(Speed, Direction, Title,
                   speed_cuts = SpeedCuts,
                   legend_title = "Wind Speed\n(m/s)",
                   legend.title.align = .5,
                   ggtheme = "bw",
                   col_pal = "GnBu")
}
#######################################
AddTidalRange <- function(dat, i) {
  Files_Tide <- list.files(path = "./TideData")
  dat_Tide <- read.table(paste("./TideData/", Files_Tide[i], sep = ""), sep = " ", header = FALSE, comment.char = '#') # Remove # 21 lines to skip (comment character)
  dat_Tide <- dat_Tide %>%
    rename(DaysSince1990 = V1, Height = V2) %>%
    mutate(Date = as.POSIXct(60*60*24*DaysSince1990, origin = '1990-01-01 00:00', tz = 'Australia/Hobart')) %>% #, tz = 'GMT') # Australia
    mutate(Date = ymd_hms(Date)) # 1 hr out of phase (Hobart time?)???
  dat_Tide$Time <- strftime(round_date(dat_Tide$Date, unit = "hour"), format="%H:%M", tz = 'UTC')
  dat_Tide$Date <- date(dat_Tide$Date) # Converts to date format
  tempory1 <- with(dat_Tide, aggregate(Height, by = list(Date), min))
  tempory2 <- with(dat_Tide, aggregate(Height, by = list(Date), max))
  tempory1$TidalRange <- tempory2$x - tempory1$x
  tempory1 <- tempory1 %>% mutate(TidalRange = tempory2$x - x) %>%
    rename(Date = Group.1) %>%
    dplyr::select(Date, TidalRange)
  dat2 <- left_join(x = dat, y = tempory1, by = c("Date"))
  return(dat2)
}
#######################################

AddTides <- function(dat, i) {
  Files_Tide <- list.files(path = "./TideData")
  dat_Tide <- read.table(paste("./TideData/", Files_Tide[i], sep = ""), sep = " ", header = FALSE, comment.char = '#') # Remove # 21 lines to skip (comment character)
  dat_Tide <- dat_Tide %>%
    rename(DaysSince1990 = V1, Height = V2) %>%
    mutate(Date = as.POSIXct(60*60*24*DaysSince1990, origin = '1990-01-01 00:00', tz = 'Australia/Hobart')) %>% #, tz = 'GMT') # Australia
    mutate(Date = ymd_hms(Date)) # 1 hr out of phase (Hobart time?)???
  
  # Separate into Incoming and Outgoing tides
  dat_Tide$Incoming[2:length(dat_Tide$Height)] <- dat_Tide$Height[2:length(dat_Tide$Height)] > dat_Tide$Height[1:(length(dat_Tide$Height)-1)]
  dat_Tide$Incoming <- as.factor(dat_Tide$Incoming)
  
  # Calculate times...
  dat_Tide$Time <- strftime(round_date(dat_Tide$Date, unit = "hour"), format="%H:%M", tz = 'UTC')
  dat$Time2 <- dat$Time
  dat_Tide$Date <- date(dat_Tide$Date)
  dat$Time <- paste(as.character(round(dat$Time)), "00", sep = ":")
  dat$Time[nchar(dat$Time) == 4] <- paste("0", dat$Time[nchar(dat$Time) == 4], sep = "") # For 7:00, 8:00, 9:00 am need to add 0 in front so 07:00 so can left_join with tide time
  dat2 <- left_join(x = dat, y = dat_Tide, by = c("Date", "Time"))
  return(dat2)
}
##################################################################

MakeData <- function(dat, Reg, WindFile) {
  ### Choose the region
  JellyDat <- dat %>%
    filter(Region2 == Reg) %>%
    group_by(Year, Month, Day) %>%
    summarise(n = n(), MnTime = getmode(Time), DOY = mean(DOY), DOYSeason = mean(DOYSeason)) %>%
    rename(Month = Month)
  
  ### Add zeros ###
  # Create dataframe containing all days from 1985 to end of 2016
  AllDays <- data.frame(Date = seq(as.Date("1985/1/1"), by = "day", length.out = 365.25*(2016-1985+1)))
  AllDays <- AllDays %>%
    mutate(Year = as.factor(year(Date)), Month = as.factor(month(Date)), Day = day(Date), Stings = 0)
  JellyDat <- left_join(x = AllDays, y = JellyDat, by = c("Year", "Month", "Day")) # Warning because no stings in August
  JellyDat <- JellyDat %>% 
    mutate(Stings = pmax(Stings, n, na.rm = TRUE)) %>%
    dplyr::select(-n) %>%
    mutate(Time = MnTime)
  # JellyDat$Time[is.na(JellyDat$MnTime)] <- 15 # Replace missing times with most frequent time
  # For MD for times generate random ones between min and max of times
  Times <- JellyDat$Time[JellyDat$MnTime != 0]
  # NumOfZeros <- length(JellyDat$Time[JellyDat$Time == 0])
  # JellyDat$Time[JellyDat$Time == 0] <- runif(min = min(Times), max = max(Times), n = NumOfZeros)
  NumOfZeros <- length(JellyDat$Time[is.na(JellyDat$Time)])
  JellyDat$Time[is.na(JellyDat$Time)] <- runif(min = min(Times, na.rm = TRUE), max = max(Times, na.rm = TRUE), n = NumOfZeros)
  rm(AllDays)
  
  ### Add winds ###
  WindDat <- read.table(paste("./WindData/", WindFile, sep = ""), sep = ",", header = FALSE, comment.char = '#') # Cairns
  WindDat <- WindDat %>% # Calculate speed and dirn
    rename(Date = V1, NS = V2, EW = V3) %>%
    mutate(Date = ymd(Date)) %>%
    mutate(Speed = sqrt(NS^2 + EW^2)) %>%
    mutate(Dirn = (atan2(NS, EW) * (180 / pi))) %>% 
    mutate(S1 = -1, S2 = -1, S3 = -1, S4 = -1, S5 = -1, S6 = -1, S7 = -1, S8 = -1, S9 = -1, S10 = -1, S11 = -1, S12 = -1, S13 = -1, S14 = -1, S15 = -1) %>% # Add columns for Speed
    mutate(D1 = -1, D2 = -1, D3 = -1, D4 = -1, D5 = -1, D6 = -1, D7 = -1, D8 = -1, D9 = -1, D10 = -1, D11 = -1, D12 = -1, D13 = -1, D14 = -1, D15 = -1) # Add columns for Dirn
  
  WindDat$Dirn[WindDat$Dirn < 0] <- WindDat$Dirn[WindDat$Dirn < 0] + 360 # This is so the atan2 angles are correct
  
  # Calculate mean winds over multiple days
  NumOfDays <- 14
  SpeedCols <- 6 # Start of Speed columns
  DirnCols <- 21 # Start of Dirn columns
  conv <- 2*pi/360 # degrees to radians for mean of angles (circular data)
  Weights <- 1/(1:15) # Weights <- 1/(1:15)^2 # Weights for each wind and dirn through time
  # Weights <- rep(1/15, 15)
  # Weights <- c(1/3, 1/3, 1/3, rep(0, 12))
  SumsW <- cumsum(Weights)
  for(j in 15:nrow(WindDat)){
    for(i in 0:NumOfDays){
      WindDat[j, i+SpeedCols] <- sum(WindDat[j:(j-i), 4] * Weights[1:(i+1)]) / SumsW[i+1]
      Dirns <- WindDat[j:(j-i), 5]
      WindDat[j, i+DirnCols] <- weighted.mean.circular(circular(Dirns*conv), Weights[1:(i+1)])/conv
      # WindDat_Cairns[j, i+DirnCols] <- mean.circular(circular(WindDat_Cairns[j:(j-i), 5]*conv))/conv
    }
  }
  WindDat[, 21:35][WindDat[, 21:35] < 0] <- WindDat[, 21:35][WindDat[, 21:35] < 0] + 360 # Add 360o to negative values to make them positive
  # write.csv(head(WindDat, 20), file = "WindDat.csv")
  
  # Check winds
  clifro::windrose(WindDat$Speed, SwitchDirn(WindDat$Dirn), "Cairns") # Winds OK!
  hist(WindDat$Dirn)
  
  # Join JellyDat with WindSpeed
  JellyDat <- left_join(x = JellyDat, y = WindDat, by = "Date") # Warning because no stings in August
  JellyDat <- JellyDat %>% 
    mutate(Month = as.factor(Month)) %>%
    mutate(Month = factor(Month, levels(Month)[c(1,5,6,7,8,9,10,11,12,2,3,4)])) %>%
    mutate(Yr = case_when(Month == "1" | Month == "2" | Month == "3" | Month == "4" | Month == "5" | Month == "6" ~ as.numeric(as.character(Year)) - 1,
                        TRUE ~ as.numeric(as.character(Year)))) %>% # Make stinger season
    mutate(Yr = as.factor(Yr)) %>%
    mutate(Mn = Month) %>%
    mutate(Mn = factor(Mn, levels(Month)[c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)])) %>%
    mutate(Angle1 = case_when( # Split Angles
      D1 >= 0 & D1 < 90 ~ "NE",
      D1 >= 90 & D1 < 180 ~ "NW",
      D1 >= 180 & D1 < 270 ~ "SW",
      D1 >= 270 & D1 < 360 ~ "SE")) %>%
    mutate(Angle1 = as.factor(Angle1)) %>%
    mutate(Angle1 = factor(Angle1, levels(Angle1)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle2 = case_when( # Split Angles
      D2 >= 0 & D2 < 90 ~ "NE",
      D2 >= 90 & D2 < 180 ~ "NW",
      D2 >= 180 & D2 < 270 ~ "SW",
      D2 >= 270 & D2 < 360 ~ "SE")) %>%
    mutate(Angle2 = as.factor(Angle2)) %>%
    mutate(Angle2 = factor(Angle2, levels(Angle2)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle3 = case_when( # Split Angles
      D3 >= 0 & D3 < 90 ~ "NE",
      D3 >= 90 & D3 < 180 ~ "NW",
      D3 >= 180 & D3 < 270 ~ "SW",
      D3 >= 270 & D3 < 360 ~ "SE")) %>%
    mutate(Angle3 = as.factor(Angle3)) %>%
    mutate(Angle3 = factor(Angle3, levels(Angle3)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle4 = case_when( # Split Angles
      D4 >= 0 & D4 < 90 ~ "NE",
      D4 >= 90 & D4 < 180 ~ "NW",
      D4 >= 180 & D4 < 270 ~ "SW",
      D4 >= 270 & D4 < 360 ~ "SE")) %>%
    mutate(Angle4 = as.factor(Angle4)) %>%
    mutate(Angle4 = factor(Angle4, levels(Angle4)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle5 = case_when( # Split Angles
      D5 >= 0 & D5 < 90 ~ "NE",
      D5 >= 90 & D5 < 180 ~ "NW",
      D5 >= 180 & D5 < 270 ~ "SW",
      D5 >= 270 & D5 < 360 ~ "SE")) %>%
    mutate(Angle5 = as.factor(Angle5)) %>%
    mutate(Angle5 = factor(Angle5, levels(Angle5)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle6 = case_when( # Split Angles
      D6 >= 0 & D6 < 90 ~ "NE", 
      D6 >= 90 & D6 < 180 ~ "NW",
      D6 >= 180 & D6 < 270 ~ "SW",
      D6 >= 270 & D6 < 360 ~ "SE")) %>%
    mutate(Angle6 = as.factor(Angle6)) %>%
    mutate(Angle6 = factor(Angle6, levels(Angle6)[c(1, 2, 4, 3)])) %>% 

    mutate(Angle7 = case_when( # Split Angles
      D7 >= 0 & D7 < 90 ~ "NE", 
      D7 >= 90 & D7 < 180 ~ "NW",
      D7 >= 180 & D7 < 270 ~ "SW",
      D7 >= 270 & D7 < 360 ~ "SE")) %>%
    mutate(Angle7 = as.factor(Angle7)) %>%
    mutate(Angle7 = factor(Angle7, levels(Angle7)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle8 = case_when( # Split Angles
      D8 >= 0 & D8 < 90 ~ "NE", 
      D8 >= 90 & D8 < 180 ~ "NW",
      D8 >= 180 & D8 < 270 ~ "SW",
      D8 >= 270 & D8 < 360 ~ "SE")) %>%
    mutate(Angle8 = as.factor(Angle8)) %>%
    mutate(Angle8 = factor(Angle8, levels(Angle8)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle9 = case_when( # Split Angles
      D9 >= 0 & D9 < 90 ~ "NE", 
      D9 >= 90 & D9 < 180 ~ "NW",
      D9 >= 180 & D9 < 270 ~ "SW",
      D9 >= 270 & D9 < 360 ~ "SE")) %>%
    mutate(Angle9 = as.factor(Angle9)) %>%
    mutate(Angle9 = factor(Angle9, levels(Angle9)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle10 = case_when( # Split Angles
      D10 >= 0 & D10 < 90 ~ "NE", 
      D10 >= 90 & D10 < 180 ~ "NW",
      D10 >= 180 & D10 < 270 ~ "SW",
      D10 >= 270 & D10 < 360 ~ "SE")) %>%
    mutate(Angle10 = as.factor(Angle10)) %>%
    mutate(Angle10 = factor(Angle10, levels(Angle10)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle11 = case_when( # Split Angles
      D11 >= 0 & D11 < 90 ~ "NE", 
      D11 >= 90 & D11 < 180 ~ "NW",
      D11 >= 180 & D11 < 270 ~ "SW",
      D11 >= 270 & D11 < 360 ~ "SE")) %>%
    mutate(Angle11 = as.factor(Angle11)) %>%
    mutate(Angle11 = factor(Angle11, levels(Angle11)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle12 = case_when( # Split Angles
      D12 >= 0 & D12 < 90 ~ "NE", 
      D12 >= 90 & D12 < 180 ~ "NW",
      D12 >= 180 & D12 < 270 ~ "SW",
      D12 >= 270 & D12 < 360 ~ "SE")) %>%
    mutate(Angle12 = as.factor(Angle12)) %>%
    mutate(Angle12 = factor(Angle12, levels(Angle12)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle13 = case_when( # Split Angles
      D13 >= 0 & D13 < 90 ~ "NE", 
      D13 >= 90 & D13 < 180 ~ "NW",
      D13 >= 180 & D13 < 270 ~ "SW",
      D13 >= 270 & D13 < 360 ~ "SE")) %>%
    mutate(Angle13 = as.factor(Angle13)) %>%
    mutate(Angle13 = factor(Angle13, levels(Angle13)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle14 = case_when( # Split Angles
      D14 >= 0 & D14 < 90 ~ "NE", 
      D14 >= 90 & D14 < 180 ~ "NW",
      D14 >= 180 & D14 < 270 ~ "SW",
      D14 >= 270 & D14 < 360 ~ "SE")) %>%
    mutate(Angle14 = as.factor(Angle14)) %>%
    mutate(Angle14 = factor(Angle14, levels(Angle14)[c(1, 2, 4, 3)])) %>% 
    
    mutate(Angle15 = case_when( # Split Angles
      D15 >= 0 & D15 < 90 ~ "NE", 
      D15 >= 90 & D15 < 180 ~ "NW",
      D15 >= 180 & D15 < 270 ~ "SW",
      D15 >= 270 & D15 < 360 ~ "SE")) %>%
    mutate(Angle15 = as.factor(Angle15)) %>%
    mutate(Angle15 = factor(Angle15, levels(Angle15)[c(1, 2, 4, 3)])) %>% 
    
    arrange(Year, Month)
  # write.csv(JellyDat, file = "JellyDat_CairnsBeach.csv")
  
  # Check
  # clifro::windrose(JellyDat$S1, SwitchDirn(JellyDat$D1), "CairnsBeaches") # Observed
  # clifro::windrose(JellyDat$S15[JellyDat$Stings > 0], SwitchDirn(JellyDat$D15[JellyDat$Stings > 0]), "CairnsBeaches") # Observed
  JellyDat <- JellyDat[-c(1:14), ] # Remove 1st 14 rows as these cannot have 15 day cumsum of winds
  return(JellyDat)
}
#########################

# PRELIMINARIES
dat <- read.csv("VJD_records_EXTRACT_20180802_QLD.csv", header = TRUE)
      
missmap(dat) # Lots of missing data
graphics.off() # need to do after the missmap

hist(hour(ymd_hm(dat$EVENT_DATE))) # Most stings at 3 pm. 0s are missing data. If estimating MD for TOD, then 3 pm is mode
table(dat$EVENT_ONSHORE_OFFSHORE) # the one bay is in Moreton Bay and gets removed; the 14 unknowns are in 1956/57, so get removed. The 139 blanks get coded later
summary(dat$REGION)
aggregate(dat$LAT, by = list(dat$REGION), range) # Remove SEQ

dat <- dat %>%
  filter(EVENT_TYPE == "sting" & YEAR >= 1985) %>% # Only stings after 1985 when data are more complete
  dplyr::select(CSIRO_ID, EVENT_DATE, YEAR, MONTH, DAY, LAT, LON, EVENT_ONSHORE_OFFSHORE, SITE_INFO) %>% # Key variables for analysis
  rename(ID = CSIRO_ID, Event_Date = EVENT_DATE, Year = YEAR, Month = MONTH, 
         Day = DAY, Latitude = LAT, Longitude = LON,
         Area = EVENT_ONSHORE_OFFSHORE, Site = SITE_INFO) %>%
  
  filter(!(Year == 1985 & Month == "JAN")) %>% # Only 1 sting before December 1985, so drop this
  filter(Year != 2017) %>% # Only 1 sting in 2017, so drop this, so 1985-2016
  
  mutate(Date = dmy_hm(Event_Date)) %>%
  mutate(Event_Date = dmy_hm(Event_Date)) %>%
  mutate(Date = date(Date)) %>% # Removes time part
  mutate(Year = as.factor(Year)) %>%
  mutate(Month = month(Event_Date, label = FALSE)) %>% # Month as a number
  mutate(Month = as.factor(Month)) %>%
  mutate(Time = hour(Event_Date) + minute(Event_Date)/60) %>% # Time part
  mutate(TimeOriginal = Time) %>% # Original time
  mutate(Time = replace(Time, Time < 2, 15)) %>% # Replace MD with 15h00, the most common swim time. Change one value at 1 am to 15h00 too
  mutate(TimeOriginal = replace(TimeOriginal, TimeOriginal < 2, values = NA)) %>% # MD code for 0s in original time variable
  mutate(Area = fct_recode(Area, "Beach" = "beach", "Island" = "island", "Reef" = "reef")) %>% # Rename areas to caps first letter (to match my file where I found areas for lat and lon)
  mutate(DOY = yday(Event_Date)) %>% # Normal DOY
  mutate(DOYSeason = (DOY <= 182) * (DOY + 182) + (DOY >= 183) * (DOY - 183)) %>% # DOY, starting in July (stinger season)
  mutate(Yr = case_when(Month == "1" | Month == "2" | Month == "3" | Month == "4" | Month == "5" | Month == "6" ~ as.numeric(as.character(Year)) - 1,
                        TRUE ~ as.numeric(as.character(Year)))) %>% # Make stinger season
  mutate(Yr = as.factor(Yr)) %>%
  mutate(Region = case_when( # Split data into 3 regions - as identified by map
    Latitude < -14 & Latitude > -17.5    ~ "Cairns",
    Latitude <= -17.5 & Latitude > -19.5 ~ "Townsville",
    Latitude <= -19.5 & Latitude > -22   ~ "Whitsundays",
  TRUE                                 ~ "Other")) %>% # Outside AOI
  mutate(Region = as.factor(Region)) %>%
  filter(Region != "Other")

missmap(dat)
table(dat$Region, dat$Area)
aggregate(dat$Longitude, by = list(dat$Region, dat$Area), mean)

### Assign MD Areas (n = 83) - assign Beach, Island or Reef from Google Earth  ###
table(dat$Area)
dat$Area[dat$Area == ""] <- NA # Replace blanks with NA
Lats <- dat$Latitude[is.na(dat$Area)] # Get their Lats, Lons and IDs
Lons <- dat$Longitude[is.na(dat$Area)]
IDs <- dat$ID[is.na(dat$Area)]
MD <- read.csv("MissingDataArea.csv", header = TRUE) # From GoogleEarth

for(i in 1:nrow(MD)){ # Put in new Areas and Sites
  Index <- dat$ID == MD$CSIRO_ID[i]
  dat$Area[Index] <- MD$Area[i]
  dat$Site[Index] <- MD$Site[i]
} 
dat$Area <- droplevels(dat$Area) # Get rid of blank levels

### Create Region2 for Cairns Beach, Cairns Reef, etc. ###
# Now thta Areas are sorted...
dat <- dat %>% 
  mutate(Region2 = case_when( # Split data into Cairns (Beach, Island, Reef), Townsville (Beach, Island+Reef), Whitsundays (Beach, Island+Reef)  
  Region == "Cairns" & Area == "Beach" ~ "CairnsBeach",
  Region == "Cairns" & Area == "Island" ~ "CairnsIsland",
  Region == "Cairns" & Area == "Reef" ~ "CairnsReef",
  Region == "Townsville" & Area == "Beach" ~ "TownsvilleBeach",
  Region == "Townsville" & Area != "Beach" ~ "TownsvilleIslandReef",
  Region == "Whitsundays" & Area == "Beach" ~ "WhitsundaysBeach",
  Region == "Whitsundays" & Area != "Beach" ~ "WhitsundaysIslandReef"))

missmap(dat) # Some missing TimeOriginal and Site - both OK

table(dat$Region2)
sum(dat$Region == "Townsville" & dat$Area == "Island")
sum(dat$Region == "Townsville" & dat$Area == "Reef")
sum(dat$Region == "Whitsundays" & dat$Area == "Island")
sum(dat$Region == "Whitsundays" & dat$Area == "Reef")

table(dat$Month[dat$Area == "Beach"]) # Most stings in Dec and Jan
table(dat$Month[dat$Region == "Townsville"]) # Most stings in December
table(dat$Month[dat$Region == "Whitsundays"]) # 

table(dat$Region)
sum(table(dat$Region))

sum(!is.na(dat$TimeOriginal))
# Seasonality
####################################
# Make a dataset for analysis
# Just do CairnsBeach first
# Function for MakeData

Files_Wind <- list.files(path = "./WindData")
Files_Wind

dat_CairnsBeach <- MakeData(dat, "CairnsBeach", "file_146.25_-17.25.txt")
dat_CairnsIsland <- MakeData(dat, "CairnsIsland", "file_146.25_-17.25.txt")
dat_CairnsReef <- MakeData(dat, "CairnsReef", "file_146.25_-17.25.txt")
dat_CairnsIslandReef <- MakeData(dat, "CairnsIslandReef", "file_146.25_-17.25.txt")
dat_TownsvilleBeach <- MakeData(dat, "TownsvilleBeach", "file_147_-18.75.txt")
dat_TownsvilleIslandReef <- MakeData(dat, "TownsvilleIslandReef", "file_147_-18.75.txt")
dat_WhitsundaysBeach <- MakeData(dat, "WhitsundaysBeach", "file_149.25_-20.25.txt")
dat_WhitsundaysIslandReef <- MakeData(dat, "WhitsundaysIslandReef", "file_149.25_-20.25.txt")

# Function to add hourly tide data
# Do it to dat_CairnsBeach first, and then put in overall functon call
Files_Tide <- list.files(path = "./TideData") # Files_Tide[2] is Cairns
Files_Tide
dat_CairnsBeach <- AddTides(dat_CairnsBeach, 2)
dat_CairnsIsland <- AddTides(dat_CairnsIsland, 2)
dat_CairnsReef <- AddTides(dat_CairnsReef, 2)
dat_CairnsIslandReef <- AddTides(dat_CairnsIslandReef, 2)
dat_TownsvilleBeach <- AddTides(dat_TownsvilleBeach, 7)
dat_TownsvilleIslandReef <- AddTides(dat_TownsvilleIslandReef, 7)
dat_WhitsundaysBeach <- AddTides(dat_WhitsundaysBeach, 7)
dat_WhitsundaysIslandReef <- AddTides(dat_WhitsundaysIslandReef, 7)

# Add tidal range
dat_CairnsBeach <- AddTidalRange(dat_CairnsBeach, 2)
dat_CairnsIsland <- AddTidalRange(dat_CairnsIsland, 2)
dat_CairnsReef <- AddTidalRange(dat_CairnsReef, 2)
dat_CairnsIslandReef <- AddTidalRange(dat_CairnsIslandReef, 2)
dat_TownsvilleBeach <- AddTidalRange(dat_TownsvilleBeach, 7)
dat_TownsvilleIslandReef <- AddTidalRange(dat_TownsvilleIslandReef, 7)
dat_WhitsundaysBeach <- AddTidalRange(dat_WhitsundaysBeach, 7)
dat_WhitsundaysIslandReef <- AddTidalRange(dat_WhitsundaysIslandReef, 7)

### ANALYSIS
# MODEL LEARNING
# 1. Models only fits if remove months with no data - 6, 7, 8
# 2. Poisson  and negative binomial fit. Negative binomial is hard to plot
# 3. Models are not overdispersed. Residuals not great though
# 4. Because wind direction (D) is so non-linear, better to use a factor (Angle)
# 5. Looks like longer cumulative winds (both S and D) give higher r2
# 6. Non-linear term (ns) for Speed works best, but only df=2
# 7. r2 ~ 0.25
# 8. r2 goes up 6% with Year in as a factor (excluding missing years). Should I include Year? Not got for prediction...
# 8. Note that forward stepwise regression doesn't work with S1, S2, etc. as you get 
# alternating declining and increasing trends. Better to try one at a time and see which
# has highest r2. Note that they are highly correlated with each other...
# mf <- stepAIC(m1, scope = list(lower = ~1, upper = ~ Month + Angle + ns(S1, df=2) + ns(S2, df=2) + ns(S3, df=2) + ns(S4, df=2) + ns(S5, df=2) + ns(S10, df=2)), direction = "forward", trace = 10)

# TO DO
# 1. Need to reorganise Month so StingerMonth - done
# 2. Package %>% more neatly for the model fitting - done
# 3. Include tides
# 4. Include mulitiple cumulative Angles - done
# 4. Stepwise model selection of Speeds and Angles - done
# dat_CairnsBeach <- dat_CairnsBeach %>%
#                    mutate(Month = factor(Month, levels(Month)[c(1,5,6,7,8,9,10,11,12,2,3,4)]))
### ****************************** ###  
# No stings in June, July, August for CairnsBeach
aggregate(dat_CairnsBeach$Stings, by = list(dat_CairnsBeach$Month), sum) 
with(dat_CairnsBeach, table(Stings, Month))
with(dat_CairnsBeach, table(Stings, Angle3))

dat_CairnsBeach_NoWinter <- dat_CairnsBeach %>%
  filter(Month !=6 & Month != 7 & Month != 8) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

with(dat_CairnsBeach_NoWinter, table(Stings))

# Presence Absence data
# dat_CairnsBeach_NoWinter$PresAbs <- dat_CairnsBeach_NoWinter$Stings
# dat_CairnsBeach_NoWinter$PresAbs[dat_CairnsBeach_NoWinter$Stings > 0] <- 1 

########## CAIRNS BEACH ####################
# Best model
# m1 <- glm(Stings ~ Mn + Angle3 + ns(S3, df=2) + Incoming * ns(Height, df=3) + TidalRange, data = dat_CairnsBeach_NoWinter, family = poisson)
# Tidal Range weak and opposite to expected - remove
m1 <- glm(Stings ~ Month + Direction3d + ns(Speed3d, df=2) + TidalRange + Incoming * ns(Height, df=3), data = dat_CairnsBeach_NoWinter, family = poisson)

1 - m1$deviance/m1$null.deviance # r2=28.2%
dropterm(m1)
m1 <- update(m1, . ~ . -TidalRange)
dropterm(m1)
quartz(width = 8, height = 6)
plot(allEffects(m1, xlevels=list(Speed3d=seq(0, 15, 5), Height=seq(-1,1.5,0.5))), main = "")
dev.copy2pdf(file = "CairnsBeaches.pdf")


# Try Presence/Absence model
m1 <- glm(PresAbs ~ Mn + ns(S3, df=2) + Angle3 + Incoming * ns(Height, df=3), data = dat_CairnsBeach_NoWinter, family = binomial)


########## CAIRNS ISLAND ####################
# No stings in June, July, August for CairnsBeach
# Only do presence/absence - not enough data for Poisson

with(dat_CairnsIsland, table(Stings, Month)) # No 6, 7, 8
with(dat_CairnsIsland, table(Stings, Angle3)) # No SW

dat_CairnsIsland_NoWinter <- dat_CairnsIsland %>%
  filter(Month !=6 & Month != 7 & Month != 8) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# Add 1  sting to NE wind the first occurrence
# Better to remove them
# dat_CairnsIsland_NoWinter$Stings[dat_CairnsIsland_NoWinter$Direction3d == "NE"][1] <- 1
# with(dat_CairnsIsland_NoWinter, table(Stings, Direction3d)) # Check

# Remove all NE Winds - infrequent
dat_CairnsIsland_NoWinter <- dat_CairnsIsland_NoWinter[dat_CairnsIsland_NoWinter$Direction3d != "NE",]

# Try coding winds as SW or NW as 1 and other winds as 0
# Does not work any better...
# dat_CairnsIsland_NoWinter$Dirn3d <- 0
# dat_CairnsIsland_NoWinter$Dirn3d[dat_CairnsIsland_NoWinter$Direction3d == "SW"] <- 1
# dat_CairnsIsland_NoWinter$Dirn3d[dat_CairnsIsland_NoWinter$Direction3d == "NW"] <- 1
# dat_CairnsIsland_NoWinter$Dirn3d <- as.factor(dat_CairnsIsland_NoWinter$Dirn3d)

# Presence Absence data
# Count is better where possible
# dat_CairnsIsland_NoWinter$PresAbs <- dat_CairnsIsland_NoWinter$Stings
# dat_CairnsIsland_NoWinter$PresAbs[dat_CairnsIsland_NoWinter$Stings > 0] <- 1 

# dat_CairnsIsland_NoWinter <- dat_CairnsIsland_NoWinter %>%
#  rename(Presence = PresAbs)

# m2 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming*ns(Height, df=3) + TidalRange, data = dat_CairnsIsland_NoWinter, family = binomial)
# m2 <- glm(Presence ~ Direction3d + Incoming*ns(Height, df=3) + TidalRange, data = dat_CairnsIsland_NoWinter, family = binomial)

m2 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming*ns(Height, df=3) + TidalRange, data = dat_CairnsIsland_NoWinter, family = poisson)
# m2 <- glm(Stings ~ Month + ns(Speed3d, df=2) + ns(D3, df = 10) + Incoming*ns(Height, df=3) + TidalRange, data = dat_CairnsIsland_NoWinter, family = poisson)
# m2 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Dirn3d + Incoming*ns(Height, df=3) + TidalRange, data = dat_CairnsIsland_NoWinter, family = poisson)

# library(pscl)
# Tried zero inflated models but deosn't work any better
# m2 <- zeroinfl(Stings ~ Month + ns(Speed3d, df=2), data = dat_CairnsIsland_NoWinter, family = poisson)
# visreg(m2, "Month", type = "conditional", rug = FALSE)
# m2 <- hurdle(Stings ~ Incoming*ns(Height, df=3) + TidalRange | Month + Direction3d + ns(Speed3d, df=2), dist = "poisson", data = dat_CairnsIsland_NoWinter)

# library(partykit)
# # Tried partykit (CARTs) but doesn't work any better
# m2 <- ctree(Presence ~ Speed3d + Direction3d + TidalRange, data = dat_CairnsIsland_NoWinter)
# m2 <- ctree(Stings ~ Direction3d, data = dat_CairnsIsland_NoWinter)
# m2 <- ctree(Presence ~ Direction3d, data = dat_CairnsIsland_NoWinter)
# m2
# plot(m2, gp = gpar(fontsize = 7))

summary(m2)
dropterm(m2)
m2 <- update(m2, . ~ . -Direction3d)
dropterm(m2)
summary(m2)
1 - m2$deviance/m2$null.deviance # r2 = 8.6%

quartz(width = 8, height = 6)
plot(allEffects(m2, xlevels=list(Speed3d=seq(0, 15, 5), Height=seq(-1,1.5,0.5))), main = "")
dev.copy2pdf(file = "CairnsIslands.pdf")

########## CAIRNS REEFS ####################
with(dat_CairnsReef, table(Stings, Month)) # No 6, 7, 8, 9
with(dat_CairnsReef, table(Stings, Angle3)) # Have alll angles with stings

dat_CairnsReef_NoWinter <- dat_CairnsReef %>%
  filter(Month !=6 & Month != 7 & Month != 8 & Month != 9) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# dat_CairnsReef_NoWinter$PresAbs <- dat_CairnsReef_NoWinter$Stings
# dat_CairnsReef_NoWinter$PresAbs[dat_CairnsReef_NoWinter$Stings > 0] <- 1 

# dat_CairnsReef_NoWinter <- dat_CairnsReef_NoWinter %>%
#  rename(Presence = PresAbs)

# m3 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_CairnsReef_NoWinter, family = binomial)
# summary(m3)
m3 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_CairnsReef_NoWinter, family = poisson)

dropterm(m3)
m3 <- update(m3, . ~ . -Incoming:ns(Height, df=3))
dropterm(m3)
m3 <- update(m3, . ~ . -Direction3d)
dropterm(m3)

1 - m3$deviance/m3$null.deviance # r2=3.8%
quartz(width = 8, height = 6)
plot(allEffects(m3, xlevels=list(Speed3d=seq(0, 15, 5), Height=seq(-1,1.5,0.5))), main = "")


########## CAIRNS ISLANDS & REEFS ####################
dat_CairnsIslandReef_NoWinter <- dat_CairnsIslandReef %>%
  filter(Month !=6 & Month != 7 & Month != 8 & Month != 9) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)
# mutate(Angle3 = fct_recode(Angle3, "NW" = "NE"))

with(dat_CairnsIslandReef_NoWinter, table(Stings, Mn))
with(dat_CairnsIslandReef_NoWinter, table(Stings, Angle3))

dat_CairnsIslandReef_NoWinter$PresAbs <- dat_CairnsIslandReef_NoWinter$Stings
dat_CairnsIslandReef_NoWinter$PresAbs[dat_CairnsIslandReef_NoWinter$Stings > 0] <- 1 

dat_CairnsIslandReef_NoWinter <- dat_CairnsIslandReef_NoWinter %>%
  rename(Presence = PresAbs)

m3.1 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_CairnsIslandReef_NoWinter, family = binomial)
summary(m3.1)
dropterm(m3.1)
m3 <- update(m3.1, . ~ . -Direction3d)
dropterm(m3.1)
m3 <- update(m3.1, . ~ . -Incoming:ns(Height, df=3))
dropterm(m3.1)
m3 <- update(m3.1, . ~ . -TidalRange)
dropterm(m3.1)
m3 <- update(m3.1, . ~ . -Incoming)
dropterm(m3.1)
summary(m3.1)
m3 <- update(m3.1, . ~ . -ns(Speed3d, df=2))
dropterm(m3.1)
1 - m3.1$deviance/m3.1$null.deviance # r2=3.8%
quartz(width = 8, height = 6)
plot(allEffects(m3.1, xlevels=list(Speed3d=seq(0, 15, 5), Height=seq(-1,1.5,0.5))), main = "")

########## TOWNSVILLE BEACH  ####################
with(dat_TownsvilleBeach, table(Stings, Month)) # No 5, 6, 7, 8, 10, 2
with(dat_TownsvilleBeach, table(Stings, Angle3)) # No NE

dat_TownsvilleBeach_NoWinter <- dat_TownsvilleBeach %>%
  filter(Month !=6 & Month != 7 & Month != 8 & Month != 10 & Month != 5 & Month != 2) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  # mutate(Angle3 = fct_recode(Angle3, "NW" = "NE")) %>%
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# Remove all NE Winds - infrequent
dat_TownsvilleBeach_NoWinter <- dat_TownsvilleBeach_NoWinter[dat_TownsvilleBeach_NoWinter$Direction3d != "NE",]

with(dat_TownsvilleBeach_NoWinter, table(Stings, Month)) # Missing Months: 2, 5-10
with(dat_TownsvilleBeach_NoWinter, table(Stings, Direction3d)) # Missing Angles: NE

# dat_TownsvilleBeach_NoWinter$PresAbs <- dat_TownsvilleBeach_NoWinter$Stings
# dat_TownsvilleBeach_NoWinter$PresAbs[dat_TownsvilleBeach_NoWinter$Stings > 0] <- 1 

# dat_TownsvilleBeach_NoWinter <- dat_TownsvilleBeach_NoWinter %>%
#  rename(Presence = PresAbs)

m4 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=2) + TidalRange, data = dat_TownsvilleBeach_NoWinter, family = binomial)
# NOTE: S4 and higher, and Angle4 and higher lead to huge SEs - Angles disappear
m4 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=2) + TidalRange, data = dat_TownsvilleBeach_NoWinter, family = poisson)
1 - m4$deviance/m4$null.deviance # 7.7%
dropterm(m4)
m4 <- update(m4, . ~ . -Incoming:ns(Height, df=2))
dropterm(m4)
m4 <- update(m4, . ~ . -ns(Height, df=2))
dropterm(m4)

1 - m4$deviance/m4$null.deviance # r2=3.8%
quartz(width = 8, height = 6)
plot(allEffects(m4, xlevels=list(Speed3d=seq(0, 15, 5), Height=seq(-1,1.5,0.5))), main = "")
summary(m4)

########## TOWNSVILLE ISLANDREEF  ####################
with(dat_TownsvilleIslandReef, table(Stings, Month)) # No 7, 8
with(dat_TownsvilleIslandReef, table(Stings, Angle3)) # No NE

dat_TownsvilleIslandReef_NoWinter <- dat_TownsvilleIslandReef %>%
  filter(Month != 7 & Month != 8) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  # mutate(Angle3 = fct_recode(Angle3, "NW" = "NE")) %>%
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# Remove all NE Winds - infrequent
dat_TownsvilleIslandReef_NoWinter <- dat_TownsvilleIslandReef_NoWinter[dat_TownsvilleIslandReef_NoWinter$Direction3d != "NE",]

# with(dat_TownsvilleIslandReef_NoWinter, table(Stings, Mn)) # Months missing: 6, 7, 8
# with(dat_TownsvilleIslandReef_NoWinter, table(Stings, Angle3)) # All Angles

# dat_TownsvilleIslandReef_NoWinter$PresAbs <- dat_TownsvilleIslandReef_NoWinter$Stings
# dat_TownsvilleIslandReef_NoWinter$PresAbs[dat_TownsvilleIslandReef_NoWinter$Stings > 0] <- 1 

# dat_TownsvilleIslandReef_NoWinter <- dat_TownsvilleIslandReef_NoWinter %>%
#  rename(Presence = PresAbs)

# m5 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_TownsvilleIslandReef_NoWinter, family = binomial)
m5 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_TownsvilleIslandReef_NoWinter, family = poisson)
1 - m5$deviance/m5$null.deviance
dropterm(m5)
m5 <- update(m5, . ~ . -ns(Speed3d, df=2))
dropterm(m5)
1 - m5$deviance/m5$null.deviance
summary(m5)
quartz(width = 8, height = 6)
plot(allEffects(m5, xlevels=list(Height=seq(-1,1.5,0.5))), main = "")

########## WHITSUNDAYS BEACH  ####################
with(dat_WhitsundaysBeach, table(Stings, Month)) # No 6-11
with(dat_WhitsundaysBeach, table(Stings, Angle3)) # No NE

dat_WhitsundaysBeach_NoWinter <- dat_WhitsundaysBeach %>%
  filter(Month !=6 & Month != 7 & Month != 8 & Month != 9 & Month != 10 & Month != 11) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  # mutate(Angle3 = fct_recode(Angle3, "NW" = "NE")) %>%
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# Remove all NE Winds - infrequent
dat_WhitsundaysBeach_NoWinter <- dat_WhitsundaysBeach_NoWinter[dat_WhitsundaysBeach_NoWinter$Direction3d != "NE",]

# dat_WhitsundaysBeach_NoWinter$PresAbs <- dat_WhitsundaysBeach_NoWinter$Stings
# dat_WhitsundaysBeach_NoWinter$PresAbs[dat_WhitsundaysBeach_NoWinter$Stings > 0] <- 1 

# dat_WhitsundaysBeach_NoWinter <- dat_WhitsundaysBeach_NoWinter %>%
#  rename(Presence = PresAbs)

# Tidal Range n.s.
# m6 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_WhitsundaysBeach_NoWinter, family = binomial)
m6 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_WhitsundaysBeach_NoWinter, family = poisson)
summary(m6)
1 - m6$deviance/m6$null.deviance
plot(allEffects(m6))
dropterm(m6)
m6 <- update(m6, . ~ . -Incoming:ns(Height, df=3))
dropterm(m6)
m6 <- update(m6, . ~ . -ns(Speed3d, df = 2))
dropterm(m6)
m6 <- update(m6, . ~ . -Incoming)
dropterm(m6)
m6 <- update(m6, . ~ . -ns(Height, df=3))
dropterm(m6)
m6 <- update(m6, . ~ . -TidalRange)
dropterm(m6)

1 - m6$deviance/m6$null.deviance # r2=4.9%
quartz(width = 8, height = 6)
plot(allEffects(m6), main = "")


########## WHITSUNDAYS ISLANDREEF  ####################
with(dat_WhitsundaysIslandReef, table(Stings, Month)) # No 8,10
with(dat_WhitsundaysIslandReef, table(Stings, Angle3)) # No NE

dat_WhitsundaysIslandReef_NoWinter <- dat_WhitsundaysIslandReef %>%
  # filter(Month != 6 & Month != 7 & Month != 8 & Month != 10) %>%
  filter(Month != 8 & Month != 10) %>%
  mutate(S3 = (S3 > 15)*15 + (S3 <= 15)*S3) %>% # Cutoff at 15
  mutate(Height = (Height > 1.5)*1.5 + (Height <= 1.5)*Height) %>% # Cutoff at 15
  # mutate(Angle3 = fct_recode(Angle3, "NW" = "NE")) %>%
  rename(Month2 = Month, Month = Mn, Direction3d = Angle3, Speed3d =  S3)

# Remove all NE Winds - infrequent
dat_WhitsundaysIslandReef_NoWinter <- dat_WhitsundaysIslandReef_NoWinter[dat_WhitsundaysIslandReef_NoWinter$Direction3d != "NE",]

# dat_WhitsundaysIslandReef_NoWinter$PresAbs <- dat_WhitsundaysIslandReef_NoWinter$Stings
# dat_WhitsundaysIslandReef_NoWinter$PresAbs[dat_WhitsundaysIslandReef_NoWinter$Stings > 0] <- 1 

# dat_WhitsundaysIslandReef_NoWinter <- dat_WhitsundaysIslandReef_NoWinter %>%
#  rename(Presence = PresAbs)

# Tidal Range n.s.
# m7 <- glm(Presence ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_WhitsundaysIslandReef_NoWinter, family = binomial)
m7 <- glm(Stings ~ Month + ns(Speed3d, df=2) + Direction3d + Incoming * ns(Height, df=3) + TidalRange, data = dat_WhitsundaysIslandReef_NoWinter, family = poisson)

summary(m7)
1 - m7$deviance/m7$null.deviance
plot(allEffects(m7))
dropterm(m7)
m7 <- update(m7, . ~ . -TidalRange)
dropterm(m7)
m7 <- update(m7, . ~ . -ns(Speed3d, df = 2))
dropterm(m7)
1 - m7$deviance/m7$null.deviance # r2=13.2%
quartz(width = 8, height = 6)
plot(allEffects(m7, xlevels=list(Height=seq(-1.75,1.5,0.5))), main = "")




# Not interpretable
m1 <- lm(Stings ~ Month + ns(D1, df = 2)*ns(S1, df= 2), data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1))

m1 <- lm(Stings ~ Month + ns(D1, df=3) + ns(D2, df=3) + ns(D3, df=3), data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1))

m1 <- lm(Stings ~ Month + Harm(D1, k=4), data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1))

m1 <- lm(Stings ~ Harm(D1, k=4), data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1))

# Try poisson - does not fit properly
m1 <- glm(Stings ~ Month + ns(S1, df = 5) + ns(D1, df = 5), data = dat_CairnsBeach, family = "poisson")
summary(m1)
plot(allEffects(m1))

# Negbin - Does not fit properly for Jun, Jul, Aug
m1 <- glm.nb(Stings ~ Month + ns(S1, df = 5) + ns(D1, df = 5), data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1))

# Negbin - Remove: Jun, Jul, Aug - fits
m1 <- glm.nb(Stings ~ Month + ns(S1, df = 5) + ns(D1, df = 5), data = dat_CairnsBeach_MM)
summary(m1)

# Not good
m1 <- glm.nb(Stings ~ Month + ns(S11, df = 2) + ns(D11, df = 2), data = dat_CairnsBeach_NoWinter)
summary(m1)
plot(m1)

# Works OK - S11 and D11 best - Month has a big effect
m1 <- glm(Stings ~ Month + ns(S13, df = 2) + Harm(D13, k = 2), data = dat_CairnsBeach_NoWinter, family = "poisson")
summary(m1)
plot(allEffects(m1))

# Without winter - pretty good model
m1 <- glm(Stings ~ Month + ns(S11, df = 2) + Harm(D11, k = 2), data = dat_CairnsBeach_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))

# Try ns for D
m1 <- glm(Stings ~ Month + ns(S10, df = 2) + ns(D10, df = 5), data = dat_CairnsBeach_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))

# Try a surface - doesn't work too well
library(mgcv)
m1 <- gam(Stings ~ Month + te(S11, D11, bs = c("cr", "cc")), data = dat_CairnsBeach_MM, family = "poisson")
summary(m1)
vis.gam(m1, c("S11", "D11"), theta = 45, phi = 10, r = 100)

table(dat_CairnsBeach_NoWinter$Angle)

# Try angle as a factor # ns(S10, df = 2) - good model. S15 best
m1 <- glm(Stings ~ Month + ns(S15, df = 2) + Angle, data = dat_CairnsBeach_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))

# Try angle as a factor # ns(S10, df = 2) and interaction
m1 <- glm(Stings ~ Month + S10*Angle, data = dat_CairnsBeach_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))

# Try a hurdle model - does not fit properly
library(pscl)
m1 <- hurdle(Stings ~ Harm(D1, k=4) + ns(S1, df = 3)| Month, dist = "negbin", na.action = na.omit, data = dat_CairnsBeach_MM)
summary(m1)

# does not fit properly
m1 <- hurdle(Stings ~ Harm(D1, k=4) + ns(S1, df = 3)| Month, dist = "poisson", na.action = na.omit, data = dat_CairnsBeach)
summary(m1)

# Drop Months 
# Negative binomial does not converge
m2 <- glm.nb(Stings ~ Year + Month, data = dat_CairnsBeach)
summary(m1)
plot(allEffects(m1, las = 2))
par(mfrow = c(1,2))

# Try angle as a factor
m1 <- glm(Stings ~ Month + ns(S15, df = 2), data = dat_CairnsIsland_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))

m1 <- glm(Stings ~ Month + Angle, data = dat_CairnsIsland_NoWinter, family = "poisson")
summary(m1)
1 - m1$deviance/m1$null.deviance
plot(allEffects(m1))


####################################
###### Plot data using ggplot ######
aus <- map_data("world2Hires")
aus <- filter(aus, region == "Australia")

ggplot() + geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = "dark grey") + 
  coord_fixed(ratio = 1, xlim = c(140, 155), ylim = c(-30, -10)) +
  geom_hex(data = dat, aes(Longitude, Latitude), binwidth = c(0.5, 0.5)) +
  scale_fill_gradient(name = "# stings", trans = "log10", 
                      low = "light blue", high = "dark blue") +
  labs(x = "Longitude", y = "Latitude") + theme_bw()
ggplot() + geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = "dark grey") + 
  coord_fixed(ratio = 1, xlim = c(140, 150), ylim = c(-25, -10)) +
  geom_hex(data = dat, aes(Longitude, Latitude), binwidth = c(0.5, 0.5)) +
  scale_fill_gradient(name = "# stings", trans = "log10", 
                      low = "light blue", high = "dark blue") +
  labs(x = "Longitude", y = "Latitude") + theme_bw()

dev.copy2pdf(file = "StingMap2.pdf", paper = "A4")
######################################

###### Plot seasonal cycle using ggplot ######
# Use DOY and do for each region
# Work out where DOY is on the axis
DaysPerMonth <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30) # July to June
Temp <- c(0, DaysPerMonth) # Add zero for start of year
Temp <- Temp[1:(length(Temp) - 1)] # Remove last element
xpos <- cumsum(Temp) + DaysPerMonth / 2
xlab <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")

gg <- ggplot(dat) + geom_histogram(aes(DOYSeason), bins = 30) +
  labs(x = "Month", y = "Number of stings") + 
  scale_x_continuous(breaks = xpos, labels = xlab, limits = c(5, 360)) +
  theme_bw()
gg
dev.copy2pdf(file = "StingSeasonalCycle_All.pdf", paper = "A4")
gg + facet_grid(Area ~ ., scales = "free")
dev.copy2pdf(file = "StingSeasonalCycle_Areas.pdf", paper = "A4")

gg + facet_grid(Region ~ ., scales = "free")
dev.copy2pdf(file = "StingSeasonalCycle_Regions.pdf", paper = "A4")
##############################################

# Plot years
test <- as.data.frame(table(dat$Yr)) # Plot Stinger Yr (i.e. based on season - Jan-Jun is part of previous year)
test$Var1 <- as.numeric(test$Var1)
gg <- ggplot(data = test, aes(x = Var1, y = Freq)) + geom_line() +
  labs(x = "Year", y = "Number of stings") +
  scale_x_continuous(breaks = seq(0, 35, 5), labels = seq(1985, 2020, 5)) +
  theme_bw()
gg
dev.copy2pdf(file = "StingInterannulVariation_StingerYrs.pdf", paper = "A4")

### Detailed map of Cairns ### Offshore islands mainly off Cairns
###### Plot data for each RegionNew using ggplot ######
ggplot() + geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = "dark grey") + 
  coord_fixed(ratio = 1, xlim = c(145, 148), ylim = c(-17.5, -14.5)) +
  geom_hex(data = dat, aes(Longitude, Latitude), binwidth = c(0.1, 0.1)) +
  scale_fill_gradient(name = "# stings", trans = "log10", 
                      low = "light blue", high = "dark blue") +
  labs(x = "Longitude", y = "Latitude") + theme_bw()
dev.copy2pdf(file = "StingMapCairns.pdf", paper = "A4")

ggplot() + geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = "dark grey") + 
  coord_fixed(ratio = 1, xlim = c(145, 148), ylim = c(-19.5, -17.5)) +
  geom_hex(data = dat, aes(Longitude, Latitude), binwidth = c(0.1, 0.1)) +
  scale_fill_gradient(name = "# stings", trans = "log10", 
                      low = "light blue", high = "dark blue") +
  labs(x = "Longitude", y = "Latitude") + theme_bw()
dev.copy2pdf(file = "StingMapTownsville.pdf", paper = "A4")

ggplot() + geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = "dark grey") + 
  coord_fixed(ratio = 1, xlim = c(147.5, 149.5), ylim = c(-22, -19.5)) +
  geom_hex(data = dat, aes(Longitude, Latitude), binwidth = c(0.1, 0.1)) +
  scale_fill_gradient(name = "# stings", trans = "log10", 
                      low = "light blue", high = "dark blue") +
  labs(x = "Longitude", y = "Latitude") + theme_bw()
dev.copy2pdf(file = "StingMapWhitsundays.pdf", paper = "A4")

######################################
# Wind locations for each RegionNew
# WindLats <- c(-11.25, -12.75, -14.25, -15.75, -17.25, -18.75, -20.25, -21.75, -21.75, -23.25, -24.75)
# Latitude < -14 & Latitude > -17.5    ~ "Cairns", # -14.25, -15.75, -17.25,
# Latitude <= -17.5 & Latitude > -19.5 ~ "Townsville", # -18.75
# Latitude <= -19.5 & Latitude > -22   ~ "Whitsundays", #  -20.25, -21.75, -21.75,
# What is the correlation between 3 Cairns winds? And for 3 Whitsundays winds?
# Open wind data and get correlation

Files_Wind <- list.files(path = "./WindData")
# NOTE: 2 files with same lat but diff lon - no data south of -20.5
WindDat_14S144E <- read.table(paste("./WindData/", Files_Wind[3], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Cairns
WindDat_15S145E <- read.table(paste("./WindData/", Files_Wind[4], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Cairns
WindDat_17S146E <- read.table(paste("./WindData/", Files_Wind[5], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Cairns
WindDat_18S147E <- read.table(paste("./WindData/", Files_Wind[6], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Townsville
WindDat_20S149E <- read.table(paste("./WindData/", Files_Wind[7], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Whitsundays
WindDat_21S150E <- read.table(paste("./WindData/", Files_Wind[8], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Whitsundays
WindDat_21S151E <- read.table(paste("./WindData/", Files_Wind[9], sep = ""), sep = ",", header = FALSE, comment.char = '#') # Whitsundays

WindDat <- cbind(WindDat_14S144E[, 1:3], WindDat_15S145E[, 2:3], WindDat_17S146E[, 2:3], WindDat_18S147E[, 2:3], WindDat_20S149E[, 2:3], WindDat_21S150E[, 2:3], WindDat_21S151E[, 2:3])
names(WindDat) <- c("Date", "N14S", "E14S", "N15S", "E15S", "N17S", "E17S", "N18S", "E18S", "N20S", "E20S", "N21S", "E21S", "N21SOff", "E21SOff")

ggpairs(WindDat[, c(2, 4, 6)]) # N-S for 3 from Cairns. r = 0.91, r = 0.84, r = 0.65 (between 14 and 17) Take wind at 17S - most data from there
ggpairs(WindDat[, c(3, 5, 7)]) # E-W for 3 from Cairns. r > 0.83. Cairns is more correlated in E-W than N-S
ggpairs(WindDat[, c(6, 8)]) # N-S Cairns 17S with Townsville. r = 0.89
ggpairs(WindDat[, c(7, 9)]) # E-W Cairns 17S with Townsville.r = 0.93
ggpairs(WindDat[, c(10, 12, 14)]) # N-S for 3 from Whitsundays. r > 0.92
ggpairs(WindDat[, c(11, 13, 15)]) # E-W for 3 from Whitsundays. r > 0.92
ggpairs(WindDat[, c(6, 8, 15)]) # r > 0.92

### Plot windroses ###
# Calculate wind speed and direction first
WindDat$WindSpeedCairns <- sqrt(WindDat[,2]^2 + WindDat[,3]^2)   
WindDat$WindDirnCairns <- (atan2(WindDat[,2], WindDat[,3]) * (180 / pi)) + 180 # +180 converts the -ve angles (to 180o away for wind) and +ve ones too angles to 

WindDat$WindSpeedTownsville <- sqrt(WindDat[,8]^2 + WindDat[,9]^2)   
WindDat$WindDirnTownsville <- (atan2(WindDat[,8], WindDat[,9]) * (180 / pi)) + 180

WindDat$WindSpeedWhitsundays <- sqrt(WindDat[,10]^2 + WindDat[,11]^2)   
WindDat$WindDirnWhitsundays <- (atan2(WindDat[,10], WindDat[,11]) * (180 / pi)) + 180

# Convert from Wide to Long format to use facet wrap in wind rose plots
tempory <- gather(WindDat, "WindSpeedCairns", "WindSpeedTownsville", "WindSpeedWhitsundays", key = "Region", value = "Speed")
tempory2 <- gather(WindDat, "WindDirnCairns", "WindDirnTownsville", "WindDirnWhitsundays", key = "Region", value = "Dirn")
WindDatAll <- cbind(tempory[, c(1, 19, 20)], tempory2[, c(18, 20)])
WindDatAll <- WindDatAll %>% mutate(Region = fct_recode(Region, "Cairns" = "WindSpeedCairns", 
                                                        "Townsville" = "WindSpeedTownsville", 
                                                        "Whitsundays" = "WindSpeedWhitsundays")) %>%
                             # dplyr::select(-WindSpeedWhit) %>%
                             mutate(Month = month(WindDatAll$Date)) %>%
                             mutate(Month = as.factor(Month))

# Plot a simple windrose using all the defaults, ignoring any facet variable
with(WindDatAll, clifro::windrose(Speed, Dirn))

# Create custom speed bins, add a legend title, and change to a B&W theme
with(WindDatAll, clifro::windrose(Speed, Dirn,
                       speed_cuts = SpeedCuts,
                       legend_title = "Wind Speed\n(m/s)",
                       legend.title.align = .5,
                       ggtheme = "bw",
                       col_pal = "GnBu")) # "Greys"

# Note that underscore-separated arguments come from the windrose method, and
# period-separated arguments come from ggplot2::theme().

# Include a facet variable with one level
with(WindDatAll, clifro::windrose(Speed, Dirn, "QLD"))
# Plot a windrose for each level of the facet variable (each station)
with(WindDatAll, clifro::windrose(Speed, Dirn, Region, n_col = 3, 
                          speed_cuts = SpeedCuts,
                          legend_title = "Wind Speed\n(m/s)",
                          legend.title.align = .5,
                          ggtheme = "bw",
                          col_pal = "GnBu"))
ggsave("QLD_windrose_region.png", dpi = 600)
with(WindDatAll, clifro::windrose(Speed, Dirn, Month, n_col = 4, 
                          speed_cuts = SpeedCuts,
                          legend_title = "Wind Speed\n(m/s)",
                          legend.title.align = .5,
                          ggtheme = "bw",
                          col_pal = "GnBu"))
ggsave("QLD_windrose_season.png", dpi = 600)
rm(WindDat, WindDat_14S144E, WindDat_15S145E, WindDat_17S146E, WindDat_18S147E, WindDat_20S149E, WindDat_21S150E, WindDat_21S151E, tempory, tempory2, WindDatAll)
### Reefs and Region
table(dat$Area, dat$Region)

### Compare windroses with and w/o stings
with(dat_CairnsBeach, clifro::windrose(S1, SwitchDirn(D1),
                                  speed_cuts = SpeedCuts,
                                  legend_title = "Wind Speed\n(m/s)",
                                  legend.title.align = .5,
                                  ggtheme = "bw",
                                  col_pal = "GnBu"))
# Stings
quartz(width = 10, height = 8)
S0 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings == 0], SwitchDirn(dat_CairnsBeach$D1[dat_CairnsBeach$Stings == 0]), "0 stings")
S1 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings == 1], SwitchDirn(dat_CairnsBeach$D1[dat_CairnsBeach$Stings == 1]), "1 sting")
S2 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings == 2], SwitchDirn(dat_CairnsBeach$D1[dat_CairnsBeach$Stings == 2]), "2 stings")
S3 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings > 2], SwitchDirn(dat_CairnsBeach$D1[dat_CairnsBeach$Stings > 2]), ">2 stings")
grid.arrange(S0, S1, S2, S3, ncol = 2)
dev.copy2pdf(file = "CairnsBeachesStingsByWind2.pdf", paper = "A4r")

# Stings - remove Winter
quartz(width = 10, height = 8)
S0 <- PlotRose(dat_CairnsBeach_NoWinter$S1[dat_CairnsBeach_NoWinter$Stings == 0], SwitchDirn(dat_CairnsBeach_NoWinter$D3[dat_CairnsBeach_NoWinter$Stings == 0]), "0 stings")
S1 <- PlotRose(dat_CairnsBeach_NoWinter$S1[dat_CairnsBeach_NoWinter$Stings == 1], SwitchDirn(dat_CairnsBeach_NoWinter$D3[dat_CairnsBeach_NoWinter$Stings == 1]), "1 sting")
S2 <- PlotRose(dat_CairnsBeach_NoWinter$S1[dat_CairnsBeach_NoWinter$Stings == 2], SwitchDirn(dat_CairnsBeach_NoWinter$D3[dat_CairnsBeach_NoWinter$Stings == 2]), "2 stings")
S3 <- PlotRose(dat_CairnsBeach_NoWinter$S1[dat_CairnsBeach_NoWinter$Stings > 2], SwitchDirn(dat_CairnsBeach_NoWinter$D3[dat_CairnsBeach_NoWinter$Stings > 2]), ">2 stings")
grid.arrange(S0, S1, S2, S3, ncol = 2)
dev.copy2pdf(file = "CairnsBeaches_NoWinter_StingsByWind2.pdf", paper = "A4r")

S0 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings == 0], SwitchDirn(dat_CairnsBeach$D10[dat_CairnsBeach$Stings == 0]), "No stings")
S1 <- PlotRose(dat_CairnsBeach$S1[dat_CairnsBeach$Stings >= 1], SwitchDirn(dat_CairnsBeach$D10[dat_CairnsBeach$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("CairnsIslandStingsByWind3.png", dpi = 600)

S0 <- PlotRose(dat_CairnsIsland$S1[dat_CairnsIsland$Stings == 0], SwitchDirn(dat_CairnsIsland$D1[dat_CairnsIsland$Stings == 0]), "No stings")
S1 <- PlotRose(dat_CairnsIsland$S1[dat_CairnsIsland$Stings >= 1], SwitchDirn(dat_CairnsIsland$D1[dat_CairnsIsland$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("CairnsIslandStingsByWind.png", dpi = 600)

S0 <- PlotRose(dat_CairnsIsland_NoWinter$S1[dat_CairnsIsland_NoWinter$Stings == 0], SwitchDirn(dat_CairnsIsland_NoWinter$D1[dat_CairnsIsland_NoWinter$Stings == 0]), "No stings")
S1 <- PlotRose(dat_CairnsIsland_NoWinter$S1[dat_CairnsIsland_NoWinter$Stings >= 1], SwitchDirn(dat_CairnsIsland_NoWinter$D1[dat_CairnsIsland_NoWinter$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)

S0 <- PlotRose(dat_CairnsReef$S1[dat_CairnsReef$Stings == 0], SwitchDirn(dat_CairnsReef$D1[dat_CairnsReef$Stings == 0]), "No stings")
S1 <- PlotRose(dat_CairnsReef$S1[dat_CairnsReef$Stings >= 1], SwitchDirn(dat_CairnsReef$D1[dat_CairnsReef$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("CairnsIslandStingsByWind.png", dpi = 600)

S0 <- PlotRose(dat_CairnsReef_NoWinter$S1[dat_CairnsReef_NoWinter$Stings == 0], SwitchDirn(dat_CairnsReef_NoWinter$D1[dat_CairnsReef_NoWinter$Stings == 0]), "No stings")
S1 <- PlotRose(dat_CairnsReef_NoWinter$S1[dat_CairnsReef_NoWinter$Stings >= 1], SwitchDirn(dat_CairnsReef_NoWinter$D1[dat_CairnsReef_NoWinter$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)

S0 <- PlotRose(dat_TownsvilleBeach$S1[dat_TownsvilleBeach$Stings == 0], SwitchDirn(dat_TownsvilleBeach$D1[dat_TownsvilleBeach$Stings == 0]), "No stings")
S1 <- PlotRose(dat_TownsvilleBeach$S1[dat_TownsvilleBeach$Stings >= 1], SwitchDirn(dat_TownsvilleBeach$D1[dat_TownsvilleBeach$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("TownsvilleBeachStingsByWind.png", dpi = 600)

S0 <- PlotRose(dat_TownsvilleIslandReef$S1[dat_TownsvilleIslandReef$Stings == 0], SwitchDirn(dat_TownsvilleIslandReef$D1[dat_TownsvilleIslandReef$Stings == 0]), "No stings")
S1 <- PlotRose(dat_TownsvilleIslandReef$S1[dat_TownsvilleIslandReef$Stings >= 1], SwitchDirn(dat_TownsvilleIslandReef$D1[dat_TownsvilleIslandReef$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("TownsvilleIslandReefStingsByWind.png", dpi = 600)


S0 <- PlotRose(dat_WhitsundaysBeach$S1[dat_WhitsundaysBeach$Stings == 0], SwitchDirn(dat_WhitsundaysBeach$D1[dat_WhitsundaysBeach$Stings == 0]), "No stings")
S1 <- PlotRose(dat_WhitsundaysBeach$S1[dat_WhitsundaysBeach$Stings >= 1], SwitchDirn(dat_WhitsundaysBeach$D1[dat_WhitsundaysBeach$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("WhitsundaysBeachStingsByWind.png", dpi = 600)

S0 <- PlotRose(dat_WhitsundaysIslandReef$S1[dat_WhitsundaysIslandReef$Stings == 0], SwitchDirn(dat_WhitsundaysIslandReef$D1[dat_WhitsundaysIslandReef$Stings == 0]), "No stings")
S1 <- PlotRose(dat_WhitsundaysIslandReef$S1[dat_WhitsundaysIslandReef$Stings >= 1], SwitchDirn(dat_WhitsundaysIslandReef$D1[dat_WhitsundaysIslandReef$Stings >= 1]), "Stings")
grid.arrange(S0, S1, ncol = 2)
ggsave("WhitsundaysIslandReefStingsByWind.png", dpi = 600)

### Is there a tide effect? Yes, but relatively weak, and in both Incoming and Height
with(dat_CairnsBeach[dat_CairnsBeach$Stings == 0,], table(Incoming))
with(dat_CairnsBeach[dat_CairnsBeach$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_CairnsBeach[dat_CairnsBeach$Stings == 0,], mean(Height))
with(dat_CairnsBeach[dat_CairnsBeach$Stings >= 1,], mean(Height)) # More on low tide
with(dat_CairnsBeach[dat_CairnsBeach$Stings == 0,], hist(Height, breaks = seq(-1.75, 1.75, 0.25)))
with(dat_CairnsBeach[dat_CairnsBeach$Stings >= 1,], hist(Height, breaks = seq(-1.75, 1.75, 0.25))) # More on low tide

with(dat_CairnsIsland[dat_CairnsIsland$Stings == 0,], table(Incoming))
with(dat_CairnsIsland[dat_CairnsIsland$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_CairnsIsland[dat_CairnsIsland$Stings == 0,], mean(Height))
with(dat_CairnsIsland[dat_CairnsIsland$Stings >= 1,], mean(Height)) # More on low tide
with(dat_CairnsIsland[dat_CairnsIsland$Stings == 0,], hist(Height, breaks = seq(-1.75, 1.75, 0.25)))
with(dat_CairnsIsland[dat_CairnsIsland$Stings >= 1,], hist(Height, breaks = seq(-1.75, 1.75, 0.25))) # More on low tide

with(dat_CairnsReef[dat_CairnsReef$Stings == 0,], table(Incoming))
with(dat_CairnsReef[dat_CairnsReef$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_CairnsReef[dat_CairnsReef$Stings == 0,], mean(Height))
with(dat_CairnsReef[dat_CairnsReef$Stings >= 1,], mean(Height)) # More on low tide
with(dat_CairnsReef[dat_CairnsReef$Stings == 0,], hist(Height, breaks = seq(-1.75, 1.75, 0.25)))
with(dat_CairnsReef[dat_CairnsReef$Stings >= 1,], hist(Height, breaks = seq(-1.75, 1.75, 0.25))) # More on low tide

with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings == 0,], table(Incoming))
with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings == 0,], mean(Height))
with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings >= 1,], mean(Height)) # More on low tide
with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings == 0,], hist(Height, breaks = seq(-2, 2, 0.25)))
with(dat_TownsvilleBeach[dat_TownsvilleBeach$Stings >= 1,], hist(Height, breaks = seq(-2, 2, 0.25))) # More on low tide

with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings == 0,], table(Incoming))
with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings == 0,], mean(Height))
with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings >= 1,], mean(Height)) # More on low tide
with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings == 0,], hist(Height, breaks = seq(-2, 2, 0.25)))
with(dat_TownsvilleIslandReef[dat_TownsvilleIslandReef$Stings >= 1,], hist(Height, breaks = seq(-2, 2, 0.25))) # More on low tide

with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings == 0,], table(Incoming))
with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings == 0,], mean(Height))
with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings >= 1,], mean(Height)) # More on low tide
with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings == 0,], hist(Height, breaks = seq(-2, 2, 0.25)))
with(dat_WhitsundaysBeach[dat_TownsvilleBeach$Stings >= 1,], hist(Height, breaks = seq(-2, 2, 0.25))) # More on low tide

with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings == 0,], table(Incoming))
with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings >= 1,], table(Incoming)) # More on incoming tide but weak
with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings == 0,], mean(Height))
with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings >= 1,], mean(Height)) # More on low tide
with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings == 0,], hist(Height, breaks = seq(-2, 2, 0.25)))
with(dat_WhitsundaysIslandReef[dat_WhitsundaysIslandReef$Stings >= 1,], hist(Height, breaks = seq(-2, 2, 0.25))) # More on low tide

### ENSO Analysis
ENSO <- read.csv("ENSOdata_Nov2018.csv")
ENSO <- gather(ENSO, 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', key = "Month", value = "SOI")
ENSO <- ENSO %>% mutate(Month = as.factor(Month)) %>% 
  mutate(Month = factor(Month, levels(Month)[c(5, 4, 8, 1, 9, 7, 6, 2, 12, 11, 10, 3)])) %>% 
  filter(Year >= 1985) %>%
  filter(Year <= 2016) %>%
  mutate(Yr = case_when(Month == "Jan" | Month == "Feb" | Month == "Mar" | Month == "Apr" | Month == "May" | Month == "Jun" ~ as.numeric(Year) - 1,
    TRUE ~ as.numeric(Year))) %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(Yr = as.factor(Yr)) %>%
  mutate(Mn = factor(Month, levels(Month)[c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)])) %>%
  arrange(Year, Month)

m3 <- lm(SOI ~ Yr + Mn, data = ENSO)
summary(m3)
plot(allEffects(m3, las = 2))
str(allEffects(m3))
YrEffect_SOI <- allEffects(m3)$Yr
MnEffect_SOI <- allEffects(m3)$Mn

# Does SOI influence jellies seasonally or year=to=year
graphics.off()
plot(MnEffect_SOI[[5]], MnEffect_Jelly[[5]])
plot(YrEffect_SOI[[5]], YrEffect_Jelly[[5]])
# No relationships
