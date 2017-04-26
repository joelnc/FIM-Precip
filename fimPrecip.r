rm(list=ls())
library(dataRetrieval)
library(zoo)
## 1) Lookup table for sample site code --> precip data source
## 2) Read FIM records
## 3) Match FIM records with precip TS
## 4) Compute metrics
## 5) SO this for opti tips

## Load the precip sites table
precipSites <- read.csv(, header=TRUE, sep=",",
                        stringsAsFactors=FALSE)



## This is coming in as UTC
## Get USGS Precip Data
##siteNumber <- "02146600"
siteNumber <- "351132080504145"

QParameterCd <- "00060"
QParameterCd <- "00045"
StartDate <- "2007-10-01"
##EndDate <- "2012-09-30"
Daily <- readNWISuv(siteNumber, QParameterCd,
                    startDate=StartDate, endDate="")

Daily <- Daily[order(Daily$dateTime),]

regTS <- data.frame(
    dateTime=seq.POSIXt(to=max(Daily$dateTime), by=300,
                  from=round(Daily$dateTime[1],"mins")-
                      (round(Daily$dateTime[1],"mins")$min%%5)*60),kp=1)
Daily <- merge(Daily, regTS, by="dateTime", all=TRUE)


## Load a WQD file and cut to FIM only
wqDF <- read.csv("LabReport.txt", sep=",", header=TRUE,
                 stringsAsFactors=FALSE)

fimDF <- wqDF[which(wqDF$Element=="ICS1.1"),]
fimDF$Coldate <- as.POSIXct(fimDF$Coldate, format="%m/%d/%Y %H:%M:%S")
rm(wqDF)

## Probably turn fimDF into a list of DFs by site...
## ... but for now just peel off one site to do the anaysls
siteDF <- fimDF[which(fimDF$Site=="MC29A1"),]
siteDF <- siteDF[,!(names(siteDF) %in% c("TransfDate", "SQL_ID",
                                         "ApprovedBy", "MilesDup",
                                         "QCDate", "CollectedBy",
                                         "Locdescr", "Acode",
                                         "DateApproved"))]

## Intialize rainfall totals, and "max 30 min precip in 6 hour prec. wind.
siteDF.DT <- data.frame(Coldate=unique(siteDF$Coldate),
                        Rain6=as.numeric(NA),
                        Rain12=as.numeric(NA),
                        Rain24=as.numeric(NA),
                        Rain36=as.numeric(NA),
                        Rain30min6hr=as.numeric(NA),
                        Rain30min12hr=as.numeric(NA),
                        DaysTenth=as.numeric(NA))


## Start with a loop over unique sampling dates
for (event in 1:nrow(siteDF.DT)) {

    ## index the rainfall record closest to sample dt
    atIndex <- which.min(
        abs(siteDF.DT$Coldate[event]-Daily$dateTime))
    ## Index rainfal record closest to minus 6, 12, 24 h
    minus6 <-  which.min(
        abs((siteDF.DT$Coldate[event]-(60*60*6))-Daily$dateTime))
    minus12 <-  which.min(
        abs((siteDF.DT$Coldate[event]-(60*60*12))-Daily$dateTime))
    minus24 <-  which.min(
        abs((siteDF.DT$Coldate[event]-(60*60*24))-Daily$dateTime))
    minus36 <-  which.min(
        abs((siteDF.DT$Coldate[event]-(60*60*36))-Daily$dateTime))

    if (abs(atIndex-minus6)>7) {
        ## Total over 6 hours
        siteDF.DT$Rain6[event] <-
            sum(Daily$X_00045_00000[minus6:atIndex], na.rm=TRUE)

        temp0 <- zoo(Daily$X_00045_00000[minus6:atIndex])
        temp <- rollapply(temp0, FUN=function(x) sum(x, na.rm=TRUE),
                          width=6)
        siteDF.DT$Rain30min6hr[event] <- max(temp,na.rm=TRUE)
        rm(temp0, temp)
    }
    if (abs(atIndex-minus12>12)) {
        siteDF.DT$Rain12[event] <-
            sum(Daily$X_00045_00000[minus12:atIndex], na.rm=TRUE)
        temp0 <- zoo(Daily$X_00045_00000[minus12:atIndex])
        temp <- rollapply(temp0, FUN=function(x) sum(x, na.rm=TRUE),
                          width=6)
        siteDF.DT$Rain30min12hr[event] <- max(temp,na.rm=TRUE)
        rm(temp0, temp)
    }
    if (atIndex!=minus24) {
        siteDF.DT$Rain24[event] <-
            sum(Daily$X_00045_00000[minus24:atIndex], na.rm=TRUE)
    }
    if (atIndex!=minus36) {
        siteDF.DT$Rain36[event] <-
            sum(Daily$X_00045_00000[minus24:atIndex], na.rm=TRUE)
    }

    ## Time since 0.10" in an hour
    temp <- 0
    tempSeries <- Daily[max((atIndex-(12*24*30)),1):atIndex,]
    ##browser()
    if (nrow(tempSeries)>11) {
        for (i in 1:(nrow(tempSeries)-12)) {
            fI <- nrow(tempSeries)-11-i
            tI <- nrow(tempSeries)+1-i
            ##browser()
            temp[i] <- sum(tempSeries$X_00045_00000[fI:tI],na.rm=TRUE)
            if (max(temp)>0.099) {
                siteDF.DT$DaysTenth[event] <- length(temp)/(12*24)
                break
            }
        }
    }


    rm(atIndex, minus6, minus12, minus24, minus36, temp, tempSeries)
}


## Merge precip totals back into df
siteDF <- merge(siteDF, siteDF.DT, by="Coldate", all=TRUE)

## Plot TSS vs....
graphics.off()
useTSS <- which(siteDF$Analyte=="Total Suspended Solids")
par(mfrow=c(3,3))
plot(siteDF$Coldate[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain6[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain12[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain24[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain36[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain30min6hr[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$Rain30min12hr[useTSS], siteDF$Result[useTSS], pch=16)
plot(siteDF$DaysTenth[useTSS], siteDF$Result[useTSS], pch=16)


##
dataSet <- Daily[10000:20300,c(1,4)]
plot(dataSet)



temp <- 0
for (i in 1:(nrow(dataSet)-12)) {
        fI <- nrow(dataSet)-11-i
        tI <- nrow(dataSet)+1-i
        temp[i] <- sum(dataSet$X_00045_00000[fI:tI],na.rm=TRUE)
        if (max(temp)>0.099)
            break
}

length(temp)/(12*24)



temp <- 0
while(max(temp)<0.1) {
    for (i in 1:10000) {
        fI <- nrow(dataSet)-11-i
        tI <- nrow(dataSet)+1-i
        temp[i] <- sum(dataSet$X_00045_00000[fI:tI],na.rm=TRUE)
    }
}


rollapply(data=dataSet$X_00045_00000, width=20, FUN=sum)
