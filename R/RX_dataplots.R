#########################################################################################################
###    Load GAMS results from cost minimization
#########################################################################################################


# this script loads the gams result file containing a complete time series of net load 
# (load - wind -pv - hydro) for all scenarios and the time series of hydro_generation, res_level,
# curtailment (fixed to zero) and hydro_spill

# check if ramping limit is OK < 
# check no generation march events

# validate reservoir level and hydro generation against real data from scen 1

options(scipen=999)
options(max.print=100)

####################################################################################################### 
#####                                       Load libraries and functions
####################################################################################################### 
library("tools")
source("MyLibraries.R")
library(readxl)
library(readr)
library(feather)
library(foreign)
library(openxlsx)
library("RColorBrewer")
library("ggpubr")
library("ggthemes")
library("cowplot")
library("tibbletime")
library("ggalt")
library("tidyquant")
library("fmsb")
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(broom)
library(purrr)
library(magrittr)
library(tibble)
library(stringr)



####################################################################################################### 
# Load GAMS result file

plt_folder <- "../Data/plots/"
wdth = 17.2

scov <- c(0,3610,7500)
scre <- c(0,880,1760,3610)
scan <- c(0,880,1760,3610)

HYDROCAPACITY <- 16155;
CF_HYDRO        = 0.85

minsimyear <- 1986
minyear <- 1986
maxyear <- 2010

minline <- (minyear - minsimyear)*8760 - 3500

readkey <- function()
{
  key <- scan(n=1)
}

calcextremevents <- function(hourdata)
{
  max_balancing_cap <- 2
  
  eetemp <- hourdata %>% 
    spread(variable,mwh) %>% 
    group_by(wf_elh, cp_elh, pr_elh, hy_rmp) %>% 
    mutate(
      extreme_load    = NET_LOAD + Electrolysis - Hydro - Thermal,
      extreme_load    = ifelse(extreme_load < max_balancing_cap, 0,extreme_load),
      balance_hydro   = Hydro,
      balance_thermal = Thermal,
      extreme_event  = ifelse(abs(extreme_load) < max_balancing_cap, 0,1 ),
      event_nr       = EventNr(extreme_event),
      event_len      = EventLen(extreme_event)) %>%
    ungroup()
  return(eetemp)
}
calcextremeventsoverview <- function(eedata)
{
  eetempoview <- eedata %>%
    filter(event_nr != 0) %>% 
    group_by(cp_elh, wf_elh, pr_elh, hy_rmp,
             event_nr,
             type = cut(extreme_load,breaks = c(-Inf,0,+Inf), labels=c("surplus","deficit"))) %>% 
    summarize(date         = first(date),
              extreme_max  = max(extreme_load),
              extreme_mean = mean(extreme_load),
              extreme_sum  = sum(extreme_load),
              extreme_min  = min(extreme_load),
              #max_ramp     = max(residual_load_ramp),
              event_len    = first(event_len)) %>% 
    ungroup()
  return(eetempoview)
}

calchourlyhydroramp <- function(hourdata)
{
  hourdata<-hourdata %>% as_tibble() 
  hydro <- hourdata %>%
    filter(variable %in% "Hydro") %>% 
    mutate(dh = (mwh - lag(mwh))) %>% 
    pull(dh)
  spill <- hourdata %>%
    filter(variable %in% "HydroSpill") %>% 
    mutate(dh = (mwh - lag(mwh))) %>% 
    pull(dh)
  hrmpt <- data.frame(hydro = hydro, spill = spill)
  hrmpt$date <- hourdata$date[hourdata$variable=="HydroSpill"]
  hrmpt$cp_elh <- hourdata[1,]$cp_elh
  hrmpt$wf_elh <- hourdata[1,]$wf_elh
  hrmpt$pr_elh <- hourdata[1,]$pr_elh
  hrmpt$hy_rmp <- hourdata[1,]$hy_rmp
  return(hrmpt)
}

calchourlyloaddur <- function(hourdata)
{
  variables <- hourdata$variable %>% unique()
  
  length <- unique(hourdata %>% group_by(cp_elh,pr_elh,wf_elh,hy_rmp,variable) %>% summarize(n=n()) %>% pull(n))
  
  ldctemp <- data.frame(x = seq(1,length))
  for (v in variables) {
    mwh <- hourdata %>%
      filter(variable == v) %>% 
      pull(mwh)
    mwh <- as.data.frame(mwh)
    mwh <- mwh[order(-mwh),]
    mwh <- data.frame(mwh)
    ldctemp[[v]] <- mwh$mwh
  }
  ldctemp$cp_elh <- hourdata[1,]$cp_elh
  ldctemp$wf_elh <- hourdata[1,]$wf_elh
  ldctemp$pr_elh <- hourdata[1,]$pr_elh
  ldctemp$hy_rmp <- hourdata[1,]$hy_rmp
  ldctemp %<>%
    gather(variable, value, -cp_elh, -wf_elh, -pr_elh, -hy_rmp, -x) %>% 
    transmute(x, cp_elh, wf_elh, pr_elh, hy_rmp,
              variable,
              mwh = as.numeric(as.character(value)))
  return(ldctemp)
}

res_folder <- "../Data/RES/"
rds_folder <- "../Data/RES/"
simulations <- setdiff(list.files(res_folder), list.dirs(res_folder, recursive = FALSE, full.names = FALSE))


hourly <- data.frame()
daily <- data.frame()
hrmp <- data.frame()
ldc <- data.frame()

for (sim in simulations) {
  simname <- file_path_sans_ext(sim)
  print(simname)
  
  tuples <- data.frame(x = 1)
  params <- strsplit(simname, "-")
  for (p in params[[1]]) {
    name <- strsplit(p, "\\.")[[1]][1]
    value <- as.integer(strsplit(p, "\\.")[[1]][2])
    tuples[name] <- value
  }
  
  #if (tuples$cp_elh %in% scan & tuples$pr_elh == 100 & tuples$hy_rmp > 1) {
  #  print("reading data")
  #}
  #else {
  #  print("skipping data")
  #  next
  #}
  
  #header <- read.table(str_c(res_folder,sim), nrows = 1, header = FALSE, sep =',', stringsAsFactors = FALSE)
  #temp <- read.table(str_c(res_folder,sim), skip = minline, header = FALSE, sep =',')
  header <- c(
    "model_stat", "solver_stat", "res_used", "scenario",
    "year", "hour", "total_cost", "NET_LOAD", 
    "Wind", "Nuclear", "Hydro", "Thermal",
    "Curtailment", "LossofLoad", "ReservoirLevel", "HydroSpill",
    "TH-B-CHP-IND", "TH-B-CHP-L", "TH-B-CHP-M", "TH-B-CHP-S",        
    "TH-BG-IC-CHP", "TH-W-CHP", "TH-R-CHP",
    "TH-NGCC-CHP-L", "TH-NGCC-CHP-M", "TH-NG-IC-CHP", "TH-OCGT",
    "TH-FST", "Electrolysis", "Export"
  )
  htemp <- read.table(str_c(res_folder,sim), header = T, sep =',')
  colnames(htemp) <- unlist(header)
  #temp[, -c(17:30)]
  
  htemp %<>%
    dplyr::select(-model_stat, -solver_stat, -res_used, -scenario, -total_cost, -Nuclear) %>%
    mutate(hour = str_sub(hour,2) %>% as.numeric(),
           year = str_sub(year,2)  %>% as.numeric(),
           date = hour2date(year, hour),
           Windmin = Wind - Curtailment)
  
  htemp %<>% filter(year <= maxyear)
  
  # check model and solver status -> 1 = OK
  mstat <- htemp %>% .$model_stat %>%  unique()
  sstat <- htemp %>% .$solver_stat %>%  unique()
  #  temp$scenario <- strsplit(simname,"-")[[1]][1]
  
  if(length(mstat) > 1) {
    cat(sprintf("model stat errors in simulation: %s, stats: %s\n", sim,mstat))
    readkey()
  }
  if(length(sstat) > 1) {
    cat(sprintf("solver stat errors in simulation: %s, stats: %s\n", sim,sstat))
    readkey()
  }
  
  params <- strsplit(simname, "-")
  for (p in params[[1]]) {
    name <- strsplit(p, "\\.")[[1]][1]
    value <- as.integer(strsplit(p, "\\.")[[1]][2])
    htemp[name] <- value
  }
  
  htemp %<>% 
    gather(variable, value, -date, -cp_elh, -wf_elh, -pr_elh, -hy_rmp, -year, -hour) %>% 
    transmute(cp_elh, wf_elh, pr_elh, hy_rmp,
              date,
              variable,
              mwh = as.numeric(as.character(value)))
  
  ### EXTREME EVENTS
  

  
  dtemp <- htemp %>% mutate(date = as.Date(date)) %>% 
    group_by(date, cp_elh, wf_elh, pr_elh, hy_rmp, variable) %>%
    summarize(sgwh=sum(mwh)/10^3, mmwh=mean(mwh), mxmwh = max(mwh), mnmwh = min(mwh))
  
  if (htemp$cp_elh %in% scan[-length(scan)] & htemp$pr_elh == 100 & htemp$hy_rmp > 1) {
    ldct <- calchourlyloaddur(htemp)
    hrmpt <- calchourlyhydroramp(htemp)
    hrmpt <- hrmpt %>% filter(yday(date)>1)
    
    ldc <- rbind(ldc, ldct)
    hrmp <- rbind(hrmp, hrmpt)
    
    remove(ldct)
    remove(hrmpt)
  }
  
  saveRDS(htemp,str_c(rds_folder,'hourly/',simname,'.RDA'))
  htemp %<>% filter(year(date) > 1986 & year(date) < 1990)
  
  hourly <- rbind(hourly, htemp)
  daily <- bind_rows(daily, dtemp)

  
  remove(htemp)
  remove(dtemp)
 
  gc()
}

#hBASE <- readRDS(str_c(rds_folder,'hBASE.RDA'))
#dBASE <- readRDS(str_c(rds_folder,'dBASE.RDA'))

#hourly <- rbind(hBASE, hourly)
#daily <- bind_rows(dBASE, daily)

saveRDS(hourly,str_c(rds_folder,'hourly.RDA'))
saveRDS(daily,str_c(rds_folder,'daily.RDA'))
saveRDS(ldc,str_c(rds_folder,'ldc.RDA'))
saveRDS(hrmp,str_c(rds_folder,'hrmp.RDA'))

# start here for read

hourly <- readRDS(str_c(rds_folder,'hourly.RDA'))
daily <- readRDS(str_c(rds_folder,'daily.RDA'))
ldc <- readRDS(str_c(rds_folder,'ldc.RDA'))
hrmp <- readRDS(str_c(rds_folder,'hrmp.RDA'))

weekly <- daily %>%
  group_by(date = cut(date, "week"), cp_elh, wf_elh, pr_elh, hy_rmp, variable) %>%
  summarize(sgwh  = sum(sgwh),
            mmwh  = mean(mmwh)) %>% 
  ungroup() %>% 
  mutate(date = as.Date(date))
saveRDS(weekly,str_c(rds_folder,'weekly.RDA'))

monthly <- daily %>%
  group_by(date = cut(date, "month"), cp_elh, wf_elh, pr_elh, hy_rmp, variable) %>%
  summarize(stwh  = sum(sgwh)/10^3,
            mmwh  = mean(mmwh)) %>% 
  ungroup() %>% 
  mutate(date = as.Date(date))
saveRDS(monthly,str_c(rds_folder,'monthly.RDA'))

annual <- daily %>% mutate(year = year(date)) %>% 
  group_by(year, cp_elh, wf_elh, pr_elh, hy_rmp, variable) %>%
  summarize(stwh  = sum(sgwh)/10^3,
            mmwh  = mean(mmwh))
saveRDS(annual,str_c(rds_folder,'annual.RDA'))

### FULL LOAD HOURS from daily ###

flh <- daily %>%
  group_by(cp_elh, wf_elh, pr_elh, hy_rmp,variable) %>%
  summarize(mean = mean(mmwh)) %>% 
  mutate(meanp = mean/cp_elh*100)

data_folder <- '../Data/HY/'
hinflow <- read.table(text = gsub(".h", "\t", readLines(str_c(data_folder,'TS_Hydro.csv'))))
colnames(hinflow) <- c('year', 'hour', 'mwh')

hinflow %<>%
  mutate(year = str_sub(year,2)  %>% as.numeric(),
  date = hour2date(year, hour))

hinflow <- left_join(hinflow, hincor, by = "date")

dinflow <- hinflow %>%
  mutate(date = as.Date(date)) %>% 
  group_by(date) %>%
  summarize(sgwh  = sum(mwh)/10^3,
            mmwh  = mean(mwh),
            mdmwh = median(mwh),
            sgwhcor  = sum(mwhcor)/10^3,
            mmwhcor  = mean(mwhcor),
            mdmwhcor = median(mwhcor))

winflow <- dinflow %>%
  group_by(date = cut(date, "week")) %>%
  summarize(sgwh  = sum(sgwh),
            mmwh  = mean(mmwh),
            mdmwh = median(mdmwh),
            sgwhcor  = sum(sgwhcor),
            mmwhcor  = mean(mmwhcor),
            mdmwhcor = median(mdmwhcor)) %>% 
  ungroup() %>% 
  mutate(date = as.Date(date))

minflow <- dinflow %>%
  group_by(date = cut(date, "month")) %>%
  summarize(stwh  = sum(sgwh)/10^3,
            mmwh  = mean(mmwh),
            mdmwh = median(mdmwh),
            stwhcor  = sum(sgwhcor)/10^3,
            mmwhcor  = mean(mmwhcor),
            mdmwhcor = median(mdmwhcor)) %>% 
  ungroup() %>% 
  mutate(date = as.Date(date))

yinflow <- dinflow %>%
  group_by(year = year(date)) %>%
  summarize(stwh  = sum(sgwh)/10^3,
            mmwh  = mean(mmwh),
            mdmwh = median(mdmwh),
            stwhcor  = sum(sgwhcor)/10^3,
            mmwhcor  = mean(mmwhcor),
            mdmwhcor = median(mdmwhcor)) %>% 
  ungroup()


p <- dinflow %>% 
  ggplot() +
  geom_line(aes(date,mmwh), color = "cornflowerblue") +
  geom_line(aes(date,mmwhcor), color = "firebrick4")
p




### FULL LOAD HOURS IN PERCENT
cpcols <- c("#0098e0", "#f5b587", "#b86377") #MID,MED,HIGH (blau, orange, rot)
hycols <- c("#4f3961", "#fc9d9d")
wfcols <- c("#0c7b93", "#c02739")
hywfcols <- c( "#527318", "#899857", "#142850", "#0c7b93") #LHF18, LHF160, HHF18, HHF160


cpflhdat <- data.frame(
  expand_grid(cp_elh = c(880,1760,3610),
              meanp = seq(0,100,by=1))
)
cpflhdat %<>% 
  mutate(meanopt = cp_elh*0.9) %>% 
  mutate(cpopt = meanopt/meanp*100)

flhdat <- flh %>% 
  filter(variable == "Electrolysis") %>% 
  filter(pr_elh == 100) %>% 
  filter(hy_rmp > 1) %>% 
  filter(cp_elh %in% c*880,1760,3610) %>% 
  mutate(meanopt = cp_elh*0.9) %>% 
  mutate(cpopt = meanopt/meanp*100)

p <- flh %>%
  filter(variable == "Electrolysis") %>% 
  filter(hy_rmp > 1) %>% 
  filter(cp_elh > 0 & cp_elh < 25000) %>% 
  filter(pr_elh == 100) %>% 
  mutate(hydro_ramp_scenario=ifelse(hy_rmp==2,"HHF","LHF")) %>% 
  mutate(meanopt = cp_elh*0.9) %>% 
  mutate(cpopt = meanopt/meanp*100) %>% 
  ggplot() +
  geom_line(aes(cp_elh, meanp, linetype = hydro_ramp_scenario, color=factor(wf_elh))) +
  geom_line(data = cpflhdat %>% filter(cp_elh ==880), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[1], size=0.3) +
  geom_vline(xintercept = 880, color = cpcols[1], size=0.3, linetype = 3) +
  geom_line(data = cpflhdat %>% filter(cp_elh ==1760), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[2], size=0.3) +
  geom_vline(xintercept = 1760, color = cpcols[2], size=0.3, linetype = 3) + 
  geom_line(data = cpflhdat %>% filter(cp_elh ==3610), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[3], size=0.3) +
  geom_vline(xintercept = 3610, color = cpcols[3], size=0.3, linetype = 3) +
  geom_text(label="LOW", y = 55, x = 1650, hjust = 0,
            color = cpcols[1], size = 3, family = "mono", fontface="plain") +
  geom_text(label="MID", y = 60, x = 2850, hjust = 0,
            color = cpcols[2], size = 3, family = "mono", fontface="plain") +
  geom_text(label="HIGH", y = 65, x = 5200, hjust = 0,
            color = cpcols[3], size = 3, family = "mono", fontface="plain") +
  scale_y_continuous(name = "Electrolysis full load operation (%)",
                     breaks = seq(0,100,by=2)) +
  scale_x_continuous(name = "Electrolysis capacity (MW)", breaks = seq(0,100000,by=2000), expand = c(0.0,0)) +
  scale_linetype_manual(name = "Hydro Flexibility", values = c("solid", "longdash", "dashed", "dotted")) +
  scale_color_manual(name = "Hydro Flexibility & Hydrogen Value", values = wfcols) +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"), legend.position = "bottom") +
  guides(col = guide_legend(nrow=1),
         linetype = guide_legend(nrow=1)) +
  coord_cartesian(ylim = c(56,100), xlim = c(0,19500))
p

ggsave(filename="PUB2020-H2_loadhours_DETAIL.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)


p <- flh %>%
  filter(variable == "Electrolysis") %>% 
  filter(hy_rmp > 1) %>% 
  filter(cp_elh > 0) %>% 
  filter(pr_elh == 100) %>% 
  mutate(hydro_ramp_scenario=ifelse(hy_rmp==2,"HHF","LHF")) %>% 
  mutate(meanopt = cp_elh*0.9) %>% 
  mutate(cpopt = meanopt/meanp*100) %>% 
  ggplot() +
  geom_line(aes(cp_elh, meanp, linetype = hydro_ramp_scenario, color=factor(wf_elh))) +
  geom_line(data = cpflhdat %>% filter(cp_elh ==880), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[1], size=0.3) +
  geom_vline(xintercept = 880, color = cpcols[1], size=0.3, linetype = 3) +
  geom_line(data = cpflhdat %>% filter(cp_elh ==1760), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[2], size=0.3) +
  geom_vline(xintercept = 1760, color = cpcols[2], size=0.3, linetype = 3) + 
  geom_line(data = cpflhdat %>% filter(cp_elh ==3610), aes(cpopt, meanp, group = factor(cp_elh)), color = cpcols[3], size=0.3) +
  geom_vline(xintercept = 3610, color = cpcols[3], size=0.3, linetype = 3) +
  geom_text(label="LOW", y = 55, x = 1650, hjust = 0,
            color = cpcols[1], size = 3, family = "mono", fontface="plain") +
  geom_text(label="MID", y = 60, x = 2850, hjust = 0,
            color = cpcols[2], size = 3, family = "mono", fontface="plain") +
  geom_text(label="HIGH", y = 65, x = 5200, hjust = 0,
            color = cpcols[3], size = 3, family = "mono", fontface="plain") +
  scale_y_continuous(name = "Electrolysis full load operation (%)",
                     breaks = seq(0,100,by=5)) +
  scale_x_continuous(name = "Electrolysis capacity (MW)", breaks = seq(0,100000,by=10000), expand = c(0.01,0)) +
  scale_linetype_manual(name = "Hydro Flexibility", values = c("solid", "longdash", "dashed", "dotted")) +
  scale_color_manual(name = "Hydro Flexibility & Hydrogen Value", values = wfcols) +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"), legend.position = "bottom") +
  guides(col = guide_legend(nrow=1),
         linetype = guide_legend(nrow=1)) +
  coord_cartesian(ylim = c(52,100), xlim = c(0,100000))
p

ggsave(filename="PUB2020-H2_loadhours.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)


### BOXPLOT H2 ###

p <- annual %>%
  filter(hy_rmp > 1) %>% ungroup(hy_rmp) %>% mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF")) %>% 
  filter(cp_elh %in% c(scre,scov,20000)) %>%
  filter(cp_elh > 0) %>% 
  filter(pr_elh == 100) %>%
  filter(variable == "Electrolysis") %>%
  mutate(hydrogen_value = factor(wf_elh)) %>% 
  ggplot() +
  geom_boxplot(aes(factor(cp_elh), mmwh/cp_elh*100,
                   group = paste(cp_elh,hydrogen_value),
                   fill = hydrogen_value,
                   color = hydrogen_value),
               size=0.3, outlier.size = 0.5) +
  facet_grid(hy_rmp~.) +
  scale_y_continuous(name = "Electrolysis full load operation (%)", breaks = seq(0,100,by=20)) +
  scale_x_discrete(name = "Electrolysis capacity (MW)") +
  scale_fill_manual(name = "Hydrogen value",
                    values = lighten(wfcols, amount = 0.3, method = "absolute"))+
  scale_color_manual(name = "Hydrogen value", values = wfcols) +
  theme(legend.position="bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"))
p

ggsave(filename="PUB2020-H2_FLH-BoxplotoverCP.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)


library(colorspace)
cfdark <- darken("cornflowerblue", amount = 0.6)

hyavail <- yinflow %>% 
  mutate(value=stwh/mean(stwh)) %>% 
  pull(value)

wnavail <- annual %>% 
  filter(cp_elh==0,wf_elh==18,pr_elh==100,hy_rmp==2) %>% 
  ungroup() %>% 
  filter(variable %in% c("Wind")) %>% 
  dplyr::select(year,variable,stwh) %>% 
  spread(variable,stwh) %>% 
  mutate("Wind availability"=Wind) %>% 
  dplyr::select(year,"Wind availability") %>% 
  gather(variable,value,-year) %>% 
  group_by(variable) %>% 
  mutate(value=value/mean(value)) %>% 
  pull(value)

hywnavail <- data.frame(
  hydro = hyavail,
  wind = wnavail
)

hywnavail %<>%
  gather() %>% 
  mutate(key = ifelse(key == "wind", "Wind availability",
                      ifelse(key == "hydro", "Hydro availability", "none")))

hywn <- c("dodgerblue4", "turquoise4")
#wind/hydro variability
p <- hywnavail %>% 
  ggplot(aes(x=key,y=value*100)) + 
  geom_boxplot(aes(fill=key,color=key)) + #, fatten = NULL) +
  ylab("Yearly availability relative to mean (%)") +
  xlab("") +
  scale_fill_manual(values = lighten(hywn,0.14,method="absolute")) +
  scale_color_manual(values = hywn) +
  theme_tq(base_size=12) +
  theme(panel.spacing = unit(0.4, "cm"), legend.position="none")
p

ggsave(filename="PUB2020-variabilityHYWND.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)

ldcsmall <- ldc %>%
  mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF")) %>% 
  filter(cp_elh %in% c(0,3610)) %>%
  filter((cp_elh == 0 & wf_elh == 160) | (cp_elh > 0))

ldcfuel <- ldcsmall %>% 
  filter(grepl("TH-", variable)) %>% 
  mutate(fuel = ifelse(variable == "TH-R-CHP", "Waste",
                       ifelse(variable == "TH-W-CHP", "Waste",
                              ifelse(variable == "TH-B-CHP-IND", "Biomass",
                                     ifelse(variable == "TH-B-CHP-L", "Biomass",
                                            ifelse(variable == "TH-B-CHP-M", "Biomass",
                                                   ifelse(variable == "TH-B-CHP-S", "Biomass",
                                                          ifelse(variable == "TH-NGCC-CHP-L", "NG",
                                                                 ifelse(variable == "TH-NGCC-CHP-M", "NG",
                                                                        ifelse(variable == "TH-NG-IC-CHP", "NG",
                                                                               ifelse(variable == "TH-BG-IC-CHP", "Biogas", "Oil (& NG)"
                                                                               )))))))))))

fuelcols <- c("goldenrod3",
              "olivedrab2",
              "yellow2",
              "tan",
              "gray50")

fuelvars <- c("Waste",
              "Biomass",
              "NG",
              "Biogas",
              "Oil (& NG)")

p <- ldcfuel %>% 
  #filter(x<5000) %>% 
  mutate(wf_elh = ifelse(cp_elh == 0, "-", wf_elh)) %>%
  mutate(wf_elh_f = factor(wf_elh, levels = c('-', '18', '160'))) %>%
  mutate(cp_elh = ifelse(cp_elh == 0, "NO", "HIGH")) %>%
  mutate(cp_elh_f = factor(cp_elh, levels = c('NO', 'HIGH'))) %>%
  ggplot() +
  geom_area(aes(x=x,y=mwh,fill=factor(fuel, levels = rev(fuelvars)))) +#, levels = rev(thvars))),color = "black", size=0.5) +
  facet_grid(cp_elh_f+wf_elh_f ~ hy_rmp, scales="fixed") +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Fuel type", values = rev(fuelcols)) +
  scale_y_continuous(name = "Thermal production (MW)") +
  scale_x_continuous(name = "hours", breaks = c(0,25000,50000,100000,200000), labels = c("0", "25k", "50k", "100k", "200k")) +
  theme(legend.position="bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"),
        panel.grid.minor.x = element_line(color="gray"))
p

ggsave(filename="PUB2020_LDC-THERMAL.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)


p <- ldc %>%
  filter(x < 500) %>% 
  filter(cp_elh %in% c(0,3610)) %>%
  mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF")) %>% 
  mutate(wf_elh = ifelse(cp_elh == 0, "-", wf_elh)) %>%
  mutate(wf_elh_f = factor(wf_elh, levels = c('-', '18', '160'))) %>%
  mutate(cp_elh = ifelse(cp_elh == 0, "NO", "HIGH")) %>%
  mutate(cp_elh_f = factor(cp_elh, levels = c('NO', 'HIGH'))) %>%
  filter(variable %in% c("LossofLoad")) %>% 
  ggplot() +
  #geom_vline(xintercept = 0, color = "darkgray") +
  geom_line(aes(x=x,y=mwh),size=0.8) +
  facet_grid(cp_elh_f+wf_elh_f ~ hy_rmp) +
  theme(legend.position = "bottom") +
  scale_y_continuous(name = "Backup Capacity (MW)") +
  scale_x_continuous(name = "hours", breaks = c(0,25,50,100,150,seq(200,800, by=100)), expand = c(0.0,0.0)) +
  theme(legend.position="bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"))
p

ggsave(filename="PUB2020_LDC-LOL-DETAIL.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)


### BOXPLOT HYDRO MIN & HYDRO RAMP ###

hflow <- daily %>% 
  filter(cp_elh %in% scan ) %>% 
  #filter(wf_elh % 18) %>% 
  filter(pr_elh == 100) %>% 
  filter(hy_rmp > 1) %>% 
  filter(variable %in% c("Hydro", "HydroSpill")) %>%
  dplyr::select(-sgwh, -mxmwh, -mnmwh) %>% 
  group_by(date, cp_elh, wf_elh, pr_elh, hy_rmp) %>% 
  spread(variable, mmwh) %>% 
  rename(hydro = Hydro) %>% 
  rename(spill = HydroSpill)
hflow$x <- "flow"
hflow$date <- NULL

hflow <- as.data.frame(hflow)

hrmp$x <- "ramp"

hydro <- rbind(hflow, hrmp %>% dplyr::select(-date))

hydro %<>% mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF"))

p <- hydro %>%
  filter(cp_elh %in% scov) %>% 
  ggplot() +
  geom_boxplot(data = . %>% filter(x == "flow"), aes(factor(cp_elh), hydro+spill, group = paste0(cp_elh,hy_rmp),
                                                     fill = factor(hy_rmp),
                                                     color = factor(hy_rmp)),
               size=0.3, outlier.size = 0.8) +
  geom_boxplot(data = . %>% filter(x == "ramp"), aes(factor(cp_elh), abs(hydro+spill)*3, group = paste0(cp_elh,hy_rmp),
                                                     fill = factor(hy_rmp),
                                                     color = factor(hy_rmp)),
               size=0.3, outlier.size = 0.4) +
  facet_wrap(x~.) +
  scale_y_continuous(name = 'Hydro flow (MW)',
                     breaks = seq(0,30000,by=3000),
                     sec.axis = sec_axis(trans = ~ ./3,
                                         name = "Hydro ramping (MW)")) +
  scale_x_discrete(name = "Electrolysis Capacity (MW)") +
  scale_fill_manual(name = "Hydro ramp scenario", values = lighten(hycols,0.4)) +
  scale_color_manual(name = "Hydro ramp scenario", values = hycols) +
  theme(legend.position="bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"))
p

ggsave(filename="PUB2020_Boxplot-HydroRamping-and-Flow.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)

p <- ldc %>%
  filter(x < 125000) %>% 
  filter(cp_elh %in% c(0,3610)) %>%
  filter((cp_elh == 0 & wf_elh == 160) | (cp_elh > 0)) %>% 
  mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF")) %>% 
  mutate(wf_elh = ifelse(cp_elh == 0, "-", wf_elh)) %>%
  mutate(wf_elh_f = factor(wf_elh, levels = c('-', '18', '160'))) %>% 
  mutate(cp_elh = ifelse(cp_elh == 0, "NO", "HIGH")) %>%
  mutate(cp_elh_f = factor(cp_elh, levels = c('NO', 'HIGH'))) %>%
  filter(variable %in% stackvars) %>% 
  ggplot() +
  geom_line(aes(x=x,y=mwh,color=factor(variable, levels = rev(stackvars)))) +
  facet_grid(cp_elh_f+wf_elh_f ~ hy_rmp) +
  theme(legend.position = "bottom") +
  scale_y_continuous(name = "Production (MW)", breaks = seq(0,20000,by=3000)) +
  scale_x_continuous(name = "hours",
                     breaks = c(0,10000,25000,50000,100000),
                     minor_breaks = c(5000,15000,20000,75000,125000)) +
  scale_color_manual(name = "Type", values = rev(stackcols), labels = rev(stacklabs)) +
  theme(legend.position="bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"),
        panel.grid.major.x = element_line(color="gray50"),
        panel.grid.minor.x = element_line(color="gray85"),
        axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(filename="PUB2020_LDC-Spill_CurtWind.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)

############################
# BARCHARTS (SPIDER REPLC) #
############################\

c <- c(0, 880, 1760, 3610)
w <- c(18,160)
h <- c(2,5)


h2flh <- flh %>%
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  filter(variable == "Electrolysis") %>%
  #replace(is.na(.), 0) %>%  
  arrange(cp_elh, hy_rmp, wf_elh) %>% 
  mutate(h2flh = meanp) %>% 
  dplyr::select(-mean,-variable,-meanp)

h2var <- annual %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  filter(variable == "Electrolysis") %>% 
  group_by(cp_elh, wf_elh, hy_rmp , variable) %>% 
  mutate(h2var = ifelse(cp_elh > 0, IQR(mmwh), NaN)) %>% 
  filter(year == 1986) %>%
  ungroup() %>% 
  arrange(cp_elh, hy_rmp, wf_elh) %>% 
  dplyr::select(-year,-variable,-stwh,-mmwh)

bkps <- annual %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  filter(variable == "LossofLoad") %>% 
  group_by(cp_elh, wf_elh, hy_rmp , variable) %>% 
  summarize(bkps = sum(stwh)) %>% 
  ungroup() %>% 
  arrange(cp_elh, hy_rmp, wf_elh) %>% 
  dplyr::select(-variable)

bkpc <- ldc %>%
  #mutate(hy_rmp = ifelse(hy_rmp == "HHF", 2, 5)) %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  filter(variable == "LossofLoad") %>% 
  group_by(cp_elh, wf_elh, hy_rmp , variable) %>% 
  slice(1) %>%
  ungroup() %>% 
  arrange(cp_elh, hy_rmp, wf_elh) %>% 
  mutate(bkpc = mwh) %>% 
  dplyr::select(-x, -variable, -mwh)

ths <- annual %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  filter(variable == "Thermal") %>% 
  group_by(cp_elh, wf_elh, hy_rmp , variable) %>% 
  summarize(ths = sum(stwh)) %>%
  ungroup() %>% 
  arrange(cp_elh, hy_rmp, wf_elh) %>% 
  dplyr::select(-variable)

hymrmp <- hrmp %>%
  dplyr::select(-date) %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  group_by(cp_elh, wf_elh, hy_rmp) %>% 
  summarize(hymrmp = mean(abs(hydro+spill))) %>%
  ungroup() %>% 
  arrange(cp_elh, hy_rmp, wf_elh)

hymin <- hflow %>% 
  filter(cp_elh %in% c) %>%
  filter(pr_elh == 100, wf_elh  %in% w, hy_rmp %in% h) %>%
  group_by(cp_elh, wf_elh, hy_rmp) %>% 
  summarize(hymin = quantile(hydro+spill)[2]) %>% 
  arrange(cp_elh, hy_rmp, wf_elh)

comptbl <- h2flh %>% left_join(h2var) %>% left_join(bkps) %>% left_join(bkpc) %>% 
  left_join(ths) %>% left_join(hymrmp) %>% left_join(hymin)

inttbl <- mutate_all(comptbl, function(x) as.integer(x))

write.csv(inttbl, paste0('comptbl.csv'))

normalize <- function(x) {
  return ((x - 0) / (max(x,na.rm = T ) - 0))
}

comp <- comptbl %>% 
  gather(variable, value, -cp_elh, -wf_elh, -pr_elh, -hy_rmp) %>% 
  group_by(variable) %>% 
  mutate(nval = normalize(value)) %>% 
  ungroup()

hywfcols <- c( "#388e3c", "#8bc34a", "#3c70a4", "#6bc5d2") #LHF18, LHF160, HHF18, HHF160

p <- comp %>%
  mutate(hy_rmp = ifelse(hy_rmp == 2, "HHF", "LHF")) %>%
  mutate(cp_elh = ifelse(cp_elh == 0, "NO",
                         ifelse(cp_elh == 880, "LOW",
                         ifelse(cp_elh == 1760, "MID","HIGH")))) %>%
  mutate(cp_elh_f = factor(cp_elh, levels = c('NO', "LOW", "MID", "HIGH"))) %>% 
  mutate(variable = if_else(variable == "ths", "Thermal Volume",
                     if_else(variable == "hymrmp", "μ Hydro Ramping",
                     if_else(variable == "hymin", "∧ Hydro Flow",
                     if_else(variable == "h2var", "var H2 production",
                     if_else(variable == "h2flh", "H2 full load production",
                     if_else(variable == "bkps", "Backup Volume",
                     if_else(variable == "bkpc", "Backup Capacity", "NONE")))))))) %>% 
  mutate(variable_f = factor(variable, levels = rev(c("H2 full load production",
                                                  "var H2 production",
                                                  "Thermal Volume",
                                                  "Backup Volume",
                                                  "Backup Capacity",
                                                  "∧ Hydro Flow",
                                                  "μ Hydro Ramping")))) %>% 
  ggplot() +
  geom_col(aes(variable_f, nval*100,
               fill = factor(paste0(hy_rmp," - ",wf_elh))),
           position = position_dodge2(padding=0.25), width=0.75) +
  coord_flip() +
  facet_wrap( ~ cp_elh_f, ncol=2) +
  scale_fill_manual(name = "", values = hywfcols) +
  scale_y_continuous(name = "", breaks = seq(0,100,by=20)) +
  scale_x_discrete(name = "") +
  theme(legend.position = "bottom") +
  theme_tq(base_size=10) +
  theme(panel.spacing = unit(0.4, "cm"),
        panel.grid.major.x = element_line(color="gray50"),
        panel.grid.minor.x = element_line(color="gray85"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p
ggsave(filename="PUB2020_BAR-Tradeoff.png", path = plt_folder, plot=p,
       width = wdth, height = wdth/4*3, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 600)
