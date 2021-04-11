##########################################################################################################
#                                                                                                        #
#  ==================================================================================================    #
#  |                                        TIME SERIES PROJECT                                      |   #     
#  |Harmonized Index of Consumer Prices: Passenger Transport by Air for European Union (28 countries)|   #
#  ==================================================================================================    #
#                                                                                                        #
#   Author :  Henri Noel Kengne                                    Instructor: Daniel Felix Ahelegbey    #
#                                                                                                        #
#                                Institution: African School of Economics                                #
#                                                                                                        #
##########################################################################################################


# Clear the console and all the variable in the memory    
cat("\014")
rm(list=ls())

# Here we specify our working directory and sub-directories
work_path <- "C:/Users/PT WORLD/Desktop/My Projects/Time_series"
setwd(work_path)
dir.create("figures") # save all the figures here



# We start by loading all the necessary libraries for our time series analysis
library("quantmod")
library("TSA")
library("tseries")
library("forecast")
library("smooth")
library("xts")
library("lubridate")
library("fpp2")
library("ggplot2")
library("gridExtra")
library("astsa")
library("cowplot")
library("latex2exp")
library("urca")
library("caschrono")
library("stargazer")
library("xtable")
library("pmisc")
library('pastecs')


" The data that we chose for our time series analysis is available on the website of Federal Bank of ST.LOUIS
here is the link St. Louis Fed: Economic Data - FRED http://research.stlouisfed.org/fred2/ "

########################
# PRELIMINARY ANALYSIS #
########################

# Load and analyze the data
Data<-read.csv("C:/Users/PT WORLD/Desktop/Academic_projects/Time_series_projects/Data/HICP.csv")
sum(is.na(Data)) # missing value in the data? No
summary(Data)
class(Data)      # Is my dataset a time series dataset? No
HICP<-ts(Data[,2], start = c(2000,12), frequency = 12) # Now the data is transformed into a time series data
start(HICP)      # the start date of the series
end(HICP)        # the end date of the series
summary(HICP)    # summary of the time series data

#Time series preliminary analysis: Visualization of the behavior of the time series
autoplot(HICP)+geom_point() +
  ggtitle("HICP- Passenger Transport by Air for European Union (28 countries)")+
  theme(plot.title = element_text(size = 8),legend.title = element_blank(),legend.position = "bottom")+
ylab("Index 2015 = 100")+
  geom_smooth(method = "lm", se = FALSE, size = 0.2, colour = "blue")# This a time plot of the time series data.
# It displays a strong positive trend and
# possible seasonality. It is also clear 
# that the time series is not stationary.
#Augmented Dickey-Fuller Test

ggsave("./figures/autoplot.pdf", width = 12, height = 8, units = "cm")
adf.test(HICP)

#Seasonality
Seasonality1<-ggseasonplot(HICP)+
  ggtitle("Seasonal plot of HICP")+
  ylab("Index 2015 = 100") # This is to investigate the seasonality in the data. The graph clearly shows that there
# are seasonal patterns in the data.
ts_plot_season <- function(x = x) {
  season <- cycle(x)
  season.factor <- factor(season)
  ggplot() + 
    geom_boxplot(mapping = aes(x = season.factor,
                               y = x)) +
    labs(x = "months", y =  "Index 2015=100")
}
seasonality2<-ts_plot_season(HICP)  # When we look at the 
# boxplot we can see that there is a certain pattern. The HICP seem to be low between
# February and July. It start increasing from AUGUST to reach it highest value around 
# September and then drop. This happens every year. 

seasonality <- grid.arrange(Seasonality1, seasonality2, ncol = 2, widths = c(4, 6))
ggsave("./figures/seasonality.pdf", seasonality, width = 30, height = 10, units = "cm")


##################################
# STATIONALIZING THE TIME SERIES #
##################################

# Fist difference analysis
D1<-HICP%>%diff(1)
stat.desc(D1)
skewness(D1)
kurtosis(D1)

autoplot(D1)+ggtitle("Time series plot of the first-difference")+
  theme(plot.title = element_text(size = 8))+
  ylab(TeX("$(1-B)X_t$}$"))+xlab("time")
adf.test(D1)
ggsave("./figures/firstdiff.pdf", width = 12, height = 8, units = "cm")

ACF<-ggAcf(D1)+ggtitle("ACF of the first difference")+
  theme(plot.title = element_text(size = 8))
PACF<-ggPacf(D1)+ggtitle("PACF of the first difference")+
  theme(plot.title = element_text(size = 8))
Correlation_diff1 <- grid.arrange(ACF, PACF, nrow = 2 )
ggsave("./figures/acf_pacf_diff1.pdf", Correlation_diff1, width = 30, height = 20, units = "cm")

# Fist difference at seasonal lag
SD1<-D1%>%diff(12)
stat.desc(SD1)
skewness(SD1)
kurtosis(SD1)

autoplot(SD1)+ggtitle("Time series plot of the first-difference at seasonal lag")+
  theme(plot.title = element_text(size = 8))+
  ylab(TeX("$(1-B)(1-B^{12})X_t$}$"))+xlab("time")
ggsave("./figures/firstdiff_seasonal_lag.pdf", width = 12, height = 8, units = "cm")
adf.test(SD1)

ACF<-ggAcf(SD1)+ggtitle("ACF of the first difference at seasonal lag")+
  theme(plot.title = element_text(size = 8))
PACF<-ggPacf(SD1)+ggtitle("PACF of the first-difference at seasonal lag")+
  theme(plot.title = element_text(size = 8))
Correlation_diff_seasonal_lag <- grid.arrange(ACF, PACF, nrow = 2 )
ggsave("./figures/acf_pacf_diff_seasonal_lag.pdf", Correlation_diff_seasonal_lag,
       width = 30, height = 20, units = "cm")


#######################################
# MODEL IDENTIFICATION AND ESTIMATION #
#######################################

dir.create("Tables") # Store the tables here

fit1 <- Arima(HICP,order = c(0,1,0), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit1
arimatex(fit1, tex = T, fixed = NULL, pr = F)


fit2 <- Arima(HICP,order = c(0,1,1), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit2
arimatex(fit2, tex = T, fixed = NULL, pr = F)


fit3 <- Arima(HICP,order = c(1,1,0), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit3
arimatex(fit3, tex = T, fixed = NULL, pr = F)

fit4 <- Arima(HICP,order = c(1,1,1), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit4
arimatex(fit4, tex = T, fixed = NULL, pr = F)

fit5 <- Arima(HICP,order = c(2,1,0), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit5
arimatex(fit5, tex = T, fixed = NULL, pr = F)

fit6 <- Arima(HICP,order = c(2,1,1), seasonal = list(order=c(1,1,1),period=12)
              ,lambda = NULL)
fit6
arimatex(fit6, tex = T, fixed = NULL, pr = F)

######################
# DIGNOSTIC CHECKING #
######################

checkresiduals(fit3, test = "LB")
checkresiduals(fit5, test = "LB")


##############
# FORCASTING #
##############

x_tr <- window(HICP,end=2020) 
# fit <- Arima(x_tr,order = c(2,1,0), seasonal = list(order=c(1,1,1),period=12),lambda = NULL)
f_fit<-forecast(fit5)

forcast <- autoplot(x_tr, series="Data") + autolayer(fit5$fitted, series="ARIMA(2,1,0)(1,1,1)[12]") +
  autolayer(f_fit, series="Prediction") +
  xlab("Year") + ylab("Index 2015 = 100") + 
  ggtitle("HICP- Passenger Transport by Air for European Union (28 countries)") + 
  theme_bw()+
  theme(plot.title = element_text(size = 8),legend.title = element_blank(),legend.position = "bottom")
ggsave("./figures/forcast.pdf", width = 12, height = 8, units = "cm")
#-------------------------------------------END-------------------------------------------------

