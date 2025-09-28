
#Clean memory log
rm(list = ls())

#Open file (open - 2. Data_for_R)
holocron<- read.table(file.choose(), header = TRUE)

#Check if ok
attach(holocron)
names(holocron)
head(holocron)

#Parameters
holocron$Spp<-as.factor(holocron$Spp)
holocron$I15N<-as.numeric(holocron$I15N)
holocron$I13C<-as.numeric(holocron$I13C)
holocron$Body_g<-as.numeric(holocron$Body_g)
holocron$B_mm<-as.numeric(holocron$B_mm)
holocron$Abd_mm<-as.numeric(holocron$Abd_mm)
holocron$Cth_mm<-as.numeric(holocron$Cth_mm)
holocron$TA_H1<-as.numeric(holocron$TA_H1)
holocron$TA_H2<-as.numeric(holocron$TA_H2)
holocron$TA_H3<-as.numeric(holocron$TA_H3)
holocron$TB_H1<-as.numeric(holocron$TB_H1)
holocron$TB_H2<-as.numeric(holocron$TB_H2)
holocron$TB_H3<-as.numeric(holocron$TB_H3)
holocron$TC_H1<-as.numeric(holocron$TC_H1)
holocron$TC_H2<-as.numeric(holocron$TC_H2)
holocron$TC_H3<-as.numeric(holocron$TC_H3)


# Load required packages
library(ecodist)
library(dplyr)  # for data manipulation
library(lmtest)  # for hypothesis testing
library(ggplot2)  # for data visualization
library(sandwich)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(gridExtra)
library(grid)
library(ggExtra)
library(patchwork)
library(tidyverse)
library(ggtext)
library(car)        #for anova
library(dplyr)      # For data manipulation
library(qgraph)
library(ape)
library(nlme)
library(MuMIn)
library(emmeans)




################################################################################
################################################################################
################################################################################
################################################################################
#                            TESTING THE VARIABLES                             #
################################################################################
################################################################################
################################################################################
################################################################################



# Because variance analysis assumes data independence, we assessed whether the
# trophic and thermal niches of individual species were spatially autocorrelated.
# To achieve this, we calculated Moran’s I statistic using the ape package
# (Paradis & Schliep 2019), which allowed us to assess the extent to which
# individuals with similar niches are spatially clustered. Specifically, we
# tested spatial autocorrelation for: (i) average temperature (the average of
# all recorded temperatures for each individual); (ii) temperature width
# (difference between the recorded minimum and maximum temperature); (iii) δ13C;
# and (iv) δ15N.  To assess whether the measured variables influenced the
# thermal and trophic niches, we then applied Generalized Least Squares (GLS)
# models using the nlme package (Pinheiro et al. 2017). Given the significant
# spatial autocorrelation detected for average temperature and δ13C, we
# incorporated proximity matrices into the models to account for spatial
# autocorrelation. These models were created with either average temperature,
# temperature width, δ13C, or δ15N as response variables, and species
# identity, body size and their interaction as predictors. If an interaction
# between predictors was found to be significant, we then tested for pairwise
# differences using Tukey-adjusted pairwise contrasts, using the
# emmeans package (Lenth 2025).




################################################################################
#                              Thermic Variables                               #
################################################################################



# Moran's Test 


# I'm using Moran's I Test in order to identify effects of spatial
# correlation in our data. Let us use our data and observe if values of isotopes
# are spatially autocorrelated. The first thing we need to
# do is to organize values of longitude ("long") and latitude ("lat")
# from our table.


geoloc<-cbind(holocron$long, holocron$lat)


# Then, let us produce a distance matrix (Euclidean) using the longitude
# and latitude values.

samples.dist <- as.matrix(dist(geoloc))




#Lets place the names on the dist matrix


# Get the names of our samples from the 'holocron' data frame
sample_names <- holocron$Spp

# Set the row and column names of the distance matrix
rownames(samples.dist) <- sample_names
colnames(samples.dist) <- sample_names



# Now, we have to divide "one by each distance value", creating an inverse
#distance matrix. Basically, instead of a "distance matrix", we'll end up
#with a "proximity matrix" (higher values indicate close sites). We'll also
#make all diagonal values equal zero, due to methodological reasons. As a 
#note, we may note that an inverse distance matrix is usually used in this
#type of spatial dependence analysis because it is better to quantify the
#idea of the first law of geography, that things that are close to each
#other are usually very similar.

samples.dist.inv <- 1/samples.dist
diag(samples.dist.inv) <- 0


# We can then test if the values of the isotopes are positively spatially autocorrelated
#(meaning that sites with similar values of isotopes are close to each
#other). To do that, let us use a Moran's test using the log of
#the isotope values and our inverse distance matrix. Note that the argument
#alternative="greater" means that we'll test if the value of Moran's I is
#greater than expected by chance (if data are more positively autocorrelated
#than expected by chance). We could, for example, use alternative="two.sided"
#if we wanted to use a two-tailed test (in order to test if spatial
#autocorrelation was either positive or negative).





################################################################################




#Average Temperature
Moran.I( holocron$T_ave , samples.dist.inv ,alternative="two.sided")


# $observed
# [1] 0.2361916
# 
# $expected
# [1] -0.02702703
# 
# $sd
# [1] 0.08667939
# 
# $p.value
# [1] 0.002391896


#SIGNIFICANT



#Temperature Width
Moran.I( holocron$T_wid , samples.dist.inv ,alternative="two.sided")


# $observed
# [1] 0.03634726
# 
# $expected
# [1] -0.02702703
# 
# $sd
# [1] 0.08726199
# 
# $p.value
# [1] 0.4676837



#NOT SIGNIFICANT






#   Accommodating spatial autocorrelation in the model      

# https://www.flutterbys.com.au/stats/tut/tut8.4a.html 






#T_ave



final.model.GLS.T_ave.Exp <- gls(T_ave ~ Spp*B_mm , data = holocron,
                                 correlation = corExp(form = ~lat + long, nugget = TRUE),
                                 method = "REML")

final.model.GLS.T_ave.Gaus <- gls(T_ave ~ Spp*B_mm , data = holocron,
                                  correlation = corGaus(form = ~lat + long, nugget = TRUE),
                                  method = "REML")


final.model.GLS.T_ave.Ratio <- gls(T_ave ~ Spp*B_mm , data = holocron,
                                   correlation = corRatio(form = ~lat + long, nugget = TRUE),
                                   method = "REML")


final.model.GLS.T_ave.Spher <- gls(T_ave ~ Spp*B_mm , data = holocron,
                                   correlation = corSpher(form = ~lat + long, nugget = TRUE),
                                   method = "REML")




AIC(final.model.GLS.T_ave.Exp, final.model.GLS.T_ave.Gaus, final.model.GLS.T_ave.Ratio, final.model.GLS.T_ave.Spher)


# df      AIC
# final.model.GLS.T_ave.Exp   11 124.7918
# final.model.GLS.T_ave.Gaus  11 122.6586
# final.model.GLS.T_ave.Ratio 11 123.5370
# final.model.GLS.T_ave.Spher 11 123.8451




#model summary
summary(final.model.GLS.T_ave.Gaus)

# Generalized least squares fit by REML
# Model: T_ave ~ Sp * B_mm 
# Data: holocron 
# AIC      BIC    logLik
# 122.6586 138.0718 -50.32931
# 
# Correlation Structure: Gaussian spatial correlation
# Formula: ~lat + long 
# Parameter estimate(s):
#   range       nugget 
# 1.072181e-05 7.808484e-02 
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)  23.426466  1.007700 23.247455  0.0000
# SpL_arg       1.606869  1.738994  0.924022  0.3628
# SpL_ery       7.468895  4.855660  1.538183  0.1345
# SpL_geo       6.096597  1.502390  4.057932  0.0003
# B_mm          0.379824  0.115535  3.287522  0.0026
# SpL_arg:B_mm -0.027372  0.279691 -0.097864  0.9227
# SpL_ery:B_mm -0.652408  0.347370 -1.878137  0.0701
# SpL_geo:B_mm -0.494856  0.192119 -2.575773  0.0152
# 
# Correlation: 
#   (Intr) SpL_rg SpL_ry SpL_ge B_mm   SpL_rg:B_ SpL_ry:B_
# SpL_arg      -0.587                                                
# SpL_ery      -0.208  0.122                                         
# SpL_geo      -0.621  0.293  0.129                                  
# B_mm         -0.946  0.554  0.196  0.606                           
# SpL_arg:B_mm  0.397 -0.945 -0.082 -0.154 -0.414                    
# SpL_ery:B_mm  0.315 -0.184 -0.983 -0.202 -0.333  0.138             
# SpL_geo:B_mm  0.534 -0.241 -0.111 -0.957 -0.582  0.142     0.194   
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -2.509544529 -0.731634789 -0.001598369  0.670956974  1.472467317 
# 
# Residual standard error: 1.020224 
# Degrees of freedom: 38 total; 30 residual



#ANOVA
anova.final.model.GLS.T_ave.Gaus <- Anova(final.model.GLS.T_ave.Gaus)


# Print the T_ave Spatial Corr ANOVA table
print(anova.final.model.GLS.T_ave.Gaus)

# Analysis of Deviance Table (Type II tests)
# 
# Response: T_ave
# Df   Chisq Pr(>Chisq)    
# Sp       3 35.4425  9.823e-08 ***
# B_mm     1  5.0311     0.0249 *  
# Sp:B_mm  3  8.7980     0.0321 *  
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1







#Contrasts

emm <- emmeans(final.model.GLS.T_ave.Gaus, ~ Spp | B_mm)

contrast(emm, method = "pairwise", adjust = "tukey")

contrast(emmeans(final.model.GLS.T_ave.Gaus, ~ Spp), method = "pairwise", adjust = "tukey")


#___________________________________________________________________________________


emm_trends <- emtrends(final.model.GLS.T_ave.Gaus, ~ Spp, var = "B_mm")
contrast(emm_trends, method = "pairwise", adjust = "tukey")

#Works

#__________________________________________________________________________________

mean_Bmm <- mean(holocron$B_mm, na.rm = TRUE)

emm <- emmeans(final.model.GLS.T_ave.Gaus, ~ Spp, at = list(B_mm = mean_Bmm))
contrast(emm, method = "pairwise", adjust = "tukey")

quantiles_Bmm <- quantile(holocron$B_mm, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

emm <- emmeans(final.model.GLS.T_ave.Gaus, ~ Spp, at = list(B_mm = quantiles_Bmm))
contrast(emm, method = "pairwise", adjust = "tukey")













################################################################################

#T_wid



final.model.GLS.TW <- gls(T_wid ~ Spp*B_mm , data = holocron, method = "REML")



#model summary
summary(final.model.GLS.TW)

# Generalized least squares fit by REML
# Model: T_wid ~ Spp * B_mm 
# Data: holocron 
# AIC      BIC    logLik
# 209.9796 222.5903 -95.98978
# 
# Coefficients:
#   Value Std.Error    t-value p-value
# (Intercept)   12.262757  4.376259  2.8021092  0.0088
# SppL_arg      -1.771893  7.264876 -0.2438985  0.8090
# SppL_ery      17.828784 19.333793  0.9221565  0.3638
# SppL_geo       8.366817  7.010692  1.1934367  0.2421
# B_mm           0.457407  0.504345  0.9069321  0.3717
# SppL_arg:B_mm  0.563758  1.166562  0.4832643  0.6324
# SppL_ery:B_mm -1.514809  1.382737 -1.0955146  0.2820
# SppL_geo:B_mm -0.476083  0.882193 -0.5396585  0.5934
# 
# Correlation: 
#   (Intr) SppL_rg SppL_ry SppL_g B_mm   SppL_rg:B_ SppL_ry:B_
# SppL_arg      -0.602                                                    
# SppL_ery      -0.226  0.136                                             
# SppL_geo      -0.624  0.376   0.141                                     
# B_mm          -0.965  0.581   0.218   0.602                             
# SppL_arg:B_mm  0.417 -0.948  -0.094  -0.260 -0.432                      
# SppL_ery:B_mm  0.352 -0.212  -0.982  -0.220 -0.365  0.158               
# SppL_geo:B_mm  0.552 -0.332  -0.125  -0.970 -0.572  0.247      0.209    
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.84731079 -0.57682134 -0.01650969  0.40801198  2.04371751 
# 
# Residual standard error: 4.145381 
# Degrees of freedom: 38 total; 30 residual


#ANOVA
anova.final.model.GLS.TW <- Anova(final.model.GLS.TW)


# Print the T_ave Spatial Corr ANOVA table
print(anova.final.model.GLS.TW)

# Analysis of Deviance Table (Type II tests)
# 
# Response: T_wid
# Df   Chisq Pr(>Chisq)  
# Spp       3 11.2530    0.01043 *
# B_mm      1  0.5706    0.45001  
# Spp:B_mm  3  1.8695    0.59993  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




################################################################################
#                                      IMAGE                                   #
################################################################################



#Argiope

A_arg <- subset(holocron, Spp == "A_arg")

#Leucauge

L_arg <- subset(holocron, Spp == "L_arg")

#Lycoside

L_ery <- subset(holocron, Spp == "L_ery")

#Latrodectus

L_geo <- subset(holocron, Spp == "L_geo")





# pivot the data to long format

#Argiope
data__A_arg <- pivot_longer(A_arg, cols = c("TA_H1", "TA_H2", "TA_H3",
                                            "TB_H1", "TB_H2", "TB_H3",
                                            "TC_H1", "TC_H2", "TC_H3"), names_to = "Time", values_to = "Temperature")


##Lycosa erythrognatha
data__L_ery <- pivot_longer(L_ery, cols = c("TA_H1", "TA_H2", "TA_H3",
                                            "TB_H1", "TB_H2", "TB_H3",
                                            "TC_H1", "TC_H2", "TC_H3"), names_to = "Time", values_to = "Temperature")



##Latrodectus geometricus
data__L_geo <- pivot_longer(L_geo, cols = c("TA_H1", "TA_H2", "TA_H3",
                                            "TB_H1", "TB_H2", "TB_H3",
                                            "TC_H1", "TC_H2", "TC_H3"), names_to = "Time", values_to = "Temperature")




##Leucauge argyra
data__L_arg <- pivot_longer(L_arg, cols = c("TA_H1", "TA_H2", "TA_H3",
                                            "TB_H1", "TB_H2", "TB_H3",
                                            "TC_H1", "TC_H2", "TC_H3"), names_to = "Time", values_to = "Temperature")



# define a custom order for Time
custom_order <- c("TA_H1", "TA_H2", "TA_H3",
                  "TB_H1", "TB_H2", "TB_H3",
                  "TC_H1", "TC_H2", "TC_H3")



# convert Time to a factor with the custom order
data__A_arg$Time <- factor(data__A_arg$Time, levels = custom_order)


# convert Time to a factor with the custom order
data__L_ery$Time <- factor(data__L_ery$Time, levels = custom_order)


# convert Time to a factor with the custom order
data__L_geo$Time <- factor(data__L_geo$Time, levels = custom_order)


# convert Time to a factor with the custom order
data__L_arg$Time <- factor(data__L_arg$Time, levels = custom_order)




# create a vector of new names for the points on the x-axis
new_names <- c("17h00", "12h00", "Day 3 | 07h00","17h00", "12h00", "Day 2 | 07h00","17h00", "12h00", "Day 1 | 07h00")



# combine the data into one data frame
combined_data <- rbind(data__A_arg, data__L_ery, data__L_geo, data__L_arg)
# set the limits for the color gradient
color_range <- c(0, 50)



# Plot1
A_arg_plot<-ggplot(data__A_arg, aes(x = Temperature, y = Time, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1.1, rel_min_height = 0.01) +
  #geom_point()+
  scale_fill_viridis(limits = color_range) +
  ggtitle("Argiope argentata")+
  theme_ipsum() +
  scale_y_discrete(labels = new_names)+
  scale_x_continuous(limits=c(0, 50))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x ="ºC")+
  theme(plot.title = element_text(color="black", size=25, face="italic"),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(color="black", size=15, face="plain"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)
  )


# Plot2
L_ery_plot<-ggplot(data__L_ery, aes(x = Temperature, y = Time, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1.1, rel_min_height = 0.01) +
  #geom_point()+
  scale_fill_viridis(limits = color_range) +
  ggtitle("Lycosa erythrognatha")+
  theme_ipsum() +
  scale_y_discrete(labels = new_names)+
  scale_x_continuous(limits=c(0, 50))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x ="ºC")+
  theme(plot.title = element_text(color="black", size=25, face="italic"),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(color="black", size=15, face="plain"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)
  )


# Plot3
L_geo_plot<-ggplot(data__L_geo, aes(x = Temperature, y = Time, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1.1, rel_min_height = 0.01) +
  #geom_point()+
  scale_fill_viridis(limits = color_range) +
  ggtitle("Latrodectus geometricus")+
  theme_ipsum() +
  scale_y_discrete(labels = new_names)+
  scale_x_continuous(limits=c(0, 50))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x ="ºC")+
  theme(plot.title = element_text(color="black", size=25, face="italic"),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(color="black", size=15, face="plain"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)
  )


# Plot4
L_arg_plot<-ggplot(data__L_arg, aes(x = Temperature, y = Time, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1.1, rel_min_height = 0.01) +
  #geom_point()+
  scale_fill_viridis(limits = color_range) +
  ggtitle("Leucauge argyra")+
  theme_ipsum() +
  scale_y_discrete(labels = new_names)+
  scale_x_continuous(limits=c(0, 50))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x ="ºC")+
  theme(plot.title = element_text(color="black", size=25, face="italic"),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(color="black", size=15, face="plain"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)
  )


image1 <- grid.arrange(A_arg_plot, L_ery_plot, L_geo_plot, L_arg_plot, ncol = 2)









################################################################################
#                              Trophic Variables                               #
################################################################################



#13C
Moran.I( holocron$I13C , samples.dist.inv ,alternative="two.sided")


# $observed
# [1] 0.2887762
# 
# $expected
# [1] -0.02702703
# 
# $sd
# [1] 0.08772295
# 
# $p.value
# [1] 0.0003182089


#SIGNIFICANT


# As a result, we'll get the following results: observed value of Moran's I, 
#expected value of Moran's I as expected by chance (under the null hypothesis),
#standard deviation of Moran's I under the null hyothesis, and p.value of the
#test of the null hypothesis against the alternative hypothesis
#(http://127.0.0.1:22465/library/ape/html/MoranI.html). Notably, a positive,
#statistically significant value indicates a positive spatial autocorrelation.




#I13C



final.model.GLS.13C.Exp <- gls(I13C ~ Spp*B_mm , data = holocron,
                               correlation = corExp(form = ~lat + long, nugget = TRUE),
                               method = "REML")

final.model.GLS.13C.Gaus <- gls(I13C ~ Spp*B_mm , data = holocron,
                                correlation = corGaus(form = ~lat + long, nugget = TRUE),
                                method = "REML")


final.model.GLS.13C.Ratio <- gls(I13C ~ Spp*B_mm , data = holocron,
                                 correlation = corRatio(form = ~lat + long, nugget = TRUE),
                                 method = "REML")


final.model.GLS.13C.Spher <- gls(I13C ~ Spp*B_mm , data = holocron,
                                 correlation = corSpher(form = ~lat + long, nugget = TRUE),
                                 method = "REML")




AIC(final.model.GLS.13C.Exp, final.model.GLS.13C.Gaus, final.model.GLS.13C.Ratio, final.model.GLS.13C.Spher)


# df      AIC
# final.model.GLS.13C.Exp   11 92.98009
# final.model.GLS.13C.Gaus  11 93.64951
# final.model.GLS.13C.Ratio 11 92.97622
# final.model.GLS.13C.Spher 11 96.43579




#model summary
summary(final.model.GLS.13C.Ratio)

# df      AIC
# final.model.GLS.13C.Exp   11 92.98009
# final.model.GLS.13C.Gaus  11 93.64951
# final.model.GLS.13C.Ratio 11 92.97622
# final.model.GLS.13C.Spher 11 96.43579
# > #model summary
#   > summary(final.model.GLS.13C.Ratio)
# Generalized least squares fit by REML
# Model: I13C ~ Spp * B_mm 
# Data: holocron 
# AIC      BIC    logLik
# 92.97622 108.3894 -35.48811
# 
# Correlation Structure: Rational quadratic spatial correlation
# Formula: ~lat + long 
# Parameter estimate(s):
#   range       nugget 
# 4.667116e-05 2.655832e-01 
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)   -23.568930  0.587446 -40.12099  0.0000
# SppL_ery       -2.653135  3.270828  -0.81115  0.4237
# SppL_arg       -0.013469  0.962850  -0.01399  0.9889
# SppA_arg        0.695327  0.875231   0.79445  0.4332
# B_mm           -0.212244  0.073994  -2.86839  0.0075
# SppL_ery:B_mm   0.245269  0.232931   1.05297  0.3008
# SppL_arg:B_mm   0.104880  0.149447   0.70179  0.4882
# SppA_arg:B_mm   0.058402  0.108891   0.53633  0.5957
# 
# Correlation: 
#   (Intr) SppL_ry SppL_rg SppA_r B_mm   SppL_ry:B_ SppL_rg:B_
# SppL_ery      -0.163                                                    
# SppL_arg      -0.602  0.111                                             
# SppA_arg      -0.564  0.037   0.262                                     
# B_mm          -0.881  0.154   0.479   0.634                             
# SppL_ery:B_mm  0.267 -0.978  -0.155  -0.145 -0.315                      
# SppL_arg:B_mm  0.507 -0.094  -0.959  -0.267 -0.447  0.146               
# SppA_arg:B_mm  0.572 -0.029  -0.235  -0.959 -0.716  0.162      0.264    
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -1.4501747 -0.1363536  0.2119646  0.7643512  2.1812928 
# 
# Residual standard error: 0.6902921 
# Degrees of freedom: 38 total; 30 residual



#ANOVA
anova.final.model.GLS.13C.Ratio <- Anova(final.model.GLS.13C.Ratio)


# Print the T_ave Spatial Corr ANOVA table
print(anova.final.model.GLS.13C.Ratio)

# Analysis of Deviance Table (Type II tests)
#
# Response: I13C
# Df   Chisq Pr(>Chisq)    
# Sp       3 23.5134  3.156e-05 ***
# B_mm     1 12.2150  0.0004741 ***
# Sp:B_mm  3  1.4744  0.6881971    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





################################################################################

#15N
# Moran.I( log(holocron$I15N+1) , samples.dist.inv ,alternative="greater")

Moran.I( holocron$I15N , samples.dist.inv ,alternative="two.sided")

# $observed
# [1] -0.03550253
# 
# $expected
# [1] -0.02702703
# 
# $sd
# [1] 0.08774368
# 
# $p.value
# [1] 0.9230489


#NOT SIGNIFICANT


final.model.GLS.15N <- gls(I15N ~ Spp*B_mm , data = holocron, method = "REML")



#model summary
summary(final.model.GLS.15N)

# Generalized least squares fit by REML
# Model: I15N ~ Sp * B_mm 
# Data: holocron 
# AIC      BIC    logLik
# 133.6373 146.2481 -57.81866
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)   11.281005  1.226082  9.200857  0.0000
# SpL_arg       -2.499933  2.035376 -1.228242  0.2289
# SpL_ery      -12.259681  5.416684 -2.263318  0.0310
# SpL_geo       -1.948452  1.964162 -0.992002  0.3291
# B_mm          -0.371178  0.141301 -2.626869  0.0134
# SpL_arg:B_mm   0.145292  0.326832  0.444547  0.6598
# SpL_ery:B_mm   0.881379  0.387397  2.275132  0.0302
# SpL_geo:B_mm   0.118277  0.247161  0.478541  0.6357
# 
# Correlation: 
#   (Intr) SpL_rg SpL_ry SpL_ge B_mm   SpL_rg:B_ SpL_ry:B_
# SpL_arg      -0.602                                                
# SpL_ery      -0.226  0.136                                         
# SpL_geo      -0.624  0.376  0.141                                  
# B_mm         -0.965  0.581  0.218  0.602                           
# SpL_arg:B_mm  0.417 -0.948 -0.094 -0.260 -0.432                    
# SpL_ery:B_mm  0.352 -0.212 -0.982 -0.220 -0.365  0.158             
# SpL_geo:B_mm  0.552 -0.332 -0.125 -0.970 -0.572  0.247     0.209   
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -2.410099085 -0.459540262 -0.007197411  0.582210936  1.941948875 
# 
# Residual standard error: 1.161397 
# Degrees of freedom: 38 total; 30 residual



#ANOVA
anova.final.model.GLS.15N <- Anova(final.model.GLS.15N)


# Print the T_ave Spatial Corr ANOVA table
print(anova.final.model.GLS.15N)

# Analysis of Deviance Table (Type II tests)
# 
# Response: I15N
# Df  Chisq Pr(>Chisq)  
# Sp       3 6.3159    0.09721 .
# B_mm     1 5.8599    0.01549 *
# Sp:B_mm  3 5.1840    0.15881  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




################################################################################
#                                      IMAGE                                   #
################################################################################




# define a custom order for Time

custom_order2 <- c("L_geo", "L_ery","L_arg", "A_arg")


# convert Time to a factor with the custom order
holocron$Spp <- factor(holocron$Spp, levels = custom_order2)




# create a vector of new names for the points on the x-axis
# new_names2 <- c("Argiope argentata", "Leucauge argyra", "Aglaoctenus lagotis", "Latrodectus geometricus")

new_names2 <- c("Latrodectus geometricus", "Aglaoctenus lagotis", "Leucauge argyra", "Argiope argentata")




# set the limits for the color gradient
color_range_ISO_C <- c(-30, -20)

color_range_ISO_N <- c(0, 20)


scale_fill_distiller(palette = "YlOrBr")


i13c<- ggplot(holocron, aes(x = I13C, y = Spp, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 0.8, rel_min_height = 0.01) +
  # geom_point()+
  scale_fill_distiller(palette = "YlOrBr", limits = color_range_ISO_C) +
  theme_ipsum() +
  scale_y_discrete(labels = new_names2)+
  scale_x_continuous(breaks = seq(-30, -10, 1))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x = expression(paste("\u03B4"^"13","C")))+
  theme(plot.title = element_blank(),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(size=40),
        axis.text.x = element_text(size=40),
        axis.text.y = element_text(color = "black", size=40, face="italic")
  )




i15n<-ggplot(holocron, aes(x = I15N, y = Spp, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 0.8, rel_min_height = 0.01) +
  # geom_point()+
  scale_fill_distiller(palette = "RdPu", limits = color_range_ISO_N) +
  theme_ipsum() +
  scale_y_discrete(labels = new_names2)+
  scale_x_continuous(breaks=seq(3,15,2))+
  removeGrid(y = TRUE, x = FALSE)+
  labs(x = expression(paste("\u03B4"^"15","N")))+
  theme(plot.title = element_text(color="black", size=40, face="bold"),
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_blank(),  #remove y axis title
        axis.title.x = element_text(size=40),
        axis.text.x = element_text(size=40),
        axis.text.y = element_text(color = "black", size=40, face="italic")
  )





image2 <- grid.arrange(i13c, i15n, ncol = 1)










################################################################################

#                         Niche Width Image                                    #

################################################################################



# Calculate temperature width and I13C range
species_summary <- holocron %>%
  group_by(Spp) %>%
  reframe(
    I13C_min = min(I13C),
    I13C_max = max(I13C),
    I15N_min = min(I15N),
    I15N_max = max(I15N),
    Temp_min = min(TA_H1, TA_H2, TA_H3, TB_H1, TB_H2, TB_H3, TC_H1, TC_H2, TC_H3),
    Temp_max = max(TA_H1, TA_H2, TA_H3, TB_H1, TB_H2, TB_H3, TC_H1, TC_H2, TC_H3),
    I13C_mean = mean(I13C),
    I15N_mean = mean(I15N),
    Temp_mean = mean(c(TA_H1, TA_H2, TA_H3, TB_H1, TB_H2, TB_H3, TC_H1, TC_H2, TC_H3))
  )




# Create a mapping of species to numbers
species_map <- c("A_arg" = "Argiope argentata", "L_arg" = "Leucauge argyra", "L_ery" = "Aglaoctenus lagotis", "L_geo" = "Latrodectus geometricus")

# Replace species names in species_summary
species_summary <- species_summary 

#Legend

# Plot horizontal and vertical lines to form a cross for each species
ggplot(species_summary) +
  # Horizontal line (I13C range)
  geom_segment(aes(x = I13C_min, xend = I13C_max, y = Temp_mean, yend = Temp_mean, color = Spp), size = 1.5) +
  # Vertical line (Temperature range)
  geom_segment(aes(x = I13C_mean, xend = I13C_mean, y = Temp_min, yend = Temp_max, color = Spp), size = 1.5) +
  # Add points for the mean values (average)
  geom_point(aes(x = I13C_mean, y = Temp_mean, color = Spp), size = 6) +
  # Define custom colors for each species
  scale_color_manual(
    breaks = c("A_arg", "L_arg", "L_ery", "L_geo"),   # <- desired order
    values = c(
      "L_geo" = "#E69F00",
      "L_ery" = "#D55E00", 
      "L_arg" = "#44AA99", 
      "A_arg" = "#F0E442"
    ),
    labels = c(
      expression(italic("Argiope argentata")),
      expression(italic("Leucauge argyra")),
      expression(italic("Aglaoctenus lagotis")),
      expression(italic("Latrodectus geometricus"))
    )
  ) +
  labs(x = "I13C", y = "Temperature", color = "Species") +
  labs(x = "\u03B413C", y = "ºC") +
  theme_ipsum() +
  removeGrid(y = FALSE, x = FALSE) +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"),
    legend.position = "bottom",
    legend.text.align = 0, 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    panel.spacing = unit(0.1, "lines"),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 30, vjust = -3),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 24, angle = 0, hjust = 1),
    legend.box.margin = margin(t = 40),   # space above legend
    legend.margin = margin(t = 10, b = 10)  # inner padding of legend box
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))




# ggplot(species_summary) +
#   # Horizontal line
#   geom_segment(aes(x = I15N_min, xend = I15N_max,
#                    y = Temp_mean, yend = Temp_mean, color = Spp),
#                linewidth = 1.2) +
#   # Vertical line
#   geom_segment(aes(x = I15N_mean, xend = I15N_mean,
#                    y = Temp_min, yend = Temp_max, color = Spp),
#                linewidth = 1.2) +
#   # Mean points
#   geom_point(aes(x = I15N_mean, y = Temp_mean, color = Spp),
#              size = 4, shape = 16) +
#   # Custom colors and legend order
#   scale_color_manual(
#     breaks = c("A_arg", "L_arg", "L_ery", "L_geo"),   # <- desired order
#     values = c(
#       "L_geo" = "#E69F00",
#       "L_ery" = "#D55E00", 
#       "L_arg" = "#44AA99", 
#       "A_arg" = "#F0E442"
#     ),
#     labels = c(
#       expression(italic("Argiope argentata")),
#       expression(italic("Leucauge argyra")),
#       expression(italic("Aglaoctenus lagotis")),
#       expression(italic("Latrodectus geometricus"))
#     )
#   ) +
#   labs(x = expression(delta^13 * "C"), y = "ºC", color = "Species") +
#   theme_minimal(base_size = 20) +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 18),
#     axis.title.x = element_text(size = 24, vjust = -2),
#     axis.title.y = element_text(size = 24, angle = 0, vjust = 0.5, hjust = 1),
#     axis.text = element_text(size = 18),
#     plot.title = element_text(size = 24, face = "bold"),
#     legend.box.margin = margin(t = 20)
#   ) +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))




# Plot horizontal and vertical lines to form a cross for each species
ggplot(species_summary) +
  # Horizontal line (I13C range)
  geom_segment(aes(x = I13C_min, xend = I13C_max, y = Temp_mean, yend = Temp_mean, color = Spp), size = 1.5) +
  # Vertical line (Temperature range)
  geom_segment(aes(x = I13C_mean, xend = I13C_mean, y = Temp_min, yend = Temp_max, color = Spp), size = 1.5) +
  # Add points for the mean values (average)
  geom_point(aes(x = I13C_mean, y = Temp_mean, color = Spp), size = 6) +
  # Define custom colors for each species
  scale_color_manual(
    values = c(
      "L_geo" = "#E69F00",
      "L_ery" = "#D55E00", 
      "L_arg" = "#44AA99", 
      "A_arg" = "#F0E442"
    ),
    labels = c(
      expression(italic("Latrodectus geometricus")),
      expression(italic("Aglaoctenus lagotis")),
      expression(italic("Leucauge argyra")),
      expression(italic("Argiope argentata"))
    )
  ) +
  labs(x = "I13C", y = "Temperature", color = "Species") +
  labs(x = expression(paste("\u03B4"^"13","C")), y = "ºC") +
  theme_ipsum() +
  removeGrid(y = FALSE, x = FALSE) +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"),
    legend.position = "none",
    legend.text.align = 0, 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    panel.spacing = unit(0.1, "lines"),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 30, vjust = -3),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 30, angle = 0, hjust = 1),
)







# Plot horizontal and vertical lines to form a cross for each species
ggplot(species_summary) +
  # Horizontal line (I15N range)
  geom_segment(aes(x = I15N_min, xend = I15N_max, y = Temp_mean, yend = Temp_mean, color = Spp), size = 1.5) +
  # Vertical line (Temperature range)
  geom_segment(aes(x = I15N_mean, xend = I15N_mean, y = Temp_min, yend = Temp_max, color = Spp), size = 1.5) +
  # Add points for the mean values (average)
  geom_point(aes(x = I15N_mean, y = Temp_mean, color = Spp), size = 6) +
  # Define custom colors for each species
  scale_color_manual(
    values = c(
      "L_geo" = "#E69F00",
      "L_ery" = "#D55E00", 
      "L_arg" = "#44AA99", 
      "A_arg" = "#F0E442"
    ),
    labels = c(
      expression(italic("Latrodectus geometricus")),
      expression(italic("Aglaoctenus lagotis")),
      expression(italic("Leucauge argyra")),
      expression(italic("Argiope argentata"))
    )
  ) +
  labs(x = "I15N", y = "Temperature", color = "Species") +
  labs(x = expression(paste("\u03B4"^"15","N")), y = "ºC") +
  theme_ipsum() +
  removeGrid(y = FALSE, x = FALSE) +
  theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"),
    legend.position = "none",
    legend.text.align = 0, 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    panel.spacing = unit(0.1, "lines"),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 30, vjust = -3),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 30, angle = 0, hjust = 1)
  )















################################################################################

#                              Mantel Test                                     #

################################################################################


# To explore the potential covariation between body temperature, stable isotope
# values, and spatial position, we conducted three different Mantel tests
# (accounting for each pairwise combination of variables) using the R package
# ecodist (Goslee & Urban 2007). Three Euclidean distance matrices were
# generated: a thermal matrix (based on individual temperature records), a
# trophic matrix (using δ13C and δ15N values for all individuals), and a
# spatial matrix (derived from the geographic coordinates of all individuals).
# The statistical significance of all Mantel tests was obtained by comparing
# the Euclidean distance matrices of each factor with 10,000 random permutations
# of the arrangement of the cells within these matrices, assessing significant
# differences from zero through Spearman correlations.


#Select the interest variables



#Transformation


# # Find the minimum value in the data
# min_val <- min(I13C)
# 
# # Shift all values to be non-negative
# holocron$I13C <- I13C - min_val
# 
# # Find the minimum value in the data
# min_val <- min(I15N)
# 
# # Shift all values to be non-negative
# holocron$I15N <- I15N - min_val


#trophic

trophic_13 <- holocron[, colnames(holocron) %in% c("I13C")]
trophic_15 <- holocron[, colnames(holocron) %in% c("I15N")]
trophic <- cbind(holocron$I13C, holocron$I15N)

#termic

thermic <- cbind(holocron$TA_H1, holocron$TA_H2, holocron$TA_H3, holocron$TB_H1,holocron$TB_H2, holocron$TB_H3,holocron$TC_H1, holocron$TC_H2, holocron$TC_H3)
thermic_ave <- holocron[, colnames(holocron) %in% c("T_ave")]

#Location (coordinates)

geo<-cbind(holocron$long, holocron$lat)




#dissimilarity index (euclidean)

#trophic
trophic.13 <- distance(trophic_13, "euclidean")
trophic.15 <- distance(trophic_15, "euclidean")
trophic.all <- distance(trophic, "euclidean")
#thermic
thermic.all <- distance(thermic, "euclidean")
thermic.ave <- distance(thermic_ave, "euclidean")
#distance
distance <- distance(geo, "euclidean")





# trophic vs location

################################################################################



# Mantel test: is the geographic distance between points
# related to the difference in the trophic signature between points?


mantel(trophic.all ~ distance, nperm = 100000, mrank = TRUE)
# mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
# 0.2133939  0.0022700  0.9977400  0.0022700  0.1576725  0.2714491 



# thermic vs location

################################################################################

# Mantel test: is the geographic distance between points
# related to the difference in the thermal signature between points?


mantel(distance ~ thermic.all, nperm = 100000, mrank = TRUE)
# mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
# -0.13867729  0.96575000  0.03426000  0.09729000 -0.17636111 -0.08834567 




################################################################################
#                                      IMAGE                                   #
################################################################################



# Flatten the matrices
elements1 <- as.vector(trophic.all)
elements2 <- as.vector(thermic.all)
elements3 <- as.vector(distance)


#join vectors in df
dis_mat <- data.frame(trophic = elements1,  thermic = elements2,  space = elements3)



#Create Individual Plots


plot2 <- ggplot(dis_mat, aes(x = space, y = trophic)) +
  geom_point(alpha = 1/4, color = "black", size=3) +
  geom_smooth(se = FALSE, method = lm, color= "black", linetype = 1, size=1)+
  # geom_rug(col="black",alpha=1/4, size=1.5)+
  labs(x = "Spatial distances",
       y = "Trophic distances") +
  theme_ipsum()+
  removeGrid(y = TRUE, x = TRUE)+
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(size=20, face="plain"),
        axis.text.y = element_text(size=20, face="plain"),
        axis.title.y = element_text(size=30, face="bold", vjust = +3),  
        axis.title.x = element_text(size=30, face="bold", vjust = -3)
  )



plot3 <- ggplot(dis_mat, aes(x = space, y = thermic)) +
  geom_point(alpha = 1/4, color = "black", size=3) +
  geom_smooth(se = FALSE, method = lm, color= "black", linetype = 1, size=1)+
  #geom_rug(col="black",alpha=1/4, size=1.5)+
  labs(x = "Spatial distances",
       y = "Thermic distances") +
  theme_ipsum()+
  removeGrid(y = TRUE, x = TRUE)+
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(size=20, face="plain"),
        axis.text.y = element_text(size=20, face="plain"),
        axis.title.y = element_text(size=30, face="bold", vjust = +3),  
        axis.title.x = element_text(size=30, face="bold", vjust = -3)
  )








ggplot(dis_mat, aes(x = space, y = trophic)) +
  # points + trend
  geom_point(alpha = 0.35, color = "grey30", size = 2.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  
  # labels
  labs(x = "Spatial distances", y = "Trophic distances") +
  
  # styling similar to your other plot
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.x = element_text(size = 24, face = "bold", vjust = -2),
    axis.title.y = element_text(size = 24, face = "bold",
                                angle = 0, hjust = 1, vjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )





# --- Space vs Trophic ---
ggplot(dis_mat, aes(x = space, y = trophic)) +
  geom_point(color = "#D55E00", size = 3, alpha = 0.6) +
  geom_smooth(se = FALSE, method = lm, color = "black", linewidth = 1.2) +
  labs(x = "Spatial distances", y = "Trophic distances") +
  scale_y_continuous(breaks = seq(1, 5, by = 1), limits = c(0, NA)) +  # <- removes 0 from axis labels
  theme_ipsum() +
  removeGrid(y = FALSE, x = FALSE) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 28, vjust = -3),
    axis.title.y = element_text(size = 28, angle = 0, hjust = 1)
  )


# --- Space vs Thermic ---
ggplot(dis_mat, aes(x = space, y = thermic)) +
  geom_point(color = "#238A8DFF", size = 3, alpha = 0.6) +
  geom_smooth(se = FALSE, method = lm, color = "black", linewidth = 1.2) +
  labs(x = "Spatial distances", y = "Thermal distances") +
  theme_ipsum() +
  removeGrid(y = FALSE, x = FALSE) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 28, vjust = -3),
    axis.title.y = element_text(size = 28, angle = 0, hjust = 1)
  )
