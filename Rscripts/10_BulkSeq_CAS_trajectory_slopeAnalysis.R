
##Here we will plot and quantify the trajectories of CAS incrase with age in both sexes, as well as split by sex
#First we load the output of the VISION score calcualtion
load('R_objects/BulkSeq_Aging_Scoretable.bin')


####We will make a coupt of the original score table. Then we duplicate the values from the CommonAgingSignatue (CAS)
#To keep it simple, we'll just assign this column the label 'score'
score_table <- BulkSeq_Aging_Scoretable 
score_table$score <- score_table[,"CommonAgingSignature"]
#First we will have to make sure that 'age'  is set as a numberic vector. These get can messed up with all the different data transformations
score_table$age <-  as.numeric(as.character(score_table$age))

#Then we plot the scores for each region along the age-axis. We also quanitfy the per-age group in a boxplot
#This creates Figure 2D and Figure S5A
(myplot <- ggplot(score_table, aes(x=age, y=score)) + geom_boxplot(notch = F, outlier.colour = NA,  mapping = aes(fill=as.factor(age))) + 
   geom_smooth(data =  score_table, mapping = aes(x=age, y=score), color='red', se = F) +
   geom_point(color='#bbbdbf', position = 'jitter') + facet_grid(~tissue) 
)

#Now we just focus on the trajectories, leaving out the faceting, boxplots and individual points
#This creates Figure 2E
(myplot <- ggplot(score_table, aes(x=age, y=score, color=tissue)) + 
    geom_smooth(data = score_table, mapping = aes(x=age, y=score), se = F,  method = 'loess') + scale_y_continuous(limits = c(-0.3,0.5),expand = c(0,0)) 
)

#Because loess models are challengnign to interpret, we'll recreate the same plot using linear models instead
#This creates Figure 2E
(myplot <- ggplot(score_table, aes(x=age, y=score, color=tissue)) + 
    geom_smooth(data = score_table, mapping = aes(x=age, y=score), se = F,  method = 'lm') + scale_y_continuous(limits = c(-0.3,0.6),expand = c(0,0))
)

#This visual inspection already gives us the impression that the CAS does not increase in every region to the same degree. Now we'll use a statistical framework to properly qunatify
#and test the differences in CAS slopes

#First, we will calculate a linear model that allows for an interaction between age and tissue (denoted with an)
m.interaction <- lm(score ~ age*tissue, data = score_table)
#We are going to setup a dataframe that contains the intercepts for each tissue for later inspection of the offset
coefficient_dataframe <- as.data.frame(m.interaction$coefficients)
colnames(coefficient_dataframe) <- 'coefficients'
coefficient_dataframe <- coefficient_dataframe[grep('^tissue', row.names(coefficient_dataframe), value = T),,drop=F]
coefficient_dataframe$tissue <- gsub('tissue', '',row.names(coefficient_dataframe))
colnames(coefficient_dataframe)[1] <- 'offset'
#This should look like this:
#           coefficients tissue
# tissuecer  0.001919534    cer
# tissueplx -0.007308881    plx
# tissuecc   0.019756150     cc
# tissueent  0.005956020    ent
# tissuehi   0.038349191     hi
# tissuehi2  0.008037825    hi2


#Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow tissue-to-tissue comparison
m.lst <- lstrends(m.interaction, "tissue", var="age")
#Let's look at that output and store it in a separte dataframe to hold the information about the slopse
slope_dataframe <- as.data.frame(m.lst)
#This should look like this. The slope is stored under the variable 'age.trend'
#     tissue   age.trend          SE  df    lower.CL    upper.CL
# 1      cp 0.015992965 0.001141669 817 0.013752015 0.018233915
# 2     cer 0.013659390 0.001113636 817 0.011473465 0.015845316
# 3     plx 0.016420592 0.001116214 817 0.014229607 0.018611578
# 4      cc 0.018581337 0.001133469 817 0.016356484 0.020806191
# 5     ent 0.005967812 0.001196206 817 0.003619813 0.008315810
#Create pairwise comparisons between slopes and test using the pais function and store that result
#This lets us know which slopes are statistically significantly different from one anothre
slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))

#Before we go deeper into the analysis and inspection of the differences between slopes, we want to examine if there is a relationship between a region's CAS increase and it's baseline
#In other words, we want to see if there is a  correlation between slope and offset
#For that, we take the coeefiecnet dataframe that we gnereate above together with the slope dataframe
coefficient_dataframe <- merge.data.frame(slope_dataframe, coefficient_dataframe, by = 'tissue')
#And test that with cor.test
cor.test(coefficient_dataframe$offset, coefficient_dataframe$age.trend)
#Doesn't look like it. 
#We will plot the actual values, which creates Figure 2F
(myplot <- ggplot(coefficient_dataframe, aes(x=offset, y=age.trend, color=tissue)) + geom_point()  + 
   geom_smooth(method = 'lm', color='black', se = F) + 
    scale_y_continuous(limits = c(0,0.02),expand = c(0,0))  + 
    scale_x_continuous(limits = c(-0.1,0.2),expand = c(0,0))
)

#This looks like we can trust the results. So we will now plot the results from lstrends as bargraphs with confidence intervals
#We'll order the rank of the tissues in ascending order along their CAS slope estimate
slope_dataframe$tissue <- factor(slope_dataframe$tissue, levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'tissue']),
                                 ordered = T)
#We will color the bars accoring to the score incrase ('age.trend')
#This will create Figure 2G and provide the information for 2H
(myplot <- ggplot(slope_dataframe, aes(x=tissue, y=age.trend,fill=age.trend)) +  
    geom_bar(position=position_dodge(), stat="identity", mapping = aes(fill=age.trend)) + xlab("") + ylab("CAS slope")  +
    scale_fill_gradient2(low=brewer.pal(11, "RdYlBu")[10], midpoint = quantile(slope_dataframe$age.trend)[3], mid=brewer.pal(11, "RdYlBu")[6], high = brewer.pal(11, "RdYlBu")[2])+
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_y_continuous(expand = c(0,0)) + 
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))
)

#We will save the slope and coefficient dataframe for better analysis by cell type in the future
save(slope_dataframe, coefficient_dataframe, file='R_objects/BulkSeq_CAS_slopes.bin')

##Now we want to turn the analysis a bit on its head and see if we can rank mice based on the CAS scores of its region
#For that we will use dplyr's summarise function to quantify the median CAS across all the regions
#So we'll establish a 'ranking' across mice based on that median
animal_order <- as.data.frame(score_table %>% dplyr::group_by(mouseID) %>% dplyr::summarise(med=median(score)))
animal_order[order(animal_order$med), 'mouseID']
#Now we turn mouseID into a vector with factors, where we use ordered factor levels 
score_table$mouseID <- factor(score_table$mouseID, levels = animal_order[order(animal_order$med), 'mouseID'], ordered = T)
#We plot this as boxplots, which creates Figure S5B
(myplot <- ggplot(score_table, aes(x=mouseID, y=score)) + geom_boxplot( outlier.colour = NA,aes(fill=as.factor(age)))+  geom_point(aes(color=tissue),  stroke=0, shape=16, size=0.8) +
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1)) 
)
#We can also plot the density of age groups along the animal oder to show how the consistence deteriorates over age
animal_order <- as.data.frame(score_table %>% dplyr::group_by(mouseID) %>% dplyr::summarise(med=median(score)))
animal_order$rank <- seq(1,59,1)
animal_order$age <- plyr::mapvalues(as.character(animal_order$mouseID), from = as.character(score_table$mouseID), to=as.character(score_table$age))
(myplot <- ggplot(animal_order, aes(x=rank, fill=age)) + geom_density(alpha=0.6, color=NA) 
)




###Now we repeat the analysis in which we'll split the groups into male and females to assess any potential differences
#For that purpose, we need to limit ourselves to the timeframe from 3 to 21 months (as we only have data for both sexes for these ages)
score_table <- BulkSeq_Aging_Scoretable 
score_table <- as.data.frame(score_table %>% filter(age %in% c(3, 12, 15, 18, 21)))
score_table$score <- score_table[,"CommonAgingSignature"]
#First we will have to make sure that 'age'  is set as a numberic vector. These get can messed up with all the different data transformations
score_table$age <-  as.numeric(as.character(score_table$age))

#Then we plot the scores for each region along the age-axis. We also quanitfy the per-age group in a boxplot
#This creates Figure 2D and Figure S5C
(myplot <- ggplot(score_table, aes(x=age, y=score, color=sex)) + geom_boxplot(notch = F, outlier.colour = NA,  mapping = aes(fill=as.factor(age))) + 
    geom_smooth(data =  score_table, mapping = aes(x=age, y=score, color=sex), se = F, method = 'lm') +
    facet_grid(~tissue) 
)


#Now,we'll investigate if there's any global difference between male and female sample independent of region identity. For that purpose, we create a model with an age-sex interaction
#and use the tissue identitvy as a covariate
m.interaction <- lm(score ~ age*sex+tissue, data = score_table)

#Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow male-to-female comparison
m.lst <- lstrends(m.interaction, "sex", var="age")
#Let's look at that output and store it in a separte dataframe to hold the information about the slope values
slope_dataframe <- as.data.frame(m.lst)
#This should look like this. The slope is stored under the variable 'age.trend'
#     sex  age.trend           SE  df    lower.CL   upper.CL
# 1 Female 0.01235224 0.0005205332 717 0.011330285 0.01337419
# 2   Male 0.01057583 0.0005162186 717 0.009562352 0.01158931ult
#This lets us know which slopes are statistically significantly different from one anothre
slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))
#This reveals that there is a significant difference present

#We plot the reult as bargraaph
slope_dataframe$sex <- factor(slope_dataframe$sex, levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'sex']),
                                 ordered = T)
#We will color the bars accoring to the score incrase ('age.trend')
#This will create Figure 2G and provide the information for 2H
(myplot <- ggplot(slope_dataframe, aes(x=sex, y=age.trend,fill=sex)) +  
    geom_bar(position=position_dodge(), stat="identity", mapping = aes(fill=sex)) + xlab("") + ylab("CAS slope")  +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_y_continuous(expand = c(0,0)) + 
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))
)

#We will save the slope and coefficient dataframe for better analysis by cell type in the future
save(slope_dataframe, coefficient_dataframe, file='R_objects/BulkSeq_CAS_slopes.bin')


#Now we can attempt to repeat the analysis for each individual region
#That is, we subset for a given region, establisha  linear model and calculate the slope for each sex and test for a potential difference
#It is important to note that the due to differences in base liine expression between regions, it is not recommended to compare the slopes across regions (i.e. male corpus callosum vs male hypothalamus)
slope_dataframe_list <- list()

for (region in unique(score_table$tissue)) {
  #Calculate linear models for the current tissue/region
  m.interaction <- lm(score ~ age*sex, data = score_table %>% filter(tissue == region) )
  #Build a dataframe that contains the intercepts for each sex
  #Obtain slopes and conduct pairwise testing
  m.lst <- lstrends(m.interaction, "sex", var=c("age"))
  slope_dataframe <- as.data.frame(m.lst)
  #Create pairwise comparisons between slopes and test
  slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))
   slope_dataframe$tissue <- region
  slope_dataframe[,'pval'] <- c(slope_pair_wiseTestResults$p.value, slope_pair_wiseTestResults$p.value)
  #Store the result in the above's list
  slope_dataframe_list[[region]] <- slope_dataframe
  
}
#Now combine the dataframes to enable side-by-side plotting
slope_dataframe_comb <- data.table::rbindlist(slope_dataframe_list)
slope_dataframe_comb <- as.data.frame(slope_dataframe_comb)


slope_dataframe_comb$tissue <- factor(slope_dataframe_comb$tissue, levels = unique(as.character(slope_dataframe_comb[order(slope_dataframe_comb$age.trend, decreasing = T),'tissue'])),
                                      ordered = T)

#Fnally, we plot these as side-by-side graphs
(myplot <- ggplot(slope_dataframe_comb, aes(x=sex, y=age.trend, fill=sex)) + 
    geom_bar(position=position_dodge(), stat="identity", mapping = aes(fill=sex)) + xlab("") +ylab("CAS slope") +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+ 
    scale_y_continuous(expand = c(0,0)) + facet_grid(~tissue) +
    
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90,hjust = 1)) 
)

