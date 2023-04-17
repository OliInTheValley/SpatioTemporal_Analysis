
##Here we will plot and quantify the effects of YMP or aDR on the CAS levels of aged mice
#First we load the output of the VISION score calcualtion
load('R_objects/BulkSeq_Rejuvenation_Scoretable.bin')


####We will make a copy of the original score table. Then we duplicate the values from the CommonAgingSignatue (CAS)
#To keep it simple, we'll just assign this column the label 'score'
score_table <- BulkSeq_Rejuvenation_Scoretable
score_table$score <- score_table[,"CommonAgingSignature"]
#First set 'treatment' as a vector of factors with ordered levels so that the controls groups are positioned left of their corresponding rejuvenation intervention
score_table$treatment <- droplevels(factor(score_table$treatment, levels = c('PBS', 'YMP', 'AL', 'aDR'),ordered = T))

#We'll also assign which of the treatments is control and which is the actual intervention
score_table$treatmentClass <- plyr::mapvalues(score_table$treatment, 
                                          from =  c('PBS','YMP','AL','aDR'),
                                          to = c('Ctrl','Rejuvenation','Ctrl','Rejuvenation'))
score_table$treatmentClass <- droplevels(factor(score_table$treatmentClass, levels = c('Ctrl', 'Rejuvenation'),ordered = T))


#Then we plot the scores for each region split by treatment and intervention-type. 
#This creates Figure 7D/H and Figure S15K/M
(myplot <- ggplot(score_table, aes(x=treatmentClass, y=score)) + geom_boxplot(notch = F, outlier.colour = NA,  mapping = aes(fill=treatmentClass)) + 
    geom_point(color='#bbbdbf', position = 'jitter') + facet_grid(experiment_group~tissue, scales = "free_x", as.table = F, drop = T) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_blank())
)

#It seems that there are some regions where the rejuvenation intervnetion causes the CAS to drop. 
#Since we are looking at categorical variables, we have little choice but assess each pair of aDR/AL or YMP/PBS per region and perform multiple testin
#rstix' pairwise_t_test function is ideal to iterate a test across multiple regions
#We filter for the PBS/YMP samples - as we do not plan to compare between interventions (they are also animals from different sexes etc.)
#Further, since we're primarily interested in regions that reduce the CAS (i.e. "rejuvenation"), we can use a one-sided test to increase the power. However, a two-sided test yields similar results
score_table %>% filter(treatment %in% c('PBS', 'YMP')) %>% group_by(tissue) %>%
  pairwise_t_test(
    score ~ treatment, paired = F,
    p.adjust.method = "BH", alternative='less'
  )
#Interestingly, YMP causes a significant reduction of the CAS in multiple regions

#Let's check if that's the same sitution for the DR-treated mice
score_table %>% filter(treatment %in% c('AL', 'aDR')) %>% group_by(tissue) %>%
  pairwise_t_test(
    score ~ treatment, paired = F,
    p.adjust.method = "BH", alternative='less'
  )
#Doesn't seem to be the case! It looks as if DR is doing something  that's orthogonal to aging. 


#However, we had identified a set of genes that were differentially regulated under aDR, and it seemed that this was happening across regions. 
#We should inspect that and see if there's something similarly happening under YMP treatment.
#We can work with the exact same analysis, we'd just have to switch out the 'score' that we've set at the beginning of the analysis
score_table <- BulkSeq_Rejuvenation_Scoretable
score_table$score <- score_table[,"aDRSignature"]
#First set 'treatment' as a vector of factors with ordered levels so that the controls groups are positioned left of their corresponding rejuvenation intervention
score_table$treatment <- droplevels(factor(score_table$treatment, levels = c('PBS', 'YMP', 'AL', 'aDR'),ordered = T))

#We'll also assign which of the treatments is control and which is the actual intervention
score_table$treatmentClass <- plyr::mapvalues(score_table$treatment, 
                                              from =  c('PBS','YMP','AL','aDR'),
                                              to = c('Ctrl','Rejuvenation','Ctrl','Rejuvenation'))
score_table$treatmentClass <- droplevels(factor(score_table$treatmentClass, levels = c('Ctrl', 'Rejuvenation'),ordered = T))


#Then we plot the scores for each region split by treatment and intervention-type. 
#This creates Figure S16A/B
(myplot <- ggplot(score_table, aes(x=treatmentClass, y=score)) + geom_boxplot(notch = F, outlier.colour = NA,  mapping = aes(fill=treatmentClass)) + 
    geom_point(color='#bbbdbf', position = 'jitter') + facet_grid(experiment_group~tissue, scales = "free_x", as.table = F, drop = T) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_blank())
)

#If we run the pairwise test for this score, we confirm that the observations made via the boxplots are statistically significant - significant effects in all aDR tissue, but no real effect with YMP
score_table %>% filter(treatment %in% c('PBS', 'YMP')) %>% group_by(tissue) %>%
  pairwise_t_test(
    score ~ treatment, paired = F,
    p.adjust.method = "BH"
  )

score_table %>% filter(treatment %in% c('AL', 'aDR')) %>% group_by(tissue) %>%
  pairwise_t_test(
    score ~ treatment, paired = F,
    p.adjust.method = "BH"
  )
