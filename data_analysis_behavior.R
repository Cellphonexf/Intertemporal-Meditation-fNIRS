#Behavioral data analysis & plotting
#This script requires one data file:"Behavioral_data.xlsx"
#Programmed by Feng XIAO (2024.6.25)

###Preparation
##Load required packages for analysis
package_list <- c('dplyr','readxl','effsize','ggpubr','ggplot2','tidyr',
                  'patchwork','cowplot','lmtest','mediation')
lapply(package_list, require, character.only = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###Data import
pre_rawdata <- read_excel("Behavioral_data.xlsx", sheet = 'pretest', na = "---")
post_rawdata <- read_excel("Behavioral_data.xlsx", sheet = 'posttest', na = "---")
p_rawdata <- read_excel("Behavioral_data.xlsx", sheet = 'present', na = "---")
e_rawdata <- read_excel("Behavioral_data.xlsx", sheet = 'end', na = "---")

###Preprocessing
participant_table <- data.frame(Ntotal = nrow(pre_rawdata),
                        Nmale = dim(filter(pre_rawdata, Gender=="m"))[1],
                        Nfemale = dim(filter(pre_rawdata, Gender=="f"))[1],
                        Mage = mean(pre_rawdata$Age),
                        Meanage_m = mean(filter(pre_rawdata, Gender=="m")$Age),
                        Sdage_m = sd(filter(pre_rawdata, Gender=="m")$Age),
                        Meanage_f = mean(filter(pre_rawdata, Gender=="f")$Age),
                        Sdage_f = sd(filter(pre_rawdata, Gender=="f")$Age))
t.test(filter(pre_rawdata, Gender=="m")$Age,filter(pre_rawdata, Gender=="f")$Age,
       paired=FALSE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
##Pretest
MDD_pre <- which(pre_rawdata$MDDHistory=="yes") 
rd_pre <- pre_rawdata[-MDD_pre,] #excluding 3 participants with MDD
pre_table <- data.frame(Ntotal = nrow(rd_pre),
                        Nmale = dim(filter(rd_pre, Gender=="m"))[1],
                        Nfemale = dim(filter(rd_pre, Gender=="f"))[1],
                        Mage = mean(rd_pre$Age),
                        SDage = sd(rd_pre$Age),
                        Meditation_yes = dim(filter(rd_pre, MeditationHistory=="yes"))[1],
                        Meditation_no = dim(filter(rd_pre, MeditationHistory=="no"))[1],
                        Righthanded = dim(filter(rd_pre, StrongHand=="r"))[1],
                        Lefthanded = dim(filter(rd_pre, StrongHand=="l"))[1])
##Posttest
MDD_post <- which(post_rawdata$MDDHistory=="yes")
rd_post <- post_rawdata[-MDD_post,] #excluding 3 participants with MDD
post_table <- data.frame(Ntotal = nrow(rd_post),
                        Nmale = dim(filter(rd_post, Gender=="m"))[1],
                        Nfemale = dim(filter(rd_post, Gender=="f"))[1],
                        Mage = mean(rd_post$Age),
                        Meditation_yes = dim(filter(rd_post, MeditationHistory=="yes"))[1],
                        Meditation_no = dim(filter(rd_post, MeditationHistory=="no"))[1],
                        Righthanded = dim(filter(rd_post, StrongHand=="r"))[1],
                        Lefthanded = dim(filter(rd_post, StrongHand=="l"))[1])
##Posttest: mindfulness meditation
MDD_present <- which(p_rawdata$MDDHistory=="yes")
CheckMed_present <- which(p_rawdata$ManiCheck=="j") #excluding 5 participants specifying the wrong meditation type
rd_present <- p_rawdata[-c(MDD_present,CheckMed_present),] #excluding 3 participants with MDD and 5 specifying the wrong meditation type
mean(rd_present$Involvement) #63.42
sd(rd_present$Involvement) #25.00
##Posttest: intertemporal meditation
MDD_end <- which(e_rawdata$MDDHistory=="yes")
CheckMed_end <- which(e_rawdata$ManiCheck=="f") #excluding 5 participants specifying the wrong meditation type
rd_end <- e_rawdata[-c(MDD_end,CheckMed_end),] #excluding 3 participants with MDD and 5 specifying the wrong meditation type
mean(rd_end$Involvement) #70.74
sd(rd_end$Involvement) #20.20

###Comparisons: pretest vs. posttest
##Time perception
time_post_p <- data.frame(SubjectNo = (filter(rd_post, Condition=="EP"))$SubjectNo,
                          Time = (filter(rd_post, Condition=="EP"))$Time)
time_post_e <- data.frame(SubjectNo = (filter(rd_post, Condition=="PE"))$SubjectNo,
                          Time = (filter(rd_post, Condition=="PE"))$Time)
time_pre_p <- data.frame()
time_pre_e <- data.frame()
for (a in time_post_p$SubjectNo) {
  temp_time <- filter(rd_pre, SubjectNo==a)
  time_pre_p <- rbind(time_pre_p, temp_time)
}
for (b in time_post_e$SubjectNo) {
  temp_time <- filter(rd_pre, SubjectNo==b)
  time_pre_e <- rbind(time_pre_e, temp_time)
}
#Mindfulness meditation (pre vs. post):
t.test(time_pre_p$Time,time_post_p$Time,paired=TRUE,
       alternative=c("two.sided"),
       var.equal=FALSE,
       conf.level=0.95) #NS
mean(time_pre_p$Time) #56.61
sd(time_pre_p$Time) #24.56
mean(time_post_p$Time) #51.18
sd(time_post_p$Time) #29.97
#Intertemporal meditation (pre vs. post):
t.test(time_pre_e$Time,time_post_e$Time,paired=TRUE,
       alternative=c("two.sided"),
       var.equal=FALSE,
       conf.level=0.95) #p < .001 
cohen.d(time_pre_e$Time, time_post_e$Time) #effect size = 0.56
mean(time_pre_e$Time) #62.94
sd(time_pre_e$Time) #22.03
mean(time_post_e$Time) #50.23
sd(time_post_e$Time) #23.47
#Mindfulness vs. Intertemporal:
t.test(time_post_p$Time,time_post_e$Time,paired=FALSE,
       alternative=c("two.sided"),
       var.equal=FALSE,
       conf.level=0.95) #NS
mean(time_post_p$Time) #51.18
sd(time_post_p$Time) #29.97
mean(time_post_e$Time) #50.23
sd(time_post_e$Time) #23.47

###Comparisons: mindfulness vs. intertemporal meditation
##Meditation involvement
t.test(rd_present$Involvement,rd_end$Involvement,paired=TRUE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
cohen.d(rd_present$Involvement,rd_end$Involvement) #effect size = 0.32
##Thrill
p_thrill <- rd_present[!is.na(rd_present$Thrill), ]
e_thrill <- rd_end[!is.na(rd_end$Thrill), ]
merge_thrill <- merge(p_thrill, e_thrill, by = "SubjectNo")
t.test(merge_thrill$Thrill.x,merge_thrill$Thrill.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
mean(merge_thrill$Thrill.x) #29.43
sd(merge_thrill$Thrill.x) #24.50
mean(merge_thrill$Thrill.y) #24.36
sd(merge_thrill$Thrill.y) #23.50
##Peace
p_peace <- rd_present[!is.na(rd_present$Peace), ]
e_peace <- rd_end[!is.na(rd_end$Peace), ]
merge_peace <- merge(p_peace, e_peace, by = "SubjectNo")
t.test(merge_peace$Peace.x,merge_peace$Peace.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
mean(merge_peace$Peace.x) #59.75
sd(merge_peace$Peace.x) #27.16
mean(merge_peace$Peace.y) #62.55
sd(merge_peace$Peace.y) #22.89
##Anxiety
p_anxious <- rd_present[!is.na(rd_present$Anxious), ]
e_anxious <- rd_end[!is.na(rd_end$Anxious), ]
merge_anxious <- merge(p_anxious, e_anxious, by = "SubjectNo")
t.test(merge_anxious$Anxious.x,merge_anxious$Anxious.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
mean(merge_anxious$Anxious.x) #27.55
sd(merge_anxious$Anxious.x) #25.90
mean(merge_anxious$Anxious.y) #28.41
sd(merge_anxious$Anxious.y) #24.75
##Joy
p_joy <- rd_present[!is.na(rd_present$Joy), ]
e_joy <- rd_end[!is.na(rd_end$Joy), ]
merge_joy <- merge(p_joy, e_joy, by = "SubjectNo")
t.test(merge_joy$Joy.x,merge_joy$Joy.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #NS
mean(merge_joy$Joy.x) #42.48
sd(merge_joy$Joy.x) #25.39
mean(merge_joy$Joy.y) #35.48
sd(merge_joy$Joy.y) #27.45
##Fear
p_fear <- rd_present[!is.na(rd_present$Fear), ]
e_fear <- rd_end[!is.na(rd_end$Fear), ]
merge_fear <- merge(p_fear, e_fear, by = "SubjectNo")
t.test(merge_fear$Fear.x,merge_fear$Fear.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #p = .008
cohen.d(merge_fear$Fear.x,merge_fear$Fear.y) #effect size = 0.54
mean(merge_fear$Fear.x) #17.53
sd(merge_fear$Fear.x) #22.25
mean(merge_fear$Fear.y) #30.08
sd(merge_fear$Fear.y) #24.15
##Sadness
p_sad <- rd_present[!is.na(rd_present$Sad), ]
e_sad <- rd_end[!is.na(rd_end$Sad), ]
merge_sad <- merge(p_sad, e_sad, by = "SubjectNo")
t.test(merge_sad$Sad.x,merge_sad$Sad.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #p = .001
cohen.d(merge_sad$Sad.x,merge_sad$Sad.y) #effect size = 0.69
mean(merge_sad$Sad.x) #21.35
sd(merge_sad$Sad.x) #23.63
mean(merge_sad$Sad.y) #39.17
sd(merge_sad$Sad.y) #27.55
##Relaxation
p_relax <- rd_present[!is.na(rd_present$Relax), ]
e_relax <- rd_end[!is.na(rd_end$Relax), ]
merge_relax <- merge(p_relax, e_relax, by = "SubjectNo")
t.test(merge_relax$Relax.x,merge_relax$Relax.y,paired=FALSE,
       alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95) #p = .033
cohen.d(merge_relax$Relax.x,merge_relax$Relax.y) #effect size = 0.43
mean(merge_relax$Relax.x) #61.33
sd(merge_relax$Relax.x) #24.51
mean(merge_relax$Relax.y) #50.08
sd(merge_relax$Relax.y) #27.95

### Comparisons of 7 emotional ratings after each meditation practices
## Intertemporal meditation
emo_intertemporal <- data.frame(thrill = rd_end$Thrill, peace = rd_end$Peace, 
                                joy = rd_end$Joy, relax = rd_end$Relax, 
                                anxiety = rd_end$Anxious, fear = rd_end$Fear,
                                sadness = rd_end$Sad)
long_emoInt <- emo_intertemporal %>%
  pivot_longer(cols = everything(), names_to = "Emotion", values_to = "Rating") %>%
  drop_na() #convert to long format and delete missing values
emoInt_result <- aov(Rating ~ Emotion, data = long_emoInt)
summary(emoInt_result) #F = 18.11, p < .001, eta square = 0.22
postHoc_emoInt <- TukeyHSD(emoInt_result)
print(postHoc_emoInt) #peace > relaxation > sadness > joy > fear > anxiety > thrill
## Mindfulness meditation
emo_mindfulness <- data.frame(thrill = rd_present$Thrill, peace = rd_present$Peace, 
                                joy = rd_present$Joy, relax = rd_present$Relax, 
                                anxiety = rd_present$Anxious, fear = rd_present$Fear,
                                sadness = rd_present$Sad)
long_emoMind <- emo_mindfulness %>%
  pivot_longer(cols = everything(), names_to = "Emotion", values_to = "Rating") %>%
  drop_na() #convert to long format and delete missing values
emoMind_result <- aov(Rating ~ Emotion, data = long_emoMind)
summary(emoMind_result) #F = 28.39, p < .001, eta square = 0.30
postHoc_emoMind <- TukeyHSD(emoMind_result)
print(postHoc_emoMind) #relaxation > peace > joy > thrill > anxiety > sadness > fear

###Multiple regression analysis for time perception on emotional ratings
##Intertemporal meditation
pcdata_e <- merge(time_post_e, rd_end, by = 'SubjectNo')
model_e <- lm(cbind(Thrill, Peace, Anxious, Joy, Fear, Sad, Relax) ~ Time, data = pcdata_e)
summary(model_e) #time perception ~ fear: p = .035, R squared = 0.140
##Mindfulness meditation
pcdata_p <- merge(time_post_p, rd_present, by = 'SubjectNo')
model_p <- lm(cbind(Thrill, Peace, Anxious, Joy, Fear, Sad, Relax) ~ Time, data = pcdata_p)
summary(model_p) #NS

###Plotting
##Time perception
viso_time <- read_excel("Behavioral_data.xlsx", sheet = 'viso_time', na = "---")
viso_time$Test <- factor(viso_time$Test, levels = c('Before', 'After'))
viso_time$Meditation <- factor(viso_time$Meditation, levels = c('Intertemporal',
                                                                'Mindfulness'))
p_time <- ggplot(data = viso_time, aes(x = Test, y = Mean, fill = Meditation,
                             Test = factor(1), Meditation = factor(1))) +
  geom_line(aes(x = Test, y = Mean, group = Meditation, color = Meditation),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD / sqrt(N), ymax = Mean + SD / sqrt(N),
                    color = Meditation), width = 0.2, position = position_dodge(0)) +
  labs(x = 'Meditation practice', y = 'Rating') +
  scale_color_manual(values = c('Intertemporal' = '#B22222',
                                'Mindfulness' = '#4169E1')) +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 100), oob = scales::squish,
                     breaks = c(1, seq(25, 100, by = 25))) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none') +
  plot_layout(nrow=1) +
  plot_annotation(title = '(A) Time perception',
                  theme = theme(plot.title = element_text(size = 7,color = 'black',
                                                          face = 'bold')))
ggsave("time_pic.pdf", plot = p_time, width = 2, height = 2)
##Emotional ratings (contrast)
viso_emo <- read_excel("Behavioral_data.xlsx", sheet = 'viso_emo', na = "---")
viso_emo$Meditation <- factor(viso_emo$Meditation, levels = c('Intertemporal',
                                                              'Mindfulness'))
filtered_emo <- viso_emo[viso_emo$Emotion %in% c('Sadness', 'Fear', 'Relaxation'), ]
colors <- c('#4169E1','#FF7F50','#B22222')
p_emo <- ggplot(data = filtered_emo, aes(x = Meditation, y = Mean, fill = Emotion,
                             Meditation = factor(1), Emotion = factor(1))) +
  geom_line(aes(x = Meditation, y = Mean, group = Emotion, color = Emotion),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD / sqrt(N), ymax = Mean + SD / sqrt(N),
                    color = Emotion), width = 0.8, position = position_dodge(0)) +
  labs(x = 'Meditation type', y = 'Rating') +
  scale_color_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 100), oob = scales::squish,
                     breaks = c(1, seq(25, 100, by = 25))) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none') +
  plot_layout(nrow=1) +
  plot_annotation(title = '(D) Emotional ratings (contrast)',
                  theme = theme(plot.title = element_text(size = 7,color = 'black',
                                                          face = 'bold')))
ggsave("contrast_emo_pic.pdf", plot = p_emo, width = 2, height = 2)
##Emotional ratings (within group)
#Intertemporal meditation
viso_intertemporal <- subset(viso_emo, Meditation == 'Intertemporal') %>%
  arrange(desc(Mean))
viso_intertemporal$Emotion <- factor(viso_intertemporal$Emotion, levels = viso_intertemporal$Emotion)
p_emo_intertemporal <- ggplot(data = viso_intertemporal, aes(x = Emotion, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", fill = "#B22222") +
  geom_errorbar(aes(ymin = Mean - SD / sqrt(N), ymax = Mean + SD / sqrt(N)), width = 0.2, position = position_dodge(0.9)) +
  labs(x = 'Emotion', y = 'Rating') +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 100), oob = scales::squish,
                     breaks = c(1, seq(25, 100, by = 25))) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none') +
  plot_annotation(title = '(B) Emotional ratings (intertemporal)',
                  theme = theme(plot.title = element_text(size = 7,color = 'black',
                                                          face = 'bold')))
ggsave("intertemporal_emo_pic.pdf", plot = p_emo_intertemporal, width = 2, height = 2)
#Mindfulness meditation
viso_mindfulness <- subset(viso_emo, Meditation == 'Mindfulness') %>%
  arrange(desc(Mean))
viso_mindfulness$Emotion <- factor(viso_mindfulness$Emotion, levels = viso_mindfulness$Emotion)
p_emo_mindfulness <- ggplot(data = viso_mindfulness, aes(x = Emotion, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", fill = "#4169E1") +
  geom_errorbar(aes(ymin = Mean - SD / sqrt(N), ymax = Mean + SD / sqrt(N)), width = 0.2, position = position_dodge(0.9)) +
  labs(x = 'Emotion', y = 'Rating') +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 100), oob = scales::squish,
                     breaks = c(1, seq(25, 100, by = 25))) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none') +
  plot_annotation(title = '(C) Emotional ratings (mindfulness)',
                  theme = theme(plot.title = element_text(size = 7,color = 'black',
                                                          face = 'bold')))
ggsave("mindfulness_emo_pic.pdf", plot = p_emo_mindfulness, width = 2, height = 2)



