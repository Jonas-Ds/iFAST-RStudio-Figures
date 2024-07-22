#### Set the directory of the data file and load dependancies####
dir <- "D:/User/iFAST_DATA"
setwd(dir)

library(ggplot2)
library(stats)
library(export)

#### loading and processing dataset####
data <- read.csv("data.csv") #load the csv file with the data 

data$Dye <- as.factor(data$Dye)      #change "Strain", "Dye" and "Time" to factors
data$Strain <- as.factor(data$Strain)
data$Time <- as.factor(data$Time)

##### making subset per dye #####
  #for each fluorogen make a separate dataset.
amber <- subset(data, Dye == "AMBER")

lime <- subset(data, Dye == "LIME")

coral <- subset(data, Dye == "CORAL")



### Figures#####
  ### defining colour palette ###
colours <- c("#05a8aa", "#bc412b") #this predefined list of colours will determine the colour assigned to the iFAST and the Wild-type strain


### Lime #####

levels(lime$Time) #check the order in which the levels of the time variable are loaded in
  
  # The for loop will calculate if there is a difference in average between the iFAST and wild-type strain. Afterwards, it will output the p-value of this calculation
for (i in levels(lime$Time)){
  set <- subset(lime, Time == i)
  result <- t.test(Fold_change ~ Strain,data=set)$p.value
  print(result)
} 
  # This part of the code calculates and stores each time point's max Fold change value. This will later be used when making the figures to display the significance
lime_max <- lime %>% 
  group_by(Time) %>%
  summarise(y_pos = max(Fold_change)+1) %>%
  pull(y_pos)
  # Following is the part of the code that makes and dispays the figure.
ggplot(lime, aes( x = Time, y = Fold_change, fill = Strain)) + # Here the dataset is loaded and the different variables are assigned to the correct axis.
      annotate("rect",xmin=1.5,xmax=5.5,ymin=-Inf,ymax=Inf,fill="#99ff94", alpha = 0.5)+ # This makes the green rectangle in the background.
      annotate('text', x= 2, y= 45, label = "Lime", fontface = 'bold')+ # Use this line to annotate the figure with the used fluorogen.
      geom_vline(xintercept = 1.5, linetype = 'dotted', color = 'black', linewidth = 0.75)+ # This line and the next are used to produce the dotted lines that demarcate the green rectangle 
      geom_vline(xintercept = 5.5, linetype = 'dotted', color = 'black', linewidth = 0.75)+
      geom_point(shape = 21, position = position_jitterdodge())+ # Here we plot all data points as dots 
      geom_boxplot(position = position_dodge(0), alpha  = 0.75, outlier.shape = 21, outlier.alpha = 1) + # Here we make a boxplot that is semi-transparent
      scale_x_discrete(name="Time (min)", limits = c('-5', '2', '5', '10', '15', '20', '25'), 
                       labels = c('-5' = '-5', '2' = '2','5' = '5','10' = '10','15' = '15','20' ='W1', '25' = 'W2'))+ # Use this line to adjust the labeling of the x-axis
      scale_fill_manual(values = colours)+ # 
      scale_colour_manual(values = colours)+
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = Strain, colour = Strain),
        linewidth = 0.75,
        position = position_dodge(0)
      )+
      ylab("Fold change fluorescence (A.U.)")+
      theme_classic(base_size = 15)+
      theme(legend.position = c(0.9,0.9),
            axis.text = element_text(face = 'bold', color = 'black'),
            axis.ticks = element_line(color = 'black'),
            axis.title = element_text(face = 'bold'))+
  annotate('text', 
           x =c('-5','2','5','10','15','20','25'),
           label = c("", "","***","***","***","***","***"),
           y = lime_max,
           size = 6)
 # Export the file with the preferred width, height and resolution
graph2png(file="Image.png", width = 6 , height = 4, dpi = 1200)

  
### Coral ####

for (i in levels(coral$Time)){
  set <- subset(coral, Time == i)
  result <- t.test(Fold_change ~ Strain,data=set)$p.value
  print(result)
} 

coral_max <- coral %>% 
  group_by(Time) %>%
  summarise(y_pos = max(Fold_change) + 50) %>%
  pull(y_pos)


ggplot(coral, aes( x = Time, y = Fold_change, fill = Strain)) +
      annotate("rect",xmin=1.5,xmax=5.5,ymin=-Inf,ymax=Inf,
               fill="coral", alpha = 0.5)+
      annotate('text', x= 2, y= 2100, label = "Coral", fontface = 'bold')+
      geom_vline(xintercept = 1.5, linetype = 'dotted', color = 'black', linewidth = 0.75)+
      geom_vline(xintercept = 5.5, linetype = 'dotted', color = 'black', linewidth = 0.75)+
      geom_point(shape = 21, position = position_jitterdodge())+
      geom_boxplot(position = position_dodge(0), alpha = 0.75, outlier.shape = 21, outlier.alpha = 1) +
  scale_x_discrete(name="Time (min)", limits = c('-5', '2', '5', '10', '15', '20', '25'), 
                   labels = c('-5' = '-5', '2' = '2','5' = '5','10' = '10','15' = '15','20' ='W1', '25' = 'W2'))+
      scale_fill_manual(values = colours)+
      scale_colour_manual(values = colours)+    
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = Strain, colour = Strain),
        size = 0.75,
        position = position_dodge(0)
      )+
      ylab("Fold change fluorescence (A.U.)")+
      theme_classic(base_size = 15)+
      theme(legend.position = c(0.9,0.9),
            axis.text = element_text(face = 'bold', color = 'black'),
            axis.ticks = element_line(color = 'black'),
            axis.title = element_text(face = 'bold'))+
  annotate('text', 
           x =c('-5','2','5','10','15','20','25'),
           label = c("","***","***","***","***","***",""),
           y = coral_max,
           size = 6)

graph2png(file="Coral reversibility combined.png", width = 6 , height = 4, dpi = 1200)
### Amber ####

for (i in levels(amber$Time)){
  set <- subset(amber, Time == i)
  result <- t.test(Fold_change ~ Strain,data=set)$p.value
  print(result)
} 

amber_max <- amber %>% 
  group_by(Time) %>%
  summarise(y_pos = max(Fold_change)+1) %>%
  pull(y_pos)

ggplot(amber, aes( x = Time, y = Fold_change, fill = Strain)) +
      annotate("rect",xmin=1.5,xmax=5.5,ymin=-Inf,ymax=Inf,
                fill="#fed186", alpha = 0.5)+
      annotate('text', x= 2.05, y= 45, label = "Amber", fontface = 'bold')+
      geom_vline(xintercept = 1.5, linetype = 'dotted', color = 'black', size = 0.75)+
      geom_vline(xintercept = 5.5, linetype = 'dotted', color = 'black', size = 0.75)+
      geom_point(shape = 21, position = position_jitterdodge())+ 
      geom_boxplot(position =  position_dodge(0), alpha = 0.75, outlier.shape = 21, outlier.alpha = 1) +
      scale_x_discrete(name="Time (min)", limits = c('-5', '2', '5', '10', '15', '20', '25'), 
                   labels = c('-5' = '-5', '2' = '2','5' = '5','10' = '10','15' = '15','20' ='W1', '25' = 'W2'))+
      scale_fill_manual(values = colours)+
      scale_colour_manual(values = colours)+   
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = Strain, colour = Strain),
        size  = 0.75,
        position = position_dodge(0)
      )+
      ylab("Fold change fluorescence (A.U.)")+
      theme_classic(base_size = 15)+
      theme(legend.position = c(0.9,0.9),
            axis.text = element_text(face = 'bold', color = 'black'),
            axis.ticks = element_line(color = 'black'),
            axis.title = element_text(face = 'bold'))+
      annotate('text', 
           x =c('-5','2','5','10','15','20','25'),
           label = c("","***","","***","***","***","***"),
           y = amber_max,
           size = 6)

graph2png(file="Amber reversibility combined.png", width = 6, height = 4, dpi = 1200)
