#### Setup ####

# Load the tidyverse package
## Usually it's a good idea to have this at the top of your code, so you and your 
## collaborators know which packages are needed to run the code
library(tidyverse)

# Set your working directory
## This can be used as an alternative to using "R projects".
## Note: the "~" symbol means "home directory", which is variable depending on your 
## username and operating system (Mac or Windows or Linux).
## You can use the `getwd()` command to see what your current working directory is 
setwd("~/Course_Materials/Day1PM-2_R_Data_Analysis")

# Clean workspace - this removes all the objects from the current environment
## Usually you don't have to do this, we are doing it to start this lesson clean 
rm(list = ls())

# Create a directory for the data
## You might already have this directory, in which case the function issues a warning
dir.create("data")

# Download the data provided by your collaborator
download.file("https://github.com/tavareshugo/data-carpentry-rnaseq/blob/master/data/fission_data.RData?raw=true",
              destfile = "data/fission_data.RData",
              mode = "wb")

# Load the data into R
## "RData" files are special R files that can contain several objects within them
## You might not use these often, tipically you would read data from files
load("data/fission_data.RData")


##### Analysis ####
# Start of your analysis


#convert matric to data.frame
trans_cts_tbl<- trans_cts %>%
  as_tibble(rownames="gene")


trans_cts_long <-trans_cts_tbl %>%
  gather(key="sample", value="cts",
         wt_0_r1:mut_180_r3)


#convert norm_cts to tibble (retain rownames)
#gather samples of norm_cts_tbl

norm_cts_tbl <- norm_cts %>%
  as_tibble(rownames="gene")

norm_cts_tbl_long <- norm_cts_tbl %>%
  gather(key="sample", value="cts",
         wt_0_r1:mut_180_r3)

trans_cts_long <-full_join(trans_cts_long,sample_info)

trans_cts_long %>%
  ggplot(aes(x=cts))+
  geom_histogram()

trans_cts_long <-full_join(trans_cts_long,sample_info,by="sample")

trans_cts_long %>%
  +   ggplot(aes(x=cts))+
  +   geom_histogram()




trans_cts_long %>%
     ggplot(aes(x=cts))+
     geom_histogram()+
  facet_grid(strain ~minute)


#scatterplot of two samples to see co-relation
wt_0_r1 vs wt_0_r2

trans_cts_tbl%>%
  ggplot(aes(x=wt_0_r1, y=wt_0_r2))+
  geom_point()+
  geom_abline(colour="brown")
  

trans_cts_tbl%>%
  ggplot(aes(x=wt_15_r1, y=wt_15_r2))+
  geom_point()+
  geom_abline(colour="brown")

trans_cts_tbl%>%
  ggplot(aes(x=wt_30_r1, y=wt_30_r2))+
  geom_point()+
  geom_abline(colour="brown")

trans_cts_tbl%>%
  ggplot(aes(x=wt_0_r1, y=wt_30_r1))+
  geom_point()+
  geom_abline(colour="brown")

#explore relationship between mean and variance of each gene

trans_cts_long%>%
  group_by(gene)%>%
  summarize(mean_cts=mean(cts),
            var_cts=var(cts))%>%
  ggplot(aes(x=mean_cts, y=var_cts))+
  geom_point()



#first need to transpose matrix
#PCA prcomp
prcomp()

sample_pca<- prcomp(t(trans_cts))
t(trans_cts[1:10, 1:10])
