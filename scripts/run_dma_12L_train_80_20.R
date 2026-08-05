#!/usr/bin/env Rscript
library(methylKit)

sample.ids <- c("T108","T117","T120","T123","T135","T140","T142","T148",
                 "T171","T176","T191","T19","T207","T20","T211","T214",
                 "T228","T231","T244","T24","T275","T284","T303","T313",
                 "T317","T357","T381","T419","T448","T449","T452","T474",
                 "T487","T493","T494","T505","T510","T515","T529","T574",
                 "T576","T578","T57","T593","T60","T61","T638","T653",
                 "T664","T674","T676","T695","T698","T706","T710","T713",
                 "T724","T730","T734","T738","T747","T749","T751","T764",
                 "T771","T781","T784","T786","T797","T799","T83","T8")

file.list <- list("T108_S18c_bismark.cov.gz", "T117_S2_bismark.cov.gz", "T120_S5b_bismark.cov.gz",
  "T123_S13c_bismark.cov.gz", "T135_S19b_bismark.cov.gz", "T140_S4b_bismark.cov.gz",
  "T142_S24_bismark.cov.gz", "T148_S1b_bismark.cov.gz", "T171_S5c_bismark.cov.gz",
  "T176_S3_bismark.cov.gz", "T191_S2c_bismark.cov.gz", "T19_S13b_bismark.cov.gz",
  "T207_S1c_bismark.cov.gz", "T20_S18b_bismark.cov.gz", "T211_S15b_bismark.cov.gz",
  "T214_S4_bismark.cov.gz", "T228_S19c_bismark.cov.gz", "T231_S20c_bismark.cov.gz",
  "T244_S5_bismark.cov.gz", "T24_S6_bismark.cov.gz", "T275_S21c_bismark.cov.gz",
  "T284_S11c_bismark.cov.gz", "T303_S3b_bismark.cov.gz", "T313_S7c_bismark.cov.gz",
  "T317_S7_bismark.cov.gz", "T357_S8_bismark.cov.gz", "T381_S6c_bismark.cov.gz",
  "T419_S9_bismark.cov.gz", "T448_S20b_bismark.cov.gz", "T449_S3c_bismark.cov.gz",
  "T452_S10_bismark.cov.gz", "T474_S11_bismark.cov.gz", "T487_S21b_bismark.cov.gz",
  "T493_S12_bismark.cov.gz", "T494_S8c_bismark.cov.gz", "T505_S22c_bismark.cov.gz",
  "T510_S13_bismark.cov.gz", "T515_S9c_bismark.cov.gz", "T529_S14_bismark.cov.gz",
  "T574_S7b_bismark.cov.gz", "T576_S10b_bismark.cov.gz", "T578_S15_bismark.cov.gz",
  "T57_S1_bismark.cov.gz", "T593_S2b_bismark.cov.gz", "T60_S15c_bismark.cov.gz",
  "T61_S17c_bismark.cov.gz", "T638_S22b_bismark.cov.gz", "T653_S16_bismark.cov.gz",
  "T664_S16b_bismark.cov.gz", "T674_S11b_bismark.cov.gz", "T676_S14b_bismark.cov.gz",
  "T695_S12b_bismark.cov.gz", "T698_S17_bismark.cov.gz", "T706_S18_bismark.cov.gz",
  "T710_S16c_bismark.cov.gz", "T713_S23b_bismark.cov.gz", "T724_S24b_bismark.cov.gz",
  "T730_S10c_bismark.cov.gz", "T734_S19_bismark.cov.gz", "T738_S23c_bismark.cov.gz",
  "T747_S24c_bismark.cov.gz", "T749_S8b_bismark.cov.gz", "T751_S20_bismark.cov.gz",
  "T764_S9b_bismark.cov.gz", "T771_S21_bismark.cov.gz", "T781_S22_bismark.cov.gz",
  "T784_S23_bismark.cov.gz", "T786_S6b_bismark.cov.gz", "T797_S14c_bismark.cov.gz",
  "T799_S12c_bismark.cov.gz", "T83_S4c_bismark.cov.gz", "T8_S17b_bismark.cov.gz")

treatment <- c(1,0,1,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,
               1,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,1,1,
               0,0,1,0,0,1,0,0,1,1,1,1)

names(file.list) <- sample.ids

#Splitting the data randomly using set.seed(42) for 80/20 split for 72 samples:


stopifnot(length(sample.ids) == 72, length(treatment) == 72, length(file.list) == 72)
table(treatment)
#treatment
#0  1
#36 36
set.seed(42)
case_ids    <- sample.ids[treatment == 1]   # 36 dead
control_ids <- sample.ids[treatment == 0]
train_cases    <- sample(case_ids, 29) #random selection here for dead group (R draws n elements from vector x without replacement)
test_cases     <- setdiff(case_ids, train_cases)
train_controls <- sample(control_ids, 29) # random selection here for alive group (R draws n elements from vector x without replacement)
test_controls  <- setdiff(control_ids, train_controls)
train_ids <- c(train_cases, train_controls)
test_ids  <- c(test_cases, test_controls)

train_treatment <- treatment[match(train_ids, sample.ids)]
test_treatment  <- treatment[match(test_ids, sample.ids)]
stopifnot(length(train_ids) == 58, sum(train_treatment) == 29, sum(train_treatment == 0) == 29)
stopifnot(length(test_ids) == 14, sum(test_treatment) == 7, sum(test_treatment == 0) == 7)
saveRDS(list(train_ids = train_ids, test_ids = test_ids,
             train_treatment = train_treatment, test_treatment = test_treatment,
             seed = 42),
        "sample_split_80_20_min12.rds")
split <- readRDS("sample_split_80_20_min12.rds")

#See what's inside
str(split)
#List of 5
 #$ train_ids      : chr [1:50] "T108" "T710" "T303" "T576" ...
 #$ test_ids       : chr [1:22] "T19" "T24" "T313" "T419" ...
 #$ train_treatment: num [1:50] 1 1 1 1 1 1 1 1 1 1 ...
 #$ test_treatment : num [1:22] 1 1 1 1 1 1 1 1 1 1 ...
 #$ seed           : num 42
length(split$train_ids)   
length(split$test_ids)    
#[1] 58
#[1] 14
table(split$train_treatment)   
table(split$test_treatment)    

# 0  1
#29 29

 #0  1
#7 7
intersect(split$train_ids, split$test_ids)
character(0)
all_ids <- c(split$train_ids, split$test_ids)
length(all_ids) == 72                         
length(unique(all_ids)) == 72                 
setdiff(sample.ids, all_ids)
#[1] TRUE
#[1] TRUE #this confirms there are no dupicate samples in each group
# character(0)
identical(
  split$train_treatment,
  treatment[match(split$train_ids, sample.ids)]
)  
identical(
  split$test_treatment,
  treatment[match(split$test_ids, sample.ids)]
)  
#[1] TRUE
#[1] TRUE
set.seed(split$seed)
case_ids_check    <- sample.ids[treatment == 1]
control_ids_check <- sample.ids[treatment == 0]
train_cases_check    <- sample(case_ids_check, 29)
train_controls_check <- sample(control_ids_check, 29)

identical(sort(train_cases_check), sort(split$train_ids[split$train_treatment == 1]))
identical(sort(train_controls_check), sort(split$train_ids[split$train_treatment == 0]))
#[1] TRUE
#[1] TRUE
train_file.list <- unname(as.list(file.list[train_ids]))
test_file.list  <- unname(as.list(file.list[test_ids]))

saveRDS(train_file.list, "train_file_list_80_20_min12.rds")
saveRDS(test_file.list,  "test_file_list_80_20_min12.rds")
split           <- readRDS("sample_split_80_20_min12.rds")
train_file.list <- readRDS("train_file_list_80_20_min12.rds")

train_ids       <- split$train_ids
train_treatment <- split$train_treatment

# Checking if data is split correctly into groups and cases/controls: 
stopifnot(length(train_file.list) == length(train_ids),
          length(train_ids) == 58,
          sum(train_treatment) == 29, sum(train_treatment == 0) == 29)
length(train_file.list)      
length(train_ids)            
sum(train_treatment)         
sum(train_treatment == 0)    
#[1] 58
#[1] 58
#[1] 29
#[1] 29
myobj_train <- methRead(
  train_file.list,
  sample.id = as.list(train_ids),
  assembly = "hg38",
  treatment = train_treatment,
  context = "CpG",
  pipeline = "bismarkCoverage",
  mincov = 10,
  dbtype = "tabix",
  dbdir = "methylDB_train_80_20_min12"
)

filtered_train   <- filterByCoverage(myobj_train, lo.count = 10, lo.perc = NULL,
                                      hi.count = NULL, hi.perc = 99.9)

normalized_train <- normalizeCoverage(filtered_train)
saveRDS(normalized_train, "normalized_train_80_20_min12.rds")

meth_train <- unite(normalized_train, destrand = FALSE, min.per.group = 12L)
saveRDS(meth_train, "meth_train_80_20_min12.rds")

beta_train <- percMethylation(meth_train) / 100
saveRDS(beta_train, "beta_matrix_train_80_20_min12.rds")

myDiff_train <- calculateDiffMeth(meth_train, mc.cores = 24)
saveRDS(myDiff_train, "myDiff_train_80_20_min12.rds")
