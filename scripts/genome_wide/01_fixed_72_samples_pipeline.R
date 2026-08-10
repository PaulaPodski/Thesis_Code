#!/usr/bin/env Rscript

library(methylKit)
file.list = list("T108_S18c_bismark.cov.gz", "T117_S2_bismark.cov.gz", "T120_S5b_bismark.cov.gz", "T123_S13c_bismark.cov.gz", "T135_S19b_bismark.cov.gz", "T140_S4b_bismark.cov.gz", "T142_S24_bismark.cov.gz", "T148_S1b_bismark.cov.gz", "T171_S5c_bismark.cov.gz", "T176_S3_bismark.cov.gz", "T191_S2c_bismark.cov.gz", "T19_S13b_bismark.cov.gz", "T207_S1c_bismark.cov.gz", "T20_S18b_bismark.cov.gz", "T211_S15b_bismark.cov.gz", "T214_S4_bismark.cov.gz", "T228_S19c_bismark.cov.gz", "T231_S20c_bismark.cov.gz", "T244_S5_bismark.cov.gz", "T24_S6_bismark.cov.gz", "T275_S21c_bismark.cov.gz", "T284_S11c_bismark.cov.gz", "T303_S3b_bismark.cov.gz", "T313_S7c_bismark.cov.gz", "T317_S7_bismark.cov.gz", "T357_S8_bismark.cov.gz", "T381_S6c_bismark.cov.gz", "T419_S9_bismark.cov.gz", "T448_S20b_bismark.cov.gz", "T449_S3c_bismark.cov.gz", "T452_S10_bismark.cov.gz", "T474_S11_bismark.cov.gz", "T487_S21b_bismark.cov.gz", "T493_S12_bismark.cov.gz", "T494_S8c_bismark.cov.gz", "T505_S22c_bismark.cov.gz", "T510_S13_bismark.cov.gz", "T515_S9c_bismark.cov.gz", "T529_S14_bismark.cov.gz", "T574_S7b_bismark.cov.gz", "T576_S10b_bismark.cov.gz", "T578_S15_bismark.cov.gz", "T57_S1_bismark.cov.gz", "T593_S2b_bismark.cov.gz", "T60_S15c_bismark.cov.gz", "T61_S17c_bismark.cov.gz", "T638_S22b_bismark.cov.gz", "T653_S16_bismark.cov.gz", "T664_S16b_bismark.cov.gz", "T674_S11b_bismark.cov.gz", "T676_S14b_bismark.cov.gz", "T695_S12b_bismark.cov.gz", "T698_S17_bismark.cov.gz", "T706_S18_bismark.cov.gz", "T710_S16c_bismark.cov.gz", "T713_S23b_bismark.cov.gz", "T724_S24b_bismark.cov.gz", "T730_S10c_bismark.cov.gz", "T734_S19_bismark.cov.gz", "T738_S23c_bismark.cov.gz", "T747_S24c_bismark.cov.gz", "T749_S8b_bismark.cov.gz", "T751_S20_bismark.cov.gz", "T764_S9b_bismark.cov.gz", "T771_S21_bismark.cov.gz", "T781_S22_bismark.cov.gz", "T784_S23_bismark.cov.gz", "T786_S6b_bismark.cov.gz", "T797_S14c_bismark.cov.gz", "T799_S12c_bismark.cov.gz", "T83_S4c_bismark.cov.gz", "T8_S17b_bismark.cov.gz")

myobj = methRead(
  file.list,
  sample.id = list(
    "T108", "T117", "T120", "T123", "T135", "T140", "T142", "T148",
    "T171", "T176", "T191", "T19", "T207", "T20", "T211", "T214",
    "T228", "T231", "T244", "T24", "T275", "T284", "T303", "T313",
    "T317", "T357", "T381", "T419", "T448", "T449", "T452", "T474",
    "T487", "T493", "T494", "T505", "T510", "T515", "T529", "T574",
    "T576", "T578", "T57", "T593", "T60", "T61", "T638", "T653",
    "T664", "T674", "T676", "T695", "T698", "T706", "T710", "T713",
    "T724", "T730", "T734", "T738", "T747", "T749", "T751", "T764",
    "T771", "T781", "T784", "T786", "T797", "T799", "T83", "T8"
  ),
  assembly = "hg38",
  treatment = c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1),
  context = "CpG",
  pipeline = "bismarkCoverage",
  mincov = 10,
  dbtype = "tabix",
  dbdir = "methylDB_fixed_3L"
)

filtered.myobj = filterByCoverage(
  myobj,
  lo.count = 10,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = 99.9
)

normalized.myobj = normalizeCoverage(filtered.myobj)
saveRDS(normalized.myobj, "normalized.myobj.rds")

meth = unite(normalized.myobj, destrand = FALSE, min.per.group = 3L)
saveRDS(meth, "meth_min3.rds")

beta_matrix = percMethylation(meth) / 100
saveRDS(beta_matrix, "beta_matrix_min3_fixed.rds")
