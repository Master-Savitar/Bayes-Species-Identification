library(tidyverse)

# read in audio clip files
clip_taxon_annotation_dat <- read.table(
  file = "clips/clip_taxon_annotations.tsv", header = TRUE, sep = "\t")
clip_taxon_annotation_dat <- subset(clip_taxon_annotation_dat, 
                                    is_extra_annotation == "False")

clip_annotation_dat <- read.table(
  file = "clips/clip_annotations.tsv", header = TRUE, sep = "\t")

# convert the annotation files into a 3-D binary array
## the number of audio clips: 3997
N1 <- clip_annotation_dat$clip_id %>% unique() %>% length()
clips <- clip_annotation_dat$clip_id %>% unique() %>% sort()
## the number of raters/users: 46
N2 <- clip_annotation_dat$user_id %>% unique() %>% length()
users <- clip_annotation_dat$user_id %>% unique() %>% sort()
## the number of bird species: 349
N3 <- clip_taxon_annotation_dat$taxon_id %>% unique() %>% length()
birds <- clip_taxon_annotation_dat$taxon_id %>% unique() %>% sort()

binary_array <- array(data = rep(NA, N1*N2*N3), 
                      dim = c(N1, N2, N3),
                      dimnames = list(clips, users, birds))

for (i in 1:nrow(clip_annotation_dat)) {
  
  # i is the row_index of the clip_annotation dataset
  # which isn't the clip_annotation_id
  # however, there relationship is one to one
  # specifically, the max of clip_annotation_id is 5972
  # the number of columns is 5932
  
  clip_id <- clip_annotation_dat$clip_id[i]
  user_id <- clip_annotation_dat$user_id[i]
  annotation_id <- clip_annotation_dat$clip_annotation_id[i]
  
  tmp <- clip_taxon_annotation_dat[clip_taxon_annotation_dat$clip_annotation_id == annotation_id, ]
  binary_array[as.character(clip_id), user_id, tmp$taxon_id] <- tmp$annotation
}
binary_array[binary_array == 2] <- NA

# save the 3-D binary data
saveRDS(binary_array, "annotate_dat.rds")      #  0.9712575 




