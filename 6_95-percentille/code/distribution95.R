setwd("/Users/brigittetombers/Desktop/Project/8.Distribution")
#setwd("/Users/ninatombers/Desktop/Project/8.Distribution")

library(tidyverse)
library(plyranges)
library(hypogen)
# this is just a factor to translate the number of haplotypes into percentages
# 100 / (n * (n-1)) 
hap_to_perc <- 100 / (49 * 48)
# we start by reading in the truffle segments
data_seg <- vroom::vroom(glue::glue("/Users/brigittetombers/Desktop/Project/8.Distribution/truffle_1k_0.5k_nofilt.segments.tsv"),
                         delim = "\t",
                         col_types = "cccciidcdci") %>%
  # I want to plot along the genome, so here I add the karyotype
  # to add the linkage group starting position to the segments
  left_join(hypogen::hypo_chrom_start) %>%
  # then the coordinates of the segments are computed
  # (both within the individual linakge group, and on the concatianted genome)
  mutate(start = POS * 10^6,                   # start position on the LG (Mb -> bp)
         end = start + (LENGTH * 10^6),        # end position on the LG
         gstart = GSTART + start,              # start position on the genome
         gend = start + (LENGTH * 10^6),       # end postion on the genome 
         ibd_hplo = str_remove(TYPE,"IBD") %>% # ibd type as integer (1 ore 2), needed for coverage weighting
           as.integer())
data_seg %>%
  dplyr::select(seqnames = CHROM,start,end,gstart,GSTART,TYPE,ibd_hplo,ID1,ID2) %>%
  arrange(gstart) %>%
  dplyr::select(-gstart) %>%
  as_granges() %>%
  GenomicRanges::coverage(weight = "ibd_hplo") %>%
  plyranges::as_ranges() %>%
  as_tibble() %>%
  dplyr::select(CHROM = seqnames, start, end, width, score) %>%
  left_join(hypogen::hypo_chrom_start) %>%
  mutate(gstart = GSTART + start, gend = GSTART + end) %>%
  dplyr::select(CHROM, gstart, gend, score) %>%
  pivot_longer(gstart:gend,values_to = "GPOS", names_to = "PART") %>%
  ggplot() +
  geom_hypo_LG() +
  geom_step(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.5)) +
  geom_hline(yintercept = 22.5, linetype = 3, color = "red") +
  ylab("IBD score (1k/0.5k 10^3 SNPs)") + xlab("Linkage Group") +
  scale_hypobg_manual(values = c("transparent",rgb(.9,.9,.9,.9),"red","blue") %>%
                        set_names(nm = c("even", "odd", "a","b")), guide = "none")+
  scale_x_hypo_LG() +
  theme_hypo()
  
cummulative_coverage <- data_seg %>%
  dplyr::select(seqnames = CHROM,start,end,gstart,GSTART,TYPE,ibd_hplo,ID1,ID2) %>%
  arrange(gstart) %>%
  dplyr::select(-gstart) %>%
  as_granges() %>%
  GenomicRanges::coverage(weight = "ibd_hplo") %>%
  plyranges::as_ranges() %>%
  as_tibble() %>%
  dplyr::select(CHROM = seqnames, start, end, width, score) %>%
  left_join(hypogen::hypo_chrom_start) %>%
  mutate(gstart = GSTART + start, gend = GSTART + end) %>%
  dplyr::select(CHROM, gstart, gend, score) %>%
  arrange(score) %>% 
  mutate(length = gend - gstart,
         cummlength = cumsum(length))

total_cov_length <- max(cummulative_coverage$cummlength)

percentile_cutoff <- .95

threshold_position <- percentile_cutoff * total_cov_length

cummulative_coverage %>% 
  mutate(below_threshold = cummlength < threshold_position) %>% 
  ggplot(aes(x =cummlength, y = score, color = below_threshold))+
  geom_hypo_LG() +
  scale_hypobg_manual(values = c("transparent",rgb(.9,.9,.9,.9),"red","blue") %>%
                        set_names(nm = c("even", "odd", "a","b")), guide = "none") +
  geom_line()+
  geom_vline(xintercept = threshold_position,
             linetype = 3, color = "red") +
  geom_hline(yintercept = 22.5 / hap_to_perc, linetype = 3, color = "red") +
  scale_x_hypo_LG() +
  theme_hypo()

cummulative_coverage %>% 
  mutate(below_threshold = cummlength < threshold_position) %>% 
  group_by(below_threshold) %>% 
  summarise(min_cov = min(score) * hap_to_perc,
            max_cov = max(score) * hap_to_perc)

filtered_segments <- cummulative_coverage %>% 
  mutate(above_threshold = cummlength > threshold_position) %>% 
  filter(above_threshold) %>% 
  arrange(gstart) %>%
  left_join(hypogen::hypo_chrom_start) %>% 
  mutate(start = gstart - GSTART,
         end = gend - GSTART)

filtered_segments %>% 
  pivot_longer(gstart:gend,values_to = "GPOS", names_to = "PART") %>%
  ggplot() +
  geom_hypo_LG() +
  geom_step(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.5)) +
  scale_hypobg_manual(values = c("transparent",rgb(.9,.9,.9,.9),"red","blue") %>%
                        set_names(nm = c("even", "odd", "a","b")), guide = "none")
  scale_x_hypo_LG() +
  theme_hypo()

 filtered_segments %>% 
    dplyr::select(CHROM, start, end) %>% 
   write_tsv("thresholds_95.bed")
 
 filtered_segments %>% 
   write_tsv("truffle1k0.5k_95.segments")
 
 
 
