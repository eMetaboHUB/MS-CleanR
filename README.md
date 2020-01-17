# mscleanr: A package for cleaning and analyzing MS data

The mscleanr package provides 3 important functions: **clean_msdial_data**, **keep_top_peaks** and **launch_msfinder_annotation**.

See the functions documentation and vignettes for more information.

Needs MSDial v4.00 or higher.

### Exemple
```R
devtools::install_github("SyrupType/mscleanr")
library(mscleanr)

choose_project_directory("/Users/ofv/Analyses/CAD-camomilla", overwrite = TRUE)

clean_msdial_data(filter_blk = TRUE,
                  filter_blk_threshold = 0.8,
                  filter_blk_ghost_peaks = TRUE,
                  filter_mz = TRUE,
                  filter_rsd = TRUE,
                  filter_rsd_threshold = 30,
                  filter_rmd = TRUE,
                  filter_rmd_range = c(50, 3000),
                  threshold_mz = 0.005,
                  threshold_rt = 0.025,
                  user_pos_adducts_refs = NA,
                  user_neg_adducts_refs = NA,
                  user_pos_neutral_refs = NA,
                  user_neg_neutral_refs = NA,
                  compute_pearson_correlation = TRUE,
                  pearson_correlation_threshold = 0.8,
                  pearson_p_value = 1)

keep_top_peaks(selection_criterion = "both",
               n = 5,
               export_filtered_peaks = TRUE)

launch_msfinder_annotation(compound_levels = c("1a", "1b"),
                           biosoc_levels = c("genre", "family"),
                           levels_scores = list("1a" = 2, "1b" = 1.5, "genre" = 2, "family" = 1.5),
                           score_only = FALSE)

convert_csv_to_msp(all = TRUE, min_score = 5)


```
