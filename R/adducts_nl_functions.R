
#' Detect possible neutral losses, isotopes and adduct links.
#'
#' @eval recurrent_params("threshold_rt", "threshold_mz", "user_neg_neutral_refs,user_pos_neutral_refs,user_pos_adducts_refs,user_neg_adducts_refs")
#' @param data_to_treat A data.frame containing peaks data.
#' @param samples_cols A list of string containing the names of samples columns.
#' @param detected_adducts A data.frame indicating which adducts were detected by MSDial.
#' @return A list of 3 data.frames: "links" containing detected links, "computed_massdiff" containing all pairs of positive/negative adducts pair mass differences and "filtered_massdiff" containing pairs of positive/negative adducts pair mass differences for adducts pair detected by MSDial.
get_adducts_nl_links <- function(data_to_treat,
                                 samples_cols,
                                 threshold_rt,
                                 threshold_mz,
                                 detected_adducts,
                                 user_pos_neutral_refs = NA,
                                 user_neg_neutral_refs = NA,
                                 user_pos_adducts_refs = NA,
                                 user_neg_adducts_refs = NA) {
    clusters_pos <- unique(data_to_treat[data_to_treat$source == "pos",]$cluster.msdial)
    clusters_neg <- unique(data_to_treat[data_to_treat$source == "neg",]$cluster.msdial)
    msc_data <- data_to_treat[, c(c("cluster.msdial", "id", "Average.Mz", "Average.Rt.min."), samples_cols)]
    msc_links <- data.frame()

    utils::data("mass_isotopes", package="mscleanr", envir=environment())
    mass_isotopes <- get("mass_isotopes", envir=environment())  # redundant, but makes the syntax checker happy

    # pos
    if (is.na(user_pos_neutral_refs)) {
        print_message("Using package neutral losses for positive mode")
        utils::data("mass_neutral_loss_pos", package="mscleanr", envir=environment())
        current_mass_nl <- get("mass_neutral_loss_pos", envir=environment())  # redundant, but makes the syntax checker happy
    } else {
        print_message("Using personalized reference file for neutral losses in positive mode")
        current_mass_nl <- user_pos_neutral_refs
        names(current_mass_nl) <- c("name", "m")
    }
    for (cpos in clusters_pos) {
        print_message("Adducts/Neutral losses detection (cluster pos ", cpos, ")")
        lks <- detect_single_polarity_combinations(msc_data[msc_data$cluster.msdial == cpos, 2:ncol(msc_data)],
                                                   threshold_rt,
                                                   threshold_mz,
                                                   mass_neutral_loss = current_mass_nl,
                                                   mass_isotopes = mass_isotopes)
        if(!is.null(lks)) msc_links <- rbind(msc_links, lks)
    }

    # neg
    if (is.na(user_neg_neutral_refs)) {
        print_message("Using package neutral losses for negative mode")
        utils::data("mass_neutral_loss_neg", package="mscleanr", envir=environment())
        current_mass_nl <- get("mass_neutral_loss_neg", envir=environment())  # redundant, but makes the syntax checker happy
    } else {
        print_message("Using personalized reference file for neutral losses in negative mode")
        current_mass_nl <- user_neg_neutral_refs
        names(current_mass_nl) <- c("name", "m")
    }
    for (cneg in clusters_neg) {
        print_message("Adducts/Neutral losses detection (cluster neg ", cneg, ")")
        lks <- detect_single_polarity_combinations(msc_data[msc_data$cluster.msdial == cneg, 2:ncol(msc_data)],
                                                   threshold_rt,
                                                   threshold_mz,
                                                   mass_neutral_loss = current_mass_nl,
                                                   mass_isotopes = mass_isotopes)
        if(!is.null(lks)) msc_links <- rbind(msc_links, lks)
    }

    # pos/neg: computing massdiff + filter on detected adducts
    all_adducts_posneg <- compute_massdiff(user_adduct_pos = user_pos_adducts_refs,
                                           user_adduct_neg = user_neg_adducts_refs,
                                           significance_mass_threshold = threshold_mz)
    selected_pos_adducts <- all_adducts_posneg$adduct.pos %in% detected_adducts[detected_adducts$source == "pos",]$adduct
    selected_neg_adducts <- all_adducts_posneg$adduct.neg %in% detected_adducts[detected_adducts$source == "neg",]$adduct
    adducts_posneg <- all_adducts_posneg[selected_pos_adducts & selected_neg_adducts,]

    for (cpos in clusters_pos) {
        for (cneg in clusters_neg) {
            # if clusters have peaks close enough for links to exist -> Adducts/Neutral losses detection
            # I.E.: if min(abs(rt1 - rt2) foreach pair of peaks) <= threshold_rt -> Adducts/Neutral losses detection
            rts <- expand.grid(msc_data[msc_data$cluster.msdial == cpos,]$Average.Rt.min.,
                               msc_data[msc_data$cluster.msdial == cneg,]$Average.Rt.min.)
            delta_min <- min(abs(rts$Var1 - rts$Var2))
            if (delta_min <= threshold_rt) {
                print_message("Adducts/Neutral losses detection (cluster pos ", cpos, " / neg ", cneg, ")")
                lks <- detect_posneg_combinations(msc_data[msc_data$cluster.msdial == cpos, 2:ncol(msc_data)],
                                                  msc_data[msc_data$cluster.msdial == cneg, 2:ncol(msc_data)],
                                                  threshold_rt,
                                                  threshold_mz,
                                                  adducts_posneg)
                if(!is.null(lks)) msc_links <- rbind(msc_links, lks)
            }
        }
    }

    msc_links$Adduct.1 <- as.character(msc_links$Adduct.1)
    msc_links$Adduct.2 <- as.character(msc_links$Adduct.2)

    return(list(links             = msc_links,
                computed_massdiff = all_adducts_posneg,
                filtered_massdiff = adducts_posneg))
}



#' Detect possible adducts and neutral losses.
#'
#' @eval recurrent_params("threshold_rt", "threshold_mz")
#' @param peaks_data A data.frame containing peaks to treat, coming from MSDial.
#' @param mass_neutral_loss An 2-column data.frame containing information about neutral losses (one column for the name and one column for the mass difference with the base compound).
#' @param mass_isotopes An 2-column data.frame containing information about isotopes (one column for the name and one column for the mass difference with the base compound).
#' @return A data.frame of possible links due to neutral losses and isotopes.
detect_single_polarity_combinations <- function(peaks_data, threshold_rt, threshold_mz, mass_neutral_loss, mass_isotopes) {
    if (nrow(peaks_data) == 0) return(NULL)

    to_detect <- list(
        list(l1 = peaks_data, l2 = peaks_data, massdiff = mass_neutral_loss[c("name", "name", "m")], nature = "neutral loss"),
        list(l1 = peaks_data, l2 = peaks_data, massdiff = mass_isotopes[c("name", "name", "m")],     nature = "isotope"),
        list(l1 = peaks_data, l2 = peaks_data, massdiff = peaks_data[c("id", "id", "Average.Mz")],   nature = "mers")
    )

    combinations <- data.frame()
    for (cdet in to_detect) {
        com <- MScombine::FindCommon(cdet$l1,
                                     cdet$l2,
                                     cdet$massdiff,
                                     Masstolerance = threshold_mz,
                                     RTtolerance = threshold_rt)
        if(nrow(com) > 0) com$simple.nature <- cdet$nature
        combinations <- rbind(combinations, com)
    }

    if (nrow(combinations) == 0) return(NULL)
    else {
        combinations <- clean_msc_combinations(combinations)
        is_mer <- combinations$simple.nature == "mers"
        if (nrow(combinations[is_mer & combinations$Adduct.1 == combinations$CpdID.2,]) > 0) {
            combinations[is_mer & combinations$Adduct.1 == combinations$CpdID.2,]$simple.nature <- "homomers"
        }
        if (nrow(combinations[is_mer & combinations$Adduct.1 != combinations$CpdID.2,]) > 0) {
            combinations[is_mer & combinations$Adduct.1 != combinations$CpdID.2,]$simple.nature <- "dimers"
        }
        return(combinations)
    }
}



#' Detect possible adducts, neutral losses and positive and negative peaks to combine.
#'
#' @eval recurrent_params("threshold_rt", "threshold_mz")
#' @param peaks_data_pos A data.frame containing positive peaks to treat, coming from MSDial.
#' @param peaks_data_neg A data.frame containing negative peaks to treat, coming from MSDial.
#' @param adducts_massdiff A data.frame containing combinations information (adduct1 adduct2 m).
#' @return A data.frame of possible links between positive and negative peaks.
detect_posneg_combinations <- function(peaks_data_pos, peaks_data_neg, threshold_rt, threshold_mz, adducts_massdiff) {
    if (nrow(peaks_data_pos) == 0)   return(NULL)
    if (nrow(peaks_data_neg) == 0)   return(NULL)
    if (nrow(adducts_massdiff) == 0) return(NULL)

    combinations <- MScombine::FindCommon(peaks_data_pos,
                                          peaks_data_neg,
                                          adducts_massdiff,
                                          Masstolerance = threshold_mz,
                                          RTtolerance = threshold_rt)

    if(nrow(combinations) == 0) return(NULL)
    else {
        combinations$simple.nature <- "pol / adduct"
        return(clean_msc_combinations(combinations))
    }
}



#' Compute mass differences between 2 lists of adducts
#'
#' @param user_adduct_pos,user_adduct_neg An optional 2-column data.frame containing information about positive adducts or negative adducts (one column for the name and one column for the mass difference with the base compound). If no data.frame is provided, the package default list is used.
#' @param significance_mass_threshold A numeric value indicating at which threshold mass differences are discarded because too small.
#' @return A 3-columns data.frame with the pairs of adducts and their mass difference values.
compute_massdiff <- function(user_adduct_pos = NA, user_adduct_neg = NA, significance_mass_threshold = NA) {
    if (is.na(user_adduct_pos)) {
        print_message("Using package positive adducts")
        # https://r.789695.n4.nabble.com/no-visible-binding-for-global-variable-for-data-sets-in-a-package-tp4696053p4696079.html
        utils::data("mass_adducts_pos", package="mscleanr", envir=environment())
        adduct_pos <- get("mass_adducts_pos", envir=environment())  # redundant, but makes the syntax checker happy
    }
    else {
        print_message("Using personalized reference file for positive adducts")
        adduct_pos <- user_adduct_pos
        names(adduct_pos) <- c("adduct", "mass")
    }

    if (is.na(user_adduct_neg)) {
        print_message("Using package negative adducts")
        utils::data("mass_adducts_neg", package="mscleanr", envir=environment())
        adduct_neg <- get("mass_adducts_neg", envir=environment())  # redundant, but makes the syntax checker happy
    }
    else {
        print_message("Using personalized reference file for negative adducts")
        adduct_neg <- user_adduct_neg
        names(adduct_pos) <- c("adduct", "mass")
    }

    diff <- expand.grid(adduct_pos$adduct, adduct_neg$adduct)
    names(diff) <- c("adduct.pos", "adduct.neg")
    diff <- merge(diff, adduct_pos, by.x = "adduct.pos", by.y = "adduct")
    names(diff)[length(names(diff))] <- "diff1"
    diff <- merge(diff, adduct_neg, by.x = "adduct.neg", by.y = "adduct")
    names(diff)[length(names(diff))] <- "diff2"
    diff$MassDiff <- diff$diff1 - diff$diff2
    diff <- diff[c("adduct.pos", "adduct.neg", "MassDiff")]

    if (!is.na(significance_mass_threshold)) diff <- diff[abs(diff$MassDiff) > significance_mass_threshold,]
    if (nrow(diff) == 0) print_warning("No mass differences remaining between positive and negative adducts, you might want to check your input data or reduce the significance mass threshold.")
    return(diff)
}



#' Cleans the columns name of data generated my MSCombine.
#'
#' @param data_to_clean A data.frame with columns to rename.
#' @return A data.frame containing the same data than \code{data_to_clean} but with new column names.
clean_msc_combinations <- function(data_to_clean) {
    data_to_clean$`Adduct+` <- as.character(data_to_clean$`Adduct+`)
    data_to_clean$`Adduct-` <- as.character(data_to_clean$`Adduct-`)
    names(data_to_clean) <- c("CpdID.1", "CpdID.2", "Adduct.1", "Adduct.2", "Mass.1", "Mass.2", "RT.1", "RT.2",
                              "Mean.1", "Mean.2", "N.1", "N.2", "Correlation", "simple.nature")
    return(data_to_clean)
}
