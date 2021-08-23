# Copyright (C) Shixiang Wang
# LICENSE: GPL3

# Other TCGA utils: https://bioconductor.org/packages/release/bioc/vignettes/TCGAutils/inst/doc/TCGAutils.html


# 过滤TCGA重名ID --------------------------------------------------------------

filterReplicates <- function(tsb, analyte_target = c("DNA", "RNA"), 
                             decreasing = TRUE, analyte_position = 20, plate = c(22, 25),
                             portion = c(18, 19), filter_FFPE = FALSE, full_barcode = FALSE) {
  ## Filter TCGA Replicate Samples
  ## Author: ShixiangWang <w_shixiang@163.com>
  ## Update: 2021-08-18, remove which sentence to speed up, use sort and index selection to make sure only unique sample left.
  ## ooooooo
  ## The filter rule following broad institute says:
  ##
  ## In many instances there is more than one aliquot for a given combination of individual, platform, and data type. However, only one aliquot may be ingested into Firehose. Therefore, a set of precedence rules are applied to select the most scientifically advantageous one among them. Two filters are applied to achieve this aim: an Analyte Replicate Filter and a Sort Replicate Filter.
  ##
  ## Analyte Replicate Filter
  ## The following precedence rules are applied when the aliquots have differing analytes. For RNA aliquots, T analytes are dropped in preference to H and R analytes, since T is the inferior extraction protocol. If H and R are encountered, H is the chosen analyte. This is somewhat arbitrary and subject to change, since it is not clear at present whether H or R is the better protocol. If there are multiple aliquots associated with the chosen RNA analyte, the aliquot with the later plate number is chosen. For DNA aliquots, D analytes (native DNA) are preferred over G, W, or X (whole-genome amplified) analytes, unless the G, W, or X analyte sample has a higher plate number.
  ##
  ## Sort Replicate Filter
  ## The following precedence rules are applied when the analyte filter still produces more than one sample. The sort filter chooses the aliquot with the highest lexicographical sort value, to ensure that the barcode with the highest portion and/or plate number is selected when all other barcode fields are identical.
  ## oooooo
  ## Ref Link: <https://confluence.broadinstitute.org/display/GDAC/FAQ#FAQ-sampleTypesQWhatTCGAsampletypesareFirehosepipelinesexecutedupon>


  # basically, user provide tsb and analyte_target is fine. If you
  # want to filter FFPE samples, please set filter_FFPE and full_barcode
  # all to TRUE, and tsb must have nchar of 28

  analyte_target <- match.arg(analyte_target)
  # Strings in R are largely lexicographic
  # see ??base::Comparison

  # filter FFPE samples
  # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html>
  if (full_barcode & filter_FFPE) {
    ffpe <- c(
      "TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01",
      "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26",
      "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08",
      "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07",
      "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01",
      "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07",
      "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01",
      "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26",
      "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08",
      "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05",
      "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07",
      "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01",
      "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26",
      "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08",
      "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05",
      "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07",
      "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01",
      "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26",
      "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08",
      "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05",
      "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07",
      "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01",
      "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26",
      "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13",
      "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01",
      "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26",
      "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13",
      "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01",
      "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26",
      "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13",
      "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07",
      "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01",
      "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26",
      "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10",
      "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05",
      "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07",
      "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01",
      "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26",
      "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10",
      "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05",
      "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07",
      "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01",
      "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26",
      "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13",
      "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01",
      "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26",
      "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10",
      "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05",
      "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07",
      "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10",
      "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05",
      "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07",
      "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10",
      "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05",
      "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13",
      "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07",
      "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09",
      "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13",
      "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07",
      "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09",
      "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05",
      "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13",
      "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01",
      "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26",
      "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13",
      "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07",
      "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10",
      "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05",
      "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07",
      "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10",
      "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05",
      "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07",
      "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10",
      "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05",
      "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07",
      "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09",
      "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05",
      "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07",
      "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09",
      "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05",
      "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13",
      "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01",
      "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26",
      "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13",
      "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01",
      "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26",
      "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13",
      "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01",
      "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26",
      "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13",
      "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05",
      "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13",
      "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01",
      "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26",
      "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13"
    )

    tsb <- setdiff(tsb, tsb[tsb %in% ffpe])
  }

  # find repeated samples
  sampleID <- substr(tsb, start = 1, stop = 15)
  dp_samples <- unique(sampleID[duplicated(sampleID)])

  if (length(dp_samples) == 0) {
    message("ooo Not find any duplicated barcodes, return original input..")
    tsb
  } else {
    uniq_tsb <- tsb[!sampleID %in% dp_samples]
    dp_tsb <- setdiff(tsb, uniq_tsb)

    add_tsb <- c()

    # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
    # if analyte_target = "DNA"
    # analyte:  D > G,W,X
    if (analyte_target == "DNA") {
      for (x in dp_samples) {
        mulaliquots <- dp_tsb[substr(dp_tsb, 1, 15) == x]
        analytes <- substr(mulaliquots,
          start = analyte_position,
          stop = analyte_position
        )
        if (any(analytes == "D") & !(all(analytes == "D"))) {
          aliquot <- mulaliquots[analytes == "D"]

          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb <- c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb <- c(add_tsb, aliquot)
            dp_tsb <- setdiff(dp_tsb, mulaliquots)
          }
        }
      }
    } else {
      # if analyte_target = "RNA"
      # analyte: H > R > T
      for (x in dp_samples) {
        mulaliquots <- dp_tsb[substr(dp_tsb, 1, 15) == x]
        analytes <- substr(mulaliquots,
          start = analyte_position,
          stop = analyte_position
        )
        if (any(analytes == "H") & !(all(analytes == "H"))) {
          aliquot <- mulaliquots[analytes == "H"]

          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb <- c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb <- c(add_tsb, aliquot)
            dp_tsb <- setdiff(dp_tsb, mulaliquots)
          }
        } else if (any(analytes == "R") & !(all(analytes == "R"))) {
          aliquot <- mulaliquots[analytes == "R"]

          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb <- c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb <- c(add_tsb, aliquot)
            dp_tsb <- setdiff(dp_tsb, mulaliquots)
          }
        } else if (any(analytes == "T") & !(all(analytes == "T"))) {
          aliquot <- mulaliquots[analytes == "T"]

          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb <- c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb <- c(add_tsb, aliquot)
            dp_tsb <- setdiff(dp_tsb, mulaliquots)
          }
        }
      }
    }


    if (length(dp_tsb) == 0) {
      message("ooo Filter barcodes successfully!")
      c(uniq_tsb, add_tsb)
    } else {
      # filter according to portion number
      sampleID_res <- substr(dp_tsb, start = 1, stop = 15)
      dp_samples_res <- unique(sampleID_res[duplicated(sampleID_res)])

      for (x in dp_samples_res) {
        mulaliquots <- dp_tsb[substr(dp_tsb, 1, 15) == x]
        portion_codes <- substr(mulaliquots,
          start = portion[1],
          stop = portion[2]
        )
        portion_keep <- sort(portion_codes, decreasing = decreasing)[1]
        if (!all(portion_codes == portion_keep)) {
          if (sum(portion_codes == portion_keep) == 1) {
            add_tsb <- c(add_tsb, mulaliquots[portion_codes == portion_keep])
            dp_tsb <- setdiff(dp_tsb, mulaliquots)
          } else {
            dp_tsb <- setdiff(dp_tsb, mulaliquots[portion_codes != portion_keep])
          }
        }
      }

      if (length(dp_tsb) == 0) {
        message("ooo Filter barcodes successfully!")
        c(uniq_tsb, add_tsb)
      } else {
        # filter according to plate number
        sampleID_res <- substr(dp_tsb, start = 1, stop = 15)
        dp_samples_res <- unique(sampleID_res[duplicated(sampleID_res)])
        for (x in dp_samples_res) {
          mulaliquots <- dp_tsb[substr(dp_tsb, 1, 15) == x]
          plate_codes <- substr(mulaliquots,
            start = plate[1],
            stop = plate[2]
          )
          plate_keep <- sort(plate_codes, decreasing = decreasing)[1]
          add_tsb <- c(add_tsb, sort(mulaliquots[plate_codes == plate_keep])[1])
          dp_tsb <- setdiff(dp_tsb, mulaliquots)
        }

        if (length(dp_tsb) == 0) {
          message("ooo Filter barcodes successfully!")
          c(uniq_tsb, add_tsb)
        } else {
          message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
          c(uniq_tsb, add_tsb)
        }
      }
    }
  }
}

# translateID ------------------------------------------------------------
translateID <- function(id, type = c("file_id", "case_id", "aliquot_ids"), legacy = FALSE) {
  ## Translcate GDC UUID to sample barcode
  ##
  ## @param a vector of UUID.
  ##
  ## @examples
  ## translateID("719fc310-5ff5-4493-a421-1533583baca1")
  .check_install("TCGAutils", bioc = TRUE)
  type <- match.arg(type)
  TCGAutils::UUIDtoBarcode(id, type, legacy)
}


# Sample mapping between PCAWG and TCGA -----------------------------------

## How the data was collected.
##
# library(UCSCXenaShiny)
# library(dplyr)
# library(IDConverter) # https://github.com/ShixiangWang/IDConverter
# pcawg_tcga_samps <- pcawg_purity %>%
#   dplyr::mutate(submitter_id = convert_pcawg(
#     icgc_specimen_id,
#     from = "icgc_specimen_id",
#     to = "submitted_specimen_id",
#     db = "simple"),
#     submitter_id = substr(submitter_id, 1, 15)) %>%
#   dplyr::filter(startsWith(submitter_id, "TCGA")) %>% 
#   dplyr::select(icgc_specimen_id, submitter_id)
# 
# pcawg_tcga_map <- pcawg_tcga_samps$submitter_id
# names(pcawg_tcga_map) <- pcawg_tcga_samps$icgc_specimen_id

PCAWG_TCGA_SAMPLE_MAP <- c(SP47708 = "TCGA-EZ-7264-01", SP85623 = "TCGA-DJ-A2Q1-01", SP22031 = "TCGA-A6-6141-01", 
  SP36498 = "TCGA-A3-3387-01", SP94588 = "TCGA-AX-A1C8-01", SP7692 = "TCGA-A2-A0EY-01", 
  SP50827 = "TCGA-05-4396-01", SP80037 = "TCGA-HC-7744-01", SP50263 = "TCGA-91-6847-01", 
  SP31046 = "TCGA-CR-6470-01", SP27201 = "TCGA-02-2485-01", SP90725 = "TCGA-AP-A0LL-01", 
  SP98082 = "TCGA-DD-A4NG-01", SP122725 = "TCGA-IW-A3M5-01", SP80367 = "TCGA-AG-3890-01", 
  SP28791 = "TCGA-06-2570-01", SP5052 = "TCGA-EW-A3U0-01", SP13206 = "TCGA-C5-A2LY-01", 
  SP82445 = "TCGA-GN-A266-06", SP83242 = "TCGA-D3-A1Q5-06", SP61703 = "TCGA-24-1544-01", 
  SP96112 = "TCGA-A6-A567-01", SP85491 = "TCGA-EL-A3TB-01", SP19215 = "TCGA-AA-A00R-01", 
  SP37516 = "TCGA-CJ-4639-01", SP13078 = "TCGA-C5-A1BN-01", SP121828 = "TCGA-DX-A3LS-01", 
  SP87582 = "TCGA-DE-A3KN-01", SP7785 = "TCGA-BH-A0BW-01", SP107607 = "TCGA-C5-A1M8-01", 
  SP123972 = "TCGA-KN-8427-01", SP96120 = "TCGA-AY-A54L-01", SP31334 = "TCGA-CN-6989-01", 
  SP23739 = "TCGA-06-0214-01", SP123850 = "TCGA-KL-8330-01", SP16886 = "TCGA-CA-6718-01", 
  SP50317 = "TCGA-91-6840-01", SP59707 = "TCGA-13-1477-01", SP54113 = "TCGA-05-4397-01", 
  SP123902 = "TCGA-KM-8639-01", SP56537 = "TCGA-56-1622-01", SP89090 = "TCGA-DJ-A13W-01", 
  SP94236 = "TCGA-EY-A1GS-01", SP1086 = "TCGA-BL-A13J-01", SP102690 = "TCGA-HC-7233-01", 
  SP85251 = "TCGA-BR-4255-01", SP50406 = "TCGA-75-7030-01", SP49334 = "TCGA-EP-A3RK-01", 
  SP109941 = "TCGA-EK-A2RE-01", SP103894 = "TCGA-EE-A3JI-06", SP57586 = "TCGA-92-8064-01", 
  SP30011 = "TCGA-CV-6433-01", SP43510 = "TCGA-B9-4115-01", SP39997 = "TCGA-A3-3363-01", 
  SP109801 = "TCGA-C5-A2LT-01", SP84408 = "TCGA-BR-7722-01", SP84491 = "TCGA-D7-6822-01", 
  SP93772 = "TCGA-D1-A17K-01", SP96124 = "TCGA-QG-A5YX-01", SP83382 = "TCGA-EE-A2GN-06", 
  SP35617 = "TCGA-CJ-6033-01", SP3415 = "TCGA-AO-A0J4-01", SP49124 = "TCGA-G3-A3CK-01", 
  SP68659 = "TCGA-09-1666-01", SP62616 = "TCGA-23-1110-01", SP59761 = "TCGA-13-1491-01", 
  SP85511 = "TCGA-EM-A22O-01", SP2799 = "TCGA-E9-A1NH-01", SP9979 = "TCGA-A2-A259-01", 
  SP82087 = "TCGA-AG-3896-01", SP82451 = "TCGA-DA-A3F5-06", SP123884 = "TCGA-KM-8438-01", 
  SP13072 = "TCGA-C5-A1MI-01", SP51824 = "TCGA-05-4389-01", SP80615 = "TCGA-F5-6814-01", 
  SP43770 = "TCGA-AL-3473-01", SP66960 = "TCGA-24-1552-01", SP84743 = "TCGA-D7-6528-01", 
  SP19295 = "TCGA-AA-3555-01", SP85952 = "TCGA-FK-A3S3-01", SP92195 = "TCGA-AP-A0L8-01", 
  SP104330 = "TCGA-ER-A2NF-06", SP55387 = "TCGA-55-8299-01", SP5393 = "TCGA-AR-A1AY-01", 
  SP122476 = "TCGA-HB-A43Z-01", SP61453 = "TCGA-13-0727-01", SP31174 = "TCGA-CV-7180-01", 
  SP83967 = "TCGA-FP-7998-01", SP85840 = "TCGA-ET-A4KN-01", SP87675 = "TCGA-EM-A2OW-01", 
  SP63796 = "TCGA-13-0723-01", SP121783 = "TCGA-DX-A1L3-01", SP32662 = "TCGA-CN-6994-01", 
  SP85836 = "TCGA-BJ-A28V-01", SP120767 = "TCGA-DJ-A13L-01", SP52779 = "TCGA-78-7158-01", 
  SP89957 = "TCGA-AP-A0LH-01", SP43651 = "TCGA-AL-3468-01", SP92268 = "TCGA-BS-A0TG-01", 
  SP31126 = "TCGA-CV-7255-01", SP127630 = "TCGA-FF-8047-01", SP11292 = "TCGA-B6-A0RT-01", 
  SP105708 = "TCGA-EM-A3FL-01", SP69155 = "TCGA-24-0982-01", SP115501 = "TCGA-BW-A5NP-01", 
  SP13036 = "TCGA-EW-A1PC-01", SP121774 = "TCGA-DX-A1L0-01", SP56541 = "TCGA-21-1082-01", 
  SP85818 = "TCGA-BJ-A28T-01", SP32958 = "TCGA-CR-6487-01", SP59950 = "TCGA-04-1347-01", 
  SP37636 = "TCGA-BP-4977-01", SP24702 = "TCGA-14-1459-01", SP123964 = "TCGA-KN-8424-01", 
  SP7421 = "TCGA-B6-A0RE-01", SP105253 = "TCGA-BR-8373-01", SP31422 = "TCGA-BA-5149-01", 
  SP56553 = "TCGA-60-2698-01", SP4875 = "TCGA-AO-A0J2-01", SP33496 = "TCGA-CN-6011-01", 
  SP8795 = "TCGA-B6-A0I1-01", SP83312 = "TCGA-GN-A26A-06", SP82614 = "TCGA-GN-A262-06", 
  SP92707 = "TCGA-A5-A0GE-01", SP55235 = "TCGA-05-4398-01", SP50713 = "TCGA-49-4486-01", 
  SP106656 = "TCGA-MH-A55W-01", SP105922 = "TCGA-EM-A3AL-01", SP56771 = "TCGA-18-3408-01", 
  SP123995 = "TCGA-KN-8434-01", SP57941 = "TCGA-66-2756-01", SP39907 = "TCGA-B2-4099-01", 
  SP25332 = "TCGA-06-1086-01", SP80213 = "TCGA-EJ-5503-01", SP4820 = "TCGA-AN-A04D-01", 
  SP4593 = "TCGA-BH-A0H0-01", SP40736 = "TCGA-CJ-4870-01", SP8085 = "TCGA-C8-A130-01", 
  SP86118 = "TCGA-FK-A3SD-01", SP98305 = "TCGA-DD-A3A7-01", SP34005 = "TCGA-CR-5250-01", 
  SP47990 = "TCGA-FG-8182-01", SP11235 = "TCGA-A2-A0D4-01", SP1003 = "TCGA-GD-A2C5-01", 
  SP89687 = "TCGA-AP-A0L9-01", SP50611 = "TCGA-44-2659-01", SP49247 = "TCGA-BC-A217-01", 
  SP84998 = "TCGA-BR-A4J4-01", SP42154 = "TCGA-CW-6087-01", SP6730 = "TCGA-A8-A08S-01", 
  SP28335 = "TCGA-06-0686-01", SP59938 = "TCGA-24-1103-01", SP123846 = "TCGA-KL-8328-01", 
  SP19670 = "TCGA-A6-2680-01", SP10635 = "TCGA-E2-A156-01", SP86306 = "TCGA-BJ-A0Z2-01", 
  SP97249 = "TCGA-A4-A57E-01", SP17905 = "TCGA-CA-6717-01", SP87099 = "TCGA-EM-A2OV-01", 
  SP84062 = "TCGA-BR-6564-01", SP106638 = "TCGA-MH-A561-01", SP64976 = "TCGA-24-1419-01", 
  SP26709 = "TCGA-06-5411-01", SP57818 = "TCGA-52-7812-01", SP1677 = "TCGA-BT-A20T-01", 
  SP49328 = "TCGA-EP-A2KB-01", SP50485 = "TCGA-55-6984-01", SP5381 = "TCGA-BH-A0DT-01", 
  SP90269 = "TCGA-B5-A11I-01", SP16934 = "TCGA-AA-A01V-01", SP107595 = "TCGA-C5-A2LV-01", 
  SP97161 = "TCGA-HE-A5NH-01", SP64296 = "TCGA-23-1118-01", SP98092 = "TCGA-DD-A4NE-01", 
  SP124003 = "TCGA-KN-8437-01", SP123876 = "TCGA-KL-8343-01", SP56460 = "TCGA-18-3415-01", 
  SP120755 = "TCGA-MR-A520-01", SP10084 = "TCGA-E2-A109-01", SP121790 = "TCGA-DX-A23R-01", 
  SP65807 = "TCGA-04-1367-01", SP56474 = "TCGA-43-5670-01", SP9251 = "TCGA-A7-A26G-01", 
  SP49417 = "TCGA-G3-A25W-01", SP25350 = "TCGA-27-2528-01", SP51037 = "TCGA-05-5429-01", 
  SP82431 = "TCGA-EE-A2M5-06", SP56704 = "TCGA-60-2724-01", SP84982 = "TCGA-BR-4280-01", 
  SP25380 = "TCGA-06-0211-02", SP24565 = "TCGA-06-0171-02", SP39594 = "TCGA-B2-4102-01", 
  SP59803 = "TCGA-25-1634-01", SP48414 = "TCGA-DU-6401-01", SP123878 = "TCGA-KL-8344-01", 
  SP55142 = "TCGA-67-3771-01", SP49175 = "TCGA-CC-5260-01", SP20993 = "TCGA-AA-3956-01", 
  SP121859 = "TCGA-DX-A3U5-01", SP48135 = "TCGA-CS-5395-01", SP84790 = "TCGA-D7-5579-01", 
  SP121808 = "TCGA-DX-A240-01", SP30113 = "TCGA-BA-6869-01", SP87434 = "TCGA-DJ-A13R-01", 
  SP82459 = "TCGA-EE-A2A0-06", SP113197 = "TCGA-ER-A19T-06", SP41453 = "TCGA-CW-5585-01", 
  SP19983 = "TCGA-AA-A02O-01", SP122372 = "TCGA-FX-A3NJ-01", SP49223 = "TCGA-G3-A25T-01", 
  SP1305 = "TCGA-FD-A3N5-01", SP3016 = "TCGA-A7-A13D-01", SP82532 = "TCGA-FS-A1ZK-06", 
  SP25833 = "TCGA-06-0745-01", SP11045 = "TCGA-A2-A04X-01", SP80271 = "TCGA-AG-4008-01", 
  SP105673 = "TCGA-D7-8570-01", SP66159 = "TCGA-13-0906-01", SP18121 = "TCGA-AD-6964-01", 
  SP94933 = "TCGA-B5-A11G-01", SP43696 = "TCGA-B3-3925-01", SP56533 = "TCGA-85-8052-01", 
  SP107640 = "TCGA-C5-A1BJ-01", SP2826 = "TCGA-BH-A0B3-01", SP56730 = "TCGA-22-5485-01", 
  SP963 = "TCGA-DK-A1A5-01", SP63356 = "TCGA-13-1487-01", SP82636 = "TCGA-GN-A264-06", 
  SP1144 = "TCGA-BT-A2LA-01", SP59245 = "TCGA-66-2766-01", SP95126 = "TCGA-AP-A0LF-01", 
  SP49187 = "TCGA-G3-A25V-01", SP49157 = "TCGA-EP-A26S-01", SP53387 = "TCGA-55-6972-01", 
  SP11171 = "TCGA-AQ-A04J-01", SP31030 = "TCGA-CV-5442-01", SP68348 = "TCGA-24-2024-01", 
  SP123842 = "TCGA-KL-8326-01", SP5784 = "TCGA-E2-A1LL-01", SP85864 = "TCGA-FE-A233-01", 
  SP3368 = "TCGA-B6-A0WX-01", SP49379 = "TCGA-CC-5262-01", SP1132 = "TCGA-DK-A3IL-01", 
  SP97113 = "TCGA-A4-A4ZT-01", SP49541 = "TCGA-EP-A2KA-01", SP47652 = "TCGA-DU-7009-01", 
  SP57024 = "TCGA-90-7767-01", SP98164 = "TCGA-FV-A2QQ-01", SP96122 = "TCGA-A6-A566-01", 
  SP95646 = "TCGA-B5-A1MY-01", SP54745 = "TCGA-49-6742-01", SP25905 = "TCGA-14-1823-01", 
  SP105807 = "TCGA-EM-A3FR-01", SP57933 = "TCGA-21-5782-01", SP56754 = "TCGA-34-5240-01", 
  SP49481 = "TCGA-DD-A1ED-01", SP121781 = "TCGA-DX-A1L2-01", SP119755 = "TCGA-QG-A5YW-01", 
  SP98078 = "TCGA-DD-A3A8-01", SP8564 = "TCGA-A2-A0YG-01", SP122489 = "TCGA-HB-A5W3-01", 
  SP26439 = "TCGA-02-2483-01", SP5844 = "TCGA-E2-A1LK-01", SP123900 = "TCGA-KM-8477-01", 
  SP33544 = "TCGA-CV-6956-01", SP56827 = "TCGA-21-1076-01", SP82796 = "TCGA-DA-A3F3-06", 
  SP106677 = "TCGA-DD-A1EG-01", SP34452 = "TCGA-CW-6093-01", SP33688 = "TCGA-CR-6472-01", 
  SP122560 = "TCGA-IE-A4EI-01", SP52284 = "TCGA-75-5147-01", SP83482 = "TCGA-EE-A2GT-06", 
  SP53618 = "TCGA-64-1678-01", SP52232 = "TCGA-73-4659-01", SP24236 = "TCGA-14-2554-01", 
  SP89291 = "TCGA-EL-A3T0-01", SP88050 = "TCGA-L6-A4ET-01", SP6519 = "TCGA-A1-A0SM-01", 
  SP42829 = "TCGA-BP-4326-01", SP104530 = "TCGA-ER-A19E-06", SP30907 = "TCGA-HD-7753-01", 
  SP62450 = "TCGA-04-1349-01", SP96126 = "TCGA-NH-A50V-01", SP89443 = "TCGA-B5-A0K8-01", 
  SP89651 = "TCGA-AX-A06B-01", SP88593 = "TCGA-DJ-A2Q2-01", SP122590 = "TCGA-IE-A4EK-01", 
  SP23925 = "TCGA-19-2629-01", SP85339 = "TCGA-IN-7806-01", SP79939 = "TCGA-CH-5763-01", 
  SP123978 = "TCGA-KN-8429-01", SP67428 = "TCGA-25-2391-01", SP109649 = "TCGA-HC-8258-01", 
  SP6766 = "TCGA-AO-A12H-01", SP50412 = "TCGA-55-6986-01", SP81172 = "TCGA-AF-2689-01", 
  SP5448 = "TCGA-A8-A07B-01", SP88158 = "TCGA-FE-A3PD-01", SP123897 = "TCGA-KM-8476-01", 
  SP87534 = "TCGA-EM-A2CN-01", SP84719 = "TCGA-CG-4474-01", SP61713 = "TCGA-24-1614-01", 
  SP106743 = "TCGA-DD-A4NB-01", SP12049 = "TCGA-GI-A2C9-01", SP88776 = "TCGA-ET-A3DV-01", 
  SP98359 = "TCGA-FV-A3I0-01", SP98192 = "TCGA-DD-A4NA-01", SP97124 = "TCGA-A4-A48D-01", 
  SP24815 = "TCGA-41-5651-01", SP19606 = "TCGA-AA-3666-01", SP105425 = "TCGA-HU-8245-01", 
  SP84056 = "TCGA-HF-7136-01", SP85725 = "TCGA-EL-A3H1-01", SP98053 = "TCGA-FV-A23B-01", 
  SP96118 = "TCGA-A6-A565-01", SP43808 = "TCGA-AL-3472-01", SP121763 = "TCGA-DX-A1KW-01", 
  SP86660 = "TCGA-E8-A418-01", SP124013 = "TCGA-KO-8405-01", SP43822 = "TCGA-AL-3466-01", 
  SP37970 = "TCGA-CZ-5453-01", SP57084 = "TCGA-85-8277-01", SP87903 = "TCGA-DJ-A4UT-01", 
  SP88098 = "TCGA-BJ-A191-01", SP82429 = "TCGA-DA-A1I8-06", SP124021 = "TCGA-KO-8407-01", 
  SP82644 = "TCGA-EE-A2MI-06", SP27957 = "TCGA-32-1970-01", SP94060 = "TCGA-D1-A16G-01", 
  SP21528 = "TCGA-AA-3534-01", SP82988 = "TCGA-FS-A1ZP-06", SP85487 = "TCGA-EL-A4K6-01", 
  SP106631 = "TCGA-ED-A4XI-01", SP50321 = "TCGA-50-5932-01", SP2714 = "TCGA-D8-A27H-01", 
  SP68725 = "TCGA-09-2045-01", SP121811 = "TCGA-DX-A2IZ-01", SP48888 = "TCGA-HW-7486-01", 
  SP48426 = "TCGA-DB-5278-01", SP34186 = "TCGA-A3-3308-01", SP53548 = "TCGA-67-6215-01", 
  SP57538 = "TCGA-98-8022-01", SP105261 = "TCGA-BR-8486-01", SP122714 = "TCGA-IW-A3M4-01", 
  SP49551 = "TCGA-CC-A1HT-01", SP49651 = "TCGA-G3-A25Y-01", SP43514 = "TCGA-B9-4117-01", 
  SP30185 = "TCGA-CR-7385-01", SP2801 = "TCGA-AR-A24Z-01", SP82461 = "TCGA-EB-A24D-01", 
  SP8891 = "TCGA-A8-A075-01", SP21400 = "TCGA-AA-A00N-01", SP9481 = "TCGA-E2-A15H-01", 
  SP85379 = "TCGA-CG-5724-01", SP5666 = "TCGA-A8-A09X-01", SP27603 = "TCGA-19-2620-01", 
  SP92787 = "TCGA-BG-A18C-01", SP98265 = "TCGA-FV-A496-01", SP5808 = "TCGA-A2-A25B-01", 
  SP98297 = "TCGA-G3-A5SL-01", SP57735 = "TCGA-66-2744-01", SP123852 = "TCGA-KL-8331-01", 
  SP105086 = "TCGA-HU-A4H0-01", SP92723 = "TCGA-BS-A0TD-01", SP33837 = "TCGA-CR-6482-01", 
  SP123998 = "TCGA-KN-8435-01", SP27825 = "TCGA-16-1063-01", SP57669 = "TCGA-68-8250-01", 
  SP6429 = "TCGA-A2-A0D0-01", SP121870 = "TCGA-DX-A3U8-01", SP121852 = "TCGA-DX-A3M1-01", 
  SP79988 = "TCGA-HC-7740-01", SP92364 = "TCGA-AP-A054-01", SP5017 = "TCGA-AN-A0G0-01", 
  SP8229 = "TCGA-E2-A14X-01", SP8157 = "TCGA-GM-A2DF-01", SP5279 = "TCGA-E2-A14P-01", 
  SP23078 = "TCGA-AA-A01T-01", SP82900 = "TCGA-DA-A1HV-06", SP80183 = "TCGA-G9-6365-01", 
  SP56566 = "TCGA-94-7943-01", SP105759 = "TCGA-BJ-A45K-01", SP97258 = "TCGA-GL-A4EM-01", 
  SP81840 = "TCGA-AG-3727-01", SP86775 = "TCGA-EM-A3FQ-01", SP88322 = "TCGA-DJ-A1QL-01", 
  SP56569 = "TCGA-60-2713-01", SP123010 = "TCGA-MO-A47R-01", SP67609 = "TCGA-09-2050-01", 
  SP58342 = "TCGA-66-2759-01", SP8660 = "TCGA-A7-A0CE-01", SP2781 = "TCGA-AR-A2LK-01", 
  SP13084 = "TCGA-EX-A1H5-01", SP24135 = "TCGA-14-1034-02", SP8532 = "TCGA-E2-A1LG-01", 
  SP1712 = "TCGA-DK-A1AG-01", SP122361 = "TCGA-FX-A2QS-01", SP80657 = "TCGA-AG-3901-01", 
  SP121837 = "TCGA-DX-A3LU-01", SP107603 = "TCGA-C5-A1BF-01", SP58349 = "TCGA-66-2793-01", 
  SP85582 = "TCGA-DE-A0Y3-01", SP84392 = "TCGA-CG-4442-01", SP31790 = "TCGA-CV-5973-01", 
  SP98060 = "TCGA-DD-A3A6-01", SP98289 = "TCGA-ED-A459-01", SP109957 = "TCGA-EK-A2RM-01", 
  SP1114 = "TCGA-DK-A1AA-01", SP62110 = "TCGA-10-0938-01", SP64546 = "TCGA-36-1574-01", 
  SP30143 = "TCGA-BB-4225-01", SP115498 = "TCGA-BW-A5NO-01", SP81494 = "TCGA-AF-2691-01", 
  SP123967 = "TCGA-KN-8425-01", SP37903 = "TCGA-CJ-4918-01", SP30071 = "TCGA-BA-5556-01", 
  SP95406 = "TCGA-BS-A0U9-01", SP91746 = "TCGA-AP-A05D-01", SP26499 = "TCGA-19-1389-02", 
  SP30493 = "TCGA-CR-7404-01", SP109953 = "TCGA-EK-A2RL-01", SP11878 = "TCGA-BH-A1FC-01", 
  SP63546 = "TCGA-24-1466-01", SP56644 = "TCGA-68-7755-01", SP98096 = "TCGA-FV-A495-01", 
  SP35849 = "TCGA-BP-4327-01", SP48480 = "TCGA-DU-5874-01", SP17430 = "TCGA-AA-3685-01", 
  SP123988 = "TCGA-KN-8432-01", SP48263 = "TCGA-HT-7695-01", SP60322 = "TCGA-13-0890-01", 
  SP9433 = "TCGA-E2-A152-01", SP32014 = "TCGA-CN-4737-01", SP79971 = "TCGA-CH-5771-01", 
  SP47628 = "TCGA-HT-7602-01", SP87446 = "TCGA-FE-A22Z-01", SP97194 = "TCGA-B9-A44B-01", 
  SP82399 = "TCGA-DA-A1HY-06", SP127634 = "TCGA-FF-8042-01", SP57629 = "TCGA-56-7582-01", 
  SP60243 = "TCGA-10-0937-01", SP96110 = "TCGA-NH-A50T-01", SP114020 = "TCGA-DG-A2KJ-01", 
  SP47808 = "TCGA-E1-5318-01", SP1431 = "TCGA-BT-A20Q-01", SP56941 = "TCGA-21-1078-01", 
  SP92332 = "TCGA-BS-A0V8-01", SP82756 = "TCGA-ER-A19L-06", SP81711 = "TCGA-AG-4007-01", 
  SP91730 = "TCGA-B5-A0JN-01", SP2881 = "TCGA-AC-A2BK-01", SP82836 = "TCGA-FS-A1ZD-06", 
  SP7378 = "TCGA-AO-A03L-01", SP109470 = "TCGA-B1-A47M-01", SP63966 = "TCGA-24-1558-01", 
  SP9816 = "TCGA-BH-A0AV-01", SP49114 = "TCGA-DD-A1EB-01", SP87337 = "TCGA-FK-A3SE-01", 
  SP25518 = "TCGA-06-0744-01", SP975 = "TCGA-DK-A1A7-01", SP92947 = "TCGA-AX-A1CI-01", 
  SP105018 = "TCGA-CD-8529-01", SP57619 = "TCGA-22-1016-01", SP123984 = "TCGA-KN-8431-01", 
  SP18946 = "TCGA-AA-3977-01", SP1059 = "TCGA-H4-A2HQ-01", SP123858 = "TCGA-KL-8334-01", 
  SP96540 = "TCGA-C5-A0TN-01", SP84858 = "TCGA-D7-6815-01", SP91666 = "TCGA-DI-A1NN-01", 
  SP60362 = "TCGA-13-0912-01", SP58245 = "TCGA-66-2789-01", SP98090 = "TCGA-DD-A3A9-01", 
  SP121824 = "TCGA-DX-A2J4-01", SP80217 = "TCGA-G9-6370-01", SP57251 = "TCGA-60-2722-01", 
  SP28659 = "TCGA-06-0125-02", SP23622 = "TCGA-14-0786-01", SP90245 = "TCGA-AJ-A23M-01", 
  SP4265 = "TCGA-B6-A0RU-01", SP1419 = "TCGA-BT-A20V-01", SP82417 = "TCGA-DA-A3F8-06", 
  SP123890 = "TCGA-KM-8441-01", SP2793 = "TCGA-AO-A124-01", SP6825 = "TCGA-A8-A09I-01", 
  SP27939 = "TCGA-27-1831-01", SP12856 = "TCGA-EW-A1J5-01", SP49286 = "TCGA-FV-A3I1-01", 
  SP90893 = "TCGA-AP-A05A-01", SP30332 = "TCGA-CV-5432-01", SP69077 = "TCGA-23-1124-01", 
  SP123950 = "TCGA-KN-8418-01", SP36586 = "TCGA-CJ-5681-01", SP34431 = "TCGA-CJ-4899-01", 
  SP49469 = "TCGA-DD-A1EJ-01", SP4523 = "TCGA-BH-A18R-01", SP55004 = "TCGA-44-6148-01", 
  SP123886 = "TCGA-KM-8439-01", SP121841 = "TCGA-DX-A3LW-01", SP9930 = "TCGA-EW-A1PB-01", 
  SP124033 = "TCGA-KO-8411-01", SP79968 = "TCGA-CH-5788-01", SP21057 = "TCGA-AA-3529-01", 
  SP121865 = "TCGA-DX-A3U7-01", SP80160 = "TCGA-EJ-5506-01", SP123892 = "TCGA-KM-8442-01", 
  SP123840 = "TCGA-KL-8325-01", SP114032 = "TCGA-EK-A2R9-01", SP56079 = "TCGA-73-4666-01", 
  SP82471 = "TCGA-D3-A1Q1-06", SP67346 = "TCGA-04-1542-01", SP89519 = "TCGA-AP-A0LO-01", 
  SP22750 = "TCGA-AA-3514-01", SP31814 = "TCGA-CN-5365-01", SP66687 = "TCGA-24-1548-01", 
  SP127636 = "TCGA-FF-8061-01", SP84086 = "TCGA-D7-6519-01", SP105213 = "TCGA-HU-8608-01", 
  SP121831 = "TCGA-DX-A3LT-01", SP49205 = "TCGA-CC-5261-01", SP29987 = "TCGA-CR-6491-01", 
  SP123969 = "TCGA-KN-8426-01", SP97278 = "TCGA-HE-A5NL-01", SP105577 = "TCGA-D7-A4YX-01", 
  SP37432 = "TCGA-BP-5010-01", SP35989 = "TCGA-BP-4807-01", SP39102 = "TCGA-AK-3455-01", 
  SP123953 = "TCGA-KN-8419-01", SP19750 = "TCGA-AA-A01S-01", SP107624 = "TCGA-C5-A1MQ-01", 
  SP23754 = "TCGA-26-1438-01", SP5820 = "TCGA-A8-A07I-01", SP56502 = "TCGA-43-3920-01", 
  SP48850 = "TCGA-IK-7675-01", SP30083 = "TCGA-CV-7090-01", SP10470 = "TCGA-A2-A04P-01", 
  SP9648 = "TCGA-B6-A0I6-01", SP29697 = "TCGA-19-2624-01", SP58668 = "TCGA-60-2726-01", 
  SP86929 = "TCGA-EL-A3MY-01", SP109544 = "TCGA-HE-A5NJ-01", SP80423 = "TCGA-AG-4015-01", 
  SP4535 = "TCGA-BH-A0H6-01", SP127628 = "TCGA-FF-8046-01", SP11808 = "TCGA-AO-A0JM-01", 
  SP34493 = "TCGA-BP-4968-01", SP80244 = "TCGA-G9-6336-01", SP32742 = "TCGA-CN-5374-01", 
  SP58101 = "TCGA-22-5492-01", SP1009 = "TCGA-DK-A1AE-01", SP1724 = "TCGA-CF-A27C-01", 
  SP58991 = "TCGA-77-7139-01", SP86130 = "TCGA-DJ-A3US-01", SP97104 = "TCGA-IA-A40Y-01", 
  SP86425 = "TCGA-BJ-A0ZB-01", SP94917 = "TCGA-BS-A0TE-01", SP122676 = "TCGA-IS-A3K7-01", 
  SP82435 = "TCGA-D3-A3MO-06", SP104056 = "TCGA-EE-A3J5-06", SP92659 = "TCGA-BS-A0TC-01", 
  SP85495 = "TCGA-EL-A3CV-01", SP106577 = "TCGA-MH-A55Z-01", SP27339 = "TCGA-06-0157-01", 
  SP63716 = "TCGA-61-2000-01", SP96163 = "TCGA-GM-A3XL-01", SP17016 = "TCGA-AA-A02Y-01", 
  SP29559 = "TCGA-06-0221-02", SP94332 = "TCGA-A5-A0GJ-01", SP107575 = "TCGA-EK-A2PK-01", 
  SP109457 = "TCGA-AL-A5DJ-01", SP36036 = "TCGA-CZ-5454-01", SP115162 = "TCGA-BA-A4IH-01", 
  SP29331 = "TCGA-06-0210-02", SP35951 = "TCGA-B0-5695-01", SP79998 = "TCGA-CH-5789-01", 
  SP16958 = "TCGA-AA-3664-01", SP96136 = "TCGA-BT-A3PJ-01", SP11948 = "TCGA-A2-A04T-01", 
  SP94661 = "TCGA-AP-A053-01", SP56607 = "TCGA-96-7545-01", SP56168 = "TCGA-55-6982-01", 
  SP17329 = "TCGA-AA-A01X-01", SP26475 = "TCGA-19-5960-01", SP29940 = "TCGA-CR-7382-01", 
  SP58326 = "TCGA-33-4586-01", SP85787 = "TCGA-EM-A2CP-01", SP98327 = "TCGA-FV-A3R2-01", 
  SP107650 = "TCGA-C5-A1M9-01", SP49449 = "TCGA-DD-A1EH-01", SP123955 = "TCGA-KN-8421-01", 
  SP48189 = "TCGA-DU-7301-01", SP80754 = "TCGA-AG-3885-01", SP54363 = "TCGA-05-4395-01", 
  SP39248 = "TCGA-DV-5566-01", SP109478 = "TCGA-B1-A47N-01", SP123882 = "TCGA-KL-8346-01", 
  SP53810 = "TCGA-50-5930-01", SP127640 = "TCGA-FF-8062-01", SP43688 = "TCGA-B9-4617-01", 
  SP50518 = "TCGA-50-6591-01", SP65376 = "TCGA-24-1557-01", SP123894 = "TCGA-KM-8443-01", 
  SP49391 = "TCGA-DD-A1EI-01", SP49531 = "TCGA-G3-A25S-01", SP23775 = "TCGA-06-0190-02", 
  SP30843 = "TCGA-CR-5249-01", SP32222 = "TCGA-CR-7391-01", SP97243 = "TCGA-HE-A5NF-01", 
  SP48073 = "TCGA-HW-7487-01", SP85431 = "TCGA-CD-5802-01", SP49322 = "TCGA-BC-A216-01", 
  SP31190 = "TCGA-CQ-6228-01", SP83083 = "TCGA-ER-A19D-06", SP28275 = "TCGA-06-2557-01", 
  SP49433 = "TCGA-FV-A3R3-01", SP124017 = "TCGA-KO-8406-01", SP109384 = "TCGA-HP-A5MZ-01", 
  SP56303 = "TCGA-38-4628-01", SP57267 = "TCGA-37-4135-01", SP97269 = "TCGA-GL-A59R-01", 
  SP121761 = "TCGA-DX-A1KU-01", SP123856 = "TCGA-KL-8333-01", SP84962 = "TCGA-BR-6456-01", 
  SP4557 = "TCGA-B6-A0X5-01", SP967 = "TCGA-CF-A3MF-01", SP96147 = "TCGA-A2-A3XX-01", 
  SP48008 = "TCGA-FG-5964-01", SP61343 = "TCGA-04-1514-01", SP123854 = "TCGA-KL-8332-01", 
  SP93652 = "TCGA-AP-A052-01", SP30077 = "TCGA-BA-6872-01", SP105006 = "TCGA-HU-A4G6-01", 
  SP86836 = "TCGA-E8-A416-01", SP84439 = "TCGA-BR-6452-01", SP60001 = "TCGA-04-1331-01", 
  SP113926 = "TCGA-B1-A47O-01", SP122412 = "TCGA-FX-A48G-01", SP121816 = "TCGA-DX-A2J0-01", 
  SP90125 = "TCGA-AX-A2H5-01", SP10150 = "TCGA-A8-A092-01", SP55711 = "TCGA-50-6597-01", 
  SP81440 = "TCGA-AG-A032-01", SP1781 = "TCGA-BT-A3PH-01", SP122392 = "TCGA-FX-A3RE-01", 
  SP123888 = "TCGA-KM-8440-01", SP106560 = "TCGA-MH-A560-01", SP34246 = "TCGA-AK-3428-01", 
  SP98065 = "TCGA-FV-A4ZQ-01", SP82433 = "TCGA-ER-A19J-06", SP82780 = "TCGA-EE-A185-06", 
  SP96116 = "TCGA-AD-A5EK-01", SP89245 = "TCGA-EL-A3CX-01", SP83019 = "TCGA-DA-A1I0-06", 
  SP28581 = "TCGA-26-5132-01", SP123870 = "TCGA-KL-8340-01", SP94540 = "TCGA-B5-A0K3-01", 
  SP57066 = "TCGA-66-2795-01", SP48504 = "TCGA-E1-5319-01", SP49591 = "TCGA-DD-A1EL-01", 
  SP57901 = "TCGA-18-4721-01", SP80216 = "TCGA-EJ-7791-01", SP32694 = "TCGA-DQ-5625-01", 
  SP96129 = "TCGA-A6-A56B-01", SP21193 = "TCGA-A6-2681-01", SP48534 = "TCGA-CS-6668-01", 
  SP35412 = "TCGA-B0-5094-01", SP83099 = "TCGA-DA-A1I2-06", SP80950 = "TCGA-AG-3582-01", 
  SP1174 = "TCGA-FD-A3N6-01", SP33774 = "TCGA-CV-5443-01", SP4472 = "TCGA-C8-A12L-01", 
  SP122634 = "TCGA-IF-A4AJ-01", SP55509 = "TCGA-55-7281-01", SP57651 = "TCGA-60-2711-01", 
  SP23639 = "TCGA-27-2523-01", SP123958 = "TCGA-KN-8422-01", SP105159 = "TCGA-BR-8381-01", 
  SP59860 = "TCGA-25-2400-01", SP86989 = "TCGA-EL-A3T9-01", SP49385 = "TCGA-ES-A2HT-01", 
  SP123844 = "TCGA-KL-8327-01", SP43792 = "TCGA-B9-4114-01", SP6813 = "TCGA-C8-A12Q-01", 
  SP85230 = "TCGA-CG-4443-01", SP122702 = "TCGA-IS-A3KA-01", SP96511 = "TCGA-A2-A3Y0-01", 
  SP18787 = "TCGA-AA-3994-01", SP39298 = "TCGA-BP-5168-01", SP98313 = "TCGA-DD-A4ND-01", 
  SP55309 = "TCGA-75-6203-01", SP5636 = "TCGA-E2-A15E-06", SP38044 = "TCGA-CJ-5682-01", 
  SP38271 = "TCGA-BP-4756-01", SP5559 = "TCGA-A8-A094-01", SP121861 = "TCGA-DX-A3U6-01", 
  SP88757 = "TCGA-EM-A3AQ-01", SP84922 = "TCGA-D7-6527-01", SP49119 = "TCGA-ES-A2HS-01", 
  SP80042 = "TCGA-HC-7737-01", SP52607 = "TCGA-64-1680-01", SP30093 = "TCGA-BA-6873-01", 
  SP79907 = "TCGA-HC-7075-01", SP3631 = "TCGA-AR-A0TX-01", SP106602 = "TCGA-MH-A562-01", 
  SP43201 = "TCGA-CJ-4878-01", SP52667 = "TCGA-05-4420-01", SP33094 = "TCGA-CV-5431-01", 
  SP95222 = "TCGA-AP-A0LI-01", SP5473 = "TCGA-EW-A1PH-01", SP121847 = "TCGA-DX-A3LY-01", 
  SP953 = "TCGA-DK-A1A6-01", SP30803 = "TCGA-CV-6961-01", SP13242 = "TCGA-DS-A0VL-01", 
  SP97145 = "TCGA-IA-A40X-01", SP79958 = "TCGA-CH-5750-01", SP80205 = "TCGA-HI-7169-01", 
  SP83146 = "TCGA-GN-A26C-01", SP30675 = "TCGA-BA-4076-01", SP109301 = "TCGA-BW-A5NQ-01", 
  SP96114 = "TCGA-QG-A5YV-01", SP31854 = "TCGA-CX-7086-01", SP8394 = "TCGA-A2-A0D1-01", 
  SP6115 = "TCGA-EW-A1P8-01", SP43664 = "TCGA-B9-4116-01", SP1377 = "TCGA-BT-A20P-01", 
  SP85733 = "TCGA-FY-A2QD-01", SP89389 = "TCGA-BK-A0CC-01", SP7291 = "TCGA-AN-A0XR-01", 
  SP57450 = "TCGA-22-5477-01", SP2997 = "TCGA-BH-A0E0-01", SP127638 = "TCGA-FF-8043-01", 
  SP123872 = "TCGA-KL-8341-01", SP57846 = "TCGA-34-2600-01", SP60610 = "TCGA-24-1562-01", 
  SP85130 = "TCGA-F1-6875-01", SP104984 = "TCGA-BR-8682-01", SP8831 = "TCGA-AR-A256-01", 
  SP90503 = "TCGA-EO-A1Y8-01", SP24391 = "TCGA-06-0155-01", SP92460 = "TCGA-A5-A0G9-01", 
  SP85222 = "TCGA-CD-5799-01", SP92931 = "TCGA-EY-A1GW-01", SP2731 = "TCGA-BH-A0DG-01", 
  SP28041 = "TCGA-06-0145-01", SP82103 = "TCGA-AG-3574-01", SP91265 = "TCGA-BK-A139-01", 
  SP105375 = "TCGA-BR-8690-01", SP6223 = "TCGA-BH-A0WA-01", SP58612 = "TCGA-60-2719-01", 
  SP31606 = "TCGA-CR-6467-01", SP30213 = "TCGA-CR-6480-01", SP56821 = "TCGA-43-3394-01", 
  SP43532 = "TCGA-B9-4113-01", SP115830 = "TCGA-PD-A5DF-01", SP1365 = "TCGA-C4-A0F7-01", 
  SP19582 = "TCGA-A6-3807-01", SP123874 = "TCGA-KL-8342-01", SP7171 = "TCGA-E2-A15K-06", 
  SP34191 = "TCGA-B2-4101-01", SP80014 = "TCGA-G9-7522-01", SP10563 = "TCGA-A8-A08B-01", 
  SP40047 = "TCGA-A3-3372-01", SP114016 = "TCGA-C5-A1ML-01", SP8987 = "TCGA-AO-A0J6-01", 
  SP89909 = "TCGA-A5-A0GG-01", SP90629 = "TCGA-E6-A1LZ-01", SP49229 = "TCGA-BC-A10Q-01", 
  SP5980 = "TCGA-AN-A0AT-01", SP12186 = "TCGA-AO-A03N-01", SP48010 = "TCGA-DU-6407-01", 
  SP25494 = "TCGA-14-1402-02", SP7456 = "TCGA-A2-A3KC-01", SP51446 = "TCGA-97-8171-01", 
  SP103866 = "TCGA-EE-A29B-06", SP39349 = "TCGA-CZ-5987-01", SP127632 = "TCGA-FF-8041-01", 
  SP50592 = "TCGA-49-4512-01", SP85737 = "TCGA-DJ-A2Q8-01", SP53073 = "TCGA-78-7535-01", 
  SP58882 = "TCGA-77-6843-01", SP123836 = "TCGA-KL-8323-01", SP81137 = "TCGA-AG-3593-01", 
  SP57189 = "TCGA-21-1083-01", SP37369 = "TCGA-B0-5693-01", SP10944 = "TCGA-BH-A18U-01", 
  SP38759 = "TCGA-BP-4781-01", SP83844 = "TCGA-DA-A1HW-06", SP17294 = "TCGA-A6-2683-01", 
  SP2766 = "TCGA-D8-A27F-01", SP85647 = "TCGA-DE-A2OL-01", SP26649 = "TCGA-06-5415-01", 
  SP95550 = "TCGA-D1-A1NU-01", SP36218 = "TCGA-AK-3454-01", SP6673 = "TCGA-A8-A08L-01", 
  SP83027 = "TCGA-FS-A1ZU-06", SP64036 = "TCGA-24-2290-01", SP32894 = "TCGA-CV-7100-01", 
  SP85287 = "TCGA-EQ-5647-01", SP60842 = "TCGA-25-1632-01", SP1491 = "TCGA-FT-A3EE-01"
)

# Utils -------------------------------------------------------------------

.check_install <- function(pkg, bioc = FALSE, ...) {
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) install_func(pkg, ...)
  message("Required package ", pkg, " has been installed.")
}
