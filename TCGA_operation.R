# 诗翔的TCGA数据操作函数

#---------------------------------------------------
# 过滤TCGA重名ID
#---------------------------------------------------
# Filter TCGA Replicate Samples
# Author: ShixiangWang <w_shixiang@163.com>
# ooooooo
# The filter rule following broad institute says:
#
# In many instances there is more than one aliquot for a given combination of individual, platform, and data type. However, only one aliquot may be ingested into Firehose. Therefore, a set of precedence rules are applied to select the most scientifically advantageous one among them. Two filters are applied to achieve this aim: an Analyte Replicate Filter and a Sort Replicate Filter.
# 
# Analyte Replicate Filter
# The following precedence rules are applied when the aliquots have differing analytes. For RNA aliquots, T analytes are dropped in preference to H and R analytes, since T is the inferior extraction protocol. If H and R are encountered, H is the chosen analyte. This is somewhat arbitrary and subject to change, since it is not clear at present whether H or R is the better protocol. If there are multiple aliquots associated with the chosen RNA analyte, the aliquot with the later plate number is chosen. For DNA aliquots, D analytes (native DNA) are preferred over G, W, or X (whole-genome amplified) analytes, unless the G, W, or X analyte sample has a higher plate number.
# 
# Sort Replicate Filter
# The following precedence rules are applied when the analyte filter still produces more than one sample. The sort filter chooses the aliquot with the highest lexicographical sort value, to ensure that the barcode with the highest portion and/or plate number is selected when all other barcode fields are identical.
# oooooo
# Ref Link: <https://confluence.broadinstitute.org/display/GDAC/FAQ#FAQ-sampleTypesQWhatTCGAsampletypesareFirehosepipelinesexecutedupon>
#-------------------------------------------------------------------

tcgaReplicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19), filter_FFPE=FALSE, full_barcode=FALSE){
    # basically, user provide tsb and analyte_target is fine. If you
    # want to filter FFPE samples, please set filter_FFPE and full_barcode
    # all to TRUE, and tsb must have nchar of 28
    
    analyte_target = match.arg(analyte_target)
    # Strings in R are largely lexicographic
    # see ??base::Comparison
    
    # filter FFPE samples
    # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
    if(full_barcode & filter_FFPE){
      ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
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
               "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
      
      tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
    }
    
    # find repeated samples
    sampleID = substr(tsb, start = 1, stop = 15)
    dp_samples = unique(sampleID[duplicated(sampleID)])
    
    if(length(dp_samples)==0){
        message("ooo Not find any duplicated barcodes, return original input..")
        tsb
    }else{
        uniq_tsb = tsb[! sampleID %in% dp_samples]
        dp_tsb = setdiff(tsb, uniq_tsb)
        
        add_tsb = c()
        
        # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
        # if analyte_target = "DNA"
        # analyte:  D > G,W,X
        if(analyte_target == "DNA"){
            for(x in dp_samples){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                analytes = substr(mulaliquots, 
                                  start = analyte_position,
                                  stop = analyte_position)
                if(any(analytes == "D") & !(all(analytes == "D"))){
                    aliquot = mulaliquots[which(analytes == "D")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
            }
        }else{
            for(x in dp_samples){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                analytes = substr(mulaliquots, 
                                  start = analyte_position,
                                  stop = analyte_position)
                if(any(analytes == "H") & !(all(analytes == "H"))){
                    aliquot = mulaliquots[which(analytes == "H")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }else if(any(analytes == "R") & !(all(analytes == "R"))){
                    aliquot = mulaliquots[which(analytes == "R")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }else if(any(analytes == "T") & !(all(analytes == "T"))){
                    aliquot = mulaliquots[which(analytes == "T")]
                    add_tsb = c(add_tsb, aliquot)
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
            }
        }
        # if analyte_target = "RNA"
        # analyte: H > R > T 
        # else{
        #     
        # }
        if(length(dp_tsb) == 0){
            message("ooo Filter barcodes successfully!")
            c(uniq_tsb, add_tsb)
        }else{
            # filter according to portion number
            sampleID_res = substr(dp_tsb, start=1, stop=15)
            dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
            
            for(x in dp_samples_res){
                mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                portion_codes = substr(mulaliquots,
                                       start = portion[1],
                                       stop = portion[2])
                portion_keep = sort(portion_codes, decreasing = decreasing)[1]
                if(!all(portion_codes == portion_keep)){
                    if(length(which(portion_codes == portion_keep)) == 1){
                        add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
                        dp_tsb = setdiff(dp_tsb, mulaliquots)
                    }else{
                        dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
                    }
                    
                }
            }
            
            if(length(dp_tsb)==0){
                message("ooo Filter barcodes successfully!")
                c(uniq_tsb, add_tsb)
            }else{
                # filter according to plate number
                sampleID_res = substr(dp_tsb, start=1, stop=15)
                dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
                for(x in dp_samples_res){
                    mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
                    plate_codes = substr(mulaliquots,
                                         start = plate[1],
                                         stop = plate[2])
                    plate_keep = sort(plate_codes, decreasing = decreasing)[1]
                    add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
                    dp_tsb = setdiff(dp_tsb, mulaliquots)
                }
                
                if(length(dp_tsb)==0){
                    message("ooo Filter barcodes successfully!")
                    c(uniq_tsb, add_tsb)
                }else{
                    message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
                    c(uniq_tsb, add_tsb)
                }
            }
        }
    }
}


#------------------------------------------------------------
# 通过manifest文件转换为TCGA ID
# Source code from biostar and mofied for easier use
tcgaTranslateID = function(file_ids, legacy = TRUE) {
    # file_ids: id of download files in TCGA manifest file
    # legacy: if TRUE, use hg19
    if(!(require(GenomicDataCommons) & require(magrittr))){
        install.packages(c("GenomicDataCommons", "magrittr"), dependencies = TRUE)
    }
    info = files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}


# manifest <- read.table("gdc_manifest_filePath", header = TRUE, stringsAsFactors = FALSE)
# file_uuids <- manifest$id
# head(file_uuids)
# res = TCGAtranslateID(file_uuids)
# head(res)
