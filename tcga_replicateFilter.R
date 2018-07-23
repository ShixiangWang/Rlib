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
# Test and Usage example include in gist <https://gist.github.com/ShixiangWang/33b2f9b49b77eaa8f773b428480f9101>


tcga_replicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19)){
    analyte_target = match.arg(analyte_target)
    # Strings in R are largely lexicographic
    # see ??base::Comparison
    
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
