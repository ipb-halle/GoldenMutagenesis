#' @import seqinr
#' @import stringr
#' @import methods
#' @importFrom stats dist
#' @importFrom utils read.csv
NULL

#' 
#' Order the replacement list
#'
#' @param replacements A list containing vectors which have a number at the first slot and the amino acid at the second one
#'
#' @return A sorted list of the input
#'
#' @examples
#' \dontrun{
#' replacements<-list(c(55, "V"), c(40, "A"))
#' GoldenMutagenesis::order_replacements(replacements)
#' }
order_replacements<-function(replacements){return(replacements[order(sapply(replacements, function(x){as.numeric(x[1])}))])}

remove_newline<-function(x){
  gsub("\r?\n|\r", "", x)
}



#' Select a Codon Usage Table
#' 
#' Get a list or an array of values for the selected Codon Usage Table
#'
#' @param name The filname of the codon usage table. You can use list_cu_table() to get the overview.
#' @param list A boolean parameter that decides if an array or a list is returned. The array output is requiered by the domesticate function
#'
#' @return An array or a list with values for the codons/amino acids.
#' @export
#'
#' @examples
#' list_cu_table()
#' \dontrun{
#' get_cu_table("e_coli_316407.csv")
#' }
#' 
get_cu_table<-function(name, list=T) {
  stopifnot(is.character(name))
  file<-system.file("cuf", name, package="GoldenMutagenesis")
  cuf<-read.csv(file)
  cuf[, "codon"]<-str_replace_all(cuf[,"codon"], "U", "T")
  cuf_vector<-as.vector(cuf$relative_frequency)
  names(cuf_vector)<-gsub("U","T",cuf$codon)
  cuf_list<-lapply(unique(cuf$amino_acid), function(x){a<-cuf[which(cuf$amino_acid==x),"relative_frequency"];a<-a*1000;names(a)<-str_to_lower(cuf[which(cuf$amino_acid==x),"codon"]);return(as.table(a))})
  names(cuf_list)<-unique(cuf$amino_acid)
  if(list==T) {
    return(cuf_list)
  }
  else{
    return(cuf_vector)
  }
}

make_fragments<-function(mutations, fsize, buffer=0, seq, start, distance=2){
  start=start
  end=NULL
  fsize<-ceiling(fsize/3)
  cm<-1
  fragments<-c()
  distance<-distance+1
  if(buffer > distance){
    stop("You can not try to use a buffer bigger than the allowed distance to integrate two mutations on the same primer!")
  }
  #Make fragments
  repeat{
    this_fragment<-fragment(start=start)
    endm<-0
    if(cm > length(mutations)){#This generates a last fragment without any mutations
      this_fragment@stop<-length(seq)
      fragments<-c(fragments, this_fragment)
      break
    }
    i<-cm
    repeat{
      if((mutations[i] + buffer - start) >= fsize){#Iterate until we have a mutation on the end of the fragment
          if(length(seq) - mutations[i] < fsize){# the last fragment would be too small 
             #next mutation should have at least 1 or buffer +1 difference
              if(((length(seq) - fsize) - start >= fsize) & !is.null(fragments)) {
                this_fragment@stop<-length(seq)-fsize
                fragments<-c(fragments, this_fragment)
                start=length(seq)-(fsize-1)
                this_fragment<-fragment(start=start)
                end=length(seq)
                endm=length(mutations)
              } else {
                end<-length(seq)
                endm=length(mutations)
              }
              break
            }
          if(i<length(mutations)){
            if((mutations[i+1] - mutations[i])<distance){#Integrate next mutation
              i<-i+1
              next
            }
          }
        endm<-i
        end<-mutations[i] + buffer
        if(i == 1) { #if there was more than one mutation, we are already safe
          if(length(seq)-end < fsize) {
            this_fragment@stop<-length(seq)-fsize-1-buffer
            fragments<-c(fragments, this_fragment)
            this_fragment@start<-length(seq)-fsize-buffer
            end<-length(seq)
            endm<-1
          }
        }
        #buffer and distance must be at least equal
        break
      }
      else{
        if(i==length(mutations)){#If there is a mutation at the end, but the end of the fragment has to be the end of the sequence
          endm<-i
          end<-length(seq)
          break
        }
      }
      i<-i+1
    }
    #Is the mutation at the end or start of a fragment?
    this_fragment@stop=end
    mid<-this_fragment@start+round((this_fragment@stop - this_fragment@start)/2)
    if(cm<=length(mutations)){
      for(j in cm:endm){
        if(mutations[j] < mid){
          this_fragment@start_mutation<-c(this_fragment@start_mutation, mutations[j])
        }
        else{
          this_fragment@stop_mutation<-c(this_fragment@stop_mutation, mutations[j])
        }
      }
    }
    fragments<-c(fragments, this_fragment)
    if(end==length(seq)){
      break
    }
    else{
      cm<-endm+1
    }
    start<-this_fragment@stop+1
  }
  #Optimize Fragments - Shift start and stop
  if(length(fragments)>1){
    for (k in 1:(length(fragments)-1)) {
      if(length(fragments[[k]]@stop_mutation)>0){
        stop_pos<-fragments[[k]]@stop_mutation[1]
        if(length(fragments[[k+1]]@start_mutation)>0){
          start_pos<-fragments[[k+1]]@start_mutation[length(fragments[[k+1]]@start_mutation)]
          mid_stop<-stop_pos+(round((start_pos-stop_pos)/2))-1
          newstop<-max(fragments[[k]]@stop, min(mid_stop, (fragments[[k+1]]@start_mutation[1]-1)))
          fragments[[k]]@stop<-newstop
          fragments[[k+1]]@start<-newstop+1
        }
        else{
          next
        }
      }
      else{
        next
        #This should never be the case
      }
    }
  }
  return(fragments)
}

# make_fragments<-function(positions_aa,codon_seq, fragment_start=1, msd=F, replacement_range=3, min_fragment=100, binding_max_length=9){
#   min_fragment<-ceiling(min_fragment/3)
#   replacement_distances<-as.matrix(dist(positions_aa))
#   #First calculate fragments
#   #then calculate primers 
#   fragments<-c()
#   i<-1
#   first<-T
#   #todo replace formular with binding_max_length
#   repeat {
#     #################################
#     #################################
#     #Creation of the first fragment/primer
#     if (first == T) {
#       #first replacement
#       #check if it is on the beginning of the first fragment
#       temp_fragment <- fragment(start = fragment_start)
#       if (positions_aa[i] <  (min_fragment+fragment_start - 1)) {
#         #if (positions_aa[i] <= fragment_start - 1 + replacement_range + 2) {
#         temp_fragment@start_mutation <- positions_aa[i]
#         if (length(positions_aa) == 1) {
#           temp_fragment@stop <- length(codon_seq)
#           fragments <- c(fragments, temp_fragment)
#           break
#         }
#       }
#       else {
#         # we will create a new fragment
#         #check if there is enough space to create a new fragment
#         if (length(positions_aa) == 1) {
#           temp_fragment@stop <-  positions_aa[i]
#           if(msd){
#             temp_fragment@stop<-temp_fragment@stop + 2
#           }
#           temp_fragment@stop_mutation <- positions_aa[i]
#           fragments <- c(fragments, temp_fragment)
#           temp_fragment <-
#             fragment(start = fragments[[1]]@stop + 1,
#                      stop = length(codon_seq))
#           fragments <- c(fragments, temp_fragment)
#           break
#         }
#         if (replacement_distances[i + 1, i] < 3) {
#           temp_fragment@stop <- positions_aa[i] - 1
#         }
#         else{
#           temp_fragment@stop <-  positions_aa[i]
#           if(msd){
#             temp_fragment@stop <- temp_fragment@stop + 2
#           }
#           temp_fragment@stop_mutation <- positions_aa[i]
#           i <- i + 1
#         }
#         fragments <- c(fragments, temp_fragment)
#       }
#       first<-F
#     }
#     #################################
#     #################################
#     #generic part for all fragments
#     mutations_in_fragment_range <-
#       as.numeric(which(
#         replacement_distances[, i] < (min_fragment) &
#           replacement_distances[, i] > 0
#       ))
#     mutations_in_fragment_range <-
#       mutations_in_fragment_range[which(mutations_in_fragment_range > i)]
#     ###distinguish between fragment start and end
#     ###if the minimal binding length is very high and the distance of the mutations is very short, you will get very long primers!
#     if (length(temp_fragment@stop) == 0) {
#       if (any(temp_fragment@start_mutation[length(temp_fragment@start_mutation)] == positions_aa[i])) {
#         if (length(mutations_in_fragment_range) != 0) {
#           #we won't create a new fragment until we have enough space between the mutations to create a forward and a reverse primer
#           for (mutation in mutations_in_fragment_range[1:(length(mutations_in_fragment_range) - 1)]) {
#             temp_fragment@start_mutation <-
#               c(temp_fragment@start_mutation, positions_aa[mutation])
#             i <-
#               mutations_in_fragment_range[length(mutations_in_fragment_range)]
#           }
#           rm(mutation)
#         }
#         else {
#           i <- i + 1
#         }
#       }
#       else {
#         #We can finish the fragment here
#         if (replacement_distances[i + 1, i] < 3) {
#           #Check if we have the next mutation on the same primer
#           temp_fragment@stop = positions_aa[i] - 1
#         }
#         else {
#           #End the fragment with the mutation
#           temp_fragment@stop <- positions_aa[i]
#           if(msd) {
#             temp_fragment@stop <- positions_aa[i] + 2
#           }
#           temp_fragment@stop_mutation <-
#             c(temp_fragment@stop_mutation, positions_aa[i])
#           i <- i + 1
#         }
#         fragments <- c(fragments, temp_fragment)
#       }
#     }  
#     else {
#       ###create a new fragment in forward direction
#       #we will always start with a mutation in forward direction (execpt the first one)
#       temp_new_fragment <- fragment(start = temp_fragment@stop + 1)
#       temp_old_fragment <- temp_fragment
#       temp_fragment <- temp_new_fragment
#       rm(temp_new_fragment)
#       #does the fragment need additional shifiting?
#       if (positions_aa[i] - temp_fragment@start < (min_fragment)) {
#         #Mutation(s) at the beginning of the fragment
#         if (length(temp_old_fragment@stop_mutation) == 0) {
#           #there is no mutation at the end, we can shift anyway
#           fragments[[length(fragments)]]@stop <- positions_aa[i]
#           if(msd){
#             fragments[[length(fragments)]]@stop <- positions_aa[i] + 2
#           }
#           fragments[[length(fragments)]]@stop_mutation<-positions_aa[i]
#           i<-i+1
#           temp_fragment<-fragments[[length(fragments)]]
#         }
#         else {
#           #Check the difference between the last and the current mutation
#           last_mutation <-
#             fragments[[length(fragments)]]@stop_mutation[length(fragments[[length(fragments)]]@stop_mutation)]
#           diff <- positions_aa[i] - last_mutation - 1
#           if (diff <= replacement_range) {
#             fragments[[length(fragments)]]@stop <- positions_aa[i]#edit here
#             if(msd){
#               fragments[[length(fragments)]]@stop <- positions_aa[i] + 2
#             }
#             fragments[[length(fragments)]]@stop_mutation <-
#               c(fragments[[length(fragments)]]@stop_mutation, positions_aa[i])
#             temp_fragment<-fragments[[length(fragments)]]
#             i <- i + 1
#           }
#           else{
#             fragments[[length(fragments)]]@stop <-
#               fragments[[length(fragments)]]@stop + replacement_range
#             temp_fragment@start <-
#               fragments[[length(fragments)]]@stop + 1
#             temp_fragment@start_mutation <-
#               c(temp_fragment@start_mutation, positions_aa[i])
#           }
#         }
#       }
#     }
#     #################################
#     #################################
#     #part for the last mutation and the last fragment
#     if (i == length(positions_aa)) {
#       if (any(temp_fragment@start_mutation[length(temp_fragment@start_mutation)] == positions_aa[i])) {
#         #we started a fragment with the mutation on the forward part
#         #just end it here
#         temp_fragment@stop <- length(codon_seq)
#         fragments <- c(fragments, temp_fragment)
#         rm(temp_fragment)
#         break
#       }
#       else{
#         if ((length(codon_seq) - positions_aa[i-1]) >= (min_fragment)) {
#           if(length(temp_fragment@stop)!=0) {
#             temp_new_fragment <- fragment(start = temp_fragment@stop + 1)
#             #temp_old_fragment <- temp_fragment
#             temp_fragment <- temp_new_fragment
#             rm(temp_new_fragment)
#           }
#           else {
#             if(length(temp_fragment@start)==0) {
#               stop("Internal error. Please give a bug report!")
#             }
#           }
#           if(positions_aa[i]-temp_fragment@start < binding_max_length+2 ) {
#             temp_fragment@start_mutation<-c(temp_fragment@start_mutation,positions_aa[i])
#             temp_fragment@stop<-length(codon_seq)
#             fragments <- c(fragments, temp_fragment)
#             break
#           }
#           else {
#             temp_fragment@stop_mutation<-positions_aa[i]
#             if ((length(codon_seq) - positions_aa[i]) > (min_fragment)) {
#               temp_fragment@stop<-temp_fragment@stop_mutation
#               if(msd){
#                 temp_fragment@stop<-temp_fragment@stop_mutation+2
#               }
#               fragments <- c(fragments, temp_fragment)
#               temp_fragment<-fragment(start = temp_fragment@stop + 1, stop = length(codon_seq))
#               fragments <- c(fragments, temp_fragment)
#               break
#             }
#             else{
#               temp_fragment@stop<-length(codon_seq)
#               fragments <- c(fragments, temp_fragment)
#               break
#             }
#           }
#         }
#         else {
#           temp_fragment@stop<-length(codon_seq)
#           temp_fragment@stop_mutation <-
#             c(temp_fragment@stop_mutation, positions_aa[i])
#           fragments[length(fragments)] <- temp_fragment
#           break
#         }
#       }
#     }
#   }
#   return(fragments)
# }

calculate_tm<-function(x, salt_concentration=50, primer_concentration=50, offset=0){
  oligo_sequence<-s2c(x)
  oligo_sequence<-oligo_sequence[offset:length(oligo_sequence)]
  #  Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
  counts<-count(s2c(x), wordsize=1, by=1, alphabet = c("A", "C", "G", "T"))
  tm<-100.5 + (41 * as.numeric(counts["G"] + counts["C"])/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) - (820/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) + 16.6*log10(salt_concentration/1000)
  return(as.numeric(tm))
}

calculate_UPAC<-function(x, func=calculate_tm, selection="max", temp=0, cuf="e_coli_316407.csv") {
  cuf_vector<-get_cu_table(cuf, list=F)
  x_matrix<-expand.grid(sapply(s2c(x), function(x){amb(x, forceToLower = TRUE, checkBase = TRUE, IUPAC = s2c("acgturymkswbdhvn"), u2t = TRUE)}))
  results<-apply(x_matrix, 1, function(x){func(str_to_upper(paste(x, collapse="")))})
  if(selection=="max") {
    candidates<-which(results==max(results))
  }
  if(selection=="min"){
    candidates<-which(results==min(results))
  }
  if(selection=="diff"){
    candidates<-which(abs(results-temp)==min(abs(results-temp)))
  }
  if(length(candidates)>1){
    #select sequence with highest probability
    sum_of_prob<-apply(x_matrix[candidates,], 1, function(x){c_x<-count(s2c(str_to_upper(paste(x, collapse=""))), wordsize = 3, by=1, alphabet = c("A", "C", "G", "T")); return(sum(c_x*cuf_vector[names(c_x)]))})
    candidate<-candidates[which(sum_of_prob==max(sum_of_prob))]
    return(str_to_upper(paste(apply(x_matrix[candidate,], 1, as.character), collapse="")))
  }
  else {
    return(str_to_upper(paste(apply(x_matrix[candidates[1],], 1, as.character), collapse="")))
  }
}

calculate_DeltaG<-function(x){
  g <- -5.0
  oligo_sequence<-s2c(x)
  counts<-count(s2c(x), wordsize=2, by=1, alphabet = c("A", "C", "G", "T"))
  g <-g + (as.numeric(counts["AA"])+as.numeric(counts["TT"]))*1.2
  g <-g + as.numeric(counts["AT"])*0.9
  g <-g + as.numeric(counts["TA"])*0.9
  g <-g + (as.numeric(counts["CA"])+as.numeric(counts["TG"]))*1.7
  g <-g + (as.numeric(counts["GT"])+as.numeric(counts["AC"]))*1.5
  g <-g + (as.numeric(counts["CT"])+as.numeric(counts["AG"]))*1.5
  g <-g + (as.numeric(counts["GA"])+as.numeric(counts["TC"]))*1.5
  g <-g + as.numeric(counts["CG"])*2.8
  g <-g + as.numeric(counts["GC"])*2.3
  g <-g + (as.numeric(counts["GG"])+as.numeric(counts["CC"]))*2.1
  return(g)
  
}

calculate_DeltaH<-function(x){
  h <- 0.0
  oligo_sequence<-s2c(x)
  counts<-count(s2c(x), wordsize=2, by=1, alphabet = c("A", "C", "G", "T"))
  h <-h + (as.numeric(counts["AA"])+as.numeric(counts["TT"]))*8
  h <-h + as.numeric(counts["AT"])*5.6
  h <-h + as.numeric(counts["TA"])*6.6
  h <-h + (as.numeric(counts["CA"])+as.numeric(counts["TG"]))*8.2
  h <-h + (as.numeric(counts["GT"])+as.numeric(counts["AC"]))*9.4
  h <-h + (as.numeric(counts["CT"])+as.numeric(counts["AG"]))*6.6
  h <-h + (as.numeric(counts["GA"])+as.numeric(counts["TC"]))*8.8
  h <-h + as.numeric(counts["CG"])*11.8
  h <-h + as.numeric(counts["GC"])*10.5
  h <-h + (as.numeric(counts["GG"])+as.numeric(counts["CC"]))*10.9
  return(h)
}

calculate_DeltaS<-function(x){
  s<-0.0
  oligo_sequence<-s2c(x)
  counts<-count(s2c(x), wordsize=2, by=1, alphabet = c("A", "C", "G", "T"))
  s <-s + (as.numeric(counts["AA"])+as.numeric(counts["TT"]))*21.9
  s <-s + as.numeric(counts["AT"])*15.2
  s <-s + as.numeric(counts["TA"])*18.4
  s <-s + (as.numeric(counts["CA"])+as.numeric(counts["TG"]))*21.0
  s <-s + (as.numeric(counts["GT"])+as.numeric(counts["AC"]))*25.5
  s <-s + (as.numeric(counts["CT"])+as.numeric(counts["AG"]))*16.4
  s <-s + (as.numeric(counts["GA"])+as.numeric(counts["TC"]))*23.5
  s <-s + as.numeric(counts["CG"])*29.0
  s <-s + as.numeric(counts["GC"])*26.4
  s <-s + (as.numeric(counts["GG"])+as.numeric(counts["CC"]))*28.4
  return(s)
}

#' Calculate melting temperature based on next neighbor calculation
#' 
#' The implementation is based on the explanations of \url{http://biotools.nubic.northwestern.edu/OligoCalc.html}.
#' 
#' More details at \url{https://doi.org/10.1093/nar/gkm234} 
#'
#' @param oligo_sequence A string containing an oligo sequence.
#' @param primer_concentration The concentration of the primer in nanomole [default: 50]
#' @param salt_concentration The concentration of Na+ in nanomole [default: 50]
#' @param offset You can skip a prefix of your oligo sequence with this parameter. The first n bases are not considered in the calculation. [default: 0] 
#' @return An array or a list with values for the codons/amino acids.
#' @return The melting temperature in \code{print('\u00B0')}C
#' 
#' @examples
#' \dontrun{
#' GoldenMutagenesis::calculate_tm_nnb("AAAAAATGGTGTGTGATGTGTCCCTCTATC")
#' }
#' 
calculate_tm_nnb<-function(oligo_sequence, primer_concentration=50, salt_concentration=50, offset=0){
  oligo_sequence_s2c<-s2c(oligo_sequence)
  oligo_sequence<-paste(oligo_sequence_s2c[offset:length(oligo_sequence_s2c)], collapse="")
  K<-1/(primer_concentration*1e-9) #Convert from nanomoles to moles
  R<-1.987
  RlnK<-R*log(K)
  result<-(1000*(calculate_DeltaH(oligo_sequence)-3.4)/(calculate_DeltaS(oligo_sequence)+RlnK)-272.9)
  result<-result+7.21*log(salt_concentration/1000)
  return(result)
}


setGeneric("sequence_length_temperature" , function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60,  gc_filter=F) {
  standardGeneric("sequence_length_temperature")
})

setMethod("sequence_length_temperature", signature(primer="Primer"),
          function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60, gc_filter=F){
            primer_seq_s2c<-s2c(primer@binding_sequence)
            temperatures<-list()
            names_i<-vector()
            sequences_i<-vector()
            for(i in (primer_min*3):length(primer_seq_s2c)){
              temperatures<-c(temperatures, temp_func(paste(primer_seq_s2c[1:i],collapse=""), offset=0))
              names_i<-c(names_i, i)
              sequences_i<-c(sequences_i, paste(primer_seq_s2c[1:i],collapse=""))
            }
            names(temperatures)<-sequences_i
            diff<-unlist(lapply(temperatures, function(x){abs(x-target_temp)}))
            #check for at least two A or T
            candidates_with_AT<-which(str_count(str_sub(names(diff), start=-5), "A|T")>=2 & str_count(str_sub(names(diff), start=-5), "A|T")<4)
            if(length(candidates_with_AT) == 0 || gc_filter==F) {
              candidate_binding_sequence<-names(diff[diff==min(diff)])
              if(gc_filter==T) {
                warning("The end (last five bases) of the binding sequence is not optimal. The primers are maybe inefficient.")
              }
            }
            else{
              diff_AT<-diff[candidates_with_AT]
              candidate_binding_sequence<-names(diff_AT[diff_AT==min(diff_AT)])
            }
            primer@binding_sequence<-candidate_binding_sequence
            primer@temperature<-as.numeric(temperatures[candidate_binding_sequence])
            primer@difference<-as.numeric(diff[candidate_binding_sequence])
            return(primer)
          }
)

setMethod("sequence_length_temperature", signature(primer="Primer_MSD"),
          function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60, gc_filter=F){
            callNextMethod()
          }
)

setMethod("sequence_length_temperature", signature(primer="Primer_SPM"),
          function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60, gc_filter=F){
            callNextMethod()
          }
          # function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60, gc_filter=F){
          #   primer_seq_s2c<-s2c(paste(primer@extra, primer@binding_sequence, sep=""))
          #   temperatures<-list()
          #   names_i<-vector()
          #   sequences_i<-vector()
          #   for(i in max((primer_min*3), nchar(primer@extra)):length(primer_seq_s2c)){
          #     temperatures<-c(temperatures, temp_func(paste(primer_seq_s2c[1:i],collapse=""), offset=0))
          #     names_i<-c(names_i, i)
          #     sequences_i<-c(sequences_i, paste(primer_seq_s2c[1:i],collapse=""))
          #   }
          #   names(temperatures)<-sequences_i
          #   diff<-unlist(lapply(temperatures, function(x){abs(x-target_temp)}))
          #   #check for at least two A or T
          #   candidates_with_AT<-which(str_count(str_sub(names(diff), start=-5), "A|T")>=2 & str_count(str_sub(names(diff), start=-5), "A|T")<4)
          #   if(length(candidates_with_AT) == 0 || gc_filter==F) {
          #     candidate_binding_sequence<-names(diff[diff==min(diff)])
          #     if(gc_filter==T) {
          #       warning("The end (last five bases) of the binding sequence is not optimal. The primers are maybe inefficient.")
          #     }
          #   }
          #   else{
          #     diff_AT<-diff[candidates_with_AT]
          #     candidate_binding_sequence<-names(diff_AT[diff_AT==min(diff_AT)])
          #   }
          #   primer@binding_sequence<-str_sub(candidate_binding_sequence, max(nchar(primer@extra)+1,0))
          #   primer@temperature<-as.numeric(temperatures[candidate_binding_sequence])
          #   primer@difference<-as.numeric(diff[candidate_binding_sequence])
          #   return(primer)
          # }
)

sequence_check<-function(input_sequence){
  input_sequence<-str_to_upper(input_sequence)
  input_sequence<-str_trim(input_sequence)
  if(nchar(input_sequence)%%3!=0) {
    stop(paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "))
  }
  
  if(str_detect(input_sequence, "^(A|C|G|T)+$") == F) {
    stop(paste("The sequence contains invalid characters that are not A|C|G|T."))
  }
  
  codon_seq<-splitseq(s2c(input_sequence))
  met<-which(str_detect(codon_seq, "ATG"))
  if(length(met) == 0) {
    stop("No Methionine in the provided sequence. Stopping here. Please check the provided sequence.")
  }
  
  if(min(met) != 1){
    warning(paste("No Methionine at first codon found! Please check the provided sequence! Took codon #", min(met), "as start.", sep=" "))
    codon_seq<-codon_seq[min(met):length(codon_seq)]
  } #else(codon_seq<-codon_seq[-1])
  
  stop<-which(str_detect(codon_seq, "(TAA)|(TGA)|(TAG)"))
  if(length(stop) == 0) {
    stop("No stop codon in the provided sequence. Stopping here. Please check the provided sequence!")
  }
  
  if(max(stop) != length(codon_seq)) {
    warning(paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "))  
    codon_seq<-codon_seq[1:max(stop)]
  }# else {
  #codon_seq <- codon_seq[-length(codon_seq)]
  #}
  return(codon_seq)
}

check_primer_overhangs<-function(primers, fragments, binding_min_length=4, target_temp=60, check_repetitive=T) {
  #ToDo: Add paramter for temperature calculation method
  overhangs<-sapply(primers, function(x){return(c(x[[1]]@overhang, x[[2]]@overhang))})
  duplicates<-table(overhangs)
  duplicates<-duplicates[names(duplicates)!="" & duplicates > 1]
  if(check_repetitive == T) {
    #Repetitive overhangs
    rep<-table(overhangs)
    rep<-rep[names(rep)!=""]
    rep_temp<-names(rep)
    rep<-str_count(names(rep), ("(^(A|T){4}$)|(^(G|C){4}$)"))
    names(rep)<-rep_temp
    rep<-rep[rep > 0]
    rm(rep_temp)
    bad_overhangs<-union(names(duplicates), names(rep))
  } else {
    bad_overhangs<-duplicates
    }
  if(length(bad_overhangs)==0) {
    return(primers)
  }
  bad_overhang<-bad_overhangs[1]
  primer_num<-which(overhangs==bad_overhang)
  primer_unlist<-unlist(primers)
  fragment_num<-ceiling(primer_num)/2
  primer_num2<-primer_num %% 2
  primer_num2[primer_num2==0]<-2
  primer_num2[primer_num2==1]<-1
  for(i in 1:length(primer_num)) {
    if(primer_num2[i]==1) {
      primer_fd_num<-primer_num[i]
      primer_rv_num<-primer_num[i] - 1
    } else {
      primer_fd_num<-primer_num[i] + 1
      primer_rv_num<-primer_num[i]
    }
    primer_fd<-primer_unlist[[primer_fd_num]]
    primer_rv<-primer_unlist[[primer_rv_num]]
    #we will move to the left direction 
    #check if primer_rv is an NDT primer
    if(class(primer_rv)=="Primer_MSD") {
      msd_mut<-sapply(c("NNN", "NNK", "NNS", "NDT", "DBK", "NRT"), FUN = function(x){paste(stringr::str_to_upper(rev(seqinr::comp(seqinr::s2c(x), ambiguous=T))), sep="", collapse="")}, USE.NAMES = F)
      if(str_sub(primer_rv@extra, 1, 3) %in% msd_mut) {
        if(i == length(primer_num)) {
          warning(paste("We can not fix overlaps or palyndromic sequences in the primers. Please consider a silent mutation at position ", fragments[[ceiling((primer_rv_num+1)/2)]]@start, ".", sep=""))
          return(primers)
        }
        else {
          next
        }
      }
      else{
        shift_base<-str_sub(primer_rv@extra, 1, 1)
        primer_rv@overhang<-paste(primer_rv@overhang,shift_base, sep="")
        primer_rv@extra<-str_sub(primer_rv@extra, 2)
        primer_rv@overhang<-str_sub(primer_rv@overhang, 2)
        primer_rv@temperature<-calculate_tm_nnb(primer_rv@binding_sequence, offset = 0)
        primer_rv@difference<-abs(primer_rv@temperature - primer_unlist[[primer_rv_num -1 ]]@temperature)
      }
    } else {
      if(nchar(primer_rv@binding_sequence) < 3 * binding_min_length) {
        if(i == length(primer_num)) {
          warning(paste("We can not fix overlaps or palyndromic sequences in the primers. Please consider a silent mutation at position ", fragments[[ceiling((primer_rv_num+1)/2)]]@start, ".", sep=""))
          return(primers)
        }
        else{
          next
        }
      }
      else{
        shift_base<-str_sub(primer_rv@binding_sequence, 1, 1)
        primer_rv@extra<-paste(primer_rv@extra, shift_base, sep="")
        primer_rv@binding_sequence<-str_sub(primer_rv@binding_sequence, 2)
        shift_base<-str_sub(primer_rv@extra, 1, 1)
        primer_rv@extra<-str_sub(primer_rv@extra, 2)
        primer_rv@overhang<-paste(primer_rv@overhang,shift_base, sep="")
        primer_rv@overhang<-str_sub(primer_rv@overhang, 2)
        primer_rv@temperature<-calculate_tm_nnb(primer_rv@binding_sequence, offset = 0)
        primer_rv@difference<-abs(primer_rv@temperature - primer_unlist[[primer_rv_num -1 ]]@temperature)
      }
    }
    #if(class(primer_fd)=="Primer_MSD") {
      primer_fd@overhang<-paste(comp(shift_base, forceToLower = F), primer_fd@overhang, sep="")
      primer_fd@extra<-paste(str_sub(primer_fd@overhang, 5), primer_fd@extra ,sep="")
      if(class(primer_fd)=="Primer_SPM") {
        primer_fd@binding_sequence<-paste(str_sub(primer_fd@extra, -1), primer_fd@binding_sequence ,sep="")
        primer_fd@extra<-str_sub(primer_fd@extra, 1, -2)
      }
      primer_fd@overhang<-str_sub(primer_fd@overhang, 1, 4)
      primer_fd@temperature<-calculate_tm_nnb(primer_fd@binding_sequence, offset = 0)
      primer_fd@difference<-abs(target_temp - primer_fd@temperature)
      primers[[ceiling(primer_fd_num/2)]][[1]]<-primer_fd
      primers[[ceiling(primer_rv_num/2)]][[2]]<-primer_rv
      break
    #}
    #else{
    #  primer_fd@overhang<-paste(comp(shift_base, forceToLower = F), primer_fd@overhang, sep="")
    #  primer_fd@binding_sequence<-paste(str_sub(primer_fd@overhang, 5), primer_fd@binding_sequence ,sep="")
    #  primer_fd@overhang<-str_sub(primer_fd@overhang, 1, 4)
    #  primer_fd@temperature<-calculate_tm_nnb(primer_fd@binding_sequence)
    #  primer_fd@difference<-abs(target_temp - primer_fd@temperature)
    #  primer_rv@difference<-abs(primer_fd@temperature - primer_rv@temperature)
    #  primers[[ceiling(primer_fd_num/2)]][[1]]<-primer_fd
    #  primers[[ceiling(primer_rv_num/2)]][[2]]<-primer_rv
    #  break
    #}
  }
  
  #overhangs<-sapply(primers, function(x){return(c(x[[1]]@overhang, x[[2]]@overhang))})
  #duplicates<-table(overhangs)
  #duplicates<-duplicates[names(duplicates)!="" & duplicates > 1]
  #if(length(duplicates)==0) {
  #  return(primers)
  #}
  #else{
  return(check_primer_overhangs(primers = primers, fragments = fragments, binding_min_length = binding_min_length, target_temp = target_temp))
  #}
} 