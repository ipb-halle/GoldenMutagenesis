#Functions for GoldenGateProject
#' @import seqinr
#' @import stringr
#' @import methods
#' @importFrom stats dist
#' @importFrom utils read.csv
NULL
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



calculate_tm<-function(x, salt=50, primer=50, offset=9){
  oligo_sequence<-s2c(x)
  oligo_sequence<-oligo_sequence[offset:length(oligo_sequence)]
  #  Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
  counts<-count(s2c(x), wordsize=1, by=1, alphabet = c("A", "C", "G", "T"))
  tm<-100.5 + (41 * as.numeric(counts["G"] + counts["C"])/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) - (820/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) + 16.6*log10(salt/1000)
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

calculate_tm_nnb<-function(oligo_sequence, primer_concentration=50, salt_concentration=50, offset=9){
  oligo_sequence_s2c<-s2c(oligo_sequence)
  oligo_sequence<-paste(oligo_sequence_s2c[offset:length(oligo_sequence_s2c)], collapse="")
  K<-1/(primer_concentration*1e-9) #Convert from nanomoles to moles
  R<-1.987
  RlnK<-R*log(K)
  result<-(1000*(calculate_DeltaH(oligo_sequence)-3.4)/(calculate_DeltaS(oligo_sequence)+RlnK)-272.9)
  result<-result+7.21*log(salt_concentration/1000)
  return(result)
}


setGeneric("sequence_length_temperature" , function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60) {
  standardGeneric("sequence_length_temperature")
})

setMethod("sequence_length_temperature", signature(primer="Primer"),
          function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60){
            primer_seq_s2c<-s2c(primer@binding_sequence)
            temperatures<-list()
            names_i<-vector()
            for(i in (primer_min*3):length(primer_seq_s2c)){
              temperatures<-c(temperatures, temp_func(paste(primer_seq_s2c[1:i],collapse=""), offset=0))
              names_i<-c(names_i, i)
            }
            names(temperatures)<-names_i
            diff<-unlist(lapply(temperatures, function(x){abs(x-target_temp)}))
            candidate<-as.numeric(names(diff[diff==min(diff)]))
            primer@binding_sequence<-paste(primer_seq_s2c[1:min(candidate)],collapse="")
            primer@temperature<-temperatures[[as.character(min(candidate))]]
            primer@difference<-as.numeric(diff[as.character(min(candidate))])
            return(primer)
          }
)

setMethod("sequence_length_temperature", signature(primer="Primer MSD"),
          function(primer, temp_func=calculate_tm_nnb, primer_min=3, target_temp=60){
            callNextMethod()
          }
)


sequence_check<-function(input_sequence){
  input_sequence<-str_to_upper(input_sequence)
  if(nchar(input_sequence)%%3!=0) {
    stop(paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "))
  }
  codon_seq<-splitseq(s2c(str_to_upper(input_sequence)))
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
    stop("No stop codon in the provided sequence.Stopping here. Please check the provided sequence!")
  }
  
  if(max(stop) != length(codon_seq)) {
    warning(paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "))  
    codon_seq<-codon_seq[1:max(stop)]
  }# else {
  #codon_seq <- codon_seq[-length(codon_seq)]
  #}
  return(codon_seq)
}

check_primers<-function(primers, fragments, binding_min_length=4, target_temp=60) {
  overhangs<-sapply(primers, function(x){return(c(x[[1]]@overhang, x[[2]]@overhang))})
  duplicates<-table(overhangs)
  duplicates<-duplicates[names(duplicates)!="" & duplicates > 1]
  if(length(duplicates)==0) {
    return(primers)
  }
  duplicate<-duplicates[1]
  primer_num<-which(overhangs==names(duplicate))
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
    if(class(primer_rv)=="Primer MSD") {
      if(str_sub(primer_rv@NDT, 1, 3) == "AHN") {
        if(i == length(primer_num)) {
          stop(paste("We can not fix overlaps in the primers. Please consider a silent mutation at position", fragments[[ceiling((primer_rv_num+1)/2)]]@start))
        }
        else {
          next
        }
      }
      else{
        shift_base<-str_sub(primer_rv@NDT, 1, 1)
        primer_rv@overhang<-paste(primer_rv@overhang,shift_base, sep="")
        primer_rv@NDT<-str_sub(primer_rv@NDT, 2)
        primer_rv@overhang<-str_sub(primer_rv@overhang, 2)
      }
    } else {
      if(nchar(primer_rv@binding_sequence) < 3 * binding_min_length) {
        if(i == length(primer_num)) {
          stop(paste("We can not fix overlaps in the primers. Please consider a silent mutation at position", fragments[[ceiling((primer_rv_num+1)/2)]]@start))
        }
        else{
          next
        }
      }
      else{
        shift_base<-str_sub(primer_rv@binding_sequence, 1, 1)
        primer_rv@overhang<-paste(primer_rv@overhang,shift_base, sep="")
        primer_rv@binding_sequence<-str_sub(primer_rv@binding_sequence, 2)
        primer_rv@overhang<-str_sub(primer_rv@overhang, 2)
        primer_rv@temperature<-calculate_tm_nnb(primer_rv@binding_sequence)
        primer_rv@difference<-abs(primer_rv@temperature - primer_unlist[[primer_rv_num -1 ]]@temperature)
      }
    }
    if(class(primer_fd)=="Primer MSD") {
      primer_fd@overhang<-paste(comp(shift_base, forceToLower = F), primer_fd@overhang, sep="")
      primer_fd@NDT<-paste(str_sub(primer_fd@overhang, 5), primer_fd@NDT ,sep="")
      primer_fd@overhang<-str_sub(primer_fd@overhang, 1, 4)
      primers[[ceiling(primer_fd_num/2)]][[1]]<-primer_fd
      primers[[ceiling(primer_rv_num/2)]][[2]]<-primer_rv
      break
    }
    else{
      primer_fd@overhang<-paste(comp(shift_base, forceToLower = F), primer_fd@overhang, sep="")
      primer_fd@binding_sequence<-paste(str_sub(primer_fd@overhang, 5), primer_fd@binding_sequence ,sep="")
      primer_fd@overhang<-str_sub(primer_fd@overhang, 1, 4)
      primer_fv@temperature<-calculate_tm_nnb(primer_fv@binding_sequence)
      primer_fv@difference<-abs(target_temp - primer_fd@temperature)
      primers[[ceiling(primer_fd_num/2)]][[1]]<-primer_fd
      primers[[ceiling(primer_rv_num/2)]][[2]]<-primer_rv
      break
    }
  }
  
  #overhangs<-sapply(primers, function(x){return(c(x[[1]]@overhang, x[[2]]@overhang))})
  #duplicates<-table(overhangs)
  #duplicates<-duplicates[names(duplicates)!="" & duplicates > 1]
  #if(length(duplicates)==0) {
  #  return(primers)
  #}
  #else{
  return(check_primers(primers = primers, fragments = fragments, binding_min_length = binding_min_length, target_temp = target_temp))
  #}
} 

