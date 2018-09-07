#Functions for GoldenGateProject
library(seqinr)
library(stringr)

list_cu_table<-function(){
  return(list.files(system.file("cuf", package="GoldenMutagenesis")))
}

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

#cuf<-read.csv("")
#cuf[, "codon"]<-str_replace_all(cuf[,"codon"], "U", "T")
#cuf_vector<-as.vector(cuf$relative_frequency)
#names(cuf_vector)<-gsub("U","T",cuf$codon)
#cuf_list<-lapply(unique(cuf$amino_acid), function(x){a<-cuf[which(cuf$amino_acid==x),"relative_frequency"];a<-a*1000;names(a)<-str_to_lower(cuf[which(cuf$amino_acid==x),"codon"]);return(as.table(a))})
#names(cuf_list)<-unique(cuf$amino_acid)


ps<-setClass("Primerset", slots=c(oldsequence="character", primers="list", newsequence="character"))
pc<-setClass("Primer", slots=c(prefix="character" ,restriction_enzyme="character", suffix="character", vector="character", overhang="character" ,binding_sequence="character", temperature="numeric", difference="numeric"))
pc_msd<-setClass("Primer MSD", contains="Primer", slots=c(NDT="character"))
fragment<-setClass("Fragment", slots=c(start="numeric", stop="numeric", start_mutation="vector", stop_mutation="vector"))

order_replacements<-function(replacements){return(replacements[order(sapply(replacements, function(x){as.numeric(x[1])}))])}


remove_newline<-function(x){
  gsub("\r?\n|\r", "", x)
}

calculate_tm<-function(x, salt=50, primer=50, offset=9){
  oligo_sequence<-s2c(x)
  oligo_sequence<-oligo_sequence[offset:length(oligo_sequence)]
  #  Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
  counts<-count(s2c(x), wordsize=1, by=1, alphabet = c("A", "C", "G", "T"))
  tm<-100.5 + (41 * as.numeric(counts["G"] + counts["C"])/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) - (820/as.numeric(counts["A"]+counts["T"]+counts["G"]+counts["C"])) + 16.6*log10(salt/1000)
  return(as.numeric(tm))
}

calculate_UPAC<-function(x, func=calculate_tm, selection="max", temp=0) {
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

setGeneric("print_primer" , function(primer) {
  standardGeneric("print_primer")
})

setMethod("print_primer", signature(primer="Primer"),
          function(primer){
            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@binding_sequence, "\n" , sep="")
            cat("Temperature of binding site: ", primer@temperature, " °C" , "\n")
            cat("Temperature difference: ", primer@difference, " K", "\n")
          }
)
setMethod("print_primer", signature(primer="Primer MSD"),
          function(primer){
            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@NDT, primer@binding_sequence, "\n", sep="")
            cat("Temperature of binding site: ", primer@temperature, " °C" , "\n")
            cat("Temperature difference: ", primer@difference, " K", "\n")
          }
)
setMethod("print_primer", signature(primer="Primerset"),
          function(primer){
            for(i in 1:length(primer@primers)){
              cat("Fragment ", i, "\n", "Forward\n", sep="")
              print_primer(primer@primers[[i]][[1]])
              cat("Reverse\n")
              print_primer(primer@primers[[i]][[2]])
              cat("\n")
            }
            cat("Input Sequence:\n", primer@oldsequence,"\n" )
            cat("\nModified Sequence:\n", primer@newsequence, "\n")
          }
)

sequence_check<-function(input_sequence){
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


mutate<-function(input_sequence, prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements, binding_min_length=4 ,primer_length=9, target_temp=60, cuf_list=get_cu_table("e_coli_316407.csv")) {#change to primer_length_max? and min?
  replacements<-order_replacements(replacements)
  sequence<-s2c(input_sequence)
  codon_seq<-sequence_check(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  prot_sequence<-translate(sequence)
  primers<-vector("list", length(replacements)+1)
  #First primer @ transcription start
  forward<-pc(prefix=prefix, restriction_enzyme = restriction_enzyme, suffix=suffix, vector=vector[1], overhang="")
  
  #start_primer_non_binding=paste(prefix, restriction_enzyme, suffix, vector[1], sep="")
  if(str_sub(vector[1], 2) == "ATG"){
    binding_start<-2
  } else{
    binding_start<-1
  }
  forward@binding_sequence<-paste(paste(codon_seq[binding_start:(binding_start+primer_length)-1], collapse=""), sep="") #Set to 2 because the sequence starts with ATG
  #Remove also trailing TAA?
  #Maybe check of the sequence before, afterwards removing of those parts?
  #Getting out shorter primer sequence based on melting temperature
  start_primer<-sequence_length_temperature(forward, primer_min=binding_min_length, target_temp=target_temp)
  rm(forward)
  primers[[1]]<-vector("list", 2)
  primers[[1]][[1]]<-start_primer
  for(i in 1:length(replacements)) {
    position_aa<-as.numeric(replacements[[i]][1])
    position<-position_aa*3
    aminoacid<-replacements[[i]][2]
    codon<-str_to_upper(names(which.max(cuf_list[[aminoacid]]))[1])
    if(codon == codon_seq[position_aa]) {
      if(length(cuf_list[[aminoacid]] == 1)) {
        stop(paste("There is no syn. codon for", aminoacid ,sep=""))
      }
      else {
        old_codon<-which.max(cuf_list[[aminoacid]])
        codon<-str_to_upper(names(which.max(cuf_list[[aminoacid]][-old_codon]))[1])
      }
    }
    codon_seq[position_aa]<-codon
    forward<-pc(prefix=prefix, restriction_enzyme = restriction_enzyme, suffix=suffix, vector="", overhang=paste(str_split(codon_seq[position_aa-1], "", simplify = T)[3], codon_seq[position_aa], sep=""))
    forward@binding_sequence<-paste(paste(codon_seq[(position_aa+1):(position_aa+min(primer_length, length(codon_seq)-position_aa))], collapse=""), sep="") 
    primer_forward<-sequence_length_temperature(forward, primer_min=binding_min_length, target_temp=target_temp)
    primers[[i+1]]<-list(primer_forward, NULL)
    reverse<-pc(prefix = prefix, restriction_enzyme = restriction_enzyme, suffix = suffix, vector = "")
    overlap_binding_sequence<-str_to_upper(paste(comp(rev(s2c(paste(codon_seq[max(position_aa-primer_length, 1):position_aa], collapse="")))), collapse =""))
    reverse@binding_sequence<-str_sub(overlap_binding_sequence, 5)
    reverse@overhang<-str_sub(overlap_binding_sequence, 1, 4)
    primer_reverse<-sequence_length_temperature(reverse, primer_min = binding_min_length, target_temp = primers[[i]][[1]]@temperature )
    primers[[i]][[2]]<-primer_reverse
    rm(forward)
    rm(reverse)
  }
  reverse<-pc(prefix=prefix, restriction_enzyme = restriction_enzyme, suffix = suffix, vector = vector[2], overhang="")
  reverse@binding_sequence<-str_to_upper(paste(comp(rev(s2c(paste(codon_seq[(length(codon_seq)-primer_length):length(codon_seq)], collapse="")))), collapse =""))
  end_primer<-sequence_length_temperature(reverse, primer_min=binding_min_length, target_temp=primers[[length(primers)]][[1]]@temperature)
  
  primers[[length(primers)]][[2]]<-end_primer
  return(ps(oldsequence=input_sequence, primers=primers, newsequence=paste(codon_seq, collapse = "")))
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


msd_mutate<-function(input_sequence, prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements, replacement_range=5, binding_min_length=4 ,primer_length=9, target_temp=60, fragment_min_size=60 ) {#change to primer_length_max? and min?
  replacements<-sort(replacements)
  sequence<-s2c(input_sequence)
  codon_seq<-sequence_check(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  min_fragment<-3*primer_length
  prot_sequence<-translate(sequence)
  primers<-vector("list")
  if(str_sub(vector[1], 2) == "ATG"){
    fragment_start<-2
  } else{
    fragment_start<-1
  }
  #next_replacement<-replacements[1]
  #current_replacement<-0
  replacement_distances<-as.matrix(dist(replacements))
  #First calculate fragments
  #then calculate primers 
  fragments<-c()
  i<-1
  #todo replace formular with primer_length
  repeat {
    #################################
    #################################
    #Creation of the first fragment/primer
    if (i == 1) {
      #first replacement
      #check if it is on the beginning of the first fragment
      temp_fragment <- fragment(start = fragment_start)
      if (replacements[i] <  (min_fragment+fragment_start - 1)) {
        #if (replacements[i] <= fragment_start - 1 + replacement_range + 2) {
        temp_fragment@start_mutation <- replacements[i]
        if (length(replacements) == 1) {
          temp_fragment@stop <- length(codon_seq)
          fragments <- c(fragments, temp_fragment)
          break
        }
      }
      else {
        # we will create a new fragment
        #check if there is enough space to create a new fragment
        if (length(replacements) == 1) {
          temp_fragment@stop <-  replacements[i] + 2
          temp_fragment@stop_mutation <- c(replacements[i])
          fragments <- c(fragments, temp_fragment)
          temp_fragment <-
            fragment(start = fragments[[1]]@start + 1,
                     stop = length(codon_seq))
          fragments <- c(fragments, temp_fragment)
          break
        }
        if (replacement_distances[i + 1, i] < 3) {
          temp_fragment@stop <- replacements[i] - 1
        }
        else{
          temp_fragment@stop <-  replacements[i] + 2
          temp_fragment@stop_mutation <- c(replacements[i])
          fragments <- c(fragments, temp_fragment)
          i <- i + 1
        }
      }
    }
    #################################
    #################################
    #generic part for all fragments
    mutations_in_fragment_range <-
      as.numeric(which(
        replacement_distances[, i] < (min_fragment) &
          replacement_distances[, i] > 0
      ))
    mutations_in_fragment_range <-
      mutations_in_fragment_range[which(mutations_in_fragment_range > i)]
    ###distinguish between fragment start and end
    ###if the minimal binding length is very high and the distance of the mutations is very short, you will get very long primers!
    if (length(temp_fragment@stop) == 0) {
      if (any(temp_fragment@start_mutation[length(temp_fragment@start_mutation)] == replacements[i])) {
        if (length(mutations_in_fragment_range) != 0) {
          #we won't create a new fragment until we have enough space between the mutations to create a forward and a reverse primer
          for (mutation in mutations_in_fragment_range[1:(length(mutations_in_fragment_range) - 1)]) {
            temp_fragment@start_mutation <-
              c(temp_fragment@start_mutation, replacements[mutation])
            i <-
              mutations_in_fragment_range[length(mutations_in_fragment_range)]
          }
          rm(mutation)
        }
        else {
          i <- i + 1
        }
      }
      else {
        #We can finish the fragment here
        if (replacement_distances[i + 1, i] < 3) {
          #Check if we have the next mutation on the same primer
          temp_fragment@stop = replacements[i] - 1
        }
        else {
          #End the fragment with the mutation
          temp_fragment@stop <- replacements[i] + 2
          temp_fragment@stop_mutation <-
            c(temp_fragment@stop_mutation, replacements[i])
          i <- i + 1
        }
        fragments <- c(fragments, temp_fragment)
      }
    }  
    else {
      ###create a new fragment in forward direction
      #we will always start with a mutation in forward direction (execpt the first one)
      temp_new_fragment <- fragment(start = temp_fragment@stop + 1)
      temp_old_fragment <- temp_fragment
      temp_fragment <- temp_new_fragment
      rm(temp_new_fragment)
      #does the fragment need additional shifiting?
      if (replacements[i] - temp_fragment@start < (min_fragment)) {
        #Mutation(s) at the beginning of the fragment
        #Check if there are any mutations after the current mutation which are too far away
        if (length(mutations_in_fragment_range[mutations_in_fragment_range > replacement_range]) > 0) {
          #Then put the current mutation to the end of the last fragment
          fragments[[length(fragments)]]@stop <- replacements[i] + 1
          fragments[[length(fragments)]]@stop_mutation <-
            c(fragments[[length(fragments)]]@stop_mutation, replacements[i])
          temp_fragment <- fragments[[length(fragments)]]
          rm(temp_new_fragment)
          rm(temp_old_fragment)
          i <- i + 1
        }
        else {
          if (length(temp_old_fragment@stop_mutation) == 0) {
            #there is no mutation at the end, we can shift anyway
            fragments[[length(fragments)]]@stop <-
              replacements[i] - 1
          }
          else {
            #Check the difference between the last and the current mutation
            last_mutation <-
              fragments[[length(fragments)]]@stop_mutation[length(fragments[[length(fragments)]]@stop_mutation)]
            diff <- replacements[i] - last_mutation - 1
            if (diff <= replacement_range) {
              fragments[[length(fragments)]]@stop <- replacements[i] - 1
            }
            else{
              fragments[[length(fragments)]]@stop <-
                fragments[[length(fragments)]]@stop + replacement_range
            }
          }
          temp_fragment@start <-
            fragments[[length(fragments)]]@stop + 1
          temp_fragment@start_mutation <-
            c(temp_fragment@start_mutation, replacements[i])
          
        }
      }
    }
    #################################
    #################################
    #part for the last mutation and the last fragment
    if (i == length(replacements)) {
      if (any(temp_fragment@start_mutation[length(temp_fragment@start_mutation)] == replacements[i])) {
        #we started a fragment with the mutation on the forward part
        #just end it here
        temp_fragment@stop <- length(codon_seq)
        fragments <- c(fragments, temp_fragment)
        rm(temp_fragment)
        break
      }
      else{
        if ((length(codon_seq) - replacements[i-1]) >= (min_fragment)) {
          temp_new_fragment <- fragment(start = temp_fragment@stop + 1)
          temp_old_fragment <- temp_fragment
          temp_fragment <- temp_new_fragment
          rm(temp_new_fragment)
          if(replacements[i]-temp_fragment@start < primer_length+2) {
            temp_fragment@start_mutation<-replacements[i]
            temp_fragment@stop<-length(codon_seq)
            fragments <- c(fragments, temp_fragment)
            break
          }
          else {
            temp_fragment@stop_mutation<-replacements[i]
            if ((length(codon_seq) - replacements[i]) > (min_fragment)) {
              temp_fragment@stop<-temp_fragment@stop_mutation+2
              fragments <- c(fragments, temp_fragment)
              temp_fragment<-fragment(start = temp_fragment@stop + 1, stop = length(codon_seq))
              fragments <- c(fragments, temp_fragment)
              break
            }
            else{
              temp_fragment@stop<-length(codon_seq)
              fragments <- c(fragments, temp_fragment)
              break
            }
          }
        }
        else {
          temp_fragment@stop<-length(codon_seq)
          temp_fragment@stop_mutation <-
            c(temp_fragment@stop_mutation, replacements[i])
          fragments[length(fragments)] <- temp_fragment
          break
        }
      }
    }
  }
  
  primers<-vector("list", length(fragments))
  for(n in 1:length(fragments)){
    cur_fragment<-fragments[[n]]
    primers[[n]]<-vector("list", 2)
    if(n==1) { #the first primer does not need any existing overlap
      vector_f=vector[1]
      vector_r<-""
      overhang_f<-""
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[cur_fragment@stop-1], start=3), codon_seq[cur_fragment@stop], sep=""))))), collapse="")
      suffix_f<-"A"
      suffix_r<-comp("A", forceToLower = F)
      stop_r<-cur_fragment@stop-2
    } 
    else if(n==length(fragments)) {
      vector_f<-""
      vector_r<-vector[2]
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-""
      suffix_f<-"A"
      suffix_r<-comp("A", forceToLower = F)
      stop_r<-cur_fragment@stop
    }
    else {
      vector_f=""
      vector_r=""
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[fragments[[n]]@stop-1], start=3), codon_seq[fragments[[n]]@stop], sep=""))))), collapse="")
      suffix_f="A"
      suffix_r=comp("A", forceToLower = F)
      stop_r<-cur_fragment@stop-2
    }
    #forward
    if(length(cur_fragment@start_mutation)==0) {
      temp_primer<-pc(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      temp_primer@binding_sequence<-paste(paste(codon_seq[cur_fragment@start:(cur_fragment@start+primer_length-1)], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    else {
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      codon_seq[cur_fragment@start_mutation]<-"NDT"
      temp_primer@NDT<-paste(paste(codon_seq[cur_fragment@start:max(cur_fragment@start_mutation)], collapse = ""), sep="")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(max(cur_fragment@start_mutation)+1):((max(cur_fragment@start_mutation)+1)+primer_length-1)], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    primers[[n]][[1]]<-temp_primer
    rm(temp_primer)
    #reverse
    if(length(cur_fragment@stop_mutation)==0){
      temp_primer<-pc(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      temp_primer@binding_sequence<-paste(paste(codon_seq[(cur_fragment@stop-2-primer_length-1):stop_r], collapse=""), sep="")
      if(n!=length(fragments))
        temp_primer@binding_sequence<-paste(temp_primer@binding_sequence, str_sub(codon_seq[cur_fragment@stop-1], end=2) ,sep="")
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    else{
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      codon_seq[cur_fragment@stop_mutation]<-"NDT"
      temp_primer@NDT<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)):stop_r], collapse=""), sep="")
      temp_primer@NDT<-paste(temp_primer@NDT, str_sub(codon_seq[cur_fragment@stop-1], end=2), sep="")
      temp_primer@NDT<-paste(comp(rev(s2c(temp_primer@NDT)), ambiguous = T,forceToLower = F), collapse = "")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)-1-primer_length-1):(min(cur_fragment@stop_mutation)-1)], collapse=""), sep="")
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    primers[[n]][[2]]<-temp_primer
    rm(temp_primer)
  }
  #check the primers with the checkprimer function:
  #Check for primers with same overlap
  #Replace them based on ? -> Primer without mutation/length of the primer in total?
  #It is easier to modify the exisiting primer 
  #If it is not possible to correct all overlaps -> return message with postion for silent mutation
  primers<-check_primers(primers, fragments, binding_min_length, target_temp)
  return(ps(oldsequence=input_sequence, primers=primers, newsequence=paste(codon_seq, collapse = "")))
}

primer_add_level<-function(primerset, prefix="TT" ,restriction_enzyme="GAAGAC", suffix="AA", vector=c("CTCA", "CTCG")){
  for(i in 1:length(primerset@primers)) {
    if(primerset@primers[[i]][[1]]@overhang=="" && primerset@primers[[i]][[1]]@vector!=""){
      primerset@primers[[i]][[1]]@overhang<-primerset@primers[[i]][[1]]@vector
    }
    primerset@primers[[i]][[1]]@vector<-vector[1]
    primerset@primers[[i]][[1]]@prefix<-prefix
    primerset@primers[[i]][[1]]@restriction_enzyme<-restriction_enzyme
    primerset@primers[[i]][[1]]@suffix<-suffix
    ############################
    
    if(primerset@primers[[i]][[2]]@overhang=="" && primerset@primers[[i]][[2]]@vector!=""){
      primerset@primers[[i]][[2]]@overhang<-primerset@primers[[i]][[2]]@vector
    }
    primerset@primers[[i]][[2]]@vector<-vector[2]
    primerset@primers[[i]][[2]]@prefix<-prefix
    primerset@primers[[i]][[2]]@restriction_enzyme<-restriction_enzyme
    primerset@primers[[i]][[2]]@suffix<-suffix
  }
  return(primerset)
}

domesticate<-function(input_sequence, restriction_enzyme="GGTCTC", cuf_vector=get_cu_table("e_coli_316407.csv", list = F)){
  sequence<-s2c(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  prot_sequence<-translate(sequence)
  matches <- do.call(rbind, str_locate_all(input_sequence, c(restriction_enzyme, restriction_enzyme_reverse))) # Returns positions of every match in a string
  if(nrow(matches) == 0) {
    print("No domestication needed.")
    return(list())
  }
  split_seq<-splitseq(sequence)
  replacements<-vector(mode = "list", length = nrow(matches))
  for(i in 1:nrow(matches)){
    start<-ceiling(matches[i,"start"]/3)
    end<-ceiling(matches[i, "end"]/3)
    protein_domest<-split_seq[start:end]
    alt_codons<-syncodons(protein_domest) #we could use synsequence for this task
    for(j in 1:length(alt_codons)){ #filter out the codon which we already know
      codons<-alt_codons[[j]]
      codons<-codons[which(codons != names(alt_codons[j]))]
      alt_codons[[j]]<-cuf_vector[str_to_upper(codons)]
    }
    max_in_list<-which.max(unlist(lapply(alt_codons, function(x) x[which.max(x)])))  
    replacements[[i]]<-c(as.numeric((start-1)+max_in_list),translate(s2c(str_to_upper(names(alt_codons[[max_in_list]][which.max(alt_codons[[max_in_list]])])))))
  }
  return(replacements)
}
