#' List the available Codon Usage Tables
#'
#' @return A list of files
#' @export
#'
#' @examples
#' list_cu_table()
list_cu_table<-function(){
  return(list.files(system.file("cuf", package="GoldenMutagenesis")))
}

#' A function to generate a slim ouput text for primers
#' 
#' An example is shown in the vignette at \url{https://github.com/ipb-halle/GoldenMutagenesis/blob/master/vignettes/Point_Mutagenesis.md}
#' 
#' @param primer An object of class Primer, Primer_MSD or Primerset
#' @return Textual output
#' @export
#' @docType methods
#' @rdname print_primer-methods
#' @examples
#' #Load results of the Point Mutation vignette and print it
#' data(Point_Mutagenesis_BbsI_result)
#' print_primer(primers)
setGeneric("print_primer" , function(primer) {
  standardGeneric("print_primer")
})
#' @rdname print_primer-methods
#' @aliases print_primer,Primer-method 
setMethod("print_primer", signature(primer="Primer"),
          function(primer){
            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@binding_sequence, "\n" , sep="")
            cat("Temperature of binding site: ", primer@temperature, " \u00b0C" , "\n")
            cat("Temperature difference: ", primer@difference, " K", "\n")
          }
)
#' @rdname print_primer-methods
setMethod("print_primer", signature(primer="Primer_MSD"),
          function(primer){
            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@NDT, primer@binding_sequence, "\n", sep="")
            cat("Temperature of binding site: ", primer@temperature, " \u00b0C" , "\n")
            cat("Temperature difference: ", primer@difference, " K", "\n")
          }
)
#' @rdname print_primer-methods
#' @aliases print_primer,Primerset-method 
setMethod("print_primer", signature(primer="Primerset"),
          function(primer){
            primerset<-primer
            for(i in 1:length(primerset@primers)){
              cat("Fragment ", i, "\n", "Forward\n", sep="")
              print_primer(primerset@primers[[i]][[1]])
              cat("Reverse\n")
              print_primer(primerset@primers[[i]][[2]])
              cat("\n")
            }
            cat("Input Sequence:\n", primerset@oldsequence,"\n" )
            cat("\nModified Sequence:\n", primerset@newsequence, "\n")
          }
)
#' @rdname print_primer-methods
#' @aliases print_primer,Extended_Primerset-method 
setMethod("print_primer", signature(primer="Extended_Primerset"),
          function(primer){
            primerset<-primer
            for(i in 1:length(primerset@fragments)){
              cat("Fragment ", i, "\n", sep="")
              cat("Start ", primerset@fragments[[i]]@start, ", ", sep="")
              cat("Stop ",  primerset@fragments[[i]]@stop, ", ", sep="")
              cat("Length ",(primerset@fragments[[i]]@stop - primerset@fragments[[i]]@start)+1, "\n", sep="")
              cat("Forward\n")
              print_primer(primerset@primers[[i]][[1]])
              cat("Reverse\n")
              print_primer(primerset@primers[[i]][[2]])
              cat("\n")
            }
            cat("Input Sequence:\n", primerset@oldsequence,"\n" )
            cat("\nModified Sequence:\n", primerset@newsequence, "\n")
          }
)

#' Domestication of the input sequence
#' 
#' The domesticate function checks for internal cleavage sites. If corresponding sites are present silent mutations are introduced to destroy the recognition sites. 
#' The functions returns a list containing the position of the choosen amino acid residue for silent mutation.
#'
#' @param input_sequence The sequence which should be modified. This is an object of type character containing the sequence. 
#' @param restriction_enzyme Recognition site sequence of the respective restriction enzyme [default: GGTCTC]
#' @param cuf The Codon Usage Table which is being used to select the codon for an exchanged amino acid (and in this case to select the codon which shoulb be replaced). [default: e_coli_316407.csv]
#'
#' @return A list with replacments: Each element has a vector with the codon number at the first slot and the amino acid of this position at the second slot.
#' @export
#' @import seqinr
#' @importFrom seqinr translate
#' @examples
#' #Load the setup of the Point Mutation vignette and run the domestication
#' data(Point_Mutagenesis_BbsI_setup)
#' domesticate(input_sequence, restriction_enzyme=recognition_site_bbsi, cuf=cuf)
domesticate<-function(input_sequence, restriction_enzyme="GGTCTC", cuf="e_coli_316407.csv"){
  cuf_vector<-get_cu_table(cuf, list=F)
  sequence<-s2c(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  prot_sequence<-seqinr::translate(sequence)
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
    replacements[[i]]<-c(as.numeric((start-1)+max_in_list),seqinr::translate(s2c(str_to_upper(names(alt_codons[[max_in_list]][which.max(alt_codons[[max_in_list]])])))))
  }
  return(replacements)
}

#' Calculate primers for a Point Mutagenesis
#' 
#' The mutate function designs the necessary set of primers for the desired mutations.
#' An example is given in the vignette at \url{https://github.com/ipb-halle/GoldenMutagenesis/blob/master/vignettes/Point_Mutagenesis.md}
#' 
#' @param input_sequence The sequence which should be modified. This is an object of type character containing the sequence. 
#' @param prefix Additional nucleobases in 5' position of the recognition site [default: TT]
#' @param restriction_enzyme Recognition site sequence of the respective restriction enzyme [default: GGTCTC]
#' @param suffix Spacer nucleotides matching the cleavage pattern of the enzyme [default: A]
#' @param vector Four basepair overhangs complementary to the created overhangs in the acceptor vector  [default: c("AATG", "AAGC")]
#' @param replacements The desired substitutions
#' @param binding_min_length The minimal threshold value of the template binding sequence [default: 4]
#' @param primer_length Maximal length of the binding sequence [default: 9]
#' @param target_temp Melting temperature of the binding sequence in \code{print('\u00B0')}C [default: 60]
#' @param cuf The Codon Usage Table which is being used to select the codon for an exchanged amino acid. [default: e_coli_316407.csv]
#'
#' @return An object of class Primerset with the designed Primers.
#' @export
#'
#' @examples
#' #Load the setup of the Point Mutation vignette and design the primers
#' data(Point_Mutagenesis_BbsI_setup)
#' primers<-mutate(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, 
#' suffix = "AA", vector=c("CTCA", "CTCG"), replacements = mutations, binding_min_length=4 ,
#' primer_length=9, target_temp=60, cuf=cuf)
#' 
mutate<-function(input_sequence, prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements, binding_min_length=4 ,primer_length=9, target_temp=60, cuf="e_coli_316407.csv") {#change to primer_length_max? and min?
  cuf_list<-get_cu_table(cuf)
  replacements<-order_replacements(replacements)
  sequence<-s2c(input_sequence)
  codon_seq<-sequence_check(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  prot_sequence<-seqinr::translate(sequence)
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


#' Calculate primers for Multiple Site Saturation Mutagenesis
#' 
#' The mutate_msd function designs the necessary set of primers for the desired mutations.
#' 
#' @param input_sequence The sequence which should be modified. This is an object of type character containing the sequence. 
#' @param codon The desired type of MSD mutation [default: NDT]
#' @param prefix Additional nucleobases in 5' position of the recognition site [default: TT]
#' @param restriction_enzyme Recognition site sequence of the respective restriction enzyme [default: GGTCTC]
#' @param suffix Spacer nucleotides matching the cleavage pattern of the enzyme [default: A]
#' @param vector Four basepair overhangs complementary to the created overhangs in the acceptor vector  [default: c("AATG", "AAGC")]
#' @param replacements The desired substitutions as a vector with positions OR a list containing vetors with position (char) and type of MSD mutation (char)
#' @param replacement_range The minimal threshold value of the template binding sequence in amino acid residues [default: 4]
#' @param binding_min_length Maximal length of the binding sequence [default: 9]
#' @param primer_length Melting temperature of the binding sequence in \code{print('\u00B0')}C [default: 60]
#' @param target_temp Maximum distance between two randomization sites to be incoporated into a single primer in amino acid residues [default: 5]
#' @param fragment_min_size Minimal size of a generated gene fragment in base pairs [default 100]
#'
#' @return An object of class Primerset with the designed Primers.
#' @export
#'
#' @examples
#' #Load the setup of the MSD vignette and design the primers
#' data(MSD_BsaI_setup_lv2)
#' print(mutations)
#' print(recognition_site_bsai)
#' primers<-msd_mutate(input_sequence, prefix="TT" ,
#' restriction_enzyme=recognition_site_bsai, suffix="A", 
#' vector=c("AATG", "AAGC"), replacements=mutations, replacement_range=5,
#' binding_min_length=4 , primer_length=9, target_temp=60,
#' fragment_min_size=60 )
msd_mutate<-function(input_sequence, codon="NDT" ,prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements, replacement_range=5, binding_min_length=4 ,primer_length=9, target_temp=60, fragment_min_size=60 ) {#change to primer_length_max? and min?
  codon<-str_to_upper(codon)
  possible_codons<-c("NNN", "NNK", "NNS", "NDT", "DBK", "NRT")
  if(!(codon %in% possible_codons)) {
    stop(paste(codon, "is not a valid codon. Please select one of the following:", paste(possible_codons, collapse = " ") ,sep=" "))
  }
  if(class(replacements)=="list"){
    replacements<-order_replacements(replacements)
    codons<-sapply(replacements, function(x){str_to_upper(as.character(x[2]))})
    if(all(is.element(codons, possible_codons))==F){
      stop(paste(codons, "contains invalid codons. Please select one of the following:", paste(possible_codons, collapse = " ") ,sep=" "))
    }
    replacements<-sapply(replacements, function(x){(as.numeric(x[1]))})
  }
  else{
    replacements<-sort(replacements)
    codons<-rep(codon, length(replacements))
  }
  sequence<-s2c(input_sequence)
  codon_seq<-sequence_check(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  min_fragment<-3*primer_length
  prot_sequence<-seqinr::translate(sequence)
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
          if(length(temp_fragment@stop)!=0) {
            temp_new_fragment <- fragment(start = temp_fragment@stop + 1)
            #temp_old_fragment <- temp_fragment
            temp_fragment <- temp_new_fragment
            rm(temp_new_fragment)
          }
          else {
            if(length(temp_fragment@start)==0) {
              stop("Internal error. Please give a bug report!")
            }
          }
          if(replacements[i]-temp_fragment@start < primer_length+2 ) {
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
      codon_seq[cur_fragment@start_mutation]<-codons[1:length(cur_fragment@start_mutation)]
      codons<-codons[-(1:length(cur_fragment@start_mutation))]
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
      codon_seq[cur_fragment@stop_mutation]<-codons[1:length(cur_fragment@stop_mutation)]
      codons<-codons[-(1:length(cur_fragment@stop_mutation))]
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
  primers<-check_primer_dupplicates(primers, fragments, binding_min_length, target_temp)
  return(eps(oldsequence=input_sequence, primers=primers, newsequence=paste(codon_seq, collapse = ""), fragments=fragments))
}

#' Add a level to exisiting Primerset
#' 
#' This function replaces the prefix, the suffix and the restriction enzyme of a given Primerset to change the design to another Level.
#' You can use this function to convert an exisiting Level 2 Primerset to a Level 0 Primerset for example.
#' Also the overhangs of the first and the last primer will be modified to match the plasmid of the new level.
#'
#' @param primerset An exisiting Primerset (in Level 2)
#' @param prefix Additional nucleobases in 5' position of the new recognition site [default: TT]
#' @param restriction_enzyme Recognition site sequence of the new restriction enzyme (Level 0) [default: GAAGAC]
#' @param suffix Spacer nucleotides matching the cleavage pattern of the enzyme (Level 0) [default: AA]
#' @param vector Four basepair overhangs complementary to the created overhangs in the acceptor vector [default: c("CTCA", "CTCG")]
#'
#' @return A Primerset in the new Level (Level 0)
#' @export
#'
#' @examples
#' #Load level 2 results of the MSD vignette
#' data(MSD_BsaI_result_lv2)
#' primer_add_level(primers,  prefix="TT", 
#' restriction_enzyme="GAAGAC", suffix="AA", vector=c("CTCA", "CTCG"))
#' 
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

#' Create a graphical evaluation of sequencing results
#' 
#' This function creates a graphical evalution of the sequencing results to determine the quality of the created library.
#'
#'The functions aligns the obtained sequencing results to the target gene sequence. 
#'It also tries to align the reverse complement of the obtained sequence. 
#'Afterwards it checks for mismatches between the sequences.
#'Mismatches are likely to be sucessfully mutated nucleotides. 
#'Positions regarded as mismatches are displayed as pie charts. 
#'The shown distributions are based on the signal intensities of the four nucleobases at the mismatch positions.
#'You can compare the pie charts with expected pattern of randomization, therefore validating the quality of the created library.
#' @importFrom dplyr slice
#' @importFrom graphics pie
#' @importFrom sangerseqR readsangerseq peakPosMatrix
#' @importFrom Biostrings pairwiseAlignment mismatchTable reverseComplement
#' @import RColorBrewer
#' @param input_sequence The sequence which was modified. This is an object of type character containing the sequence. 
#' @param ab1file The path to the ab1file which was provided by the sequencer/sequencing service
#' @param replacements The mutations which were desired.
#' @param trace_cutoff The minimal sum of signals (4 nucleotides) for a position in the sequence. [default: 80]
#'
#' @return Plots on the active/default graphics device.
#' @export
#'
#' @examples
#' \dontrun{
#' data(MSD_BsaI_setup_lv2)
#' abfile<-"activesite_for_200718.ab1"
#' base_distribution(input_sequence=input_sequence, ab1file=abfile, replacements=mutations)
#' }
base_distribution<-function(input_sequence, ab1file, replacements, trace_cutoff=80){
  sanger_seq<-sangerseqR::readsangerseq(ab1file) #reading in the data
  global_Align<-Biostrings::pairwiseAlignment(input_sequence, sanger_seq@primarySeq)
  global_Align_rev<-Biostrings::pairwiseAlignment(input_sequence, Biostrings::reverseComplement(sanger_seq@primarySeq))
  reverse=F
  if(global_Align_rev@score > global_Align@score) {
    reverse=T
    global_Align<-global_Align_rev
    print("Reverse sequence detected!")
  }
  mismatches<-Biostrings::mismatchTable(global_Align)
  replacements_basepairs<-as.vector(sapply(replacements, FUN<-function(x){return(c(x*3-2, x*3-1, x*3))}))
  candidates<-unlist(sapply(replacements_basepairs, FUN = function(x){which(mismatches[,"PatternStart"]==x)}, simplify = array))
  mismatches_candidates<-mismatches[candidates, ]
  mismatches_candidates$pos<-mismatches_candidates[,"PatternStart"]%%3
  mismatches_candidates[mismatches_candidates["pos"]==0, "pos"]<-3
  subject_pos<-vector()
  pattern_pos<-vector()
  for (i in 1:nrow(mismatches_candidates)) {
    subject_start<-mismatches_candidates[i, "SubjectStart"]
    pos<-mismatches_candidates[i, "pos"]
    pattern_start<-mismatches_candidates[i, "PatternStart"]
    if(pos==1) {
      subject_pos<-c(subject_pos, subject_start, subject_start+1, subject_start+2)
      pattern_pos<-c(pattern_pos, pattern_start, pattern_start+1, pattern_start+2)
      
    }
    if(pos==2) {
      subject_pos<-c(subject_pos, subject_start-1, subject_start, subject_start+1)
      pattern_pos<-c(pattern_pos, pattern_start-1, pattern_start, pattern_start+1)
      
    }
    if(pos==3) {
      subject_pos<-c(subject_pos, subject_start-2, subject_start-1, subject_start)
      pattern_pos<-c(pattern_pos, pattern_start-2, pattern_start-1, pattern_start)
      
    }
  }
  subject_pos<-unique(subject_pos)
  pattern_pos<-unique(pattern_pos)
  
  if(reverse==T) {
    subject_pos<-length(sanger_seq@primarySeq)-subject_pos+1
  }
  tracematrix_subject<-sangerseqR::traceMatrix(sanger_seq)[sangerseqR::peakPosMatrix(sanger_seq)[subject_pos],]
  sums_row<-which(rowSums(tracematrix_subject)>=trace_cutoff)
  tracematrix_subject<-as.data.frame(tracematrix_subject[sums_row,])
  for(element in sums_row) {
    # plotting as pie chart
    sliceit <- dplyr::slice (tracematrix_subject,element)
    slices <- as.numeric(sliceit)
    lbls <- c("Adenine", "Cytosine", "Guanine", "Thymine")
    if(reverse==T) {
      lbls <- c("Thymine", "Guanine", "Cytosine", "Adenine")
    }
    pct <- round(slices/sum(slices)*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    pie(slices,labels = lbls, col=brewer.pal(4,"Spectral"),main = paste("Peak intensity distribution for \nPosition", pattern_pos[element], "(Template) -", subject_pos[element], "(Sequencing)", sep=" ")) 
    }
}
