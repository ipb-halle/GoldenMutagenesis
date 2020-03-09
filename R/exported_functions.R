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

#' A function to generate an output text for primers
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
            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@extra, primer@binding_sequence, "\n", sep="")
            cat("Temperature of binding site: ", primer@temperature, " \u00b0C" , "\n")
            cat("Temperature difference: ", primer@difference, " K", "\n")
          }
)
#' @rdname print_primer-methods
#setMethod("print_primer", signature(primer="Primer_MSD"),
#          function(primer){
#            cat(primer@prefix, primer@restriction_enzyme, primer@suffix, primer@vector, primer@overhang, primer@extra, primer@binding_sequence, "\n", sep="")
#            cat("Temperature of binding site: ", primer@temperature, " \u00b0C" , "\n")
#            cat("Temperature difference: ", primer@difference, " K", "\n")
#          }
#)
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
  restriction_enzyme<-str_to_upper(restriction_enzyme)
  input_sequence<-str_to_upper(input_sequence)
  sequence<-s2c(input_sequence)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  prot_sequence<-translate(sequence)
  matches <- do.call(rbind, str_locate_all(input_sequence, c(restriction_enzyme, restriction_enzyme_reverse))) # Returns positions of every match in a string
  if(nrow(matches) == 0) {
    message("No domestication needed.")
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

#' Calculate Primers for Point Mutagenesis
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
#' @param replacement_range  Maximum distance in amino acid residues between two randomization sites to be incoporated into a single primer (reverse, end of the fragment) - has a cascading effect for following mutations [default: 2]
#' @param binding_min_length The minimal threshold value of the length of the template binding sequence in amino acid residues [default: 4]
#' @param binding_max_length Maximal length of the binding sequence in amino acid residues [default: 9]
#' @param target_temp Melting temperature of the binding sequence in \code{print('\u00B0')}C [default: 60]
#' @param fragment_min_size Minimal size of a generated gene fragment in base pairs [default 100]
#' @param cuf The Codon Usage Table which is being used to select the codon for an exchanged amino acid. [default: e_coli_316407.csv]
#'
#' @return An object of class Primerset with the designed Primers.
#' @export
#'
#' @examples
#' #Load the setup of the Point Mutation vignette and design the primers
#' data(Point_Mutagenesis_BbsI_setup)
#' primers<-mutate_spm(input_sequence, prefix="TT", restriction_enzyme = recognition_site_bbsi, 
#' suffix = "AA", vector=c("CTCA", "CTCG"), replacements = mutations, binding_min_length=4 ,
#' binding_max_length=9, target_temp=60, cuf=cuf)
#' 
mutate_spm<-function(input_sequence, prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements,replacement_range=2, binding_min_length=4, binding_max_length=9, target_temp=60, cuf="e_coli_316407.csv", fragment_min_size=100) {#change to binding_max_length_max? and min?
  cuf_list<-get_cu_table(cuf)
  prefix<-str_to_upper(prefix)
  vector<-str_to_upper(vector)
  suffix<-str_to_upper(suffix)
  replacements<-order_replacements(replacements)
  sequence<-s2c(input_sequence)
  codon_seq<-sequence_check(input_sequence)
  restriction_enzyme<-str_to_upper(restriction_enzyme)
  restriction_enzyme_s2c<-s2c(restriction_enzyme)
  restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
  restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
  restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
  min_fragment<-3*binding_max_length
  prot_sequence<-seqinr::translate(sequence)
  primers<-vector("list")
  
  check_offtargets<-function(codon_seq, codons, positions_aa, restriction_enzyme="GGTCTC", cuf="e_coli_316407.csv") {
    cuf_vector<-get_cu_table(cuf, list=F)
    seq<-paste(codon_seq, collapse = "")
    restriction_enzyme_s2c<-s2c(restriction_enzyme)
    restriction_enzyme_s2c_reverse<-comp(restriction_enzyme_s2c)
    restriction_enzyme_s2c_reverse<-rev(restriction_enzyme_s2c_reverse)
    restriction_enzyme_reverse<-str_to_upper(paste(restriction_enzyme_s2c_reverse, collapse = ""))
    matches <- do.call(rbind, str_locate_all(seq, c(restriction_enzyme, restriction_enzyme_reverse))) # Returns positions of every match in a string
    if(nrow(matches) > 0) {
      for(i in 1:nrow(matches)){
        start<-ceiling(matches[i,"start"]/3)
        end<-ceiling(matches[i, "end"]/3)
        positions<-positions_aa[positions_aa >= start & positions_aa <= end]
        if(length(positions > 0)){
          alt_codons<-syncodons(codon_seq_tmp[positions])
          for(j in 1:length(alt_codons)){ #filter out the codon which we already know
            codons_tmp<-alt_codons[[j]]
            codons_tmp<-codons_tmp[which(codons_tmp != names(alt_codons[j]))]
            alt_codons[[j]]<-cuf_vector[str_to_upper(codons_tmp)]
          }
          max_in_list<-which.max(unlist(lapply(alt_codons, function(x) x[which.max(x)]))) 
          codons[which(positions_aa==positions[max_in_list])]<-str_to_upper(names(alt_codons[[max_in_list]][which.max(alt_codons[[max_in_list]])]))
        }
      }
    }
    return(codons)
  }
  
  if(str_sub(vector[1], 2) == "ATG"){
    fragment_start<-2
  } else{
    fragment_start<-1
  }

  positions<-vector()
  positions_aa<-vector()
  aminoacids<-vector()
  codons<-vector()
  codon_seq_tmp<-codon_seq
  
  for(i in 1:length(replacements)) {
      position_aa<-as.numeric(replacements[[i]][1])
      positions_aa<-c(positions_aa, position_aa)
      position<-position_aa*3
      positions<-c(positions, position)
      aminoacid<-replacements[[i]][2]
      aminoacids<-c(aminoacids, aminoacid)
      codon<-str_to_upper(names(which.max(cuf_list[[aminoacid]]))[1])
      if(codon == codon_seq[position_aa]) {
        if(length(cuf_list[[aminoacid]] == 1)) {
          stop(paste("There is no syn. codon for ", aminoacid ,sep=""))
        }
        else {
          old_codon<-which.max(cuf_list[[aminoacid]])
          codon<-str_to_upper(names(which.max(cuf_list[[aminoacid]][-old_codon]))[1])
        }
      }
      codons<-c(codons, codon)
      codon_seq_tmp[position_aa]<-codon
  }
  
  codons<-check_offtargets(codon_seq = codon_seq_tmp, codons = codons, positions_aa = positions_aa, restriction_enzyme = "GGTCTC", cuf = cuf)
  codons<-check_offtargets(codon_seq = codon_seq_tmp, codons = codons, positions_aa = positions_aa, restriction_enzyme = "GAAGAC", cuf = cuf)
  if(restriction_enzyme != "GGTCTC" | restriction_enzyme != "GAAGAC") {
    codons<-check_offtargets(codon_seq = codon_seq_tmp, codons = codons, positions_aa = positions_aa, restriction_enzyme = restriction_enzyme, cuf = cuf)
  }
  
  fragments<-make_fragments(mutations=positions_aa,seq=codon_seq,start=fragment_start, distance=replacement_range, fsize=fragment_min_size, buffer=0)

  primers<-vector("list", length(fragments))
  for(n in 1:length(fragments)){
    cur_fragment<-fragments[[n]]
    primers[[n]]<-vector("list", 2)
    if(length(cur_fragment@start_mutation) != 0) {
      codon_seq[cur_fragment@start_mutation]<-codons[1:length(cur_fragment@start_mutation)]
      codons<-codons[-(1:length(cur_fragment@start_mutation))]
    }
    if(length(cur_fragment@stop_mutation) != 0) {
      codon_seq[cur_fragment@stop_mutation]<-codons[1:length(cur_fragment@stop_mutation)]
      codons<-codons[-(1:length(cur_fragment@stop_mutation))]
      }
    if(n==1) { #the first primer does not need any existing overlap
      vector_f=vector[1]
      vector_r<-""
      overhang_f<-""
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[cur_fragment@stop-1], start=3), codon_seq[cur_fragment@stop], sep=""))))), collapse="")
      suffix_f<-suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop-2
      if(n==length(fragments)){#if there is just one single fragment - e.g. one mutation at the start or end
        overhang_r<-""
        stop_r<-cur_fragment@stop
        vector_r<-vector[2]
      }
    } 
    else if(n==length(fragments)) {
      vector_f<-""
      vector_r<-vector[2]
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-""
      suffix_f<-suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop
    }
    else {
      vector_f=""
      vector_r=""
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[fragments[[n]]@stop-1], start=3), codon_seq[fragments[[n]]@stop], sep=""))))), collapse="")
      suffix_f=suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop-2
    }
    #forward
    if(length(cur_fragment@start_mutation)==0) {
      temp_primer<-pc_spm(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      temp_primer@binding_sequence<-paste(paste(codon_seq[cur_fragment@start:(cur_fragment@start+binding_max_length-1)], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    else {
      temp_primer<-pc_spm(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      temp_primer@extra<-paste(paste(codon_seq[cur_fragment@start:max(cur_fragment@start_mutation)], collapse = ""), sep="")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(max(cur_fragment@start_mutation)+1):((max(cur_fragment@start_mutation)+1)+binding_max_length-1)], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    primers[[n]][[1]]<-temp_primer
    rm(temp_primer)
    #reverse
    if(length(cur_fragment@stop_mutation)<=1){
      temp_primer<-pc_spm(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      temp_primer@binding_sequence<-paste(paste(codon_seq[(cur_fragment@stop-2-binding_max_length-1):stop_r], collapse=""), sep="")
      if(n!=length(fragments))
        temp_primer@binding_sequence<-paste(temp_primer@binding_sequence, str_sub(codon_seq[cur_fragment@stop-1], end=2) ,sep="")
      #check if the mutation is also the end
      if(length(cur_fragment@stop_mutation)==1){
          if(cur_fragment@stop_mutation[1] != cur_fragment@stop) {
            if(cur_fragment@stop_mutation[1]<=stop_r) # if it is bigger than stop_r, it is partly in the overlap...
              temp_primer@extra<-paste(paste(codon_seq[cur_fragment@stop_mutation[1]:stop_r], collapse = ""), sep="")
            if(n!=length(fragments)) #...and handeled here, for the last fragment stop_r is stop
              temp_primer@extra<-paste(temp_primer@extra, str_sub(codon_seq[cur_fragment@stop-1], end=2) ,sep="")
            temp_primer@extra<-paste(comp(rev(s2c(temp_primer@extra)), ambiguous = T,forceToLower = F), collapse = "")
            temp_primer@binding_sequence<-paste(paste(codon_seq[(cur_fragment@stop_mutation[1]-1-binding_max_length-1):(cur_fragment@stop_mutation[1]-1)], collapse=""), sep="")
          }
      }
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    else{
      temp_primer<-pc_spm(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      temp_primer@extra<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)):stop_r], collapse=""), sep="")
      temp_primer@extra<-paste(temp_primer@extra, str_sub(codon_seq[cur_fragment@stop-1], end=2), sep="")
      temp_primer@extra<-paste(comp(rev(s2c(temp_primer@extra)), ambiguous = T,forceToLower = F), collapse = "")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)-1-binding_max_length-1):(min(cur_fragment@stop_mutation)-1)], collapse=""), sep="")
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    primers[[n]][[2]]<-temp_primer
    rm(temp_primer)
  }
  #Check the primers with the checkprimer function:
  #Check for primers with same overlap
  #Replace them based on ? -> Primer without mutation/length of the primer in total?
  #It is easier to modify the exisiting primer 
  #If it is not possible to correct all overlaps -> return message with position for silent mutation
  primers<-check_primer_overhangs(primers, fragments, binding_min_length, target_temp)
  return(eps(oldsequence=input_sequence, primers=primers, newsequence=paste(codon_seq, collapse = ""), fragments=fragments))
}


#' Calculate Primers for Multiple Site Saturation Mutagenesis
#' 
#' The mutate_msd function designs the necessary set of primers for the desired mutations.
#' Note that you can also select TGG in saturation mutagenesis to apply the 22c trick.
#' 
#' @param input_sequence The sequence which should be modified. This is an object of type character containing the sequence. 
#' @param codon The desired type of MSD mutation [default: NDT]
#' @param prefix Additional nucleobases in 5' position of the recognition site [default: TT]
#' @param restriction_enzyme Recognition site sequence of the respective restriction enzyme [default: GGTCTC]
#' @param suffix Spacer nucleotides matching the cleavage pattern of the enzyme [default: A]
#' @param vector Four basepair overhangs complementary to the created overhangs in the acceptor vector  [default: c("AATG", "AAGC")]
#' @param replacements The desired substitutions. Can be a numeric vector or a list of character vectors with position and codon. See \href{articles/MSD3.html}{\code{vignette("MSD3", package = "goldenmutagenesis")}} for examples.
#' @param replacement_range  Maximum distance in amino acid residues between two randomization sites to be incoporated into a single primer (reverse, end of the fragment) - has a cascading effect for following mutations [default: 3]
#' @param binding_min_length The minimal threshold value of the length of the template binding sequence in amino acid residues [default: 4]
#' @param binding_max_length Maximal length of the binding sequence in amino acid residues [default: 9]
#' @param target_temp Melting temperature of the binding sequence in \code{print('\u00B0')}C [default: 60]
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
#' primers<-mutate_msd(input_sequence, prefix="TT" ,
#' restriction_enzyme=recognition_site_bsai, suffix="A", 
#' vector=c("AATG", "AAGC"), replacements=mutations, replacement_range=5,
#' binding_min_length=4 , binding_max_length=9, target_temp=60,
#' fragment_min_size=60 )
mutate_msd<-function(input_sequence, codon="NDT" ,prefix="TT" ,restriction_enzyme="GGTCTC", suffix="A", vector=c("AATG", "AAGC"), replacements, replacement_range=3, binding_min_length=4 ,binding_max_length=9, target_temp=60, fragment_min_size=100 ) {
  codon<-str_to_upper(codon)
  prefix<-str_to_upper(prefix)
  suffix<-str_to_upper(suffix)
  vector<-str_to_upper(vector)
  restriction_enzyme<-str_to_upper(restriction_enzyme)
  input_sequence<-str_to_upper(input_sequence)
  possible_codons<-c("NNN", "NNK", "NNS", "NDT", "DBK", "NRT", "VHG", "VRK", "NYC", "KST", "TGG")
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
  min_fragment<-3*binding_max_length
  prot_sequence<-seqinr::translate(sequence)
  primers<-vector("list")
  if(str_sub(vector[1], 2) == "ATG"){
    fragment_start<-2
  } else{
    fragment_start<-1
  }
  fragments<-make_fragments(mutations=replacements,seq=codon_seq,start=fragment_start, distance=replacement_range, fsize=fragment_min_size, buffer=2)
  
  primers<-vector("list", length(fragments))
  for(n in 1:length(fragments)){
    cur_fragment<-fragments[[n]]
    primers[[n]]<-vector("list", 2)
    if(n==1) { #the first primer does not need any existing overlap
      vector_f=vector[1]
      vector_r<-""
      overhang_f<-""
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[cur_fragment@stop-1], start=3), codon_seq[cur_fragment@stop], sep=""))))), collapse="")
      suffix_f<-suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop-2
      if(n==length(fragments)){#if there is just one single fragment - e.g. one mutation at the start or end
        overhang_r<-""
        stop_r<-cur_fragment@stop
        vector_r<-vector[2]
      }
    } 
    else if(n==length(fragments)) {
      vector_f<-""
      vector_r<-vector[2]
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-""
      suffix_f<-suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop
    }
    else {
      vector_f=""
      vector_r=""
      overhang_f<-paste(str_to_upper(s2c(paste(str_sub(codon_seq[fragments[[n-1]]@stop-1], start=3), codon_seq[fragments[[n-1]]@stop], sep=""))), collapse="")
      overhang_r<-paste(str_to_upper(comp(rev(s2c(paste(str_sub(codon_seq[fragments[[n]]@stop-1], start=3), codon_seq[fragments[[n]]@stop], sep=""))))), collapse="")
      suffix_f<-suffix
      suffix_r<-paste(comp(s2c(suffix), forceToLower = F), sep="", collapse = "")
      stop_r<-cur_fragment@stop-2
    }
    #forward
    if(length(cur_fragment@start_mutation)==0) {
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      temp_primer@binding_sequence<-paste(paste(codon_seq[cur_fragment@start:(cur_fragment@start+(binding_max_length-1))], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    else {
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_f, vector=vector_f, overhang=overhang_f)
      codon_seq[cur_fragment@start_mutation]<-codons[1:length(cur_fragment@start_mutation)]
      codons<-codons[-(1:length(cur_fragment@start_mutation))]
      temp_primer@extra<-paste(paste(codon_seq[cur_fragment@start:max(cur_fragment@start_mutation)], collapse = ""), sep="")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(max(cur_fragment@start_mutation)+1):((max(cur_fragment@start_mutation)+1)+binding_max_length-1)], collapse=""), sep="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=target_temp)
    }
    primers[[n]][[1]]<-temp_primer
    rm(temp_primer)
    #reverse
    if(length(cur_fragment@stop_mutation)==0){
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      temp_primer@binding_sequence<-paste(paste(codon_seq[(cur_fragment@stop-2-binding_max_length-1):stop_r], collapse=""), sep="")
      if(n!=length(fragments))
        temp_primer@binding_sequence<-paste(temp_primer@binding_sequence, str_sub(codon_seq[cur_fragment@stop-1], end=2) ,sep="")
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    else{
      temp_primer<-pc_msd(prefix=prefix ,restriction_enzyme=restriction_enzyme, suffix=suffix_r, overhang=overhang_r, vector=vector_r)
      codon_seq[cur_fragment@stop_mutation]<-codons[1:length(cur_fragment@stop_mutation)]
      codons<-codons[-(1:length(cur_fragment@stop_mutation))]
      temp_primer@extra<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)):stop_r], collapse=""), sep="")
      temp_primer@extra<-paste(temp_primer@extra, str_sub(codon_seq[cur_fragment@stop-1], end=2), sep="")
      temp_primer@extra<-paste(comp(rev(s2c(temp_primer@extra)), ambiguous = T,forceToLower = F), collapse = "")
      temp_primer@binding_sequence<-paste(paste(codon_seq[(min(cur_fragment@stop_mutation)-1-binding_max_length-1):(min(cur_fragment@stop_mutation)-1)], collapse=""), sep="")
      temp_primer@binding_sequence<-paste(str_to_upper(comp(rev(s2c(temp_primer@binding_sequence)))), collapse="")
      temp_primer<-sequence_length_temperature(temp_primer, primer_min=binding_min_length, target_temp=primers[[n]][[1]]@temperature)
    }
    primers[[n]][[2]]<-temp_primer
    rm(temp_primer)
  }
  #Check the primers with the checkprimer function:
  #Check for primers with same overlap
  #Replace them based on ? -> Primer without mutation/length of the primer in total?
  #It is easier to modify the exisiting primer 
  #If it is not possible to correct all overlaps -> return message with position for silent mutation
  primers<-check_primer_overhangs(primers, fragments, binding_min_length, target_temp)
  return(eps(oldsequence=input_sequence, primers=primers, newsequence=paste(codon_seq, collapse = ""), fragments=fragments))
}

#' Add a level to exisiting Primerset
#' 
#' This function replaces the prefix, the suffix and the restriction enzyme of a given Primerset to change the design to another level.
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
  prefix<-str_to_upper(prefix)
  restriction_enzyme<-str_to_upper(restriction_enzyme)
  suffix<-str_to_upper(suffix)
  vector<-str_to_upper(vector)
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

#' Prepare a Single Point Mutation Primerset to be used in Level 2
#' 
#' This function adds definied vector overhangs to Level 0 Primersets to express them in a Level 2 vector.
#' #'
#' @param primerset An exisiting Primerset (in Level 0)
#' @param vector Four basepair overhangs complementary to the created overhangs in the acceptor vector [default: c("AATG", "AAGC")]
#'
#' @return A Primerset prepared for expression in Level 2
#' @export
#'
#' @examples
#' #Load level 0 results of the SPM vignette
#' data(Point_Mutagenesis_BbsI_result)
#' primer_prepare_level(primers)
#' 
primer_prepare_level<-function(primerset, vector=c("AATG", "AAGC")){
  vector<-str_to_upper(vector)
  if(str_sub(vector[1], 2, 4)!="ATG") {
    warning("Working with unknown sequence for Level 2 vector")
    primerset@primers[[1]][[1]]@extra<-paste(primerset@primers[[1]][[1]]@extra, vector[1], sep="")
  }
  else if(str_sub(paste(primerset@primers[[1]][[1]]@extra,primerset@primers[[1]][[1]]@binding_sequence, sep=""), 1, 3) == "ATG") {
    primerset@primers[[1]][[1]]@extra<-paste(str_sub(vector[1], 1, 1),primerset@primers[[1]][[1]]@extra, sep="")
  } else {
    stop("The extra+binding_sequence of the primer did not start with ATG. Something went wrong. Please send a bug report to us.")
  }
  primerset@primers[[length(primerset@primers)]][[2]]@extra<-paste(vector[2], primerset@primers[[length(primerset@primers)]][[2]]@extra,  sep="")
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
