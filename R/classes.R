ps<-setClass("Primerset", slots=c(oldsequence="character", primers="list", newsequence="character"))
pc<-setClass("Primer", slots=c(prefix="character" ,restriction_enzyme="character", suffix="character", vector="character", overhang="character" ,binding_sequence="character", temperature="numeric", difference="numeric"))
pc_msd<-setClass("Primer MSD", contains="Primer", slots=c(NDT="character"))
fragment<-setClass("Fragment", slots=c(start="numeric", stop="numeric", start_mutation="vector", stop_mutation="vector"))