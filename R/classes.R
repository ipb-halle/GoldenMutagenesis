ps<-setClass("Primerset", slots=c(oldsequence="character", primers="list", newsequence="character"))
pc<-setClass("Primer", slots=c(prefix="character" ,restriction_enzyme="character", suffix="character", vector="character", overhang="character", extra="character" ,binding_sequence="character", temperature="numeric", difference="numeric"))
setMethod("initialize", "Primer",
          function(.Object, prefix="", restriction_enzyme="", suffix="", vector=c("", ""), overhang="", extra="", binding_sequence="", temperature=60, difference=0,...) {
            .Object<-callNextMethod(.Object, ...)
            .Object@prefix<-prefix
            .Object@restriction_enzyme<-restriction_enzyme
            .Object@suffix<-suffix
            .Object@vector<-vector
            .Object@overhang<-overhang
            .Object@extra<-extra
            .Object@binding_sequence<-binding_sequence
            .Object@temperature<-temperature
            .Object@difference<-difference
            .Object
          }
          )
pc_msd<-setClass("Primer_MSD", contains="Primer")
pc_spm<-setClass("Primer_SPM", contains="Primer")
fragment<-setClass("Fragment", slots=c(start="numeric", stop="numeric", start_mutation="vector", stop_mutation="vector"))
eps<-setClass("Extended_Primerset", contains="Primerset", slots=c(fragments="list"))
