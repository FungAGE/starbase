devtools::install_github("zachcp/msaR")
library(msaR)

prot<-"shoot/YRsuperfamRefs/tree/YRsuperfamRefs.mafft.faa"

# loading AA with ape. Can also use Biostrings
proteins <- ape::read.FASTA(prot, type="AA")
                            
# note the seqlogo will show up in your widget but 
# not the vignette static output
prot.html<-msaR(proteins, menu=F, overviewbox = F,  colorscheme = "clustal")

htmlwidgets::saveWidget(prot.html, app_sys("captain.prot.html"))

nuc<-"shoot/YRsuperfamRefs/tree/YRsuperfamRefs.mafft.fa"

# display the MSA.
nuc.html<-msaR(nuc, menu=F, overviewbox = T)
htmlwidgets::saveWidget(nuc.html, app_sys("captain.nuc.html"))
