# classifcation pipeline that should be used:
# - classify query sequences for blast page, display results
# - classify submitted sequences from submission page, assign an accession, and input into submissionsdatabase

# accession format:
# - normal ship accession: SBS123456
# - updated ship accession: SBS123456.1

########################################################
# assigning accessions
########################################################
# part 1: check if the sequence is already in the database
# part 2: if it is completely identical, almost identical, or contained within an existing sequence (also flag for review), assign the normal accession
# part 3: if it is novel, assign a new accession

# workflow:
# first check for exact matches
# then check for contained within matches
# then check for almost identical matches
# if no matches, assign a new accession

# if a sequence is a truncated version of a longer sequence, assign the longer sequence accession, flag for review

# fast method for checking for exact matches:
# md5 hash of sequence

# method for checking for contained/highly similar sequences:
# k-mers

########################################################
# classification pipeline
########################################################
# part 1: sequence similarity -> haplotype assignment
# part 2: cargo jaccard scores -> navis assignment
# part 3: captain gene similarity (hmmsearch) -> family assignment

# use methods from starfish pipeline