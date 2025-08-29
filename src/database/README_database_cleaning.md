"""
# Database cleaning tasks
## correction of `accession_tag`
### checking reverse compliment of sequence
- for each sequence:
1. generate md5sums for both normal and reverse compliment
2. generate list of accessions with identical md5sums for normal and reverse complimented sequences
3. check for existing sequences in the database that are nested, and therefore should be grouped under the same accession
6. if the list is longer than 1, consolidate the accessions based on the lowest number

### checking of nesting and updating accession version
- for each sequence:
1. look for existing seqeuences that are nested in other sequences
2. create a list of these nested sequences and the accessions of both seqeuences
3. consolidate nested sequences under the accession of the longer sequence (i.e. the sequence that they are nested within)

- once accessions have been corrected, proceed with other database cleaning tasks

## check genome table
- we need to add genome information where it is missing

## check taxonomic table
- check that taxonomic information is internally consistent

## check genomic features
- check that coordinates make sense and are linked to genomes

## after all previous steps have been completed
- update foreign keys and check consistency

## general schema check?
- can we check if database entries violate schema?
"""