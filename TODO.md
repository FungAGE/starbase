# TODO

- [ ] change name? existing tools/databases named "starbase"... one which is mRNA, etc. database
  - "starfleet", "armada"
- [x] "SSA" why did we change?
  - I didn't really consult with anyone before changing this

## data integrity
- [x] sequences without accessions? accessions without sequences? check zenodo version, Andres found inconsistencies...
- [x] wiki table has issues with metadata
  - [x] assembly accessions should not have the zenodo doi in there...
  - [x] genome information is usually not visible at all
- [ ] check accessions that have lots of ships collapsed in them
  - for example, asp fum, lots under a single accession
    - these are missing taxonomy ids
    - we have strains for these taxonomy entries, they should not have the same tax id anymore
- [x] accession version tags, we should keep them
- [x] unique identifiers should be a new level of accessions

SSAXXXXXX.X : collapsed navis/haplotype identifiers
SSBXXXXXX.X : individual ships


- [x] re-organize modals to reflect the new organization of the multiple accessions
  - have a table to metadata with rows for each individual ships with their specific genome, taxonomy info
- synteny viewer should show the genome accession, instead of the SSA twice
- flag one unique entry in each SSA accession as the "reference" 
  - this should also be designated at the "SSB" accession level, arbitrarily chosen if there are multiple same ships in the same genome
- switch to more direct comparisons for synteny viewer, the table is overwhelming for a user
- blast page
  - does a sequence with many N's result in a truncated blast result(s)?
- fix horizontal scroll on the modals: it's making the close button hard to click
- submission page
  - instead of a text field for annotation method, have radio buttons, and a custom field
  - discourage people from just submitting raw starfish output



- [x] handle operations with accession tags, when multiple sequences are returned by the accession tag

- synteny
  - should viz all sequences, even if they are identical (it's still useful)

- blast
  - really only care about one representative from each accession

- wiki
  - it should be clear that the table is displaying a list of accessions
  - from the table we can get information about the individual ships within that accession

- [x] sourmash signatures should be created at the same time as other databases?
  - [x] use existing sourmash signatures at the time of similarity comparison, rather than creating it fresh