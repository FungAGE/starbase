# -b int
# --breaklen 	Distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)
# -c int
# --mincluster 	Minimum cluster length (default 65)
# -l int
# --minmatch 	Minimum length of an maximal exact match (default 20)

nucmer -b 1000 -c 25 -l 10 --prefix=ref_qry SQL/data/fna/ships/mycodb/mycodb.final.starships.part_altals1_s00058__+.fna SQL/data/fna/ships/mycodb/mycodb.final.starships.part_altalt7_s00064__+.fna -o ref_qry.coords

# lastz SQL/data/fna/ships/mycodb/mycodb.final.starships.part_altals1_s00058__+.fna SQL/data/fna/ships/mycodb/mycodb.final.starships.part_altalt7_s00064__+.fna