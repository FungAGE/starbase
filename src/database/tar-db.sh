#!/bin/bash
cd src/database || exit
tar -cf db.tar db/captain/tyr/faa/blastdb db/captain/tyr/faa/hmm/combined.hmm db/captain/tyr/fna/hmm/combined.hmm db/ships/fna/blastdb db/starbase.sqlite db/submissions.sqlite db/telemetry.sqlite
gzip -f db.tar