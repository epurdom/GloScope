#!/bin/bash
# Run this file in bash with this command:  ./filename
cd inst/extdata
HOST=ftp.ebi.ac.uk
USER=anonymous
ftp -pinv $HOST <<EOF
user $USER
cd biostudies/nfs/E-MTAB-/026/E-MTAB-10026/Files
binary
mget "covid_portal_210320_with_raw.h5ad"
disconnect
bye
EOF
