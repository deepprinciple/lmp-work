hostname

# set variables to hold names of:
set OD = $DQS_O_WORKDIR                   # Original Directory
set LD = /scratch/bill/$JOB_ID            # Local Directory
set JN = $JOB_NAME:r                      # Job name
set JN = $JN.sos
set AN = $JN.tar                          # Archive name

cd $OD/$JN                               # Create the archive, compress it
tar cf ../$AN *
cd ..
gzip -9 $AN

mkdir $LD                                 # Create local directory, copy archive
/bin/cp $AN.gz $LD                       #    to it, remove original
/bin/rm $AN.gz

cd $LD                                    # Go to local directory, unpack and
gzip -cd $AN.gz | tar xof -               #    remove archive
/bin/rm $AN.gz

csh soscmd21                               # run
grep DelSOS f*sum >>delgsum
/home/bill/boss/miscexec/sumdelg
rm delgsum
#

tar cf ../$AN f* sumout                   # Create output archive, compress
cd ..                                      #   it and move to Original
gzip -9 $AN                               #   directory, remove Local
/bin/rm -rf $LD:t                         #   directory.
/bin/cp $AN.gz $OD/$JN
/bin/rm $AN.gz
cd $OD/$JN
gzip -cd $AN.gz | tar xf -
/bin/rm $AN.gz
