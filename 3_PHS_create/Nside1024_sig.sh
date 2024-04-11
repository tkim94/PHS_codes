#!/bin/csh
 
i=1


for (( j=1; j<=20; j+=1  ))
do

# go and run MadGraph
#cd /afs/crc.nd.edu/user/t/tkim12/Work/CMB_ML/Eventcode/single_pair_sigbkg_gen/eta160/g6

cat Nside1024_sig.sub | sed     -e "s/ number/ $j/" > Nside1024_sig_$j.sub

qsub Nside1024_sig_$j.sub

rm Nside1024_sig_$j.sub

((i+=1))
done