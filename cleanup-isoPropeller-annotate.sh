#!/bin/sh

# 14.08.2023 14:39:28 EDT

for i in ????-??-??_freeze*postfilt*/
do
   if [ -d "${i}/01_sqanti3/" ] 
   then
      echo "Starting cleanup of '${i}'"
   
      rm -f ${i}/*/core.*
      
      # Cleanup 01_sqanti3
      rm -rf ${i}/01_sqanti3/GMST/
      rm -rf ${i}/01_sqanti3/RTS/
      for file in ${i}/01_sqanti3/*_corrected.* ${i}/01_sqanti3/*_junctions.txt
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
      
      # Cleanup/compress 05_cpatv3
      for file in ${i}/05_cpatv3/*.gtf ${i}/05_cpatv3/*ORF_seqs.fa
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
      
      # Cleanup/compress 06_interpro
      for file in ${i}/06_interpro/*_corrected.*
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
      
      # Cleanup/compress 09_niap_asef
      for file in ${i}/09_niap_asef/*.gtf
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
   
      # Cleanup/compress 11_transdecoder
      rm -rf ${i}/11_transdecoder/transdecoder_dir/
      for file in ${i}/11_transdecoder/*.transdecoder.*
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
      
      # Cleanup/compress 12_tracks
      for file in ${i}/12_tracks/*.gtf
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
     
      # Cleanup/compress pipeline logs
      for file in ${i}/pipeline*.log 
      do
         submitjob 1 -c 2 -m 5 -q express -P acc_pintod02c gzip ${file}
      done
   else
      echo "Skipping folder '${i}': Path does not exist."
   fi
done
