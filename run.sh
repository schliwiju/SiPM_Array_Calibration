#!/bin/bash
#usage: ./run.sh runNr
#       where runNr is a run number from runlist

here=`pwd`
if [ ! -d "$here/runs" ]; then
  mkdir $here/runs
fi


runNr=$1
echo $runNr
while read line
do
    [[ $line = \#* ]] && continue
    lineArr=($line)
    if [ "${lineArr[0]}" = "$runNr" ]; then 
      runName=${lineArr[1]}
      mkdir $here/runs/$runName
      if [ ! -e $here/runs/$runName/$runName.list ]; then
        ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
      fi
      time $here/read $here/runs/$runName/$runName.list $here/data/$runName/ $here/runs/$runName/out.root ${lineArr[0]} ${lineArr[2]} ${lineArr[3]} ${lineArr[4]} ${lineArr[5]}
    fi
done < ./runlist



# #valgrind --trace-children=yes --tool=massif time ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root  
# #run for memory check
# #HEAPCHECK=normal LD_PRELOAD=/usr/lib/libtcmalloc.so ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root
