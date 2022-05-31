#!/bin/bash 

begin=11
end=20
for i in $(seq $begin $end); do
    cat 03_eq2smd.sh | sed -e "s/#{ID}/${i}/g" > job.sh
    grep "id=" job.sh 
    pjsub job.sh
done
rm job.sh
