#!/bin/bash
for i in {1..100}
do
sbatch run_1rep_264nodes.sh $i
done
