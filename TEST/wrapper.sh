#!/bin/bash -l 

nvidia-cuda-mps-control -d

IFS=',' read -r -a gpus <<< "$CUDA_VISIBLE_DEVICES"

echo $gpus
let num_gpus=${#gpus[@]}
echo "Number of visible gpus:  $num_gpus"

let mod_rank=$((OMPI_COMM_WORLD_LOCAL_RANK % $num_gpus))
$echo $mod_rank

export CUDA_VISIBLE_DEVICES=$mod_rank
$echo $CUDA_VISIBLE_DEVICES

$*
