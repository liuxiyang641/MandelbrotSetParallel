#!/bin/bash
# usage: 
# chmod +x run.sh
# ./run.sh <二进制文件>

file=$1
proc=3
thd_per_proc=4
real_lower=-4
real_upper=4
imag_lower=-4
imag_upper=4
w=4000
h=4000
output_path="./${file}.png"
if [ "$file" == "sequential" ]; then
    ./${file} ${thd_per_proc} ${real_lower} ${real_upper} ${imag_lower} ${imag_upper} ${w} ${h} ${output_path}
else
    mpirun -n ${proc} ./${file} ${thd_per_proc} ${real_lower} ${real_upper} ${imag_lower} ${imag_upper} ${w} ${h} ${output_path}
fi