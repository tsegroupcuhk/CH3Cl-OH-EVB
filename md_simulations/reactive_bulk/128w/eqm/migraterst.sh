#!/bin/bash

for i in {106..305}
do
filename=restart/nvtrun-$i.rst
cp $filename ../config_bank/$((i-105)).rst

input_file="../config_bank/$((i-105)).rst"
# Remove the line with <Parameters>
sed -i '/<Parameters/d' "$input_file"
# Reset stepCount to 0
sed -i 's/stepCount="[^"]*"/stepCount="0"/' "$input_file"

done
