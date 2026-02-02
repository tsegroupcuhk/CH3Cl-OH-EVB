#!/bin/bash

count=150

echo "Generating folders for $count pulls ..."

for ((i=1; i<=count; i++))
do
  idx=$i
  pull_folder="pull-$idx"
  mkdir "$pull_folder"
  cp -r template/* "$pull_folder"
  sed -i "s/1.rst/${idx}.rst/g" "$pull_folder/runSMD.py"
done

