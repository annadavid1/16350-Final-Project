#!/bin/bash

g++ -DNUMBALLS=$1 planner.cpp -o planner.out \
-I/opt/homebrew/opt/sfml@2/include \
-L/opt/homebrew/opt/sfml@2/lib \
-lsfml-graphics -lsfml-window -lsfml-system \
-pthread

max=-1
balls=$1

for file in data/plans${balls}_*.csv; do
  [[ -e "$file" ]] || continue  # handles no matches case

  if [[ $file =~ data/plans${balls}_([0-9]+)\.csv ]]; then
    num=${BASH_REMATCH[1]}
    (( num > max )) && max=$num
  fi
done

next=$((max + 1))
planfilename="data/plans${balls}_${next}.csv"
exfilename="data/executed${balls}_${next}.csv"

./planner.out "$planfilename" "$exfilename"