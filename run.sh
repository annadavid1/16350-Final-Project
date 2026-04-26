#!/bin/bash

g++ -DNUMBALLS=$1 planner.cpp -o planner.out \
-I/opt/homebrew/opt/sfml@2/include \
-L/opt/homebrew/opt/sfml@2/lib \
-lsfml-graphics -lsfml-window -lsfml-system \
-pthread

max=-1
balls=$1

mkdir -p data/ball${balls}
mkdir -p data/ball${balls}

for file in data/ball${balls}/results_*.csv; do
  [[ -e "$file" ]] || continue  # handles no matches case

  if [[ $file =~ data/ball${balls}/results_([0-9]+)\.csv ]]; then
    num=${BASH_REMATCH[1]}
    (( num > max )) && max=$num
  fi
done

next=$((max + 1))
filename="data/ball${balls}/results_${next}.csv"

./planner.out "$filename"