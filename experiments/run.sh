#!/usr/bin/bash

set -euo pipefail

for nn in casa pointrcnn pvrcnn; do
  for thresh in {00..80..10}; do
    ./Point_Cloud.jl "$nn" "$thresh" >> log.txt
  done
  for thresh in {50..70..1}; do
    ./Point_Cloud.jl "$nn" "$thresh" >> log.txt
  done
done
