#!/bin/bash

g++-9.1.0 -Wall -g -O3 -fopenmp -std=c++17 deliverable1.cpp -o deliverable1

module load perf

# Parametri di test
SCHEDULES=("static" "dynamic" "guided")
CHUNKS=(10 100 1000)
THREADS=(1 2 4 8 16 32 64)

# Numero di ripetizioni per calcolare la media
REPEAT=2

OUTFILE="results.csv"


echo "threads,schedule,chunk,avg_time(s),l1_cache_load,l1_cache_miss,l1_miss_rate,ll_cache_load,ll_cache_miss,ll_miss_rate" > "$OUTFILE"

for t in "${THREADS[@]}"; do
  export OMP_NUM_THREADS=$t
  
  for s in "${SCHEDULES[@]}"; do
    for c in "${CHUNKS[@]}"; do
        
        total_time=0.0
        total_l1_load=0.0
        total_l1_miss=0.0
        total_ll_load=0.0
        total_ll_miss=0.0 
             
        for ((i=1; i<=REPEAT; i++)); do

          result=$( (perf stat -e L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses ./deliverable1 $s $c) 2>&1 )
          
          #the printf of time has to be written exactly like that: printf("Time %e",time)
          time_val=$(echo "$result" | grep "Time" | awk '{printf "%f", $2}')
          
          l1_cache_load=$(echo "$result" | grep "L1-dcache-loads" | awk '{print $1}' | tr -d ',')
          l1_cache_miss=$(echo "$result" | grep "L1-dcache-load-misses" | awk '{print $1}' | tr -d ',')
          ll_cache_load=$(echo "$result" | grep "LLC-loads" | awk '{print $1}' | tr -d ',')
          ll_cache_miss=$(echo "$result" | grep "LLC-load-misses" | awk '{print $1}' | tr -d ',')
          
        
          # Se una delle variabili è vuota, metti 0 per sicurezza
          time_val=${time_val:-0}
          l1_cache_load=${l1_cache_load:-0}
          l1_cache_miss=${l1_cache_miss:-0}
          ll_cache_load=${ll_cache_load:-0}
          ll_cache_miss=${ll_cache_miss:-0}
          
          #echo "DEBUG [$t $s $ c $i]: $l1_cache_load, $l1_cache_miss, $ll_cache_load, $ll_cache_miss"
  

          #total data
          total_time=$(echo "$total_time + $time_val" | bc -l)
          total_l1_load=$(echo "$total_l1_load + $l1_cache_load" | bc -l)
          total_l1_miss=$(echo "$total_l1_miss + $l1_cache_miss" | bc -l)
          total_ll_load=$(echo "$total_ll_load + $ll_cache_load" | bc -l)
          total_ll_miss=$(echo "$total_ll_miss + $ll_cache_miss" | bc -l)
          
          
        done
        
        #average data
        avg_time=$(printf "%.10f" "$(echo "$total_time / $REPEAT" | bc -l)")
        avg_l1_load=$(printf "%.2f" "$(echo "$total_l1_load / $REPEAT" | bc -l)")
        avg_l1_miss=$(printf "%.2f" "$(echo "$total_l1_miss / $REPEAT" | bc -l)")
        avg_ll_load=$(printf "%.2f" "$(echo "$total_ll_load / $REPEAT" | bc -l)")
        avg_ll_miss=$(printf "%.2f" "$(echo "$total_ll_miss / $REPEAT" | bc -l)")
        
        
        avg_l1_miss=${avg_l1_miss:-0}
        avg_l1_load=${avg_l1_load:-0}
        avg_ll_miss=${avg_ll_miss:-0}
        avg_ll_load=${avg_ll_load:-0}
        
        #echo "DEBUG: "$avg_l1_miss" "$avg_l1_load" "$avg_ll_miss" "$avg_ll_load"
    
         # calcolo percentuali di cache miss
        if (( $(echo "$avg_l1_load > 0" | bc -l) )); then
            l1_miss_rate=$(printf "%.2f" "$(echo "($avg_l1_miss / $avg_l1_load)*100" | bc -l)")
        else
            l1_miss_rate=0
        fi
        
        if (( $(echo "$avg_ll_load > 0" | bc -l) )); then
            ll_miss_rate=$(printf "%.2f" "$(echo "($avg_ll_miss / $avg_ll_load)*100" | bc -l)")
        else
            ll_miss_rate=0
        fi
        
        #l1_miss_rate=0
        #ll_miss_rate=0
        
        
        #save on csv
         echo "$t, $s, $c, $avg_time, $avg_l1_load, $avg_l1_miss, $l1_miss_rate, $avg_ll_load, $avg_ll_miss, $ll_miss_rate" >> "$OUTFILE"
                
    done
  done
done