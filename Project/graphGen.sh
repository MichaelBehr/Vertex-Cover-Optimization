#!/bin/bash

# Generate a bunch of graphs, save them in a graph file.
for i in $(seq 5 1 20)
do
    rm GeneratedGraphs/$i.txt
    for j in $(seq 1 10)
    do
        GENERATEDGRAPH=$(./graphGen $i)
        for k in $(seq 1 10)
        do
            echo "$GENERATEDGRAPH" >> GeneratedGraphs/$i.txt 
        done
    done
done

for i in $(seq 25 5 50)
do
    rm GeneratedGraphs/$i.txt
    for j in $(seq 1 10)
    do
        GENERATEDGRAPH=$(./graphGen $i)
        for k in $(seq 1 10)
        do
            echo "$GENERATEDGRAPH" >> GeneratedGraphs/$i.txt 
        done
    done
done
