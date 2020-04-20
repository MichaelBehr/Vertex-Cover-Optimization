rm GeneratedGraphs/alternate15.txt 
for j in $(seq 1 10)
do
    GENERATEDGRAPH=$(./graphGen 15)
    for k in $(seq 1 10)
    do
        echo "$GENERATEDGRAPH" >> GeneratedGraphs/alternate15.txt 
    done
done
