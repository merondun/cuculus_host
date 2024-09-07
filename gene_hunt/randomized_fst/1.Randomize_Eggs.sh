#randomize
awk '{print $2}' W7.popfile > temp2

for i in {1..10}
do
    # Shuffle lines in the file
    awk '{print $1}' W7.popfile | shuf > temp.list

    # Add $p1 to first 5 lines and $p2 to the next 5
    paste temp.list temp2 > "rf${i}.popfile"
done

