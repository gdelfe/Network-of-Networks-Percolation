#!/bin/sh
rm p2.txt
#
ind=9135
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
ind=9165
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
ind=9565
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
ind=9565
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
ind=9735
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
ind=9765
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt

for ind in $(seq 7000 100 9900)
do
./perc J_th_0.$ind.txt $1 0 | awk -v ind="$ind" '{print "0."ind, $0}' >>  p2.txt
done