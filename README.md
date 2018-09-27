# scalemb

## Info:
Implementation of [1]: "Billion-scale Network Embedding with Iterative Random Projection", Ziwei Zhang, Peng Cui, Haoyang Li, Xiao Wang and Wenwu Zhu, ICDM2018.

NOTE THAT THIS IS NOT THE IMPLEMENTATION OF THE AUTHORS AND IT IS NOT GUARANTEED TO WORK.
Implementation of the authors is here: https://github.com/ZW-ZHANG/RandNE

## to compile:

gcc scalemb.c -o scalemb -O9 -lm

## to execute:

./scalemb net.txt emb.txt d q a0 a1 ... aq
- net.txt should contain on each line: "i j\n" that is the input directed graph.
- emb.txt contains the resulting embedding (d floats on each line)
- d is the dimention of the embeding
- q is the order of the embbeding
- ak's are the coeficient of A^k matrix such as defined in [1]

## Reference:
[1] : https://papers-gamma.link/paper/110

## First contributor:
Maximilien Danisch  
Septembre 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com

