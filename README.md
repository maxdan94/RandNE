
## Info:
Implementation of the scalable graph embeding method suggested in [1]

NOTE THAT THIS IS NOT THE IMPLEMENTATION OF THE AUTHORS AND IT IS NOT GUARANTEED TO WORK.
Implementation of the authors is/will be here: https://github.com/ZW-ZHANG/RandNE.

## to compile:

gcc RandNE.c -o RandNE -O9 -lm

## to execute:

./RandNE net.txt emb.txt d q a0 a1 ... aq
- net.txt should contain on each line: "u v\n" that is the input undirected graph.
- emb.txt contains the resulting embedding (d floats on each line)
- d is the dimention of the embeding
- q is the order of the embbeding
- ak's are the coeficient of A^k matrix such as defined in [1]

## Reference:
[1] : "Billion-scale Network Embedding with Iterative Random Projection", Ziwei Zhang, Peng Cui, Haoyang Li, Xiao Wang and Wenwu Zhu, ICDM2018. https://papers-gamma.link/paper/110

## First contributor:
Maximilien Danisch  
Septembre 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com

