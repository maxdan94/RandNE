
## Info:
Implementation of the scalable graph embedding method proposed in [1]

NOTE THAT THIS IS NOT THE IMPLEMENTATION OF THE AUTHORS OF [1].  
The implementation of the authors of [1] is here: https://github.com/ZW-ZHANG/RandNE.

## to compile:

gcc RandNE.c -o RandNE -O9 -lm

## to execute:

./RandNE net.txt emb.txt d q a0 a1 ... aq
- net.txt should contain on each line: "u v\n" that is the input undirected graph.
- emb.txt contains the resulting embedding (d floats on each line)
- d is the dimension of the embedding
- q is the order of the embedding
- ak's are the coefficient of A^k matrix such as defined in [1]

## running time on a commodity machine:

On http://snap.stanford.edu/data/com-Orkut.html
- d=32 and q=3: 5 minutes
- d=32 and q=5: 10 minutes
- d=128 and q=3: 20 minutes


## Reference:
[1] : "Billion-scale Network Embedding with Iterative Random Projection", Ziwei Zhang, Peng Cui, Haoyang Li, Xiao Wang and Wenwu Zhu, ICDM2018. https://papers-gamma.link/paper/110

## First contributor:
Maximilien Danisch  
Septembre 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com


