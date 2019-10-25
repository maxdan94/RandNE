/*
Maximilien Danisch
Septembre 2018
http://bit.ly/danisch
maximilien.danisch@gmail.com

Info:
Implementation of [1]: "Billion-scale Network Embedding with Iterative
Random Projection" ICDM2018.

NOTE THAT THIS IS NOT THE IMPLEMENTATION OF THE AUTHORS.

to compile:
gcc RandNE.c -o RandNE -O9 -lm

to execute:
./RandNE net.txt emb.txt d q a0 a1 ... aq
- net.txt should contain on each line: "u v\n" that is the input graph.
- emb.txt contains the resulting embedding (d floats on each line)
- d is the dimention of the embeding
- q is the order of the embbeding
- ak's are the coeficient of A^k matrix such as defined in [1]

Reference:
[1]: https://papers-gamma.link/paper/110

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <time.h>

#define NLINKS 100000000 //maximum number of links, will increase if needed

// graph datastructure:
typedef struct {
	unsigned long s;//source node
	unsigned long t;//target node
} edge;//not directed!

//sparse graphe structure
typedef struct {
	unsigned long n;//number of nodes
	unsigned long long e;//number of edges
	edge *el;//edge list
} sparse;

//compute the maximum of three unsigned long ints
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the graph
sparse* readedgelist(char* edgelist){
	unsigned long long e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	g->el=malloc(e1*sizeof(edge));
	FILE *file;
	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	while (fscanf(file,"%lu %lu\n", &(g->el[g->e].s), &(g->el[g->e].t))==2) {
		g->n=max3(g->n,g->el[g->e].s,g->el[g->e].t);
		if (++(g->e)==e1) {
			e1+=NLINKS;
			g->el=realloc(g->el,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->el=realloc(g->el,g->e*sizeof(edge));

	return g;
}

//free the graph stucture
void freegraph(sparse *g){
	free(g->el);
	free(g);
}

//normalize the vector vect such that ||vect||_2=1.
void normalize(unsigned long n, double* vect){
	unsigned long i;
	double s=0;
	for (i=0;i<n;i++){
		s+=vect[i]*vect[i];
	}
	s=sqrt(s);
	for (i=0;i<n;i++){
		vect[i]/=s;
	}
}

//compute the scalar product v1.v2
double scallarproduct(unsigned long n, double *v1, double *v2){
	unsigned long i;
	double s=0;
	for (i=0;i<n;i++){
		s+=v1[i]*v2[i];
	}
	return s;
}

//Return d orthonormal vectors of dimension n using https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
//The d n-vecors are concatenated into a single vector of length n*d.
double *GramSchmidt(unsigned long n, unsigned d){
	double s;
	unsigned long long n_ull=(unsigned long long)n,	m=n_ull * d;//n_ull*d can exeed 4,294,967,295 (largest number stored with unsigned long long)
	unsigned long long m1,m2,i;
	unsigned long j;
	unsigned d1,d2;
	double *vect=malloc(m*sizeof(double));

	srand(time(NULL));

	for (i=0;i<m;i++){
		vect[i]= ((double)rand()/(RAND_MAX));//CAREFUL!!!!!! According to the version of the authors, this is supposed to be the normal distribution with sigma=1/d
	}
	for (d1=0;d1<d;d1++){
		m1=n_ull * d1;
		for (d2=0;d2<d1;d2++){
			m2=n_ull * d2;
			s=scallarproduct(n,vect+m1,vect+m2);
			for (j=0;j<n;j++){
				vect[j+m1]-=vect[j+m2]*s;
			}
		}
		normalize(n,vect+m1);
	}
	return vect;
}

//let A be the adjacency matrix of graph g. Puts in v2 the product of the matrix A by v2. That is: v2=A*v1.
void prod(sparse* g, double* v1, double* v2){
	unsigned long long i;
	bzero(v2,sizeof(double)*g->n);
	for (i=0;i<g->e;i++){//better using a classical martix vector multiplication to make it parallel: v2[i]=sum_{(i,j)} v1[i].
		v2[g->el[i].s]+=v1[g->el[i].t];
		v2[g->el[i].t]+=v1[g->el[i].s];
	}
}

//add a*v_in to v_out. V_in and v_out can be very larde. n is their size.
void add(unsigned long long n, double a, double* v_in, double* v_out){
	unsigned long long i;
	for (i=0;i<n;i++){
		v_out[i]+=a*v_in[i];
	}
}

//cf algo 1 of [1]
double *RandNE(sparse *g,unsigned d,unsigned q, double* a){
	unsigned long n=g->n;//number of nodes
	unsigned long long n_ull=(unsigned long long)n;
	unsigned i,j;
	double *r1=GramSchmidt(n,d),*r2=malloc(n*d*sizeof(double)),*r3;
	double* emb=calloc(n_ull*d,sizeof(double));//resulting embedding
	if (a[0]>0.){
		add(n_ull*d,a[0],r1,emb);
	}

	for (i=1;i<q+1;i++){
		for (j=0;j<d;j++){
			prod(g,r1+n_ull*j,r2+n_ull*j);
		}
		r3=r1;
		r1=r2;
		r2=r3;
		if (a[i]>0.){
			add(n_ull*d,a[i],r1,emb);
		}
	}

	free(r1);
	free(r2);

	return emb;
}

//printing the result
void printres(FILE* file,unsigned long n,unsigned d,double *emb){
	unsigned long i;
	unsigned j;
	unsigned long long n_ull=(unsigned long long)n;
	for (i=0;i<n;i++){
		fprintf(file,"%le",emb[i]);
		for (j=1;j<d;j++){
			fprintf(file," %le",emb[i+j*n_ull]);
		}
		fprintf(file,"\n");
	}
}


//./RandNE net.txt emb.txt d q a0 a1 ... aq

int main(int argc,char** argv){
	sparse* g;//input graph
	unsigned d,q;//dimension and order
	double *a;//coefficients
	unsigned i;
	double* emb;//resulting embedding

	FILE* file;

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Number of nodes = %lu\n",g->n);
	printf("Number of edges = %llu\n",g->e);

	printf("Dimensions of the embedding = %s\n",argv[3]);
	d=atoi(argv[3]);

	printf("Order of the embedding = %s\n",argv[4]);
	q=atoi(argv[4]);

	a=malloc((q+1)*sizeof(double));
	printf("Coefficients:\n");
	for (i=5;i<6+q;i++){
		printf("a_%d = %s\n",i-5,argv[i]);
		a[i-5]=atof(argv[i]);
	}

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing the embeding\n");

	emb=RandNE(g, d, q, a);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Printing embedding in file %s\n",argv[2]);
	file=fopen(argv[2],"w");
	printres(file,g->n,d,emb);
	fclose(file);

	freegraph(g);
	free(emb);
	free(a);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}

