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
gcc scalemb.c -o scalemb -O9 -lm

to execute:
./scalemb net.txt emb.txt d q a0 a1 ... aq
- net.txt should contain on each line: "i j\n" that is the input graph.
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

#define NLINKS 10000000 //maximum number of links, will increase if needed


// graph datastructure:
typedef struct {
	unsigned long s;//source node
	unsigned long t;//target node
} edge;//not directed!

//sparse graphe structure
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	edge *el;//edge list
} sparse;

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the graph
sparse* readedgelist(char* edgelist){
	unsigned long i, e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	g->el=malloc(e1*sizeof(edge));
	FILE *file;
	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	while (fscanf(file,"%lu %lu\n", &(g->el[g->e].s), &(g->el[g->e].t))==2) {
		g->n=max3(g->n,g->el[g->e].s,g->el[g->e].t);
		if (g->e++==e1) {
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

double scallarproduct(unsigned long n, double *v1, double *v2){
	unsigned long i;
	double s=0;
	for (i=0;i<n;i++){
		s+=v1[i]*v2[i];
	}
	return s;
}

double *GramSchmidt(unsigned long n, unsigned d){
	double s;
	unsigned long m=n*d,d1,d2,m1,m2,i;
	double *vect=malloc(m*sizeof(double));

	srand(time(NULL));
	for (i=0;i<m;i++){
		vect[i]= ((double)rand()/(RAND_MAX));//supposed to be the normal distribution with sigma=1/d
	}
	for (d1=0;d1<d;d1++){
		m1=n*d1;
		for (d2=0;d2<d1;d2++){
			m2=n*d2;
			s=scallarproduct(n,vect+m1,vect+m2);
			for (i=0;i<n;i++){
				vect[i+m1]-=vect[i+m2]*s;
			}
		}
		normalize(n,vect+m1);
	}
	return vect;
}

void prod(sparse* g, double* v1, double* v2){
	unsigned long i;
	bzero(v2,sizeof(double)*g->n);
	for (i=0;i<g->e;i++){
		v2[g->el[i].s]+=v1[g->el[i].t];
		v2[g->el[i].t]+=v1[g->el[i].s];
	}
}

void add(unsigned long n, double a, double* v_in, double* v_out){
	unsigned long i;
	for (i=0;i<n;i++){
		v_out[i]+=a*v_in[i];
	}
}


//cf algo 1 of [1]
double *RandNE(sparse *g,unsigned d,unsigned q, double* a){
	unsigned long n=g->n;//number of nodes
	unsigned i,j;
	double *r1=GramSchmidt(n,d),*r2=malloc(n*d*sizeof(double)),*r3;
	double* emb=calloc(n*d,sizeof(double));//resulting embedding
	if (a[0]>0.){
		add(n*d,a[0],r1,emb);
	}

	for (i=1;i<q+1;i++){
		for (j=0;j<d;j++){
			prod(g,r1+n*j,r2+n*j);
		}
		r3=r1;
		r1=r2;
		r2=r3;
		if (a[i]>0.){
			add(n*d,a[i],r1,emb);
		}
	}

	free(r1);
	free(r2);

	return emb;
}

void printres(FILE* file,unsigned long n,unsigned d,double *emb){
	unsigned long i,j;
	for (i=0;i<n;i++){
		fprintf(file,"%le",emb[i]);
		for (j=1;j<d;j++){
			fprintf(file," %le",emb[i+j*n]);
		}
		fprintf(file,"\n");
	}
}


//./scalemb net.txt emb.txt d q a0 a1 ... aq

int main(int argc,char** argv){
	sparse* g;//input graph
	unsigned d,q;//dimension and order
	double *a;//coefficients
	unsigned i;
	double* emb;//resulting embedding

	FILE* file;

	time_t t1,t2,t3;
	t1=time(NULL);

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Number of nodes = %lu\n",g->n);
	printf("Number of edges = %lu\n",g->e);

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

	printf("Computing the embeding\n");

	emb=RandNE(g, d, q, a);

	printf("Printing embedding in file %s\n",argv[2]);
	file=fopen(argv[2],"w");
	printres(file,g->n,d,emb);
	fclose(file);

	freegraph(g);
	free(emb);
	return 0;
}

