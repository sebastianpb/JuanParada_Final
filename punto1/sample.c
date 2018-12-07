#include <stdio.h>
#include <stdlib.h>

/*Incluye modulo math*/
#include <math.h>
#define PI 3.14159265f

double randn (double mu, double sigma){

	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1)
	{
	call = !call;
	return (mu + sigma *(double) X2);
	}

	do
	{
	U1 = -1 + ((double)rand() / RAND_MAX) * 2;
	U2 = -1 + ((double)rand() / RAND_MAX) * 2;
	W = pow(U1, 2) + pow(U2, 2);
	}
	while (W >= 1 || W == 0);

	mult = sqrt ((-2 * log(W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;
	
	return (mu + sigma * (double)X1);
}

// Devuelve gaussiana centrada en mu y con desviacion sigma
double gauss (double x, double mu, double sigma){
  double c;
  c=(1.0/(sigma*sqrt(2*PI)));
	return c*exp(-0.5*(pow(((x-mu)/sigma),2)));
}



int main(int argc, char *argv[]){

	int i;
	double aleatorio1;
	double mu,sigma,rand_num;
	char str_N;
  //str_N=argv[1];
  //printf("%s\n", argv[1]);

	char *p;
	int N=atoi(argv[1]);
	mu=atof(argv[2]);
	sigma=atof(argv[3]);
  //i<nt N = strtol(argv[1], &p, 10); //tamaño del puntero con los numeros de la distribucion

	int N_iter = 40;
	double lista [N_iter];
	double * plista;
	plista=lista;
	*(plista)=((double)rand()/ (double)RAND_MAX);
	double sigma_delta = 1.0;


	for(i=1;i<N_iter;i++){
      double propuesta  = *(plista+i-1) + gauss(1-2*((double)rand()/ (double)RAND_MAX),0,sigma_delta);
      double f1=gauss(propuesta,mu,sigma);
      double f2=gauss(*(plista+i-1),mu,sigma);

      double r=f1/f2;
      if(r>1){
        r=1;
      }
      
      double alpha = (double)rand()/(double)RAND_MAX;
      if(alpha<r){
        *(plista+i)=propuesta;
      }
          
      else{
        *(plista+i)=*(plista+i-1);
      }
      
	}


	double myarray [N];
	double * numeros;

	numeros=myarray;  //Puntero que almacena los numeros

	for(i=0;i<N;i++){
    rand_num=randn(mu,sigma);
		*(numeros+i)=rand_num;
	}	


  	// open file for writing 
	FILE *outfile; 
	outfile = fopen ("sample.dat", "w+"); 
	//fwrite(numeros , 1 , N , outfile );
	fwrite(myarray, 1, sizeof(myarray), outfile);
	//fwrite (numeros, sizeof(struct person), 1, outfile); 
	fclose(outfile);

/*
  	// open file for writing 
	FILE *outfile; 
  outfile = fopen ("sample.dat", "w+"); 
	//fwrite(numeros , 1 , N , outfile );
  fwrite(lista, 1, sizeof(lista), outfile);
	//fwrite (numeros, sizeof(struct person), 1, outfile); 
  fclose(outfile);
*/

  return 0;
}
