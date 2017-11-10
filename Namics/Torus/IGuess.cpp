#include <math.h>
#include <stdio.h>

double u;
double X,Y,Z;
double pi=3.141593;
int countA=0,countB=0;
int n_layers;
double valueA=0,valueB=0;

char fname[150];
main() {
	printf("Number of lattice layers in the system: "); scanf("%d", &n_layers);
	printf("Enter file name: ");  scanf("%s", &fname);
	FILE *out_file = fopen(fname,"w");
	fprintf(out_file,"gradients \n3 \n%i \n%i \n%i \nmolecule \nall \nstate \nA \n", n_layers, n_layers, n_layers);
	for (int x=1; x<=n_layers; x++)
	for (int y=1; y<=n_layers; y++)
	for (int z=1; z<=n_layers; z++){
		X=1.0*x/(1.0*n_layers); Y=1.0*y/(1.0*n_layers); Z=1.0*z/(1.0*n_layers);
		if ( cos(pi*X) + cos(pi*Y) + cos(pi*Z) <0){u = 1.5; countA++;} else {u = 1.0; countB++;}
		fprintf(out_file,"%1f \n",u);
	}
	fprintf(out_file,"molecule \n all \n state \n B \n");
	for (int x=1; x<=n_layers; x++)
	for (int y=1; y<=n_layers; y++)
	for (int z=1; z<=n_layers; z++){
		X=1.0*x/(1.0*n_layers); Y=1.0*y/(1.0*n_layers); Z=1.0*z/(1.0*n_layers);
		if ( cos(pi*X) + cos(pi*Y) + cos(pi*Z) <0){u = 0.68; } else {u = 1.0; }
		fprintf(out_file,"%1f \n",u);
	}
	fprintf(out_file,"phibulk solvent \n0.8147400257 \nalphabulk \nA \n1\nalphabulk \nB \n1 \n ");
	fclose(out_file);
	//intf("count A = %i \n count B = %i \n", countA,countB);
	//intf("valueA = %1f \n valueB = %1f \n ",valueA/(countA+countB),valueB/(countA+countB));
	return(0);
};
