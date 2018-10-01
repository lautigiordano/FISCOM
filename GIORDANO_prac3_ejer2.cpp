#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////Estructuras/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Ranq1 {
    typedef unsigned long long int Ullong;
    typedef unsigned int Uint;
    typedef double Doub;

    Ullong v;
    Ranq1(Ullong j) : v(4101842887655102017LL) {
    v^=j;
    v=int64();
	}
    inline Ullong int64() {
        v^=v>>21;v^=v<<35;v^=v>>4;
        return v * 2685821657736338717LL;
    }
    inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
    inline Uint int32() { return (Uint)int64(); }
};

struct system{					//Solo voy a comentar las cosas diferentes al ejercicio 1. Si hay alguna duda mirar GIORDANO_prac3_ejer1.cpp

	system(int,double,double);

	int **	AllocoMatriz			 (int L);
	void    Calculo_PESO_DE_BOLTZMAN (double T);
	int 	Map 					 (double DeltaE);
	int 	Propongo_Nuevo_Spin      (int i, int j);
	double 	Calculo_deltaE			 (int i, int j, int nova);
	void	update					 (int midvalues);
	void	FreeSystem				 ();

	int L;
	double   kb,J,W[100],E,M,T,Delta,Emid,E2mid,Mmid,M2mid;    
	int **s, *plus, *minus;
	Ranq1 randf;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////Funciones de la practica///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

system::system(int L0, double T0, double delta): randf(777){

		L=L0;    J=1;     kb=1;		T=T0;	 Delta=delta*J;

		plus = (int *) calloc(L,sizeof(int));
		minus = (int *) calloc(L,sizeof(int));
		
		for(int i=0;i<L-1;i++){
			plus[i]=i+1;
			minus[i+1]=i;
		}

		plus[L-1]=0;
		minus[0]=L-1;

		s=AllocoMatriz(L);
		if(s==NULL) printf("Matrix NULL!!!\n");
		
		for(int i=0;i<L;i++)
			for(int j=0;j<L;j++)
				s[i][j]=1;

		Calculo_PESO_DE_BOLTZMAN(T0);

		Emid = E = -0.5*J*4*pow(L,2) + Delta*pow(L,2)*1;	
		Mmid = M = 1*pow(L,2);
		E2mid = pow(Emid,2);
		M2mid = pow(Mmid,2);
}

int ** system::AllocoMatriz(int L){

	int **s = (int **) calloc(L,sizeof(int *));

	for(int i=0;i<L;i++)
		s[i] = (int *) calloc(L,sizeof(int));

	return s;
}

int system::Map(double DeltaE){ //Funcion que devuelve indice para el peso de boltzmann. No resuelvo colisiones porque los indices son todos distintos (hice la prueba)
	double num = 5137*DeltaE;
 	return ((int)num%100);
}	

void system::Calculo_PESO_DE_BOLTZMAN(double T){
	double  DeltaE = 0;		//3 casos: (1) S!=0, S cambia de signo || (2) S=+-1 a S=0 || (3)  S=0 a S=+-1;

	for(size_t i=0;i<100;i++)	W[i]=0;

	for(int neighsum=-4;neighsum<5;neighsum++){	//Caso (1). neighsum es la suma de vecinos. Se contemplan los casos S=1 y S=-1 porque recorro neighsum de -4 a 4
		DeltaE= 2*J*(neighsum);
		if(DeltaE>0)
			W[Map(DeltaE)] = exp(-DeltaE/(kb*T));
	}

	for(int neighsum=-4;neighsum<5;neighsum++){	//Caso (2)
		DeltaE= J*(neighsum) - Delta;
		if(DeltaE>0)
			W[Map(DeltaE)] = exp(-DeltaE/(kb*T));
	}

	for(int neighsum=-4;neighsum<5;neighsum++){	//Caso (3)
		DeltaE= -J*(neighsum) + Delta;
		if(DeltaE>0)
			W[Map(DeltaE)] = exp(-DeltaE/(kb*T));
	}
}

int system::Propongo_Nuevo_Spin(int i, int j){		//Devuelve un numero random distinto al espin actual para proponer el cambio
	int nova, retries = 100;		//Esta variable es para que no proponga el mismo spin al anterior
	double rand;
	
	while(retries--){	//El retry es para ahorrarme que el nuevo indice sea igual al anterior. La probabilidad de que retry llegue a 0 es bajisima. 
		rand = 3*randf.doub();
		nova = (int)rand - 1;
		if(s[i][j]!=nova)	return nova;
	}

	printf("Error in Propongo_Nuevo_Spin\n");
	return 0;
}

double system::Calculo_deltaE (int i, int j, int nova){				//Calculo DeltaE con la nueva definicion de energia.
	int neighsum = s[minus[i]][j] + s[i][minus[j]] + s[i][plus[j]] + s[plus[i]][j];
	double DeltaE = -J*(nova-s[i][j])*neighsum + Delta*(pow(nova,2)-pow(s[i][j],2));
	return DeltaE;
}

void system::update(int midvalues){
	double DeltaE=0;
	int nova;

	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){

			nova = Propongo_Nuevo_Spin(i,j);		//Ahora en vez de dar vuelta el spin propongo uno y lo guardo en la variable nova. Es el potencial nuevo espin.
			DeltaE = Calculo_deltaE(i,j,nova);
						
			if(DeltaE<=0 || randf.doub()<W[Map(DeltaE)]) {
				E+=DeltaE;
				M+= (nova - s[i][j]);
				s[i][j] = nova;				//Hago el cambio.
			}
						
			if(midvalues){
				Emid+=E;
				E2mid+=pow(E,2);
				Mmid+=M;
				M2mid+=pow(M,2);
			}
		}
	}
}


void system::FreeSystem(){
	for(int i=0;i<L;i++)	free (*(s+i));
	free(s);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////Resolucion////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Ej_2a(){

	size_t MCS=10000;
	int L0; 
	double T; 
	char buf[50];

	printf("T Eja: ");
	scanf("%lf",&T);
	printf("L Eja: ");
	scanf("%d",&L0);

	sprintf(buf,"datos/2a/2aT%.2lfL%d.dat",T,L0);
	FILE *f= fopen(buf,"w");

	struct system a(L0,T,1.975);

	for(int i=0;i<MCS;i++){
		fprintf(f, "%d\t%lf\t%lf\n",i,a.E,a.M);
		a.update(0);
	}

	fclose(f);
	a.FreeSystem();
	return;
}

void Ej_2b(){

	size_t MCS=10000;
	int L0; 
	double T,Delta,Cmid=0,Xmid=0; 
	char buf[50];

	printf("L Ejb: ");
	scanf("%d",&L0);
	printf("Delta Ejb: ");
	scanf("%lf",&Delta);

	sprintf(buf,"datos/2b/2bL%dDelta%.3lf.dat",L0,Delta);
	FILE *f= fopen(buf,"w");
	
	double Ntry=MCS*L0*L0;

	for(size_t j=0;j<100;j++){
		struct system a(L0,.1+.05*j,Delta);

		for(size_t i=0;i<MCS;i++)
			a.update(1);

		a.Emid  /= Ntry;	a.E2mid /= Ntry;	a.Mmid  /= Ntry;	a.M2mid /= Ntry;
		Cmid  = (a.E2mid - pow(a.Emid,2))/(pow(a.T,2));
		Xmid  = (a.M2mid - pow(a.Mmid,2))/a.T;

		fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", a.T, a.Emid, a.Mmid, Cmid, Xmid);

		a.FreeSystem();
	}

	fclose(f);
	return;
}

void Ej_2c(int L0){

	size_t MCS=10000;
	int nits=321; 
	double T=.2, Cmid=0, Xmid=0, Delta; 
	char buf[50];
	clock_t start, end;
    double cpu_time_used;

	printf("Delta Ejc: ");
	scanf("%lf",&Delta);

	sprintf(buf,"datos/2c/2cL%dDelta%.3lf.dat",L0,Delta);
	FILE *f= fopen(buf,"w");
	double Ntry=MCS*L0*L0;

	struct system a(L0,T,Delta);						//Inicializo el sistema una sola vez para que el resto continue desde el esatdo anterior

	for(int j=0;j<nits;j++){							//Barrido en temperaturas ida y vuelta para ver la histeresis
		start = clock();

		a.Calculo_PESO_DE_BOLTZMAN(T);

		a.Emid = a.Mmid = a.E2mid = a.M2mid = 0;		//Reinicio los valores de EMid, Mmid, etc. del sistema con los actuales en cada iteracion
		
		if(j<=20){										// Voy desde 0.2 a 0.4 con step 0.01
			T=0.2+0.01*j;
			a.T = T;									//Aca cambio la temperatura del sistema en cada iteracion

			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		if(j>20 && j<=141){								// Voy desde 0.4 a 0.8 con step 0.003333
			T=0.4+0.00333333*(j-21);
			a.T = T;

			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		if(j>141 && j<=161){							// Voy desde 0.8 a 1
			T=0.8+0.01*(j-141);
			a.T = T;
			
			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		if(j>161 && j<=181){							// Vuelvo desde 1 a 0.8
			T=1-0.01*(j-161);
			a.T = T;

			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		if(j>181 && j<=301){							// Vuelvo desde 0.8 a 0.4
			T=0.8-0.003333*(j-181);
			a.T = T;

			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		if(j>301 && j<=321){							// Vuelvo desde 0.4 a 0.2
			T=0.4-0.01*(j-301);
			a.T = T;

			for(size_t i=0;i<MCS;i++)	
				a.update(1);			
		}

		a.Emid  /= Ntry;	a.E2mid /= Ntry;	a.Mmid  /= Ntry;	a.M2mid /= Ntry;
		Cmid  = (a.E2mid - pow(a.Emid,2))/(pow(a.T,2));
		Xmid  = (a.M2mid - pow(a.Mmid,2))/a.T;

		fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", a.T, a.Emid, a.Mmid, Cmid, Xmid);

		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf("It: %d of %d, Temperature: %.3lf,  It time:  %.1lf s, Time left: %.1lf min.\n",j,nits-1,T,cpu_time_used,(nits-j)*cpu_time_used/60);
		printf("--------------------------------------------------------------------------------\n");
	}

	a.FreeSystem();

	fclose(f);
	return;
}

void Ej_2d(){

	size_t MCS=10000,i;
	int L0=20, idx; 
	double T,Delta,Ehist[101]; 
	char buf[50];

	for(i=0;i<101;i++)	
		Ehist[i]=0;

	printf("T Ejd: ");
	scanf("%lf",&T);
	printf("Delta Ejc: ");
	scanf("%lf",&Delta);

	struct system a(L0,T,Delta);

	sprintf(buf,"datos/2d/2dT%.3lfL%dD%.3lf.dat",T,L0,Delta);
	FILE *f= fopen(buf,"w");

	for(i=0;i<100*MCS;i++){
		a.update(0);
//		fprintf(f,"%d\t%lf\n",i,a.E/(20*20));

		idx = (int)((100/1.3)*a.E/(L0*L0)+(100/1.3));	//Centro el histograma. La energia empieza en ~ -0.05. Emax/L**2 ~ 0.45 // 0.45*220 + 1 = 100
		Ehist[idx]++;
	}

	for(i=0;i<101;i++)
		fprintf(f, "%lf\t%lf\n",(double)i*(1.3/100) - 1, Ehist[i]/(100*MCS));		//Printeo la densidad de probabilidad

	a.FreeSystem();
	fclose(f);
	return;
}

int main()
{
	int run,retry=1, L0;
	Ranq1(43229);	

	while(retry){
		printf("--------------------------------------------------------------------------------\n");
		printf("Gimme the exercise number (a=1, b=2 ...):  ");
		scanf("%d",&run);

		switch(run){
			case(1):
				Ej_2a();
				break;		
			
			case(2):
				Ej_2b();
				break;

			case(3):
				printf("L0?: ");
				scanf("%d",&L0);
				Ej_2c(L0);
				break;

			case(4):
				Ej_2d();
				break;
		}
		printf("Retry?: ");
		scanf("%d",&retry);
	}
	printf("\nTermine!\n");
   	return 0;
}