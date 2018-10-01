#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////Estructuras/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Ranq1 {	//Estructura generadora de random de NR 3rd edition.
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

struct system{	//Estructura/clase con toda la informacion del sistema

	system(int,double);

	int **	AllocoMatriz			 (int L);
	void    Calculo_PESO_DE_BOLTZMAN (double T);
	int     Selecciona_IndexW 		 (int  DeltaE);
	void 	Acepto_Nuevo_Spin        (int i, int j);
	double 	Calculo_deltaE			 (int i, int j);
	void	update					 (int plot, int midvalues);
	void	Plot_XY					 ();
	void	FreeSystem				 ();

	int L;
	double   kb,J,W[2],E,M,T,Emid,E2mid,Mmid,M2mid;
	int **s, *plus, *minus;
	Ranq1 randf;						//Variable del tipo que devuelve el random
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////Funciones de la practica///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

system::system(int L0, double T0): randf(777){

		L=L0;    J=1;     kb=1;		T=T0;

		plus = (int *) calloc(L,sizeof(int));		//Vectores para las condiciones de contorno periodicas
		minus = (int *) calloc(L,sizeof(int));
		
		for(int i=0;i<L-1;i++){
			plus[i]=i+1;
			minus[i+1]=i;
		}

		plus[L-1]=0;
		minus[0]=L-1;

		s=AllocoMatriz(L);							//s es la matriz que me guarda los valores de los espines
		if(s==NULL) printf("Matrix NULL!!!\n");
		
		for(int i=0;i<L;i++)
			for(int j=0;j<L;j++)
				s[i][j]=1;							//Condicion inicial s[i][j]=1

		Calculo_PESO_DE_BOLTZMAN(T0);	

		Emid = E = -0.5*J*4*pow(L,2);	//Energia del estado ferromagnetico
		Mmid = M = 1*pow(L,2);			//Magnetizacion incial (= Magnetizacion de saturacion)
		E2mid = pow(Emid,2);
		M2mid = pow(Mmid,2);
}

int **system::AllocoMatriz(int L){

	int **s = (int **) calloc(L,sizeof(int *));

	for(int i=0;i<L;i++)
		s[i] = (int *) calloc(L,sizeof(int));

	return s;
}

void system::Calculo_PESO_DE_BOLTZMAN(double T){		//Funcion que llena el vector con el peso de boltzmann. Se calcula una vez porque los resultados de DeltaE son finitos (y pocos)
	
	double  DeltaE =0;
	DeltaE =8*J;               W[0]= exp(-DeltaE/(kb*T));
	DeltaE =4*J;               W[1]= exp(-DeltaE/(kb*T));
}

int system::Selecciona_IndexW(int DeltaE){				//Segun el valor de DeltaE devuelve el indice que necesita el vector de probabilidad
	if(DeltaE ==8)    return  0;
	if(DeltaE ==4)    return  1;
}

void system::Acepto_Nuevo_Spin(int i, int j){ s[i][j]*=-1;}		//Si acepto nuevo espin lo doy vuelta

double system::Calculo_deltaE (int i, int j){ return (2*J*s[i][j]*(s[plus[i]][j] + s[minus[i]][j] + s[i][plus[j]] + s[i][minus[j]]));}	//Calculo de DeltaE visto en la teoria

void system::update(int plot, int midvalues){			//Funcion que actualiza el sistema. Cada vez que la llamo se cumple un paso de montecarlo

	double DeltaE=0;

	for(int i=0;i<L;i++){								//Recorro SECUENCIALMENTE la red. Esto no es puramente random pero agiliza la ejecucion un monton. Seleccionar el indice al azar no tiene mucha incidencia en el resultado de la simulacion
		for(int j=0;j<L;j++){

			DeltaE = Calculo_deltaE(i,j);
						
			if(DeltaE<=0 || randf.doub()<W[Selecciona_IndexW(DeltaE)]) {	//Si DeltaE negativo o el numero random es manor que el vector de boltzmann acepto el cambio.
				Acepto_Nuevo_Spin(i,j);	
				E+=DeltaE;					//Si hago cambios actualizo los valores de E y M
				M+=2*s[i][j];				//Esta magnetizacion no esta normalizada. Lo hago cuando grafico.
			}
			
			if(plot && j==0) Plot_XY();		//Funcion que plotea la red en el estado actual
			
			if(midvalues){					//Calculo los valores medios solo si se lo pido
				Emid+=E;
				E2mid+=pow(E,2);
				Mmid+=fabs(M);
				M2mid+=pow(M,2);
			}
		}
	}
}


void system::FreeSystem(){					//Funcion que libera la memoria de la red.
	for(int i=0;i<L;i++) free (*(s+i));
	free(s);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para graficar
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

FILE *fid;
void Plot_ini(int L){

    fid = popen( "gnuplot -persist", "w");

        fprintf(fid, "set title 'Le Ferromagnetique'\n");
        fprintf(fid, "set xrange [-.1*%d:1.1*%d-1]\n",L,L);
        fprintf(fid, "set yrange [-.1*%d:1.1*%d-1]\n",L,L);
        fprintf(fid, "set palette maxcolors 2\n");
        fprintf(fid, "set palette defined(-1 'red', 1 'blue')\n");
        fprintf(fid, "unset colorbox\n");
 //       fprintf(fid, "set terminal gif animate delay 1 size 600,600\n");
 //       fprintf(fid, "set output 'GIORDANO_prac3_Ej_1a'\n");
}

void system::Plot_XY(){
      
    fprintf(fid,"plot '-' u 1:2:3 not ps %lf pt 5 lt palette\n",1.0*60/L);

    for(int i=0;i<L;i++)
    	for(int j=0;j<L;j++)
        	fprintf(fid,"%d\t%d\t%d\n",j,i,s[i][j]);

    fprintf(fid,"e\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////Resolucion////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Ej_1a(){

	size_t MCS=10000;
	int L0; double T; int plot=0; 
	char buf[50];

	printf("T Eja: ");
	scanf("%lf",&T);
	printf("L Eja: ");
	scanf("%d",&L0);
	printf("Plot?: ");
	scanf("%d",&plot);

	sprintf(buf,"datos/1a/1aT%.1lfL%d.dat",T,L0);
	FILE *f= fopen(buf,"w");

	struct system a(L0,T);							//Construyo el sistema

	if(plot) Plot_ini(a.L);

	for(int i=0;i<MCS;i++){						//Itero MCS veces actualizando el sistema
		fprintf(f, "%d\t%lf\t%lf\n",i,a.E,a.M);
		a.update(plot,0);
	}

	fclose(f);
	a.FreeSystem();
	return;
}

void Ej_1b(){

	size_t MCS=10000;
	int L0; 
	double T,Cmid=0,Xmid=0; 
	char buf[50];

	printf("L Ejb: ");
	scanf("%d",&L0);

	sprintf(buf,"datos/1b/1bL%d.dat",L0);
	FILE *f= fopen(buf,"w");
	double Ntry=MCS*L0*L0;

	for(size_t j=0;j<75;j++){
		struct system a(L0,1.5+.02*j);

		for(size_t i=0;i<MCS;i++)	
			a.update(0,1);			//Al updatear el 0 significa que no plotea y el 1 que calcula valores medios.
		
		a.Emid  /= Ntry;	a.E2mid /= Ntry;	a.Mmid  /= Ntry;	a.M2mid /= Ntry;		//Calculo los valores medios segun lo visto en la teoria
		Cmid  = (a.E2mid - pow(a.Emid,2))/(pow(a.T,2));
		Xmid  = (a.M2mid - pow(a.Mmid,2))/a.T;

		fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", a.T, a.Emid, a.Mmid, Cmid, Xmid);

		a.FreeSystem();
	}

	fclose(f);
	return;
}

void Ej_1c(int L0){

	size_t MCS=10000;
	int nits=215; 
	double T, Cmid=0, Xmid=0; 
	char buf[50];
	clock_t start, end;
    double cpu_time_used;

	sprintf(buf,"datos/1c/1cL%d.dat",L0);
	FILE *f= fopen(buf,"w");
	double Ntry=MCS*L0*L0;

	for(int j=0;j<nits;j++){	//3 casos : T<<Tc.  T~Tc.   T>>Tc 		//Tc=2.2692  		//Hago un barrido en temperatura con distinto step, siendo mas fino cerca de la transicion
	
		start = clock();
		
		if(j<=40 || j>160){
			if(j<=40)
				T=1.0+0.031*j;	//de 1 a 2,24

			if(j>160)
				T=2.3+(0.05)*(j-160); //de 2.35 a 5

			struct system a(L0,T);

			for(size_t i=0;i<MCS;i++)	
				a.update(0,1);		
	
			a.Emid  /= Ntry;	a.E2mid /= Ntry;	a.Mmid  /= Ntry;	a.M2mid /= Ntry;
			Cmid  = (a.E2mid - pow(a.Emid,2))/(pow(a.T,2));
			Xmid  = (a.M2mid - pow(a.Mmid,2))/a.T;

			fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", a.T, a.Emid, a.Mmid, Cmid, Xmid);
	
			a.FreeSystem();
		}

		if(j>40 && j<161){
			T = 2.24+.0005*(j-40);		//de 2.2405 a 2.3

			struct system a(L0,T);
		
			for(size_t i=0;i<MCS;i++)		//Itero en la zona cercana a la critica.
				a.update(0,1);			
	
			a.Emid  /= Ntry;	a.E2mid /= Ntry;	a.Mmid  /= Ntry;	a.M2mid /= Ntry;
			Cmid  = (a.E2mid - pow(a.Emid,2))/(pow(a.T,2));
			Xmid  = (a.M2mid - pow(a.Mmid,2))/a.T;

			fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", a.T, a.Emid, a.Mmid, Cmid, Xmid);

			a.FreeSystem();
		}

		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf("It: %d of %d, Temperature: %.3lf,  It time:  %.1lf s, Time left: %.1lf min.\n",j,nits-1,T,cpu_time_used,(nits-j)*cpu_time_used/60);
		printf("--------------------------------------------------------------------------------\n");
	}

	fclose(f);
	return;
}

void Ej_1d(){

	size_t MCS=10000;
	int L0, idx; 
	double T,mhist[101],tot=0; 
	char buf[50];
	size_t i;

	for(i=0;i<101;i++)		//Hago un vector que guarda el histograma de 101 bins entre -1 y 1 
		mhist[i]=0;

	printf("L Ejd: ");
	scanf("%d",&L0);
	printf("T Ejd: ");
	scanf("%lf",&T);

	struct system a(L0,T);

	sprintf(buf,"datos/1d/1dT%.1lfL%d.dat",T,L0);
	FILE *f= fopen(buf,"w");

	for(i=0;i<1000*MCS;i++){						//Podria tirar los primeros pasos hasta que entre en equilibrio pero no tienen mucha peso dentro de los 10^8 pasos de montecarlo
		a.update(0,0);
		idx = (int)(a.M*50/(L0*L0)+50);			//Calculo en que bin caeria el valor de M en cada update.
		mhist[idx]++;							//Sumo el bin donde cae idx.
	}

		tot+=1000*MCS;							//Esta es la constante de normalizacion de los bins/dx	No multiplico por dx porque grafico la densidad de probabilidad	

	for(i=0;i<101;i++)
		fprintf(f, "%lf\t%lf\n",(double)i/50 - 1, mhist[i]/tot);		//Printeo la densidad de probabilidad

	a.FreeSystem();
	fclose(f);
	return;
}

int main()
{
	int run,retry=1, L0;
	Ranq1(12934);		//Inicializo ranq con una seed cualquiera

	while(retry){
		printf("\n--------------------------------------------------------------------------------\n");
		printf("Gimme the exercise number (a=1 b=2 ...):  ");
		scanf("%d",&run);

		switch(run){
			case(1):
				Ej_1a();
				break;		
			
			case(2):
				Ej_1b();
				break;

			case(3):{
				printf("L0?: ");
				scanf("%d",&L0);
				Ej_1c(L0);
				break;
			}

			case(4):
				Ej_1d();
				break;
		}
		
		printf("Retry?: ");
		scanf("%d",&retry);
	}

	printf("\nTermine!\n");
   	return 0;
}