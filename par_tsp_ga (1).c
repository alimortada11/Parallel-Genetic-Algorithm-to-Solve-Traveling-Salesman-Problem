#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include "mpi.h"
#include <omp.h>
#include<time.h>

#define POPULATION_SIZE 8
#define MAX_NB_GENERATIONS 100
int CHROMOSOME_LENGTH;


// code run only on wi10.tsp and on population of size 8 on 4 Nodes


int	**genesMatrix;
float **matrix; //distance matrix
float *x , *y ;
int *city_id;

typedef struct Chromosome {
	int *genes; // a sequence cities or path
	float fitness;
}Chromosome;

typedef struct Generation {
	Chromosome *population;
	int index_of_fittest; // useless if population is sorted
}Generation;

Chromosome *population;
float *fitness_array;
void swap( Chromosome *pop , int src , int dest){
	Chromosome chrom;
	chrom  = pop[src];
	pop[src] = pop[dest];
	pop[dest] = chrom;
}

int getRandomNumberHalf(){
	int seed = (unsigned)(time(NULL)+rand());
	srand(seed);
	return rand()%CHROMOSOME_LENGTH/2;
}

void selection(Chromosome *pop){
	int n = POPULATION_SIZE;
	int random=0.0;
	n = (40*n)/100;
	for(int i = 0 ; i < ((10*POPULATION_SIZE)/100) ;i++ ){
		random =(POPULATION_SIZE/2) + getRandomNumberHalf(); 
		swap(pop , (n+i) , random);
	}	
}
void sort(Chromosome *pop){

	for(int i = 0 ;i < POPULATION_SIZE ; i++){
		for(int j = i+1  ; j < POPULATION_SIZE   ; j++){
			if(pop[i].fitness > pop[j].fitness ){
				swap(pop , i , j);
			}
		}
	}
}

int getRandomNumber(){
	int seed = (unsigned)(time(NULL)+rand());
	srand(seed);
	return rand()%CHROMOSOME_LENGTH;
}


void print_matrix(){																// print MATRIX 
	for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){
		for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
			printf(" %f " , matrix[i][j]);
		}	
	printf("\n");
	}
}

float calculate_distance(float x, float y , float x1 , float y1){						// calculate distance 
	return sqrt(pow((x-x1) , 2) + pow((y-y1) , 2));
}

void fill_distances(){																// fill matrix data
	//float result = calculate_distance(x[0] , y[0] , x[1] , y[1]);
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){
		for(int j = 0 ;j < CHROMOSOME_LENGTH ;j++){
			matrix[i][j] =calculate_distance(x[i],y[i] , x[j],y[j]); 
		}
	}
}
void fillGeneMatrix(){
	for(int i = 0 ; i < POPULATION_SIZE ; i++){
		for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
			genesMatrix[i][j] = population[i].genes[j];
		}
	}
}
void printGeneMatrix(){
	printf("Matrix");
	for(int i = 0 ; i < POPULATION_SIZE ; i++){
		printf("\n");
		for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
			printf("%d , ",genesMatrix[i][j]);
		}
	}
	printf("\n");
}
// code run only on wi10.tsp and on population of size 8 on 4 Nodes


void parse_arguments(int argc, char** argv)										// INITIALIZE matrix cities coordinates
{
	char buffer[100],c;  // to escape first 4 lines
	char *filename , *str;
	int read=0;

	FILE *file;
	if(argc == 1){
		printf("Please pass the file name");
		exit(0);
	}

	filename = argv[1];
	
	file = fopen( filename, "r");		
	int X =0 ;
	for(int i = 0 ; i < 7 ;i++){
		if(i == 4){
			fscanf(file,"%s %s %d",buffer, buffer,&CHROMOSOME_LENGTH);
		}
		fgets(buffer, 100, file);
	}
	genesMatrix = (int **)malloc( POPULATION_SIZE* sizeof(int *));		 // attach space for matrix in heap 
        for (int i=0; i < CHROMOSOME_LENGTH; i++)
			genesMatrix[i] = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int )); 

	x = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));				// allocate cities and their coordinates in heap memory
	y= (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
	city_id = (int *)malloc(CHROMOSOME_LENGTH*sizeof(int));

	//printf("use fscanf %d" , CHROMOSOME_LENGTH);
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){   							//load cities and their coordinates
		fscanf(file,"%d %f %f",&city_id[i], &x[i],&y[i]);
	}
	

	matrix = (float **)malloc( CHROMOSOME_LENGTH* sizeof(float *));		 // attach space for matrix in heap 
       for (int i=0; i < CHROMOSOME_LENGTH; i++)
		matrix[i] = (float *)malloc( CHROMOSOME_LENGTH* sizeof(float ));  

	for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){							// Initialize matrix
		for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
			matrix[i][j] = 0.0;
		}	
	}

	
													
}
void calculate_fitness(Chromosome *chro){
	float fit = 0.0;
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		if(i != CHROMOSOME_LENGTH-1){
			fit =fit + matrix[chro->genes[i] -1][chro->genes[i+1] -1];
		}
		else{
			fit =fit + matrix[chro->genes[i] -1][chro->genes[0]-1];	
		}
	}
	chro -> fitness = fit;
}
void calculate_fitness_i(Chromosome *chro ,int i){
	float fit = 0.0;
	//printf("iiiiiiiiiiii : %d",i);
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		if(i != CHROMOSOME_LENGTH-1){
			fit =fit + matrix[chro->genes[i] -1][chro->genes[i+1] -1];
		}
		else{
			fit =fit + matrix[chro->genes[i] -1][chro->genes[0]-1];	
		}
	}
	fitness_array[i] = fit;
	chro -> fitness = fit;
}
// code run only on wi10.tsp and on population of size 8 on 4 Nodes


void getRandomChromosome(Chromosome *chro ){
	int array[CHROMOSOME_LENGTH] ;
	int temp = 0,nb = 0 ;
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++ ){
		array[i] = i+1;
			
	}
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		nb = getRandomNumber()%(CHROMOSOME_LENGTH - i);
		temp = array[nb] ;
		array[nb] = array[CHROMOSOME_LENGTH - i -1];
		array[CHROMOSOME_LENGTH - i -1 ] = temp;
		chro -> genes[i] = temp;
		
	}
	calculate_fitness(chro ); 
}

void initPopulation(Chromosome *population)
{	
	for (int i = 0; i < POPULATION_SIZE ; i++){
		population[i].genes = (int *)malloc(CHROMOSOME_LENGTH*sizeof(int));
		getRandomChromosome(&population[i]);
	}
}


void printChro(){
	printf("\nChromosomes\n");
	for(int i = 0 ; i < POPULATION_SIZE ; i++){
		for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
			printf(" %d , " , population[i].genes[j]);
		}
		printf("\n");
	}
}

int if_exist(Chromosome *chrom , int x){
	//printf("chrom : %d , %d %d \n" ,chrom.genes[0] , chrom.genes[1] , x );
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		if(x == chrom->genes[i]){
			return 1;
		}
	}
	return 0;
}

void create_Child(Chromosome father  , Chromosome mother , Chromosome *Chro){
	//Chro->genes[0]=1;	
	//printf("B0 %d ",Chro->genes[0]);
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){
		Chro->genes[i]=0;
	}
	int clock =0 ,father_i =0 , mother_i = 0 ; 
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){
		if(clock == 0 ){
		//printf("hello clock 0\n");
			clock=1;
			while(mother_i < CHROMOSOME_LENGTH){
				if(if_exist(Chro , mother.genes[mother_i]) == 1){
					//printf("\n mother node %d " ,mother.genes[mother_i] );
					
					mother_i++;	
				}
				else{
					Chro->genes[i] = mother.genes[mother_i];
					break;
				}
			}
		}
		i++;
		if(clock == 1 ){
		//printf("hello clock 1\n");
			clock =0;
			while(father_i < CHROMOSOME_LENGTH){
				if(if_exist(Chro , father.genes[father_i]) == 1){
					father_i++;	
					
				}
				else{
					Chro->genes[i] = father.genes[father_i];
					break;
				}
			}
		}
	}
}
void print_fitness(Chromosome *pop){
	printf("\n");
	for(int i = 0 ; i < POPULATION_SIZE ;i++){		
		printf("%f    " , pop[i].fitness);
	}
	printf("\n");
}
float percentage_of_difference(Chromosome chro1 , Chromosome chro2){
	float sum = 0;  
	for(int i = 0 ; i < CHROMOSOME_LENGTH ;i++){
		if(chro1.genes[i] != chro2.genes[i]){
			sum++;
		}
	}	
	//printf("%f" , (sum*100)/CHROMOSOME_LENGTH);
	return (sum*100)/CHROMOSOME_LENGTH;
}
void crossover(Chromosome *pop){
	int j=0 , nb=0; 
	for(int i = 0 ; i <( POPULATION_SIZE/2) -1 ; i++){
	
		do{
			nb= getRandomNumber()%(POPULATION_SIZE/2);
		}while(nb == i && percentage_of_difference(pop[i] , pop[nb]) < 70);
		//printf("(%d ,%d)" , i , nb);
		create_Child(pop[i] , pop[nb] , &pop[(POPULATION_SIZE/2) +i]);
	}
	create_Child(pop[(POPULATION_SIZE/2)-1] , pop[0] , &pop[POPULATION_SIZE-1]);
}
void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro){
	//printf("Here");
	int n=0 , i=0,z=1;
	n = getRandomNumber()%(CHROMOSOME_LENGTH );
	//printf("{%d , %d}", n , n+((CHROMOSOME_LENGTH*30)/100));

	for(i=0 ; i < CHROMOSOME_LENGTH ;i++){
		Chro->genes[i] =0 ;
	}

	for( i = n ; i < n+((CHROMOSOME_LENGTH*30)/100);i++){
		z=i%CHROMOSOME_LENGTH;
		Chro->genes[z]=p.genes[z];
	}
	int c=0;
	i=(z+1)%CHROMOSOME_LENGTH;
	while( i!=z){
		c = c%CHROMOSOME_LENGTH;
		//printf("\nin loop %d\n",i);
		if(if_exist(Chro , m.genes[c]) != 1){
			Chro->genes[i] = m.genes[c];
		}
		else{
			if(Chro->genes[i] == 0){
				while(if_exist(Chro , m.genes[c]) == 1){
					c++;
				}
				Chro->genes[i] = m.genes[c];
			}
		}
		c++;
		i++;
		i=i%CHROMOSOME_LENGTH;
	}
	//rintf("z %d",z);
}
void crossoverV2(Chromosome *pop){
	int j=0 , nb=0; 
	for(int i = 0 ; i <( POPULATION_SIZE/2) -1 ; i++){
		do{
			nb= getRandomNumber()%(POPULATION_SIZE/2);
		}while(nb == i && percentage_of_difference(pop[i] , pop[nb]) < 70);
		create_ChildV2(pop[i] , pop[nb] , &pop[(POPULATION_SIZE/2) +i]);
	}
	//create_ChildV2(pop[(POPULATION_SIZE/2)-1] , pop[0] , &pop[POPULATION_SIZE-1]);
}
void mutation(Chromosome *pop){
	int i,j,k; 	
	for(int z =0 ; z < 5 ; z++){
		i = getRandomNumber()%(CHROMOSOME_LENGTH );
		j = getRandomNumber()%(CHROMOSOME_LENGTH );
		k = getRandomNumber()%(POPULATION_SIZE -(20*POPULATION_SIZE/100));
		
		//printf("%d ,%d , %d" , i , j , (20*POPULATION_SIZE/100)+k);
		
		int temp = pop[(20*POPULATION_SIZE/100)+k].genes[j];
		pop[(20*POPULATION_SIZE/100)+k].genes[j] = pop[(20*POPULATION_SIZE/100)+k].genes[i];
		pop[(20*POPULATION_SIZE/100)+k].genes[i] = temp;
	}
}

void main(int argc , char **argv){
	int   numtasks, taskid, tag1=0, tag2=1,tag3=2,tag4=3,tag5=4,tag6=5,tag7=7,tag8=8,tag9=9,tag10=10,tag11 =11,tag12 =12,offset=2 , numberOfChromosomesForEachNode ,extra;
	MPI_Status status;
	float *lastDistanceLine;
	float *BeforeLastLineDistance;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	//printf("\n task id  : %d \n" , taskid);
	numberOfChromosomesForEachNode = POPULATION_SIZE/numtasks;
	extra = POPULATION_SIZE%numtasks;
	fitness_array = (float *)malloc(POPULATION_SIZE*sizeof(float));
	for(int i = 0; i < POPULATION_SIZE ;i++){
		fitness_array[i] = 0.0;
	}
	int chromosomesNumberToBeSend;
	int tid = 0,nthreads = 0;
	if(taskid ==0 ){
	
		parse_arguments(argc , argv);
		fill_distances();
		population = (Chromosome *)malloc(POPULATION_SIZE*sizeof(Chromosome));
		initPopulation(population);
		sort(population);
		fillGeneMatrix();
		printChro();
		print_fitness(population);
		printf("\n");
		lastDistanceLine = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
		BeforeLastLineDistance = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
		for(int i =0 ;i < CHROMOSOME_LENGTH ;i++){
			lastDistanceLine[i] = matrix[CHROMOSOME_LENGTH-1][i];
			BeforeLastLineDistance[i] = matrix[CHROMOSOME_LENGTH-2][i];
		}
		for(int i = 1 ; i < numtasks ;i++  ){
			MPI_Send(&CHROMOSOME_LENGTH, 1, MPI_INT, i, tag6, MPI_COMM_WORLD);
		}
		int k =0,check=0;
		
		while(k < MAX_NB_GENERATIONS){
			k++;
			selection(population);
			crossoverV2(population);
			mutation(population);
			fillGeneMatrix();
			offset = 0;
			for(int i =0 ;i < numberOfChromosomesForEachNode ;i++){
				calculate_fitness_i(&population[i] ,i  );
			}
			
			for(int i=1 ; i < numtasks ;i++){
				offset = i*numberOfChromosomesForEachNode;
				if(i  <= extra)
					chromosomesNumberToBeSend = numberOfChromosomesForEachNode+1;
				else
					chromosomesNumberToBeSend = numberOfChromosomesForEachNode;
				
				// FItness Array contians a sorded Fitness value of Chromosomes in population 
				MPI_Send(&offset, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
				MPI_Send(&chromosomesNumberToBeSend, 1, MPI_INT, i, tag2, MPI_COMM_WORLD);
				MPI_Send(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, i, tag3, MPI_COMM_WORLD);				
				MPI_Send(&genesMatrix[offset][0], chromosomesNumberToBeSend*CHROMOSOME_LENGTH, MPI_INT, i, tag4,MPI_COMM_WORLD);
				MPI_Send(&genesMatrix[offset+1][CHROMOSOME_LENGTH-1], 1, MPI_INT, i, tag7, MPI_COMM_WORLD);
				MPI_Send(&genesMatrix[offset+1][CHROMOSOME_LENGTH-2], 1, MPI_INT, i, tag8, MPI_COMM_WORLD);
				MPI_Send(&matrix[0][0], CHROMOSOME_LENGTH*CHROMOSOME_LENGTH, MPI_FLOAT, i, tag9,MPI_COMM_WORLD);
				MPI_Send(&BeforeLastLineDistance[0], CHROMOSOME_LENGTH, MPI_FLOAT, i, tag11, MPI_COMM_WORLD);
				MPI_Send(&matrix[CHROMOSOME_LENGTH-1][0], CHROMOSOME_LENGTH, MPI_FLOAT, i, tag10, MPI_COMM_WORLD);	
					

				//printf("\nSender Offset : %d\n " ,  offset);
			}
			for(int i=1 ; i < numtasks ;i++){
				MPI_Recv(&offset, 1, MPI_INT, i, tag1, MPI_COMM_WORLD, &status);
				MPI_Recv(&chromosomesNumberToBeSend, 1, MPI_INT, i, tag2, MPI_COMM_WORLD, &status);
				MPI_Recv(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, i, tag3, MPI_COMM_WORLD, &status);
			}
			for(int i = 0 ; i < POPULATION_SIZE ; i++){
				population[i].fitness = fitness_array[i];
			}
			sort(population);
			
			printf("\nnew population : ");
			printf("Generation number : %d\n" , k);
			printChro();
			printf("\n -----------------------------------FITNESS OF CHROMOSOMES----------------------------------\n");
			print_fitness(population);
			printf("\n ------------------------------------------------------------------------------------------- \n");	
			printf("\n\n\n\n");
			if(k == MAX_NB_GENERATIONS){
				check=1;
				for(int i =1 ;i < numtasks ;i++){
					MPI_Send(&check, 1, MPI_INT, i, tag12, MPI_COMM_WORLD);
				}
				break;
			}
			else{
				check=0;
				for(int i =1;i < numtasks ;i++){
					MPI_Send(&check, 1, MPI_INT, i, tag12, MPI_COMM_WORLD);
				}
				
			}
			
		}
	}
	else{
		int check  = 0;
		MPI_Recv(&CHROMOSOME_LENGTH, 1, MPI_INT, 0, tag6, MPI_COMM_WORLD, &status);
		while(check == 0 ){
			
			genesMatrix = (int **)malloc( POPULATION_SIZE* sizeof(int *));		 // attach space for matrix in heap 
		        for (int i=0; i < CHROMOSOME_LENGTH; i++)
					genesMatrix[i] = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int )); 
			for(int i = 0 ; i < POPULATION_SIZE ; i++){
				for(int j = 0 ; j  < CHROMOSOME_LENGTH ;j++ ){
					genesMatrix[i][j] = 0;
				}
			}

			matrix = (float **)malloc( CHROMOSOME_LENGTH* sizeof(float *));		 // attach space for matrix in heap 
	       	for (int i=0; i < CHROMOSOME_LENGTH; i++)
				matrix[i] = (float *)malloc( CHROMOSOME_LENGTH* sizeof(float ));  

			for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){							// Initialize matrix
				for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
					matrix[i][j] = 0.0;
				}	
			}
			int last, prelast;
			MPI_Recv(&offset, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);
			MPI_Recv(&chromosomesNumberToBeSend, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD, &status);
			MPI_Recv(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, 0 , tag3, MPI_COMM_WORLD, &status);
			MPI_Recv(&genesMatrix[offset][0], chromosomesNumberToBeSend*CHROMOSOME_LENGTH, MPI_INT , 0 , tag4, MPI_COMM_WORLD, &status);
			MPI_Recv(&last, 1, MPI_INT, 0, tag7, MPI_COMM_WORLD, &status);
			MPI_Recv(&prelast, 1, MPI_INT, 0, tag8, MPI_COMM_WORLD, &status);
			MPI_Recv(&matrix[0][0], CHROMOSOME_LENGTH*CHROMOSOME_LENGTH, MPI_FLOAT , 0 , tag9, MPI_COMM_WORLD, &status);
			MPI_Recv(&BeforeLastLineDistance[0], CHROMOSOME_LENGTH, MPI_FLOAT, 0 , tag11, MPI_COMM_WORLD, &status);
			MPI_Recv(&matrix[CHROMOSOME_LENGTH-1][0], CHROMOSOME_LENGTH, MPI_FLOAT, 0 , tag10, MPI_COMM_WORLD, &status);

			for(int i =0 ; i < CHROMOSOME_LENGTH ;i++){
				matrix[CHROMOSOME_LENGTH-2][i] = BeforeLastLineDistance[i];
			}


			Chromosome chrom; 
			chrom.genes = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int ));
			for(int i =0 ; i < CHROMOSOME_LENGTH ; i++ ){
				chrom.genes[i] = 0;
			}
			genesMatrix[offset+1][CHROMOSOME_LENGTH-1] = last;
			genesMatrix[offset+1][CHROMOSOME_LENGTH-2] = prelast;
			for(int i = offset ; i < offset + chromosomesNumberToBeSend ; i++ ){
				//printf("\n matrix%d",taskid);
				for(int j =0 ;j < CHROMOSOME_LENGTH ;j++){
					chrom.genes[j] = genesMatrix[i][j];
					//printf(" %d -" , chrom.genes[j]);	
				}
				calculate_fitness_i(&chrom , i);
			}
			//printf("fitness array value %f %f \n",fitness_array[offset],fitness_array[offset+1]);

			MPI_Send(&offset, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD);
			MPI_Send(&chromosomesNumberToBeSend, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD);
			MPI_Send(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, 0, tag3, MPI_COMM_WORLD);
			
			MPI_Recv(&check, 1, MPI_INT, 0, tag12, MPI_COMM_WORLD, &status);
			if(check == 1){
				break;
			}
		}
	}
	MPI_Finalize();
	
}

//THE ABOVE RUN WITH MPI ONLY IT RUNS ONLY ON WI10.TXT on (CHROMOSOMES OF SIZE 10 only 10 cities ) 


//CODE WITH OPENMP IS IN COMMENT BUT the below code needs repair so i will upload it as soon as possible 
/*
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include "mpi.h"
#include <omp.h>
#include<time.h>

#define POPULATION_SIZE 8
#define MAX_NB_GENERATIONS 100
int CHROMOSOME_LENGTH;


int	**genesMatrix;
float **matrix; //distance matrix
float *x , *y ;
int *city_id;

typedef struct Chromosome {
	int *genes; // a sequence cities or path
	float fitness;
}Chromosome;

typedef struct Generation {
	Chromosome *population;
	int index_of_fittest; // useless if population is sorted
}Generation;

Chromosome *population;
float *fitness_array;
void swap( Chromosome *pop , int src , int dest){
	Chromosome chrom;
	chrom  = pop[src];
	pop[src] = pop[dest];
	pop[dest] = chrom;
}

int getRandomNumberHalf(){
	int seed = (unsigned)(time(NULL)+rand());
	srand(seed);
	return rand()%CHROMOSOME_LENGTH/2;
}

void selection(Chromosome *pop){
	int n = POPULATION_SIZE;
	int random=0.0;
	n = (40*n)/100;
	#pragma omp parallel shared(pop)
	{
		#pragma omp for schedule( dynamic, 2)
		for(int i = 0 ; i < ((10*POPULATION_SIZE)/100) ;i++ ){
			random =(POPULATION_SIZE/2) + getRandomNumberHalf(); 
			swap(pop , (n+i) , random);
		}
	}
	
}


void sort(Chromosome *pop){
	#pragma omp parallel shared(pop)
	{
		#pragma omp for schedule( dynamic, 2)
		for(int i = 0 ;i < POPULATION_SIZE ; i++){
			for(int j = i+1  ; j < POPULATION_SIZE   ; j++){
				if(pop[i].fitness > pop[j].fitness ){
					swap(pop , i , j);
				}
			}
		}
	}
}

int getRandomNumber(){
	int seed = (unsigned)(time(NULL)+rand());
	srand(seed);
	return rand()%CHROMOSOME_LENGTH;
}


void print_matrix(){																// print MATRIX 
	for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){
		for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
			printf(" %f " , matrix[i][j]);
		}	
	printf("\n");
	}
}

float calculate_distance(float x, float y , float x1 , float y1){						// calculate distance 
	return sqrt(pow((x-x1) , 2) + pow((y-y1) , 2));
}

void fill_distances(){																// fill matrix data
	//float result = calculate_distance(x[0] , y[0] , x[1] , y[1]);
	#pragma omp parallel 
	{
		#pragma omp for schedule( dynamic, 2)	
		for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){
			for(int j = 0 ;j < CHROMOSOME_LENGTH ;j++){
				matrix[i][j] =calculate_distance(x[i],y[i] , x[j],y[j]); 
			}
		}
	}
}
void fillGeneMatrix(){
	#pragma omp parallel 
	{
		#pragma omp for schedule( dynamic, 2)	
		for(int i = 0 ; i < POPULATION_SIZE ; i++){
			for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
				genesMatrix[i][j] = population[i].genes[j];
			}
		}
	}
}
void printGeneMatrix(){
	printf("Matrix");
	for(int i = 0 ; i < POPULATION_SIZE ; i++){
		printf("\n");
		for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
			printf("%d , ",genesMatrix[i][j]);
		}
	}
	printf("\n");
}


void parse_arguments(int argc, char** argv)										// INITIALIZE matrix cities coordinates
{
	char buffer[100],c;  // to escape first 4 lines
	char *filename , *str;
	int read=0;

	FILE *file;
	if(argc == 1){
		printf("Please pass the file name");
		exit(0);
	}

	filename = argv[1];
	
	file = fopen( filename, "r");		
	int X =0 ;
	for(int i = 0 ; i < 7 ;i++){
		if(i == 4){
			fscanf(file,"%s %s %d",buffer, buffer,&CHROMOSOME_LENGTH);
		}
		fgets(buffer, 100, file);
	}
	genesMatrix = (int **)malloc( POPULATION_SIZE* sizeof(int *));		 // attach space for matrix in heap 
        for (int i=0; i < CHROMOSOME_LENGTH; i++)
			genesMatrix[i] = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int )); 

	x = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));				// allocate cities and their coordinates in heap memory
	y= (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
	city_id = (int *)malloc(CHROMOSOME_LENGTH*sizeof(int));

	//printf("use fscanf %d" , CHROMOSOME_LENGTH);
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++){   							//load cities and their coordinates
		fscanf(file,"%d %f %f",&city_id[i], &x[i],&y[i]);
	}
	

	matrix = (float **)malloc( CHROMOSOME_LENGTH* sizeof(float *));		 // attach space for matrix in heap 
       for (int i=0; i < CHROMOSOME_LENGTH; i++)
		matrix[i] = (float *)malloc( CHROMOSOME_LENGTH* sizeof(float ));  

	for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){							// Initialize matrix
		for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
			matrix[i][j] = 0.0;
		}	
	}

	
													
}
void calculate_fitness(Chromosome *chro){
	float fit = 0.0;
	#pragma omp parallel shared(chro,fit)
	{
		#pragma omp for reduction(+:fit)
		for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
			if(i != CHROMOSOME_LENGTH-1){
				fit =fit + matrix[chro->genes[i] -1][chro->genes[i+1] -1];
			}
			else{
				fit =fit + matrix[chro->genes[i] -1][chro->genes[0]-1];	
			}
		}
	}
	chro -> fitness = fit;
}
void calculate_fitness_i(Chromosome *chro ,int i){
	float fit = 0.0;
	//printf("iiiiiiiiiiii : %d",i);
	#pragma omp parallel shared(chro)
	{
		#pragma omp for reduction(+:fit)
		for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
			if(i != CHROMOSOME_LENGTH-1){
				fit =fit + matrix[chro->genes[i] -1][chro->genes[i+1] -1];
			}
			else{
				fit =fit + matrix[chro->genes[i] -1][chro->genes[0]-1];	
			}
		}
	}
	fitness_array[i] = fit;
	chro -> fitness = fit;
}

void getRandomChromosome(Chromosome *chro ){
	int array[CHROMOSOME_LENGTH] ;
	int temp = 0,nb = 0 ;
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++ ){
		array[i] = i+1;
			
	}
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		nb = getRandomNumber()%(CHROMOSOME_LENGTH - i);
		temp = array[nb] ;
		array[nb] = array[CHROMOSOME_LENGTH - i -1];
		array[CHROMOSOME_LENGTH - i -1 ] = temp;
		chro -> genes[i] = temp;
		
	}
	calculate_fitness(chro ); 
}

void initPopulation(Chromosome *population)
{	
	#pragma omp parallel shared(population)
	{
		#pragma omp for schedule(dynamic,2)
		for (int i = 0; i < POPULATION_SIZE ; i++){
			population[i].genes = (int *)malloc(CHROMOSOME_LENGTH*sizeof(int));
			getRandomChromosome(&population[i]);
		}
	}
}


void printChro(){
	printf("\nChromosomes\n");
	for(int i = 0 ; i < POPULATION_SIZE ; i++){
		for(int j = 0 ; j < CHROMOSOME_LENGTH ; j++){
			printf(" %d , " , population[i].genes[j]);
		}
		printf("\n");
	}
}

int if_exist(Chromosome *chrom , int x){
	//printf("chrom : %d , %d %d \n" ,chrom.genes[0] , chrom.genes[1] , x );

		for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
			if(x == chrom->genes[i]){
				return 1;
			}
		}
	
	return 0;
}

void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro){
	//printf("Here");
	int n=0 , i=0,z=1;
	n = getRandomNumber()%(CHROMOSOME_LENGTH );
	//printf("{%d , %d}", n , n+((CHROMOSOME_LENGTH*30)/100));
	#pragma omp parallel shared(Chro)
	{
		#pragma omp for schedule(dynamic,2)
		for(i=0 ; i < CHROMOSOME_LENGTH ;i++){
			Chro->genes[i] =0 ;
		}
	}
	#pragma omp parallel shared(Chro)
	{
		#pragma omp for schedule(dynamic,2)
		for( i = n ; i < n+((CHROMOSOME_LENGTH*30)/100);i++){
			z=i%CHROMOSOME_LENGTH;
			Chro->genes[z]=p.genes[z];
		}
	}
	int c=0;
	i=(z+1)%CHROMOSOME_LENGTH;
	while( i!=z){
		c = c%CHROMOSOME_LENGTH;
		//printf("\nin loop %d\n",i);
		if(if_exist(Chro , m.genes[c]) != 1){
			Chro->genes[i] = m.genes[c];
		}
		else{
			if(Chro->genes[i] == 0){
				while(if_exist(Chro , m.genes[c]) == 1){
					c++;
				}
				Chro->genes[i] = m.genes[c];
			}
		}
		c++;
		i++;
		i=i%CHROMOSOME_LENGTH;
	}
	//rintf("z %d",z);
}
void crossoverV2(Chromosome *pop){
	int j=0 , nb=0; 
	#pragma omp parallel shared(pop)
	{
		#pragma omp for schedule(dynamic,2)
		for(int i = 0 ; i <( POPULATION_SIZE/2) -1 ; i++){
			do{
				nb= getRandomNumber()%(POPULATION_SIZE/2);
			}while(nb == i && percentage_of_difference(pop[i] , pop[nb]) < 70);
			create_ChildV2(pop[i] , pop[nb] , &pop[(POPULATION_SIZE/2) +i]);
		}
	}	
}
void mutation(Chromosome *pop){
	int i,j,k; 	
	#pragma omp parallel shared(pop)
	{
		#pragma omp for schedule(dynamic,2)
		for(int z =0 ; z < 5 ; z++){
			i = getRandomNumber()%(CHROMOSOME_LENGTH );
			j = getRandomNumber()%(CHROMOSOME_LENGTH );
			k = getRandomNumber()%(POPULATION_SIZE -(20*POPULATION_SIZE/100));
			
			//printf("%d ,%d , %d" , i , j , (20*POPULATION_SIZE/100)+k);
			
			int temp = pop[(20*POPULATION_SIZE/100)+k].genes[j];
			pop[(20*POPULATION_SIZE/100)+k].genes[j] = pop[(20*POPULATION_SIZE/100)+k].genes[i];
			pop[(20*POPULATION_SIZE/100)+k].genes[i] = temp;
		}
	}
}

void main(int argc , char **argv){


	int   numtasks, taskid, tag1=0, tag2=1,tag3=2,tag4=3,tag5=4,tag6=5,tag7=7,tag8=8,tag9=9,tag10=10,tag11 =11,tag12 =12,offset=2 , numberOfChromosomesForEachNode ,extra;
	MPI_Status status;
	float *lastDistanceLine;
	float *BeforeLastLineDistance;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	//printf("\n task id  : %d \n" , taskid);
	numberOfChromosomesForEachNode = POPULATION_SIZE/numtasks;
	extra = POPULATION_SIZE%numtasks;
	fitness_array = (float *)malloc(POPULATION_SIZE*sizeof(float));
	for(int i = 0; i < POPULATION_SIZE ;i++){
		fitness_array[i] = 0.0;
	}
	int chromosomesNumberToBeSend;
	int tid = 0,nthreads = 0;
	if(taskid ==0 ){
	
		parse_arguments(argc , argv);
		fill_distances();
		population = (Chromosome *)malloc(POPULATION_SIZE*sizeof(Chromosome));
		initPopulation(population);
		sort(population);
		fillGeneMatrix();
		printChro();
		print_fitness(population);
		printf("\n");
		lastDistanceLine = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
		BeforeLastLineDistance = (float *)malloc(CHROMOSOME_LENGTH*sizeof(float));
		for(int i =0 ;i < CHROMOSOME_LENGTH ;i++){
			lastDistanceLine[i] = matrix[CHROMOSOME_LENGTH-1][i];
			BeforeLastLineDistance[i] = matrix[CHROMOSOME_LENGTH-2][i];
		}
		for(int i = 1 ; i < numtasks ;i++  ){
			MPI_Send(&CHROMOSOME_LENGTH, 1, MPI_INT, i, tag6, MPI_COMM_WORLD);
		}
		int k =0,check=0;
		
		while(k < MAX_NB_GENERATIONS){
			k++;
			selection(population);
			crossoverV2(population);
			mutation(population);
			fillGeneMatrix();
			offset = 0;
			for(int i =0 ;i < numberOfChromosomesForEachNode ;i++){
				calculate_fitness_i(&population[i] ,i  );
			}
			
			for(int i=1 ; i < numtasks ;i++){
				offset = i*numberOfChromosomesForEachNode;
				if(i  <= extra)
					chromosomesNumberToBeSend = numberOfChromosomesForEachNode+1;
				else
					chromosomesNumberToBeSend = numberOfChromosomesForEachNode;
				
				MPI_Send(&offset, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
				MPI_Send(&chromosomesNumberToBeSend, 1, MPI_INT, i, tag2, MPI_COMM_WORLD);
				MPI_Send(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, i, tag3, MPI_COMM_WORLD);				
				MPI_Send(&genesMatrix[offset][0], chromosomesNumberToBeSend*CHROMOSOME_LENGTH, MPI_INT, i, tag4,MPI_COMM_WORLD);
				MPI_Send(&genesMatrix[offset+1][CHROMOSOME_LENGTH-1], 1, MPI_INT, i, tag7, MPI_COMM_WORLD);
				MPI_Send(&genesMatrix[offset+1][CHROMOSOME_LENGTH-2], 1, MPI_INT, i, tag8, MPI_COMM_WORLD);
				MPI_Send(&matrix[0][0], CHROMOSOME_LENGTH*CHROMOSOME_LENGTH, MPI_FLOAT, i, tag9,MPI_COMM_WORLD);
				MPI_Send(&BeforeLastLineDistance[0], CHROMOSOME_LENGTH, MPI_FLOAT, i, tag11, MPI_COMM_WORLD);
				MPI_Send(&matrix[CHROMOSOME_LENGTH-1][0], CHROMOSOME_LENGTH, MPI_FLOAT, i, tag10, MPI_COMM_WORLD);	
					

				//printf("\nSender Offset : %d\n " ,  offset);
			}
			for(int i=1 ; i < numtasks ;i++){
				MPI_Recv(&offset, 1, MPI_INT, i, tag1, MPI_COMM_WORLD, &status);
				MPI_Recv(&chromosomesNumberToBeSend, 1, MPI_INT, i, tag2, MPI_COMM_WORLD, &status);
				MPI_Recv(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, i, tag3, MPI_COMM_WORLD, &status);
			}
			for(int i = 0 ; i < POPULATION_SIZE ; i++){
				population[i].fitness = fitness_array[i];
			}
			sort(population);
			
			printf("\nnew population : ");
			printf("Generation number : %d\n" , k);
			printChro();
			printf("\n -----------------------------------FITNESS OF CHROMOSOMES----------------------------------\n");
			print_fitness(population);
			printf("\n ------------------------------------------------------------------------------------------- \n");	
			printf("\n\n\n\n");
			if(k == MAX_NB_GENERATIONS){        // if counter k == generation the parent send check =1 and the child when read check =0 he exit from loop
				check=1;
				for(int i =1 ;i < numtasks ;i++){
					MPI_Send(&check, 1, MPI_INT, i, tag12, MPI_COMM_WORLD);
				}
				break;
			}
			else{
				check=0;
				for(int i =1;i < numtasks ;i++){
					MPI_Send(&check, 1, MPI_INT, i, tag12, MPI_COMM_WORLD);
				}
				
			}
			
		}
	}
	else{
		int check  = 0;
		MPI_Recv(&CHROMOSOME_LENGTH, 1, MPI_INT, 0, tag6, MPI_COMM_WORLD, &status);
		while(check == 0 ){
			
			genesMatrix = (int **)malloc( POPULATION_SIZE* sizeof(int *));		 // attach space for matrix in heap 
		        for (int i=0; i < CHROMOSOME_LENGTH; i++)
					genesMatrix[i] = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int )); 
			for(int i = 0 ; i < POPULATION_SIZE ; i++){
				for(int j = 0 ; j  < CHROMOSOME_LENGTH ;j++ ){
					genesMatrix[i][j] = 0;
				}
			}

			matrix = (float **)malloc( CHROMOSOME_LENGTH* sizeof(float *));		 // attach space for matrix in heap 
	       	for (int i=0; i < CHROMOSOME_LENGTH; i++)
				matrix[i] = (float *)malloc( CHROMOSOME_LENGTH* sizeof(float ));  

			for(int i =0  ; i < CHROMOSOME_LENGTH ; i++){							// Initialize matrix
				for(int j = 0 ; j < CHROMOSOME_LENGTH ;j++ ){
					matrix[i][j] = 0.0;
				}	
			}
			int last, prelast;
			MPI_Recv(&offset, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status); // recieve offset 
			MPI_Recv(&chromosomesNumberToBeSend, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD, &status); // recieve number of rows 
			MPI_Recv(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, 0 , tag3, MPI_COMM_WORLD, &status);  //recieve array to fill the fitness by childs
			MPI_Recv(&genesMatrix[offset][0], chromosomesNumberToBeSend*CHROMOSOME_LENGTH, MPI_INT , 0 , tag4, MPI_COMM_WORLD, &status);// send matrix of genes
			MPI_Recv(&last, 1, MPI_INT, 0, tag7, MPI_COMM_WORLD, &status); // send last cell of the second row in matrix
			MPI_Recv(&prelast, 1, MPI_INT, 0, tag8, MPI_COMM_WORLD, &status);
			MPI_Recv(&matrix[0][0], CHROMOSOME_LENGTH*CHROMOSOME_LENGTH, MPI_FLOAT , 0 , tag9, MPI_COMM_WORLD, &status); // recieve matrix of distances 
			MPI_Recv(&BeforeLastLineDistance[0], CHROMOSOME_LENGTH, MPI_FLOAT, 0 , tag11, MPI_COMM_WORLD, &status);
			MPI_Recv(&matrix[CHROMOSOME_LENGTH-1][0], CHROMOSOME_LENGTH, MPI_FLOAT, 0 , tag10, MPI_COMM_WORLD, &status);

			for(int i =0 ; i < CHROMOSOME_LENGTH ;i++){
				matrix[CHROMOSOME_LENGTH-2][i] = BeforeLastLineDistance[i];
			}


			Chromosome chrom; 
			chrom.genes = (int *)malloc( CHROMOSOME_LENGTH* sizeof(int ));
			for(int i =0 ; i < CHROMOSOME_LENGTH ; i++ ){
				chrom.genes[i] = 0;
			}
			genesMatrix[offset+1][CHROMOSOME_LENGTH-1] = last;
			genesMatrix[offset+1][CHROMOSOME_LENGTH-2] = prelast;
			for(int i = offset ; i < offset + chromosomesNumberToBeSend ; i++ ){
				//printf("\n matrix%d",taskid);
				for(int j =0 ;j < CHROMOSOME_LENGTH ;j++){
					chrom.genes[j] = genesMatrix[i][j];
					//printf(" %d -" , chrom.genes[j]);	
				}
				calculate_fitness_i(&chrom , i);
			}
			//printf("fitness array value %f %f \n",fitness_array[offset],fitness_array[offset+1]);

			MPI_Send(&offset, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD);
			MPI_Send(&chromosomesNumberToBeSend, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD);
			MPI_Send(&fitness_array[offset], chromosomesNumberToBeSend, MPI_FLOAT, 0, tag3, MPI_COMM_WORLD);
			
			MPI_Recv(&check, 1, MPI_INT, 0, tag12, MPI_COMM_WORLD, &status);
			if(check == 1){
				break;
			}
		}
	}
	MPI_Finalize();
}*/
