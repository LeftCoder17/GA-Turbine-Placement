#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#define L 10 //square side length
#define Rr  40 //Rotor radius
#define N 100 //number of possible turbines
#define z 60  // hub height
#define z0 0.3 //ground roughness
#define CT 0.88 //Coefficient thurst
#define POP 600  //number of individuals
#define NSUBPOP 20 // number of subpopulations
#define SUBPOP (POP / NSUBPOP) // number of individuals per subpopulation
#define e 2.718281828
#define PI 3.14159265
#define NGEN 1200 //number of generations
#define a 0.3267949 //axial induction factor
#define R1 27.88071  //downstream radius
#define alpha 0.0943695 //entrainment constant

// Define structures
typedef struct {
    char genes[N];
    double fitness;
} Individual;

typedef struct {
    Individual individuals[SUBPOP];
} Subpopulation;

typedef struct {
    Subpopulation subpopulations[NSUBPOP];
    Individual best_individual;
} Population;

// Define functions
double compute_fitness_a(char *);
double compute_fitness_b(char *);
double compute_fitness_c(char *, float *);
void sort_subpopulation(Subpopulation *);
void sort_population(Population *);
void migration(Population *, int );
void crossover(char *parent1, char *parent2, char *child1, char *child2);
double power(double *u);
void speeds_a(char * indi, double u0, double *u); //genes of the individual as input
void speeds_bc(char * indi, double u0, double *u, double theta); //genes of the individual as input
double cost(char * indi);


int main(int argv, char **argc) {
    // Wind variables
    char c;
    c = 'c'; // Specify a, b or c
    int nturbines=0;
    float wind_prob[108];
    FILE *wind_file = fopen("wind_prob.txt", "r");
    FILE * fitnessfile;
    if (wind_file == NULL) {
        printf("Error opening file\n");
        exit(1);
    }
    if (c == 'c')
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 36; j++)
                fscanf(wind_file, "%f", &wind_prob[i*36 + j]);
    fclose(wind_file);

    // GA parameters
    float mutation_probability = 0.05;
    int migration_interval = 50;
    int n_migrants = 2; // 10% of the subpopulation
    int n_parents = 5;
    srand(1);

    // GA variables
    int generation;
    Population population, new_population;
    double totalFitnesses, runningTotal;
    int randomSelectPoint;
    int non_existed;
    int new_2_ind, parent1, parent2;
    
    // Others
    int subpop, ind, gene;
    int i,j;

    // File to store the evolution of the fitness. Choose the case at hand
    fitnessfile = fopen("fitness_evolution_case_c_migration_best","w");


    // 1. Initiallize subpopulations
    for (subpop=0; subpop<NSUBPOP; subpop++)
        for (ind=0; ind < SUBPOP; ind++)
            for (gene=0; gene<N; gene++) {
                i = rand() % 2; // random bit value
                population.subpopulations[subpop].individuals[ind].genes[gene] = i + '0';
            }
    population.best_individual = population.subpopulations[0].individuals[0];
    if (c == 'a') {
        population.best_individual.fitness = compute_fitness_a(population.best_individual.genes);
    } else if (c == 'b') {
        population.best_individual.fitness = compute_fitness_b(population.best_individual.genes);
    } else {
        population.best_individual.fitness = compute_fitness_c(population.best_individual.genes, wind_prob);
    }

    // 2. Start Generations
    for (generation=0; generation<NGEN; generation++) {
        // 2.1. Compute the fitness of every individual
        for (subpop=0; subpop<NSUBPOP; subpop++){
            for (ind=0; ind < SUBPOP; ind++) {
                if (c == 'a') {
                    population.subpopulations[subpop].individuals[ind].fitness = compute_fitness_a(population.subpopulations[subpop].individuals[ind].genes);
                } else if (c == 'b') {
                    population.subpopulations[subpop].individuals[ind].fitness = compute_fitness_b(population.subpopulations[subpop].individuals[ind].genes);
                } else {
                    population.subpopulations[subpop].individuals[ind].fitness = compute_fitness_c(population.subpopulations[subpop].individuals[ind].genes, wind_prob);

                }
            }
        }

        // 2.2. Sort every subpopulation
        for (subpop=0; subpop<NSUBPOP; subpop++){
                sort_subpopulation(&population.subpopulations[subpop]);
        }

        // 2.3. Check if there is a better new best_individual
        for (subpop=0; subpop<NSUBPOP; subpop++)
            if (population.best_individual.fitness > population.subpopulations[subpop].individuals[0].fitness){
                population.best_individual = population.subpopulations[subpop].individuals[0];
            }
        if (generation == NGEN - 1) break;

        fprintf(fitnessfile,"%.8f\n",population.best_individual.fitness);

        // 2.4. Migration
        if ((generation % migration_interval) == migration_interval - 1) {
            migration(&population, n_migrants);
            for (subpop=0; subpop<NSUBPOP; subpop++){
                sort_subpopulation(&population.subpopulations[subpop]);
            }
        }

        // 2.5. Start breeding new generation
        for (subpop=0; subpop<NSUBPOP; subpop++) {
            totalFitnesses = 0.0;
            for (ind=0; ind<n_parents; ind++){
                totalFitnesses += 1.0 / population.subpopulations[subpop].individuals[ind].fitness;
                }
            for (new_2_ind=0; new_2_ind<SUBPOP/2; new_2_ind++) {
                // 2.5.1. Roulette to choose parents
                runningTotal=0.0;
                randomSelectPoint = rand() % (int)(totalFitnesses + 1);
                for (ind=0; ind<n_parents; ind++) {
                    runningTotal += 1.0 / population.subpopulations[subpop].individuals[ind].fitness;
                    if (runningTotal >= randomSelectPoint){
                        parent1 = ind;
                        break;
                    }
                }
                runningTotal=0.0;
                randomSelectPoint = rand() % (int)(totalFitnesses + 1);
                for (ind=0; ind<n_parents; ind++) {
                    runningTotal += 1.0 / population.subpopulations[subpop].individuals[ind].fitness;
                    if (runningTotal >= randomSelectPoint){
                        parent2 = ind;
                        break;
                    }
                }

                // 2.5.2. Crossover
                crossover(population.subpopulations[subpop].individuals[parent1].genes, population.subpopulations[subpop].individuals[parent2].genes, new_population.subpopulations[subpop].individuals[2*new_2_ind].genes, new_population.subpopulations[subpop].individuals[2*new_2_ind + 1].genes);

                // 2.5.3. Mutation
                for(i=0; i<2; i++) {
                    j = rand() % (int)(1/mutation_probability);
                    if (j==0) {
                        gene = rand() % N;
                        if (new_population.subpopulations[subpop].individuals[2*new_2_ind + i].genes[gene] == '0') new_population.subpopulations[subpop].individuals[2*new_2_ind + i].genes[gene] = '1';
                        else new_population.subpopulations[subpop].individuals[2*new_2_ind + i].genes[gene] = '0';
                    }
                }

                // 2.5.4. Check whether the children existed in the last generation or not
                for(i=0; i<2; i++) {
                    for (ind=0; ind<SUBPOP; ind++) {
                        non_existed = strcasecmp(new_population.subpopulations[subpop].individuals[2*new_2_ind + i].genes, population.subpopulations[subpop].individuals[ind].genes);
                        if (!non_existed) break;
                    }
                    if (!non_existed) {
                        new_2_ind--;
                        break;
                    }
                }
            }
        } // Breeding new generation loop

        // 2.6. Copy new individuals to individuals
        for (subpop=0; subpop<NSUBPOP; subpop++)
            for (ind=0; ind<SUBPOP; ind++) 
                strcpy(population.subpopulations[subpop].individuals[ind].genes, new_population.subpopulations[subpop].individuals[ind].genes);

        // 2.7. Print some results
        if (generation % 100 == 0) printf("Generation %d completed. Best fitness is %f\n", generation, population.best_individual.fitness);

    } // Finish generations loop

    // 3. Print final results
    printf("Generations Finished!\nBest individual is:");
    for (gene=0; gene<N; gene++) {
        if (gene%10 == 0) printf("\n");
        printf("%c", population.best_individual.genes[gene]);
    }
    printf("\nBest fitness is: %.8f\n", population.best_individual.fitness);

    for (gene=0; gene<N; gene++)
        if (population.best_individual.genes[gene] == '1')
            nturbines++;
    
    printf("Number of turbines: %d\n",nturbines);  

    fclose(fitnessfile);
    return 0;
}  // Finish main

// Auxiliar functions
double compute_fitness_a(char *genes) {
    double u[N]; // turbine speeds
    double u0 = 12.0; //case a

    memset(u, 0, N * sizeof(double));
    speeds_a(genes,u0,u);
    double fit = cost(genes) / power(u);
    return fit;
}

double compute_fitness_b(char *genes) {
    double u[N]; // turbine speeds
    double u0 = 12.0; // case b
    double theta;
    double inv_power_sum = 0.0;
    double inv_power;

    for (int i=0; i<36; i++) {
        theta=i*10.0*PI/180.0; //0, 10, 20, 30, ... and to rad
        memset(u, 0, N * sizeof(double));
        speeds_bc(genes,u0,u, theta);
        inv_power = 1.0 / power(u);
        inv_power_sum += inv_power;
    }
    double fit = cost(genes) * inv_power_sum / 36.0;
    return fit;
}

double compute_fitness_c(char *genes, float *wind_prob) {
    double u[N]; // turbine speeds
    double theta, weight;
    double inv_power_sum = 0.0;
    double fit;
    double v[3]; //wind speeds case c
    v[0]=8.0;
    v[1]=12.0;
    v[2]=17.0;

    for (int i=0; i<3; i++)
        for (int j=0; j<36; j++) {
            theta=j*10.0*PI/180.0; //0, 10, 20, 30, ... and to rad
            memset(u, 0, N * sizeof(double));
            speeds_bc(genes,v[i],u, theta);
            weight = wind_prob[i*36 + j];
            inv_power_sum += (1.0 / (power(u))) * weight;
    }
    fit = cost(genes) * inv_power_sum;
    return fit;
}

void sort_subpopulation(Subpopulation *subpopulation) {
    Individual aux_individual;
    int ind, ind2;
    for (ind=0; ind<SUBPOP; ind++) {
        //individual = subpopulation->individuals[ind];
        for (ind2=ind+1; ind2<SUBPOP; ind2++) {
            //individual2 = subpopulation->individuals[ind2];
            if (subpopulation->individuals[ind].fitness > subpopulation->individuals[ind2].fitness) {
                aux_individual = subpopulation->individuals[ind];
                subpopulation->individuals[ind] =subpopulation->individuals[ind2];
                subpopulation->individuals[ind2] = aux_individual;
            }
        }
    }
}

void migration(Population *population, int n_migrants) {
    int subpop, ind;
    for (subpop=0; subpop<NSUBPOP; subpop++) {
        if (subpop != NSUBPOP - 1) {
            for (ind = 0; ind<n_migrants; ind++)
                population->subpopulations[subpop+1].individuals[N - 1 - ind] = population->subpopulations[subpop].individuals[ind];
        } else {
            for (ind = 0; ind<n_migrants; ind++)
                population->subpopulations[0].individuals[N - 1 - ind] = population->subpopulations[subpop].individuals[ind];
        }
    }
}

double cost(char * indi){

    double result;
    double Ntur2;
    int i;
    double Ntur;

    Ntur=0.0;
    for(i=0; i<N;i++){
        if(indi[i]=='1'){
            Ntur++;
        }
    }

    Ntur2 = Ntur*Ntur;
    result = Ntur * (2.0/3.0 + 1.0/3.0 * exp(-0.00174*Ntur2));

    return result;
}

double power(double *u){

    int i;
    double result = 0.0;

    for(i=0; i<N; i++){
        result = result + 0.3*pow(u[i], 3.0);
    }

    return result;
}


void speeds_a(char * indi, double u0, double *u){
    int i, j;
    double xtur, ytur, xtur2, ytur2; // turbine position
    double dist;
    double uwake, sum;

    for(i=0; i<N; i++){

        if(indi[i]=='0'){ // If there is no turbine, speed 0 as it doesn't produce power and pass the the next
            u[i]=0;
            continue;
        }

        xtur=200 * (i%10) + 100; // It is located at the center of the cell, at (100,100) 
        ytur=200 *(i/10) + 100;  
        if(xtur==100){
            u[i]=u0;
            continue;
        }

        sum=0.0;
        for(j=0; j<N; j++){

            xtur2=200 * (j%10) + 100;
            ytur2=200 *(j/10) + 100;
            if(ytur2 != ytur){ // ignore it if it is not behind, in the wind's direction
                
                continue;

            }else{
                if(xtur2 < xtur && indi[j] != '0'){ // if it is ahead of the turbine of interest
                    dist = xtur - xtur2;
                    uwake = u0 * (1 - 2*a/((1 + alpha *(dist/R1))*(1 + alpha *(dist/R1))));
                    sum = sum + pow(1 - uwake/u0, 2.0);
                }
            }
        }
        u[i] = u0* (1-sqrt(sum));
    }
    return;
}

void speeds_bc(char * indi, double u0, double *u, double theta){
    int i, j;
    double xtur, ytur, xtur2, ytur2; // turbine position
    double xrot, yrot, xrot2, yrot2; // rotated positions
    double dist;
    double uwake, sum;

    for(i=0; i<N; i++){

        if(indi[i]=='0'){ // If there is no turbine, speed 0 as it doesn't produce power and pass the the next
            u[i]=0;
            continue;
        } 

        xtur=200 * (i%10) + 100; // It is located at the center of the cell, at (100,100) 
        ytur=200 *(i/10) + 100;

        // clockwise rotation (axis are non-standard)
        xrot = xtur*cos(theta) - ytur*sin(theta);
        yrot = xtur*sin(theta) + ytur*cos(theta);

        sum=0.0;
        for(j=0; j<N; j++){

            xtur2=200 * (j%10) + 100;
            ytur2=200 *(j/10) + 100;

            //clockwise rotation
            xrot2 = xtur2*cos(theta) - ytur2*sin(theta);
            yrot2 = xtur2*sin(theta) + ytur2*cos(theta);

            if(indi[j]=='0'){ //if there is no turbine at indi[j], no influence
                continue;
            }

            if(yrot > alpha*(xrot-xrot2) + Rr/2 + yrot2 || yrot < -alpha*(xrot-xrot2) -Rr/2 + yrot2 ){ // If it is not located inside the wake, ignore
                continue;
            }else{
                if(xrot2 < xrot ){ // wind from left to right, so only turbines to the left can affect others
                    
                    dist = xrot - xrot2; // formula does NOT use euclidean distance, only x-distance. 
                                         // wind speed "isoclines" are parallel to the y axis
                    uwake = u0 * (1 - 2*a/((1 + alpha *(dist/R1))*(1 + alpha *(dist/R1))));
                    sum = sum + pow(1 - uwake/u0, 2.0);
                }
            }
        }
        u[i] = u0* (1-sqrt(sum));
    }

    return;

}

void crossover(char *parent1, char *parent2, char *child1, char *child2) {
    int length = strlen(parent1);

    // Crossover points
    int crossoverPoint1 = rand() % length;
    int crossoverPoint2 = rand() % length;

    if (crossoverPoint1 > crossoverPoint2) {
        int temp = crossoverPoint1;
        crossoverPoint1 = crossoverPoint2;
        crossoverPoint2 = temp;
    }

    memcpy(child1, parent1, crossoverPoint1);
    memcpy(child2, parent2, crossoverPoint1);
    memcpy(child1 + crossoverPoint1, parent2 + crossoverPoint1, crossoverPoint2 - crossoverPoint1);
    memcpy(child2 + crossoverPoint1, parent1 + crossoverPoint1, crossoverPoint2 - crossoverPoint1);
    strcpy(child1 + crossoverPoint2, parent1 + crossoverPoint2);
    strcpy(child2 + crossoverPoint2, parent2 + crossoverPoint2);
}