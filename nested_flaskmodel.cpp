#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <bitset>
#include <string>
#include <sstream>

/*******************************************
        CONSTANT PARAMETERS
*******************************************/

   // MICROBE PARAMETERS
#define initial_population 100         // seed population size
#define genome_length 8                // length of genome bitstring
#define reproduction_thresh 120        // biomass threashold to reproduce
#define starve_thresh 50               // if a microbe's biomass drops below this, it dies
#define initial_biomass 80             // starting biomass for every microbe
#define max_consumption 10             // the maximum rate at which micrcomes can comsume
#define nutrient_conversion_eff 0.6 // efficiency of microbe conversion
#define maintainence_cost 1            // how much biomass it costs per timestep to live
#define p_mut 0.01                  // probability per gene of mutation in reproduction event
#define p_kill 0.002                //probability of random death per timestep (not starvation)
#define prefered_abiotic 150.0      // ideal temperature for microbes (universal)
#define abiotic_scaling 0.0         // temperature sensitivity of microbes (universal)
    
   // FLASK PARAMETERS
#define num_nutrients 4
#define num_abiotic 1  // if we change this we need to change other code...
#define num_flasks 4   // the number of flasks

   //FLOW PARAMETERS
#define percentage_outflow 0.2     // to calculate outflow
#define max_nutrient_inflow 150.0
#define abiotic_start 100.0
#define abiotic_end 100.0
#define abiotic_influx 0.2    // the percentage that temp will change each timestep
#define mini_flask_inflow 0.2
#define main_flask_scale 1.0 

    //TIME PARAMETERS
#define max_timesteps 10*pow(10,4)   // max number of generations per experiment
#define init_period 5000             // initialisation period (timesteps)


using namespace std;


/*******************************************
		OBJECTS
*******************************************/

//********* microbe ****************
class microbe {  // class for our microbe species
    
    public:
    int genome;
    int population;
    int nutrient; // amount of unconverted food 
    int biomass;
    int waste;  // amount of un excreted waste
    
};

//********* flask ****************
class flask {

     public:
     vector < microbe > species;    // vector containing all species in this flask
     vector < double > environment; // vector containing nutrient stocks in this flask
     double temperature;            // value of abiotic parameter of each individual flask
     int iteration;		    // number of iterations per timestep for flask (= total population at start of timestep)
     string file;

};

//********* large flask ****************
class large_flask {

     public:
     int mini_flasks;  		     // number of flasks contained within main flask
     vector < double > environment;  // vector containing nutrient stocks in main flask
     double temperature;             // value of abiotic parameter of main flask

};


/*******************************************
		FUNCTIONS
*******************************************/

//********** sorting rule ****************
struct myclass {  // rank species by population, largest first
    bool operator() ( microbe i, microbe j) { return (i.population > j.population);}
} sorting_rule;


//********* choose individual *********
int chooseAgent(vector< microbe > &species) {
    
  //    vector<double> prob_array (species.size());
    int p = 0;
    for (int s = 0; s < species.size(); s++) { p += species[s].population; }
    double r = 0;
    
    double p_num = drand48();
    
    for (int j = 0; j < species.size(); j++){
        r+= species[j].population;
        if (p*p_num <= r){
            return j;
        }
        
    }
    return -1; // this should never happen   
}

//********* greater common denomenator ***************
int GCD(int a, int b) {
    while( 1 ) {
        a = a % b;
	if( a == 0 ) { return b; }
	b = b % a;
        if( b == 0 ) { return a; }
    }
}

//******** nutrient genome interactions (metabolisms) **************
vector < vector <int> > nutrient_genome_interactions(default_random_engine &generator){
  
   vector <int> cons_ex_vector (8, 0);  // consumption excretion vector
   vector <int> temp(4, 0);
   vector < vector <int> > all_metabolisms;
    
   for (int i = 0; i < pow(2,genome_length); i++){
	int pos = 0;  // number nutrients consumed
	int neg = 0;  // number nutrients excreted

	for (int j = 0; j < temp.size(); j++) {
	   double a = 11*(2.0*drand48() - 1); //random number between (11,-11) actual range with floor / ceil is [10,-10]
	   double b = 11*(2.0*drand48() - 1);
	   if (a > 0) { a = floor(a); }
	   else       { a = ceil(a);  }
	   if (b > 0) { b = floor(b); }
	   else       { b = ceil(b);  }
	   temp[j] = a+b; 

	   if      (temp[j] > 0) { pos += 1;    }    
	   else if (temp[j] < 0) { neg += 1;    } 
	   else                  { temp[j] = 0; }
	}

	if ( pos > 2 && neg > 2) { cout << "metabolism setup problem"; }  // bug check

	if (pos == 1) {   // we eat one nutrient
	   for ( int k = 0; k < temp.size(); k++ ) { if (temp[k] > 0) { temp[k] = 1; } }
	}
      
	if (pos == 2) {   // we eat two nutrients
	   int loc1 = -1;  // find the location of the positive numbers in the vector
	   int loc2 = -1;
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0             ) { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0 && k != loc1) { loc2 = k; } }
	   temp[loc1] = temp[loc1]/GCD(temp[loc1],temp[loc2]);
	   temp[loc2] = temp[loc2]/GCD(temp[loc1],temp[loc2]);
	}

	if (pos == 3) {   // we eat 3 nutrients and must therefore excrete the other
	   int loc1 = -1;
	   int loc2 = -1;
	   int loc3 = -1;
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0                         ) { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0 && k != loc1            ) { loc2 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0 && k != loc1 && k!= loc2) { loc3 = k; } }
	
	   int gcd_1 = GCD(temp[loc1],temp[loc2]);
	   int gcd_2 = GCD(temp[loc1],temp[loc3]);
	   int gcd_3 = GCD(temp[loc2],temp[loc3]);

	   if (gcd_1 == gcd_2 && gcd_1 > 0) {
		temp[loc1] = temp[loc1]/gcd_1;
		temp[loc2] = temp[loc2]/gcd_1;
		temp[loc3] = temp[loc3]/gcd_1;

	   } else if (gcd_1 == gcd_3 && gcd_1 > 0) {
		temp[loc1] = temp[loc1]/gcd_1;
		temp[loc2] = temp[loc2]/gcd_1;
		temp[loc3] = temp[loc3]/gcd_1;
	   }
	}

	if (neg == 1) { // we excrete 1 nutrient
	   for ( int k = 0; k < temp.size(); k++ ) { if (temp[k] < 0) { temp[k] = -1; } }
	}

	if (neg == 2) {  // we excrete 2 nutrients
	   int loc1 = -1;
	   int loc2 = -1;
	
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] < 0)              { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] < 0 && k != loc1) { loc2 = k; } }
	   temp[loc1] = temp[loc1]/GCD(-1*temp[loc1],-1*temp[loc2]);
	   temp[loc2] = temp[loc2]/GCD(-1*temp[loc1],-1*temp[loc2]);
	}

	if (neg == 3) {  // we excrete 3 nutrients

	   int loc1 = -1; 
	   int loc2 = -1;
	   int loc3 = -1;
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] < 0)                          { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] < 0 && k != loc1)             { loc2 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] < 0 && k != loc1 && k!= loc2) { loc3 = k; } }
	
	   int gcd_1 = GCD(-1*temp[loc1],-1*temp[loc2]);
	   int gcd_2 = GCD(-1*temp[loc1],-1*temp[loc3]);
	   int gcd_3 = GCD(-1*temp[loc2],-1*temp[loc3]);

	   if (gcd_1 == gcd_2 && gcd_1 > 0) {
		temp[loc1] = temp[loc1]/gcd_1;
		temp[loc2] = temp[loc2]/gcd_1;
		temp[loc3] = temp[loc3]/gcd_1;

	   } else if (gcd_1 == gcd_3 && gcd_1 > 0) {
		temp[loc1] = temp[loc1]/gcd_1;
		temp[loc2] = temp[loc2]/gcd_1;
		temp[loc3] = temp[loc3]/gcd_1;
	   }
	}

	if ( pos > 0 && neg > 0 && pos < 4 && neg < 4 ) {
	   for (int m = 0; m < num_nutrients; m ++ ) {
		if ( temp[m]  > 0 ) {  cons_ex_vector[m] = temp[m];     cons_ex_vector[m+num_nutrients] = 0;          }
		if ( temp[m]  < 0 ) {  cons_ex_vector[m] = 0;           cons_ex_vector[m+num_nutrients] = -1*temp[m]; }
		if ( temp[m] == 0 ) {  cons_ex_vector[m] = 0;           cons_ex_vector[m+num_nutrients] = 0;          }
	   }
	} else { // metabolism not viable - must both consume and excrete!
	   for (int m = 0; m < 2*num_nutrients; m++) { cons_ex_vector[m] = 0; }
	}
	all_metabolisms.push_back(cons_ex_vector);        
   }
   return all_metabolisms;  
}

//************* metaboblism temperature impacts ***********
vector < double > abiotic_genome_interactions() {
    //determines the temperature impact due to metabolism for every genome
    
    vector <double> abiotic_ints (pow(2,genome_length),0);
    for (int i = 0; i < pow(2,genome_length); i++){ abiotic_ints[i] = 2.0*drand48() - 1.0; }
    return abiotic_ints;  
}

/**************** initialise flasks ****************************/
int initialise_flasks (large_flask &main_flask, vector <flask> &flask_list) {

   vector < double > large_nutrient_init (num_nutrients, 0.0);
   main_flask.temperature = 0.0;
   main_flask.mini_flasks = num_flasks;
   main_flask.environment = large_nutrient_init;
   for (int k = 0; k < num_flasks; k++) {  // set up each flask
	flask new_flask;
	microbe new_microbe;
	vector < microbe > species_init; //(1, new_microbe);
	for (int j = 0; j < initial_population; j++) { // currently set up with diverse population
	   int genome_new = floor(drand48()*pow(2,genome_length)); // randomly generate microbes
	   int already_exists = 0;
	   for (int k = 0; k < species_init.size(); k++){

		if (genome_new == species_init[k].genome) {
		   species_init[k].population++;
		   species_init[k].biomass += initial_biomass;
		   already_exists = 1;
		}
	   }
	   if (already_exists == 0) {
		new_microbe.population = 1;
		new_microbe.genome = genome_new;
		new_microbe.nutrient = 0;
		new_microbe.biomass = initial_biomass;
		new_microbe.waste = 0;
		species_init.push_back(new_microbe);
	   }
	}
	new_flask.species = species_init;
	vector < double > nutrient_init (num_nutrients, 0.0);	
	new_flask.environment = nutrient_init; // start all flasks empty
	new_flask.temperature = 0.0; 	 // each flask starts at 0.0
	flask_list.push_back(new_flask);
	// filenames are initialised in the main code
   }
   return 0;
}

//*************** update all flask environments **********************************//
int update_all_flasks (large_flask &main_flask, vector < flask > &flask_list, double abiotic_env) {
   // outflow and inflow
   for (int j = 0; j < num_nutrients; j++){ main_flask.environment[j] = main_flask.environment[j]*(1.0-abiotic_influx);}     
   for (int j = 0; j < num_nutrients; j++){ main_flask.environment[j] += max_nutrient_inflow; } 

   abiotic_env -= (abiotic_start - abiotic_end)/max_timesteps;
	    
   main_flask.temperature = main_flask.temperature*(1.0-abiotic_influx) + abiotic_env*abiotic_influx; // dilute with fresh inflow

   double current_temperature_main = main_flask.temperature;
   double temperature_changes = 0.0;
   vector < double > current_nutrients_main = main_flask.environment;
   vector < double > current_nutrients_mini;
   double current_temperature_mini;

   for (int j = 0; j < num_flasks; j++) {
	current_temperature_mini = flask_list[j].temperature;
	flask_list[j].temperature = flask_list[j].temperature*(1.0-mini_flask_inflow) + current_temperature_main*mini_flask_inflow;
	temperature_changes += (current_temperature_mini*mini_flask_inflow + current_temperature_main*(1.0-mini_flask_inflow))/(num_flasks);

	for (int k = 0; k < num_nutrients; k++) {
	   double nutrient_k_transfer = -1.0*mini_flask_inflow*flask_list[j].environment[k] + mini_flask_inflow*current_nutrients_main[k]/(1.0*num_flasks);
	   flask_list[j].environment[k] += nutrient_k_transfer;
	   main_flask.environment[k] -= nutrient_k_transfer;
	}
   } 
   main_flask.temperature = temperature_changes;

   return 0;

}

/**************** kill event ****************************/
int kill_event (vector<microbe> &species, default_random_engine &generator){

   int i = chooseAgent(species); // choses a microbe at random
   if (i > -1) {
	double average_biomass = species[i].biomass/(1.0*species[i].population);
	normal_distribution<double> biomass_dist( average_biomass, average_biomass*0.1 ); // distribution of biomass in pop
	int i_biomass = floor(biomass_dist(generator));
	double nutrient_avg = 1.0*species[i].nutrient/species[i].population;
	normal_distribution<double> nutrient_species_dist( nutrient_avg, nutrient_avg*0.1);
	int i_nutrient = floor(nutrient_species_dist(generator)); // we'll just round down as can't use half a nutrient

	if (i_biomass <= starve_thresh) { // death event via starvation
	   species[i].population--;
	   species[i].biomass -= i_biomass; // remove biomass of dead microbe 
	   species[i].nutrient -= i_nutrient; // remove unconverted nutrients 
	   if (species[i].nutrient < 0) { species[i].nutrient = 0; }
	   if (species[i].biomass < 1) { species[i].biomass = 0; species[i].population = 0; }
	   if (species[i].population == 0){ species.erase(species.begin() + i); } // remove from list if extinct
	} else if (drand48() <= p_kill && species[i].population > 0) {  // death event via p_kill probability
	   species[i].population--;
	   species[i].biomass -= i_biomass; // remove biomass of dead microbe 
	   species[i].nutrient -= i_nutrient; // remove unconverted nutrients
	   if (species[i].nutrient < 0) { species[i].nutrient = 0; }
	   if (species[i].biomass < 1) { species[i].biomass = 0; species[i].population = 0; }
	   if (species[i].population == 0){ species.erase(species.begin() + i);} // remove from list if extinct
	}
   }
   return 0;
}

/**************** maintenance event ****************************/
int maintenance_event (vector<microbe> &species) {

   int i = chooseAgent(species);
   if (i > -1 && species[i].biomass > 0) { species[i].biomass--; }
   return 0;
}

/**************** metabolism event ****************************/
int metabolism_event(flask &current, vector<vector<int>> metabolism) {

   int i = chooseAgent(current.species);
   if( i > -1) {
	double factor_i = abiotic_scaling*sqrt(pow(current.temperature - prefered_abiotic, 2.0));
	double satisfaction = exp (-1.0*pow(factor_i,2.0));
	int total_count_eat = floor(max_consumption * satisfaction);
	int minimum_count_eat = 0; // the minumum total number of nutrients microbe can intake
	int max_count_eat = 0; 
	for (int k = 0; k < num_nutrients; k++) { minimum_count_eat += metabolism[current.species[i].genome][k];}
	max_count_eat = minimum_count_eat;
	bool enough_nutrient = true;
	while ( max_count_eat <= total_count_eat && minimum_count_eat > 0 && total_count_eat > 0 ) {

	   for (int k = 0; k < num_nutrients; k++) {
		if (floor(current.environment[k]) < metabolism[current.species[i].genome][k]) { enough_nutrient = false; }
		// if there is not enough nutrient to satisfy metabolism in correct proportions, cannot eat
	   }
	   if (enough_nutrient == false) { break; } // exit if cannot eat
	   for (int k = 0; k < num_nutrients; k++) {
		current.environment[k] -= metabolism[current.species[i].genome][k];
		current.species[i].nutrient += metabolism[current.species[i].genome][k];
		if( current.environment[k] < 0) { cout << "metabolism event problem" << endl; }
	   }
	   max_count_eat += minimum_count_eat; // can we have more than one eating event?
	}
   }
   return 0;
}

/**************** biomass creation event ****************************/
int biomassCreation_event (vector<microbe> &species, double &temperature, vector<double> T_impact, default_random_engine &generator) {

   int i = chooseAgent(species);
   if (i > -1){
	double species_nutrient_avg = 1.0*species[i].nutrient/species[i].population;
	normal_distribution<double> nutrient_species_dist( species_nutrient_avg, species_nutrient_avg*0.1);
	int nutrient_available = floor(nutrient_species_dist(generator)); // we'll just round down as can't use half a nutrient
	while ( nutrient_available >= 5) {     
	   species[i].nutrient -= 5;
	   nutrient_available -= 5;
	   species[i].biomass += int(5.0*nutrient_conversion_eff);
	   species[i].waste += int(5*(1.0 - nutrient_conversion_eff));
	   temperature += 5.0*nutrient_conversion_eff*T_impact[species[i].genome]; // the effect on abiotic environment
	}
   }
   return 0;
}

/**************** waste event ****************************/
int waste_event (vector<microbe> &species, vector<double> &environment, vector<vector<int>> metabolism, default_random_engine &generator){

   int i = chooseAgent(species);
   if (i > -1){
	double species_waste_avg = 1.0*species[i].waste/species[i].population;
	normal_distribution<double> waste_species_dist( species_waste_avg, species_waste_avg*0.1);
	int waste_available = round(waste_species_dist(generator));
	if (waste_available > species[i].waste) { waste_available = species[i].waste; } 
	int waste_count = 0;
	for (int k = 0; k < num_nutrients; k++) { waste_count += metabolism[species[i].genome][k+num_nutrients]; }
	while (species[i].waste >= waste_count && species[i].waste > 0) {
	   for (int k = 0; k < num_nutrients; k++) {
		species[i].waste -=  metabolism[species[i].genome][k+num_nutrients];
		environment[k] +=  metabolism[species[i].genome][k+num_nutrients];
	   }
	}
   }
   return 0;
}

/**************** reproduction event ****************************/
int reproduction_event(vector<microbe> &species, default_random_engine &generator) {

   int i = chooseAgent(species);
   if (i > -1){
	double average_biomass = (1.0*species[i].biomass)/species[i].population; // average biomass per indiviual
	normal_distribution<double> biomass_species_dist( average_biomass, average_biomass*0.1 );
	int i_biomass = floor(biomass_species_dist(generator));
	if (i_biomass >= reproduction_thresh) {  // reproduction occurs
 	   bitset<genome_length> mutant_genome(species[i].genome);
	   int did_we_mutate = 0; // do we mutate?
 	   for (int j = 0; j < genome_length; j++){
		if (drand48() <= p_mut){ // probability of mutation per gene
		   did_we_mutate = 1;
		   if (mutant_genome[j] == 1) { mutant_genome[j] = 0; } // flip gene if mutated
		   else { mutant_genome[j] = 1; }
		}
	   }
	   if (did_we_mutate == 1) {
		int mutant_number = int(mutant_genome.to_ulong());
		int species_exists = 0;
		for (int q = 0; q < species.size(); q++){ // check to see if species exists
		   if (species[q].genome == mutant_number){
			species[q].population++;
			species[q].biomass += int(i_biomass / 2.0); // half biomass goes to new mutant
			species[i].biomass -= int(i_biomass / 2.0);
			species_exists = 1;
			break;
		   }
		}
		if (species_exists == 0){ // add species if it doesn't exist
		   microbe temp_mutant;
		   temp_mutant.genome = mutant_number;
		   temp_mutant.nutrient = 0; //  no nutrient count to begin with
		   temp_mutant.biomass = int(i_biomass / 2.0);  // half biomass goes to new mutant
		   species[i].biomass -= int(i_biomass / 2.0);
		   temp_mutant.population = 1; // initial population of 1
		   temp_mutant.waste = 0;
		   species.push_back(temp_mutant);
		}
	   } else { species[i].population++; } // no mutation
	}
   }
   return 0;
}


/****************************************************************************************
**********		               MAIN CODE                               **********
****************************************************************************************/

int main(int argc, char **argv) {

   // INITIALISE INFLOW TEMPERATURE
   double abiotic_env = abiotic_start;

   // RANDOM NUMBER GENERATORS
   int t = atoi(argv[1]); // seed for randomness
   srand48 (t);
   mt19937 rng(t);
   default_random_engine generator;    // used for our distributions
   generator.seed(t);                  //provide seed for randomness in different runs
    
   // TIME VARIABLES
   int number_gens = 0;                // track the number of generations / timesteps
   int init_counter = 0;
   int num_alive_flasks = 0;           // track number of flasks hosting life throughout experiment
    
   // DATA FILES NUMBER
   int file_num = atoi(argv[2]);

   vector < vector< int > >  n_g_interacts = nutrient_genome_interactions(generator); 
   // initialises which metabolism each genome codes for

   vector < double > a_g_interacts = abiotic_genome_interactions(); 
   // how each genome affects the local temperature

   // FLASK SET UP
   large_flask main_flask;
   vector < flask > flask_list;
   initialise_flasks(main_flask, flask_list);

   // FILE SET UP
   ofstream macro_data ("main_flask_macro_data_"+to_string(file_num)+".txt");
   ofstream nutrient_data ("main_flask_nutrient_data_"+to_string(file_num)+".txt");
   for (int f = 0; f < num_flasks; f++){
	ostringstream filename;
	filename << "mini_flask_"+to_string(f)+"_run_"+to_string(file_num)+".txt";
	flask_list[f].file = filename.str();
	ofstream temp_file (flask_list[f].file);
	temp_file.close();
   }
    

   /****************************************************************************************
					INITIALISE ENVIRONMENT 
   *****************************************************************************************/

    while (init_counter < init_period) {
	update_all_flasks(main_flask, flask_list, abiotic_env);
	init_counter++;
    }

   /****************************************************************************************
					MAIN EXPERIMENT 
   *****************************************************************************************/
    
    while (number_gens < max_timesteps) {

	// sorts species by order of population, largest first
	for (int f = 0; f < num_flasks; f++) { 
	   stable_sort (flask_list[f].species.begin(), flask_list[f].species.end(), sorting_rule);
	}
            
	/* ********************************************************************************
                                    RECORD DATA
	********************************************************************************/

	macro_data << number_gens << " " << num_alive_flasks << " " << main_flask.temperature << " " << abiotic_env << endl;
	num_alive_flasks = 0;

	for (int f = 0; f < num_flasks; f++) {
	   int local_pop = 0;
	   for (int s = 0; s < flask_list[f].species.size(); s++) { local_pop += flask_list[f].species[s].population; }
	   if (local_pop > 0) {
		ofstream temp_file;
		temp_file.open(flask_list[f].file, ofstream::app);
		temp_file << number_gens << " " << flask_list[f].species.size() << " " << local_pop << " " << flask_list[f].temperature << endl;
		temp_file.close();
	   } 
	}

	/************************************************************************************
	                          NUTRIENT FLOW
	**************************************************************************************/

	update_all_flasks(main_flask, flask_list, abiotic_env);

	/********************************************************************************
				LOOP OVER MINI FLASKS
	*********************************************************************************/

	for (int f = 0; f < num_flasks; f++) {

	   int iterations = 0;
	   for (int k = 0; k < flask_list[f].species.size(); k++) { 
		iterations += flask_list[f].species[k].population; 
	   }
	     
	   for (int l = 0; l < iterations; l++) {

        	/* ********************************************************************************
                                	    KILL EVENT
        	 ********************************************************************************/

		kill_event(flask_list[f].species, generator);

        	/* ********************************************************************************
             	                       MAINTENANCE COST EVENT
        	 ********************************************************************************/

		maintenance_event(flask_list[f].species);

        	/* ********************************************************************************
                         	         METABOLISM EVENT
         	********************************************************************************/
        
		metabolism_event(flask_list[f], n_g_interacts);

        	/* ********************************************************************************
            	                      BIOMASS CREATION EVENT
        	 ********************************************************************************/

		biomassCreation_event(flask_list[f].species, flask_list[f].temperature, a_g_interacts, generator);

		/*********************************************************************************
                    	                  WASTE EVENT
	 	**********************************************************************************/

		waste_event(flask_list[f].species, flask_list[f].environment, n_g_interacts, generator);
	
		/* ********************************************************************************
               		               REPRODUCTION EVENT
         	********************************************************************************/

		reproduction_event(flask_list[f].species, generator);

		/* END OF MICROBE EVENTS */
	   }

	   if (flask_list[f].species.size() > 0) { num_alive_flasks++; }

	} /* END OF LOOPING OVER MINI FLASKS */

	//if (num_alive_flasks == 0) { break; } // end code when all flasks are dead
        number_gens++;
        
    }
    
    /* CLOSE FILES */
    macro_data.close();
    nutrient_data.close();

    return 0;
}
