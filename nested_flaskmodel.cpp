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
#define initial_population        100    // seed population size
#define genome_length             8      // length of genome bitstring
#define reproduction_thresh       120    // biomass threashold to reproduce
#define starve_thresh             50     // if a microbe's biomass drops below this, it dies
#define initial_biomass           80     // starting biomass for every microbe
#define max_consumption           10     // the maximum rate at which micrcomes can comsume
#define nutrient_conversion_eff   0.6    // efficiency of microbe conversion
#define maintainence_cost         1      // how much biomass it costs per timestep to live
#define p_mut                     0.01   // probability per gene of mutation in reproduction event
#define p_kill                    0.002  // probability of random death per timestep (not starvation)
#define prefered_abiotic          150.0  // ideal temperature for microbes (universal)
#define tau           		  0.02   // temperature sensitivity of microbes (universal)
    
   // FLASK PARAMETERS
#define num_flasks                4      // the number of mini flasks contained within the large flask
#define num_nutrients             4      // number of different nutrient types within experiment

   //FLOW PARAMETERS
#define nutrient_inflow	          1200.0 // number of units of each nutrient type inflowing per timestep to large flask
#define inflow_T_start            100.0  // temperature of inflow medium at start of experiment
#define abiotic_influx            0.2    // the percentage of the main flask liquid (by volume) exchanged for fresh inflow each timestep
#define mini_flask_exchange       0.2    // the percentage of the mini flask liquid (by volume) swapped for main flask liquid each timestep
#define main_flask_scale          1.0    // defines how large the main flask is in relation to the mini flask sizes
					 // total volume of liquid in large flask = num_flasks*main_flask_scale
					 // (excludes the volume taken up by the mini flasks)
					 // cannot have a value < mini_flask_exchange !!!
#define large_flask_exchange      mini_flask_exchange/main_flask_scale  // percentage of main flask liquid (by volume) that is  
									// swapped for mini flask liquid each timestep.
									// Between the main flask and each mini flask
									// large_flask_exchange/num_flasks % of the large flask liquid
									// is exchanged with the mini flask
									// If main_flask_scale = 1.0
									// then 20% of mini flask liquid = (20/num_flasks)% of main flask
									// liquid. If main_flask_scale = 0.5, then 20% of 
									// mini flask liquid = (40/num_flasks)% of main flask liquid
   //TEMPERATURE PERTURBATION PARAMETERS
#define t_perturbation = true						//Is there temperature perturbation? true or false values only
#define t_perturb_type = "smooth"					//Allowed values =  ["smooth", "step"]
#define inflow_T_end   = 1000						//If "smooth" chosen choose temperature inflow value at end of run
#define max_step_size  = 20						//If "step" chosen, choose the max temperature change during perturb
#define step_freq      = 200						//If "step" chosen, choose how often perturbation occurs (timesteps)

   //CULLING PARAMETERS
#define cull 	        	  true					//Do we cull? true or false values only
#define cull_percent    	  0.5					//How much do we cull each flask by? values from range [0,1]
#define cull_freq       	  100					//How frequently do we cull? (timesteps)
#define cull_number     	  1					//Number of flasks to cull. Ints from range [0,num_flasks]

   //TIME PARAMETERS
#define max_timesteps             1*pow(10,3)  // max number of generations per experiment
#define init_period               500          // initialisation period (timesteps)

   // DATA FILE NAMES
#define macro_filename(file_num)           "main_flask_macro_data_"+to_string(file_num)+".txt"
#define nutrient_filename(file_num)        "main_flask_nutrient_data_"+to_string(file_num)+".txt"
#define mini_filename(flask_num, filenum)  "mini_flask_"+to_string(flask_num)+"_run_"+to_string(file_num)+".txt"


using namespace std;


/*******************************************
		OBJECTS
*******************************************/

//********* microbe ****************
class microbe {    
    
    public:
    int genome;      // decimal representation of binary genome
    int population;  // population of this species
    int nutrient;    // amount of unconverted food 
    int biomass;     // amount of biomass within the species
    int waste;       // amount of un-excreted waste
    
};

//********* flask ****************
class flask {

     public:
     vector < microbe >  species;     // vector containing all species in this flask
     vector < double >   environment; // vector containing nutrient stocks in this flask
     double              temperature; // value of abiotic parameter of each individual flask

};

//********* large flask ****************
class large_flask {

     public:
     int                mini_flasks;  // number of flasks contained within main flask
     vector < double >  environment;  // vector containing nutrient stocks in main flask
     double             temperature;  // value of abiotic parameter of main flask

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


microbe generate_individual(int i, vector<microbe> &species, default_random_engine &generator) {

   double biomass_avg = species[i].biomass/(1.0*species[i].population);
   normal_distribution<double> biomass_dist( biomass_avg, biomass_avg*0.1 ); // distribution of biomass in pop
   int i_biomass = floor(biomass_dist(generator));

   double nutrient_avg = species[i].nutrient/(1.0*species[i].population);
   normal_distribution<double> nutrient_dist( nutrient_avg, nutrient_avg*0.1);
   int i_nutrient = floor(nutrient_dist(generator)); // we'll just round down as can't use half a nutrient

   double waste_avg = species[i].waste/(1.0*species[i].population);
   normal_distribution<double> waste_dist( waste_avg, waste_avg*0.1);
   int i_waste = floor(waste_dist(generator)); // we'll just round down as can't use half a nutrient

   if (i_biomass > species[i].biomass)   { i_biomass = species[i].biomass;   }
   if (i_biomass < 0)                    { i_biomass = 0;                    }
   if (i_nutrient > species[i].nutrient) { i_nutrient = species[i].nutrient; }
   if (i_nutrient < 0)                   { i_nutrient = 0;                   }
   if (i_waste > species[i].waste)       { i_waste = species[i].waste;       }
   if (i_waste < 0)                      { i_waste = 0;                      }

   microbe chosen_microbe;
   chosen_microbe.genome =     i;
   chosen_microbe.population = 1;
   chosen_microbe.nutrient =   1.0*i_nutrient;
   chosen_microbe.biomass =    1.0*i_biomass;
   chosen_microbe.waste =      1.0*i_waste;

   return chosen_microbe;

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
  
   vector <int> cons_ex_vector (8, 0);           // consumption excretion vector
   vector <int> temp(4, 0);
   vector < vector <int> > all_metabolisms;
    
   for (int i = 0; i < pow(2,genome_length); i++){
	int pos = 0;  // number nutrients consumed
	int neg = 0;  // number nutrients excreted

	for (int j = 0; j < temp.size(); j++) {
	   double a = 11*(2.0*drand48() - 1);       //random number between (-11,+11) actual range with floor / ceil is [-10,+10]
	   double b = 11*(2.0*drand48() - 1);
	   if (a >= 0) { a = floor(a);  }	   // integer values
	   else        { a = ceil(a);   }      
	   if (b >= 0) { b = floor(b);  }
	   else        { b = ceil(b);   }
	   temp[j] = a+b;                          // sum to get random number between [-10,+10]

	   if      (temp[j] > 0) { pos += 1;    }  // if positive - species consumes this nutrient   
	   else if (temp[j] < 0) { neg += 1;    }  // if negative - species excretes this nutrient
	   else                  { temp[j] = 0; }  // if 0 - species does not interact with this nutrient
	}

	if ( pos > 2 && neg > 2) { cout << "metabolism setup problem"; }  // bug check

	if (pos == 1) {   // we eat one nutrient
	   for ( int k = 0; k < temp.size(); k++ ) { if (temp[k] > 0) { temp[k] = 1; } }
	}
      
	if (pos == 2) {   // we eat two nutrients
	   int loc1 = -1;  // find the location of the positive numbers in the vector
	   int loc2 = -1;
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0)              { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0 && k != loc1) { loc2 = k; } }
	   temp[loc1] = temp[loc1]/GCD(temp[loc1],temp[loc2]);
	   temp[loc2] = temp[loc2]/GCD(temp[loc1],temp[loc2]);
	}

	if (pos == 3) {   // we eat 3 nutrients and must therefore excrete the other
	   int loc1 = -1;
	   int loc2 = -1;
	   int loc3 = -1;
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0)                          { loc1 = k; } }
	   for (int k = 0; k < temp.size(); k++) { if (temp[k] > 0 && k != loc1)             { loc2 = k; } }
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

   vector < double > large_nutrient_init (num_nutrients, 0.0);     // initially nutrient empty environment
   main_flask.temperature = 0.0;                                   // initial temperature at 0
   main_flask.mini_flasks = num_flasks;                            // set number of mini flasks
   main_flask.environment = large_nutrient_init;
   for (int k = 0; k < num_flasks; k++) {                          // set up each flask
	flask new_flask;
	microbe new_microbe;
	vector < microbe > species_init;
	for (int j = 0; j < initial_population; j++) {             // currently set up with diverse population
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
	new_flask.environment = nutrient_init; // start all flasks empty of nutrients
	new_flask.temperature = 0.0; 	       // each flask starts at 0.0
	flask_list.push_back(new_flask);       // append new flask to flask list
					       // filenames are initialised in the main code
   }
   return 0;
}

//*************** update all flask environments **********************************//
int update_all_flasks (large_flask &main_flask, vector < flask > &flask_list, double &inflow_T, double T_switch, int number_gens) {
   // outflow and inflow to large flasks
   for (int j = 0; j < num_nutrients; j++){ main_flask.environment[j]  = main_flask.environment[j]*(1.0-abiotic_influx);}     
   for (int j = 0; j < num_nutrients; j++){ main_flask.environment[j] += nutrient_inflow;            		        } 

   if (T_switch == 1.0 && t_perturbation ) { // update the T of inflow medium if perturbations are switched on
	if (t_perturb_type == "smooth")    { inflow_T -= (inflow_T_start - inflow_T_end)/(max_timesteps*1.0);} 
	else if (t_perturb_type == "step" && number_gens > 0 && fmod(number_gens, step_freq) == 0 ) { inflow_T += max_step_size*(1.0-2.0*drand48()); }
	else if (t_perturb_type != "step" && t_perturb_type != "smooth") { cout << "Invalid temperature perturbation option! \n"; }
   }       
	    
   main_flask.temperature = main_flask.temperature*(1.0-abiotic_influx) + inflow_T*abiotic_influx; // dilute main flask with fresh inflow

   double current_temperature_main = main_flask.temperature;   					   // record the temperature of the main flask 
   double temperature_changes = 0.0;                           					   // T change of main flask due to mini flask exchange
   vector < double > current_nutrients_main = main_flask.environment;
   double current_temperature_mini;

   for (int j = 0; j < num_flasks; j++) {
	current_temperature_mini  = flask_list[j].temperature;
	flask_list[j].temperature = current_temperature_mini*(1.0-mini_flask_exchange)  + current_temperature_main*mini_flask_exchange;
	temperature_changes      += current_temperature_main*(1.0-large_flask_exchange) + current_temperature_mini*large_flask_exchange;
	for (int k = 0; k < num_nutrients; k++) {
	   // net movement of nutrients from miniflask -> main flask
	   double nutrient_k_transfer = -1.0*mini_flask_exchange*flask_list[j].environment[k] + large_flask_exchange*current_nutrients_main[k]/(1.0*num_flasks);  
	   flask_list[j].environment[k] += nutrient_k_transfer;
	   main_flask.environment[k]    -= nutrient_k_transfer;
	}
   } 
   main_flask.temperature = temperature_changes/num_flasks;

   return 0;

}

/**************** kill event ****************************/
int kill_event (vector<microbe> &species, default_random_engine &generator){

   int i = chooseAgent(species); 								  // choses a microbe at random
   if (i > -1) {
	microbe chosen_microbe = generate_individual(i, species, generator); 			  // generate the chosen microbe
	if (chosen_microbe.biomass <= starve_thresh) { 						  // death event via starvation
	   species[i].population--;
	   species[i].biomass -= chosen_microbe.biomass; 					  // remove biomass of dead microbe 
	   species[i].nutrient -= chosen_microbe.nutrient; 					  // remove unconverted nutrients 
	   if (species[i].nutrient < 0)    { species[i].nutrient = 0;                           }
	   if (species[i].biomass < 1)     { species[i].biomass = 0; species[i].population = 0; }
	   if (species[i].population == 0) { species.erase(species.begin() + i);                } // remove from list if extinct
	} else if (drand48() <= p_kill && species[i].population > 0) {                            // death event via p_kill probability
	   species[i].population--;
	   species[i].biomass -= chosen_microbe.biomass; 					  // remove biomass of dead microbe 
	   species[i].nutrient -= chosen_microbe.nutrient; 					  // remove unconverted nutrients
	   if (species[i].nutrient < 0)    { species[i].nutrient = 0;                           }
	   if (species[i].biomass < 1)     { species[i].biomass = 0; species[i].population = 0; }
	   if (species[i].population == 0) { species.erase(species.begin() + i);                } // remove from list if extinct
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
	double factor_i = tau*sqrt(pow(current.temperature - prefered_abiotic, 2.0));
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
	microbe chosen_microbe = generate_individual(i, species, generator);       // generate the chosen microbe
	while ( chosen_microbe.nutrient >= 5) {     
	   species[i].nutrient -= 5;
	   chosen_microbe.nutrient -= 5;
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
	microbe chosen_microbe = generate_individual(i, species, generator); // generate the chosen microbe
	int waste_count = 0;
	for (int k = 0; k < num_nutrients; k++) { waste_count += metabolism[species[i].genome][k+num_nutrients]; }
	while (chosen_microbe.waste >= waste_count && chosen_microbe.waste > 0) {
	   for (int k = 0; k < num_nutrients; k++) {
		species[i].waste -=  metabolism[species[i].genome][k+num_nutrients];
		environment[k] +=  metabolism[species[i].genome][k+num_nutrients];
	   }
	   chosen_microbe.waste -= waste_count;
	}
   }
   return 0;
}

/**************** reproduction event ****************************/
int reproduction_event(vector<microbe> &species, default_random_engine &generator) {

   int i = chooseAgent(species);
   if (i > -1){
	microbe chosen_microbe = generate_individual(i, species, generator); // generate the chosen microbe
	if (chosen_microbe.biomass >= reproduction_thresh) { 	             // reproduction occurs
 	   bitset<genome_length> mutant_genome(species[i].genome);
	   int did_we_mutate = 0; 					     // do we mutate?
 	   for (int j = 0; j < genome_length; j++){
		if (drand48() <= p_mut){ 				     // probability of mutation per gene
		   did_we_mutate = 1;
		   if (mutant_genome[j] == 1) { mutant_genome[j] = 0; }      // flip gene if mutated
		   else { mutant_genome[j] = 1; }
		}
	   }
	   if (did_we_mutate == 1) {
		int mutant_number = int(mutant_genome.to_ulong());
		int species_exists = 0;
		for (int q = 0; q < species.size(); q++){ 			 // check to see if species exists
		   if (species[q].genome == mutant_number){
			species[q].population++;
			species[q].biomass += int(chosen_microbe.biomass / 2.0); // half biomass goes to new mutant
			species[i].biomass -= int(chosen_microbe.biomass / 2.0);
			species_exists = 1;
			break;
		   }
		}
		if (species_exists == 0){ // add species if it doesn't exist
		   microbe temp_mutant;
		   temp_mutant.genome = mutant_number;
		   temp_mutant.nutrient = 0; 				     // no nutrient count to begin with
		   temp_mutant.biomass = int(chosen_microbe.biomass / 2.0);  // half biomass goes to new mutant
		   species[i].biomass -= int(chosen_microbe.biomass / 2.0);
		   temp_mutant.population = 1; 				     // initial population of 1
		   temp_mutant.waste = 0;
		   species.push_back(temp_mutant);
		}
	   } else { species[i].population++; } // no mutation
	}
   }
   return 0;
}

/************************ cull event *********************************/
int cull_event(vector <flask> &flask_list, default_random_engine &generator) {

	vector < int > cull_list;
	int c; 	      // records the flask randomly selected for culling
	int cull_pop; // records how many individuals get culled in a particular flask
	int f; 	      // flask marker
	
	if      ( cull_number <= 0 || cull_number > num_flasks ) { cout << "Invalid cull set up \n"; return 0; }
	else if ( cull_number == num_flasks ) { for (int f = 0; f < num_flasks; f++) { cull_list.push_back(f); } }
	else    {
	   while (cull_list.size() < cull_number) {  // generate list of flasks to cull
		c = floor( drand48()*num_flasks );
		if ( find(cull_list.begin(), cull_list.end(), c) == cull_list.end() ) { cull_list.push_back(c); }
	   }
	}

	for (int v = 0; v < cull_list.size(); v++) {  // cull the selected flasks
	   f = cull_list[v];
	   cull_pop = 0;
	   for (int s = 0; s < flask_list[f].species.size(); s++) { cull_pop += flask_list[f].species[s].population; }
	   cull_pop = floor(cull_pop*cull_percent);  // calculate population to cull
	   
	   for (int p = 0; p < cull_pop; p++) {      // kill individuals until cull requirement satisfied
  		int cull_i = chooseAgent(flask_list[f].species);  // choses a microbe at random to cull
   		if (cull_i > -1) {
		    microbe chosen_microbe = generate_individual(cull_i, flask_list[f].species, generator);  // generate the chosen microbe
	   	    flask_list[f].species[cull_i].population--;
	   	    flask_list[f].species[cull_i].biomass -= chosen_microbe.biomass; 		          // remove biomass of dead microbe 
	   	    flask_list[f].species[cull_i].nutrient -= chosen_microbe.nutrient; 		  // remove unconverted nutrients 
	   	    if (flask_list[f].species[cull_i].nutrient < 0) { flask_list[f].species[cull_i].nutrient = 0;   }
	   	    if (flask_list[f].species[cull_i].biomass < 1)  { flask_list[f].species[cull_i].biomass = 0;
								      flask_list[f].species[cull_i].population = 0; }
	   	    if (flask_list[f].species[cull_i].population == 0) { flask_list[f].species.erase(flask_list[f].species.begin() + cull_i); } //remove extinct species from list
		}
	   }
	}

	return 0;

}


/****************************************************************************************
**********		               MAIN CODE                               **********
****************************************************************************************/

int main(int argc, char **argv) {

   // INITIALISE INFLOW TEMPERATURE
   double inflow_T = inflow_T_start;

   // RANDOM NUMBER GENERATORS
   int t = atoi(argv[1]); 	       // seed for randomness
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
   ofstream macro_data (macro_filename(file_num));
   ofstream nutrient_data (nutrient_filename(file_num));
   vector < ofstream > mini_flask_data;
   for (int f = 0; f < num_flasks; f++){
	mini_flask_data.emplace_back(ofstream{ mini_filename(f, filenum) });
   }    

   /****************************************************************************************
					INITIALISE ENVIRONMENT 
   *****************************************************************************************/

    while (init_counter < init_period) {
	update_all_flasks(main_flask, flask_list, inflow_T, 0.0, 0);
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

	// record data from main flask
	macro_data << number_gens << " " << num_alive_flasks << " " << main_flask.temperature << " " << inflow_T << endl;
	nutrient_data << number_gens;
	for (int n = 0; n < num_nutrients; n++){
	   nutrient_data << " " << main_flask.environment[n];
	}
	nutrient_data << endl;
	
	num_alive_flasks = 0;  // reset counter for number of mini flasks hosting live biospheres for this timestep

	//record data from the mini flasks to individual files
	for (int f = 0; f < num_flasks; f++) {
	   int local_pop = 0;
	   for (int s = 0; s < flask_list[f].species.size(); s++) { local_pop += flask_list[f].species[s].population; }
	   mini_flask_data[f] << number_gens << " " << flask_list[f].species.size() << " " << local_pop << " " << flask_list[f].temperature; 
	   for (int n = 0; n < num_nutrients; n++) {
		mini_flask_data[f] << " " << flask_list[f].environment[n];
	   }
	   mini_flask_data[f] << endl;
	}

	/************************************************************************************
	                          NUTRIENT FLOW
	**************************************************************************************/

	update_all_flasks(main_flask, flask_list, inflow_T, 1.0, number_gens);


	/*****************************************************************************************
				CULL EVENT
	*****************************************************************************************/

	if (cull && fmod(number_gens,cull_freq) == 0) { cull_event(flask_list, generator); }

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
    for (int f = 0; f < num_flasks; f++) {
	mini_flask_data[f].close();
    }

    return 0;
}
