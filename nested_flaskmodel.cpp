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
#include <boost/shared_ptr.hpp>
#include <sstream>

using namespace std;


//********* MICROBES ****************

class microbe {  // class for our microbe species
    
    public:
    int genome;
    int population;
    int nutrient; // amount of unconverted food 
    int biomass;
    int waste;  // amount of un excreted waste
    
};

class flask {

     public:
     vector < microbe > species;  // vector containing all species in this flask
     vector < double > nutrient;  // vector containing nutrient stocks in this flask
     double temperature;          // value of abiotic parameter of each individual flask
     int iteration;		  // number of iterations per timestep for flask (same as total population at start of timestep)
     //shared_ptr<ofstream> file;	  // file handle for mini flask data
     string file;

};

class large_flask {

     public:
     int mini_flasks;  		  // number of flasks contained within main flask
     vector < double > nutrient;  // vector containing nutrient stocks in main flask
     double temperature;          // value of abiotic parameter of main flask

};

//********** SORTING ****************

struct myclass {  // rank species by population, largest first
    bool operator() ( microbe i, microbe j) { return (i.population > j.population);}
} sorting_rule;


//********* CHOOSE INDIVIDUAL *********

// code to chose a random individual from the system
int chooseAgent(vector< microbe > &species, int p) {
    
  //    vector<double> prob_array (species.size());
    double r = 0;
    
    double p_num = drand48();
    
    for (int j = 0; j < species.size(); j++){
        r+= species[j].population;
        if (p*p_num <= r){
            return j;
        }
        
    }
    
    // cout << "\n SHOULD NEVER SEE THIS! \n";
    return -1; // this should never happen 
    
}

//********* GREATER COMMON DENOMENATOR ***************

int GCD(int a, int b)
{
    while( 1 )
    {
        a = a % b;
		if( a == 0 )
			return b;
		b = b % a;

        if( b == 0 )
			return a;
    }
}

//******** NUTRIENT GENOME INTERACTIONS **************

//determines the metabolism for every genome

vector < vector <int> > nutrient_genome_interactions(int genome_length, int num_nutrients, default_random_engine &generator){
    // this maps which geomones will code for which metabolims
  
    vector <int> cons_ex_vector (8, 0);  // consumption excretion vector
    vector <int> temp(4, 0);
    
    vector < vector <int> > all_metabolisms;
    //std::normal_distribution<double> distribution(3,1.0); // mean 2, variance 1
    
    for (int i = 0; i < pow(2,genome_length); i++){
      int pos = 0;
      int neg = 0;

      for (int j = 0; j < temp.size(); j++) {
	double a = 11*(2.0*drand48() - 1); //random number between (11,-11) actual range with floor / ceil is [10,-10]
	double b = 11*(2.0*drand48() - 1);
	if (a > 0) { a = floor(a); }
	else { a = ceil(a); }
	if (b > 0) { b = floor(b); }
	else { b = ceil(b); }
	temp[j] = a+b; 

	if (temp[j] > 0) { pos += 1; }    
	else if (temp[j] < 0) { neg += 1; } 
	else { temp[j] = 0; }
	
      }

      if ( pos > 2 && neg > 2) { cout << "we have a problem"; }  // bug check

      if (pos == 1) {   // we eat one nutrient
	for ( int k = 0; k < temp.size(); k++ ) {
	  if (temp[k] > 0) { temp[k] = 1; }
	}
      }

      if (neg == 1) {   // we excrete one nutrient
	for ( int k = 0; k < temp.size(); k++ ) {
	  if (temp[k] < 0) { temp[k] = -1; }

	}

      }
      
      if (pos == 2) {   // we eat two nutrients
	int loc1 = -1;   // find the location of the positive numbers in the vector
	int loc2 = -1;
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] > 0) { loc1 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] > 0 && k != loc1) { loc2 = k; }
	  
	}

	temp[loc1] = temp[loc1]/GCD(temp[loc1],temp[loc2]);
	temp[loc2] = temp[loc2]/GCD(temp[loc1],temp[loc2]);
      }

      if (pos == 3) {   // we eat 3 nutrients and must therefore excrete the other
	int loc1 = -1;
	int loc2 = -1;
	int loc3 = -1;
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] > 0) { loc1 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] > 0 && k != loc1) { loc2 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] > 0 && k != loc1 && k!= loc2) { loc3 = k; }
	  
	}
	
	int gcd_1 = GCD(temp[loc1],temp[loc2]);
	int gcd_2 = GCD(temp[loc1],temp[loc3]);

	  if (gcd_1 == gcd_2 && gcd_1 > 0) {
	    temp[loc1] = temp[loc1]/gcd_1;
	    temp[loc2] = temp[loc2]/gcd_1;
	    temp[loc3] = temp[loc3]/gcd_1;
	  }

      }

      if (neg == 1) {
	for ( int k = 0; k < temp.size(); k++ ) {
	  if (temp[k] < 0) { temp[k] = -1; }
	}
      }

      if (neg == 2) {
	int loc1 = -1;
	int loc2 = -1;
	
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] < 0) { loc1 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] < 0 && k != loc1) { loc2 = k; }
	  
	}
	temp[loc1] = temp[loc1]/GCD(-1*temp[loc1],-1*temp[loc2]);
	temp[loc2] = temp[loc2]/GCD(-1*temp[loc1],-1*temp[loc2]);

      }
      if (neg == 3) {

	int loc1 = -1; 
	int loc2 = -1;
	int loc3 = -1;
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] < 0) { loc1 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] < 0 && k != loc1) { loc2 = k; }
	  
	}
	for (int k = 0; k < temp.size(); k++) {
	  if (temp[k] < 0 && k != loc1 && k!= loc2) { loc3 = k; }
	  
	}
	
	int gcd_1 = GCD(-1*temp[loc1],-1*temp[loc2]);
	int gcd_2 = GCD(-1*temp[loc1],-1*temp[loc3]);


	  if (gcd_1 == gcd_2 && gcd_1 > 0) {
	    temp[loc1] = temp[loc1]/gcd_1;
	    temp[loc2] = temp[loc2]/gcd_1;
	    temp[loc3] = temp[loc3]/gcd_1;
	  }

      }
      if ( pos > 0 && neg > 0 && pos < 4 && neg < 4 ) {

	for (int m = 0; m < num_nutrients; m ++ ) {
	  
	  if ( temp[m] > 0 ) {   cons_ex_vector[m] = temp[m];     cons_ex_vector[m+num_nutrients] = 0; }
	  if ( temp[m] < 0 ) {   cons_ex_vector[m] = 0;           cons_ex_vector[m+num_nutrients] = -1*temp[m]; }
	  if ( temp[m] == 0 ) {  cons_ex_vector[m] = 0;           cons_ex_vector[m+num_nutrients] = 0; }

	}

      } else { // metabolism not viable - must both consume and excrete!

	for (int m = 0; m < 2*num_nutrients; m++) {
	  
	  cons_ex_vector[m] = 0;
	}

      }
      
        all_metabolisms.push_back(cons_ex_vector);
       
        
    }

    return all_metabolisms;
    
}

//************* GENOME ABIOTIC FACTORS INTERACTIONS ***********

//determines the abiotic impact due to metabolism for every genome

vector < double > abiotic_genome_interactions(int genome_length) {
    
    vector <double> abiotic_ints (pow(2,genome_length),0);
    
    for (int i = 0; i < pow(2,genome_length); i++){
        
      abiotic_ints[i] = 2.0*drand48() - 1.0; // now our range is between -1.0 and 1.0

    }
    
    return abiotic_ints;
    
}



//************* CALCULATE NUTRIENT DEMAND **********

vector <int> nutrient_demand_calc (vector <microbe> &species, double abiotic_T, double prefered_abiotic, int num_nutrients, double abiotic_scaling, vector< vector<int> > &n_g_interacts, int max_consumption) {

  
  double factor_i = abiotic_scaling*sqrt(pow(abiotic_T - prefered_abiotic, 2.0));
  
  double satisfaction = exp (-1.0*pow(factor_i,2.0));

  int total_demand;
  int i_demand;

  vector <int> nutrient_dem(num_nutrients,0);

  for (int j = 0; j < species.size(); j++) {
    total_demand = 0;

    for (int k = 0; k < num_nutrients; k++) {
      total_demand += n_g_interacts[species[j].genome][k];

    }

    if (total_demand <= floor(max_consumption * satisfaction) && total_demand > 0 ){

      i_demand = total_demand;
      
      while (total_demand + i_demand <= floor(max_consumption*satisfaction)){
	
	total_demand += i_demand;

      }
	
      while (total_demand > 0) {
	nutrient_dem[0] += n_g_interacts[species[j].genome][0] * species[j].population;
	nutrient_dem[1] += n_g_interacts[species[j].genome][1] * species[j].population;
	nutrient_dem[2] += n_g_interacts[species[j].genome][2] * species[j].population;
	nutrient_dem[3] += n_g_interacts[species[j].genome][3] * species[j].population;
	total_demand -= i_demand;
      }
      
    }
  }
  
  return nutrient_dem;
  
}


//*************** UPDATE ALL FLASKS **********************************//


int update_all_flasks (large_flask &main_flask, vector < flask > &flask_list, double abiotic_influx, int max_nutrient_inflow, int num_nutrients, double abiotic_start, double abiotic_end, double max_timesteps, double abiotic_env, int num_flasks, double mini_flask_inflow) {

      for (int j = 0; j < num_nutrients; j++ ){
                
	main_flask.nutrient[j] = main_flask.nutrient[j]*(1.0-abiotic_influx);
      }
            
      //nutrient inflow
            
      for (int j = 0; j < num_nutrients; j++){
  
	main_flask.nutrient[j] += max_nutrient_inflow;
            
      }

      double current_temperature_main = main_flask.temperature;
      abiotic_env -= (abiotic_start - abiotic_end)/max_timesteps;
	    
      main_flask.temperature = current_temperature_main*(1.0-abiotic_influx) + abiotic_env*abiotic_influx; // dilute with fresh inflow

      current_temperature_main = main_flask.temperature;
      double temperature_changes = 0.0;
      vector < double > current_nutrients_main = main_flask.nutrient;
      vector < double > current_nutrients_mini;
      double current_temperature_mini;

     for (int j = 0; j < num_flasks; j++) {

	current_temperature_mini = flask_list[j].temperature;
	flask_list[j].temperature = flask_list[j].temperature*(1.0-mini_flask_inflow) + current_temperature_main*mini_flask_inflow;
	temperature_changes += (current_temperature_mini*mini_flask_inflow + current_temperature_main*(1.0-mini_flask_inflow))/(num_flasks);

	for (int k = 0; k < num_nutrients; k++) {
		double nutrient_k_transfer = -1.0*mini_flask_inflow*flask_list[j].nutrient[k] + mini_flask_inflow*current_nutrients_main[k]/(num_flasks);
		flask_list[j].nutrient[k] += nutrient_k_transfer;
		main_flask.nutrient[k] -= nutrient_k_transfer;
	}
	

     } 

      main_flask.temperature = temperature_changes;
      return 0;

}

// ********************* MAIN CODE ***********************************//



int main(int argc, char **argv) {

    
    // MICROBE PARAMETERS
    int i;  // this is a marker for chosing and individual
    const int initial_population = 100;
    const int genome_length = 8;
    const int reproduction_thresh = 120; // biomass threashold to reproduce
    const int starve_thresh = 50; // if a microbe's biomass drops below this, it dies
    const int initial_biomass = 80; // starting biomass for every microbe
    const int max_consumption = 10; // the maximum rate at which micrcomes can comsume
    const double nutrient_conversion_eff = 0.6; // efficiency of microbe conversion
    const int maintainence_cost = 1; // how much biomass it costs per timestep to live
    const double p_mut = 0.01; // probability per gene of mutation in reproduction event
    const double p_kill = 0.002;  //probability of death per timestep from causes other than starvation
    const double prefered_abiotic = 150.0;
    
    
    // FLOW PARAMETERS
    const int num_nutrients = 4;
    const int num_abiotic = 1;  // if we change this we need to change other code...
    const int num_flasks = 4; // the number of flasks


    const double percentage_outflow = 0.2; // to calculate outflow
    const double max_nutrient_inflow = 150.0;
    int temp_nutrient_count;
    const double abiotic_start = 100.0;//190.0;
    const double abiotic_end = 100.0;
    double abiotic_env = abiotic_start;
    const double abiotic_influx = 0.2; // the percentage that temp will change each timestep
    const double mini_flask_inflow = 0.2;
    double abiotic_trickle = 0;
    const double main_flask_scale = 1.0; // this determines how large the main flask is with respect to 
					 // the smaller flasks. A value of 1.0 would mean that for each  
					 // mini flask there is a mini flasks worth of medium in the main flask
					 // Cannot be less than the value of mini_flask_inflow


    // RANDOM NUMBER GENERATORS

    int t = atoi(argv[1]); // seed for randomness
    srand48 (t);
    mt19937 rng(t);
    default_random_engine generator; // used for our distributions
    generator.seed(t);  //provide seed for randomness in different runs

    
    // VARIABLES FOR METABOLISM CALCULATIONS ETC
    int average_biomass;
    int i_biomass;
    int average_nutrient;
    vector <int> nutrient_demand(num_nutrients,0);    // Nutrient demand vector
    int demand_population;
    int count_eaten;
    double species_nutrient_avg;
    double species_biomass_avg;
    int nutrient_available;
    const double abiotic_scaling = 0.0; //denoted as tau in the paper
    double satisfaction;
    double factor_i;
    //double abiotic_T_sum = 0;
    
    
    // VARIABLES FOR KEEPING TRACK OF TIME
    int loop_variable;
    int timestep_length = initial_population; // num of iterations per timestep determined by total population
    int number_gens = 0;
    int timestep_counter = 0;
    int max_timesteps = 10*pow(10,4);
    int init_period = 5000;
    int init_counter = 0;
    int num_alive_flasks = 0;
    
    
    // MUTATION VARIABLES USED IN MUTATION EVENTS
    microbe temp_mutant;
    int did_we_mutate = 0;
    
    
    // DATA FILES NUMBER
    int file_num = atoi(argv[2]);


    // main flask set up
    vector < double > large_nutrient_init (num_nutrients, 0.0);
    large_flask main_flask;
    main_flask.temperature = 0.0;
    main_flask.mini_flasks = num_flasks;
    main_flask.nutrient = large_nutrient_init;

    // ENVIRONMENT SET UP -- NOW DONE LATER
    double abiotic_T = abiotic_start;
    
    vector<double> environment(num_nutrients, 0);
    vector<double> nutrient_trickle(num_nutrients,0);
    
    // METABOLISM SET UP
    vector < vector< int > >  n_g_interacts = nutrient_genome_interactions(genome_length,num_nutrients, generator); // initialises what metabolisms each genome codes for
    vector < double > a_g_interacts = abiotic_genome_interactions(genome_length); // how each genome affects the temperature


    // MINI FLASK SET UP
    flask new_flask;
    vector < flask > flask_list;

    for (int k = 0; k < num_flasks; k++) {  // set up each flask

	    microbe new_microbe;


	    vector < microbe > species_init; //(1, new_microbe);
	    int total_population = 0;
    
    
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

		      total_population++;
      
	    }

	    if (total_population != initial_population) {cout << "POPULATION  BUG" << endl;}

	new_flask.species = species_init;
	
	vector < double > nutrient_init (num_nutrients, 0.0);	
	
	new_flask.nutrient = nutrient_init; // start all flasks empty
	new_flask.temperature = 0.0; 	 // each flask starts at 0.0
	ostringstream filename;
	filename << "mini_flask_"+to_string(k)+"_run_"+to_string(file_num)+".txt";
	new_flask.file = filename.str();
        flask_list.push_back(new_flask);
	ofstream temp_file (new_flask.file);
	temp_file.close();
    }

    ofstream macro_data ("main_flask_macro_data_"+to_string(file_num)+".txt");
    ofstream nutrient_data ("main_flask_nutrient_data_"+to_string(file_num)+".txt");

    
    // BUG CHECKING VARIABLES
    int total_pop_bug = 0;


    while (init_counter < init_period) {   // INITIALISE OUR ENVIRONMENT

      /* ********************************************************************************
                                    NUTRIENT FLOW
      ********************************************************************************/

            
      // only update the flow once every time step       
      // nutrient outflow    
      /*for (int j = 0; j < num_nutrients; j++ ){ main_flask.nutrient[j] = main_flask.nutrient[j]*(1.0-percentage_outflow); }
            
      //nutrient inflow   
      for (int j = 0; j < num_nutrients; j++){ main_flask.nutrient[j] += max_nutrient_inflow; }

      double current_temperature_main = main_flask.temperature;
      abiotic_env -= (abiotic_start - abiotic_end)/max_timesteps;
	    
      main_flask.temperature = current_temperature_main*(1.0-abiotic_influx) + abiotic_env*abiotic_influx; // dilute with fresh inflow

      current_temperature_main = main_flask.temperature;
      double temperature_changes = 0.0;
      vector < double > current_nutrients_main = main_flask.nutrient;
      double current_temperature_mini;

     for (int j = 0; j < num_flasks; j++) {

	current_temperature_mini = flask_list[j].temperature;
	flask_list[j].temperature = flask_list[j].temperature*(1.0-mini_flask_inflow) + current_temperature_main*mini_flask_inflow;

	temperature_changes += (current_temperature_mini*mini_flask_inflow + current_temperature_main*(1.0-mini_flask_inflow))/(num_flasks);

	for (int k = 0; k < num_nutrients; k++) {
		double nutrient_k_transfer = -1.0*mini_flask_inflow*flask_list[i].nutrient[k] + mini_flask_inflow*current_nutrients_main[k]/(num_flasks);
		flask_list[i].nutrient[k] += nutrient_k_transfer;
		main_flask.nutrient[k] -= nutrient_k_transfer;
	}
	

     } 

      main_flask.temperature = temperature_changes; */

	update_all_flasks(main_flask, flask_list, abiotic_influx, max_nutrient_inflow, num_nutrients, abiotic_start, abiotic_end, max_timesteps, abiotic_env, num_flasks, mini_flask_inflow);
	init_counter++;

    }

    
    while (number_gens < max_timesteps) {
          
      
	/* ********************************************************************************
                                RECORD DATA
	********************************************************************************/

	// sorts species by order of population, largest first
	for (int f = 0; f < num_flasks; f++) { stable_sort (flask_list[f].species.begin(), flask_list[f].species.end(), sorting_rule);}

	int total_population = 0; // UPDATE THIS
	// Record data here!!!!!
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
	                       CALCULATE NUTRIENT TRICKLE
	**************************************************************************************/

	/*for (int j = 0; j < num_nutrients; j++ ){
	   nutrient_trickle[j] =  max_nutrient_inflow - environment[j]*percentage_outflow;
	   if (timestep_length > 0) { nutrient_trickle[j] = nutrient_trickle[j] / (1.0*timestep_length); }
	}

	abiotic_env -= (abiotic_start - abiotic_end)/max_timesteps;
	abiotic_trickle = (abiotic_T*(1.0-abiotic_influx) + abiotic_env*abiotic_influx) - abiotic_T;
	     
	if (timestep_length > 0) { abiotic_trickle = abiotic_trickle / (1.0*timestep_length); } */

	update_all_flasks(main_flask, flask_list, abiotic_influx, max_nutrient_inflow, num_nutrients, abiotic_start, abiotic_end, max_timesteps, abiotic_env, num_flasks, mini_flask_inflow);
	    
	//}

	/********************************************************************************
				LOOP OVER MINI FLASKS
	*********************************************************************************/

	for (int f = 0; f < num_flasks; f++) {

	   vector <microbe> species = flask_list[f].species;
	   vector <double> environment = flask_list[f].nutrient;
	   double abiotic_T = flask_list[f].temperature;

	   int iterations = 0;
	   for (int k = 0; k < species.size(); k++) { iterations += species[k].population; }
	     
	   for (int l = 0; l < iterations; l++) {
        
               /* ********************************************************************************
                                       NUTRIENT FLOW
                ********************************************************************************/

            
               // have a trickle every iteration adding up to the alloted count per timestep     
               // nutrient outflow

		for (int j = 0; j < num_nutrients; j++ ){
		   environment[j] += nutrient_trickle[j];
		   if (environment[j] < 0) { environment[j] = 0; }
		   // just in case the net nutrient influx is negative not possibl though I don't think...
		}

		//  double current_env = abiotic_T;
		abiotic_T += abiotic_trickle;
        
        	/* ********************************************************************************
                                	    KILL
        	 ********************************************************************************/
        
        	// NEED TO REMOVE BIOMASS WHEN AN INDIVIDUAL DIES
        	// USE SAME DISTRIBUTION TO CALCULATE BIOMASS TO REMOVE
        	// AS USING WHEN CALCUATING IF INDIVIDUALS REPRODUCE

        	i = chooseAgent(species, total_population); // choses a microbe at random

        	if (i > -1) {

	  	   // death event starvation
	  	   average_biomass = species[i].biomass/(1.0*species[i].population);
	  	   normal_distribution<double> biomass_dist( average_biomass, species_biomass_avg*0.01 ); // distribution of biomass in pop
	  	   i_biomass = floor(biomass_dist(generator));

	  	   if (i_biomass <= starve_thresh) {
            		// so if the biomass count is low, the above will be high so high prob of death
            		species[i].population--;
            		species[i].biomass -= i_biomass; // remove biomass of dead microbe 
            		total_population--;
            		if (species[i].biomass < 1) { species[i].biomass = 0; species[i].population = 0; }
            		if (species[i].population == 0){ species.erase(species.begin() + i); } // remove from list if extinct
	  	   } else if (drand48() <= p_kill && species[i].population > 0) {
            
			species[i].population--;
			species[i].biomass -= i_biomass; // remove biomass of dead microbe 
            		total_population--;
            		if (species[i].population == 0){ species.erase(species.begin() + i);} // remove from list if extinct
	  	   }
		}
        
        	/* ********************************************************************************
             	                   MAINTENANCE COST
        	 ********************************************************************************/

        	i = chooseAgent(species, total_population);
	
        	if (i > -1) { species[i].biomass--; }

        	/* ********************************************************************************
                         	       METABOLISM
         	********************************************************************************/
        
        	// metabolism event
        	i = chooseAgent(species, total_population);
	
		if( i > -1) {

		   factor_i = abiotic_scaling*sqrt(pow(abiotic_T - prefered_abiotic, 2.0));
	  	   satisfaction = exp (-1.0*pow(factor_i,2.0));

		   int total_count_eat = floor(max_consumption * satisfaction);
	  	   int minimum_count_eat = 0; // the minumum total number of nutrients microbe can intake
	  	   int max_count_eat = 0;
	  
	  	   for (int k = 0; k < num_nutrients; k++) {  minimum_count_eat += n_g_interacts[species[i].genome][k]; }

	  	   max_count_eat = minimum_count_eat;
	  	   int enough_nutrient = true;
	  
	  	   while ( max_count_eat <= total_count_eat && minimum_count_eat > 0 && total_count_eat > 0 ) {

	    		for (int k = 0; k < num_nutrients; k++) {
	      		   if (floor(environment[k]) < n_g_interacts[species[i].genome][k]) { enough_nutrient = false; }
	  		   // if there is not enough nutrient to satisfy metabolism in correct proportions, cannot eat
	    		}
	    		if (enough_nutrient == false) { break; } // exit if cannot eat
	    		for (int k = 0; k < num_nutrients; k++) {
	     		   environment[k] -= n_g_interacts[species[i].genome][k];
	     		   species[i].nutrient += n_g_interacts[species[i].genome][k];
	     		   if( environment[k] < 0) { cout << "NUTRIENT EATING PROBLEM" << endl; }
	    		}
	     		max_count_eat += minimum_count_eat; // can we have more than one eating event?
	  	   }
		}

        	/* ********************************************************************************
            	                    BIOMASS CREATION
        	 ********************************************************************************/

        	i = chooseAgent(species, total_population);
		if (i > -1){

	  	   species_nutrient_avg = 1.0*species[i].nutrient/species[i].population;

	  	   normal_distribution<double> nutrient_species_dist( species_nutrient_avg, species_nutrient_avg*0.1);
	  	   nutrient_available = floor(nutrient_species_dist(generator)); // we'll just round down as can't use half a biomass

	  	   while ( nutrient_available >= 5) { 
           
			species[i].nutrient -= 5;
			nutrient_available -= 5;
			species[i].biomass += int(5.0*nutrient_conversion_eff);
			species[i].waste += int(5*(1.0 - nutrient_conversion_eff));
			abiotic_T += 5.0*nutrient_conversion_eff*a_g_interacts[species[i].genome]; // the effect on abiotic environment ie pH has an effect per biomass created
		   }
		}

		/*********************************************************************************
                    	              WASTE 
	 	**********************************************************************************/

        	i = chooseAgent(species, total_population);
		if (i > -1){

	  	   double species_waste_avg = 1.0*species[i].waste/species[i].population;

	  	   normal_distribution<double> waste_species_dist( species_waste_avg, species_waste_avg*0.1);
	  	   int waste_available = round(waste_species_dist(generator));

	  	   if (waste_available > species[i].waste) { waste_available = species[i].waste; }
	  
	  	   int waste_count = 0;
	  	   for (int k = 0; k < num_nutrients; k++) {
	    		waste_count += n_g_interacts[species[i].genome][k+num_nutrients];
	  	   }
	  
	  	   while (species[i].waste >= waste_count && species[i].waste > 0) {

	    		for (int k = 0; k < num_nutrients; k++) {

	      		   species[i].waste -=  n_g_interacts[species[i].genome][k+num_nutrients];
	      		   environment[k] +=  n_g_interacts[species[i].genome][k+num_nutrients];
	    		}
	  	   }
		}
	
		/* ********************************************************************************
               		                 REPRODUCTION
         	********************************************************************************/

        	i = chooseAgent(species, total_population);
        	if (i > -1){
        
	  	   species_biomass_avg = (1.0*species[i].biomass)/species[i].population; // average biomass per indiviual
        
	  	   normal_distribution<double> biomass_species_dist( species_biomass_avg, species_biomass_avg*0.01 );
	  	   i_biomass = floor(biomass_species_dist(generator));
        
	  	   if (i_biomass >= reproduction_thresh) {
            		// DO WE MUTATE?
            		bitset<genome_length> mutant_genome(species[i].genome);
            		did_we_mutate = 0;
            
            		for (int j = 0; j < genome_length; j++){
	      		   if (drand48() <= p_mut){
				did_we_mutate = 1;
				if (mutant_genome[j] == 1) { mutant_genome[j] = 0; }
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
			   temp_mutant.genome = mutant_number;
			   temp_mutant.nutrient = 0; //  no nutrient count to begin with
			   temp_mutant.biomass = int(i_biomass / 2.0);  // half biomass goes to new mutant
			   species[i].biomass -= int(i_biomass / 2.0);
			   temp_mutant.population = 1; // initial population of 1
			   temp_mutant.waste = 0;
			   species.push_back(temp_mutant);
	      		}
	      		total_population++;
            		} else { // no mutation takes place, we add one to the population
	      		   species[i].population++;
	      		   total_population++;
            		}

        
	  	   }
        	} /* END OF MICROBE ACTIONS */

	   }
	   flask_list[f].species = species;
	   flask_list[f].nutrient = environment;
	   flask_list[f].temperature = abiotic_T;
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
