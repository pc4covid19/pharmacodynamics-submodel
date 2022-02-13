#include "./pharmacodynamics_submodel.h" 

using namespace PhysiCell; 

std::string pharmacodynamics_submodel_version = "0.1.0"; 

Submodel_Information pharmacodynamics_submodel_info; 

void pharmacodynamics_submodel_setup( void )
{
	// set up the model 
		// set version info 
	pharmacodynamics_submodel_info.name = "pharmacodynamics on virus"; 
	pharmacodynamics_submodel_info.version = pharmacodynamics_submodel_version; 
		// set functions 
	pharmacodynamics_submodel_info.main_function = NULL; 
	pharmacodynamics_submodel_info.phenotype_function = pharmacodynamics_response; 
	pharmacodynamics_submodel_info.mechanics_function = NULL; 
	
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	
	pharmacodynamics_submodel_info.cell_variables.push_back( "Uptaken_drug" ); 
	pharmacodynamics_submodel_info.cell_variables.push_back( "Intracellular_drug_concentration" ); 
	pharmacodynamics_submodel_info.cell_variables.push_back( "Intracellular_drug_effect" ); 
	pharmacodynamics_submodel_info.cell_variables.push_back( "Extracellular_drug_effect" ); 

	pharmacodynamics_submodel_info.microenvironment_variables.push_back( "drug" ); 
	
		// register the submodel  
	pharmacodynamics_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
		
//	pCD = find_cell_definition( "lung epithelium" ); 
//	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
//	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	
	return; 
}


void advance_internalization( Cell* pCell, Phenotype& phenotype , double dt )
{
	
	// static double drug_saturation = parameters.doubles( "drug_saturation" );  // 6000
	static double internalization_rate = parameters.doubles( "internalization_rate" );  // r = 0.0058
	
	static int cell_drug_i = pCell->custom_data.find_variable_index("Uptaken_drug");
	static int intra_drug_concentration_i = pCell->custom_data.find_variable_index("Intracellular_drug_concentration");
       static int drug_substrate_ID = microenvironment.find_density_index( "drug" );
	
	// obtain the total uptaken drug
	pCell->custom_data[ cell_drug_i ] = pCell->phenotype.molecular.internalized_total_substrates[ drug_substrate_ID ]; 

    // set the uptake rate based on the endocytosis drug
    /*
	double temp = pCell->custom_data[ cell_drug_i ]; // n(i)
	temp /= drug_saturation;  // n(i)/n*
	phenotype.secretion.uptake_rates[drug_substrate_ID] = 1.0; // 1.0 
	phenotype.secretion.uptake_rates[drug_substrate_ID] -= temp; // 1 - n(i)/n* 	
	phenotype.secretion.uptake_rates[drug_substrate_ID] *= internalization_rate; // r*( 1- n(i)/n* )
	*/

	// In this project, we use passive drug internalization dynamics
	double extra_drug_concentration = pCell->nearest_density_vector()[drug_substrate_ID];

	// net_export_rate = - r*( rho*Vcell - n(i) )
	// "-" is because of actually "uptake"
	phenotype.secretion.net_export_rates[drug_substrate_ID] = 
	   - internalization_rate * ( extra_drug_concentration * phenotype.volume.total - pCell->custom_data[ cell_drug_i ] );


	// the concentration of drug inside one cell
	pCell->custom_data[ intra_drug_concentration_i ] = pCell->custom_data[ cell_drug_i ]; // C
	pCell->custom_data[ intra_drug_concentration_i ] /= phenotype.volume.total; // C/Vcell
	
	// update AUC 
	// pCell->custom_data[ Drug_concentration_AUC_i ] += dt* pCell->custom_data[ cell_Drug_concentraion_i ];
	
	return; 

}


// simple drug concentration model
void simple_effect_model( Cell* pCell, Phenotype& phenotype, double dt )
{	
	static int intra_drug_concentration_i = pCell->custom_data.find_variable_index("Intracellular_drug_concentration");    
	static int intra_drug_effect_i = pCell->custom_data.find_variable_index("Intracellular_drug_effect"); 
	static int extra_drug_effect_i = pCell->custom_data.find_variable_index("Extracellular_drug_effect"); 
	static int drug_substrate_ID = microenvironment.find_density_index( "drug" );
	static double EC_50 = parameters.doubles( "EC_50" ); // 1.2
	static double Hill_power =  parameters.doubles( "Hill_power" ); // 4
	static double S = parameters.doubles( "drug_efficiency" );  // 100%
	
	pCell->custom_data[ intra_drug_effect_i ] = Hill_function( pCell->custom_data[ intra_drug_concentration_i ] , Hill_power , EC_50 ); 
	pCell->custom_data[ intra_drug_effect_i ] *= S;

	double extra_drug_concentration = pCell->nearest_density_vector()[drug_substrate_ID];
	pCell->custom_data[ extra_drug_effect_i ] = Hill_function( extra_drug_concentration , Hill_power , EC_50 ); 
	pCell->custom_data[ extra_drug_effect_i ] *= S;


	return;
}


double Hill_function( double input, double Hill_power , double EC_50 )
{
	double temp = input; // c
	temp /= EC_50; // c/c0
	temp = std::pow( temp, Hill_power ); // (c/c0)^n 
	double output = temp; // (c/c0)^n
	temp += 1.0; // 1 + (c/c0)^n 
	output /= temp; // // (c/c0)^n / ( 1 + (c/c0)^n ) 
	
	return output; 	
}


// output = (1-param)*input1 + param*input2 ( 0 <= param <= 1 )
double interpolate( double input1, double input2, double param )
{
	 double Output = input2; // input2
	 Output -= input1; // input2-input1
	 Output *= param; // param*(input2-input1)
	 Output += input1; // input1 + param*(input2-input1) 
		// == (1-param)*input1 + param*input2 
		
	 return Output; 
}


void pharmacodynamics_response( Cell* pCell, Phenotype& phenotype, double dt )
{
    // drug effect !!!!!
	advance_internalization( pCell, phenotype , dt ); 
       simple_effect_model( pCell, phenotype , dt );

	static int drug_effect_extra = pCell->custom_data.find_variable_index("Extracellular_drug_effect"); 
	static bool endo_export_enabled = parameters.bools( "drug_endo_export" );  
	static int drug_effect_intra = pCell->custom_data.find_variable_index("Intracellular_drug_effect"); 
	static bool replication_enabled = parameters.bools( "drug_replication" );  

	static int RNA_index =  pCell->custom_data.find_variable_index( "viral_RNA" ); 
	
    // uptake new production rate based on drug effect 
	if ( endo_export_enabled )
	{
	   if( pCell->custom_data[RNA_index] >= parameters.doubles("RNA_threshold") )
		{
                 pCell->custom_data[ "ACE2_endocytosis_rate" ] = interpolate( 
        	   parameters.doubles("ACE2_endocytosis_rate_feedback"),  
                 parameters.doubles("drug_ACE2_endocytosis_rate"), 
                 pCell->custom_data[ drug_effect_intra]); 
		}
	   else 
	    {
	    	pCell->custom_data[ "ACE2_endocytosis_rate" ] = interpolate( 
        	   parameters.doubles("ACE2_endocytosis_rate_original"),  
               parameters.doubles("drug_ACE2_endocytosis_rate"), 
               pCell->custom_data[ drug_effect_intra]); 
	    }

		pCell->custom_data[ "virion_export_rate"] = interpolate( 
               parameters.doubles("virion_export_rate_original"),  
               parameters.doubles("drug_virion_export"), 
               pCell->custom_data[ drug_effect_intra]);   
    }

    // uptake new production rate based on drug effect
	if ( replication_enabled )
	{

	// Maraviroc doesn't impact viral binding process accoding to Aarith experiental results
	// Let's only focus on uncoating RNA and protein synthesis rate 
	/*	
		if( pCell->custom_data[RNA_index] >= parameters.doubles("RNA_threshold") )
		{

	          pCell->custom_data[ "ACE2_binding_rate" ] = interpolate( 
	 	    	parameters.doubles("ACE2_binding_rate_feedback"),  
	        	parameters.doubles("drug_ACE2_binding_rate"), 
	        	pCell->custom_data[ drug_effect_intra]); 
		}
		else
		{
		   pCell->custom_data[ "ACE2_binding_rate" ] = interpolate( 
	 	   	parameters.doubles("ACE2_binding_rate_original"),  
	        	parameters.doubles("drug_ACE2_binding_rate"), 
	        	pCell->custom_data[ drug_effect_intra]); 
		}
 	*/

	/*
        pCell->custom_data[ "virion_uncoating_rate"] = interpolate( 
               parameters.doubles("virion_uncoating_rate_original"),  
               parameters.doubles("drug_virion_uncoated"), 
               pCell->custom_data[ drug_effect_intra]);   
	*/

	if( pCell->custom_data[RNA_index] >= parameters.doubles("RNA_threshold") )
	{
		pCell->custom_data[ "uncoated_to_RNA_rate" ] = interpolate( 
        		parameters.doubles("uncoated_to_RNA_rate_feedback"),  
               	parameters.doubles("drug_uncoated_to_RNA_rate"), 
               	pCell->custom_data[ drug_effect_intra]);
	}
	else
	{
		pCell->custom_data[ "uncoated_to_RNA_rate" ] = interpolate( 
        	   	parameters.doubles("uncoated_to_RNA_rate_original"),  
                 	parameters.doubles("drug_uncoated_to_RNA_rate"), 
                 	pCell->custom_data[ drug_effect_intra]);
	}

        pCell->custom_data[ "protein_synthesis_rate" ] = interpolate( 
        	   parameters.doubles("protein_synthesis_rate_original"),  
               parameters.doubles("drug_protein_synthesis_rate"), 
               pCell->custom_data[ drug_effect_intra]); 

        /*
        pCell->custom_data[ "virion_assembly_rate" ] = interpolate( 
        	   parameters.doubles("virion_assembly_rate_original"),  
               parameters.doubles("drug_virion_assembly_rate"), 
               pCell->custom_data[ drug_effect_intra]);  
        */ 
    }  
}


// change drug does, inject delay and injection duration
void apply_therapies( void )  
{
	
	static double doses = parameters.doubles( "drug_doses" ); // 5.0

	static double therapy_dt = 3; // update therapy every 3 minutes
	static double next_therapy_time = 0.0; 
	
	static double tolerance = 0.01 * diffusion_dt; 
	
	static int therapy_update_time_multiplier = 0; 

	// DON'T CHECK FOR THERAPY EVERY 0.01 MINUTES! 
	// Let's check it every 3 minutes instead
	
	if( PhysiCell_globals.current_time <= next_therapy_time - tolerance )
	{ return; }
	else
	{ 
		therapy_update_time_multiplier++;
		next_therapy_time = therapy_update_time_multiplier*therapy_dt;     
	}

	static int dose_dt = parameters.ints( "dose_duration" ); // 24*60   inject drug every 24 hours
	static int start_inject = parameters.ints( "start_inject" );  // 0  the first time injecting drug
	static int next_dose_time = 0;
	static int dose_update = 0;
	static int index = microenvironment.find_density_index( "drug" );  

	if ( PhysiCell_globals.current_time <= start_inject - tolerance )
	{ return; }
    else
    {   
	    if ( PhysiCell_globals.current_time > next_dose_time - tolerance )
		{

			#pragma omp parallel for 
			for( int i=0 ; i < microenvironment.mesh.voxels.size() ;i++ )
			{
				// if( microenvironment.mesh.voxels[i].is_Dirichlet == true )
				// { microenvironment.update_dirichlet_node( i, index, dose ); }
				// microenvironment(i)[index] += doses;
				microenvironment.set_substrate_dirichlet_activation(index, i, true); 
				microenvironment.update_dirichlet_node( i, index, doses ); 

			}	
		
			dose_update++;
			next_dose_time = start_inject + dose_update*dose_dt; 
		}

    }
		
	
	return; 

}






























