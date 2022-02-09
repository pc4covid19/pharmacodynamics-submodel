#include "./epithelium_submodel.h" 

using namespace PhysiCell; 

std::string epithelium_submodel_version = "0.1.0"; 

Submodel_Information epithelium_submodel_info; 

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions 
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );
	
	return; 
}

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// immune stuff, not need for PBPD model
	// static int debris_index = microenvironment.find_density_index( "debris");
	
	// receptor dynamics 
	// requires faster time scale - done in main function 
	
	// viral dynamics model 
	internal_viral_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_model(pCell,phenotype,dt);
	
	// viral response model 
	internal_virus_response_model_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_response_model(pCell,phenotype,dt);	

	// PBPD submodel
	pharmacodynamics_submodel_info.phenotype_function(pCell,phenotype,dt); 

	// cell fusion model
	epithelium_fusion(pCell,phenotype,dt); 
	

 /*  // immune stuff, not need for PBPD model
	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt ); 
	
	// if I am dead, remove all adhesions 
	static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
	standard_elastic_contact_functionif( phenotype.death.dead == true )
	{
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
	}
*/
	
/*
	// cell secretion belongs in viral response 
	
	// if I am dead, make sure to still secrete the chemokine 
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	static int nP = pCell->custom_data.find_variable_index( "viral_protein"); 
	double P = pCell->custom_data[nP];
	
	static int nAV = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	double AV = pCell->custom_data[nAV]; 

	static int nR = pCell->custom_data.find_variable_index( "viral_RNA");
	double R = pCell->custom_data[nR];
	
	if( R >= 1.00 - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = AV; // P; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0; 
	}
*/	
	
	// if I am dead, don't bother executing this function again 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}



void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	
	pCell->is_movable = false; 
	
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		// remove_all_adhesions( pCell ); 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		pCell->functions.contact_function = NULL; 

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		return; 
	}	
	
	// this is now part of contact_function 
	/*
	// if I'm adhered to something ... 
	if( pCell->state.neighbors.size() > 0 )
	{
		// add the elastic forces 
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
	}
	*/
	return; 
}


// simulate cell fusion, YW 2022 !!!
void epithelium_fusion( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	static int drug_effect_intra = pCell->custom_data.find_variable_index("Intracellular_drug_effect"); 
	static int time_ingest_index = pCell->custom_data.find_variable_index("time_to_next_phagocytosis"); 
	static int ingest_rate_index = pCell->custom_data.find_variable_index("phagocytosis_rate"); 
	static int n_fusion = pCell->custom_data.find_variable_index("cell_fusion_number"); 

	static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_ACE2" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_ACE2" ); 
	static int nR_IU = pCell->custom_data.find_variable_index( "unbound_internal_ACE2" ); 
	static int nR_IB = pCell->custom_data.find_variable_index( "bound_internal_ACE2" ); 

	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	static int nUV = pCell->custom_data.find_variable_index( "uncoated_virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral_RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral_protein" ); 

	static bool drug_cell_fusion_enabled = parameters.bools( "drug_cell_fusion" );  

	// follow COVID19 model: (Adrianne) check if still phagocytosing something, added if statement to say that if cell is still 
	// internalising current material not to phagocytose anything else

	if( pCell->custom_data[ time_ingest_index ] > PhysiCell_globals.current_time )
	{ return; }	

	// https://github.com/MathCancer/PhysiCell/blob/81e85ce1dcb7156864a4460666e60878b82d3e97/sample_projects/virus_macrophage/custom_modules/custom.cpp

	Cell* pTestCell = NULL; 
	std::vector<Cell*> neighbors = pCell->cells_in_my_container(); 
	
	for( int n=0; n < neighbors.size() ; n++ )
	{
		pTestCell = neighbors[n]; 
		// if it is not me and not a macrophage 
		if( pTestCell != pCell && pTestCell->phenotype.death.dead == false )
		{
			// calculate distance to the cell 
			std::vector<double> displacement = pTestCell->position;
			displacement -= pCell->position;
			double distance = norm( displacement ); 
			
			double max_distance = pCell->phenotype.geometry.radius + pTestCell->phenotype.geometry.radius; 
			max_distance *= 1.5;  // 1.1

			// caculate the mean location
			std::vector<double> mean_position = 0.5*(pTestCell->position + pCell->position);

			// assume fusion probabbility is positive related with intra-viral load
			double assembled_virion_pop = pTestCell->phenotype.molecular.internalized_total_substrates[nA_external]; 
			double prob_cell_fusion = Hill_function(assembled_virion_pop, 1, parameters.doubles("half_max_cell_fusion"))
					*pCell->custom_data[ ingest_rate_index ]*dt;  

			// caculate the probability of ingest based on drug effect (the larger it is, the smaller probability of cell fusion)
			// nonliear interplotation, assume r_min = 0
			// prob = (1-E)^n*r_0 *dt + E^n*r_min *dt = (1-E)^n*r_0 *dt 
			double temp = 1 - pCell->custom_data[drug_effect_intra];
			if( drug_cell_fusion_enabled )
			{
				prob_cell_fusion *= std::pow(temp, parameters.doubles("drug_fusion_power"));  // 3, 5
			}
			
			// test for viral load, distance, probability and cell fusion times
			if( assembled_virion_pop > 10 && distance < max_distance && 
				   UniformRandom() < prob_cell_fusion &&  pCell->custom_data[ n_fusion ] < 16 )
			{
				// std::cout << "\t\tnom nom nom" << std::endl; 

				// what does the ingest_cell mean in core function 
				// https://github.com/MathCancer/PhysiCell/blob/d45d290dd6e956bcc4c6899ff193748aad54bab4/core/PhysiCell_cell.cpp#L1163

				pCell->ingest_cell( pTestCell );

				pCell->position = mean_position;      // update the new position as of average 
				pCell->update_voxel_in_container();   // update in the data structure 

				// cell fusion number adding one and TestCell, as TestCell may have already fused !!!!!!!
				pCell->custom_data[ n_fusion ] += ( 1+ pTestCell->custom_data[n_fusion] ); 

				// update ACE2 stuff
				pCell->custom_data[ nR_EU ] += pTestCell->custom_data[ nR_EU ]; 
				pCell->custom_data[ nR_EB ] += pTestCell->custom_data[ nR_EB ]; 
				pCell->custom_data[ nR_IU] += pTestCell->custom_data[ nR_IU ]; 
				pCell->custom_data[ nR_IB ] += pTestCell->custom_data[ nR_IB ]; 

				// update intracellular viral stuff
				pCell->custom_data[ nV_internal ] += pTestCell->custom_data[ nV_internal ]; 
				pCell->custom_data[ nA_internal ] += pTestCell->custom_data[ nA_internal ]; 
				pCell->custom_data[ nUV ] += pTestCell->custom_data[ nUV ]; 
				pCell->custom_data[ nR ] += pTestCell->custom_data[ nR ]; 
				pCell->custom_data[ nP ] += pTestCell->custom_data[ nP ]; 

				// if cell fusion, virion replicate faster
				// pCell->custom_data["uncoated_to_RNA_rate"] += parameters.doubles("uncoated_to_RNA_rate_original");
				// pCell->custom_data["protein_synthesis_rate"] += parameters.doubles("protein_synthesis_rate_original");
				pCell->custom_data["uncoated_to_RNA_rate"] += pTestCell->custom_data["uncoated_to_RNA_rate"]; 
				pCell->custom_data["protein_synthesis_rate"] += pTestCell->custom_data["protein_synthesis_rate"]; 

				// if cell fusion, cell can tolerate more virus
				// pCell->custom_data["max_apoptosis_half_max"] += parameters.doubles("max_apoptosis_half_max_orginal"); 
				pCell->custom_data["max_apoptosis_half_max"] += pTestCell->custom_data["max_apoptosis_half_max"];

				double time_to_ingest = pTestCell->phenotype.volume.total / pCell->custom_data["material_internalisation_rate"];
				pCell->custom_data[ time_ingest_index ] = PhysiCell_globals.current_time + time_to_ingest;	
				
			}
		}
	}

	return; 
}



void epithelium_submodel_setup( void )
{
	Cell_Definition* pCD;
	
	// set up any submodels you need 
	// viral replication 


    // PBPD submodel
    // pharmacodynamics_submodel_setup(); 


	// receptor trafficking 
	receptor_dynamics_model_setup(); // done 
	// viral replication 
	internal_virus_model_setup();	
	// single-cell response 
	internal_virus_response_model_setup(); 
 	
	// set up epithelial cells
		// set version info 
	epithelium_submodel_info.name = "epithelium model"; 
	epithelium_submodel_info.version = epithelium_submodel_version; 
		// set functions 
	epithelium_submodel_info.main_function = NULL; 
	epithelium_submodel_info.phenotype_function = epithelium_phenotype; 

	// immune stuff, not need for PBPD model
	// epithelium_submodel_info.mechanics_function = epithelium_mechanics; 
	epithelium_submodel_info.mechanics_function = NULL;
	
		// what microenvironment variables do you expect? 
	epithelium_submodel_info.microenvironment_variables.push_back( "virion" ); 

    // immune stuff, not need for PBPD model
	// epithelium_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	// epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	// epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
		// what custom data do I need? 
	//epithelium_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	epithelium_submodel_info.register_model();	

	
	// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "lung epithelium" ); 
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;

	// immune stuff, not need for PBPD model
	// pCD->functions.contact_function = epithelium_contact_function; 
	pCD->functions.contact_function = NULL; 
	
	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	
	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		#pragma omp critical
		{
		std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl; 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

