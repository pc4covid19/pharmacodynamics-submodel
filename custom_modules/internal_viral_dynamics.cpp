#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

std::string internal_virus_replication_version = "0.3.0"; 

Submodel_Information internal_viral_dynamics_info; 


void internal_virus_model_setup( void )
{
		// set version
	internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
	internal_viral_dynamics_info.version = internal_virus_replication_version; 
		// set functions 
	internal_viral_dynamics_info.main_function = NULL; 
	internal_viral_dynamics_info.phenotype_function = internal_virus_model; 
	internal_viral_dynamics_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled_virion" ); 

	internal_viral_dynamics_info.cell_variables.push_back( "virion_uncoating_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_to_RNA_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "protein_synthesis_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_assembly_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_export_rate" ); 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	
	return; 
}

void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index("virion");
	static int nA_external = microenvironment.find_density_index("assembled virion");

	static int nV_internal = pCell->custom_data.find_variable_index("virion");
	static int nA_internal = pCell->custom_data.find_variable_index("assembled_virion");

	static int nUV = pCell->custom_data.find_variable_index("uncoated_virion");
	static int nR = pCell->custom_data.find_variable_index("viral_RNA");
	static int nP = pCell->custom_data.find_variable_index("viral_protein");

	double sVol = 2494;
	double Vol = pCell->phenotype.volume.total;
	if (phenotype.death.dead == true)
	{
		return;
	}

	static int n_fusion = pCell->custom_data.find_variable_index("cell_fusion_number"); 

	pCell->custom_data[nA_internal] = phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 
	double x[4][5] = { {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0} };//initialize x
	double f[4][5] = { {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0} };//initialize f
	int j;//initialize counter

	//initial values for RK4
	x[0][0] = pCell->custom_data[nV_internal];
	x[0][0] = x[0][0] / Vol;
	x[0][1] = pCell->custom_data[nUV];
	x[0][1] = x[0][1] / Vol;
	x[0][2] = pCell->custom_data[nR];
	x[0][2] = x[0][2] / Vol;
	x[0][3] = pCell->custom_data[nP];
	x[0][3] = x[0][3] / Vol;
	x[0][4] = pCell->custom_data[nA_internal];
	x[0][4] = x[0][4] / Vol;

	double hRNA = pCell->custom_data["RNA_replication_half"] / sVol;

	// update the feedback with including the cell fusion problem, YW 2022
    // add negative feedback for replication
	
	// testing change to a MM like inhibition
	
	double kmod = pCell->custom_data["uncoated_to_RNA_rate"] * (1-pCell->custom_data[nR]/((pCell->custom_data[n_fusion] +1) *parameters.doubles("RNA_threshold") + pCell->custom_data[nR]));
	//possibly put this in the loop

	for (j = 0; j < 4; j++) {
		f[j][0] = { -pCell->custom_data["virion_uncoating_rate"] * x[j][0] }; // uncoat endocytosed virus
		f[j][1] = { pCell->custom_data["virion_uncoating_rate"] * x[j][0] - pCell->custom_data["uncoated_to_RNA_rate"] * x[j][1] }; // convert uncoated virus to usable mRNA 
		f[j][2] = { pCell->custom_data["uncoated_to_RNA_rate"] * x[j][1] + pCell->custom_data["max_RNA_replication_rate"] / sVol * x[j][2] / (x[j][2] + hRNA) - pCell->custom_data["basal_RNA_degradation_rate"] * x[j][2] }; // RNA replication post uncoated to RNA calc
		f[j][3] = { pCell->custom_data["protein_synthesis_rate"] * x[j][2] - pCell->custom_data["virion_assembly_rate"] * x[j][3] }; // use mRNA to create viral protein 
		f[j][4] = { pCell->custom_data["virion_assembly_rate"] * x[j][3] - pCell->custom_data["virion_export_rate"] * x[j][4] }; // assemble virus
		if (j == 0 || j == 1) {
			x[j + 1][0] = x[0][0] + dt / 2 * f[j][0]; //first and second x approximations
			x[j + 1][1] = x[0][1] + dt / 2 * f[j][1]; //first and second x approximations
			x[j + 1][2] = x[0][2] + dt / 2 * f[j][2]; //first and second x approximations
			x[j + 1][3] = x[0][3] + dt / 2 * f[j][3]; //first and second x approximations
			x[j + 1][4] = x[0][4] + dt / 2 * f[j][4]; //first and second x approximations
		}
		if (j == 2) {
			x[j + 1][0] = x[0][0] + dt * f[j][0]; //third approximation
			x[j + 1][1] = x[0][1] + dt * f[j][1]; //third approximation
			x[j + 1][2] = x[0][2] + dt * f[j][2]; //third approximation
			x[j + 1][3] = x[0][3] + dt * f[j][3]; //third approximation
			x[j + 1][4] = x[0][4] + dt * f[j][4]; //third approximation
		}
	}

	pCell->custom_data[nV_internal] = x[0][0] + dt * (f[0][0] / 6 + f[1][0] / 3 + f[2][0] / 3 + f[3][0] / 6);
	pCell->custom_data[nV_internal] = pCell->custom_data[nV_internal] * Vol;

	pCell->custom_data[nUV] = x[0][1] + dt * (f[0][1] / 6 + f[1][1] / 3 + f[2][1] / 3 + f[3][1] / 6);
	pCell->custom_data[nUV] = pCell->custom_data[nUV] * Vol;

	pCell->custom_data[nR] = x[0][2] + dt * (f[0][2] / 6 + f[1][2] / 3 + f[2][2] / 3 + f[3][2] / 6); //detirmine n+1
	pCell->custom_data[nR] = pCell->custom_data[nR] * Vol;

	pCell->custom_data[nP] = x[0][3] + dt * (f[0][3] / 6 + f[1][3] / 3 + f[2][3] / 3 + f[3][3] / 6); //detirmine n+1
	pCell->custom_data[nP] = pCell->custom_data[nP] * Vol;

	pCell->custom_data[nA_internal] = x[0][4] + dt * (f[0][4] / 6 + f[1][4] / 3 + f[2][4] / 3 + f[3][4] / 6); //detirmine n+1
	pCell->custom_data[nA_internal] = pCell->custom_data[nA_internal] * Vol;
 
	phenotype.secretion.net_export_rates[nA_external] =
		pCell->custom_data["virion_export_rate"] * pCell->custom_data[nA_internal];

	// copy data from custom variables to "internalized variables" 

	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];	
	
	return; 
}
	
