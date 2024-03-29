#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __internal_viral_dynamics__
#define __internal_viral_dynamics__
	
void internal_virus_model_setup( void );
void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt );

double interpolate( double input1, double input2, double param ); 

#endif 