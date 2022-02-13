#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __pharmacodynamics_submodel__
#define __pharmacodynamics_submodel__
	
// inputs: 

// outputs: 

extern Submodel_Information pharmacodynamics_submodel_info; 

void pharmacodynamics_submodel_setup( void );

void pharmacodynamics_response( Cell* pCell, Phenotype& phenotype, double dt ); 

void advance_internalization( Cell* pCell, Phenotype& phenotype , double dt );
void simple_effect_model( Cell* pCell, Phenotype& phenotype, double dt );
double Hill_function( double input, double Hill_power , double EC_50 );
double interpolate( double input1, double input2, double param );

void apply_therapies( void );

#endif 