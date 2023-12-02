#ifndef LOCALPROJECTIONS_H
#define LOCALPROJECTIONS_H

#include <string>

/**
 * Enumeration for the different field views.
 */

 namespace Proj 
 {

// TODO: OMG SO MUCH!!! 
enum Proj_Type {
    svd, 
    gauss_newton_L2_rank1,
    gauss_newton_L4_rank1,
    gauss_newton_L2_L4_rank1,
    
    gauss_newton_L2_rank2,
    gauss_newton_L4_rank2,
    gauss_newton_L2_L4_rank2,

    gauss_newton_L2_rank3,
    gauss_newton_L4_rank3,
    gauss_newton_L2_L4_rank3,

// TODO later

    gauss_newton_L2_rankN,
    gauss_newton_L4_rankN,
    gauss_newton_L2_L4_rankN,

    Element_COUNT
};

/**
 * Converts a field view enumeration to a string.
 * 
 * @param view The field view to convert.
 * @return The name of the field view as a string.
 */
// std::string fieldViewToString(Field_View view);
// std::string fieldViewToFileStub(Field_View view);


 };

#endif // LOCALPROJECTIONS_H
