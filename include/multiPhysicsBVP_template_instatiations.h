/*
 * matrixFreePDE_template_instantiations.h
 *
 *  Created on: August 1, 2024
 *      Author: david_montiel_t
 */

#ifndef INCLUDE_MULTIPHYSICSBVP_TEMPLATE_INSTANTIATIONS_H_
#define INCLUDE_MULTIPHYSICSBVP_TEMPLATE_INSTANTIATIONS_H_

#ifndef MULTIPHYSICSBVP_TEMPLATE_INSTANTIATION
#define MULTIPHYSICSBVP_TEMPLATE_INSTANTIATION
template class MultiPhysicsBVP<2,1>;
template class MultiPhysicsBVP<3,1>;
template class MultiPhysicsBVP<2,2>;
template class MultiPhysicsBVP<3,2>;
template class MultiPhysicsBVP<3,3>;
template class MultiPhysicsBVP<2,3>;
template class MultiPhysicsBVP<3,4>;
template class MultiPhysicsBVP<2,4>;
#endif

#endif /* INCLUDE_MULTIPHYSICSBVP_TEMPLATE_INSTANTIATIONS_H_ */
