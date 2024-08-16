//utility functions for the MultiPhysicsBVP class

#include "../../include/multiPhysicsBVP.h"

//return index of given field name if exists, else throw error
template <int dim, int degree>
unsigned int MultiPhysicsBVP<dim,degree>::getFieldIndex(std::string _name) {
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
     if (it->name.compare(_name)==0) return it->index;
   }
   pcout << "\nutilities.h: field '" << _name.c_str() << "' not initialized\n";
   exit(-1);
}

#include "../../include/multiPhysicsBVP_template_instantiations.h"

