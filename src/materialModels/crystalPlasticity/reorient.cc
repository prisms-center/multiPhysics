#include "../../../include/crystalPlasticity.h"
#include <deal.II/base/symmetric_tensor.h>

template <int dim>
void crystalPlasticity<dim>::reorient() {
    //Update the history variables

    SymmetricTensor<2, dim, double> C_old;
    SymmetricTensor<2, dim, double> C_new;

    Vector<double> eigenvals(dim);
    FullMatrix<double> C_old_temp(dim,dim),C_new_temp(dim,dim),Fe_old(dim,dim), Fe_new(dim,dim), eigenvecs(dim,dim),Lambda(dim,dim),U_old(dim,dim),U_new(dim,dim),R_old(dim,dim),R_new(dim,dim),Omega(dim,dim),temp(dim,dim);
    Lambda=IdentityMatrix(dim);
    Omega=0.0;
    Vector<double> rot1(dim),Omega_vec(dim),rold(dim),dr(dim),rnew(dim);
    FullMatrix<double> rotmat(dim,dim);
    //int itgno;
    unsigned int num_local_cells = this->triangulation_cp.n_locally_owned_active_cells();

    for (unsigned int i=0; i<num_local_cells; ++i) {
        for(unsigned int j=0;j<N_qpts;j++){

            C_old_temp=0.0;
            C_old.clear();

            Fe_old=Fe_conv[i][j];
            Fe_new=Fe_iter[i][j];
            
            /////In the case of twin reorientation, do not update the orientation for one loading increment
            if (twin_conv[i][j]!=twin_iter[i][j]){
              Fe_old=Fe_new;
            }
            ///////////////////////////////////////////////////////////////////////////////////
            
            Fe_old.Tmmult(C_old_temp,Fe_old);

            // Copy only the symmetric part of C_old_temp into C_old
            for (unsigned int i = 0; i < dim; ++i) {
                for (unsigned int j = i; j < dim; ++j) {
                    C_old[i][j] = C_old_temp(i, j);
                    if (i != j)
                        C_old[j][i] = C_old_temp(i, j);  // Ensure symmetry
                }
            }         

            //Compute eigenvalues and eigenvectors for C_old
            auto eigenpairs_old = eigenvectors(C_old);

            // Extract eigenvalues and eigenvectors
            for (unsigned int i = 0; i < dim; ++i) {
                eigenvals[i] = eigenpairs_old[i].first;  // Store eigenvalues
                for (unsigned int j = 0; j < dim; ++j) {
                        eigenvecs(j, i) = eigenpairs_old[i].second[j];  // Store eigenvectors column-wise
                }
            }

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvals(k));
            }

            eigenvecs.mmult(U_old,Lambda);
            temp=U_old; temp.mTmult(U_old,eigenvecs);
            R_old.invert(U_old);
            temp=R_old; Fe_old.mmult(R_old,temp);
            Fe_new.Tmmult(C_new_temp,Fe_new);

            C_new.clear();
            // Copy only the symmetric part of C_new_temp into C_new
            for (unsigned int i = 0; i < dim; ++i) {
                for (unsigned int j = i; j < dim; ++j) {
                    C_new[i][j] = C_new_temp(i, j);
                    if (i != j)
                        C_new[j][i] = C_new_temp(i, j);  // Ensure symmetry
                }
            }
            //Compute eigenvalues and eigenvectors for C_new
            auto eigenpairs_new = eigenvectors(C_new);

            // Extract eigenvalues and eigenvectors
            for (unsigned int i = 0; i < dim; ++i) {
                eigenvals[i] = eigenpairs_new[i].first;  // Store eigenvalues
                for (unsigned int j = 0; j < dim; ++j) {
                        eigenvecs(j, i) = eigenpairs_new[i].second[j];  // Store eigenvectors column-wise
                }
            }

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvals(k));
            }

            eigenvecs.mmult(U_new,Lambda);
            temp=U_new; temp.mTmult(U_new,eigenvecs);
            R_new.invert(U_new);
            temp=R_new; Fe_new.mmult(R_new,temp);

            Omega=0.0; Omega.add(1.0,R_new); Omega.add(-1.0,R_old);
            temp=Omega; temp.mTmult(Omega,R_new);

            rold=rotnew_conv[i][j];


            Omega_vec(0)=-0.5*(Omega(1,2)-Omega(2,1));Omega_vec(1)=0.5*(Omega(0,2)-Omega(2,0));Omega_vec(2)=-0.5*(Omega(0,1)-Omega(1,0));

            double dot;
            dot=Omega_vec(0)*rold(0)+Omega_vec(1)*rold(1)+Omega_vec(2)*rold(2);
            Vector<double> cross(dim),dot_term(dim);

            dot_term(0)=dot*rold(0);dot_term(1)=dot*rold(1);dot_term(2)=dot*rold(2);

            cross(0)=Omega_vec(1)*rold(2)-Omega_vec(2)*rold(1);
            cross(1)=Omega_vec(2)*rold(0)-Omega_vec(0)*rold(2);
            cross(2)=Omega_vec(0)*rold(1)-Omega_vec(1)*rold(0);

            dr=0.0;	dr.add(1.0, Omega_vec); dr.add(1.0,dot_term); dr.add(1.0,cross); dr.equ(0.5,dr);

            rnew=0.0; rnew.add(1.0,rold); rnew.add(1.0,dr);
            
            ///////Very large Rodrigues vector norm leads to Nan or Inf for updated rnew. Accordingly, we keep the maximum norm to max_rnew_Norm=10000.
            double rnew_Norm, max_rnew_Norm;
            max_rnew_Norm=10000;
            rnew_Norm=sqrt(rnew(0)*rnew(0)+rnew(1)*rnew(1)+rnew(2)*rnew(2));

            if (rnew_Norm>max_rnew_Norm){
                rnew(0)=rnew(0)*max_rnew_Norm/rnew_Norm;
                rnew(1)=rnew(1)*max_rnew_Norm/rnew_Norm;
                rnew(2)=rnew(2)*max_rnew_Norm/rnew_Norm;
            }
            ///////////////////////

            rotnew_conv[i][j]=rnew;


        }
    }


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
