#include "../../../include/crystalPlasticity.h"
#include <deal.II/base/symmetric_tensor.h>

template <int dim>
void crystalPlasticity<dim>::reorient2(Vector<double> &rnew, Vector<double> rold, FullMatrix<double> FE_tau, FullMatrix<double> FE_t) {
    //Update the history variables

    SymmetricTensor<2, dim, double> C_old;
    SymmetricTensor<2, dim, double> C_new;

    Vector<double> eigenvals(dim);
    FullMatrix<double> C_old_temp(dim,dim),C_new_temp(dim,dim), Fe_old(dim, dim), Fe_new(dim, dim), eigenvecs(dim,dim),Lambda(dim,dim),U_old(dim,dim),U_new(dim,dim),R_old(dim,dim),R_new(dim,dim),Omega(dim,dim),temp(dim,dim);
    Lambda=IdentityMatrix(dim);
    Omega=0.0;
    Vector<double> rot1(dim),Omega_vec(dim);
    FullMatrix<double> rotmat(dim,dim);

    C_old_temp=0.0;
    C_old.clear();

    Fe_old= FE_t;
    Fe_new= FE_tau;
    FE_t.Tmmult(C_old_temp, FE_t);
    
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
    temp=R_old; FE_t.mmult(R_old,temp);
    FE_tau.Tmmult(C_new_temp, FE_tau);

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
    temp=R_new; FE_tau.mmult(R_new,temp);

    Omega=0.0; Omega.add(1.0,R_new); Omega.add(-1.0,R_old);
    temp=Omega; temp.mTmult(Omega,R_new);

    Omega_vec(0)=-0.5*(Omega(1,2)-Omega(2,1));Omega_vec(1)=0.5*(Omega(0,2)-Omega(2,0));Omega_vec(2)=-0.5*(Omega(0,1)-Omega(1,0));

    double dot;
    dot=Omega_vec(0)*rold(0)+Omega_vec(1)*rold(1)+Omega_vec(2)*rold(2);
    Vector<double> cross(dim),dot_term(dim);

    dot_term(0)=dot*rold(0);dot_term(1)=dot*rold(1);dot_term(2)=dot*rold(2);

    cross(0)=Omega_vec(1)*rold(2)-Omega_vec(2)*rold(1);
    cross(1)=Omega_vec(2)*rold(0)-Omega_vec(0)*rold(2);
    cross(2)=Omega_vec(0)*rold(1)-Omega_vec(1)*rold(0);

    Vector<double> dr(dim);

    dr=0.0;	dr.add(1.0, Omega_vec); dr.add(1.0,dot_term); dr.add(1.0,cross); dr.equ(0.5,dr);

    rnew=0.0; rnew.add(1.0,rold); rnew.add(1.0,dr);
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
