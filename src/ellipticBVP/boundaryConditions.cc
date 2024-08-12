//methods to apply Dirichlet boundary conditons
#include "../../include/multiPhysicsBVP.h"

//Specify Dirichlet boundary conditions
template <int dim>
void multiPhysicsBVP<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  unsigned int i,dof_1, dof_2 ;

  Vector<double> externalMeshParameterBCs(dim),nodalDisplacementBCsToleranceVector(dim);
  double alpha_Torsion;
  for (unsigned int i=0; i<dim; ++i) {
    externalMeshParameterBCs(i)=userInputs_cp.externalMeshParameter*userInputs_cp.span[i];
    nodalDisplacementBCsToleranceVector(i)=userInputs_cp.nodalDisplacementBCsTolerance*userInputs_cp.span[i];
  }

  if(userInputs_cp.enableCyclicLoading){
    if(dof==(userInputs_cp.cyclicLoadingDOF-1)){
      switch (userInputs_cp.cyclicLoadingFace){
        case 1:
        if (node[0] <= externalMeshParameterBCs(0))
        {
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
	}
      	break;

        case 2:
        if (node[0] >= (userInputs_cp.span[0]-externalMeshParameterBCs(0)))
        {
          //pcout<<"Positive12"<<std::endl;
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            //pcout<<"Positive"<<std::endl;
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            //pcout<<"negative"<<std::endl;
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{//pcout<<"Positive12"<<std::endl;
          flag=true;value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
        }
        break;
        case 3:
        if (node[1] <= externalMeshParameterBCs(1))
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
        break;
        case 4:
        if (node[1] >= (userInputs_cp.span[1]-externalMeshParameterBCs(1)))
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
        break;
        case 5:
        if (node[2] <= externalMeshParameterBCs(2))
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
        break;
        case 6:
        if (node[2] >= (userInputs_cp.span[2]-externalMeshParameterBCs(2)))
          if(fmod((currentIncrement_cp*delT),cycleTime)<userInputs_cp.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement_cp*delT),cycleTime)<3*userInputs_cp.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs_cp.cyclicLoadingFace-1][dof];return;}
      }
    }
  }

  if((userInputs_cp.enableSimpleBCs)||(userInputs_cp.enableCyclicLoading)){
    for (i=0;i<2*dim;i++){
      if(faceDOFConstrained[i][dof])
        switch (i+1){
          case 1:
          if (node[0] <= externalMeshParameterBCs(0))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 2:
          if (node[0] >= (userInputs_cp.span[0]-externalMeshParameterBCs(0)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 3:
          if (node[1] <= externalMeshParameterBCs(1))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 4:
          if (node[1] >= (userInputs_cp.span[1]-externalMeshParameterBCs(1)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 5:
          if (node[2] <= externalMeshParameterBCs(2))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 6:
          if (node[2] >= (userInputs_cp.span[2]-externalMeshParameterBCs(2)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
        }
    }
  }

  if(userInputs_cp.enableTorsionBCs){
    currentTime_cp=delT*(currentIncrement_cp+1);

    if (currentIncrement_cp==0){
      timeCounter=1;
    }
    if (currentTime_cp>userInputs_cp.tabularTimeInputTorsion[timeCounter]){
      timeCounter=timeCounter+1;
    }
    alpha_Torsion=userInputs_cp.tabularTorsionBCsInput[timeCounter]*delT;

    if (userInputs_cp.torsionAxis==2){ //z-axis is torsion axis
      dof_1=0; dof_2=1;
    }
    else if (userInputs_cp.torsionAxis==0){ //x-axis is torsion axis
      dof_1=1; dof_2=2;
    }
    else{ //y-axis is torsion axis
      dof_1=2; dof_2=0;
    }

    if ((node[userInputs_cp.torsionAxis] >= (userInputs_cp.span[userInputs_cp.torsionAxis]-externalMeshParameterBCs(userInputs_cp.torsionAxis)))&&(dof==dof_1)){
                flag=true; value=-(node[dof_2]-userInputs_cp.centerTorsion[1])*alpha_Torsion;}
    if ((node[userInputs_cp.torsionAxis] >= (userInputs_cp.span[userInputs_cp.torsionAxis]-externalMeshParameterBCs(userInputs_cp.torsionAxis)))&&(dof==dof_2)){
                flag=true; value=(node[dof_1]-userInputs_cp.centerTorsion[0])*alpha_Torsion;}
  }

  if(userInputs_cp.enableNodalDisplacementBCs){
    for (i=0;i<userInputs_cp.numberOfNodalBCs;i++){
      if ((fabs(node[0]-nodalDisplacement[i][0])<= nodalDisplacementBCsToleranceVector(0))&&(fabs(node[1]-nodalDisplacement[i][1])<= nodalDisplacementBCsToleranceVector(1))&&(fabs(node[2]-nodalDisplacement[i][2])<= nodalDisplacementBCsToleranceVector(2))){
        if (dof==dofNodalDisplacement[i]){
          flag=true; value=deluNodalDisplacement[i];return;
        }
      }
    }
  }

  if(userInputs_cp.enableTabularBCs){

    for (i=0;i<2*dim;i++){
      if(faceDOFConstrained[i][dof])
        switch (i+1){
          case 1:
          if (node[0] <= externalMeshParameterBCs(0))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true;

                currentTime_cp=delT*(currentIncrement_cp+1);

                if (currentIncrement_cp==0){
                  timeCounter=1;
                }
                if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                  timeCounter=timeCounter+1;
                }

                value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

                return;}
          break;
          case 2:
          if (node[0] >= (userInputs_cp.span[0]-externalMeshParameterBCs(0)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime_cp=delT*(currentIncrement_cp+1);

              if (currentIncrement_cp==0){
                timeCounter=1;
              }
              if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 3:
          if (node[1] <= externalMeshParameterBCs(1))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime_cp=delT*(currentIncrement_cp+1);

              if (currentIncrement_cp==0){
                timeCounter=1;
              }
              if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 4:
          if (node[1] >= (userInputs_cp.span[1]-externalMeshParameterBCs(1)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime_cp=delT*(currentIncrement_cp+1);

              if (currentIncrement_cp==0){
                timeCounter=1;
              }
              if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 5:
          if (node[2] <= externalMeshParameterBCs(2))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime_cp=delT*(currentIncrement_cp+1);

              if (currentIncrement_cp==0){
                timeCounter=1;
              }
              if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 6:
          if (node[2] >= (userInputs_cp.span[2]-externalMeshParameterBCs(2)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime_cp=delT*(currentIncrement_cp+1);

              if (currentIncrement_cp==0){
                timeCounter=1;
              }
              if (currentTime_cp>userInputs_cp.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs_cp.tabularTimeInput[timeCounter]-userInputs_cp.tabularTimeInput[timeCounter-1])*delT ;

              return;}
        }
    }
  }

  if(userInputs_cp.useVelocityGrad){
    value = 0;
    for(i=0;i<dim;i++)
      value+=deltaF[dof][i]*node[i];
    flag=true;
    return;
  }

  if(userInputs_cp.enableDICpipeline){

  //back boundary
  if (node[0] <= externalMeshParameterBCs(0)){
    double value_x, value_y;
    bcFunction1(node[1],value_x,value_y, currentIncrement_cp);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}
  }


  //front boundary
  if (node[0] >= (userInputs_cp.span[0]-externalMeshParameterBCs(0))){
    double value_x, value_y;
    bcFunction2(node[1],value_x,value_y, currentIncrement_cp);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}
  }


  //left boundary
  if (node[1] <= externalMeshParameterBCs(1)){
    double value_x, value_y;
    bcFunction3(node[0],value_x,value_y, currentIncrement_cp);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}

  }

  //left boundary
  if (node[1] >= (userInputs_cp.span[1]-externalMeshParameterBCs(1))){
    double value_x, value_y;
    bcFunction4(node[0],value_x,value_y, currentIncrement_cp);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}

  }


  //bottom boundary: u_z=0
  if (node[2] <= externalMeshParameterBCs(2)){
    if(node[0]<= externalMeshParameterBCs(0) && node[1] <= externalMeshParameterBCs(1)){
      if (dof==2) {flag=true; value=0.0;}
    }
  }
  }

}


template <int dim>
void ellipticBVP<dim>::bcFunction1(double yval, double &value_x, double &value_y, double currentIncr){

  unsigned int n_div=userInputs_cp.Y_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(yval-bc_new1[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime_cp=delT*(currentIncrement_cp+1);

  if (currentIncrement_cp==0){
    timeCounter=1;
  }
  if (currentTime_cp>userInputs_cp.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new1[coord_pos][1+(timeCounter-1)*2]+bc_new1[coord_pos][1+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new1[coord_pos][2+(timeCounter-1)*2]+bc_new1[coord_pos][2+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;
  //return value1 ; // displacement along X-Direction

}

template <int dim>
void multiPhysicsBVP<dim>::bcFunction2(double yval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs_cp.Y_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(yval-bc_new2[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime_cp=delT*(currentIncrement_cp+1);

  if (currentIncrement_cp==0){
    timeCounter=1;
  }
  if (currentTime_cp>userInputs_cp.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new2[coord_pos][1+(timeCounter-1)*2]+bc_new2[coord_pos][1+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new2[coord_pos][2+(timeCounter-1)*2]+bc_new2[coord_pos][2+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;



  //return value1 ; // displacement along X-Direction

}

template <int dim>
void multiPhysicsBVP<dim>::bcFunction3(double xval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs_cp.X_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(xval-bc_new3[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime_cp=delT*(currentIncrement_cp+1);

  if (currentIncrement_cp==0){
    timeCounter=1;
  }
  if (currentTime_cp>userInputs_cp.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new3[coord_pos][1+(timeCounter-1)*2]+bc_new3[coord_pos][1+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new3[coord_pos][2+(timeCounter-1)*2]+bc_new3[coord_pos][2+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;


}

template <int dim>
void multiPhysicsBVP<dim>::bcFunction4(double xval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs_cp.X_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(xval-bc_new4[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);
  currentTime_cp=delT*(currentIncrement_cp+1);

  if (currentIncrement_cp==0){
    timeCounter=1;
  }
  if (currentTime_cp>userInputs_cp.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new4[coord_pos][1+(timeCounter-1)*2]+bc_new4[coord_pos][1+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new4[coord_pos][2+(timeCounter-1)*2]+bc_new4[coord_pos][2+timeCounter*2])/(userInputs_cp.timeInputDIC[timeCounter]-userInputs_cp.timeInputDIC[timeCounter-1])*delT ;


}



//methods to apply dirichlet BC's
template <int dim>
void multiPhysicsBVP<dim>::applyDirichletBCs_cp(){
    //pcout<<"debug setDirichlet 0\n";
  if(!userInputs_cp.enablePeriodicBCs){
     // pcout<<"debug setDirichlet b1\n";
    constraints.clear();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dofHandler, constraints);
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
     // pcout<<"debug setDirichlet b2\n";
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs_cp.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs_cp.quadOrder), update_values);
     // pcout<<"debug setDirichlet b3\n";
    //parallel loop over all elements
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
      if (cell->is_locally_owned()){
        cell->get_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        for (unsigned int faceID=0; faceID<GeometryInfo<dim>::faces_per_cell; faceID++){
          if (cell->face(faceID)->at_boundary()){
            fe_face_values.reinit (cell, faceID);
            for (unsigned int i=0; i<dofs_per_cell; ++i) {
              if (fe_face_values.shape_value(i, 0)!=0){
                const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                unsigned int globalDOF=local_dof_indices[i];
                bool flag=false;
                double value=0;
                node=supportPoints[globalDOF];
                setBoundaryValues(node, dof, flag, value);
                if (flag){
                  constraints.add_line (globalDOF);
                  if (currentIteration==0){
                    value*=loadFactorSetByModel;
                  }
                  else{
                    value=0.0;
                  }
                  constraints.set_inhomogeneity(globalDOF, value);
                }
              }
            }
          }
        }
      }
    }
  }
  else{
    setPeriodicityConstraints_cp();
  }
  if (userInputs_cp.enableIndentationBCs){
      //pcout<<"debug setDirichlet c1\n";
      setIndentationConstraints();
      //pcout<<"debug setDirichlet c2\n";
      constraints.merge(indentation_constraints,dealii::AffineConstraints<>::right_object_wins);
      //pcout<<"debug setDirichlet c3\n";
      //indentation conditions will overwrite dirichlet conditions from other sources (as we prefer)
  }
  constraints.close ();
  //pcout<<"debug setDirichlet 4\n";

}
#include "../../include/multiPhysicsBVP_template_instatiations.h"
