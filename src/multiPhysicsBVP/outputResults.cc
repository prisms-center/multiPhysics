// outputResults() method for MatrixFreePDE class
#include <deal.II/numerics/data_out.h>

#include "../../include/matrixFreePDE.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// helper functions
bool
directory_exists(const std::string &path)
{
  struct stat info;
  return stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR);
}

bool
create_directory_if_not_exists(const std::string &path)
{
  if (!directory_exists(path))
    {
      return mkdir(path.c_str(), 0755) == 0;
    }
  return true;
}

// output results
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::outputResults()
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: output");

  // create DataOut object
  DataOut<dim> data_out;

  // loop over fields

  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      // mark field as scalar/vector
      std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType(
        fields[fieldIndex].numComponents,
        (fields[fieldIndex].type == SCALAR
           ? DataComponentInterpretation::component_is_scalar
           : DataComponentInterpretation::component_is_part_of_vector));
      // add field to data_out
      std::vector<std::string> solutionNames(fields[fieldIndex].numComponents,
                                             fields[fieldIndex].name.c_str());
      data_out.add_data_vector(*dofHandlersSet[fieldIndex],
                               *solutionSet[fieldIndex],
                               solutionNames,
                               dataType);
    }

  // Test section for outputting postprocessed fields
  // Currently there are hacks in place, using the matrixFreeObject, invM, constraints,
  // and DoFHandler as the primary variables
  if (userInputs.postProcessingRequired)
    {
      std::vector<vectorType_pf *> postProcessedSet;
      computePostProcessedFields(postProcessedSet);

      unsigned int invM_size = invM.locally_owned_size();
      for (unsigned int fieldIndex = 0; fieldIndex < postProcessedSet.size();
           fieldIndex++)
        {
          for (unsigned int dof = 0;
               dof < postProcessedSet[fieldIndex]->locally_owned_size();
               ++dof)
            {
              postProcessedSet[fieldIndex]->local_element(dof) =
                invM.local_element(dof % invM_size) *
                postProcessedSet[fieldIndex]->local_element(dof);
            }
          constraintsOtherSet[0]->distribute(*postProcessedSet[fieldIndex]);
          postProcessedSet[fieldIndex]->update_ghost_values();
        }

      // Integrate over selected post-processed fields and output them to the screen and a
      // text file
      std::ofstream output_file;

      if (userInputs.num_integrated_fields > 0)
        {
          if (first_integrated_var_output_complete)
            {
              output_file.open("integratedFields.txt", std::ios::app);
            }
          else
            {
              output_file.open("integratedFields.txt", std::ios::out);
            }
          output_file.precision(10);

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              output_file << currentTime;
            }

          for (unsigned int i = 0; i < userInputs.pp_number_of_variables; i++)
            {
              if (userInputs.pp_calc_integral[i])
                {
                  double integrated_field;
                  computeIntegral(integrated_field, i, postProcessedSet);
                  pcout << "Integrated value of "
                        << userInputs.pp_var_name[userInputs.integrated_field_indices[i]]
                        << ": " << integrated_field << std::endl;
                  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                    {
                      output_file
                        << "\t"
                        << userInputs.pp_var_name[userInputs.integrated_field_indices[i]]
                        << "\t" << integrated_field;
                    }
                  integrated_postprocessed_fields.at(
                    userInputs.integrated_field_indices[i]) = integrated_field;
                }
            }
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              output_file << std::endl;
            }
          output_file.close();
          first_integrated_var_output_complete = true;
        }

      // Add the postprocessed fields to data_out
      for (unsigned int fieldIndex = 0; fieldIndex < userInputs.pp_number_of_variables;
           fieldIndex++)
        {
          // mark field as scalar/vector
          unsigned int components;
          if (userInputs.pp_varInfoList[fieldIndex].is_scalar)
            {
              components = 1;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                dataType(components, DataComponentInterpretation::component_is_scalar);
              std::vector<std::string> solutionNames(
                components,
                userInputs.pp_var_name[fieldIndex].c_str());
              // add field to data_out
              data_out.add_data_vector(*dofHandlersSet[0],
                                       *postProcessedSet[fieldIndex],
                                       solutionNames,
                                       dataType);
            }
          else
            {
              components = dim;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                                       dataType(components,
                         DataComponentInterpretation::component_is_part_of_vector);
              std::vector<std::string> solutionNames(
                components,
                userInputs.pp_var_name[fieldIndex].c_str());
              // add field to data_out
              // data_out.add_data_vector(*vector_dofHandler,
              // *postProcessedSet[fieldIndex], solutionNames, dataType);
              data_out.add_data_vector(*dofHandlersSet[0],
                                       *postProcessedSet[fieldIndex],
                                       solutionNames,
                                       dataType);
            }
        }
    }

  data_out.build_patches(degree);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      if (!create_directory_if_not_exists(userInputs.output_directory_pf))
        {
          std::cerr << "Error creating output directory: " << userInputs.output_directory_pf
                    << std::endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

  // Wait for rank 0 to create the directory
  MPI_Barrier(MPI_COMM_WORLD);

  // Write to results file
  std::ostringstream cycleAsString;
  cycleAsString << std::setw(std::floor(std::log10(userInputs.totalIncrements_pf)) + 1)
                << std::setfill('0') << currentIncrement;

  // Create directory prefix (with trailing slash if needed)
  std::string prefix = userInputs.output_directory_pf;
  if (!prefix.empty() && prefix.back() != '/')
    prefix += "/";

  char baseFileName[200], vtuFileName[200];
  snprintf(baseFileName,
           sizeof(baseFileName),
           "%s%s-%s",
           prefix.c_str(),
           userInputs.output_file_name.c_str(),
           cycleAsString.str().c_str());

  snprintf(vtuFileName,
           sizeof(vtuFileName),
           "%s.%u.%s",
           baseFileName,
           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
           userInputs.output_file_type.c_str());

  // Write to file in either vtu or vtk format
  if (userInputs.output_file_type == "vtu")
    {
      dealii::DataOutBase::VtkFlags flags;
      flags.time                = currentTime;
      flags.cycle               = currentIncrement;
      flags.print_date_and_time = true;
      data_out.set_flags(flags);

      if (userInputs.output_vtu_per_process)
        {
          std::ofstream output(vtuFileName);
          data_out.write_vtu(output);

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              std::vector<std::string> filenames;
              for (unsigned int i = 0;
                   i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
                   ++i)
                {
                  char vtuProcFileName[200];
                  snprintf(vtuProcFileName,
                           sizeof(vtuProcFileName),
                           "%s%s-%s.%u.%s",
                           prefix.c_str(),
                           userInputs.output_file_name.c_str(),
                           cycleAsString.str().c_str(),
                           i,
                           userInputs.output_file_type.c_str());
                  filenames.push_back(vtuProcFileName);
                }

              char pvtuFileName[200];
              snprintf(pvtuFileName,
                       sizeof(pvtuFileName),
                       "%s.p%s",
                       baseFileName,
                       userInputs.output_file_type.c_str());

              std::ofstream master_output(pvtuFileName);
              data_out.write_pvtu_record(master_output, filenames);
              pcout << "PF Output written to: " << pvtuFileName << "\n\n";
            }
        }
      else
        {
          char svtuFileName[200];
          snprintf(svtuFileName,
                   sizeof(svtuFileName),
                   "%s.%s",
                   baseFileName,
                   userInputs.output_file_type.c_str());
          data_out.write_vtu_in_parallel(svtuFileName, MPI_COMM_WORLD);
          pcout << "PF Output written to: " << svtuFileName << "\n\n";
        }
    }
  else if (userInputs.output_file_type == "vtk")
    {
      std::ofstream output(vtuFileName);
      data_out.write_vtk(output);
      pcout << "PF Output written to: " << vtuFileName << "\n\n";
    }
  else
    {
      std::cerr << "PRISMS-PF Error: The parameter 'outputFileType' must be either "
                   "\"vtu\" or \"vtk\""
                << std::endl;
      abort();
    }
  // log time
  computing_timer.leave_subsection("matrixFreePDE: output");
}

#include "../../include/matrixFreePDE_template_instantiations.h"