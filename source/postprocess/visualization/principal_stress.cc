/*
  Copyright (C) 2020 - 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/postprocess/visualization/principal_stress.h>

#include <deal.II/base/symmetric_tensor.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      PrincipalStress<dim>::
      PrincipalStress ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("Pa")
      {}



      template <int dim>
      std::vector<std::string>
      PrincipalStress<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        // dim principal stress values
        for (unsigned int i=0; i<dim; ++i)
          solution_names.push_back("principal_stress_" + std::to_string(i+1));

        // dim principal stress directions
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            solution_names.push_back("principal_stress_direction_" + std::to_string(i+1));

        return solution_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      PrincipalStress<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> solution_components;

        // dim principal stress values
        for (unsigned int i=0; i<dim; ++i)
          solution_components.push_back(DataComponentInterpretation::component_is_scalar);

        // dim principal stress directions
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            solution_components.push_back(DataComponentInterpretation::component_is_part_of_vector);

        return solution_components;
      }



      template <int dim>
      UpdateFlags
      PrincipalStress<dim>::
      get_needed_update_flags () const
      {
        return update_values | update_gradients | update_quadrature_points;
      }



      template <int dim>
      void
      PrincipalStress<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == dim*dim + dim, ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // We do not need to compute anything but the viscosity and ElasticAdditionalOutputs
        in.requested_properties = MaterialModel::MaterialProperties::viscosity | MaterialModel::MaterialProperties::additional_outputs;

        this->get_material_model().create_additional_named_outputs(out);

        // Compute the viscosity and additional outputs
        this->get_material_model().evaluate(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            SymmetricTensor<2,dim> stress;

            const double eta = out.viscosities[q];

            const SymmetricTensor<2, dim> deviatoric_strain_rate = (this->get_material_model().is_compressible()
                                                                    ? in.strain_rate[q] - 1. / 3. * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>()
                                                                    : in.strain_rate[q]);

            // Add elastic stresses if existent
            if (this->get_parameters().enable_elasticity == false)
              {
                // Compressive stress is positive in geoscience applications
                stress = -2. * eta * deviatoric_strain_rate;
              }
            else
              {
                // Visco-elastic stresses are stored on the fields
                SymmetricTensor<2, dim> stress_0, stress_old;
                stress_0[0][0] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
                stress_0[1][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
                stress_0[0][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

                stress_old[0][0] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx_old")];
                stress_old[1][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy_old")];
                stress_old[0][1] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy_old")];

                if (dim == 3)
                  {
                    stress_0[2][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                    stress_0[0][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                    stress_0[1][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];

                    stress_old[2][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz_old")];
                    stress_old[0][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz_old")];
                    stress_old[1][2] = in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz_old")];
                  }

                const MaterialModel::ElasticAdditionalOutputs<dim> *elastic_additional_out = out.template get_additional_output<MaterialModel::ElasticAdditionalOutputs<dim>>();

                Assert(elastic_additional_out != nullptr, ExcMessage("Elastic Additional Outputs are needed for the 'principal stress' postprocessor, but they have not been created."));
                const double shear_modulus = elastic_additional_out->elastic_shear_moduli[q];

                // Retrieve the timestep ratio dtc/dte and the elastic viscosity,
                // only two material models support elasticity.
                double timestep_ratio = 0.;
                double elastic_viscosity = 0.;
                if (Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim>>(this->get_material_model()))
                  {
                    const MaterialModel::ViscoPlastic<dim> &vp = Plugins::get_plugin_as_type<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model());
                    elastic_viscosity = vp.get_elastic_viscosity(shear_modulus);
                    timestep_ratio = vp.get_timestep_ratio();
                  }
                else if (Plugins::plugin_type_matches<MaterialModel::Viscoelastic<dim>>(this->get_material_model()))
                  {
                    const MaterialModel::Viscoelastic<dim> &ve = Plugins::get_plugin_as_type<const MaterialModel::Viscoelastic<dim>>(this->get_material_model());
                    elastic_viscosity = ve.get_elastic_viscosity(shear_modulus);
                    timestep_ratio = ve.get_timestep_ratio();
                  }
                else
                  AssertThrow(false, ExcMessage("The principle stress visualization plugin cannot be used with elasticity for material models other than ViscoPlastic and Viscoelastic."));

                // Apply the stress update to get the total stress of timestep t.
                // Both eta and the elastic viscosity have been scaled with the timestep ratio.
                stress = in.pressure[q] * unit_symmetric_tensor<dim>() -
                         (2. * eta * deviatoric_strain_rate + eta / elastic_viscosity * stress_0 + (1. - timestep_ratio) * (1. - eta / elastic_viscosity) * stress_old);
              }

            if (use_deviatoric_stress == false)
              stress += in.pressure[q] * unit_symmetric_tensor<dim>();

            const std::array<std::pair<double, Tensor<1,dim>>, dim> principal_stresses_and_directions = eigenvectors(stress);

            // dim principal stress values
            for (unsigned int i=0; i<dim; ++i)
              computed_quantities[q][i] = principal_stresses_and_directions[i].first;

            // dim principal stress directions
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                computed_quantities[q][dim + i*dim + j] = principal_stresses_and_directions[i].second[j];
          }
      }



      template <int dim>
      void
      PrincipalStress<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Principal stress");
            {
              prm.declare_entry("Use deviatoric stress", "false",
                                Patterns::Bool(),
                                "Whether to use the deviatoric stress tensor "
                                "instead of the full stress tensor to compute "
                                "principal stress directions and values.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      PrincipalStress<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Principal stress");
            {
              use_deviatoric_stress = prm.get_bool("Use deviatoric stress");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(PrincipalStress,
                                                  "principal stress",
                                                  "A visualization output object that outputs the "
                                                  "principal stresses and directions, i.e., the "
                                                  "eigenvalues and eigenvectors of the stress tensor. "
                                                  "Wikipedia defines principal stresses as follows: "
                                                  "At every point in a stressed body there are at least "
                                                  "three planes, called principal planes, with normal "
                                                  "vectors, called principal directions, where the "
                                                  "corresponding stress vector is perpendicular to the "
                                                  "plane, and where there are no normal "
                                                  "shear stresses. The three stresses normal to these "
                                                  "principal planes are called principal stresses. "
                                                  "This postprocessor can either operate on the full "
                                                  "stress tensor or only on the deviatoric stress tensor, "
                                                  "depending on what run-time parameters are set."
                                                  "\n\n"
                                                  "Physical units: \\si{\\pascal}.")
    }
  }
}
