/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#include <aspect/particle/property/heat_flux.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      HeatFlux<dim>::HeatFlux ()
        :
        material_inputs(1,0),
        material_outputs(1,0)
      {}



      template <int dim>
      void
      HeatFlux<dim>::initialize ()
      {
        material_inputs = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());

        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

      }



      template <int dim>
      void
      HeatFlux<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                           std::vector<double> &data) const
      {
        
        data.push_back(0.);

      }



      template <int dim>
      void
      HeatFlux<dim>::update_particle_property(const unsigned int data_position,
                                                   const Vector<double> &solution,
                                                   const std::vector<Tensor<1,dim>> &gradients,
                                                   typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        material_inputs.position[0] = particle->get_location();

#if DEAL_II_VERSION_GTE(9,4,0)
        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                      &(this->get_dof_handler()));
#else
        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(this->get_triangulation()),
                                                                                      &(this->get_dof_handler()));
#endif
        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];


        Tensor<1,dim> temperature_gradient;
	for (unsigned int d = 0; d < dim; ++d)
            temperature_gradient[d] = gradients[this->introspection().component_indices.temperature][d];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        this->get_material_model().evaluate (material_inputs,material_outputs);

        const Tensor<1,dim> advective_flux = material_inputs.velocity[0] * material_inputs.temperature[0] *
                                             material_outputs.densities[0]*material_outputs.specific_heat[0];
 
        const Tensor<1,dim> conductive_flux = -temperature_gradient*
                                              material_outputs.thermal_conductivities[0];

        for (unsigned int i=0; i<dim; ++i)
        particle->get_properties()[data_position + i] = advective_flux[i] + conductive_flux[i];

        }



      template <int dim>
      UpdateTimeFlags
      HeatFlux<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      HeatFlux<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      HeatFlux<dim>::get_property_information() const
      {
	const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("heat_flux",dim));

        return property_information;
      }
  }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(HeatFlux,
                                        "heat flux",
                                        "A plugin in which the particle property tensor is "
                                        "defined as the total elastic stress a particle has "
                                        "accumulated. See the viscoelastic material model "
                                        "documentation for more detailed information.")

    }
  }
}
