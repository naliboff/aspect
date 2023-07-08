/*
  Copyright (C) 2023 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_iterative_viscosity_dampening_h
#define _aspect_material_model_rheology_iterative_viscosity_dampening_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {

      template <int dim>
      class IterativeViscosityDampening : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * A function that fills the reaction terms for the finite strain invariant(s) in
           * MaterialModelOutputs object that is handed over.
           */
          void fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                      const int i,
                                      const double old_viscosity,
                                      MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * A function that calculate the dampened viscosity between successive nonlinear iterations
           */
          double calculate_viscosity (const double old_viscosity,
                                      const double new_viscosity) const;

          /**
           * A function that returns a ComponentMask, which indicates that components
           * associated with viscosity should be excluded during the volume fraction computation.
           */
          ComponentMask get_viscosity_composition_mask() const;

          /**
           * A dampening factor for the viscosity that controls the rate of change
           * between the viscosity calculated on the previous and current nonlinear
           * itearation.
           */
          double iterative_viscosity_dampening_factor;

        private:

      };
    }
  }
}
#endif
