/*
  Copyright (C) 2020 by the authors of the ASPECT code.
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


#include <aspect/material_model/rheology/viscous_damper.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      ViscousDamper<dim>::ViscousDamper ()
      {}



      template <int dim>
      double
      ViscousDamper<dim>::compute_viscosity (const double yield_stress,
                                             const double strain_rate,
                                             const double pre_yield_viscosity) const
      {
        const double total_stress = (yield_stress + 2.*damper_viscosity*strain_rate)/
                                    (1. + damper_viscosity/pre_yield_viscosity);

        const double plastic_strain_rate = strain_rate - total_stress/(2.*pre_yield_viscosity);

        const double plastic_viscosity = yield_stress/(2.*plastic_strain_rate) + damper_viscosity;

        return plastic_viscosity;
      }



      template <int dim>
      void
      ViscousDamper<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Damper viscosity", "1e21",
                           Patterns::Double(0),
                           "Viscosity of a damper used to stabalize the plasticity formulation and "
                           "provide an internal length scale. Decreasing the viscosity of the damper "
                           "will lead to narrower plastic shear bands. Units: Pa s.");
      }



      template <int dim>
      void
      ViscousDamper<dim>::parse_parameters (ParameterHandler &prm)
      {
        damper_viscosity = prm.get_double("Damper viscosity");
      }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class ViscousDamper<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
