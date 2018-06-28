/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/melt_visco_plastic.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    MeltViscoPlastic<dim>::initialize ()
    {
      base_model->initialize();
    }


    template <int dim>
    double
    MeltViscoPlastic<dim>::
    reference_viscosity () const
    {
      return base_model->reference_viscosity();
    }

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    MeltViscoPlastic<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    std::vector<double>
    MeltViscoPlastic<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions( compositional_fields.size()+1);
      double sum_composition = 0.0;

      std::vector<double> x_comp = compositional_fields;
      for (unsigned int i=0; i < x_comp.size(); ++i)
        {
          // The first element represents the background mantle, which is
          // not included in x_comp
          if (field_used_in_viscosity_averaging[i+1] == true)
            {
              // clip the compositional fields so they are between zero and one
              // and sum them for normalization purposes
              x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);
              sum_composition += x_comp[i];
            }
          else
            x_comp[i] = 0.0;
        }

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  //background material
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background material
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;
    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.temperature.size(); ++q)
        melt_fractions[q] = melt_fraction(in.temperature[q],
                                          std::max(0.0, in.pressure[q]));
      return;
    }


    template <int dim>
    void
    MeltViscoPlastic<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // 1) Get the background viscosities

      // get all material properties from the visco-plastic model
      base_model->evaluate(in, out);

      // Modify viscosity if not using viscosity values from the base model
      if (use_linear_viscosities == true)
        {
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              // Compute volume fractions
              const std::vector<double> volume_fractions = compute_volume_fractions(in.composition[i]);
              double linear_viscosity = 0.;
              for (unsigned int c=0; c< volume_fractions.size(); ++c)
                {
                  linear_viscosity += volume_fractions[c] * linear_viscosities[c];
                }
              out.viscosities[i] = linear_viscosity;
            }
        }

      // Store the intrinsic viscosities for computing the compaction viscosities later on
      // (Keller et al. eq. 51).
      const std::vector<double> xis = out.viscosities;

      // 2) Retrieve porosity, depletion, fluid pressure and volumetric strain rate
      // for timesteps bigger than 0.

      std::vector<double> maximum_melt_fractions(in.position.size());
      std::vector<double> old_porosity(in.position.size());
      std::vector<double> fluid_pressures(in.position.size());
      std::vector<double> volumetric_strain_rates(in.position.size());
      std::vector<double> volumetric_yield_strength(in.position.size());

      // The porosity (if no melt transport is used, it is set to zero)
      double porosity = 0.0;

      // we want to get the porosity and the peridotite field (the depletion) from the old
      // solution here, so we can compute differences between the equilibrium in the current
      // time step and the melt generated up to the previous time step
      if (this->include_melt_transport() )
        {
          // get peridotite and porosity field indices
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

          if ( in.current_cell.state() == IteratorState::valid
               && this->get_timestep_number() > 0)
            {
              // Prepare the field function
              Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
              fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());


              fe_value.set_active_cell(in.current_cell);
              fe_value.value_list(in.position,
                                  maximum_melt_fractions,
                                  this->introspection().component_indices.compositional_fields[peridotite_idx]);

              fe_value.value_list(in.position,
                                  old_porosity,
                                  this->introspection().component_indices.compositional_fields[porosity_idx]);

              // get fluid pressure from the current solution
              Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
              fe_value_current(this->get_dof_handler(), this->get_solution(), this->get_mapping());
              fe_value_current.set_active_cell(in.current_cell);

              fe_value_current.value_list(in.position,
                                          fluid_pressures,
                                          this->introspection().variable("fluid pressure").first_component_index);

              // get volumetric strain rate
              // see Keller et al. eq. 11.
              std::vector<Tensor<1,dim> > velocity_gradients(in.position.size());
              for (unsigned int d=0; d<dim; ++d)
                {
                  fe_value_current.gradient_list(in.position,
                                                 velocity_gradients,
                                                 this->introspection().component_indices.velocities[d]);
                  for (unsigned int i=0; i<in.position.size(); ++i)
                    volumetric_strain_rates[i] += velocity_gradients[i][d];
                }
            }


          // 3) Get porosity, melt density and update melt reaction terms
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              // we can not use the densities from the visco-plastic model
              // as they assume compositional fields between 0 and 1, and the depletion field can become negative
              const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                       ?
                                       depletion_density_change * in.composition[i][peridotite_idx]
                                       :
                                       0.0;
              out.densities[i] += delta_rho;

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              double melting = 0.0;
              if (this->get_timestep_number() > 0)
                {
                  // batch melting
                  melting = melt_fraction(in.temperature[i], in.pressure[i])
                            - std::max(maximum_melt_fractions[i], 0.0);
                }
              // freezing of melt below the solidus
              {
                const double freezing_potential = melt_fraction(in.temperature[i], in.pressure[i]) - old_porosity[i];
                const double freezing = freezing_rate * this->get_timestep() / year_in_seconds * 0.5 * (freezing_potential - std::abs(freezing_potential));
                melting += freezing;
              }

              // do not allow negative porosity
              if (old_porosity[i] + melting < 0)
                melting = -old_porosity[i];

              // because depletion is a volume-based, and not a mass-based property that is advected,
              // additional scaling factors on the right hand side apply
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  if (c == peridotite_idx && this->get_timestep_number() > 0 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting * (1 - maximum_melt_fractions[i])
                                               / (1 - maximum_melt_fractions[i]);
                  else if (c == porosity_idx && this->get_timestep_number() > 0 && (in.strain_rate.size()))
                    out.reaction_terms[i][c] = melting
                                               * out.densities[i] / this->get_timestep();
                  else
                    out.reaction_terms[i][c] = 0.0;
                }

              porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // 4) Reduce shear viscosity due to melt presence
              if ( in.strain_rate.size() )
                {
                  // melt dependence of shear viscosity
                  out.viscosities[i] *= exp(- alpha_phi * porosity);
                }
            }
        }

      if ( in.strain_rate.size() )
        {

          // 5) Compute plastic weakening of visco(elastic) viscosity
          for (unsigned int i=0; i<in.position.size(); ++i)
            {

              // Compute volume fractions
              const std::vector<double> volume_fractions = compute_volume_fractions(in.composition[i]);

              // calculate deviatoric strain rate (Keller et al. eq. 13)
              const double edot_ii = ( (this->get_timestep_number() == 0 && in.strain_rate[i].norm() <= std::numeric_limits<double>::min())
                                       ?
                                       ref_strain_rate
                                       :
                                       std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                                                min_strain_rate) );

              // compute visco(elastic) stress
              const double viscous_stress = 2. * out.viscosities[i] * edot_ii * (1.0 - porosity);

              // overwrite the reaction terms, which is needed to track the finite strain invariant.
              // TODO clarify what field tracks strain
              // TODO I removed the condition time_step > 0, because there can be initial strain present, correct?
              // TODO so far, strain weakening is not implemented anywhere else in this material model (e.g. compute_volume_fractions,
              // or to actually change the cohesions and friction angles)
              if (use_strain_weakening && in.strain_rate.size())
                {
                  // New strain invariant is old strain invariant plus the strain invariant of the current time step
                  out.reaction_terms[i][0] = edot_ii*this->get_timestep();
                }
              // In case porosity lies above the melt transport threshold
              // P_effective = P_bulk - P_f = (1-porosity) * P_s + porosity * P_f - P_f = (1-porosity) * (P_s - P_f)
              // otherwise,
              // P_effective = P_bulk, which equals P_solid (which is given by in.pressure[i])
              const double effective_pressure = ((this->include_melt_transport() && this->get_melt_handler().is_melt_cell(in.current_cell))
                                                 ?
                                                 (1. - porosity) * (in.pressure[i] - fluid_pressures[i])
                                                 :
                                                 in.pressure[i]);

              double yield_strength = 0.0;
              double tensile_strength = 0.0;

              for (unsigned int c=0; c< volume_fractions.size(); ++c)
                {
                  const double tensile_strength_c = cohesions[c]/strength_reductions[c];

                  // Convert friction angle from degrees to radians
                  double phi = angles_internal_friction[c] * numbers::PI/180.0;
                  const double transition_pressure = (cohesions[c] * std::cos(phi) - tensile_strength_c) / (1.0 -  sin(phi));

                  double yield_strength_c = 0.0;
                  // In case we're not using the Keller et al. formulation,
                  // or the effective pressure is bigger than the transition pressure, use
                  // the normal yield strength formulation
                  if (effective_pressure > transition_pressure || !this->include_melt_transport())
                    yield_strength_c = ( (dim==3)
                                         ?
                                         ( 6.0 * cohesions[c] * std::cos(phi) + 2.0 * effective_pressure * std::sin(phi) )
                                         / ( std::sqrt(3.0) * (3.0 + std::sin(phi) ) )
                                         :
                                         cohesions[c] * std::cos(phi) + effective_pressure * std::sin(phi) );
                  else
                    // Note typo in Keller et al. paper eq. (37) (minus sign)
                    yield_strength_c = tensile_strength_c + effective_pressure;

                  // TODO add different averagings?
                  yield_strength += volume_fractions[c]*yield_strength_c;
                  tensile_strength += volume_fractions[c]*tensile_strength_c;
                }

              // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface
              // TODO add in shear stress evolution parameter! (or in viscous_stress) (see Keller et al. eq (29) and eq (40))
              if (viscous_stress >= yield_strength)
                out.viscosities[i] = yield_strength / (2.0 * edot_ii);

              // Limit the viscosity with specified minimum and maximum bounds
              out.viscosities[i] = std::min(std::max(out.viscosities[i], min_visc), max_visc);

              // Compute the volumetric yield strength (Keller et al. eq (38))
              volumetric_yield_strength[i] = viscous_stress - tensile_strength;
            }
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      if (melt_out != NULL)
        {
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = (this->get_melt_handler().is_melt_cell(in.current_cell)
                                             ?
                                             std::max(reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2),0.0)
                                             :
                                             0.0);

              melt_out->fluid_densities[i] = out.densities[i] + melt_density_change;
              melt_out->fluid_density_gradients[i] = 0.0;

              const double compaction_pressure = (1.0 - porosity) * (in.pressure[i] - fluid_pressures[i]);

              const double phi_0 = 0.05;
              porosity = std::max(std::min(porosity,0.995),1.e-7);
              // compaction viscosities (Keller et al. eq (51)
              melt_out->compaction_viscosities[i] = xis[i] * phi_0 / porosity;

              // visco(elastic) compaction viscosity
              // Keller et al. eq (36)
              // TODO include elastic part
              melt_out->compaction_viscosities[i] *= (1. - porosity);

              // TODO compaction stress evolution parameter

              // effective compaction viscosity (Keller et al. eq (43) )
              // NB: I've added a minus sign as according to eq 43
              if (in.strain_rate.size() && compaction_pressure > volumetric_yield_strength[i])
                melt_out->compaction_viscosities[i] =  -volumetric_yield_strength[i] / std::max(volumetric_strain_rates[i], min_strain_rate);

              // Limit the viscosity with specified minimum and maximum bounds
              melt_out->compaction_viscosities[i] = std::min(std::max(melt_out->compaction_viscosities[i], min_visc), max_visc);
            }
        }


    }


    template <int dim>
    void
    MeltViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {
          prm.declare_entry("Base model","diffusion dislocation",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a visco-plastic  "
                            "rheology and melt migration outputs. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Melt density change", "-500",
                             Patterns::Double (),
                             "Difference between solid density $\\rho_{s}$ and melt/fluid$\\rho_{f}$. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Freezing rate", "0.0",
                             Patterns::Double (0),
                             "Freezing rate of melt when in subsolidus regions."
                             "Units: $1/yr$.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of"
                             "depleted material. Depletion is indicated by the compositional"
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model."
                             "Units: $kg/m^3$.");

          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");

          prm.declare_entry ("Linear viscosities", "1.e22",
                             Patterns::List(Patterns::Double(0)),
                             "List of linear viscosities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The values can be used instead of the viscosities derived from the "
                             "base material model. Units: Pa s.");
          prm.declare_entry ("Use linear viscosities", "false",
                             Patterns::Bool (),
                             "Use user-specified linear visocsity values in replace of viscosity "
                             "values derived from the DiffusionDislocation material model. "
                             "Units: None");

          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");

          // Strain weakening parameters
          prm.declare_entry ("Use strain weakening", "false",
                             Patterns::Bool (),
                             "Apply strain weakening to viscosity, cohesion and internal angle "
                             "of friction based on accumulated finite strain.  Units: None");
          prm.declare_entry ("Start strain weakening intervals", "0.",
                             Patterns::List(Patterns::Double(0)),
                             "List of strain weakening interval initial strains "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None");
          prm.declare_entry ("End strain weakening intervals", "1.",
                             Patterns::List(Patterns::Double(0)),
                             "List of strain weakening interval final strains "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None");
          prm.declare_entry ("Cohesion strain weakening factors", "1.",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesion strain weakening factors "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None");
          prm.declare_entry ("Friction strain weakening factors", "1.",
                             Patterns::List(Patterns::Double(0)),
                             "List of friction strain weakening factors "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None");

          // Plasticity parameters
          prm.declare_entry ("Angles of internal friction", "0",
                             Patterns::List(Patterns::Double(0)),
                             "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesions", "1e20",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesions, $C$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: $Pa$.");
          prm.declare_entry ("Host rock strength reductions", "4",
                             Patterns::List(Patterns::Double(0)),
                             "List of reduction factors of strength of the host rock under tensile stress, $R$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "Units: none.");
          prm.declare_entry ("Use compositional field for viscosity averaging", "1",
                             Patterns::List(Patterns::Integer(0,1)),
                             "List of integers, detailing for each compositional field if it should be included in the "
                             "averaging scheme when the yield strength is computed (if 1) or not (if 0).");
          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");

          // Elasticity parameters
          prm.declare_entry ("Elastic shear moduli", "75.0e9",
                             Patterns::List(Patterns::Double(0)),
                             "List of elastic shear moduli, $G$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The default value of 75 GPa is representative of mantle rocks. Units: none.");
          prm.declare_entry ("Elastic time step", "1.e3",
                             Patterns::Double (0),
                             "The elastic time step $elastic_time_step$. Units: $yr$.");
          prm.declare_entry ("Model is viscoelastic", "false",
                             Patterns::Bool (),
                             "Viscosity is modified according elastic parameters. Units: None ");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {
          AssertThrow( prm.get("Base model") != "melt visco plastic",
                       ExcMessage("You may not use ``melt visco plastic'' itself as the base model for "
                                  "the melt visco plastic model.") );

          // check if the applicable compositional fields exist
          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt peridotite eclogite only works if there is a "
                                 "compositional field called peridotite."));

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt peridotite eclogite with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          melt_density_change        = prm.get_double ("Melt density change");
          xi_0                       = prm.get_double ("Reference bulk viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
          freezing_rate              = prm.get_double ("Freezing rate");
          depletion_density_change   = prm.get_double ("Depletion density change");

          ref_strain_rate = prm.get_double("Reference strain rate");

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          M_cpx           = prm.get_double ("Mass fraction cpx");

          // Plasticity parameters
          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");
          strength_reductions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Host rock strength reductions"))),
                                                                        n_fields,
                                                                        "Host rock strength reductions");
          field_used_in_viscosity_averaging = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_int(Utilities::split_string_list(prm.get("Use compositional field for viscosity averaging"))),
                                                                                      n_fields,
                                                                                      "Use compositional field for viscosity averaging");
          min_strain_rate = prm.get_double("Minimum strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");

          use_linear_viscosities = prm.get_bool ("Use linear viscosities");
          linear_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear viscosities"))),
                                                                       n_fields,
                                                                       "Linear viscosities");

          // Strain weakening parameters
          use_strain_weakening             = prm.get_bool ("Use strain weakening");
          start_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Start strain weakening intervals"))),
                                                                                     n_fields,
                                                                                     "Start strain weakening intervals");
          end_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("End strain weakening intervals"))),
                                                                                   n_fields,
                                                                                   "End strain weakening intervals");
          cohesion_strain_weakening_factors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesion strain weakening factors"))),
                                                                                      n_fields,
                                                                                      "Cohesion strain weakening factors");
          friction_strain_weakening_factors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Friction strain weakening factors"))),
                                                                                      n_fields,
                                                                                      "Friction strain weakening factors");

          // Elastic parameters
          elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                                                         n_fields,
                                                                         "Elastic shear moduli");
          elastic_time_step = prm.get_double ("Elastic time step");
          if (this->convert_output_to_years ())
            elastic_time_step *= year_in_seconds;

          model_is_viscoelastic = prm.get_bool ("Model is viscoelastic");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for depth dependent, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltViscoPlastic,
                                   "melt visco plastic",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to the melting "
                                   "model for dry peridotite of \\cite{KSL2003}. All other material "
                                   "properties are taken from the visco-plastic model.")
  }
}
