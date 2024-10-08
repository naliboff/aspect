#!/bin/bash

# Run a sequence of benchmark models with varying resolution
# and for different advection methods. Individual advection methods
# can be activated/deactivated with the variables below.

run_analytical_density=true
run_dg_compositional_field=true
run_continuous_compositional_field=true
run_particle_models=true

processes=8
ASPECT_EXEC="../plugin/aspect"

for refinement in 2 3 4; do
  echo "subsection Mesh refinement" > resolution.prm
  echo "  set Initial global refinement = $refinement" >> resolution.prm
  echo "end" >> resolution.prm

  if [ $run_analytical_density = true ]; then
    echo "subsection Annulus benchmark" > analytic.prm
    echo "  set Use analytical density = true" >> analytic.prm
    echo "end" >> analytic.prm

    # Reference model at this resolution
    echo "set Output directory = refinement_${refinement}_analytical_density" > output.prm
    echo "Starting refinement_${refinement}_analytical_density" > log
    cat transient_annulus.prm resolution.prm analytic.prm output.prm | mpirun -np $processes $ASPECT_EXEC --
  fi

  if [ $run_dg_compositional_field = true ]; then
    echo "subsection Compositional fields" > field.prm
    echo "  set Number of fields = 1" >> field.prm
    echo "  set Compositional field methods = field" >> field.prm
    echo "end" >> field.prm

    # Reference model at this resolution
    echo "set Output directory = refinement_${refinement}_compositional_field" > output.prm
    echo "Starting refinement_${refinement}_compositional_field" > log
    cat transient_annulus.prm resolution.prm field.prm output.prm | mpirun -np $processes $ASPECT_EXEC --
  fi

  if [ $run_continuous_compositional_field = true ]; then
    echo "subsection Compositional fields" > continuous_field.prm
    echo "  set Number of fields = 1" >> continuous_field.prm
    echo "  set Compositional field methods = field" >> continuous_field.prm
    echo "end" >> continuous_field.prm

    echo "subsection Discretization" >> continuous_field.prm
    echo "set Use discontinuous composition discretization = false" >> continuous_field.prm
    echo "end" >> continuous_field.prm

    # Reference model at this resolution
    echo "set Output directory = refinement_${refinement}_continuous_compositional_field" > output.prm
    echo "Starting refinement_${refinement}_continuous_compositional_field" > log
    cat transient_annulus.prm resolution.prm continuous_field.prm output.prm | mpirun -np $processes $ASPECT_EXEC --
  fi

  if [ $run_particle_models = true ]; then
    particles_per_cell=64
    number_of_particles=`echo "12 * 4^${refinement} * $particles_per_cell" | bc -l`

    for higher_order_time in 'true' 'false'; do

        echo "subsection Compositional fields" > particles.prm
        echo "  set Number of fields = 1" >> particles.prm
        echo "  set Compositional field methods = particles" >> particles.prm
        echo "  set Mapped particle properties = density_field:particle_density" >> particles.prm
        echo "end" >> particles.prm

        echo "subsection Postprocess" >> particles.prm
        echo "  set List of postprocessors = visualization, velocity statistics, AnnulusPostprocessor, rotation statistics, particles, particle count statistics" >> particles.prm
        echo "end" >> particles.prm

        echo "subsection Particles" >> particles.prm
        echo "  subsection Generator" >> particles.prm
        echo "    subsection Probability density function" >> particles.prm
        echo "      set Number of particles = $number_of_particles" >> particles.prm
        echo "    end" >> particles.prm
        echo "  end" >> particles.prm
        echo "  set Integration scheme = rk2" >> particles.prm
        echo "  subsection Integrator" >> particles.prm
        echo "    subsection RK2" >> particles.prm
        echo "      set Higher order accurate in time = $higher_order_time" >> particles.prm
        echo "    end" >> particles.prm
        echo "  end" >> particles.prm
        echo "end" >> particles.prm

        echo "set Output directory = refinement_${refinement}_higher_order_${higher_order_time}" > output.prm
        echo "Starting refinement_${refinement}_higher_order_${higher_order_time}" > log
        cat transient_annulus.prm resolution.prm particles.prm output.prm | mpirun -np $processes $ASPECT_EXEC --
    done
  fi
done
