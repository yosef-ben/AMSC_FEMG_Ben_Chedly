#ifndef HEAT_HPP
#define HEAT_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class representing the non-linear diffusion problem.
class Heat
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 1;


  class DiffusionCoefficient : public Function<dim>
  {
  public:
  // Constructor.
    DiffusionCoefficient()
    {}

    //Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

  class TransportCoefficient : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
        values[i] = 0.0;

      values[dim - 1] = -g;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      if (component == dim - 1)
        return -g;
      else
        return 0.0;
    }

  protected:
    const double g = 0.0;
  };

  // Forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
  // Constructor.
    ForcingTerm()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  ///// Boundary conditions://///

  class FunctionG : public Function<dim>
  {
  public:
    // Constructor.
    FunctionG()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return p[0] + p[1];
    }
  };


  // Exact solution.
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return std::sin(5 * M_PI * get_time()) * std::sin(2 * M_PI * p[0]) *
             std::sin(3 * M_PI * p[1]) * std::sin(4 * M_PI * p[2]);
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      Tensor<1, dim> result;

      result[0] = 2.0 * std::exp(2.0 * (p[0] - get_time()));

      result[1] = 0;

      return result;
    }
  };

  //Initial condition:
  class FunctionU0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]) * p[2] * (1.0 - p[2]);
    }
  };
  

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  Heat(const unsigned int  &N_,
       const unsigned int &r_,
       const double       &T_,
       const double       &deltat_,
       const double       &theta_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , N(N_)
    , r(r_)
    , deltat(deltat_)
    , theta(theta_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

  // Compute the error.
  double
  compute_error(const VectorTools::NormType &norm_type);

protected:
  // Assemble the mass and stiffness matrices.
  void
  assemble_matrices();

  // Assemble the right-hand side of the problem.
  void
  assemble_rhs(const double &time);

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output.
  void
  output(const unsigned int &time_step) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////



  DiffusionCoefficient diffusion_coefficient;

  TransportCoefficient transport_coefficient;

  // Forcing term.
  ForcingTerm forcing_term;

  FunctionG function_g;


  ExactSolution UEx;

  // Initial condition.
  FunctionU0 u_0;

  // Final time.
  const double T;

  double       time      = 0.0;

  // Discretization. ///////////////////////////////////////////////////////////

  // N+1 is the number of elements.
  const unsigned int N;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Theta parameter of the theta method.
  const double theta;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Mass matrix M / deltat.
  TrilinosWrappers::SparseMatrix mass_matrix;

  // Stiffness matrix A.
  TrilinosWrappers::SparseMatrix stiffness_matrix;

  // Matrix on the left-hand side (M / deltat + theta A).
  TrilinosWrappers::SparseMatrix lhs_matrix;

  // Matrix on the right-hand side (M / deltat - (1 - theta) A).
  TrilinosWrappers::SparseMatrix rhs_matrix;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;
};

#endif