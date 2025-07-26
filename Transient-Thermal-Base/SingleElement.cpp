#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <filesystem>

using namespace dealii;

class TimeDependentTemperature : public Function<3>
{
public:
    TimeDependentTemperature(const std::string &filename)
        : Function<3>(), current_time(0.0)
    {
        std::ifstream input(filename);
        AssertThrow(input, ExcMessage("Could not open temperature input file"));

        std::string line;
        while (std::getline(input, line))
        {
            if (line.empty() || line[0] == '#')
                continue;  // Skip header or empty lines

            std::istringstream iss(line);
            double t, temp;
            if (!(iss >> t >> temp))
                continue;  // Skip bad lines

            temperature_data.emplace_back(t, temp+273.15);
        }

        AssertThrow(temperature_data.size() >= 2,
                    ExcMessage("Temperature data must contain at least two points"));

        // Sort the data by time to ensure lower_bound works correctly
        std::sort(temperature_data.begin(), temperature_data.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
    }

    void update_time(const double t)
    {
        current_time = t;
    }

    virtual double value(const Point<3> &, const unsigned int component = 0) const override
    {
        AssertIndexRange(component, 1);
        return interpolate_temperature(current_time);
    }

private:
    double interpolate_temperature(double t) const
    {
        // Find the first element with time >= t
        auto it = std::lower_bound(temperature_data.begin(), temperature_data.end(), std::make_pair(t, 0.0),
                                   [](const auto &a, const auto &b) { return a.first < b.first; });

        if (it == temperature_data.begin())
            return it->second;  // t before data range, clamp to first value

        if (it == temperature_data.end())
            return (it - 1)->second;  // t after data range, clamp to last value

        auto [t2, T2] = *it;
        auto [t1, T1] = *(it - 1);

        double alpha = (t - t1) / (t2 - t1);
        return T1 * (1 - alpha) + T2 * alpha;
    }

    std::vector<std::pair<double, double>> temperature_data;
    double current_time;
};

class TransientThermalProblem
{
public:
	TransientThermalProblem();
	void run();

private:
	void setup_system();
	void assemble_system();
	void solve_time_step();
	void load_material_tables();
	double get_cp(double T);
	double get_k(double T);
	void output_results(unsigned int timestep) const;

	Triangulation<3> triangulation;
	FE_Q<3>		  fe;
	DoFHandler<3>	dof_handler;

	SparsityPattern	  sparsity_pattern;
	SparseMatrix<double> mass_matrix, stiffness_matrix, system_matrix;
	Vector<double>	   solution, old_solution, system_rhs;

	double time, time_step;
	unsigned int timestep_number;
	const double total_time = 555.0;

	std::vector<std::pair<double, double>> cp_table, k_table;
	const double density = 4620; // example: steel
};

double linear_interpolate(const std::vector<std::pair<double, double>> &table, double T)
{
	Assert(table.size() >= 2, ExcInternalError());
	auto it = std::lower_bound(table.begin(), table.end(), std::make_pair(T, 0.0),
							   [](const auto &a, const auto &b) { return a.first < b.first; });

	if (it == table.begin())
		return it->second;
	if (it == table.end())
		return (table.end() - 1)->second;

	auto [T1, val1] = *(it - 1);
	auto [T2, val2] = *it;

	double alpha = (T - T1) / (T2 - T1);
	return val1 * (1 - alpha) + val2 * alpha;
}

TransientThermalProblem::TransientThermalProblem()
	: fe(2), dof_handler(triangulation), time(0), time_step(0.1), timestep_number(0) {}

void TransientThermalProblem::setup_system()
{
	GridGenerator::subdivided_hyper_rectangle(triangulation, {2,2,2},
		Point<3>(0,0,0), Point<3>(1,1,1));
	
	dof_handler.distribute_dofs(fe);
	solution.reinit(dof_handler.n_dofs());
	old_solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	mass_matrix.reinit(sparsity_pattern);
	stiffness_matrix.reinit(sparsity_pattern);
	system_matrix.reinit(sparsity_pattern);
}

void TransientThermalProblem::load_material_tables()
{
	auto load_table = [](const std::string &filename) {
		std::ifstream file(filename);
		AssertThrow(file, ExcMessage("Can't open " + filename));
		std::vector<std::pair<double, double>> table;
		double T, value;
		while (file >> T >> value)
			table.emplace_back(T, value);
		std::sort(table.begin(), table.end());
		return table;
	};

	cp_table = load_table("cp_table.dat");
	k_table  = load_table("k_table.dat");
}

double TransientThermalProblem::get_cp(double T)
{
	return linear_interpolate(cp_table, T);
}

double TransientThermalProblem::get_k(double T)
{
	return linear_interpolate(k_table, T);
}

void TransientThermalProblem::assemble_system()
{
	mass_matrix = 0;
	stiffness_matrix = 0;

	QGauss<3> quadrature_formula(fe.degree + 1);
	FEValues<3> fe_values(fe, quadrature_formula,
						  update_values | update_gradients |
						  update_JxW_values);

	const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
	const unsigned int n_q_points	= quadrature_formula.size();

	FullMatrix<double> cell_mass(dofs_per_cell), cell_stiffness(dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		cell_mass = 0;
		cell_stiffness = 0;

		double T_avg = 0;
		for (unsigned int q = 0; q < n_q_points; ++q)
		{
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				cell->get_dof_indices(local_dof_indices);
				for (unsigned int i = 0; i < dofs_per_cell; ++i)
					T_avg += old_solution(local_dof_indices[i]) * fe_values.shape_value(i, q);
			}
		}
		T_avg /= (n_q_points * dofs_per_cell);

		double cp = get_cp(T_avg);
		double k  = get_k(T_avg);

		for (unsigned int q = 0; q < n_q_points; ++q)
		{
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
				{
					cell_mass(i,j) += density * cp *
						fe_values.shape_value(i,q) *
						fe_values.shape_value(j,q) *
						fe_values.JxW(q);

					cell_stiffness(i,j) += k *
						fe_values.shape_grad(i,q) *
						fe_values.shape_grad(j,q) *
						fe_values.JxW(q);
				}
			}
		}

		cell->get_dof_indices(local_dof_indices);
		mass_matrix.add(local_dof_indices, cell_mass);
		stiffness_matrix.add(local_dof_indices, cell_stiffness);
	}

	// theta method (implicit Euler)
	system_matrix.copy_from(mass_matrix);
	system_matrix.add(time_step, stiffness_matrix);

	// RHS
	system_rhs = 0;
	mass_matrix.vmult(system_rhs, old_solution);
}

void TransientThermalProblem::solve_time_step()
{
	SparseDirectUMFPACK solver;
	solver.initialize(system_matrix);
	solver.vmult(solution, system_rhs);
}

void TransientThermalProblem::output_results(unsigned int timestep) const
{
	const std::string output_dir = "Solution";
	// Create directory if it doesn't exist
    std::filesystem::create_directory(output_dir);

    DataOut<3> data_out;
    data_out.attach_dof_handler(dof_handler);
	if (timestep == 0) 
		data_out.add_data_vector(old_solution, "temperature");
	else
		data_out.add_data_vector(solution, "temperature");
    data_out.build_patches();

    const std::string filename = output_dir + "/solution-" + std::to_string(timestep);
    std::ofstream out(filename + ".vtu");
    data_out.write_vtu(out);
    std::ofstream gnuplot_file(filename + ".gpl");
    data_out.write_gnuplot(gnuplot_file);
}

void TransientThermalProblem::run()
{
    setup_system();
    load_material_tables();

    TimeDependentTemperature boundary_temp("temperature_input.dat");

    // Set initial temperature to 293.15 K everywhere
    VectorTools::interpolate(dof_handler, Functions::ConstantFunction<3>(293.15), old_solution);

    timestep_number = 0;
    time = 0.0;

    // Output initial condition
    output_results(timestep_number);  // Make sure it outputs old_solution, not solution

    while (time < total_time)
    {
        time += time_step;
        ++timestep_number;

        assemble_system();

        boundary_temp.update_time(time);

        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,  // boundary_id
                                                 boundary_temp,
                                                 boundary_values);

        MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);

        solve_time_step();

        old_solution = solution;

        output_results(timestep_number);

        std::cout << "Time step " << timestep_number << " at time " << time << std::endl;
    }
}

int main()
{
    TransientThermalProblem thermal_problem;
    thermal_problem.run();
    return 0;
}
