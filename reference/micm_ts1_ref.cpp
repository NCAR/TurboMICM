// MICM reference run of the TS1 benchmark for the MICM-to-C experiment.
//
//   micm_ts1_ref [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV
//   micm_ts1_ref --write-defaults FILE
//
// Builds TS1 from ROOT/micm/test/benchmark/ts1_mechanism.hpp with MICM's
// three-stage Rosenbrock solver in micm_bench's default configuration (standard
// ordering, in-place Mozart LU). Every cell starts from TS1's InitState and then
// takes the values in its CELLS_CSV row: temperature, pressure, species
// concentrations, and rate parameters. CELLS uses the first CELLS rows.
//
// One untimed warm-up step is followed by STEPS timed steps of DT seconds. The
// program prints elapsed_ms and the solver statistics summed over the timed
// steps, and writes the final state to ts1_final.csv in the working directory.
//
// Cells are solved in groups, each group a MICM state with its own step-size
// control. By default there is one group per thread: with N threads (default:
// all cores) the cells are split into N contiguous groups. With --group-size G
// the cells are split into groups of G, and each thread takes a contiguous block
// of groups, as an OpenMP static schedule does, solving each group's STEPS steps
// in turn with one solver per thread.
//
// --write-defaults writes TS1's InitState for one cell, as a CSV with the same
// columns that make_cells_csv.py writes.

#include "ts1_mechanism.hpp"

#include <algorithm>
#include <barrier>
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{
  using DenseMatrix = micm::Matrix<micm::Real>;
  using SparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>;
  using SolverBuilder = micm::CpuSolverBuilderInPlace<
      micm::RosenbrockSolverParameters,
      DenseMatrix,
      SparseMatrix,
      micm::LuDecompositionMozartInPlace<SparseMatrix>>;

  auto BuildSolver()
  {
    return bench::Ts1::Build(SolverBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()));
  }

  struct Table
  {
    std::vector<std::string> names_;
    std::vector<std::vector<double>> rows_;
  };

  std::vector<std::string> SplitCsvLine(std::string line)
  {
    if (!line.empty() && line.back() == '\r')
    {
      line.pop_back();
    }
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
    {
      fields.push_back(field);
    }
    return fields;
  }

  Table ReadCsv(const std::string& path)
  {
    std::ifstream file(path);
    if (!file)
    {
      throw std::runtime_error("cannot open " + path);
    }
    Table table;
    std::string line;
    if (!std::getline(file, line))
    {
      throw std::runtime_error(path + " is empty");
    }
    table.names_ = SplitCsvLine(line);
    while (std::getline(file, line))
    {
      if (line.empty() || line == "\r")
      {
        continue;
      }
      auto fields = SplitCsvLine(line);
      if (fields.size() != table.names_.size())
      {
        throw std::runtime_error(path + ": row " + std::to_string(table.rows_.size() + 1) + " has " +
                                 std::to_string(fields.size()) + " fields, the header has " +
                                 std::to_string(table.names_.size()));
      }
      std::vector<double> row;
      row.reserve(fields.size());
      for (const auto& value : fields)
      {
        row.push_back(std::stod(value));
      }
      table.rows_.push_back(std::move(row));
    }
    return table;
  }

  /// @brief Names indexed by position, from a name-to-index map.
  template<class Map>
  std::vector<std::string> NamesByIndex(const Map& map)
  {
    std::vector<std::string> names(map.size());
    for (const auto& [name, index] : map)
    {
      names.at(index) = name;
    }
    return names;
  }

  /// @brief Throw unless every column names temperature, pressure, a species, or a rate parameter.
  template<class State>
  void CheckColumns(const State& state, const Table& table)
  {
    for (const auto& name : table.names_)
    {
      if (name != "temperature" && name != "pressure" && !state.variable_map_.contains(name) &&
          !state.custom_rate_parameter_map_.contains(name))
      {
        throw std::runtime_error("CELLS_CSV column '" + name + "' is not a TS1 species or rate parameter");
      }
    }
  }

  /// @brief TS1's InitState, then the CSV rows from `first` on, one per cell of the state.
  template<class State>
  void Initialize(State& state, micm::Index num_cells, const Table& table, std::size_t first)
  {
    bench::Ts1::InitState(state, num_cells);
    for (std::size_t column = 0; column < table.names_.size(); ++column)
    {
      const auto& name = table.names_[column];
      for (micm::Index cell = 0; cell < num_cells; ++cell)
      {
        const double value = table.rows_[first + cell][column];
        if (name == "temperature")
        {
          state.conditions_[cell].temperature_ = value;
        }
        else if (name == "pressure")
        {
          state.conditions_[cell].pressure_ = value;
        }
        else if (auto species = state.variable_map_.find(name); species != state.variable_map_.end())
        {
          state.variables_[cell][species->second] = value;
        }
        else
        {
          state.custom_rate_parameters_[cell][state.custom_rate_parameter_map_.at(name)] = value;
        }
      }
    }
    for (micm::Index cell = 0; cell < num_cells; ++cell)
    {
      state.conditions_[cell].CalculateIdealAirDensity();
    }
  }

  struct Statistics
  {
    std::size_t internal_steps_{};
    std::size_t accepted_{};
    std::size_t rejected_{};
    std::size_t function_calls_{};
    std::size_t jacobian_updates_{};
    std::size_t decompositions_{};
    std::size_t solves_{};

    void Add(const micm::SolverStats& stats)
    {
      internal_steps_ += stats.number_of_steps_;
      accepted_ += stats.accepted_;
      rejected_ += stats.rejected_;
      function_calls_ += stats.function_calls_;
      jacobian_updates_ += stats.jacobian_updates_;
      decompositions_ += stats.decompositions_;
      solves_ += stats.solves_;
    }

    void Add(const Statistics& other)
    {
      internal_steps_ += other.internal_steps_;
      accepted_ += other.accepted_;
      rejected_ += other.rejected_;
      function_calls_ += other.function_calls_;
      jacobian_updates_ += other.jacobian_updates_;
      decompositions_ += other.decompositions_;
      solves_ += other.solves_;
    }
  };

  /// @brief One group of cells: one MICM state.
  struct Group
  {
    std::size_t first_cell_{};
    micm::Index num_cells_{};
    std::vector<double> final_state_{};  // cell-major: num_cells_ x species
  };

  /// @brief One thread and its contiguous block of groups.
  struct Worker
  {
    std::size_t first_group_{};
    std::size_t num_groups_{};
    Statistics statistics_{};
    std::string error_{};
  };

  std::string Format(double value)
  {
    char buffer[32];
    std::snprintf(buffer, sizeof buffer, "%.17g", value);
    return buffer;
  }

  int WriteDefaults(const std::string& path)
  {
    auto solver = BuildSolver();
    auto state = solver.GetState(1);
    bench::Ts1::InitState(state, 1);
    const auto species = NamesByIndex(state.variable_map_);
    const auto parameters = NamesByIndex(state.custom_rate_parameter_map_);
    std::ofstream file(path);
    file << "temperature,pressure";
    for (const auto& name : species)
    {
      file << ',' << name;
    }
    for (const auto& name : parameters)
    {
      file << ',' << name;
    }
    file << '\n' << Format(state.conditions_[0].temperature_) << ',' << Format(state.conditions_[0].pressure_);
    for (std::size_t i = 0; i < species.size(); ++i)
    {
      file << ',' << Format(state.variables_[0][i]);
    }
    for (std::size_t i = 0; i < parameters.size(); ++i)
    {
      file << ',' << Format(state.custom_rate_parameters_[0][i]);
    }
    file << '\n';
    std::cout << "wrote " << path << ": " << species.size() << " species, " << parameters.size() << " rate parameters"
              << std::endl;
    return 0;
  }

  constexpr const char* kUsage =
      "usage: micm_ts1_ref [--threads N] [--group-size G] CELLS STEPS DT CELLS_CSV\n"
      "       micm_ts1_ref --write-defaults FILE";
}  // namespace

int main(int argc, char** argv)
{
  try
  {
    std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
    std::size_t group_size = 0;  // 0: one group per thread
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg = argv[i];
      if (arg == "--write-defaults" && i + 1 < argc)
      {
        return WriteDefaults(argv[i + 1]);
      }
      if (arg == "--threads" && i + 1 < argc)
      {
        threads = std::stoul(argv[++i]);
      }
      else if (arg == "--group-size" && i + 1 < argc)
      {
        group_size = std::stoul(argv[++i]);
        if (group_size < 1)
        {
          std::cerr << kUsage << std::endl;
          return 2;
        }
      }
      else
      {
        positional.push_back(arg);
      }
    }
    if (positional.size() != 4 || threads < 1)
    {
      std::cerr << kUsage << std::endl;
      return 2;
    }
    const std::size_t num_cells = std::stoul(positional[0]);
    const std::size_t num_steps = std::stoul(positional[1]);
    const double dt = std::stod(positional[2]);
    const Table table = ReadCsv(positional[3]);
    if (num_cells < 1 || table.rows_.size() < num_cells)
    {
      throw std::runtime_error(
          "CELLS is " + std::to_string(num_cells) + " but CELLS_CSV has " + std::to_string(table.rows_.size()) + " rows");
    }

    std::vector<std::string> species;
    {
      auto solver = BuildSolver();
      auto state = solver.GetState(1);
      CheckColumns(state, table);
      species = NamesByIndex(state.variable_map_);
    }

    // Groups of cells: one per thread by default, or groups of group_size.
    std::vector<Group> groups;
    const std::size_t num_groups = group_size ? (num_cells + group_size - 1) / group_size : std::min(threads, num_cells);
    for (std::size_t g = 0, first = 0; g < num_groups; ++g)
    {
      const std::size_t count =
          group_size ? std::min(group_size, num_cells - first) : num_cells / num_groups + (g < num_cells % num_groups ? 1 : 0);
      groups.push_back(Group{ first, count, {} });
      first += count;
    }

    // Each thread takes a contiguous block of groups, as an OpenMP static schedule does.
    threads = std::min(threads, num_groups);
    std::vector<Worker> workers(threads);
    for (std::size_t t = 0, first = 0; t < threads; ++t)
    {
      workers[t].first_group_ = first;
      workers[t].num_groups_ = num_groups / threads + (t < num_groups % threads ? 1 : 0);
      first += workers[t].num_groups_;
    }

    // The main thread and every worker meet here twice: after the warm-up, and
    // after the timed steps. The time between the two meetings is elapsed_ms.
    std::barrier sync(static_cast<std::ptrdiff_t>(threads + 1));
    auto run_worker = [&](Worker& worker)
    {
      try
      {
        auto solver = BuildSolver();
        std::vector<decltype(solver.GetState(1))> states;
        states.reserve(worker.num_groups_);
        for (std::size_t k = 0; k < worker.num_groups_ && worker.error_.empty(); ++k)
        {
          const auto& group = groups[worker.first_group_ + k];
          states.push_back(solver.GetState(group.num_cells_));
          auto& state = states.back();
          Initialize(state, group.num_cells_, table, group.first_cell_);
          solver.UpdateStateParameters(state);
          auto warmup = solver.Solve(dt, state);
          if (warmup.state_ != micm::SolverState::Converged)
          {
            worker.error_ = "warm-up solve: " + micm::SolverStateToString(warmup.state_);
          }
        }
        sync.arrive_and_wait();
        for (std::size_t k = 0; k < states.size() && worker.error_.empty(); ++k)
        {
          for (std::size_t step = 0; step < num_steps && worker.error_.empty(); ++step)
          {
            auto result = solver.Solve(dt, states[k]);
            worker.statistics_.Add(result.stats_);
            if (result.state_ != micm::SolverState::Converged)
            {
              worker.error_ = "timed step " + std::to_string(step) + ": " + micm::SolverStateToString(result.state_);
            }
          }
        }
        sync.arrive_and_wait();
        for (std::size_t k = 0; k < states.size(); ++k)
        {
          auto& group = groups[worker.first_group_ + k];
          group.final_state_.resize(group.num_cells_ * species.size());
          for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
          {
            for (std::size_t i = 0; i < species.size(); ++i)
            {
              group.final_state_[cell * species.size() + i] = states[k].variables_[cell][i];
            }
          }
        }
      }
      catch (const std::exception& e)
      {
        // Keep the meetings balanced so the other threads can finish.
        worker.error_ = e.what();
        sync.arrive_and_drop();
      }
    };

    std::vector<std::thread> threads_running;
    for (auto& worker : workers)
    {
      threads_running.emplace_back(run_worker, std::ref(worker));
    }
    sync.arrive_and_wait();
    const auto t0 = std::chrono::steady_clock::now();
    sync.arrive_and_wait();
    const auto t1 = std::chrono::steady_clock::now();
    for (auto& thread : threads_running)
    {
      thread.join();
    }

    Statistics statistics;
    for (const auto& worker : workers)
    {
      if (!worker.error_.empty())
      {
        const auto& first = groups[worker.first_group_];
        const auto& last = groups[worker.first_group_ + worker.num_groups_ - 1];
        std::cerr << "cells " << first.first_cell_ << "-" << last.first_cell_ + last.num_cells_ - 1
                  << " did not converge: " << worker.error_ << std::endl;
        return 1;
      }
      statistics.Add(worker.statistics_);
    }

    std::ofstream file("ts1_final.csv");
    file << "temperature,pressure";
    for (const auto& name : species)
    {
      file << ',' << name;
    }
    file << '\n';
    for (const auto& group : groups)
    {
      for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
      {
        const auto& row = table.rows_[group.first_cell_ + cell];
        double temperature = 287.45, pressure = 101319.9;  // InitState's conditions
        for (std::size_t column = 0; column < table.names_.size(); ++column)
        {
          if (table.names_[column] == "temperature")
          {
            temperature = row[column];
          }
          else if (table.names_[column] == "pressure")
          {
            pressure = row[column];
          }
        }
        file << Format(temperature) << ',' << Format(pressure);
        for (std::size_t i = 0; i < species.size(); ++i)
        {
          file << ',' << Format(group.final_state_[cell * species.size() + i]);
        }
        file << '\n';
      }
    }

    const double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << std::fixed << std::setprecision(3) << "elapsed_ms=" << elapsed_ms << '\n'
              << "internal_steps=" << statistics.internal_steps_ << " accepted=" << statistics.accepted_
              << " rejected=" << statistics.rejected_ << " function_calls=" << statistics.function_calls_
              << " jacobian_updates=" << statistics.jacobian_updates_
              << " decompositions=" << statistics.decompositions_ << " solves=" << statistics.solves_ << std::endl;
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "micm_ts1_ref: " << e.what() << std::endl;
    return 1;
  }
}
