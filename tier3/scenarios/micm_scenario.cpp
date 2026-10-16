// MICM run of a Tier 3 TS1 scenario, for reference solutions and baselines.
//
//   micm_scenario [--threads N] [--group-size G] [--method NAME] [--rtol R]
//                 [--atol A] [--max-steps M] [--no-output] SCENARIO
//
// SCENARIO is a scenario file (see README.md): a CELLS_CSV, a step length DT, a
// number of steps, a photolysis-factor file f(t), and the snapshot interval. Every
// cell starts from TS1's InitState and then takes its CELLS_CSV row, as in Tier 1.
// At the start of each step, each cell's photolysis rates (the rate parameters
// whose labels begin with "j") are set to its CELLS_CSV values times f at that
// time, and the rate constants are updated; the step is then solved with MICM.
//
// --method picks MICM's Rosenbrock parameters: rosenbrock2, rosenbrock3 (the
// default, Tier 1's solver), rosenbrock4, rosenbrock4dae, or rosenbrock6dae.
// --rtol and --atol override MICM's default tolerances (relative 1e-6; absolute
// 1e-3 mol m-3 for every TS1 species). --max-steps raises MICM's limit of 1000
// internal steps per step. Cells are solved in groups as micm_ts1_ref does: one
// group per thread by default, or groups of G cells, each thread taking a
// contiguous block of groups.
//
// The program times all the steps, including the rate updates, but not reading
// the input, setup, or writing output, and prints elapsed_ms and the solver
// statistics summed over all steps and groups. Unless --no-output is given, it
// writes ts1_final.csv and, if the scenario asks for snapshots, ts1_snapshots.csv.

#include "ts1_mechanism.hpp"

#include <algorithm>
#include <barrier>
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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

  micm::RosenbrockSolverParameters Parameters(const std::string& method, std::size_t max_steps)
  {
    const auto pick = [&]
    {
      using P = micm::RosenbrockSolverParameters;
      if (method == "rosenbrock2")
        return P::TwoStageRosenbrockParameters();
      if (method == "rosenbrock3")
        return P::ThreeStageRosenbrockParameters();
      if (method == "rosenbrock4")
        return P::FourStageRosenbrockParameters();
      if (method == "rosenbrock4dae")
        return P::FourStageDifferentialAlgebraicRosenbrockParameters();
      if (method == "rosenbrock6dae")
        return P::SixStageDifferentialAlgebraicRosenbrockParameters();
      throw std::runtime_error("unknown --method " + method);
    };
    auto parameters = pick();
    if (max_steps > 0)
    {
      parameters.max_number_of_steps_ = max_steps;
    }
    return parameters;
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
        throw std::runtime_error(path + ": a row has " + std::to_string(fields.size()) + " fields, the header has " +
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

  std::string Trim(const std::string& text)
  {
    const auto first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos)
    {
      return "";
    }
    const auto last = text.find_last_not_of(" \t\r");
    return text.substr(first, last - first + 1);
  }

  /// @brief A scenario file: `key = value` lines; `#` starts a comment.
  struct Scenario
  {
    std::string cells_;
    double dt_{};
    std::size_t steps_{};
    std::string photolysis_;
    std::size_t snapshot_every_{};
  };

  Scenario ReadScenario(const std::string& path)
  {
    std::ifstream file(path);
    if (!file)
    {
      throw std::runtime_error("cannot open " + path);
    }
    const auto slash = path.find_last_of('/');
    const std::string directory = slash == std::string::npos ? "" : path.substr(0, slash + 1);
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(file, line))
    {
      line = Trim(line.substr(0, line.find('#')));
      if (line.empty())
      {
        continue;
      }
      const auto equals = line.find('=');
      if (equals == std::string::npos)
      {
        throw std::runtime_error(path + ": expected key = value, got: " + line);
      }
      values[Trim(line.substr(0, equals))] = Trim(line.substr(equals + 1));
    }
    for (const char* key : { "cells", "dt", "steps", "photolysis" })
    {
      if (!values.count(key))
      {
        throw std::runtime_error(path + ": missing " + key);
      }
    }
    Scenario scenario;
    scenario.cells_ = directory + values["cells"];
    scenario.dt_ = std::stod(values["dt"]);
    scenario.steps_ = std::stoul(values["steps"]);
    scenario.photolysis_ = directory + values["photolysis"];
    scenario.snapshot_every_ = values.count("snapshot_every") ? std::stoul(values["snapshot_every"]) : 0;
    return scenario;
  }

  /// @brief f(t), piecewise linear between the points, constant beyond them.
  struct Photolysis
  {
    std::vector<double> times_, factors_;

    double operator()(double time) const
    {
      if (time <= times_.front())
        return factors_.front();
      if (time >= times_.back())
        return factors_.back();
      const auto upper = std::upper_bound(times_.begin(), times_.end(), time) - times_.begin();
      const double t0 = times_[upper - 1], t1 = times_[upper];
      const double w = (time - t0) / (t1 - t0);
      return (1.0 - w) * factors_[upper - 1] + w * factors_[upper];
    }
  };

  Photolysis ReadPhotolysis(const std::string& path)
  {
    const Table table = ReadCsv(path);
    if (table.names_.size() != 2 || table.rows_.empty())
    {
      throw std::runtime_error(path + ": expected a header and rows of time_s,factor");
    }
    Photolysis photolysis;
    for (const auto& row : table.rows_)
    {
      if (!photolysis.times_.empty() && row[0] <= photolysis.times_.back())
      {
        throw std::runtime_error(path + ": times must increase");
      }
      photolysis.times_.push_back(row[0]);
      photolysis.factors_.push_back(row[1]);
    }
    return photolysis;
  }

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
          state.conditions_[cell].temperature_ = value;
        else if (name == "pressure")
          state.conditions_[cell].pressure_ = value;
        else if (auto species = state.variable_map_.find(name); species != state.variable_map_.end())
          state.variables_[cell][species->second] = value;
        else
          state.custom_rate_parameters_[cell][state.custom_rate_parameter_map_.at(name)] = value;
      }
    }
    for (micm::Index cell = 0; cell < num_cells; ++cell)
    {
      state.conditions_[cell].CalculateIdealAirDensity();
    }
  }

  struct Statistics
  {
    std::size_t internal_steps_{}, accepted_{}, rejected_{}, function_calls_{}, jacobian_updates_{}, decompositions_{},
        solves_{};

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

  struct Group
  {
    std::size_t first_cell_{};
    micm::Index num_cells_{};
    std::vector<std::vector<double>> snapshots_{};  // per snapshot: num_cells_ x species
  };

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

  constexpr const char* kUsage =
      "usage: micm_scenario [--threads N] [--group-size G] [--method NAME] [--rtol R] [--atol A]\n"
      "                     [--max-steps M] [--no-output] SCENARIO";
}  // namespace

int main(int argc, char** argv)
{
  try
  {
    std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
    std::size_t group_size = 0, max_steps = 0;
    std::string method = "rosenbrock3";
    double rtol = -1.0, atol = -1.0;
    bool output = true;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg = argv[i];
      const bool has_value = i + 1 < argc;
      if (arg == "--threads" && has_value)
        threads = std::stoul(argv[++i]);
      else if (arg == "--group-size" && has_value)
        group_size = std::stoul(argv[++i]);
      else if (arg == "--method" && has_value)
        method = argv[++i];
      else if (arg == "--rtol" && has_value)
        rtol = std::stod(argv[++i]);
      else if (arg == "--atol" && has_value)
        atol = std::stod(argv[++i]);
      else if (arg == "--max-steps" && has_value)
        max_steps = std::stoul(argv[++i]);
      else if (arg == "--no-output")
        output = false;
      else
        positional.push_back(arg);
    }
    if (positional.size() != 1 || threads < 1)
    {
      std::cerr << kUsage << std::endl;
      return 2;
    }
    const Scenario scenario = ReadScenario(positional[0]);
    const Table table = ReadCsv(scenario.cells_);
    const Photolysis photolysis = ReadPhotolysis(scenario.photolysis_);
    const std::size_t num_cells = table.rows_.size();
    if (num_cells < 1)
    {
      throw std::runtime_error(scenario.cells_ + " has no cells");
    }
    const auto parameters = Parameters(method, max_steps);
    const auto build = [&] { return bench::Ts1::Build(SolverBuilder(parameters)); };

    std::vector<std::string> species;
    std::vector<micm::Index> photolysis_index;  // custom rate parameters whose labels begin with "j"
    {
      auto solver = build();
      auto state = solver.GetState(1);
      CheckColumns(state, table);
      species = NamesByIndex(state.variable_map_);
      for (const auto& [label, index] : state.custom_rate_parameter_map_)
      {
        if (!label.empty() && label[0] == 'j')
        {
          photolysis_index.push_back(index);
        }
      }
      std::sort(photolysis_index.begin(), photolysis_index.end());
    }

    std::vector<Group> groups;
    const std::size_t num_groups = group_size ? (num_cells + group_size - 1) / group_size : std::min(threads, num_cells);
    for (std::size_t g = 0, first = 0; g < num_groups; ++g)
    {
      const std::size_t count =
          group_size ? std::min(group_size, num_cells - first) : num_cells / num_groups + (g < num_cells % num_groups ? 1 : 0);
      groups.push_back(Group{ first, count, {} });
      first += count;
    }
    threads = std::min(threads, num_groups);
    std::vector<Worker> workers(threads);
    for (std::size_t t = 0, first = 0; t < threads; ++t)
    {
      workers[t].first_group_ = first;
      workers[t].num_groups_ = num_groups / threads + (t < num_groups % threads ? 1 : 0);
      first += workers[t].num_groups_;
    }

    // Snapshots are taken after every snapshot_every steps; the final state is always kept.
    std::vector<std::size_t> snapshot_steps;
    if (scenario.snapshot_every_ > 0)
    {
      for (std::size_t step = scenario.snapshot_every_; step <= scenario.steps_; step += scenario.snapshot_every_)
      {
        snapshot_steps.push_back(step);
      }
    }
    if (snapshot_steps.empty() || snapshot_steps.back() != scenario.steps_)
    {
      snapshot_steps.push_back(scenario.steps_);
    }

    std::barrier sync(static_cast<std::ptrdiff_t>(threads + 1));
    auto run_worker = [&](Worker& worker)
    {
      try
      {
        auto solver = build();
        std::vector<decltype(solver.GetState(1))> states;
        std::vector<std::vector<double>> base_photolysis;  // per group: cell-major, photolysis_index order
        states.reserve(worker.num_groups_);
        for (std::size_t k = 0; k < worker.num_groups_; ++k)
        {
          const auto& group = groups[worker.first_group_ + k];
          states.push_back(solver.GetState(group.num_cells_));
          auto& state = states.back();
          Initialize(state, group.num_cells_, table, group.first_cell_);
          if (rtol > 0)
            state.SetRelativeTolerance(rtol);
          if (atol > 0)
            state.SetAbsoluteTolerances(std::vector<micm::Real>(species.size(), atol));
          std::vector<double> base(group.num_cells_ * photolysis_index.size());
          for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
            for (std::size_t j = 0; j < photolysis_index.size(); ++j)
              base[cell * photolysis_index.size() + j] = state.custom_rate_parameters_[cell][photolysis_index[j]];
          base_photolysis.push_back(std::move(base));
        }
        sync.arrive_and_wait();
        for (std::size_t k = 0; k < states.size() && worker.error_.empty(); ++k)
        {
          auto& state = states[k];
          auto& group = groups[worker.first_group_ + k];
          const auto& base = base_photolysis[k];
          std::size_t next_snapshot = 0;
          for (std::size_t step = 0; step < scenario.steps_ && worker.error_.empty(); ++step)
          {
            const double factor = photolysis(static_cast<double>(step) * scenario.dt_);
            for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
              for (std::size_t j = 0; j < photolysis_index.size(); ++j)
                state.custom_rate_parameters_[cell][photolysis_index[j]] = base[cell * photolysis_index.size() + j] * factor;
            solver.UpdateStateParameters(state);
            auto result = solver.Solve(scenario.dt_, state);
            worker.statistics_.Add(result.stats_);
            if (result.state_ != micm::SolverState::Converged)
            {
              worker.error_ = "step " + std::to_string(step) + ": " + micm::SolverStateToString(result.state_);
            }
            if (next_snapshot < snapshot_steps.size() && step + 1 == snapshot_steps[next_snapshot])
            {
              std::vector<double> snapshot(group.num_cells_ * species.size());
              for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
                for (std::size_t i = 0; i < species.size(); ++i)
                  snapshot[cell * species.size() + i] = state.variables_[cell][i];
              group.snapshots_.push_back(std::move(snapshot));
              ++next_snapshot;
            }
          }
        }
        sync.arrive_and_wait();
      }
      catch (const std::exception& e)
      {
        worker.error_ = e.what();
        sync.arrive_and_drop();
      }
    };

    std::vector<std::thread> running;
    for (auto& worker : workers)
    {
      running.emplace_back(run_worker, std::ref(worker));
    }
    sync.arrive_and_wait();
    const auto t0 = std::chrono::steady_clock::now();
    sync.arrive_and_wait();
    const auto t1 = std::chrono::steady_clock::now();
    for (auto& thread : running)
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
                  << " failed: " << worker.error_ << std::endl;
        return 1;
      }
      statistics.Add(worker.statistics_);
    }

    if (output)
    {
      std::vector<double> temperature(num_cells, 287.45), pressure(num_cells, 101319.9);  // InitState's conditions
      for (std::size_t column = 0; column < table.names_.size(); ++column)
      {
        for (std::size_t cell = 0; cell < num_cells; ++cell)
        {
          if (table.names_[column] == "temperature")
            temperature[cell] = table.rows_[cell][column];
          else if (table.names_[column] == "pressure")
            pressure[cell] = table.rows_[cell][column];
        }
      }
      auto write = [&](std::ostream& file, std::size_t snapshot, bool with_step)
      {
        for (const auto& group : groups)
        {
          for (micm::Index cell = 0; cell < group.num_cells_; ++cell)
          {
            const std::size_t c = group.first_cell_ + cell;
            if (with_step)
              file << snapshot_steps[snapshot] << ',' << c << ',';
            file << Format(temperature[c]) << ',' << Format(pressure[c]);
            for (std::size_t i = 0; i < species.size(); ++i)
              file << ',' << Format(group.snapshots_[snapshot][cell * species.size() + i]);
            file << '\n';
          }
        }
      };
      std::ofstream final_file("ts1_final.csv");
      final_file << "temperature,pressure";
      for (const auto& name : species)
        final_file << ',' << name;
      final_file << '\n';
      write(final_file, snapshot_steps.size() - 1, false);
      if (scenario.snapshot_every_ > 0)
      {
        std::ofstream snapshot_file("ts1_snapshots.csv");
        snapshot_file << "step,cell,temperature,pressure";
        for (const auto& name : species)
          snapshot_file << ',' << name;
        snapshot_file << '\n';
        for (std::size_t s = 0; s < snapshot_steps.size(); ++s)
          write(snapshot_file, s, true);
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
    std::cerr << "micm_scenario: " << e.what() << std::endl;
    return 1;
  }
}
