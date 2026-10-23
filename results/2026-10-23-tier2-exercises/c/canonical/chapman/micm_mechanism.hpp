// GENERATED FILE by TurboMICM admin/micm_header.py; do not edit.
// chapman: MUSICA configs/v0/chapman, resolved by admin/v0_to_canonical.py.
// 5 species and 7 reactions.
#pragma once

#include <micm/CPU.hpp>

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace bench
{
  namespace ts1
  {
    inline const micm::Species& SpeciesByName(const micm::Phase& phase, std::string_view name)
    {
      for (const auto& phase_species : phase.phase_species_)
        if (phase_species.species_.name_ == name)
          return phase_species.species_;
      throw std::runtime_error("no species named " + std::string(name));
    }

    inline micm::Phase CreateGasPhase()
    {
      std::vector<micm::PhaseSpecies> phase_species;
      {
        auto species = micm::Species("M");
        species.SetThirdBody();
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O1D");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O3");
        phase_species.emplace_back(species);
      }
      return micm::Phase{ "gas", phase_species };
    }

    inline std::vector<micm::Process> CreateProcesses(const micm::Phase& gas_phase)
    {
      std::vector<micm::Process> processes;
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 2.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_jO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 48.28552461368, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_jO3_to_O", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 949.089383776, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_jO3_to_O1D", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 42817420.803600006, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 72265689.12, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      return processes;
    }
  }  // namespace ts1

  struct Ts1
  {
    static constexpr std::string_view kName = "chapman";

    template<class Builder>
    static auto Build(Builder builder)
    {
      auto gas_phase = ts1::CreateGasPhase();
      return builder.SetSystem(micm::System(gas_phase))
          .SetReactions(ts1::CreateProcesses(gas_phase))
          .SetIgnoreUnusedSpecies(true)
          .Build();
    }

    template<class State>
    static void InitState(State& state, micm::Index num_cells)
    {
      std::vector<micm::Real> cells(num_cells);
      auto concentration = [&](const char* name, micm::Real value)
      {
        std::fill(cells.begin(), cells.end(), value);
        state[name] = cells;
      };
      auto parameter = [&](const char* label, micm::Real value)
      {
        std::fill(cells.begin(), cells.end(), value);
        state.SetCustomRateParameter(label, cells);
      };
      concentration("O2", 0.75);
      concentration("O", 0.0);
      concentration("O1D", 0.0);
      concentration("O3", 8.1e-06);
      parameter("PHOTO_jO2", 2.42e-17);
      parameter("PHOTO_jO3_to_O", 1.15e-05);
      parameter("PHOTO_jO3_to_O1D", 6.61e-09);
      for (micm::Index c = 0; c < num_cells; ++c)
      {
        state.conditions_[c].temperature_ = 206.6374207;
        state.conditions_[c].pressure_ = 6152.049805;
        state.conditions_[c].CalculateIdealAirDensity();
      }
    }
  };
}  // namespace bench
