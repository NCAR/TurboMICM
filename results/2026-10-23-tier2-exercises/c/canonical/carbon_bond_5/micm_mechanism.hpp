// GENERATED FILE by TurboMICM admin/micm_header.py; do not edit.
// carbon_bond_5: MUSICA configs/v0/carbon_bond_5, resolved by admin/v0_to_canonical.py.
// 67 species and 200 reactions.
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
        auto species = micm::Species("AACD");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ALD2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ALDX");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZENE");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZRO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2O3");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH4");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CL2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CLO");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CO");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CRES");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CRO");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CXO3");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ETH");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ETHA");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ETOH");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("FACD");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("FMCL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("FORM");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2O");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2O2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCO3");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HNO3");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HOCL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HONO");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("IOLE");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOP");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISPD");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("M");
        species.SetThirdBody();
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEOH");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEPX");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MGLY");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("N2O5");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO3");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NTR");
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
      {
        auto species = micm::Species("OH");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("OLE");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("OPEN");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PACD");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PAN");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PANX");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PAR");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PNA");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ROOH");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ROR");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SESQ");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SULF");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERP");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOLRO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XO2");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XO2N");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYL");
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLRO2");
        phase_species.emplace_back(species);
      }
      return micm::Phase{ "gas", phase_species };
    }

    inline std::vector<micm::Process> CreateProcesses(const micm::Phase& gas_phase)
    {
      std::vector<micm::Process> processes;
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_NO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "O2"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 217.59707599952029, .B_ = -2.4, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1806642.2280000001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3372398.8256, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 90665.44833313345, .k0_B_ = -1.8, .k0_C_ = 0.0, .kinf_A_ = 13248709.672, .kinf_B_ = -0.7, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 32639.56139992804, .k0_B_ = -1.5, .k0_C_ = 0.0, .kinf_A_ = 18066422.28, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 72265.68912, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_O3_to_O3P", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_O3_to_O1D", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 12646495.595999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 132487096.72, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1023763.9292, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6022.14076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_NO3_to_NO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_NO3_to_NO", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 9033211.14, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 27099.63342, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2O5"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 725323.5866650676, .k0_B_ = -4.4, .k0_C_ = 0.0, .kinf_A_ = 843099.7064, .kinf_B_ = -0.7, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O5"), SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.00015055351899999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O5"), SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.0006527912279985608, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 602214076000000.0, .k0_B_ = -3.5, .k0_C_ = -11000.0, .kinf_A_ = 970000000000000.0, .kinf_B_ = 0.1, .kinf_C_ = -11080.0, .Fc_ = 0.45, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.0011967839179973614, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HONO"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.00018133089666626687, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HONO"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 253863.25533277364, .k0_B_ = -2.6, .k0_C_ = 0.0, .kinf_A_ = 21679706.735999998, .kinf_B_ = -0.1, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HONO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_HONO", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HONO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 10839853.367999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HONO"), SpeciesByName(gas_phase, "HONO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.00602214076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 725323.5866650676, .k0_B_ = -3.0, .k0_C_ = 0.0, .kinf_A_ = 15055351.9, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HNO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 14453.137824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HNO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 235.73016566614695, .k0_B_ = 0.0, .k0_C_ = 1335.0, .kinf_A_ = 16.259780052, .kinf_B_ = 0.0, .kinf_C_ = 2199.0, .Fc_ = 1.0, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2107749.266, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PNA"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 65279.12279985608, .k0_B_ = -3.2, .k0_C_ = 0.0, .kinf_A_ = 2830406.1572, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PNA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 24690777116000.0, .k0_B_ = 0.0, .k0_C_ = -10650.0, .kinf_A_ = 4800000000000000.0, .kinf_B_ = 0.0, .kinf_C_ = -11170.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "PNA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 782878.2988, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 138509.23747999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 616.5250486653074, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 116.77709745307588, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.5197920875179664, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 2.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_H2O2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "H2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1746420.8204, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 66243548.36, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3312177.418, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 13248709.672, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2529299.1192, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 250236.63739944831, .k0_B_ = -1.0, .k0_C_ = 0.0, .kinf_A_ = 15657565.976, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 18066422.28, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H2O2"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 843099.7064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6022140.76, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 13248709.672, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2107749.266, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.02214076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 511881.9646, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PNA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.61 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.61 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.39 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 0.39 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_PNA", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HNO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_HNO3", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_N2O5", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1565756.5976, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2N"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1565756.5976, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ROOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2N"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ROOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NTR"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -0.66 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 355306.30484, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NTR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -0.66 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_NTR", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ROOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1812664.36876, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ROOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.5 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_ROOH", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 86718.826944, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CO"), SpeciesByName(gas_phase, "M") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1243.9299511305908, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH4") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1475424.4862, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1686199.4128, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEPX"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 246907.77116, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEO2"), SpeciesByName(gas_phase, "MEO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.37 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.74 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MEOH"), 0.63 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 57210.33722, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEPX"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.3 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2288413.4888, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEPX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_MEPX", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4396162.7548, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 5419926.683999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_FORM_to_HO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_FORM_to_CO", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 20475278.584, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 349.28416408, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FORM"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 5841.4765372, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2400000000000.0, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCO3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FACD"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3372398.8256, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEPX"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3372.3988256000002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "FACD"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 240885.6304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALD2"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 10839853.367999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALD2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3372398.8256, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALD2"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 843099.7064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALD2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_ALD2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4877934.0156, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PAN"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 97918684.19978411, .k0_B_ = -7.1, .k0_C_ = 0.0, .kinf_A_ = 7226568.9120000005, .kinf_B_ = -0.9, .kinf_C_ = 0.0, .Fc_ = 0.3, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PAN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 2950848972400000.0, .k0_B_ = 0.0, .k0_C_ = -12100.0, .kinf_A_ = 5.4e+16, .kinf_B_ = 0.0, .kinf_C_ = -13830.0, .Fc_ = 0.3, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PAN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_PAN", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PACD"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 258952.05268, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "MEO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1204428.152, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "XO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 264974.19344, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2O3"), SpeciesByName(gas_phase, "C2O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1746420.8204, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PACD"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 240885.6304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PACD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_PACD", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "AACD"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 240885.6304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALDX"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 7828782.988, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALDX"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3071291.7876, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALDX"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3914.391494, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALDX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_ALDX", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4034834.3092, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PANX"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 97918684.19978411, .k0_B_ = -7.1, .k0_C_ = 0.0, .kinf_A_ = 7226568.9120000005, .kinf_B_ = -0.9, .kinf_C_ = 0.0, .Fc_ = 0.3, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PANX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 2950848972400000.0, .k0_B_ = 0.0, .k0_C_ = -12100.0, .kinf_A_ = 5.4e+16, .kinf_B_ = 0.0, .kinf_C_ = -13830.0, .Fc_ = 0.3, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PANX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_PANX", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PANX"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 180664.2228, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PACD"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 258952.05268, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "MEO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1204428.152, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "XO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "AACD"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 264974.19344, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "CXO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1746420.8204, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CXO3"), SpeciesByName(gas_phase, "C2O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1746420.8204, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PAR"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.87 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.13 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.11 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.06 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -0.11 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ROR"), 0.76 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.05 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 487793.40155999997, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ROR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.96 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.6 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.94 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -2.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.04 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ROR"), 0.02 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1000000000000000.0, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ROR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1600.0, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ROR"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 9033211.14, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "OLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.01 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6022140.76, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.62 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.95 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -0.7 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 19270850.432, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "OLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.18 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.74 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.32 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.22 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.44 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3914.391494, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "OLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.91 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.09 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.56 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.35 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 421549.8532, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "ETH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.3 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6263026.3904, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ETH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.56 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.22 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 36266179.333253376, .k0_B_ = -0.8, .k0_C_ = 0.0, .kinf_A_ = 5299483.868799999, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "ETH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.63 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.13 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.13 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FACD"), 0.37 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 7226.568912, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "ETH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1987306.4508, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IOLE"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.24 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.66 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 13850923.748000002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IOLE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6022140.76, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IOLE"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.65 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.35 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 5058.5982384, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IOLE"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.18 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.64 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 578125.51296, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.44 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.08 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CRES"), 0.36 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "TO2"), 0.56 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLRO2"), 0.765 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1083985.3368, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OPEN"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4877934.0156, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CRES"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4.2, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CRES") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CRO"), 0.4 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.6 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.6 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OPEN"), 0.3 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 24690777.116, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CRES"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CRO"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 13248709.672, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CRO"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CRO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CRES"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3312177.418, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OPEN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_OPEN", .scaling_factor_ = 9.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OPEN"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 18066422.28, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OPEN"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.03 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 0.62 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.03 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.69 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.08 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.76 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MGLY"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 32.519560104, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XYL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CRES"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MGLY"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 1.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "TO2"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLRO2"), 0.804 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 10237639.292, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MGLY") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 10839853.367999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MGLY") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_MGLY", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 0.75 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.5 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.25 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 21679706.735999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 0.912 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.629 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.991 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.912 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.088 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 15296237.5304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 0.65 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.6 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.066 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.266 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.15 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.35 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.066 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4733.402637360001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 2.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1824708.65028, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ISPD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 1.565 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.167 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.713 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.503 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.334 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MGLY"), 0.168 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.252 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 0.21 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.12 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 20234392.9536, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "ISPD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 0.114 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.15 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MGLY"), 0.85 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.154 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.268 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.064 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.02 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.36 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.225 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4.2757199396, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "ISPD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.357 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.282 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 1.282 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.925 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.643 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 0.85 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 0.075 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.075 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 0.15 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 602.2140760000001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISPD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.333 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.067 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.832 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.033 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 0.967 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_ISPD", .scaling_factor_ = 0.0036 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.15 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 5.12 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 21679706.735999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.75 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.28 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 1.66 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.47 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 9033211.14, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.57 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.07 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.76 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.18 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.24 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.001 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 7.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.21 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 0.39 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 722.6568912, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.47 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.28 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.03 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.25 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.47 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 0.53 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2228192.0812, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SULF"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 108798.53799976013, .k0_B_ = -3.3, .k0_C_ = 0.0, .kinf_A_ = 903321.1140000001, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.6, .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ETOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.9 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.05 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 0.1 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 4155277.1244, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ETHA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.991 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.991 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.009 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 5239262.4612, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NTR"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 0.2 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.8 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 2.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 0.0903321114, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_CL2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 13850923.748000002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 9816.0894388, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3854170.0864, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "FMCL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 301107.038, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "CH4") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "MEO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3974612.9016, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "PAR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.87 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.13 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.11 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.06 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -0.11 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ROR"), 0.76 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.05 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 30110703.8, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ETHA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.991 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.991 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2N"), 0.009 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 49983768.308, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ETH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FMCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 64436906.13200001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "OLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FMCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.33 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.67 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 2.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), -1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 150553519.0, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "IOLE") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FMCL"), 0.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 0.45 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALDX"), 0.55 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "OLE"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 0.3 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.7 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 210774926.6, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ISOP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 0.15 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FMCL"), 0.85 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ISPD"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 258952052.68, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "FORM") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 49381554.232, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ALD2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "C2O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 47574912.004, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ALDX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "CXO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 78287829.88, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "MEOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 33121774.18, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "ETOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 49381554.232, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 396256.862008, .B_ = 1.16, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLRO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLRO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 114420.67444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLRO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLRO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 114420.67444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 }, micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZRO2"), 0.764 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1487468.7677200001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZRO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZRO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 114420.67444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SESQ"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6985.683281600001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SESQ"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 118636172.97199999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SESQ"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 11442067.444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 2.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "PHOTO_O2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_NO", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_NO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_CO", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_SO2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "FORM"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_FORM", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEOH"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_MEOH", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALD2"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_ALD2", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PAR"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_PAR", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ETH"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_ETH", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OLE"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_OLE", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "IOLE"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_IOLE", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TOL"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_TOL", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYL"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_XYL", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({  })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOP"), 1.0 } })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "EMIS_ISOP", .scaling_factor_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      return processes;
    }
  }  // namespace ts1

  struct Ts1
  {
    static constexpr std::string_view kName = "carbon_bond_5";

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
      concentration("AACD", 2e-09);
      concentration("ALD2", 1e-08);
      concentration("ALDX", 0.0);
      concentration("BENZENE", 0.0);
      concentration("BENZRO2", 0.0);
      concentration("C2O3", 0.0);
      concentration("CH4", 9e-05);
      concentration("CL", 0.0);
      concentration("CL2", 0.0);
      concentration("CLO", 0.0);
      concentration("CO", 6e-06);
      concentration("CRES", 0.0);
      concentration("CRO", 0.0);
      concentration("CXO3", 0.0);
      concentration("ETH", 2e-09);
      concentration("ETHA", 1e-08);
      concentration("ETOH", 0.0);
      concentration("FACD", 0.0);
      concentration("FMCL", 0.0);
      concentration("FORM", 9e-08);
      concentration("H2", 3e-05);
      concentration("H2O", 128.0);
      concentration("H2O2", 5e-08);
      concentration("HCL", 9e-08);
      concentration("HCO3", 0.0);
      concentration("HNO3", 1e-08);
      concentration("HO2", 0.0);
      concentration("HOCL", 0.0);
      concentration("HONO", 0.0);
      concentration("IOLE", 3e-11);
      concentration("ISOP", 2e-07);
      concentration("ISPD", 0.0);
      concentration("MEO2", 0.0);
      concentration("MEOH", 9e-09);
      concentration("MEPX", 2e-08);
      concentration("MGLY", 0.0);
      concentration("N2O5", 0.0);
      concentration("NO", 1e-09);
      concentration("NO2", 1e-08);
      concentration("NO3", 0.0);
      concentration("NTR", 1e-09);
      concentration("O", 0.0);
      concentration("O1D", 0.0);
      concentration("O3", 2e-06);
      concentration("OH", 0.0);
      concentration("O2", 56.0);
      concentration("OLE", 4e-10);
      concentration("OPEN", 0.0);
      concentration("PACD", 0.0);
      concentration("PAN", 3e-08);
      concentration("PANX", 0.0);
      concentration("PAR", 2e-08);
      concentration("PNA", 0.0);
      concentration("ROOH", 1e-09);
      concentration("ROR", 0.0);
      concentration("SESQ", 0.0);
      concentration("SO2", 3e-08);
      concentration("SULF", 0.0);
      concentration("TERP", 0.0);
      concentration("TO2", 0.0);
      concentration("TOL", 1e-09);
      concentration("TOLRO2", 0.0);
      concentration("XO2", 0.0);
      concentration("XO2N", 0.0);
      concentration("XYL", 1e-09);
      concentration("XYLRO2", 0.0);
      parameter("PHOTO_NO2", 0.00477);
      parameter("PHOTO_O3_to_O3P", 0.000252999999999999);
      parameter("PHOTO_O3_to_O1D", 2.26e-06);
      parameter("PHOTO_NO3_to_NO2", 0.116999999999999);
      parameter("PHOTO_NO3_to_NO", 0.0144);
      parameter("PHOTO_HONO", 0.000918);
      parameter("PHOTO_H2O2", 2.59e-06);
      parameter("PHOTO_PNA", 1.89e-06);
      parameter("PHOTO_HNO3", 8.61e-08);
      parameter("PHOTO_N2O5", 0.0);
      parameter("PHOTO_NTR", 4.77e-07);
      parameter("PHOTO_ROOH", 1.81e-06);
      parameter("PHOTO_MEPX", 1.81e-06);
      parameter("PHOTO_FORM_to_HO2", 7.93e-06);
      parameter("PHOTO_FORM_to_CO", 2.2e-05);
      parameter("PHOTO_ALD2", 2.2e-06);
      parameter("PHOTO_PAN", 0.0);
      parameter("PHOTO_PACD", 1.81e-06);
      parameter("PHOTO_ALDX", 2.2e-06);
      parameter("PHOTO_PANX", 0.0);
      parameter("PHOTO_OPEN", 0.000645);
      parameter("PHOTO_MGLY", 7.64e-05);
      parameter("PHOTO_ISPD", 1.98e-09);
      parameter("PHOTO_CL2", 0.0);
      parameter("PHOTO_O2", 0.0);
      parameter("EMIS_NO", 1.44e-10);
      parameter("EMIS_NO2", 7.56e-12);
      parameter("EMIS_CO", 1.96e-09);
      parameter("EMIS_SO2", 1.06e-09);
      parameter("EMIS_FORM", 1.02e-11);
      parameter("EMIS_MEOH", 5.92e-13);
      parameter("EMIS_ALD2", 4.25e-12);
      parameter("EMIS_PAR", 4.27e-10);
      parameter("EMIS_ETH", 4.62e-11);
      parameter("EMIS_OLE", 1.49e-11);
      parameter("EMIS_IOLE", 1.49e-11);
      parameter("EMIS_TOL", 1.53e-11);
      parameter("EMIS_XYL", 1.4e-11);
      parameter("EMIS_ISOP", 6.03e-12);
      for (micm::Index c = 0; c < num_cells; ++c)
      {
        state.conditions_[c].temperature_ = 287.45;
        state.conditions_[c].pressure_ = 101319.9;
        state.conditions_[c].CalculateIdealAirDensity();
      }
    }
  };
}  // namespace bench
