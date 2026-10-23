#!/usr/bin/env python3
"""Write a MICM mechanism header from a canonical mechanism, for MICM's reference run.

    admin/micm_header.py CANONICAL.json OUT_DIR

canonical.json comes from admin/v0_to_canonical.py. OUT_DIR/ts1_mechanism.hpp
has the interface of MICM's test/benchmark/ts1_mechanism.hpp (struct
bench::Ts1 with Build() and InitState()), so reference/micm_ts1_ref.cpp builds
unchanged against it when OUT_DIR comes first on the include path:

    clang++ -std=c++20 -O3 -DNDEBUG -I OUT_DIR -I micm/include -I micm/test/benchmark \\
        reference/micm_ts1_ref.cpp -o micm_ref

The third body is a parameterized species, as in TS1; user-defined rates take
their labels; InitState sets the canonical initial state in every cell.
Needs only the Python standard library.
"""

import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("canonical", type=Path)
parser.add_argument("out", type=Path)
args = parser.parse_args()
mechanism = json.loads(args.canonical.read_text())


def number(value):
    """A C++ double literal that round-trips."""
    text = repr(float(value))
    return text if any(c in text for c in ".en") else text + ".0"


def species_list(names):
    return "{ " + ", ".join(f'SpeciesByName(gas_phase, "{name}")' for name in names) + " }"


def products_list(products):
    return "{ " + ", ".join(f'micm::StoichSpecies{{ SpeciesByName(gas_phase, "{name}"), {number(y)} }}'
                            for name, y in products) + " }"


def rate(reaction):
    kind = reaction["type"]
    if kind == "ARRHENIUS":
        fields = ", ".join(f".{key}_ = {number(reaction[key])}" for key in ("A", "B", "C", "D", "E"))
        return f"micm::ArrheniusRateConstantParameters{{ {fields} }}"
    if kind == "TROE":
        fields = ", ".join(f".{key}_ = {number(reaction[key])}"
                           for key in ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N"))
        return f"micm::TroeRateConstantParameters{{ {fields} }}"
    if kind == "USER_DEFINED":
        return (f'micm::UserDefinedRateConstantParameters{{ .label_ = "{reaction["label"]}", '
                f'.scaling_factor_ = {number(reaction["scaling_factor"])} }}')
    raise SystemExit(f"reaction type {kind} is not handled")


lines = [
    "// GENERATED FILE by TurboMICM admin/micm_header.py; do not edit.",
    f"// {mechanism['name']}: {mechanism['source']}.",
    f"// {len(mechanism['species'])} species and {len(mechanism['reactions'])} reactions.",
    "#pragma once",
    "",
    "#include <micm/CPU.hpp>",
    "",
    "#include <stdexcept>",
    "#include <string>",
    "#include <string_view>",
    "#include <vector>",
    "",
    "namespace bench",
    "{",
    "  namespace ts1",
    "  {",
    "    inline const micm::Species& SpeciesByName(const micm::Phase& phase, std::string_view name)",
    "    {",
    "      for (const auto& phase_species : phase.phase_species_)",
    "        if (phase_species.species_.name_ == name)",
    "          return phase_species.species_;",
    '      throw std::runtime_error("no species named " + std::string(name));',
    "    }",
    "",
    "    inline micm::Phase CreateGasPhase()",
    "    {",
    "      std::vector<micm::PhaseSpecies> phase_species;",
]
for species in mechanism["species"]:
    lines.append("      {")
    lines.append(f'        auto species = micm::Species("{species["name"]}");')
    if species["third_body"]:
        lines.append("        species.SetThirdBody();")
    lines.append("        phase_species.emplace_back(species);")
    lines.append("      }")
lines += [
    '      return micm::Phase{ "gas", phase_species };',
    "    }",
    "",
    "    inline std::vector<micm::Process> CreateProcesses(const micm::Phase& gas_phase)",
    "    {",
    "      std::vector<micm::Process> processes;",
]
for reaction in mechanism["reactions"]:
    lines.append("      processes.push_back(micm::ChemicalReactionBuilder()")
    lines.append(f"                              .SetReactants({species_list(reaction['reactants'])})")
    lines.append(f"                              .SetProducts({products_list(reaction['products'])})")
    lines.append(f"                              .SetRateConstant({rate(reaction)})")
    lines.append("                              .SetPhase(gas_phase)")
    lines.append("                              .Build());")
lines += [
    "      return processes;",
    "    }",
    "  }  // namespace ts1",
    "",
    "  struct Ts1",
    "  {",
    f'    static constexpr std::string_view kName = "{mechanism["name"]}";',
    "",
    "    template<class Builder>",
    "    static auto Build(Builder builder)",
    "    {",
    "      auto gas_phase = ts1::CreateGasPhase();",
    "      return builder.SetSystem(micm::System(gas_phase))",
    "          .SetReactions(ts1::CreateProcesses(gas_phase))",
    "          .SetIgnoreUnusedSpecies(true)",
    "          .Build();",
    "    }",
    "",
    "    template<class State>",
    "    static void InitState(State& state, micm::Index num_cells)",
    "    {",
    "      std::vector<micm::Real> cells(num_cells);",
    "      auto concentration = [&](const char* name, micm::Real value)",
    "      {",
    "        std::fill(cells.begin(), cells.end(), value);",
    "        state[name] = cells;",
    "      };",
    "      auto parameter = [&](const char* label, micm::Real value)",
    "      {",
    "        std::fill(cells.begin(), cells.end(), value);",
    "        state.SetCustomRateParameter(label, cells);",
    "      };",
]
initial = mechanism["initial"]
for species in mechanism["species"]:
    if not species["third_body"]:
        lines.append(f'      concentration("{species["name"]}", {number(initial["concentrations"].get(species["name"], 0.0))});')
for label in mechanism["parameters"]:
    lines.append(f'      parameter("{label}", {number(initial["parameters"][label])});')
lines += [
    "      for (micm::Index c = 0; c < num_cells; ++c)",
    "      {",
    f"        state.conditions_[c].temperature_ = {number(initial['temperature'])};",
    f"        state.conditions_[c].pressure_ = {number(initial['pressure'])};",
    "        state.conditions_[c].CalculateIdealAirDensity();",
    "      }",
    "    }",
    "  };",
    "}  // namespace bench",
    "",
]
args.out.mkdir(parents=True, exist_ok=True)
(args.out / "ts1_mechanism.hpp").write_text("\n".join(lines))
print(f"{args.out / 'ts1_mechanism.hpp'}: {len(lines)} lines")
