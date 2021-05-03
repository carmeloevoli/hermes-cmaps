// Copyright HERMES Team 2021
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

const std::string dragonFile = "./data/DRAGON_Fornieri2020.fits.gz";

const auto minFrequency = 1_MHz;
const auto maxFrequency = 100_GHz;
const int sizeFrequency = 20;

void computeSynchroSpectrum(int nside, std::string filename) {
  // magnetic field
  auto mfield = std::make_shared<magneticfields::JF12>();
  mfield->randomTurbulent(1);
  mfield->randomStriated(1);
  // mfield->setUseTurbulent(false);

  // cosmic ray density models
  std::vector<PID> particletypes = {Electron, Positron};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  auto integrator = std::make_shared<SynchroIntegrator>(mfield, dragonModel);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
  integrator->setSunPosition(sun_pos);

  // skymap
  auto skymaps = std::make_shared<RadioSkymapRange>(RadioSkymapRange(nside, minFrequency, maxFrequency, sizeFrequency));
  skymaps->setIntegrator(integrator);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeSynchroSpectrumAbsorption(int nside, std::string filename) {
  // magnetic field
  auto mfield = std::make_shared<magneticfields::JF12>();
  mfield->randomTurbulent(1);
  mfield->randomStriated(1);
  // mfield->setUseTurbulent(false);

  // ionized gas
  auto gas = std::make_shared<ionizedgas::YMW16>();

  // cosmic ray density models
  std::vector<PID> particletypes = {Electron, Positron};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  auto integrator = std::make_shared<SynchroAbsorptionIntegrator>(mfield, dragonModel, gas);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
  integrator->setSunPosition(sun_pos);

  // skymap
  auto skymaps = std::make_shared<RadioSkymapRange>(RadioSkymapRange(nside, minFrequency, maxFrequency, sizeFrequency));
  skymaps->setIntegrator(integrator);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeSynchroSpectrumNoTurb(int nside, std::string filename) {
  // magnetic field
  auto mfield = std::make_shared<magneticfields::JF12>();
  mfield->setUseTurbulent(false);
  mfield->setUseStriated(false);

  // cosmic ray density models
  std::vector<PID> particletypes = {Electron, Positron};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  auto integrator = std::make_shared<SynchroIntegrator>(mfield, dragonModel);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
  integrator->setSunPosition(sun_pos);

  // skymap
  auto skymaps = std::make_shared<RadioSkymapRange>(RadioSkymapRange(nside, minFrequency, maxFrequency, sizeFrequency));
  skymaps->setIntegrator(integrator);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeFreeFreeSpectrum(int nside, std::string filename) {
  // ionized gas
  auto gas = std::make_shared<ionizedgas::YMW16>();

  // integrator
  auto integrator = std::make_shared<FreeFreeIntegrator>(gas);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
  integrator->setSunPosition(sun_pos);

  // skymap
  auto skymaps = std::make_shared<RadioSkymapRange>(RadioSkymapRange(nside, minFrequency, maxFrequency, sizeFrequency));
  skymaps->setIntegrator(integrator);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

}  // namespace hermes

int main() {
  hermes::computeSynchroSpectrumNoTurb(64, "!spectrum-Synchro-noturb-64.fits.gz");
  hermes::computeSynchroSpectrum(64, "!spectrum-Synchro-64.fits.gz");
  hermes::computeFreeFreeSpectrum(64, "!spectrum-FreeFree-64.fits.gz");
  hermes::computeSynchroSpectrumAbsorption(64, "!spectrum-Synchro-absorption-64.fits.gz");

  return 0;
}
