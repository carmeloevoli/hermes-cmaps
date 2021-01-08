// Copyright HERMES Team 2021
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

const std::string dragonFile = "./data/DRAGON_Fornieri2020.fits.gz";

const auto minEnergy = 100_MeV;
const auto maxEnergy = 1_TeV;
const int sizeEnergy = 32;

void computePionDecaySpectrum(int nside, neutralgas::GasType gasType, std::string filename) {
  // cosmic ray density models
  std::vector<PID> particletypes = {Proton, Helium};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  // interaction
  auto kamae_crosssection = std::make_shared<interactions::Kamae06Gamma>();

  // target gas
  auto gas = std::make_shared<neutralgas::RingModel>(gasType);

  // integrator
  auto intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, kamae_crosssection);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intPion->setSunPosition(sun_pos);
  intPion->setupCacheTable(100, 100, 50);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, minEnergy, maxEnergy, sizeEnergy);
  skymaps->setIntegrator(intPion);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeInverseComptonSpectrum(int nside, std::string filename) {
  // photon field
  auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());

  // cosmic ray density models
  std::vector<PID> particletypes = {Electron, Positron};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  // interaction
  auto kleinnishina = std::make_shared<interactions::KleinNishina>();

  // integrator
  auto intIC = std::make_shared<InverseComptonIntegrator>(dragonModel, isrf, kleinnishina);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intIC->setSunPosition(sun_pos);
  intIC->setupCacheTable(50, 50, 30);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, minEnergy, maxEnergy, sizeEnergy);

  skymaps->setIntegrator(intIC);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeBremsstrahlungSpectrum(int nside, neutralgas::GasType gasType, std::string filename) {
  // cosmic ray density models
  std::vector<PID> leptons = {Electron, Positron};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, leptons);

  // interaction
  auto crosssection = std::make_shared<interactions::BremsstrahlungTsai74>();

  // target gas
  auto gas = std::make_shared<neutralgas::RingModel>(gasType);

  // integrator
  auto intBremss = std::make_shared<BremsstrahlungIntegrator>(dragonModel, gas, crosssection);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intBremss->setSunPosition(sun_pos);
  intBremss->setupCacheTable(100, 100, 50);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, minEnergy, maxEnergy, sizeEnergy);
  skymaps->setIntegrator(intBremss);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

}  // namespace hermes

using GasType = hermes::neutralgas::GasType;

int main(void) {
  hermes::computePionDecaySpectrum(64, GasType::HI, "!spectrum-PionDecay-HI-100MeV-1TeV-64.fits.gz");
  hermes::computePionDecaySpectrum(64, GasType::H2, "!spectrum-PionDecay-H2-100MeV-1TeV-64.fits.gz");
  hermes::computeInverseComptonSpectrum(64, "!spectrum-InverseCompton-100MeV-1TeV-64.fits.gz");
  hermes::computeBremsstrahlungSpectrum(64, GasType::HI, "!spectrum-Bremsstrahlung-HI-100MeV-1TeV-64.fits.gz");
  hermes::computeBremsstrahlungSpectrum(64, GasType::H2, "!spectrum-Bremsstrahlung-H2-100MeV-1TeV-64.fits.gz");
  return 0;
}
