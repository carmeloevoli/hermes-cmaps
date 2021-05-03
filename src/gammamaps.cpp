// Copyright HERMES Team 2021
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

const std::string dragonFile = "./data/DRAGON_Fornieri2020.fits.gz";

void computePionDecayMap(int nside, double E_gamma, neutralgas::GasType gasType, std::string filename) {
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
  intPion->setupCacheTable(200, 200, 100);

  // skymap
  auto skymaps = std::make_shared<GammaSkymap>(nside, E_gamma * 1_GeV);
  skymaps->setIntegrator(intPion);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computePionDecayRingMap(int iring, int nside, double E_gamma, neutralgas::GasType gasType, std::string filename) {
  // cosmic ray density models
  std::vector<PID> particletypes = {Proton, Helium};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  // interaction
  auto kamae_crosssection = std::make_shared<interactions::Kamae06Gamma>();

  // target gas
  auto gas = std::make_shared<neutralgas::RingModel>(gasType);
  for (size_t i = 0; i < 12; ++i)
     if (i != iring) gas->disableRingNo(i);

  // integrator
  auto intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, kamae_crosssection);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intPion->setSunPosition(sun_pos);
  intPion->setupCacheTable(200, 200, 100);

  // skymap
  auto skymaps = std::make_shared<GammaSkymap>(nside, E_gamma * 1_GeV);
  skymaps->setIntegrator(intPion);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeInverseComptonMap(int nside, double E_gamma, std::string filename) {
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
  intIC->setupCacheTable(100, 100, 50);

  // skymap
  auto skymaps = std::make_shared<GammaSkymap>(nside, E_gamma * 1_GeV);
  skymaps->setIntegrator(intIC);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

}  // namespace hermes

using GasType = hermes::neutralgas::GasType;

int main(void) {
  //hermes::computePionDecayMap(256, 10, GasType::HI, "!map-Pi0-HI-10GeV-256.fits.gz");
  //hermes::computePionDecayMap(256, 10, GasType::H2, "!map-Pi0-H2-10GeV-256.fits.gz");
  //hermes::computeInverseComptonMap(128, 10, "!map-IC-10GeV-128.fits.gz");
  
  for (int i = 0; i < 11; i++) {
     hermes::computePionDecayRingMap(i, 256, 5, GasType::H2, "!map-ring-Pi0-H2-5GeV-256-" + std::to_string(i) + ".fits.gz");
     hermes::computePionDecayRingMap(i, 256, 5, GasType::HI, "!map-ring-Pi0-HI-5GeV-256-" + std::to_string(i) + ".fits.gz");
  }
    
  return 0;
}
