// Copyright HERMES Team 2021
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

const std::string dragonFile = "./data/DRAGON_Fornieri2020.fits.gz";

void computePionDecayHESpectrum(int nside, neutralgas::GasType gasType, std::string filename, bool doAbsorption) {
  // cosmic ray density models
  std::vector<PID> particletypes = {Proton, Helium};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  // interaction
  //auto crosssection = std::make_shared<interactions::KelnerAharonianGamma>();
  auto crosssection = std::make_shared<interactions::AAfragGamma>(interactions::AAfragParticle::GAMMA);

  std::cout << "cross section read\n";
  exit(1);

  // target gas
  auto gas = std::make_shared<neutralgas::RingModel>(gasType);

  // mask
  const std::array<QAngle, 2> longitude{0_deg, 360_deg};
  const std::array<QAngle, 2> latitude{-8_deg, 8_deg};
  const auto mask = std::make_shared<RectangularWindow>(latitude, longitude);

  // integrator
  std::shared_ptr<PiZeroIntegrator> intPion;
  if (doAbsorption)
      intPion = std::make_shared<PiZeroAbsorptionIntegrator>(dragonModel, gas, crosssection);
  else
      intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, crosssection);
    
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intPion->setObsPosition(sun_pos);
  intPion->setupCacheTable(100, 100, 50);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 8);
  skymaps->setMask(mask);
  skymaps->setIntegrator(intPion);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

//void computePionDecayAbsorptionSpectrum(int nside, neutralgas::GasType gasType, std::string filename) {
//  // cosmic ray density models
//  std::vector<PID> particletypes = {Proton, Helium};
//  const auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);
//
//  // photon background
//  const auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());
//
//  // interaction
//  const auto crosssection = std::make_shared<interactions::KelnerAharonianGamma>();
//
//  // target gas
//  const auto gas = std::make_shared<neutralgas::RingModel>(gasType);
//
//  // mask
//  const std::array<QAngle, 2> longitude{0_deg, 360_deg};
//  const std::array<QAngle, 2> latitude{-8_deg, 8_deg};
//  const auto mask = std::make_shared<RectangularWindow>(latitude, longitude);
//
//  // integrator
//  auto intPion = std::make_shared<PiZeroAbsorptionIntegrator>(dragonModel, gas, isrf, crosssection);
//  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
//  intPion->setObsPosition(sun_pos);
//  intPion->setupCacheTable(100, 100, 50);
//
//  // skymap
//  auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 8);
//  skymaps->setMask(mask);
//  skymaps->setIntegrator(intPion);
//
//  auto output = std::make_shared<outputs::HEALPixFormat>(filename);
//
//  skymaps->compute();
//  skymaps->save(output);
//}

void computePionDecayNeutrinoSpectrum(int nside, neutralgas::GasType gasType, std::string filename) {
  // cosmic ray density models
  std::vector<PID> particletypes = {Proton, Helium};
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

  // interaction
  auto crosssection = std::make_shared<interactions::KelnerAharonianNeutrino>();

  // target gas
  auto gas = std::make_shared<neutralgas::RingModel>(gasType);

  // integrator
  auto intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, crosssection);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intPion->setObsPosition(sun_pos);
  intPion->setupCacheTable(100, 100, 50);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 8);
  skymaps->setIntegrator(intPion);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeDarkMatterSpectrum(int nside, std::string filename) {
  // photon field
  auto gamma_slope = 1.0;
  auto concentration = 20;
  auto M_200 = 1.4 * 0.7 * 8e11 * units::sun_mass;  // check

  auto dmProfile = std::make_shared<darkmatter::NFWGProfile>(gamma_slope, concentration, M_200);

  darkmatter::Channel dmChannel = darkmatter::Channel::W;
  darkmatter::Mass dmMass = darkmatter::Mass::m30TeV;

  auto dmSpectrum = std::make_shared<darkmatter::PPPC4DMIDSpectrum>(dmChannel, dmMass);

  // integrator
  auto intDM = std::make_shared<DarkMatterIntegrator>(dmSpectrum, dmProfile);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intDM->setObsPosition(sun_pos);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 30_TeV, 50);

  skymaps->setIntegrator(intDM);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

void computeDarkMatterSecondarySpectrum(int nside, std::string filename) {
  // photon field
  auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());

  // cosmic ray density models
  std::vector<PID> particletypes = {Electron, Positron};

  auto datapath = "./data/run_2D_DM_30000.fits.gz";
  auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(datapath, particletypes);

  // interaction
  auto kleinnishina = std::make_shared<interactions::KleinNishina>();

  // integrator
  auto intIC = std::make_shared<InverseComptonIntegrator>(dragonModel, isrf, kleinnishina);
  auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
  intIC->setObsPosition(sun_pos);
  intIC->setupCacheTable(40, 40, 20);

  // skymap
  auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 30_TeV, 30);

  skymaps->setIntegrator(intIC);

  auto output = std::make_shared<outputs::HEALPixFormat>(filename);

  skymaps->compute();
  skymaps->save(output);
}

}  // namespace hermes

using GasType = hermes::neutralgas::GasType;

int main() {
    hermes::computePionDecayHESpectrum(64, GasType::HI, "!spectrum-PionDecay-HI-1TeV-1PeV-64.fits.gz", false);
  //hermes::computePionDecayHESpectrum(64, GasType::HI, "!spectrum-PionDecayWa-HI-1TeV-1PeV-64.fits.gz", true);

  // hermes::computePionDecayHESpectrum(64, GasType::HI, "!spectrum-PionDecay-HI-1TeV-1PeV-64.fits.gz");
  // hermes::computePionDecayHESpectrum(64, GasType::H2, "!spectrum-PionDecay-H2-1TeV-1PeV-64.fits.gz");
  // hermes::computePionDecayNeutrinoSpectrum(64, GasType::HI, "!spectrum-PionDecayNu-HI-1TeV-1PeV-64.fits.gz");
  // hermes::computePionDecayNeutrinoSpectrum(64, GasType::H2, "!spectrum-PionDecayNu-H2-1TeV-1PeV-64.fits.gz");
  // hermes::computeDarkMatterSpectrum(256, "!spectrum-DarkMatter-1TeV-30TeV-256.fits.gz");
  // hermes::computeDarkMatterSecondarySpectrum(64, "!spectrum-DarkMatterSecondary-1TeV-30TeV-64.fits.gz");
  // hermes::computePionDecayAbsorptionSpectrum(16, GasType::HI,
  // "!spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.fits.gz",
  //                                          true);
  // hermes::computePionDecayAbsorptionSpectrum(16, GasType::HI, "!spectrum-PionDecay-HI-1TeV-1PeV-16.fits.gz", false);

//  hermes::computePionDecayAbsorptionSpectrum(16, GasType::H2, "!spectrum-PionDecayAbsorption-H2-1TeV-1PeV-16.fits.gz",
//                                             true);
//  hermes::computePionDecayAbsorptionSpectrum(16, GasType::H2, "!spectrum-PionDecay-H2-1TeV-1PeV-16.fits.gz", false);
}
