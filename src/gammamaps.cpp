// Copyright
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

void computePionDecayMap(int nside, double E_gamma, neutralgas::GasType gasType, std::string filename) {
    // cosmic ray density models
    std::vector<PID> particletypes = {Proton, Helium};
    auto datapath = getDataPath("CosmicRays/Gaggero17/run_2D.fits.gz");
    auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(datapath, particletypes);
    
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

void computePionDecayAbsorptionMap(int nside, double E_gamma, neutralgas::GasType gasType, std::string filename) {
    // cosmic ray density models
    std::vector<PID> particletypes = {Proton, Helium};
    auto datapath = getDataPath("CosmicRays/Gaggero17/run_2D.fits.gz");
    auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(datapath, particletypes);
    
    // photon background
    auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());

    // interaction
    auto kamae_crosssection = std::make_shared<interactions::Kamae06Gamma>();
    
    // target gas
    auto gas = std::make_shared<neutralgas::RingModel>(gasType);
    
    // integrator
    auto intPion = std::make_shared<PiZeroAbsorptionIntegrator>(dragonModel, gas, isrf, kamae_crosssection);
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

void computeInverseComptonMap(int nside, std::string filename) {
    // photon field
    auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());
    
    // cosmic ray density models
    std::vector<PID> particletypes = {Electron, Positron};
    
    auto datapath = getDataPath("CosmicRays/Gaggero17/run_2D.fits.gz");
    
    auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(datapath, particletypes);
    
    // interaction
    auto kleinnishina = std::make_shared<interactions::KleinNishina>();
    
    // integrator
    auto intIC = std::make_shared<InverseComptonIntegrator>(dragonModel, isrf, kleinnishina);
    auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
    intIC->setSunPosition(sun_pos);
    intIC->setupCacheTable(100, 100, 50);
    
    // skymap
    auto skymaps = std::make_shared<GammaSkymap>(nside, 10_GeV);
    skymaps->setIntegrator(intIC);
    
    auto output = std::make_shared<outputs::HEALPixFormat>(filename);
    
    skymaps->compute();
    skymaps->save(output);
}

}  // namespace hermes

using GasType = hermes::neutralgas::GasType;

int main(void) {
//    hermes::computePionDecayMap(256, GasType::HI, "!map-PionDecay-HI-10GeV-256.fits.gz");
//    hermes::computePionDecayMap(256, GasType::H2, "!map-PionDecay-H2-10GeV-256.fits.gz");
//    hermes::computeInverseComptonMap(256, "!map-IC-10GeV-256.fits.gz");
    hermes::computePionDecayAbsorptionMap(8, 1e3, GasType::HI, "!map-PionDecayAbsorption-HI-1TeV-8.fits.gz");
    hermes::computePionDecayMap(8, 1e3, GasType::HI, "!map-PionDecay-HI-1TeV-8.fits.gz");

    return 0;
}
