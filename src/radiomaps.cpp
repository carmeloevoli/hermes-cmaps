// Copyright HERMES Team 2021
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "hermes.h"

namespace hermes {

void computeDMmap(int nside, std::string filename) {
    // ionized gas
    auto gas = std::make_shared<chargedgas::YMW16>();

    // integrator
    auto integrator = std::make_shared<DispersionMeasureIntegrator>(gas);
    auto sun_pos = Vector3QLength(8.3_kpc, 0_kpc, 6_pc);
    integrator->setSunPosition(sun_pos);

    // skymap
    auto skymaps = std::make_shared<DispersionMeasureSkymap>(nside);
    skymaps->setIntegrator(integrator);

    auto output = std::make_shared<outputs::HEALPixFormat>(filename);

    skymaps->compute();
    skymaps->save(output);
}

void computeRMmap(int nside, bool doTurb, std::string filename) {
    // ionized gas
    auto gas = std::make_shared<chargedgas::NE2001Simple>();

    auto mfield = std::make_shared<magneticfields::JF12>();
    if (doTurb) {
        mfield->randomTurbulent(1);
        mfield->randomStriated(1);
    }

    auto integrator = std::make_shared<RotationMeasureIntegrator>(mfield, gas);
    auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
    integrator->setSunPosition(sun_pos);

    // skymap
    auto skymaps = std::make_shared<RotationMeasureSkymap>(nside);
    skymaps->setIntegrator(integrator);

    auto output = std::make_shared<outputs::HEALPixFormat>(filename);

    skymaps->compute();
    skymaps->save(output);
}

void computeSynchroMap(int nside, double freq, std::string filename) {
    // magnetic field
    auto mfield = std::make_shared<magneticfields::JF12>();
    mfield->randomTurbulent(1);
    mfield->randomStriated(1);

    // cosmic ray density models
    std::vector<PID> particletypes = {Electron, Positron};
    auto datapath = "./data/DRAGON_Fornieri2020.fits.gz";
    auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(datapath, particletypes);

    auto integrator = std::make_shared<SynchroIntegrator>(mfield, dragonModel);
    auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_pc);
    integrator->setSunPosition(sun_pos);

    // skymap
    auto skymaps = std::make_shared<RadioSkymap>(nside, freq * 1_MHz);
    skymaps->setIntegrator(integrator);

    auto output = std::make_shared<outputs::HEALPixFormat>(filename);

    skymaps->compute();
    skymaps->save(output);
}

}  // namespace hermes

int main() {
    hermes::computeDMmap(128, "!map-DM-128.fits.gz");
    hermes::computeRMmap(128, true, "!map-RM-128.fits.gz");
    hermes::computeRMmap(128, false, "!map-RM-noturb-128.fits.gz");
    hermes::computeSynchroMap(128, 408, "!map-Synchro-408MHz-128.fits.gz");
    return 0;
}
