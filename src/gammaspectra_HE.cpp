// void computePionDecayHESpectrum(int nside, neutralgas::GasType gasType, std::string filename) {
//     // cosmic ray density models
//     std::vector<PID> particletypes = {Proton, Helium};
//     auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

//     // interaction
//     auto crosssection = std::make_shared<interactions::KelnerAharonianGamma>();

//     // target gas
//     auto gas = std::make_shared<neutralgas::RingModel>(gasType);

//     // integrator
//     auto intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, crosssection);
//     auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
//     intPion->setSunPosition(sun_pos);
//     intPion->setupCacheTable(100, 100, 50);

//     // skymap
//     auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 8);
//     skymaps->setIntegrator(intPion);

//     auto output = std::make_shared<outputs::HEALPixFormat>(filename);

//     skymaps->compute();
//     skymaps->save(output);
// }

// void computePionDecayNeutrinoSpectrum(int nside, neutralgas::GasType gasType,
//                                       std::string filename) {
//     // cosmic ray density models
//     std::vector<PID> particletypes = {Proton, Helium};
//     auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

//     // interaction
//     auto crosssection = std::make_shared<interactions::KelnerAharonianNeutrino>();

//     // target gas
//     auto gas = std::make_shared<neutralgas::RingModel>(gasType);

//     // integrator
//     auto intPion = std::make_shared<PiZeroIntegrator>(dragonModel, gas, crosssection);
//     auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
//     intPion->setSunPosition(sun_pos);
//     intPion->setupCacheTable(100, 100, 50);

//     // skymap
//     auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 8);
//     skymaps->setIntegrator(intPion);

//     auto output = std::make_shared<outputs::HEALPixFormat>(filename);

//     skymaps->compute();
//     skymaps->save(output);
// }

// void computeDarkMatterSpectrum(int nside, std::string filename) {
//     // photon field
//     auto gamma_slope = 1.0;
//     auto concentration = 20;
//     auto M_200 = 1.4 * 0.7 * 8e11 * units::sun_mass;  // check

//     auto dmProfile = std::make_shared<darkmatter::NFWGProfile>(gamma_slope, concentration,
//     M_200);

//     //    auto units = 1.6726219e-27_kg / 0.938 / (1_cm * 1_cm * 1_cm);
//     //    std::cout << dmProfile->getMassDensity(8.5_kpc) / units << "\n";

//     darkmatter::Channel dmChannel = darkmatter::Channel::W;
//     darkmatter::Mass dmMass = darkmatter::Mass::m30TeV;

//     auto dmSpectrum = std::make_shared<darkmatter::PPPC4DMIDSpectrum>(dmChannel, dmMass);

//     // integrator
//     auto intDM = std::make_shared<DarkMatterIntegrator>(dmSpectrum, dmProfile);
//     auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
//     intDM->setSunPosition(sun_pos);

//     // skymap
//     auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 30_TeV, 50);

//     skymaps->setIntegrator(intDM);

//     auto output = std::make_shared<outputs::HEALPixFormat>(filename);

//     skymaps->compute();
//     skymaps->save(output);
// }

// // void computePionDecayAbsorptionSpectrum(int nside, neutralgas::GasType gasType,
// //                                         std::string filename) {
// //     // cosmic ray density models
// //     std::vector<PID> particletypes = {Proton, Helium};
// //     auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

// //     // photon background
// //     auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());

// //     // interaction
// //     auto crosssection = std::make_shared<interactions::KelnerAharonianGamma>();

// //     // target gas
// //     auto gas = std::make_shared<neutralgas::RingModel>(gasType);

// //     // mask
// //     const std::array<QAngle, 2> longitude{0_deg, 360_deg};
// //     const std::array<QAngle, 2> latitude{-8_deg, 8_deg};
// //     auto mask = std::make_shared<RectangularWindow>(latitude, longitude);

// //     // integrator
// //     auto intPion =
// //         std::make_shared<PiZeroAbsorptionIntegrator>(dragonModel, gas, isrf, crosssection);
// //     auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
// //     intPion->setSunPosition(sun_pos);
// //     intPion->setupCacheTable(100, 100, 50);

// //     // skymap
// //     auto skymaps = std::make_shared<GammaSkymapRange>(nside, 1_TeV, 1_PeV, 3 * 4);
// //     skymaps->setMask(mask);
// //     skymaps->setIntegrator(intPion);

// //     auto output = std::make_shared<outputs::HEALPixFormat>(filename);

// //     skymaps->compute();
// //     skymaps->save(output);
// // }

//     // HE PLOT
//     // hermes::computePionDecayHESpectrum(64, GasType::HI,
//     //                                    "!spectrum-PionDecay-HI-1TeV-1PeV-64.fits.gz");
//     // hermes::computePionDecayHESpectrum(64, GasType::H2,
//     //                                    "!spectrum-PionDecay-H2-1TeV-1PeV-64.fits.gz");
//     // hermes::computePionDecayNeutrinoSpectrum(64, GasType::HI,
//     //                                          "!spectrum-PionDecayNu-HI-1TeV-1PeV-64.fits.gz");
//     // hermes::computePionDecayNeutrinoSpectrum(64, GasType::H2,
//     //                                          "!spectrum-PionDecayNu-H2-1TeV-1PeV-64.fits.gz");
//     // hermes::computeDarkMatterSpectrum(64, "!spectrum-DarkMatter-1TeV-30TeV-64.fits.gz");

//     // hermes::computePionDecayAbsorptionSpectrum(16, GasType::HI,
//     // "!spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.fits.gz");

// void computePionDecayAbsorptionMap(int nside, double E_gamma, neutralgas::GasType gasType, std::string filename) {
//   // cosmic ray density models
//   std::vector<PID> particletypes = {Proton, Helium};
//   auto dragonModel = std::make_shared<cosmicrays::Dragon2D>(dragonFile, particletypes);

//   // photon background
//   auto isrf = std::make_shared<photonfields::ISRF>(photonfields::ISRF());

//   // interaction
//   auto kamae_crosssection = std::make_shared<interactions::Kamae06Gamma>();

//   // target gas
//   auto gas = std::make_shared<neutralgas::RingModel>(gasType);

//   // integrator
//   auto intPion = std::make_shared<PiZeroAbsorptionIntegrator>(dragonModel, gas, isrf, kamae_crosssection);
//   auto sun_pos = Vector3QLength(8.5_kpc, 0_kpc, 0_kpc);
//   intPion->setSunPosition(sun_pos);
//   intPion->setupCacheTable(200, 200, 100);

//   // skymap
//   auto skymaps = std::make_shared<GammaSkymap>(nside, E_gamma * 1_GeV);
//   skymaps->setIntegrator(intPion);

//   auto output = std::make_shared<outputs::HEALPixFormat>(filename);

//   skymaps->compute();
//   skymaps->save(output);
// }