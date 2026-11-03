/**
 * @file configureAutoPas.h
 * @date 03.11.2026
 */
#pragma once

#include <array>
#include <string>

#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"

/**
 * Forwards all AutoPas related options of a configuration to an AutoPas object which is not yet initialized.
 * This is shared between md-flexible's Simulation and the Python interface.
 * @param autoPas The object to configure. AutoPas::init() is left to the caller.
 * @param configuration The source of all options.
 * @param boxMin Lower corner of the domain of this AutoPas object.
 * @param boxMax Upper corner of the domain of this AutoPas object.
 * @param outputSuffix Identifier for all files which AutoPas creates.
 */
void configureAutoPas(autopas::AutoPas<ParticleType> &autoPas, const MDFlexConfig &configuration,
                      const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                      const std::string &outputSuffix);
