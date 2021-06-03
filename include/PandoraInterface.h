/**
 *  @file   LArReco/include/PandoraInterface.h
 *
 *  @brief  Header file for PandoraInterface.
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_INTERFACE_H
#define PANDORA_ND_INTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

using namespace pandora;

namespace pandora {
class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_reco {

/**
 *  @brief  Parameters class
 */
class Parameters {
public:
  /**
   *  @brief Default constructor
   */
  Parameters();

  std::string m_settingsFile;  ///< The path to the pandora settings file
                               ///< (mandatory parameter)
  std::string m_inputFileName; ///< The path to the input file containing events
                               ///< and/or geometry information

  int m_nEventsToProcess; ///< The number of events to process (default all
                          ///< events in file)
  bool m_shouldDisplayEventNumber; ///< Whether event numbers should be
                                   ///< displayed (default false)

  pandora::InputInt m_nEventsToSkip; ///< The number of events to skip
};

/**
 *  @brief  LArVoxel class
 */
class LArVoxel {
public:
  /**
   *  @brief Default constructor
   */
  LArVoxel(int voxelID, double energyInVoxel, CartesianVector voxelPosVect);

  LArVoxel(const LArVoxel &rhs);

  int m_voxelID; ///< The ID of the voxel

  double m_energyInVoxel; ///< The energy in the voxel

  CartesianVector m_voxelPosVect; ///< A vector containing the x, y, z, position
                                  ///< of the voxel
};

/**
 *  @brief  Create the detector geometry based on the C++ root file
 *
 *  @param  parameters the application parameters
 *  @param  pPrimaryPandora the address of the primary pandora instance
 */
void CreateGeometry(const Parameters &parameters,
                    const pandora::Pandora *const pPrimaryPandora);

/**
 *  @brief  Process events using the supplied pandora instances
 *
 *  @param  parameters the application parameters
 *  @param  pPrimaryPandora the address of the primary pandora instance
 */
void ProcessEvents(const Parameters &parameters,
                   const pandora::Pandora *const pPrimaryPandora);

/**
 *  @brief  Make voxels from g4hits
 *
 *  @param  the g4Hits
 */
std::vector<LArVoxel> makeVoxels(const TG4HitSegment &g4Hit);

/**
 *  @brief  Find the crossings between a given box and hit segment
 *
 *  @param  coordinates of the top of box
 *  @param  coordinates of the bottom of box
 *  @param  start coordinates of the hit segment
 *  @param  stop coordinates of the hit segment
 *  @param  entry point in the box
 *  @param  exit point in the box
 */
int Intersections(const Double_t *const boxBottom, const Double_t *const boxTop,
                  const CartesianVector start, const CartesianVector stop,
                  CartesianVector &pt0, CartesianVector &pt1);

/**
 *  @brief  Find the crossings between a given box and hit segment, using a ray
 *
 *  @param  coordinates of the bottom of box
 *  @param  coordinates of the top of box
 *  @param  start coordinates of the hit segment
 *  @param  sign indicator from the inverse direction vector
 *  @param  inverse direction vector
 *  @param  cross point start
 *  @param  cross point stop
 */
int BoxCrossings(const Double_t *const boxBottom, const Double_t *const boxTop,
                 const CartesianVector start, int *const sign,
                 const CartesianVector invdir, double &t0, double &t1);

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 *
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], Parameters &parameters);

/**
 *  @brief  Print the list of configurable options
 *
 *  @return false, to force abort
 */
bool PrintOptions();

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline Parameters::Parameters()
    : m_settingsFile(""), m_inputFileName(""), m_nEventsToProcess(-1),
      m_shouldDisplayEventNumber(false) {}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(int voxelID, double energyInVoxel,
                          CartesianVector voxelPosVect)
    : m_voxelID(voxelID), m_energyInVoxel(energyInVoxel),
      m_voxelPosVect(voxelPosVect) {}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(const LArVoxel &rhs)
    : m_voxelID(rhs.m_voxelID), m_energyInVoxel(rhs.m_energyInVoxel),
      m_voxelPosVect(rhs.m_voxelPosVect) {}

} // namespace lar_nd_reco

#endif // #ifndef PANDORA_ND_INTERFACE_H
