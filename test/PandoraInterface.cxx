/**
 *  @file   LArRecoMP/test/PandoraInterface.cc
 *
 *  @brief  Implementation of the lar reco mp application
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"

#include "TG4Event.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#include "Api/PandoraApi.h"
#include "Helpers/XmlHelper.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "PandoraInterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <getopt.h>
#include <iostream>
#include <random>
#include <string>

using namespace pandora;
using namespace lar_nd_reco;

int main(int argc, char *argv[]) {
  int errorNo(0);
  const Pandora *pPrimaryPandora(nullptr);

  try {
    Parameters parameters;

    if (!ParseCommandLine(argc, argv, parameters))
      return 1;

#ifdef MONITORING
    TApplication *pTApplication = new TApplication("LArReco", &argc, argv);
    pTApplication->SetReturnFromRun(kTRUE);
#endif

    pPrimaryPandora = new pandora::Pandora();

    if (!pPrimaryPandora)
      throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            LArContent::RegisterAlgorithms(*pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            LArContent::RegisterBasicPlugins(*pPrimaryPandora));

    MultiPandoraApi::AddPrimaryPandoraInstance(pPrimaryPandora);

    CreateGeometry(parameters, pPrimaryPandora);

    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetPseudoLayerPlugin(
            *pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetLArTransformationPlugin(
            *pPrimaryPandora,
            new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraApi::ReadSettings(*pPrimaryPandora, parameters.m_settingsFile));

    ProcessEvents(parameters, pPrimaryPandora);
  } catch (const StatusCodeException &statusCodeException) {
    std::cerr << "Pandora StatusCodeException: "
              << statusCodeException.ToString()
              << statusCodeException.GetBackTrace() << std::endl;
    errorNo = 1;
  } catch (...) {
    std::cerr << "Unknown exception: " << std::endl;
    errorNo = 1;
  }

  MultiPandoraApi::DeletePandoraInstances(pPrimaryPandora);
  return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_reco {

void CreateGeometry(const Parameters &parameters,
                    const Pandora *const pPrimaryPandora) {
  // Get the geometry info from the input root file
  TFile fileSource(parameters.m_inputFileName.c_str(), "READ");
  TGeoManager *pEDepSimGeo = (TGeoManager *)fileSource.Get("EDepSimGeometry");

  if (!pEDepSimGeo) {
    std::cout << "Missing the geometry info" << std::endl;
    return;
  }

  // Start by looking at the top level volume and move down to the one we need
  std::string name;
  std::string needednode("volArgonCubeDetector_PV_0");
  TGeoNode *currentnode = pEDepSimGeo->GetCurrentNode();
  bool foundnode(false);

  TGeoVolume *pMasterVol = pEDepSimGeo->GetMasterVolume();
  for (int i = 0; i < pMasterVol->GetNdaughters(); i++) {
    pEDepSimGeo->CdDown(i);
  }

  while (foundnode == false) {
    currentnode = pEDepSimGeo->GetCurrentNode();
    name = currentnode->GetName();

    int i1 = 0;
    for (int i = 0; i < currentnode->GetNdaughters(); i++) {
      pEDepSimGeo->CdDown(i1);
      TGeoNode *node = pEDepSimGeo->GetCurrentNode();
      name = node->GetName();
      if (name == needednode) {
        foundnode = true;
        break;
      } else if (i + 1 != currentnode->GetNdaughters()) {
        pEDepSimGeo->CdUp();
        i1++;
      }
    }

    if (foundnode == true) {
      break;
    }
  }

  // This should now be the ArgonCube volume
  currentnode = pEDepSimGeo->GetCurrentNode();
  name = currentnode->GetName();
  std::cout << "Current Node: " << name << std::endl;
  std::cout << "Current N daughters: "
            << currentnode->GetVolume()->GetNdaughters() << std::endl;
  std::cout << "  " << std::endl;

  // Get the BBox for the ArgonCube
  TGeoVolume *pCurrentVol = currentnode->GetVolume();
  TGeoShape *pCurrentShape = pCurrentVol->GetShape();
  pCurrentShape->InspectShape();
  TGeoBBox *pBox = dynamic_cast<TGeoBBox *>(pCurrentShape);

  // Now can get origin/width data from the BBox
  double dx = pBox->GetDX(); // Note these are the half widths
  double dy = pBox->GetDY();
  double dz = pBox->GetDZ();
  const double *origin = pBox->GetOrigin();

  // Translate the origin coordinates from the 4th level to the first.
  // Doesn't seem to change anything. Needed?
  Double_t level3[3] = {0., 0., 0.};
  currentnode->LocalToMasterVect(origin, level3);
  Double_t level2[3] = {0., 0., 0.};
  currentnode->LocalToMasterVect(level3, level2);
  Double_t level1[3] = {0., 0., 0.};
  currentnode->LocalToMasterVect(level2, level1);

  // Can now create a geometry using the found parameters
  PandoraApi::Geometry::LArTPC::Parameters geoparameters;

  try {
    geoparameters.m_centerX = level1[0];
    // ATTN: offsets taken by visual comparison with edep-disp
    geoparameters.m_centerY = level1[1] - 675;
    geoparameters.m_centerZ = level1[2] + 6660;
    geoparameters.m_widthX = dx * 2;
    geoparameters.m_widthY = dy * 2;
    geoparameters.m_widthZ = dz * 2;
    // ATTN: parameters past here taken from uboone
    geoparameters.m_larTPCVolumeId = 0;
    geoparameters.m_wirePitchU = 0.300000011921;
    geoparameters.m_wirePitchV = 0.300000011921;
    geoparameters.m_wirePitchW = 0.300000011921;
    geoparameters.m_wireAngleU = 1.04719758034;
    geoparameters.m_wireAngleV = -1.04719758034;
    geoparameters.m_wireAngleW = 0.;
    geoparameters.m_sigmaUVW = 1;
    geoparameters.m_isDriftInPositiveX = 0;
  } catch (const pandora::StatusCodeException &) {
    std::cout << "CreatePandoraLArTPCs - invalid tpc parameter provided"
              << std::endl;
  }

  try {
    PANDORA_THROW_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraApi::Geometry::LArTPC::Create(*pPrimaryPandora, geoparameters));
  } catch (const pandora::StatusCodeException &) {
    std::cout << "CreatePandoraLArTPCs - unable to create tpc, insufficient or "
                 "invalid information supplied"
              << std::endl;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEvents(const Parameters &parameters,
                   const Pandora *const pPrimaryPandora) {
  int nEvents(0);

  TFile fileSource(parameters.m_inputFileName.c_str(), "READ");
  TTree *pEDepSimTree = (TTree *)fileSource.Get("EDepSimEvents");

  if (!pEDepSimTree) {
    std::cout << "Missing the event tree" << std::endl;
    return;
  }

  TG4Event *pEDepSimEvent(nullptr);
  pEDepSimTree->SetBranchAddress("Event", &pEDepSimEvent);

  while ((nEvents < parameters.m_nEventsToProcess) ||
         (0 > parameters.m_nEventsToProcess)) {
    if (parameters.m_shouldDisplayEventNumber)
      std::cout << std::endl
                << "   PROCESSING EVENT: " << nEvents << std::endl
                << std::endl;

    pEDepSimTree->GetEntry(nEvents++);

    if (!pEDepSimEvent)
      return;

    int hitCounter(0);

    for (TG4HitSegmentDetectors::iterator detector =
             pEDepSimEvent->SegmentDetectors.begin();
         detector != pEDepSimEvent->SegmentDetectors.end(); ++detector) {
      std::cout << "Show hits for " << detector->first << " ("
                << detector->second.size() << " hits)" << std::endl;
      std::cout << "                                 " << std::endl;

      std::vector<LArVoxel> voxelList;

      for (TG4HitSegment &g4Hit : detector->second) {
        std::vector<LArVoxel> currentVoxelList = makeVoxels(g4Hit);

        for (LArVoxel &voxel : currentVoxelList) {
          voxelList.push_back(voxel);
        }
      }

      std::cout << "Voxels produced: " << voxelList.size() << std::endl;

      // ATTN: Here we might need to add something to check if there are
      // multiple energy deposits from the same particle into one voxel. How can
      // we tell if they are from the same particle though?

      // Loop over the voxels and make them into caloHits
      for (int i = 0; i < voxelList.size(); i++) {
        PandoraApi::CaloHit::Parameters caloHitParameters;
        caloHitParameters.m_positionVector = voxelList[i].m_voxelPosVect;
        caloHitParameters.m_expectedDirection =
            pandora::CartesianVector(0.f, 0.f, 1.f);
        caloHitParameters.m_cellNormalVector =
            pandora::CartesianVector(0.f, 0.f, 1.f);
        caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
        // ATTN: need to make flexible? Currently just using the fixed size
        caloHitParameters.m_cellSize0 = 0.4;
        caloHitParameters.m_cellSize1 = 0.4;
        caloHitParameters.m_cellThickness = 0.4;
        caloHitParameters.m_nCellRadiationLengths = 1.f;
        caloHitParameters.m_nCellInteractionLengths = 1.f;
        caloHitParameters.m_time = 0.f;
        caloHitParameters.m_inputEnergy = voxelList[i].m_energyInVoxel;
        caloHitParameters.m_mipEquivalentEnergy = voxelList[i].m_energyInVoxel;
        caloHitParameters.m_electromagneticEnergy =
            voxelList[i].m_energyInVoxel;
        caloHitParameters.m_hadronicEnergy = voxelList[i].m_energyInVoxel;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_hitType = pandora::TPC_3D;
        caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_pParentAddress =
            (void *)(static_cast<uintptr_t>(++hitCounter));

        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=,
            PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters));
      } // end voxel loop

      // ATTN: the voxelisation only works with ArgonCube
      break;

    } // end segment detector loop

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraApi::ProcessEvent(*pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraApi::Reset(*pPrimaryPandora));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
std::vector<LArVoxel> makeVoxels(TG4HitSegment &g4Hit) {

  std::vector<LArVoxel> currentVoxelList;

  // Set the variables for the AV bounding box and voxel size
  Double_t voxelSize[3] = {0.4, 0.4, 0.4};
  Double_t boxTop[3] = {370, 160, 930};
  Double_t boxBottom[3] = {-370, -160, 400};
  Double_t boxLength[3] = {740, 320, 530};

  const double epsilon = 1.e-3;
  double energy_deposit = 0.;

  // Get start and stop points for what we want to voxelise
  CartesianVector start(g4Hit.GetStart().X(), g4Hit.GetStart().Y(),
                        g4Hit.GetStart().Z());
  CartesianVector stop(g4Hit.GetStop().X(), g4Hit.GetStop().Y(),
                       g4Hit.GetStop().Z());
  // ATTN: need to match up with geometry before visulisation
  start *= 0.1f; // convert unit to cm
  stop *= 0.1f;

  // Eliminate any tracks which are just points
  if ((stop - start).GetMagnitude() == 0) {
    std::cout << "Cannot have 0cm track length." << std::endl;
    std::cout << "                             " << std::endl;
    return currentVoxelList;
  }

  // Vectors for the intersection points of the hit and the test box
  CartesianVector pt0(0, 0, 0);
  CartesianVector pt1(0, 0, 0);

  // Find intersections (check hit contained within AV)
  int crossings = Intersections(boxBottom, boxTop, start, stop, pt0, pt1);

  if (crossings == 0) {
    std::cout << "No crossing point found..." << std::endl;
    // ATTN: In ML they return voxel list here if no crossing points
    // found. Think I'm doing this a slightly different way, so maybe
    // remove
  }

  std::cout << "   Intersects with bounding box at"
            << " (" << pt0.GetX() << "," << pt0.GetY() << "," << pt0.GetZ()
            << ")"
            << " and (" << pt1.GetX() << "," << pt1.GetY() << "," << pt1.GetZ()
            << ")" << std::endl;

  // Get a unit vector in the direction of the hit segment
  CartesianVector dir = pt1 - pt0;
  double length = dir.GetMagnitude();
  CartesianVector dirnorm = dir.GetUnitVector();

  // Need this to check the inverse vector doesn't end up with a divide by
  // 0
  float val1, val2, val3;
  if (dir.GetX() != 0) {
    val1 = 1 / dir.GetX();
  } else {
    val1 = std::numeric_limits<float>::max();
  }

  if (dir.GetY() != 0) {
    val2 = 1 / dir.GetY();
  } else {
    val2 = std::numeric_limits<float>::max();
  }

  if (dir.GetZ() != 0) {
    val3 = 1 / dir.GetZ();
  } else {
    val3 = std::numeric_limits<float>::max();
  }

  CartesianVector invdirnorm(val1, val2, val3);
  int sign[3];
  sign[0] = (invdirnorm.GetX() < 0);
  sign[1] = (invdirnorm.GetY() < 0);
  sign[2] = (invdirnorm.GetZ() < 0);

  double t0(0);
  double t1(0);
  // size_t nx, ny, nz;

  // Shuffle along this hit segment
  while (true) {
    // Get the position vector that we are going to check out
    CartesianVector pt = pt0 + dirnorm * (t1 + epsilon);
    std::cout << "    New point: " << pt << std::endl;

    // Find which voxel this lives in
    //-----voxel Iding------------
    if (pt.GetX() > boxTop[0] || pt.GetX() < boxBottom[0] ||
        pt.GetY() > boxTop[1] || pt.GetY() < boxBottom[1] ||
        pt.GetZ() > boxTop[2] || pt.GetZ() < boxBottom[2]) {
      std::cout << "Invalid voxel! Out of Geometry!" << std::endl;
      std::cout << "                               " << std::endl;
      return currentVoxelList;
    }

    double xnum = boxLength[0] / voxelSize[0];
    double ynum = boxLength[1] / voxelSize[1];
    double znum = boxLength[2] / voxelSize[2];

    // double xlen = boxLength[0]/xnum;

    double xindex = (pt.GetX() - boxBottom[0]) / voxelSize[0];
    double yindex = (pt.GetY() - boxBottom[1]) / voxelSize[1];
    double zindex = (pt.GetZ() - boxBottom[2]) / voxelSize[2];

    if (xindex == xnum)
      xindex -= 1;
    if (yindex == ynum)
      yindex -= 1;
    if (zindex == znum)
      zindex -= 1;

    int voxelID = (zindex * (xnum * ynum) + yindex * xnum + xindex);
    //---------------------------

    // ATTN: This wasn't working for x or y, but was for z. I couldn't
    // work out why, but seemed to work okay if I just used the indexes
    // straight out
    //--------------id_to_xyz_index---------------
    /*
    nz = voxelID / (xnum * ynum);
    voxelID -= nz * (xnum * ynum);
    ny = voxelID / xnum;
    nx = (voxelID - ny * xnum);

    std::cout << "xindex : " << xindex << "  yindex : " << yindex << "    zindex
    : " << zindex << std::endl; std::cout << "nx : " << nx << "  ny : " << ny <<
    "    nz : " << nz << std::endl;
    */
    //-----------------------------------------------------------------------
    // move the box bounds such that they are moved a voxel along
    // z in good....x and y aren't....but can just use index straight here
    // boxBottom[0] = boxBottom[0] + nx * voxelSize[0];
    // boxBottom[1] = boxBottom[1] + ny * voxelSize[1];
    // boxBottom[2] = boxBottom[2] + nz * voxelSize[2];

    // Define an updated test box
    boxBottom[0] = boxBottom[0] + (xindex * voxelSize[0] - 1e-3);
    boxBottom[1] = boxBottom[1] + (yindex * voxelSize[1] - 1e-3);
    boxBottom[2] = boxBottom[2] + (zindex * voxelSize[2] - 1e-3);
    boxTop[0] = boxBottom[0] + voxelSize[0];
    boxTop[1] = boxBottom[1] + voxelSize[1];
    boxTop[2] = boxBottom[2] + voxelSize[2];

    std::cout << "    Inspecting a voxel id " << voxelID << " ... "
              << std::endl;

    double t1before = t1;
    int cross = BoxCrossings(boxBottom, boxTop, pt0, sign, invdirnorm, t0, t1);
    double t1after = t1;

    // ATTN: This is here so that it can not get stuck in an infinite loop
    // where t1 doesn't shuffle along. Probably should think of a better
    // way to do this
    if (t1before == t1after)
      return currentVoxelList;

    // Consider crossings with the test box
    if (cross == 0) {
      // Test box should have been set up to contain a section of this hit
      std::cout << "      No crossing (not expected) ... breaking" << std::endl;
      return currentVoxelList;
    }

    double dx;
    if (cross == 1) {
      std::cout << "      One crossing: " << pt0 + dir * t1 << std::endl;
      dx = std::min(t1, length);
    } else {
      std::cout << "      Two crossing" << pt0 + dir * t0 << " => "
                << pt0 + dir * t1 << std::endl;
      if (t0 > length)
        dx = length;
      else if (t1 > length)
        dx = length - t0;
      else
        dx = t1 - t0;
    }

    // Find the fraction of energy contained in voxel from the fraction of
    // track in voxel
    double energyInVoxel = dx / length * g4Hit.GetEnergyDeposit();

    if (energyInVoxel < 0) {
      std::cout << "Voxel with negative energy deposited!" << std::endl
                << "  ID = " << voxelID << std::endl
                << "  edep computed from:" << std::endl
                << "      dx = " << dx << ", length = " << length
                << ", TG4HitSegment edep = " << g4Hit.GetEnergyDeposit()
                << std::endl;
      // ATTN: ML throw an error here. Guess we should throw one too?
    }

    energy_deposit += energyInVoxel;

    std::cout << "      Registering voxel id " << voxelID << " t1 =" << t1
              << " (total length = " << length << ")" << std::endl;

    // Push back voxels back into a list

    // Multiply by 10 to match with the detector geometry
    CartesianVector voxelPosVect(boxBottom[0] * 10, boxBottom[1] * 10,
                                 boxBottom[2] * 10);
    LArVoxel currentVoxel(voxelID, energyInVoxel, voxelPosVect);
    currentVoxelList.push_back(currentVoxel);

    // Once t1 is longer than the voxel, break out of the loop
    if (t1 > length) {
      std::cout << "      Reached the segment end (t1 = " << t1
                << " fractional length " << t1 / length << ") ... breaking"
                << std::endl;
      std::cout << "                      " << std::endl;
      return currentVoxelList;
    }

    std::cout << "      Updated t1 = " << t1 << " (fractional length "
              << t1 / length << ")" << std::endl;
    std::cout << "                      " << std::endl;

  } // end while true

  std::cout << "current num of voxels: " << currentVoxelList.size()
            << std::endl;
  std::cout << "                      " << std::endl;

  return currentVoxelList;
}

//------------------------------------------------------------------------------------------------------------------------------------------
int Intersections(const Double_t *const boxBottom, const Double_t *const boxTop,
                  CartesianVector start, CartesianVector stop,
                  CartesianVector &pt0, CartesianVector &pt1) {
  CartesianVector displVec = (stop - start);
  CartesianVector dir = displVec.GetUnitVector();
  float val1, val2, val3;
  if (dir.GetX() != 0) {
    val1 = 1 / dir.GetX();
  } else {
    val1 = std::numeric_limits<float>::max();
  }

  if (dir.GetY() != 0) {
    val2 = 1 / dir.GetY();
  } else {
    val2 = std::numeric_limits<float>::max();
  }

  if (dir.GetZ() != 0) {
    val3 = 1 / dir.GetZ();
  } else {
    val3 = std::numeric_limits<float>::max();
  }

  CartesianVector invdir(val1, val2, val3);
  int sign[3];
  sign[0] = (invdir.GetX() < 0);
  sign[1] = (invdir.GetY() < 0);
  sign[2] = (invdir.GetZ() < 0);

  // Consider if start and stop points are contained
  bool startContained =
      (start.GetX() > boxBottom[0] && start.GetX() < boxTop[0]) &&
      (start.GetY() > boxBottom[1] && start.GetY() < boxTop[1]) &&
      (start.GetZ() > boxBottom[2] && start.GetZ() < boxTop[2]);
  bool stopContained =
      (stop.GetX() > boxBottom[0] && stop.GetX() < boxTop[0]) &&
      (stop.GetY() > boxBottom[1] && stop.GetY() < boxTop[1]) &&
      (stop.GetZ() > boxBottom[2] && stop.GetZ() < boxTop[2]);

  // If the start and stop points are contained, we know entry/exit points
  if (startContained) {
    pt0 = start;
  }
  if (stopContained) {
    pt1 = stop;
  }

  if (!startContained || !stopContained) {
    double t0, t1;

    int cross = BoxCrossings(boxBottom, boxTop, start, sign, invdir, t0, t1);

    if (cross > 0) {
      if ((!startContained && t0 < 0) || t0 > displVec.GetMagnitude())
        cross--;
      if (t1 < 0 || t1 > displVec.GetMagnitude())
        cross--;
    }

    if (cross > 0) {
      const float epsilon = 0.0001;
      if (!startContained)
        pt0 = start + (dir * (t0 + epsilon));

      if (!stopContained)
        pt1 = start + (dir * (t1 - epsilon));
    }

    std::cout << "Number of crossings=" << cross << " for bounding box "
              << boxBottom << "-" << boxTop << " and ray between "
              << "(" << start.GetX() << "," << start.GetY() << ","
              << start.GetZ() << ")"
              << " and (" << stop.GetX() << "," << stop.GetY() << ","
              << stop.GetZ() << ")" << std::endl;

    if (cross > 0) {
      std::cout << "Start point contained?: " << startContained << std::endl;
      if (!startContained)
        std::cout << "  entry point: " << pt0 << "; t0=" << t0 << std::endl;
      std::cout << "Stop point contained?: " << stopContained << std::endl;
      if (!stopContained)
        std::cout << "  exit point: " << pt1 << "; t1=" << t1 << std::endl;
    }

    if (cross == 1 && startContained == stopContained) {
      std::cout << "Unexpected number of crossings (" << cross << ")"
                << " for bounding box and ray between "
                << "(" << start.GetX() << "," << start.GetY() << ","
                << start.GetZ() << ")"
                << " and (" << stop.GetX() << "," << stop.GetY() << ","
                << stop.GetZ() << ")" << std::endl;
      std::cout << "Start point contained?: " << startContained
                << ".  Stop point contained?: " << stopContained << std::endl;
    }

    return cross;
  }

  return 2;
}

//------------------------------------------------------------------------------------------------------------------------------------------
int BoxCrossings(const Double_t *const boxBottom, const Double_t *const boxTop,
                 CartesianVector start, int *const sign, CartesianVector invdir,
                 double &t0, double &t1) {

  float tmin(0), tmax(0), tymin(0), tymax(0), tzmin(0), tzmax(0);

  if (sign[0] == 0) {
    tmin = (boxBottom[0] - start.GetX()) * invdir.GetX();
    tmax = (boxTop[0] - start.GetX()) * invdir.GetX();
  } else if (sign[0] == 1) {
    tmin = (boxTop[0] - start.GetX()) * invdir.GetX();
    tmax = (boxBottom[0] - start.GetX()) * invdir.GetX();
  }

  if (sign[1] == 0) {
    tymin = (boxBottom[1] - start.GetY()) * invdir.GetY();
    tymax = (boxTop[1] - start.GetY()) * invdir.GetY();
  } else if (sign[1] == 1) {
    tymin = (boxTop[1] - start.GetY()) * invdir.GetY();
    tymax = (boxBottom[1] - start.GetY()) * invdir.GetY();
  }

  if ((tmin > tymax) || (tymin > tmax))
    return 0;

  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;

  if (sign[2] == 0) {
    tzmin = (boxBottom[2] - start.GetZ()) * invdir.GetZ();
    tzmax = (boxTop[2] - start.GetZ()) * invdir.GetZ();
  } else if (sign[2] == 1) {
    tzmin = (boxTop[2] - start.GetZ()) * invdir.GetZ();
    tzmax = (boxBottom[2] - start.GetZ()) * invdir.GetZ();
  }

  if ((tmin > tzmax) || (tzmin > tmax))
    return 0;

  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;

  t0 = tmin;
  t1 = tmax;

  if (t1 <= 0)
    return 0;
  if (t0 <= 0)
    return 1;
  return 2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters) {
  if (1 == argc)
    return PrintOptions();

  int c(0);

  while ((c = getopt(argc, argv, ":i:e:n:s:Nh")) != -1) {
    switch (c) {
    case 'i':
      parameters.m_settingsFile = optarg;
      break;
    case 'e':
      parameters.m_inputFileName = optarg;
      break;
    case 'n':
      parameters.m_nEventsToProcess = atoi(optarg);
      break;
    case 's':
      parameters.m_nEventsToSkip = atoi(optarg);
      break;
    case 'N':
      parameters.m_shouldDisplayEventNumber = true;
      break;
    case 'h':
    default:
      return PrintOptions();
    }
  }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions() {
  std::cout
      << std::endl
      << "./bin/PandoraInterface " << std::endl
      << "    -i Settings            (required) [algorithm description: xml]"
      << std::endl
      << "    -e EventsFile          (required) [input edep-sim file, "
         "typically containing events and geometry]"
      << std::endl
      << "    -n NEventsToProcess    (optional) [no. of events to process]"
      << std::endl
      << "    -s NEventsToSkip       (optional) [no. of events to skip in "
         "first file]"
      << std::endl
      << "    -p                     (optional) [print status]" << std::endl
      << "    -N                     (optional) [print event numbers]"
      << std::endl
      << std::endl;

  return false;
}

} // namespace lar_nd_reco
