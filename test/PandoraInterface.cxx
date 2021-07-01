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
#include "Geometry/LArTPC.h"
#include "Helpers/XmlHelper.h"
#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
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

  // mm to cm conversion
  const double mm2cm(0.1);

  // Now can get origin/width data from the BBox
  double dx = pBox->GetDX() * mm2cm; // Note these are the half widths
  double dy = pBox->GetDY() * mm2cm;
  double dz = pBox->GetDZ() * mm2cm;
  const double *origin = pBox->GetOrigin();

  // Translate the origin coordinates from the 4th level to the first.
  // Doesn't seem to change anything. Needed?
  double level3[3] = {0.0, 0.0, 0.0};
  currentnode->LocalToMasterVect(origin, level3);
  double level2[3] = {0.0, 0.0, 0.0};
  currentnode->LocalToMasterVect(level3, level2);
  double level1[3] = {0.0, 0.0, 0.0};
  currentnode->LocalToMasterVect(level2, level1);

  // Can now create a geometry using the found parameters
  PandoraApi::Geometry::LArTPC::Parameters geoparameters;

  try {
    geoparameters.m_centerX = level1[0] * mm2cm;
    // ATTN: offsets taken by visual comparison with edep-disp
    geoparameters.m_centerY = (level1[1] - 675.0) * mm2cm;
    geoparameters.m_centerZ = (level1[2] + 6660.0) * mm2cm;
    geoparameters.m_widthX = dx * 2.0;
    geoparameters.m_widthY = dy * 2.0;
    geoparameters.m_widthZ = dz * 2.0;

    // ATTN: parameters past here taken from uboone
    geoparameters.m_larTPCVolumeId = 0;
    geoparameters.m_wirePitchU = 0.300000011921;
    geoparameters.m_wirePitchV = 0.300000011921;
    geoparameters.m_wirePitchW = 0.300000011921;
    geoparameters.m_wireAngleU = 1.04719758034;
    geoparameters.m_wireAngleV = -1.04719758034;
    geoparameters.m_wireAngleW = 0.0;
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

  // Factory for creating LArCaloHits
  lar_content::LArCaloHitFactory m_larCaloHitFactory;

  // Detector volume for voxelising the hits
  const GeometryManager *geom = pPrimaryPandora->GetGeometry();
  const LArTPC &tpc = geom->GetLArTPC();

  double botX = tpc.GetCenterX() - 0.5 * tpc.GetWidthX();
  double botY = tpc.GetCenterY() - 0.5 * tpc.GetWidthY();
  double botZ = tpc.GetCenterZ() - 0.5 * tpc.GetWidthZ();
  double topX = botX + tpc.GetWidthX();
  double topY = botY + tpc.GetWidthY();
  double topZ = botZ + tpc.GetWidthZ();
  const double voxelWidth(parameters.m_voxelWidth);
  LArGrid grid(pandora::CartesianVector(botX, botY, botZ),
               pandora::CartesianVector(topX, topY, topZ),
               pandora::CartesianVector(voxelWidth, voxelWidth, voxelWidth));

  std::cout << "Total grid volume: bot = " << grid.m_bottom
            << "\n top = " << grid.m_top << std::endl;
  std::cout << "Making voxels with size " << grid.m_binWidths << std::endl;

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
        std::vector<LArVoxel> currentVoxelList = makeVoxels(g4Hit, grid);

        std::vector<int> idList;
        std::vector<int> indexList;

        // Check for double matches in voxel ID - if found, combine energy
        for (LArVoxel &voxel : currentVoxelList) {
          idList.push_back(voxel.m_voxelID);
        }

        int prev_id = -1;
        int prev_num_counts = -1;

        for (int &id : idList) {
          // Find number of matches to a given VoxelID
          int num_counts = std::count(idList.begin(), idList.end(), id);

          if (num_counts > 1) {
            // If a match is found, push the index back into a list
            for (int i = 0; i < num_counts; i++) {
              if (i == 0) {
                auto it = find(idList.begin(), idList.end(), id);
                int index = it - idList.begin();
                indexList.push_back(index);
              } else {
                std::vector<int>::iterator iter =
                    idList.begin() + indexList.back();
                auto it = find(iter, idList.end(), id);
                int index = it - iter + indexList.back() + 1;
                indexList.push_back(index);
              }
            }
          }
          // Erasing of extra index counts
          if (num_counts > 1 && num_counts == prev_num_counts &&
              id == prev_id) {
            for (int i = 0; i < num_counts; i++) {
              indexList.pop_back();
            }
          }

          if (num_counts < 2) {
            prev_id = -1;
            prev_num_counts = -1;
          } else {
            prev_id = id;
            prev_num_counts = num_counts;
          }
        }
        while (indexList.size() > 1 && currentVoxelList.size() > 0) {
          int first = indexList[0];
          int second = indexList[1];
          LArVoxel voxel1 = currentVoxelList[first];
          LArVoxel voxel2 = currentVoxelList[second];

          // If voxels have the same ID, add their energies
          if (voxel1.m_voxelID == voxel2.m_voxelID) {
            LArVoxel newVoxel(voxel1.m_voxelID,
                              voxel1.m_energyInVoxel + voxel2.m_energyInVoxel,
                              voxel1.m_voxelPosVect);

            currentVoxelList.insert(currentVoxelList.begin() + first, newVoxel);

            currentVoxelList.erase(currentVoxelList.begin() + second + 1);
            currentVoxelList.erase(currentVoxelList.begin() + first + 1);
            indexList.erase(indexList.begin());
            for (int &index : indexList) {
              index = index - 1;
            }
          } else {
            indexList.erase(indexList.begin());
          }
        }

        for (LArVoxel &voxel : currentVoxelList) {
          voxelList.push_back(voxel);
        }
      }

      std::cout << "Produced " << voxelList.size() << " voxels from "
                << detector->second.size() << " hit segments." << std::endl;

      // Loop over the voxels and make them into caloHits
      for (int i = 0; i < voxelList.size(); i++) {
        const LArVoxel voxel = voxelList[i];
        const CartesianVector voxelPos(voxel.m_voxelPosVect);
        const double voxelE = voxel.m_energyInVoxel;

        lar_content::LArCaloHitParameters caloHitParameters;
        caloHitParameters.m_positionVector = voxelPos;
        caloHitParameters.m_expectedDirection =
            pandora::CartesianVector(0.f, 0.f, 1.f);
        caloHitParameters.m_cellNormalVector =
            pandora::CartesianVector(0.f, 0.f, 1.f);
        caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
        caloHitParameters.m_cellSize0 = voxelWidth;
        caloHitParameters.m_cellSize1 = voxelWidth;
        caloHitParameters.m_cellThickness = voxelWidth;
        caloHitParameters.m_nCellRadiationLengths = 1.f;
        caloHitParameters.m_nCellInteractionLengths = 1.f;
        caloHitParameters.m_time = 0.f;
        caloHitParameters.m_inputEnergy = voxelE;
        caloHitParameters.m_mipEquivalentEnergy = voxelE;
        caloHitParameters.m_electromagneticEnergy = voxelE;
        caloHitParameters.m_hadronicEnergy = voxelE;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_hitType = pandora::TPC_3D;
        caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_pParentAddress =
            (void *)(static_cast<uintptr_t>(++hitCounter));
        caloHitParameters.m_larTPCVolumeId = 0;
        caloHitParameters.m_daughterVolumeId = 0;

        if (parameters.m_use3D) {
          PANDORA_THROW_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters,
                                          m_larCaloHitFactory));
        }

        if (parameters.m_useLArTPC) {
          // Create U, V and W views assuming x is the common drift coordinate
          const double x0_cm = voxelPos.GetX();
          const double y0_cm = voxelPos.GetY();
          const double z0_cm = voxelPos.GetZ();

          lar_content::LArCaloHitParameters caloHitPars_UView(
              caloHitParameters);
          caloHitPars_UView.m_hitType = pandora::TPC_VIEW_U;
          const double upos_cm(pPrimaryPandora->GetPlugins()
                                   ->GetLArTransformationPlugin()
                                   ->YZtoU(y0_cm, z0_cm));
          caloHitPars_UView.m_positionVector =
              CartesianVector(x0_cm, 0.0, upos_cm);

          lar_content::LArCaloHitParameters caloHitPars_VView(
              caloHitParameters);
          caloHitPars_VView.m_hitType = pandora::TPC_VIEW_V;
          const double vpos_cm(pPrimaryPandora->GetPlugins()
                                   ->GetLArTransformationPlugin()
                                   ->YZtoV(y0_cm, z0_cm));
          caloHitPars_VView.m_positionVector =
              CartesianVector(x0_cm, 0.0, vpos_cm);

          lar_content::LArCaloHitParameters caloHitPars_WView(
              caloHitParameters);
          caloHitPars_WView.m_hitType = pandora::TPC_VIEW_W;
          const double wpos_cm(pPrimaryPandora->GetPlugins()
                                   ->GetLArTransformationPlugin()
                                   ->YZtoW(y0_cm, z0_cm));
          caloHitPars_WView.m_positionVector =
              CartesianVector(x0_cm, 0.0, wpos_cm);

          // Create LArCaloHits for U, V and W views
          PANDORA_THROW_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_UView,
                                          m_larCaloHitFactory));
          PANDORA_THROW_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_VView,
                                          m_larCaloHitFactory));
          PANDORA_THROW_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_WView,
                                          m_larCaloHitFactory));
        }

      } // end voxel loop

      // The voxelisation only works with ArgonCube
      break;

    } // end segment detector loop

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraApi::ProcessEvent(*pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraApi::Reset(*pPrimaryPandora));
  }
}

std::vector<LArVoxel> makeVoxels(const TG4HitSegment &g4Hit,
                                 const LArGrid &grid) {
  std::vector<LArVoxel> currentVoxelList;

  // mm to cm conversion
  const double mm2cm(0.1);

  // Start and end positions
  const TLorentzVector &hitStart = g4Hit.GetStart() * mm2cm;
  const TLorentzVector &hitStop = g4Hit.GetStop() * mm2cm;

  const CartesianVector start(hitStart.X(), hitStart.Y(), hitStart.Z());
  const CartesianVector stop(hitStop.X(), hitStop.Y(), hitStop.Z());

  // Direction vector and hit segment length
  const CartesianVector dir = stop - start;
  const double hitLength = dir.GetMagnitude();

  // Hit segment total energy
  const double g4HitEnergy = g4Hit.GetEnergyDeposit();

  std::cout << "Hit start: " << start << std::endl;
  std::cout << "Hit stop: " << stop << std::endl;
  std::cout << "Hit direction: " << dir << std::endl;
  std::cout << "Hit energy: " << g4HitEnergy << std::endl;
  if (hitLength < 1e-10) {
    std::cout << "Cannot have zero track length" << std::endl;
    return currentVoxelList;
  }

  // Define ray trajectory
  const CartesianVector dirNorm = dir.GetUnitVector();
  LArRay ray(start, dirNorm);

  // We need to shuffle along the hit segment path and create voxels as we go.
  // There are 4 cases for the start and end points inside the voxelisation
  // region. Case 1: start & stop are both inside the voxelisation boundary Case
  // 2: start & stop are both outside, but path direction intersects boundary
  // Case 3: start is inside boundary, stop point = intersection at region
  // boundary Case 4: end is inside boundary, start point = intersection at
  // region boundary

  double t0(0.0), t1(0.0);
  pandora::CartesianVector point1(0.0, 0.0, 0.0), point2(0.0, 0.0, 0.0);

  // Check if the start and end points are inside the voxelisation region
  const bool inStart = grid.inside(start);
  const bool inStop = grid.inside(stop);

  if (inStart && inStop) {
    std::cout << "Case 1: Start and end points are inside boundary"
              << std::endl;
    point1 = start;
    point2 = stop;
  } else if (!inStart && !inStop) {
    std::cout << "Case 2: Start and end points are outside boundary"
              << std::endl;
    bool ok = grid.intersect(ray, t0, t1);
    if (ok) {
      point1 = ray.getPoint(t0);
      point2 = ray.getPoint(t1);
      std::cout << "Boundary points: " << point1 << " and " << point2
                << std::endl;
    } else {
      std::cout << "Voxel ray does not intersect boundary limits." << std::endl;
      return currentVoxelList;
    }
  } else if (inStart && !inStop) {
    std::cout << "Case 3: Start inside boundary" << std::endl;
    point1 = start;
    bool ok = grid.intersect(ray, t0, t1);
    if (ok) {
      point2 = ray.getPoint(t1);
    } else {
      std::cout << "Voxel ray does not intersect boundary" << std::endl;
      return currentVoxelList;
    }
  } else if (!inStart && inStop) {
    std::cout << "Case 4: End inside boundary" << std::endl;
    point2 = stop;
    bool ok = grid.intersect(ray, t0, t1);
    if (ok) {
      point1 = ray.getPoint(t0);
    } else {
      std::cout << "Voxel ray does not intersect boundary" << std::endl;
      return currentVoxelList;
    }
  }

  // Now create voxels between point1 and point2.
  // Ray direction will be the same, but update starting point
  ray.updateOrigin(point1);

  const double epsilon(1e-3);
  bool shuffle(true);

  // Keep track of total voxel path length and energy (so far)
  double totalPath(0.0), ETot(0.0);
  int loop(0);

  while (shuffle) {

    // Get point along path to define voxel bin (bottom corner)
    pandora::CartesianVector voxelPoint = ray.getPoint(epsilon);

    // Grid 3d bin containing this point; 4th element is the total bin number
    std::array<int, 4> gridBins = grid.getBinIndices(voxelPoint);
    int voxelID = gridBins[3];

    int xBin = gridBins[0];
    int yBin = gridBins[1];
    int zBin = gridBins[2];

    // Voxel bottom and top corners
    const pandora::CartesianVector voxBot = grid.getPoint(xBin, yBin, zBin);
    const pandora::CartesianVector voxTop =
        grid.getPoint(xBin + 1, yBin + 1, zBin + 1);

    // Voxel box
    const LArBox vBox(voxBot, voxTop);
    // Get ray intersections with this box
    bool result = vBox.intersect(ray, t0, t1);

    if (!result)
      shuffle = false;

    // Voxel extent = intersection path difference
    double dL(t1 - t0);
    // For the first path length, use the distance from the
    // starting ray point to the 2nd intersection t1
    if (loop == 0)
      dL = t1;

    // Stop processing if we are not moving the path along
    if (dL < epsilon)
      shuffle = false;

    totalPath += dL;

    // Stop adding voxels if we have enough
    if (totalPath > hitLength) {
      shuffle = false;
      // Adjust final path according to hit segment total length
      dL = hitLength - totalPath + dL;
    }

    // Voxel energy using path length fraction w.r.t hit length
    const double voxelEnergy = g4HitEnergy * dL / hitLength;
    ETot += voxelEnergy;

    const double fracTotE = g4HitEnergy > 0.0 ? ETot / g4HitEnergy : 0.0;
    std::cout << "Voxel " << voxelID << ": pos =" << voxBot
              << ", E = " << voxelEnergy << ", fracTotE = " << fracTotE
              << std::endl;

    // Store voxel object in vector
    const LArVoxel voxel(voxelID, voxelEnergy, voxBot);
    currentVoxelList.push_back(voxel);

    // Update ray starting position using intersection path difference
    const pandora::CartesianVector newStart = ray.getPoint(dL);
    ray.updateOrigin(newStart);
    loop++;
  }

  std::cout << "Returning vector of " << currentVoxelList.size()
            << " voxels for G4HitSegment\n"
            << std::endl;

  return currentVoxelList;
}

bool LArBox::intersect(const LArRay &ray, double &t0, double &t1) const {
  // Brian Smits ray-box algorithm with improvements from Amy Williams et al
  double x(ray.m_origin.GetX());
  double y(ray.m_origin.GetY());

  double invDirX(ray.m_invDir.GetX());
  double invDirY(ray.m_invDir.GetY());

  double tMin(0.0), tMax(0.0), tyMin(0.0), tyMax(0.0), tzMin(0.0), tzMax(0.0);

  if (ray.m_sign[0] == 0) {
    tMin = (m_bottom.GetX() - x) * invDirX;
    tMax = (m_top.GetX() - x) * invDirX;
  } else {
    tMin = (m_top.GetX() - x) * invDirX;
    tMax = (m_bottom.GetX() - x) * invDirX;
  }

  if (ray.m_sign[1] == 0) {
    tyMin = (m_bottom.GetY() - y) * invDirY;
    tyMax = (m_top.GetY() - y) * invDirY;
  } else {
    tyMin = (m_top.GetY() - y) * invDirY;
    tyMax = (m_bottom.GetY() - y) * invDirY;
  }

  if ((tMin > tyMax) || (tyMin > tMax)) {
    return false;
  }

  if (tyMin > tMin) {
    tMin = tyMin;
  }

  if (tyMax < tMax) {
    tMax = tyMax;
  }

  double z(ray.m_origin.GetZ());
  double invDirZ(ray.m_invDir.GetZ());
  if (ray.m_sign[2] == 0) {
    tzMin = (m_bottom.GetZ() - z) * invDirZ;
    tzMax = (m_top.GetZ() - z) * invDirZ;
  } else {
    tzMin = (m_top.GetZ() - z) * invDirZ;
    tzMax = (m_bottom.GetZ() - z) * invDirZ;
  }

  if ((tMin > tzMax) || (tzMin > tMax)) {
    return false;
  }

  if (tzMin > tMin) {
    tMin = tzMin;
  }

  if (tzMax < tMax) {
    tMax = tzMax;
  }

  // First and second intersection path lengths
  t0 = tMin;
  t1 = tMax;
  return true;
}

bool LArBox::inside(const pandora::CartesianVector &point) const {
  bool result(false);

  float x = point.GetX();
  float y = point.GetY();
  float z = point.GetZ();

  if (x > m_bottom.GetX() && x < m_top.GetX() && y > m_bottom.GetY() &&
      y < m_top.GetY() && z > m_bottom.GetZ() && z < m_top.GetZ()) {
    result = true;
  }

  return result;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters) {
  if (1 == argc)
    return PrintOptions();

  int c(0);

  std::string viewOption("3d");

  while ((c = getopt(argc, argv, ":i:e:n:s:j:Nh")) != -1) {
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
    case 'j':
      viewOption = optarg;
      break;
    case 'N':
      parameters.m_shouldDisplayEventNumber = true;
      break;
    case 'h':
    default:
      return PrintOptions();
    }
  }

  ProcessViewOption(viewOption, parameters);
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
      << "    -j Projection          (optional) [3D (default), LArTPC, Both]"
      << std::endl
      << std::endl;

  return false;
}

void ProcessViewOption(const std::string &viewOption, Parameters &parameters) {
  std::string chosenViewOption(viewOption);
  std::transform(chosenViewOption.begin(), chosenViewOption.end(),
                 chosenViewOption.begin(), ::tolower);

  if (chosenViewOption == "3d") {
    // 3D hits only
    std::cout << "Using 3D hits" << std::endl;
    parameters.m_useLArTPC = false;
    parameters.m_use3D = true;
  } else if (chosenViewOption == "both") {
    // LArTPC and 3D hits
    std::cout << "Using LArTPC projections _and_ 3D hits" << std::endl;
    parameters.m_useLArTPC = true;
    parameters.m_use3D = true;
  } else {
    // LArTPC hits only
    std::cout << "Using LArTPC projections" << std::endl;
    parameters.m_useLArTPC = true;
    parameters.m_use3D = false;
  }
}

} // namespace lar_nd_reco
