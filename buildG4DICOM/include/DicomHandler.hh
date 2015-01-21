//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  GAMOS software  is  copyright of the Copyright  Holders  of *
// * the GAMOS Collaboration.  It is provided  under  the  terms  and *
// * conditions of the GAMOS Software License,  included in the  file *
// * LICENSE and available at  http://fismed.ciemat.es/GAMOS/license .*
// * These include a list of copyright holders.                       *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GAMOS collaboration.                       *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the GAMOS Software license.           *
// ********************************************************************
//
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************//

//*******************************************************
//
// DicomHandler.hh :
//	- Handling of DICM images
//	- Transforming *.dcm to *.g4 ( pixels->density )
//	- Reading headers and pixels
//	- Transforming pixel to density and creating *.g4
//	  files
//	- Functions are in DicomHandler.cc
//
//
// Base on previous code by :
//	Dragan Tubic <tdragan@gel.ulaval.ca>
//*******************************************************

#ifndef DicomHandler_h
#define DicomHandler_h 1

#include <cstdio>
#include <map>
#include <vector>
#include <fstream>

#include "globals.hh"

class DicomPhantomZSliceHeader;

class DicomHandler
{
public:

  DicomHandler();
  
  ~DicomHandler();

  void CheckFileFormat();
    
  G4int ReadFile(FILE *);
  G4int ReadData(FILE *); // note: always use readHeader 
  // before readData
  
  // use ImageMagick to display the image
  //G4int displayImage(char[500]);
  
private:
  void ReadMaterialIndices( std::ifstream& finData);
  void GetInformation(G4int &, char *);
  void WriteTextZSlice(char *); // note: always use readHeader 
  void WriteTextZSliceMaterialIndices(std::ofstream& foutG4DCM, G4int** HUs);
  void WriteTextZSliceDensities(std::ofstream& foutG4DCM, G4int** HUs);

  unsigned int GetMaterialIndex( G4float density );

  void WriteBinZSlice(char *); // note: always use readHeader 
  void WriteBinZSliceMaterialIndices( FILE* fileOut, G4int** HUs );
  void WriteBinZSliceDensities( FILE* fileOut, G4int** HUs );
  
  void WriteTextAllSlices( DicomPhantomZSliceHeader* ZSliceHeaderMerged );
  void WriteBinAllSlices( DicomPhantomZSliceHeader* ZSliceHeaderMerged );

  G4float Pixel2density(G4int pixel);
  void WriteZSliceData(std::ofstream& foutG4DCM);
  G4int read_defined_nested(FILE *, G4int);
  void read_undefined_nested(FILE *);
  void read_undefined_item(FILE *);

  template <class Type> void GetValue(char *, Type &);

  void ReadPixel2Density();

private:
  const int DATABUFFSIZE;
  const int LINEBUFFSIZE;
  const int FILENAMESIZE;
    
  short compression;
  G4int nFiles;
  short rows;
  short columns;
  short bitAllocated;
  G4int maxPixelValue, minPixelValue;
  
  G4double pixelSpacingX, pixelSpacingY;
  G4double sliceThickness;
  G4double sliceLocation;
  
  G4int rescaleIntercept, rescaleSlope;
  
  G4bool littleEndian, implicitEndian;
  short pixelRepresentation;
  
  std::vector<G4int**> theHUs_ALL; // the Hounsfield number of all slices
  G4int** theHUs; // the Hounsfield numbers of 1 slice
  std::map<G4float,G4String> fMaterialIndices;

  G4bool bMergeZSlices;

  G4int theCT2D_NValues;
  G4double * theCT2D_HU;
  G4double * theCT2D_Density;
  std::ifstream* theCalibrationFile;

};
#endif

