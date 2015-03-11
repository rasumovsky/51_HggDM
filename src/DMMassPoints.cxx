////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMassPoints.cxx                                                    //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This class uses a ROOT file with a TTree to implement the analysis event  //
//  categorization and produce text files or RooDataSets.                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMassPoints.h"

/**
   Initialize the MassPoint class.
*/
DMMassPoints::DMMassPoints() {
  std::cout << std::endl << "DMMassPoints::Initializing..." << std::endl;
  
  nEvents = 0;
    
  int currLineIndex = 0;
  std::string currText;
  char *outFileNameC = (char*)outFileName.c_str();
  // Load output file, configure output TTree:
  outputT3MAPS = new TFile(outFileNameC,"recreate");
  treeT3MAPS = new TTree("TreeT3MAPS","TreeT3MAPS");
  treeT3MAPS->Branch("nHits", &nHits, "nHits/I");
  treeT3MAPS->Branch("timestamp_start", &timestamp_start, "timestamp_start/D");
  treeT3MAPS->Branch("timestamp_stop", &timestamp_stop, "timestamp_stop/D");
  treeT3MAPS->Branch("hit_row", "std::vector<int>", &hit_row);
  treeT3MAPS->Branch("hit_column", "std::vector<int>", &hit_column);
  
  // Open input text file from T3MAPS run:
  char *inFileNameC = (char*)inFileName.c_str();
  ifstream historyFile(inFileNameC);
  if (historyFile.is_open()) {
    while (getline(historyFile, currText) ) {
      
      // Start counting the line numbers (one run is 0-23)
      std::size_t foundText = currText.find("BEGIN SCAN");
      if (foundText!=std::string::npos) {
	currLineIndex = 0;
	hit_row.clear();
	hit_column.clear();
	nHits = 0;
	if (nEvents % 100 == 0) { std::cout << currText << std::endl; }
      }
      
      // start time recorded:
      if (currLineIndex == 2) { 
	timestamp_start = atoi((char*)currText.c_str());
      }
      // stop time recorded:
      else if (currLineIndex == 4) {
	timestamp_stop = atoi((char*)currText.c_str());
      }
      
      // get hit table information:
      else if (currLineIndex > 4 && currLineIndex < 23) {
	int currColumn = currLineIndex - 5;
	
	std::vector<std::string> hitColumns = delimString(currText, " ");
	
	// iterate over the columns that were hit in each row:
	for (std::vector<std::string>::iterator it = hitColumns.begin(); 
	     it != hitColumns.end(); ++it) {
	  int currRow = atoi(it->c_str()) + 1;
	  hit_row.push_back(currRow);
	  hit_column.push_back(currColumn);
	  nHits++;
	}
      }
      
      // end scan, save event information:
      else if (currLineIndex == 23) {
	treeT3MAPS->Fill();
	nEvents++;
      }
            
      // increment the line number:
      currLineIndex++;
    }
  }
  
  historyFile.close();
  treeT3MAPS->Write();
  outputT3MAPS->Close();
  
  return;
}

/**
   Description.
*/
int DMMassPoints::getNEvents() {
  return nEvents;
}
