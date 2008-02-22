// -*- C++ -*-
//
// Package:    CaloMiscalibTools
// Class:      CaloMiscalibTools
// 
/**\class CaloMiscalibTools CaloMiscalibTools.h CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibTools.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Lorenzo AGOSTINO
//         Created:  Wed May 31 10:37:45 CEST 2006
// $Id: CaloMiscalibTools.cc,v 1.1 2007/07/14 21:25:58 malgeri Exp $
//
// Modified       : Luca Malgeri 
// Date:          : 11/09/2006 
// Reason         : split class definition (.h) from source code (.cc)
//
//


// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibTools.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalBarrel.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalEndcap.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibMapEcal.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

//
// constructors and destructor
//
CaloMiscalibTools::CaloMiscalibTools(const edm::ParameterSet& iConfig)
{
   //the following line is needed to tell the framework what
   // data is being produced
  map_.prefillMap();

  barrelfileinpath_=iConfig.getUntrackedParameter<std::string> ("fileNameBarrel","");
  endcapfileinpath_=iConfig.getUntrackedParameter<std::string> ("fileNameEndcap","");

  edm::FileInPath barrelfiletmp("CalibCalorimetry/CaloMiscalibTools/data/"+barrelfileinpath_);
  edm::FileInPath endcapfiletmp("CalibCalorimetry/CaloMiscalibTools/data/"+endcapfileinpath_);
  
  
  barrelfile_=barrelfiletmp.fullPath();
  endcapfile_=endcapfiletmp.fullPath();

   // added by Zhen (changed since 1_2_0)
   setWhatProduced(this,&CaloMiscalibTools::produce);
   findingRecord<EcalIntercalibConstantsRcd>();
   //now do what ever other initialization is needed
}


CaloMiscalibTools::~CaloMiscalibTools()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
CaloMiscalibTools::ReturnType
CaloMiscalibTools::produce(const EcalIntercalibConstantsRcd& iRecord)
{
    map_.prefillMap();
    MiscalibReaderFromXMLEcalBarrel barrelreader_(map_);
    MiscalibReaderFromXMLEcalEndcap endcapreader_(map_);
    if(!barrelfile_.empty()) barrelreader_.parseXMLMiscalibFile(barrelfile_);
    if(!endcapfile_.empty())endcapreader_.parseXMLMiscalibFile(endcapfile_);
    map_.print();
    // Added by Zhen, need a new object so to not be deleted at exit
    //    std::cout<<"about to copy"<<std::endl;
    EcalIntercalibConstants* mydata=new EcalIntercalibConstants(map_.get());
    //    std::cout<<"mydata "<<mydata<<std::endl;
    return mydata;
}

 void CaloMiscalibTools::setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&,
 edm::ValidityInterval & oValidity)
 {
 oValidity = edm::ValidityInterval(edm::IOVSyncValue::beginOfTime(),edm::IOVSyncValue::endOfTime());
 
 }

