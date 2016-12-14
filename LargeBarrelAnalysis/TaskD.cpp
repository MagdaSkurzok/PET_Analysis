/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file TaskD.cpp
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cctype>
#include <JPetWriter/JPetWriter.h>
#include "TaskD.h"
#include "TF1.h"
#include "TString.h"
#include <TDirectory.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

TaskD::TaskD(const char * name, const char * description):JPetTask(name, description){}

void TaskD::init(const JPetTaskInterface::Options& opts){

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fBarrelMap.buildMappings(getParamBank());
// create histograms for time differences at each slot and each threshold

	for(auto & scin : getParamBank().getScintillators()){//over slots

		for (int thr=1;thr<=4;thr++){//over threshold

			const char * histo_name = formatUniqueSlotDescription(scin.second->getBarrelSlot(), thr, "timeDiffAB_");
			getStatistics().createHistogram( new TH1F(histo_name, histo_name, 100, -20., 20.) );
		}
	}




// create histograms for time diffrerence vs slot ID

	for(auto & layer : getParamBank().getLayers()){

		for (int thr=1;thr<=4;thr++){

			const char * histo_name = Form("TimeDiffVsID_layer_%d_thr_%d", fBarrelMap.getLayerNumber(*layer.second), thr);
			const char * histo_titile = Form("%s;Slot ID; TimeDiffAB [ns]", histo_name); 

			int n_slots_in_layer = fBarrelMap.getNumberOfSlots(*layer.second);

			getStatistics().createHistogram( new TH2F(histo_name, histo_titile, n_slots_in_layer, 0.5, n_slots_in_layer+0.5,
								  120, -20., 20.) );
		}
	}


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TaskD::exec(){
	//getting the data from event in propriate format
	if(auto hit =dynamic_cast<const JPetHit*const>(getEvent())){
		fillHistosForHit(*hit);
		fWriter->write(*hit);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TaskD::terminate(){

//getting Layer and Slots numbers from analysing files name

//gDirectory->GetFile()->GetName(); //path and name of output file example: 
//zwraca np. "/home/mskurzok/j-pet-framework-examples/build/data_ref/Layer1_Slot16.phys.hit.means.root"

std::vector<int> vectorOfNumbers;
std::string str(gDirectory->GetFile()->GetName());//getting directory of analysing file
int number=0;

for (unsigned int i=0; i < str.size(); i++){

	if (isdigit(str[i])){

          std::stringstream ss;
          ss<<str[i];
          ss>>number; //convert string into int and store it in "asInt"

//std::cout << "number" << number << std::endl;

         vectorOfNumbers.push_back(number);

	}
}

 //int a = vectorOfNumbers[0];
 //int b = vectorOfNumbers[1];

 //std::cout << a << std::endl;  
 //std::cout << b << std::endl; 

//////////////////////////////////////////////////////////////////////////////////////////////////////

// save timeDiffAB mean values for each slot and each threshold in a JPetAuxilliaryData object
	// so that they are available to the consecutive modules
	getAuxilliaryData().createMap("timeDiffAB mean values");


//create output txt file with calibration parameters 
 
	std::ofstream results_fit;
	results_fit.open("results.txt", std::ios::app); //file will be overwrote

	for(auto & slot : getParamBank().getBarrelSlots()){

		for (int thr=1;thr<=4;thr++){

			const char * histo_name = formatUniqueSlotDescription(*(slot.second), thr, "timeDiffAB_");
			double mean = getStatistics().getHisto1D(histo_name).GetMean();
			getAuxilliaryData().setValue("timeDiffAB mean values", histo_name, mean);

			TH1F* histoToSave = &(getStatistics().getHisto1D(histo_name));

//non zero histos 
// slot.first - ID
// slot.second - wskaznik na JPetBarrelSlot
//save fit parameters only for layerX and SlotY
//fit just for proper slot

			if(histoToSave->GetEntries() != 0 && (slot.second)->getLayer().getId()== vectorOfNumbers[0] && (slot.first==(10*vectorOfNumbers[1]+vectorOfNumbers[2]))){

			//if(histoToSave->GetEntries() != 0 && (slot.second)->getLayer().getId()== 1 && (slot.first==20)){
			int highestBin = histoToSave->GetBinCenter(histoToSave->GetMaximumBin());
			histoToSave->Fit("gaus","","", highestBin-5, highestBin+5);
			TCanvas* c = new TCanvas();
			histoToSave->Draw();
			std::string sHistoName = histo_name; sHistoName+=".png";
			//c->SaveAs( sHistoName.c_str());
		
			TF1 *fit = histoToSave->GetFunction("gaus");

			double position_peak = fit->GetParameter(1);
   			double position_peak_error=fit->GetParError(1);
			double sigma_peak =fit->GetParameter(2);
			double chi2_ndf = fit->GetChisquare()/fit->GetNDF();

				results_fit <<  (slot.second)->getLayer().getId() << "\t" <<  slot.first << "\t" << thr << "\t" << position_peak << "\t" << position_peak_error << "\t" << sigma_peak << "\t" << chi2_ndf << "\t"  << std::endl;
			}
			
		}
	}

results_fit.close();
	
}

//const char * histo_name = Form("TimeDiff_layer_%d_slot_%d_thr_%d", LayerNR, SlotNR, thr);
//getStatistics().createHistogram( new TH1F(histo_name, histo_name, 2000, -20., 20.));

//////////////////////////////////

void TaskD::fillHistosForHit(const JPetHit & hit){

	auto lead_times_A = hit.getSignalA().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);



	//if(hit.getBarrelSlot().getID()!=193){
	auto lead_times_B = hit.getSignalB().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
	//}

	//if(hit.getBarrelSlot().getID()==193){
	//auto lead_times_C = hit.getSignalB().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
	//}


	for(auto & thr_time_pair : lead_times_A){
		int thr = thr_time_pair.first;

		if( lead_times_B.count(thr) > 0 ){ // if there was leading time at the same threshold at opposite side

			double timeDiffAB = lead_times_A[thr] - lead_times_B[thr];
			timeDiffAB /= 1000.; // we want the plots in ns instead of ps

			// fill the appropriate histogram
			const char * histo_name = formatUniqueSlotDescription(hit.getBarrelSlot(), thr, "timeDiffAB_");
			getStatistics().getHisto1D(histo_name).Fill( timeDiffAB );

			// fill the timeDiffAB vs slot ID histogram
			int layer_number = fBarrelMap.getLayerNumber( hit.getBarrelSlot().getLayer() );
			int slot_number = fBarrelMap.getSlotNumber( hit.getBarrelSlot() );
			getStatistics().getHisto2D(Form("TimeDiffVsID_layer_%d_thr_%d", layer_number, thr)).Fill( slot_number,
														  timeDiffAB);
		}
	}




}


const char * TaskD::formatUniqueSlotDescription(const JPetBarrelSlot & slot, int threshold, const char * prefix = ""){

	int slot_number = fBarrelMap.getSlotNumber(slot);
	int layer_number = fBarrelMap.getLayerNumber(slot.getLayer()); 

	return Form("%slayer_%d_slot_%d_thr_%d",prefix,layer_number,slot_number,threshold);

}
void TaskD::setWriter(JPetWriter* writer){fWriter =writer;}
