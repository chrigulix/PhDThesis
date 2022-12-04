// C++ headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

//ROOT header
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPaletteAxis.h"
#include "TAxis.h"

const std::string FilePath = "./BNBFlux/numode_bnb_470m_r200.root";
const std::string UncertaintyFile = "./BNBFlux/bnb_sys_error_uboone.txt";

// Read Beam systematics for nu_mu, anti-nu_mu, nu_e, and anti-nu_e into a TGraph
std::vector<TGraph*> ReadFluxSystematics(const std::string &PathToFile);

void DrawBNBFlux()
{
    // Create uncertainty graphs
    std::vector<TGraph*> UncertaintyGraphs = ReadFluxSystematics(UncertaintyFile);
    
    // Open file
    TFile* BNBFluxFile = new TFile(FilePath.c_str());
    
    // Initialize Vectors
    std::vector< std::pair<std::string,int> > Flavors;
    Flavors.push_back(std::make_pair("numu",kRed+2));
    Flavors.push_back(std::make_pair("numubar",kRed-4));
    Flavors.push_back(std::make_pair("nue",kBlue+2));
    Flavors.push_back(std::make_pair("nuebar",kBlue-4));
//     Flavors.push_back(std::make_pair("numu",46));
//     Flavors.push_back(std::make_pair("numubar",30));
//     Flavors.push_back(std::make_pair("nue",38));
//     Flavors.push_back(std::make_pair("nuebar",40));
    
    std::vector<TH1D*> FluxHistogramsError;
    std::vector<TH1D*> FluxHistograms;
    std::vector<double> TotalFlux;
    
    double FluxAllNeutrinos = 0.0;
    
    double iter = 0;
    
    // Fill histograms and handle their content
    for(auto nu : Flavors)
    {
        // Fill histograms from file into memory
        FluxHistogramsError.push_back( (TH1D*)BNBFluxFile -> Get(nu.first.c_str()) );
        FluxHistograms.push_back( (TH1D*)BNBFluxFile -> Get(nu.first.c_str()) );
        
        // Scale Histos to bin width and per POT
        double BinWidth = FluxHistogramsError.back() -> GetBinWidth(2);
        FluxHistogramsError.back() -> Scale(1/BinWidth/1e20);
        FluxHistograms.back() -> Scale(1/BinWidth/1e20);
        
        // Fill total flux
        TotalFlux.push_back(FluxHistogramsError.back() -> Integral()*BinWidth);
        FluxAllNeutrinos += TotalFlux.back();
        
        // Integrated flux output
        std::cout << nu.first + " flux: " << TotalFlux.back() << " cm^-2 POT^-1" << std::endl;
        
        // Histogram drawing style and options
        FluxHistogramsError.back() -> SetLineColorAlpha(nu.second,1);
        FluxHistogramsError.back() -> SetFillColorAlpha(nu.second,0.35);
        FluxHistogramsError.back() -> SetMarkerColorAlpha(nu.second,0);
        FluxHistogramsError.back() -> SetLineWidth(2.);
        
        FluxHistograms.back() -> SetLineColor(nu.second);
        FluxHistograms.back() -> SetLineWidth(2.);
        
        double BinContent;
        double BinCenter;
        double BinError;
        
        // Loop though bins
        for(unsigned int bin = 1; bin <= FluxHistogramsError.back()->GetNbinsX(); bin++)
        {
            // Get Histogram information
            BinContent = FluxHistogramsError.back() -> GetBinContent(bin);
            BinCenter = FluxHistogramsError.back() -> GetBinCenter(bin);
            
            // Calculate bin uncertainty
            BinError = FluxHistogramsError.back()->GetBinError(bin,BinError) + BinContent * UncertaintyGraphs.at(iter)->Eval(BinCenter);
            
            FluxHistogramsError.back() -> SetBinError(bin,BinError);
            FluxHistograms.back() -> SetBinError(bin,0);
        }
        
        iter++;
    }

    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "Total Flux: " << FluxAllNeutrinos << " cm^-2 POT^-1" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    
    unsigned int index = 0;
    // Make flux ratios
    for(auto nuFlux : TotalFlux)
    {
        // Flux ratio output
        std::cout << Flavors.at(index).first + " ratio: " << nuFlux/FluxAllNeutrinos*100 << " %" << std::endl;
        index++;
    }
    
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "Mean Energy: " << FluxHistograms.front() -> GetMean() << " GeV" << std::endl;
    std::cout << "Maximum Flux: " << FluxHistograms.front() -> GetBinContent(FluxHistograms.front()->GetMaximumBin()) << " cm^-2 GeV^-1 POT^-1" << std::endl;
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.60,0.70,0.80,0.80);
    Legend->SetNColumns(2);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.5 );
    Legend->SetEntrySeparation(0.3);
    Legend->SetTextSize(0.03);
    Legend->SetHeader("Neutrino Flavour","L");
    Legend->SetTextAlign(12);
    
    // Change header size
    TLegendEntry *header = (TLegendEntry*)Legend->GetListOfPrimitives()->First();
    header->SetTextSize(0.04);
    
    // Fill Legend 
    Legend -> AddEntry(FluxHistogramsError.at(0), " #nu_{#mu}", "LF") ;
    Legend -> AddEntry(FluxHistogramsError.at(1), " #bar{#nu}_{#mu}", "LF") ;
    Legend -> AddEntry(FluxHistogramsError.at(2), " #nu_{e}", "LF") ;
    Legend -> AddEntry(FluxHistogramsError.at(3), " #bar{#nu}_{e}", "LF") ;
    
    // This is needed because root is stupid
    TH1D* GhostHist = new TH1D("GhostHist", "GhostHist", 200,-0.2, 7);
    GhostHist -> SetTitle("BNB Flux in Neutrino Mode");
    GhostHist -> GetXaxis() -> SetTitle("Neutrino Energy E_{#nu} [GeV]");
    GhostHist -> GetYaxis() -> SetTitle("BNB Flux #Phi(E_{#nu}) [cm^{-2} GeV^{-1} POT^{-1}]");
    GhostHist -> GetXaxis() -> SetTitleOffset(0.9);
    GhostHist -> GetYaxis() -> SetTitleOffset(1.1);
    GhostHist -> GetXaxis() -> SetRangeUser(-0.2,7.0);
    GhostHist -> GetYaxis() -> SetRangeUser(2e-16,2e-9);
    
    // Draw histograms
    TCanvas* C0 = new TCanvas("C0","BNB Flux",0,0,1400,1000);
    C0 -> SetLogy();
    C0 -> SetLeftMargin(0.11);
    C0 -> SetRightMargin(0.9);
    GhostHist -> Draw();
    for(unsigned int hist = 0; hist < FluxHistograms.size(); hist++)
    {
        FluxHistogramsError.at(hist) -> Draw("E2 SAME");
        FluxHistograms.at(hist) -> Draw("][ SAME");
    }
    gPad->RedrawAxis();
    Legend -> Draw();
    C0 -> SaveAs("../images/MicroBooNE/BNBFluxMicroBooNE.pdf");
    
}

std::vector<TGraph*> ReadFluxSystematics(const std::string &PathToFile)
{
    std::ifstream InputFile(PathToFile);
    
    // Readout variables
    double Energy;
    double nu_mu;
    double nu_mubar;
    double nu_e;
    double nu_ebar;
    std::string Line;
    
    // Initialize data vectors
    std::vector<double> EnergyData;
    std::vector<double> nu_muData;
    std::vector<double> nu_mubarData;
    std::vector<double> nu_eData;
    std::vector<double> nu_ebarData;
    
    // Skip header
    std::getline(InputFile,Line);
    
    // Read variables
    while(InputFile >> Energy >> nu_mu >> nu_mubar >> nu_e >> nu_ebar)
    {
        EnergyData.push_back(Energy);
        nu_muData.push_back(nu_mu);
        nu_mubarData.push_back(nu_mubar);
        nu_eData.push_back(nu_e);
        nu_ebarData.push_back(nu_ebar);
    }
    
    // Convert, because root sucks
    double* EnergyEntry = EnergyData.data();
    double* nu_muEntry =nu_muData.data();
    double* nu_mubarEntry = nu_mubarData.data();
    double* nu_eEntry = nu_eData.data();
    double* nu_ebarEntry =nu_ebarData.data();
    
    std::vector<TGraph*> Graphs;
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_muEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_mubarEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_eEntry));
    Graphs.push_back(new TGraph(EnergyData.size(),EnergyEntry,nu_ebarEntry));
    
    return Graphs;
}
