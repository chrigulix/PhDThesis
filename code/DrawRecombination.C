// Recombination in liquid argon

// C++ headers
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const double alpha = 0.70;
const double A =  0.8;
const double k = 0.0486; // kV/cm g/cm^2/MeV
const double rho = 1.396; // g/cm^3
const double dEdx_mip = 2.105; // MeV/cm

const unsigned long int SampleNo = 1000;
const double x_min = 0.00000001;
const double x_max = 2;
const double y_min = 0;
const double y_max = 1;



double ChargeFunction(double EField, double dEdx)
{
    double R_c = A/(1+k/EField*dEdx);
    return R_c;
}

double LumiFunction(double EField, double dEdx)
{   
    double R_l = 1-alpha*ChargeFunction(EField,dEdx);
    return R_l;
}

void DrawRecombination()
{
    // Create dEdx vector with x*dE/dx_mip
    std::vector<double> dEdxVec = {1.,2.,5.,10.,20};
    std::vector<unsigned int> Color = {46,38,30,28,kViolet-8};
    
    // Create vector containing different functions for different particles
    std::vector<TGraph*> ChargeGraphs;
    std::vector<TGraph*> LumiGraphs;
    
    
    // Create tick length
    double XAxisTic = (x_max - x_min)/(double)(SampleNo-1);
    
    // Create temporary arrays
    double XValue[SampleNo];
    double YValueCharge[SampleNo];
    double YValueLumi[SampleNo];
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.5,0.15,0.85,0.35);
    Legend->SetNColumns(2);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.5);
    Legend->SetTextSize(0.03);
    Legend->SetHeader("#LT-dE/dx#GT =","C");
    Legend->SetTextAlign(12);
    
    TLegendEntry *header = (TLegendEntry*)Legend->GetListOfPrimitives()->First();
    header->SetTextSize(0.04);
    
    unsigned long int index = 0;
    // Loop over all dE/dx entries and fill Graphs
    for(auto dEdx : dEdxVec)
    {
        // Scale entries to mip dE/dx and MeV cm^2/g
        dEdx *= dEdx_mip/rho;
        
        for(unsigned long int i = 0; i < SampleNo; i++)
        {
            // Calculate x-values
            XValue[i] = x_min + (double)i*XAxisTic;
        
            // Use Function to calculate y-values
            YValueCharge[i] = ChargeFunction(XValue[i],dEdx);
            YValueLumi[i] = LumiFunction(XValue[i],dEdx);
        }
        // Fill graphs
        ChargeGraphs.push_back( new TGraph(SampleNo,XValue,YValueCharge) );
        ChargeGraphs.back() -> SetLineColor(Color.at(index));
        
        LumiGraphs.push_back( new TGraph(SampleNo,XValue,YValueLumi) );
        LumiGraphs.back() -> SetLineColor(Color.at(index));
        LumiGraphs.back() -> SetLineStyle(7);
        
        // Fill legend
//         if(std::lround(dEdxVec.at(index)) < 10)
//         {
//             Legend -> AddEntry( ChargeGraphs.back(),("   " + std::to_string( std::lround(dEdxVec.at(index))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
//         }
//         else
//         {
//             Legend -> AddEntry( ChargeGraphs.back(),(" " + std::to_string( std::lround(dEdxVec.at(index))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
//         }
        
        index++;
    }
    
    // Fill legend manually because root is stupid
    Legend -> AddEntry( ChargeGraphs.at(0),(" " + std::to_string( std::lround(dEdxVec.at(0))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
    Legend -> AddEntry( ChargeGraphs.at(3),(" " + std::to_string( std::lround(dEdxVec.at(3))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
    Legend -> AddEntry( ChargeGraphs.at(1),(" " + std::to_string( std::lround(dEdxVec.at(1))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
    Legend -> AddEntry( ChargeGraphs.at(4),(" " + std::to_string( std::lround(dEdxVec.at(4))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
    Legend -> AddEntry( ChargeGraphs.at(2),(" " + std::to_string( std::lround(dEdxVec.at(2))) + " #times #LT-dE/dx#GT_{MIP}" ).c_str(),"L" );
    
    // Setup 2nd Legend
    TLegend* Legend2 = new TLegend(0.5,0.75,0.85,0.85);
    Legend2->SetNColumns(1);
    Legend2->SetLineStyle ( 0 );
    Legend2->SetLineColorAlpha ( 0,0 );
    Legend2->SetFillStyle ( 0 );
    Legend2->SetMargin ( 0.2 );
    Legend2->SetEntrySeparation(0.5);
    Legend2->SetTextSize(0.04);
//     Legend2->SetHeader("#LT-dE/dx#GT =","C");
    Legend2->SetTextAlign(12);
    
    // Empty TGraph for legend
    TGraph* PhotonicSignal = new TGraph();
    PhotonicSignal -> SetLineColor(12);
    PhotonicSignal -> SetLineStyle(7);
    TGraph* ChargeSignal = new TGraph();
    ChargeSignal -> SetLineColor(12);
    Legend2 -> AddEntry(PhotonicSignal," Luminescence fraction L/L_{0}", "L");
    Legend2 -> AddEntry(ChargeSignal," Charge fraction Q/Q_{0}", "L");
    
    // Add Legend entries
    
    // Create and fill multi graph
    TMultiGraph* MultiGraph = new TMultiGraph();
    
    // Loop over lumi graphs
    for(auto Graph : LumiGraphs)
    {
        // Fill the multi graph
        MultiGraph->Add(Graph);
    }
    
    // Loop over charge graphs
    for(auto Graph : ChargeGraphs)
    {
        // Fill the multi graph
        MultiGraph->Add(Graph);
    }
    
    // Multigraph settings
    MultiGraph -> SetTitle("Recombination in LAr");
    MultiGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(0,x_max);
    MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    MultiGraph -> GetXaxis() -> SetTitle("Drift Electric Field Strength [kV cm^{-1}]");
    MultiGraph -> GetXaxis() -> SetTitleOffset(0.9);
    MultiGraph -> GetYaxis() -> SetTitle("Q/Q_{0} and L/L_{0} [ ]");
    MultiGraph -> GetYaxis() -> SetTitleOffset(0.9);

    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Recombination",0,0,1400,1000);
    
    MultiGraph -> Draw("AC");
    Legend -> Draw();
    Legend2 -> Draw();
    
    C0 -> SaveAs("../images/Detector/Recombination.pdf");
    
//     TCanvas* C1 = new TCanvas("C1","Drift Velocity",0,0,1400,1000);
    
//     DiffusionGraphs.at(1) -> Draw("AC");
//     Legend -> Draw("SAME");
    
//     C1 -> SaveAs("../images/Detector/Diffusion2.pdf");
}
