#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <deque>
#include <iostream>
#include <cmath>

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>
#include <TSpline.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

TSpline3* KEvsRSpline; // Global spline for momentum calculation
const double Density = 1.396; // g/cm^3

// Momentum calculation
void MomentumSplinePreparation();

// Get Momentum
float GetMomentum(float TrackLength);

void DrawLengthToMomentum()
{
    // Fill momentum calculation spline
    MomentumSplinePreparation();
    
    // TGraph track length to energy conversion
    TGraph* LengthToEnergy = new TGraph();
    
    const unsigned int MaxBin = 10000;
    float MaxLength = 1000.0; //cm
    float CurrentLength = 0; //cm
    
    float x_value[MaxBin];
    float y_value[MaxBin];
    
    // Loop over all graph bins
    for(unsigned int bin_no = 0; bin_no < MaxBin; bin_no++)
    {        
        // Fill Points
        x_value[bin_no] = CurrentLength;
        y_value[bin_no] = GetMomentum(CurrentLength);
        
        // Add track length increments
        CurrentLength += MaxLength/(double)MaxBin;
    }
    
    // TGraph track length to energy conversion and fill it
    TGraph* LengthToMomentum = new TGraph(MaxBin,x_value,y_value);
    TGraph* LengthToMomentumMini = new TGraph(MaxBin,x_value,y_value);
    
    // TGraph momentum to track length conversion and fill it
    TGraph* MomentumToLength = new TGraph(MaxBin,y_value,x_value);
    TGraph* MomentumToLengthMini = new TGraph(MaxBin,y_value,x_value);
    
    TCanvas *C0 = new TCanvas("C0", "C0", 1400, 1000);
    LengthToMomentum -> GetXaxis() -> SetTitle("Muon Track Length [cm]");
    LengthToMomentum -> GetYaxis() -> SetTitle("Muon Momentum [GeV/c]");
    LengthToMomentum -> SetTitle("Track Length to Momentum Conversion");
    LengthToMomentum -> GetHistogram() -> GetXaxis() -> SetRangeUser(0,1000);
    LengthToMomentum -> GetHistogram() -> GetYaxis() -> SetRangeUser(0,2.5);
    LengthToMomentum -> SetLineColor(39);
    LengthToMomentum -> SetLineWidth(4);
    LengthToMomentum -> Draw("AC");
    gPad->RedrawAxis();
    TPad *subpad0 = new TPad("subpad0","",0.12,0.52,0.52,0.88);
    subpad0 -> SetTopMargin(0.05);
    subpad0 -> SetBottomMargin(0.05);
    subpad0 -> SetLeftMargin(0.07);
    subpad0 -> SetRightMargin(0.05);
    subpad0 -> Draw();
    subpad0 -> cd();
    LengthToMomentumMini -> SetLineColor(39);
    LengthToMomentumMini -> SetLineWidth(4);
    LengthToMomentumMini -> SetTitle("");
    LengthToMomentumMini -> GetHistogram() -> GetXaxis() -> SetRangeUser(0,20);
    LengthToMomentumMini -> GetHistogram() -> GetYaxis() -> SetRangeUser(0,0.15);
    LengthToMomentumMini -> Draw("AC");
    subpad0->RedrawAxis();
    C0->SaveAs("../images/FirstCCInclusive/TrackLengthMomentumRelation.pdf");
    
    TCanvas *C1 = new TCanvas("C1", "C1", 1400, 1000);
    C1 -> SetLogy();
    MomentumToLength -> GetXaxis() -> SetTitle("Muon Momentum [GeV/c]");
    MomentumToLength -> GetYaxis() -> SetTitle("Muon Track Length [cm]");
    MomentumToLength -> SetTitle("Momentum to Track Length Conversion");
    MomentumToLength -> GetHistogram() -> GetXaxis() -> SetRangeUser(0,2.5);
    MomentumToLength -> GetHistogram() -> GetYaxis() -> SetRangeUser(0.1,1000);
    MomentumToLength -> SetLineColor(8);
    MomentumToLength -> SetLineWidth(4);
    MomentumToLength -> Draw();
    gPad->RedrawAxis();
//     TPad *subpad1 = new TPad("subpad1","",0.52,0.12,0.88,0.52);
//     subpad1 -> SetTopMargin(0.05);
//     subpad1 -> SetBottomMargin(0.05);
//     subpad1 -> SetLeftMargin(0.07);
//     subpad1 -> SetRightMargin(0.05);
//     subpad1 -> Draw();
//     subpad1 -> cd();
//     MomentumToLengthMini -> SetLineColor(39);
//     MomentumToLengthMini -> SetLineWidth(4);
//     MomentumToLengthMini -> SetTitle("");
//     MomentumToLengthMini -> GetHistogram() -> GetXaxis() -> SetRangeUser(0,0.2);
//     MomentumToLengthMini -> GetHistogram() -> GetYaxis() -> SetRangeUser(0,40);
//     MomentumToLengthMini -> Draw("AC");
//     subpad1->RedrawAxis();
//     C1->SaveAs("../images/FirstCCInclusive/TrackLengthMomentumRelation.pdf");
    
}

void MomentumSplinePreparation()
{
    float RangeGramPerCM[29] = {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
                                2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3, 4.359E3, 5.354E3, 7.298E3,
                                1.013E4, 1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5
                               };

    float KEMeV[29] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000, 4000,
                       8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000
                      };

    // convert to cm
    for(auto & RangePoint : RangeGramPerCM)
    {
        RangePoint /= Density;
    }

    TGraph* KEvsR = new TGraph(29, RangeGramPerCM, KEMeV);
//     KEvsR -> Draw();

    KEvsRSpline = new TSpline3("KEvsRS",KEvsR);

    delete KEvsR;
}

float GetMomentum(float TrackLength)
{
    float MuonMass = 105.7; //MeV

    // Change Track length to kinetic energy
    TrackLength = KEvsRSpline->Eval(TrackLength);

    // Convert kinetic energy to momentum
    TrackLength = std::sqrt( std::pow(TrackLength,2) + 2*TrackLength*MuonMass );

    // Convert MeV to GeV
    TrackLength /= 1000;

    return TrackLength;
}
