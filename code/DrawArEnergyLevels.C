// Ion mobility in liquid argon

// C++ headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const double x_min = 0;
const double x_max = 1;
const double y_min = 0;
const double y_max = 16;

const double Opacity = 0.5;

const double h = 6.62607015e-34; // Js
const double c = 299792458; // ms^-1
const double e = 1.602176634e-19; // C

const double I_LAr = 13.4; // eV

const double LaserPhotonEnergy = 4.66; // eV

void DrawArEnergyLevels() 
{
    // Create vector containing relevant particle information
    std::vector<std::pair<double, double>> EnergyLevels;
    
    // Open file with levels
    std::ifstream InputFile("ArgonEnergyLevels");
    std::string Line;
    
    // Read file in loop line by line and fill the data into container
    double temp1;
    double temp2;
    
    // Band width in liquid
    double BWidth = 0.08;
    
    // Scale factor for cm^-1 -> eV unit conversion
    double ScaleFactor = c*h/e*100; // scale factor for cm^-1 -> eV unit transition
    
    while(InputFile >> temp1 >> temp2)
    {
        EnergyLevels.push_back(std::make_pair(temp1,temp2*ScaleFactor));
    }
    // close input file
    InputFile.close();
    
    // Debug file reading
//     for(auto Level : EnergyLevels)
//     {
//         std::cout.precision(10);
//         std::cout << Level.first << " " << Level.second << std::endl;
//     }
    
    // Calculate Ionisation potential difference of gas and liquid argon
    double DeltaE = I_LAr - EnergyLevels.back().second;
    
    // Create root TLine containers
    std::vector<TLine*> GasLines;
    std::vector<TLine*> LiquidLines;
    std::vector<TBox*> LiquidBands;
    std::vector<TBox*> LaserBands;
    
    // Iterator for writing results to console
    unsigned int iterator = 0;
    
    for(auto Level : EnergyLevels)
    {
        GasLines.push_back(new TLine(0.1,Level.second,0.95,Level.second));
        
        if (Level.second > 0.1) 
        {
            LiquidLines.push_back(new TLine(0.1,Level.second+DeltaE,0.95,Level.second+DeltaE));
            LiquidBands.push_back(new TBox(0.1,Level.second+DeltaE-BWidth,0.95,Level.second+DeltaE+BWidth));
            LaserBands.push_back(new TBox(0.1,Level.second+DeltaE-BWidth,0.95,Level.second+DeltaE+BWidth));
        }
        else
        {
            LiquidLines.push_back(new TLine(0.1,Level.second,0.95,Level.second));
            LiquidBands.push_back(new TBox(0.1,Level.second-0.2,0.95,Level.second));
            LaserBands.push_back(new TBox(0.1,Level.second-0.2,0.95,Level.second));
        }
        // Colour even and odd spins differently
        if((int)Level.first % 2)
        {
            // odd spins in red
            GasLines.back() -> SetLineColor(46);
            LiquidLines.back() -> SetLineColor(46);
            LiquidBands.back() -> SetFillColorAlpha(46,Opacity);
            
            // Remove all odd spin excitation bands
            LaserBands.pop_back();
        }
        else
        {
            // even spins in blue
            GasLines.back() -> SetLineColor(38);
            LiquidLines.back() -> SetLineColor(38);
            LiquidBands.back() -> SetFillColorAlpha(38,Opacity);
            
            // Eliminate not allowed two photon transitions with J > 2, colour the others
            if(Level.first > 2.1)
            {
//                 LaserBands.back() -> SetFillColorAlpha(4,Opacity);
                LaserBands.pop_back();
            }
            else if(Level.first > 0.1)
            {
                LaserBands.back() -> SetFillColorAlpha(31,Opacity);
                
                // Print values for laser bands
                if(iterator)
                {
                    std::cout << "Excitation band " << iterator << " with J = " << Level.first << ": " << Level.second+DeltaE-BWidth << " eV to " << Level.second+DeltaE+BWidth << " eV" << std::endl;
                }
            }
            else
            {
                LaserBands.back() -> SetFillColorAlpha(38,Opacity);
                
                // Print values for laser bands
                if(iterator)
                {
                    std::cout << "Excitation band " << iterator << " with J = " << Level.first << ": " << Level.second+DeltaE-BWidth << " eV to " << Level.second+DeltaE+BWidth << " eV" << std::endl;
                }
            }
        }
        iterator++;
    }
    
    // Pop back ionisation state band since it is a threshold
    LiquidBands.pop_back();
    LaserBands.pop_back();
    
    // Colour first and last line black
    GasLines.front() -> SetLineColor(kBlack);
    GasLines.back() -> SetLineColor(kBlack);    
    LiquidLines.front() -> SetLineColor(kBlack);
    LiquidLines.back() -> SetLineColor(kBlack);
    LiquidBands.front() -> SetFillColorAlpha(kBlack,Opacity);
    LaserBands.front() -> SetFillColorAlpha(kBlack,Opacity);
    
    // Create laser (virtual) energy levels and arrows
    std::vector<TLine*> LaserLine;
    std::vector<TArrow*> LaserArrows;
    
    // Fill lines and arrows
    for(unsigned int i = 0; i < 3; i++)
    {
        // Create three horizontal lines for multiphoton ionisation levels
        LaserLine.push_back(new TLine(0.05,LaserPhotonEnergy*(i+1),1,LaserPhotonEnergy*(i+1)));
        LaserLine.back() -> SetLineColor(kViolet-6);
        LaserLine.back() -> SetLineStyle(2);
        
        // Create arrows between the levels
        LaserArrows.push_back(new TArrow((0.05+1)/4*(i+1),LaserPhotonEnergy*i,(0.05+1)/4*(i+1),LaserPhotonEnergy*(i+1)-0.07,0.02,"|>"));
        LaserArrows.back() -> SetAngle(30);
        LaserArrows.back() -> SetLineWidth(2);
        LaserArrows.back() -> SetFillColor(kViolet-6);
        LaserArrows.back() -> SetLineColor(kViolet-6);
    }
    
    // Create an empty histogram to generate a default TPad with axis labels etc.
    TH1D* EmptyHistogramGas = new TH1D("h0","Gaseous Argon Energy Levels",10,x_min,x_max);
    EmptyHistogramGas -> GetYaxis() -> SetRangeUser(y_min,y_max);
    EmptyHistogramGas -> GetXaxis() -> SetLabelOffset(999);
    EmptyHistogramGas -> GetXaxis() -> SetLabelSize(0);
    EmptyHistogramGas -> GetXaxis() -> SetTickLength(0);
    EmptyHistogramGas -> GetXaxis() -> SetAxisColor(kWhite);
    EmptyHistogramGas -> GetYaxis() -> SetTitle("Energy [eV]");
    EmptyHistogramGas -> GetYaxis() -> SetTitleOffset(0.95);
    
    TH1D* EmptyHistogramLiquid = new TH1D("h1","Liquid Argon Energy Bands",10,x_min,x_max);
    EmptyHistogramLiquid -> GetYaxis() -> SetRangeUser(y_min,y_max);
    EmptyHistogramLiquid -> GetXaxis() -> SetLabelOffset(999);
    EmptyHistogramLiquid -> GetXaxis() -> SetLabelSize(0);
    EmptyHistogramLiquid -> GetXaxis() -> SetTickLength(0);
    EmptyHistogramLiquid -> GetXaxis() -> SetAxisColor(kWhite);
//     EmptyHistogramLiquid -> GetYaxis() -> SetTitle("Energy [eV]");
    
    TH1D* EmptyHistogramLaser = new TH1D("h2","Three-photon Ionisation",10,x_min,x_max);
    EmptyHistogramLaser -> GetYaxis() -> SetRangeUser(y_min,y_max);
    EmptyHistogramLaser -> GetXaxis() -> SetLabelOffset(999);
    EmptyHistogramLaser -> GetXaxis() -> SetLabelSize(0);
    EmptyHistogramLaser -> GetXaxis() -> SetTickLength(0);
    EmptyHistogramLaser -> GetXaxis() -> SetAxisColor(kWhite);
    EmptyHistogramLaser -> GetYaxis() -> SetTitle("Energy [eV]");
    
    // Change Style
    gStyle->SetPadTickY(0);
    gStyle->SetPadTickX(0);
    gStyle->SetPadBottomMargin(0.02);
//     gStyle->SetPadLeftMargin(0.11);
    gStyle->SetPadRightMargin(0.2);
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","",0,0,900,1600);
    C0 -> SetFrameLineColor(kWhite);    
    EmptyHistogramGas -> Draw();
    // Draw all lines
    for(auto Line : GasLines)
    {
        Line -> Draw("SAME");
    }
    // Create Labels
    TLatex LabelGas;
    LabelGas.SetTextAlign(12); // Centred 
    LabelGas.DrawLatex(1,0,"ground state");
    LabelGas.DrawLatex(1,EnergyLevels.back().second, "ionised state");
//     Legend -> Draw("SAME");
    C0 -> SaveAs("../images/Detector/ArEnergyLevelsGas.pdf");
    
    //-------------------------------------------------------------------------
    TCanvas* C1 = new TCanvas("C1","",0,0,900,1600);
    C1 -> SetFrameLineColor(kWhite);    
    EmptyHistogramLiquid -> Draw();
//     for(auto Line : LiquidLines)
//     {
//         Line -> Draw("SAME");
//     }
    LiquidLines.back() -> Draw("SAME");
    for(auto Band : LiquidBands)
    {
        Band -> Draw("SAME");
    }
    TLatex LabelLiquid;
    LabelLiquid.SetTextAlign(12); // Centred 
    LabelLiquid.DrawLatex(1,-0.1,"ground state");
    LabelLiquid.DrawLatex(1,I_LAr, "ionised state");
    
    C1 -> SaveAs("../images/Detector/ArEnergyLevelsLiquid.pdf");
    
    //-------------------------------------------------------------------------
    TCanvas* C2 = new TCanvas("C2","",0,0,1000,1500);
    C2 -> SetFrameLineColor(kWhite);    
    EmptyHistogramLaser -> Draw();
    LiquidLines.back() -> Draw("SAME");
    for(auto Band : LaserBands)
    {
        Band -> Draw("SAME");
    }
    for(auto Line : LaserLine)
    {
        Line -> Draw("SAME");
    }
    for(auto Arrow : LaserArrows)
    {
        Arrow -> Draw();
    }
    TLatex LabelLaser;
    LabelLaser.SetTextAlign(12); // Centred 
    LabelLaser.DrawLatex(1,-0.1,"ground state");
    LabelLaser.DrawLatex(1,I_LAr, "ionised state");
    LabelLaser.DrawLatex(1.02,LaserPhotonEnergy,"#color[873]{4.66 eV}");
    LabelLaser.DrawLatex(1.02,LaserPhotonEnergy*2,"#color[873]{9.32 eV}");
    LabelLaser.DrawLatex(1.02,LaserPhotonEnergy*3,"#color[873]{13.98 eV}");
    
    C2 -> SaveAs("../images/Detector/ArEnergyLevelsLaser.pdf");
}
