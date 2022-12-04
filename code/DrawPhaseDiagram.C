// Phase diagram drawing program for argon
// All constants and formulas formulas are found in https://doi.org/10.1063/1.556037 or https://lar.bnl.gov/properties/basic.html

// C++ headers
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const double pi = 3.14159265358979323846;
const double T_c = 150.687; // critical temperature in K
const double p_c = 4.863*10; // critical pressure in bar (*10 conversion from MPa)
const double T_t = 83.8058; // triple point temperature K
const double p_t = 0.0688909*10; // triple point pressure bar
const double p_ref = 0.1*10; // atmospheric pressure in bar
const double T_b = 87.28; // Boiling temperature at 1 bar
const double T_m = 83.8; // Melting temperature at 1 bar
const double rho_c = 0.5356; // g/cm^3

// Set these variables to your liking
const double T_min = 60; // K
const double T_max = 180; // K
const double p_min = 1e-3*10; // bar
const double p_max = 1e3*10; // bar
const unsigned long int SampleNo = 200;

double BoilingLine(double T)
{   
    double phi = (1-T/T_c);
    
    const double a1 = -5.9409785;
    const double a2 = 1.3553888;
    const double a3 = -0.46497607;
    const double a4 = -1.5399043;
    
    // Function Values
    double FunctionValue = p_c*std::exp( T_c/T*(a1*phi + a2*std::pow(phi,1.5) + a3*std::pow(phi,2) + a4*std::pow(phi,4.5)) );
        
    return FunctionValue;
}

double MeltingLine(double T)
{
    const double a1 = -7476.2665;
    const double a2 = 9959.0613;
    
    double FunctionValue = p_t*( 1 + a1*(std::pow(T/T_t,1.05)-1) + a2*(std::pow(T/T_t,1.275)-1) );
    
    return FunctionValue;
}

double SublimationLine(double T)
{
    const double a1 = -11.391604;
    const double a2 = -0.39513431;
    
    double FunctionValue = p_t*std::exp( T_t/T*(a1*(1-T/T_t) + a2*std::pow(1-T/T_t,2.7)) );
    
    return FunctionValue;
}

double DensityCalculator(double T)
{
    const double a1 = 1.5004262; 
    const double a2 = -0.31381290; 
    const double a3 = 0.086461622; 
    const double a4 = -0.041477525;
    
    double theta = 1 - T/T_c; 
    
    double Density = rho_c * std::exp( a1*std::pow(theta,1./3.) + a2*std::pow(theta,2./3.) + a3*std::pow(theta,7./3.) + a4*std::pow(theta,4) );
    
    return Density;
}

double MinMaxFunction(double Temp)
{
    double x_min = T_min;
    double x_max = T_max;
    
    double y_min = SublimationLine(x_min);
    double y_max = MeltingLine(x_max);
    
    double FunctionValue = y_min + (y_max-y_min)/(x_max-x_min)*(Temp-x_min);
    
    return FunctionValue;
}

void DrawPhaseDiagram() 
{
    // Print density at standard boiling temperature
    std::cout << "Density at standard boiling temperature: " <<  DensityCalculator(T_b) << " g/cm^3" << std::endl; 
    
    // Create point markers
    TMarker* TriplePoint = new TMarker(T_t, p_t, 8);
    TMarker* CriticalPoint = new TMarker(T_c, p_c, 8);
    TMarker* BoilingPoint = new TMarker(T_b,p_ref,8);
    TMarker* MeltingPoint = new TMarker(T_m,p_ref,8);
    
    TriplePoint -> SetMarkerSize(2);
    CriticalPoint -> SetMarkerSize(2);
    BoilingPoint -> SetMarkerSize(2);
    MeltingPoint -> SetMarkerSize(2);
    
    // Create 1 atmosphere line at 0.1325
    TLine* AtmosphereLine = new TLine(T_min,p_ref,T_max,p_ref);
    
    AtmosphereLine -> SetLineWidth(2);
    AtmosphereLine -> SetLineColor(1);
    AtmosphereLine -> SetLineStyle(2);
    
    // Create text describing stuff
    TText* TriplePointText = new TText(T_t-1,p_t/2, "TP");
    TText* CriticalPointText = new TText(T_c+1,p_c, "CP");
    TText* BoilingPointText = new TText(T_b+1,p_ref/1.7, "BP");
    TText* MeltingPointText = new TText(T_m-6.5,p_ref*1.2, "MP");
    
    TText* AtmosphereText = new TText(T_max-20,p_ref*1.3, "p = 1 bar");

    // Use bold font and large size
    TText* SolidText = new TText((T_t+T_min)/2.-7.5,p_ref*10, "Solid");
    SolidText -> SetTextFont(22);
    SolidText -> SetTextSize(0.08);
    SolidText -> SetTextColor(38);
    
    TText* LiquidText = new TText(T_t+20,p_c*2, "Liquid");
    LiquidText -> SetTextFont(22);
    LiquidText -> SetTextSize(0.08);
    LiquidText -> SetTextColor(40);
    
    TText* GasText = new TText(T_t+30,p_ref/9, "Gas");
    GasText -> SetTextFont(22);
    GasText -> SetTextSize(0.08);
    GasText -> SetTextColor(46);
    
    
    // Create Function Arrays
    double T_Boil[SampleNo], T_Melt[SampleNo], T_Subl[SampleNo]; // X-Values
    double p_Boil[SampleNo], p_Melt[SampleNo], p_Subl[SampleNo]; // Y-Values
    
    
    // Initialise X-Axis tics
    double BoilTic = (T_c -T_t)/(double)(SampleNo-1);
    double MeltTic = (T_max -T_t)/(double)(SampleNo-1);
    double SublTic = (T_t - T_min)/(double)(SampleNo-1);
    
    // Fill Function Arrays
    for(unsigned long int i = 0; i < SampleNo; i++)
    {
        T_Boil[i] = T_t + (double)i*BoilTic;
        T_Melt[i] = T_t + (double)i*MeltTic;
        T_Subl[i] = T_min + (double)i*SublTic;
        
        p_Boil[i] = BoilingLine(T_Boil[i]);
        p_Melt[i] = MeltingLine(T_Melt[i]);
        p_Subl[i] = SublimationLine(T_Subl[i]);
    }
    
    // Create graphs using filled data arrays
    TGraph* BoilingGraph = new TGraph(SampleNo,T_Boil,p_Boil);
    TGraph* MeltingGraph = new TGraph(SampleNo,T_Melt,p_Melt);
    TGraph* SublimationGraph = new TGraph(SampleNo,T_Subl,p_Subl);
    
    // Change line colour
    BoilingGraph -> SetLineColor(46);
    MeltingGraph -> SetLineColor(40);
    SublimationGraph -> SetLineColor(38);
    
    // Change line width
    BoilingGraph -> SetLineWidth(4);
    MeltingGraph -> SetLineWidth(4);
    SublimationGraph -> SetLineWidth(4);
    
    
//     double x, y;
//     SublimationGraph -> GetPoint(0,x,y);
//     std::cout << x << " " << y << std::endl;
    
    // Create MultiGraph and add all graphs to it 
    TMultiGraph* MultiGraph = new TMultiGraph();
    MultiGraph -> Add(BoilingGraph);
    MultiGraph -> Add(MeltingGraph);
    MultiGraph -> Add(SublimationGraph);
    
    // Create canvas and draw functions
    TCanvas* C0 = new TCanvas("C0","Phase Diagram of Argon",0,0,1400,1000);
    C0 -> cd();
    C0 -> SetLogy();
    
    MultiGraph -> SetTitle("Phase Diagram of Argon");
    MultiGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(T_min,T_max);
    MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(p_min,p_max);
//     MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(MinMaxFunction(T_min)/1.5,MinMaxFunction(T_max)*1.5);
    MultiGraph -> GetXaxis() -> SetTitle("Temperature T [K]");
    MultiGraph -> GetYaxis() -> SetTitle("Pressure p [bar]");
    
    MultiGraph -> Draw("AC");
    AtmosphereLine -> Draw("SAME");
    AtmosphereText -> Draw("SAME");
    TriplePoint -> Draw("SAME");
    TriplePointText -> Draw("SAME");
    CriticalPoint -> Draw("SAME");
    CriticalPointText -> Draw("SAME");
    BoilingPoint -> Draw("SAME");
    BoilingPointText -> Draw("SAME");
    MeltingPoint -> Draw("SAME");
    MeltingPointText -> Draw("SAME");
    
    SolidText -> Draw("SAME");
    LiquidText -> Draw("SAME");
    GasText -> Draw("SAME");
    
    C0 -> SaveAs("../images/Detector/PhaseDiagramArgon.pdf");
}
