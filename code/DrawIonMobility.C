// Ion mobility in liquid argon

// C++ headers
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const double pi = 3.14159265358979323846;
const double m_e = 9.1093837015e-31; // kg 
const double m_e_MeV = 0.510998918; // MeV
const double alpha = 1/137.03599911;
const double N_A = 6.0221415e23; // mol^-1
const double m_mu = 105.6583755; // MeV
const double c = 299792458; // m/s
const double h_bar = 1.054571596e-34; // J s
const double e = 1.602176634e-19; // C
const double epsilon_0 = 8.85418782e-12;
const double rho = 1.396; // g/cm^3
const double X_0 = 19.55; // g/cm^2

const double Z = 18; // of argon
const double A =  39.948; // g/mol of argon
const double I = 188.0e-6; // MeV of argon

const unsigned long int SampleNo = 84;
const double x_min = 0;
const double x_max = 4100;
const double y_min = 0;
const double y_max = 1.5e-3;

const double y_min_v = 0;
const double y_max_v = 5;

double IonMobilityFunction(double E)
{
    // Set step function values (error = 5%)
    double mu_1 = 6.00e-4;
    double mu_2 = 9.75e-4;
    double mu_3 = 8.50e-4;
    double mu_4 = 7.75e-4;
    double mu_5 = 7.25e-4;
    
    // Set Step function thresholds
    double ThresholdStart1 = 0.00;
    double ThresholdEnd1 = 238.42;
    double ThresholdStart2 = 502.07;
    double ThresholdEnd2 = 1012.94;
    double ThresholdStart3 = 1237.91;
    double ThresholdEnd3 = 2232.82;
    double ThresholdStart4 = 2461.22;
    double ThresholdEnd4 = 3445.87;
    double ThresholdStart5 = 3695.37;
    double ThresholdEnd5 = 4103.05;
    
    double mu = 0;
    
    // Step function of mobility
    if(E <= ThresholdEnd1) mu = mu_1;
    else if(E < ThresholdStart2) mu = mu_1 + (E-ThresholdEnd1)/(ThresholdStart2-ThresholdEnd1)*(mu_2-mu_1);
    else if(E < ThresholdEnd2) mu = mu_2;
    else if(E < ThresholdStart3) mu = mu_2 + (E-ThresholdEnd2)/(ThresholdStart3-ThresholdEnd2)*(mu_3-mu_2);
    else if(E < ThresholdEnd3) mu = mu_3;
    else if(E < ThresholdStart4) mu = mu_3 + (E-ThresholdEnd3)/(ThresholdStart4-ThresholdEnd3)*(mu_4-mu_3);
    else if(E < ThresholdEnd4) mu = mu_4;
    else if(E < ThresholdStart5) mu = mu_4 + (E-ThresholdEnd4)/(ThresholdStart5-ThresholdEnd4)*(mu_5-mu_4);
    else if(E < ThresholdEnd5) mu = mu_5;
    
    return mu;
}

double IonDriftVelocityFunction(double E, double mu)
{
    double v_Drift = mu*E;
    
    return v_Drift;
}

void DrawIonMobility() 
{
    // Create vector containing relevant particle information: name, mass, charge, critical energy (e and mu only, pi estimated, rest irrelevant), and colour in plot
    std::vector<std::tuple<std::string,double, double, double,unsigned int>> Particle;
    
    
    // Create vector containing different functions for different particles
    std::vector<TGraph*> IonMobilityGraph;
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.15,0.78,0.35,0.85);
    Legend->SetNColumns(1);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.2);
//     Legend->SetHeader("Particles");
    
    // Prepare Graph fill variables
    double XValue[SampleNo], YValue[SampleNo], YYValue[SampleNo];
    
    // Create tick length
    double XAxisTic = (x_max - x_min)/(double)(SampleNo-1);
    
    // Create Log tick length
//     double XAxisTic = (std::log10(x_max) - std::log10(x_min))/(double)(SampleNo-1);
    
    // Debug
//     std::cout << XAxisTic << std::endl;
    std::cout << IonMobilityFunction(250) << std::endl;
    
    double temp = 0;
    // Loop over all Graph samples
    for(unsigned long int i = 0; i < SampleNo; i++)
    {
        // Calculate x-values
        XValue[i] = x_min + (double)i*XAxisTic;
        // Calculate log spaced x-values
//         XValue[i] = std::pow(10, std::log10(x_min) + i*XAxisTic);
            
        // Use Function to calculate y-values
        YValue[i] = IonMobilityFunction(XValue[i]);
        YYValue[i] = IonDriftVelocityFunction(XValue[i],YValue[i]);
        
        if(temp <= YYValue[i]) temp = YYValue[i];
        else std::cout << i << " " << temp << " " << YYValue[i] << std::endl;
    }
    
    // Loop over different particles
//     for(unsigned long int i = 0; i < MiniSampleNo; i++)
//     {
        // Calculate x-values
//         MiniXValue[i] = mini_x_min + (double)i*MiniXAxisTic;
        // Calculate log spaced x-values
//         XValue[i] = std::pow(10, std::log10(x_min) + i*XAxisTic);
            
        // Use Function to calculate y-values
//         MiniYValue[i] = IonMobilityFunction(MiniXValue[i]);
//     }
    
    // Fill and configure graphs
    IonMobilityGraph.push_back( new TGraph(SampleNo,XValue,YValue) );
    IonMobilityGraph.back()->SetLineColor(38);
    IonMobilityGraph.back()->SetLineWidth(4);
    IonMobilityGraph.back() -> SetTitle("Argon Ion Mobility");
    IonMobilityGraph.back() -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    IonMobilityGraph.back() -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    IonMobilityGraph.back() -> GetXaxis() -> SetTitle("Electric Field E [V cm^{-1}]");
    IonMobilityGraph.back() -> GetYaxis() -> SetTitle("Argon Ion Mobility  #mu_{i} [cm^{2} V^{-1} s^{-1}]");
    IonMobilityGraph.back() -> GetYaxis() -> SetMaxDigits(2);
    
    IonMobilityGraph.push_back( new TGraph(SampleNo,XValue,YYValue) );
    IonMobilityGraph.back()->SetLineColor(38);
    IonMobilityGraph.back()->SetLineWidth(4);
    IonMobilityGraph.back() -> SetTitle("Argon Ion Drift Velocity");
    IonMobilityGraph.back() -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    IonMobilityGraph.back() -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min_v,y_max_v);
    IonMobilityGraph.back() -> GetXaxis() -> SetTitle("Electric Field E [V cm^{-1}]");
    IonMobilityGraph.back() -> GetYaxis() -> SetTitle("Argon Ion Velocity  v_{i} [cm s^{-1}]");
    IonMobilityGraph.back() -> GetYaxis() -> SetMaxDigits(2);
    
    // Add legend entry
    Legend -> AddEntry( IonMobilityGraph.at(0), "T = 87.303 K","" ); // "L" for line
    Legend -> AddEntry( IonMobilityGraph.at(0), "p = 1 bar","" ); // "L" for line
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Ion Mobility",0,0,1400,1000);
    
    IonMobilityGraph.at(0) -> Draw("AC");
    Legend -> Draw("SAME");
    
    C0 -> SaveAs("../images/Detector/IonMobility.pdf");
    
    TCanvas* C1 = new TCanvas("C1","Drift Velocity",0,0,1400,1000);
    
    IonMobilityGraph.at(1) -> Draw("AC");
    Legend -> Draw("SAME");
    
    C1 -> SaveAs("../images/Detector/IonVelocity.pdf");
}
