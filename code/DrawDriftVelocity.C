// Drift velocity of electrons in liquid argon.

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

const unsigned long int SampleNo = 101;
const double x_min = 0;
const double x_max = 1e4;
const double y_min = 0;
const double y_max = 5e5;

const unsigned long int MiniSampleNo = 11;
const double mini_x_min = 0;
const double mini_x_max = 1e3;
const double mini_y_min = 0;
const double mini_y_max = 2.5e5;

const double c_s = 84290; // cm/s

// Atlas fit
const double p1 = -0.01481;
const double p2 = -0.0075;
const double p3 = 0.141;
const double p4 = 12.4;
const double p5 = 1.627;
const double p6 = 0.317;
const double T_0 = 90.371; // °K

// Icarus Fit
const double p1l = -0.0462553;
const double p2l = 0.0148508;
const double p3l = 1.64156;
const double p4l = 1.273;
const double p5l = 0.0086608;
const double p6l = 4.71489;
const double T_0l = 104.326; // °K
// const double T_0l = 104.0; // °K

// Mobility
const double mu = 5.16;

const double mu_err = 0.27;

const double p1_err = 0.00095;
const double p2_err = 0.0028;
const double p3_err = 0.023;
const double p4_err = 2.7;
const double p5_err = 0.078;
const double p6_err = 0.021;

double DriftFunction(double E)
{   
    // Temperature (variable)
    double T = 87.303; // K
//     double T = 89.3; // K
    
    // Change units: mm -> cm, kV -> V, and us-> s
    E /= 1000.;
    
    // Drift velocity
    double v_Drift;
    
    // Drift function
    if(E < 0.2) v_Drift = mu*E;
    else if(E >= 0.2 && E < 0.69) v_Drift = ( 1 + p1l*(T - T_0l) ) * ( p3l*E*std::log(1+p4l/E) + p5l*std::pow(E,p6l) ) + p2l*(T - T_0l);
//     if(E == 0) v_Drift = 0.0;
//     else if(E > 0. && E <= 0.7) v_Drift = ( 1 + p1l*(T - T_0l) ) * ( p3l*E*std::log(1+p4l/E) + p5l*std::pow(E,p6l) ) + p2l*(T - T_0l);
    else v_Drift = ( 1 + p1*(T - T_0) ) * ( p3*E*std::log(1+p4/E) + p5*std::pow(E,p6) ) + p2*(T - T_0);
    
    // Change units: mm -> cm, kV -> V, and us-> s
    v_Drift *= 1e5;
    
    return v_Drift;
}

void DrawDriftVelocity() 
{
    // Create vector containing relevant particle information: name, mass, charge, critical energy (e and mu only, pi estimated, rest irrelevant), and colour in plot
    std::vector<std::tuple<std::string,double, double, double,unsigned int>> Particle;
    
    
    // Create vector containing different functions for different particles
    std::vector<TGraph*> DrifVelocityGraph;
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.15,0.78,0.35,0.85);
    Legend->SetNColumns(2);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.2);
//     Legend->SetHeader("Particles");
    
    // Prepare Graph fill variables
    double XValue[SampleNo], YValue[SampleNo], MiniXValue[MiniSampleNo], MiniYValue[MiniSampleNo];
    
    // Create tick length
    double XAxisTic = (x_max - x_min)/(double)(SampleNo-1);
    double MiniXAxisTic = (mini_x_max - mini_x_min)/(double)(MiniSampleNo-1);
    
    // Create Log tick length
//     double XAxisTic = (std::log10(x_max) - std::log10(x_min))/(double)(SampleNo-1);
    
    // Debug
//     std::cout << XAxisTic << std::endl;
    std::cout << "MicroBooNE drift velocity: " << DriftFunction(273.4) << " cm/s" << std::endl;
    std::cout << "MicroBooNE full drift time: " << 256/DriftFunction(273.4) << " s" << std::endl;
//     std::cout << DriftFunction(0.5) << std::endl;
    
    // Loop over all Graph samples
    for(unsigned long int i = 0; i < SampleNo; i++)
    {
        // Calculate x-values
        XValue[i] = x_min + (double)i*XAxisTic;
        // Calculate log spaced x-values
//         XValue[i] = std::pow(10, std::log10(x_min) + i*XAxisTic);
            
        // Use Function to calculate y-values
        YValue[i] = DriftFunction(XValue[i]);
    }
    
    // Loop over different particles
    for(unsigned long int i = 0; i < MiniSampleNo; i++)
    {
        // Calculate x-values
        MiniXValue[i] = mini_x_min + (double)i*MiniXAxisTic;
        // Calculate log spaced x-values
//         XValue[i] = std::pow(10, std::log10(x_min) + i*XAxisTic);
            
        // Use Function to calculate y-values
        MiniYValue[i] = DriftFunction(MiniXValue[i]);
    }
    
    // Fill and configure graphs
    DrifVelocityGraph.push_back( new TGraph(SampleNo,XValue,YValue) );
    DrifVelocityGraph.back()->SetLineColor(46);
    DrifVelocityGraph.back()->SetLineWidth(4);
    DrifVelocityGraph.back() -> SetTitle("Drift Velocity of Electrons");
    DrifVelocityGraph.back() -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    DrifVelocityGraph.back() -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    DrifVelocityGraph.back() -> GetXaxis() -> SetTitle("Electric Field E [V cm^{-1}]");
    DrifVelocityGraph.back() -> GetYaxis() -> SetTitle("Electron Drift Velocity v_{e} [cm s^{-1}]");
    DrifVelocityGraph.back() -> GetYaxis() -> SetMaxDigits(1);
//     DrifVelocityGraph.back() -> GetYaxis() -> SetDecimals();
    
    DrifVelocityGraph.push_back( new TGraph(MiniSampleNo,MiniXValue,MiniYValue) );
    DrifVelocityGraph.back()->SetLineColor(46);
    DrifVelocityGraph.back()->SetLineWidth(2);
    DrifVelocityGraph.back() -> SetTitle("");
    DrifVelocityGraph.back() -> GetHistogram() -> GetXaxis() -> SetRangeUser(mini_x_min,mini_x_max);
    DrifVelocityGraph.back() -> GetHistogram() -> GetYaxis() -> SetRangeUser(mini_y_min,mini_y_max);
    DrifVelocityGraph.back() -> GetYaxis() -> SetMaxDigits(1);
    DrifVelocityGraph.back() -> GetYaxis() -> SetDecimals();
    
    // Add legend entry
    Legend -> AddEntry( DrifVelocityGraph.at(0), " T = 87.303 K","" ); // "L" for lie
    
    // Create sound velocity line and text
    TLine* ConstMobility = new TLine(x_min,mu*0.0,0.2,mu*0.2);
    ConstMobility -> SetLineWidth(2);
    ConstMobility -> SetLineColor(1);
    ConstMobility -> SetLineStyle(0);
    
    
    // Create sound velocity line and text
    TLine* CsLine = new TLine(x_min,c_s,x_max/2,c_s);
    CsLine -> SetLineWidth(2);
    CsLine -> SetLineColor(1);
    CsLine -> SetLineStyle(0);
    TLatex* CsText = new TLatex(x_min+1000,c_s + 1e4,"v_{e} = c_{s}");
    
    TLine* CsMiniLine = new TLine(mini_x_min,c_s,mini_x_max,c_s);
    CsMiniLine -> SetLineWidth(2);
    CsMiniLine -> SetLineColor(1);
    CsMiniLine -> SetLineStyle(0);
    TLatex* CsMiniText = new TLatex(mini_x_max-200,c_s + 5e3,"v_{e} = c_{s}");
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Drift Velocity",0,0,1400,1000);
//     C0 -> SetLogx();
//     C0 -> SetLogy();
    
    DrifVelocityGraph.at(0) -> Draw("AC");
    CsLine -> Draw("SAME");
    CsText -> Draw("SAME");
//     ConstMobility -> Draw("SAME");
    Legend -> Draw("SAME");
    
    TPad *subpad = new TPad("subpad","",0.35,0.13,0.87,0.60);
    subpad -> SetTopMargin(0.05);
    subpad -> SetBottomMargin(0.05);
    subpad -> SetLeftMargin(0.07);
    subpad -> SetRightMargin(0.05);
    subpad -> Draw();
    subpad -> cd();
    DrifVelocityGraph.at(1) -> Draw("AC");
    CsMiniLine -> Draw("SAME");
    CsMiniText -> Draw("SAME");
    
    C0 -> SaveAs("../images/Detector/DriftVelocityElectron.pdf");
}
