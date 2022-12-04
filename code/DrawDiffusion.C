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

const double Z = 18; // of argon
const double A =  39.948; // g/mol of argon
const double I = 188.0e-6; // MeV of argon

const double D_L = 3.6; // cm^2 s^-1 former 5.4 microboone document factor 2/3 missing
const double D_T = 8.9; // cm^2 S^-1 former 13.9 microboone document factor 2/3 missing
const double epsilon_L = 0.013; // eV
const double epsilon_T = 0.026; // eV

const double E_uB = 273.4; // V cm^-1
const double mu = 516; // cm^2 V^-1
const double x_drift = 256; // cm

const unsigned long int SampleNo = 10000;
// const double x_min = 254;
// const double x_max = 258;
const double x_min = -5;
const double x_max = 5;
const double y_min = -2;
const double y_max = 25;

const double erg = 1e7;

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

double DiffusionFunction(double x, double y, double z, double t)
{
//     double N = 1/(4*pi*D_T*t*std::sqrt(4*pi*D_L*t) ) * std::exp( -(std::pow(y,2)+std::pow(z,2))/(4*D_T*t) ) * std::exp( -std::pow(x-DriftFunction(E_uB)*t,2)/(4*D_L*t) );
    // Distribution not including drift
//     double N = 1/(4*pi*D_T*t*std::sqrt(4*pi*D_L*t) ) * std::exp( -(std::pow(y,2)+std::pow(z,2))/(4*D_T*t) ) * std::exp( -std::pow(x,2)/(4*D_L*t) );
    
    // Distribution only in the drift-coordinate without drift shift for overlay purposes
    double N = 1/( std::sqrt(4*pi*D_L*t) ) * std::exp( -std::pow(x,2)/(4*D_L*t) );
    return N;
}

void DrawDiffusion() 
{
    // Create vector containing relevant particle information: name, mass, charge, critical energy (e and mu only, pi estimated, rest irrelevant), and colour in plot
    std::vector<std::tuple<std::string,double, double, double,unsigned int>> Particle;
    
    std::cout << "Calculated D_L = " << 2./3.*epsilon_L*DriftFunction(E_uB)/E_uB << " Mobility mu = " << DriftFunction(E_uB)/E_uB << std::endl;
    std::cout << "Calculated D_T = " << 2./3.*epsilon_T*mu << " Mobility mu_0 = " << mu << std::endl;
    
    // Create vector containing different functions for different particles
    std::vector<TGraph*> DiffusionGraphs;
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.15,0.70,0.45,0.85);
    Legend->SetNColumns(1);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.5);
//     Legend->SetHeader("Particles");
    
    // Prepare Graph fill variables
    double XValue[SampleNo], YValue[SampleNo], YYValue[SampleNo];
    
    std::vector< std::tuple<double,unsigned int,std::string> > TimeContainer;
    std::vector<double> t = {5e-5,5e-4,x_drift/DriftFunction(E_uB)};
    
    std::cout << x_drift/DriftFunction(E_uB) << std::endl;
    
    // Create tick length
    double XAxisTic = (x_max - x_min)/(double)(SampleNo-1);
    
    // Create Log tick length
//     double XAxisTic = (std::log10(x_max) - std::log10(x_min))/(double)(SampleNo-1);
    
//     t.push_back(std::make_tuple(1e-6,39));
    TimeContainer.push_back(std::make_tuple(t.at(0),46,"after 0.05 ms"));
    TimeContainer.push_back(std::make_tuple(t.at(1),39,"after 0.50 ms"));
    TimeContainer.push_back(std::make_tuple(t.at(2),38,"after 2.25 ms"));
    
    // Loop over all times
    for(auto time : TimeContainer)
    {
        std::cout << "Time: " << std::get<0>(time) << "  peak value: " << DiffusionFunction(0,0,0,std::get<0>(time)) << " sigma_L " << std::sqrt(2*D_L*std::get<0>(time)) << " sigma_T " << std::sqrt(2*D_T*std::get<0>(time)) <<std::endl;
        // Loop over all Graph samples
        for(unsigned long int i = 0; i < SampleNo; i++)
        {
            // Calculate x-values
            XValue[i] = x_min + (double)i*XAxisTic;
            
            // Use Function to calculate y-values
            YValue[i] = DiffusionFunction(XValue[i]/10,0,0,std::get<0>(time));
        }
        
        // Fill Graphs
        DiffusionGraphs.push_back( new TGraph(SampleNo,XValue,YValue) );
        DiffusionGraphs.back()->SetLineColor(std::get<1>(time));
        DiffusionGraphs.back()->SetLineWidth(4);
        
        Legend -> AddEntry( DiffusionGraphs.back(),(" "+std::get<2>(time)).c_str(),"L" );
    }
    
    // Create and fill multi graph
    TMultiGraph* MultiGraph = new TMultiGraph();
    
    // Loop over all graphs
    for(auto Graph : DiffusionGraphs)
    {
        // Sanity check if integrals are the same
        std::cout << "Integral: " << Graph -> Integral() << std::endl;
        
        // Fill the multi graph
        MultiGraph->Add(Graph);
    }
    
    MultiGraph -> SetTitle("Longitudinal Electron Diffusion in LAr");
    MultiGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    MultiGraph -> GetXaxis() -> SetTitle("Drift Coordinate x = x' + v_{d}t [mm]");
    MultiGraph -> GetYaxis() -> SetTitle("Longitudinal e^{-} Distribution n(t,z)/N_{0} [cm^{-1}]");
    MultiGraph -> GetYaxis() -> SetMaxDigits(2);
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Electron Diffusion",0,0,1400,1000);
    
    MultiGraph -> Draw("AC");
    Legend -> Draw("SAME");
    
    C0 -> SaveAs("../images/Detector/Diffusion.pdf");
    
//     TCanvas* C1 = new TCanvas("C1","Drift Velocity",0,0,1400,1000);
    
//     DiffusionGraphs.at(1) -> Draw("AC");
//     Legend -> Draw("SAME");
    
//     C1 -> SaveAs("../images/Detector/Diffusion2.pdf");
}
