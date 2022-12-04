// Space charge model in liquid argon
// WARNING: DOES NOT WORK FOR LArTPCs!!!
// They are too big 

// C++ headers
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const unsigned long int SampleNo = 10000;
const double x_min = 0;
const double x_max = 2.56; // m
// const double y_min = 0;
// const double y_max = 1;

const double epsilon_0 = 8.85418782e-12;
const double epsilon = 1.504;
const double mu_i = 1.5e-7; // m^2/(Vs) (2e-4 cm^2/(Vs))
const double J = 1.6e-10; // C/(m^3 s)
const double E_0 = 27340; // V/m

// Palestini model
double alpha = x_max/E_0*std::sqrt(J/epsilon_0/epsilon/mu_i);
double E_a = E_0*(1 - std::pow(alpha,2)/6 - std::pow(alpha,4)/180);

// UBooNE approach
const double eta = J/E_0/mu_i; // C/m^4
double E_aUB = E_0 - eta/(6*epsilon*epsilon_0)*std::pow(x_max,2);

double ElectricField(double x)
{
    double EField = E_0 * std::sqrt( std::pow(E_a/E_0,2) + std::pow(alpha*x/x_max,2) );
    return EField;
}

double ChargeDensity(double x)
{   
    double CDensity = std::pow(alpha*E_0/x_max,2)*epsilon*epsilon_0*x / ElectricField(x);  //J*x / ElectricField(x) / mu_i;
    return CDensity;
}

double ElectricFieldUB(double x)
{
    double EField = E_aUB + eta/(2*epsilon*epsilon_0)*std::pow(x,2);
    return EField;
}

double ChargeDensityUB(double x)
{
    double CDensity = eta*x;
    return CDensity;
}

void DrawSpaceCharge()
{
    std::cout << "alpha = " << alpha << " if larger than 2, model invalid" << std::endl;
    std::cout << "Paper:  E_a = " << E_a/100 << " V/cm" << std::endl;
    std::cout << "eta = " << eta << " C/m^4" << std::endl;
    std::cout << "uBooNE: E_a = " << E_aUB/100 << " V/cm" << std::endl;
    
    // Create vector containing different functions for different particles
    TGraph* EFieldGraphPaper;
    TGraph* EFieldGraphConst;
    TGraph* EFieldGraphUB;
    
    TGraph* CDensityGraphPaper;
    TGraph* CDensityGraphUB;
    
    // Create tick length
    double XAxisTic = (x_max - x_min)/(double)(SampleNo-1);
    
    // Create temporary arrays
    double XValue[SampleNo];
    double YValueEField[SampleNo];
    double YValueEFieldConst[SampleNo];
    double YValueEFieldUB[SampleNo];
    double YValueCDensity[SampleNo];
    double YValueCDensityUB[SampleNo];
    
    // Sample fill loop
    for(unsigned long int i = 0; i < SampleNo; i++)
    {
        // Calculate x-values
        XValue[i] = x_min + (double)i*XAxisTic;
        
        
        // Use Function to calculate y-values
        YValueEField[i] = ElectricField(XValue[i])/100; // V/cm
        YValueEFieldConst[i] = E_0/100; // V/cm
        YValueEFieldUB[i] = ElectricFieldUB(XValue[i])/100; // V/cm
        YValueCDensity[i] = ChargeDensity(XValue[i])*1e9; // nC/m^3
        YValueCDensityUB[i] = ChargeDensityUB(XValue[i])*1e9; // nC/m^3

        // Unit conversion m -> cm
        XValue[i] *= 100;
    }
    
    // Fill graphs
    EFieldGraphPaper = new TGraph(SampleNo,XValue,YValueEField);
    EFieldGraphPaper -> SetLineColor(38);
    EFieldGraphConst = new TGraph(SampleNo,XValue,YValueEFieldConst);
    EFieldGraphConst -> SetLineColor(46);
    EFieldGraphUB = new TGraph(SampleNo,XValue,YValueEFieldUB);
    EFieldGraphUB -> SetLineColor(30);
    
    CDensityGraphPaper = new TGraph(SampleNo,XValue,YValueCDensity);
    CDensityGraphPaper -> SetLineColor(38);
    CDensityGraphUB = new TGraph(SampleNo,XValue,YValueCDensityUB);
    CDensityGraphUB -> SetLineColor(30);
    
    // Create and fill multi graph
    TMultiGraph* EFieldGraph = new TMultiGraph();
    EFieldGraph -> Add(EFieldGraphConst);
    EFieldGraph -> Add(EFieldGraphUB);
    EFieldGraph -> Add(EFieldGraphPaper);
    
    TMultiGraph* CDensityGraph = new TMultiGraph();
    CDensityGraph -> Add(CDensityGraphUB);
    CDensityGraph -> Add(CDensityGraphPaper);
    
    // Setup Legend
    TLegend* FieldLegend = new TLegend(0.15,0.70,0.50,0.85);
    FieldLegend->SetNColumns( 1 );
    FieldLegend->SetLineStyle ( 0 );
    FieldLegend->SetLineColorAlpha ( 0,0 );
    FieldLegend->SetFillStyle ( 0 );
    FieldLegend->SetMargin ( 0.2 );
    FieldLegend->SetEntrySeparation(0.5);
    FieldLegend->SetTextSize(0.03);
    
    // Add Legend Entries
    FieldLegend -> AddEntry(EFieldGraphConst," Constant field E_{0} (#rho_{i} = 0)", "L");
    FieldLegend -> AddEntry(EFieldGraphPaper," Palestini model", "L");
    FieldLegend -> AddEntry(EFieldGraphUB," MicrBooNE approximation", "L");
    
    // Setup Legend
    TLegend* DensityLegend = new TLegend(0.15,0.75,0.50,0.85);
    DensityLegend->SetNColumns( 1 );
    DensityLegend->SetLineStyle ( 0 );
    DensityLegend->SetLineColorAlpha ( 0,0 );
    DensityLegend->SetFillStyle ( 0 );
    DensityLegend->SetMargin ( 0.2 );
    DensityLegend->SetEntrySeparation(0.5);
    DensityLegend->SetTextSize(0.03);
    
    // Add Legend Entries
    DensityLegend -> AddEntry(CDensityGraphPaper," Palestini model", "L");
    DensityLegend -> AddEntry(CDensityGraphUB," MicrBooNE approximation", "L");
    
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","SpaceChargeField",0,0,1400,1000);
    
    EFieldGraph -> SetTitle("Electric Field Distortion");
    EFieldGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min*100,x_max*100);
    EFieldGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(200,400);
    EFieldGraph -> GetYaxis() -> SetMaxDigits(3);
    EFieldGraph -> GetXaxis() -> SetTitle("Drift Axis x [cm]");
    EFieldGraph -> GetYaxis() -> SetTitle("Drift Electric Field E_{x}(x) [V cm^{-1}]");

    EFieldGraph -> Draw("AC");
    FieldLegend -> Draw();
    
    C0 -> SaveAs("../images/Detector/SpaceChargeField.pdf");
    
    // Create canvas, set it up, and draw functions
    TCanvas* C1 = new TCanvas("C1","SpaceChargeDensity",0,0,1400,1000);
    
    CDensityGraph -> SetTitle("Space Charge Density");
    CDensityGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min*100,x_max*100);
    CDensityGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(0,100);
    CDensityGraph -> GetXaxis() -> SetTitle("Drift Axis x [cm]");
    CDensityGraph -> GetYaxis() -> SetTitle("Space Charge Density #rho_{i}(x) [nC m^{-3}]");
    CDensityGraph -> GetYaxis() -> SetTitleOffset(0.9);
    
    CDensityGraph -> Draw("AC");
    DensityLegend -> Draw();
    
    C1 -> SaveAs("../images/Detector/SpaceChargeDensity.pdf");
    
}
