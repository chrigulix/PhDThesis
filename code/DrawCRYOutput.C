// Particle analysis Program

// C++ Headers
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include <typeinfo>

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

// Define globals
void analysis(std::string FileName);
const double PI = 4.0*atan(1.0);
const double e = 1.60217733e-19;
const double epsilon_0 = 8.85418782e-12;
const double epsilon_r = 1.53;

const double SubBoxArea = 1e8; // cm^2

// Input file name and output path
const std::string FileName = "/home/christoph/Programming/CosmicAnalysis/Cosmic/CRY-Output.root";
const std::string StoragePath = "../images/CosmicGammaBackground/";

void DrawCRYOutput()
{

    // Names
    std::string HitParticleName;
    std::string DecayParticleName;
    std::string HitFileName;
    std::string DecayFileName;
    
    // Numbers
    double SinPhi;
    double CosPhi;
    double PhiAngle;
    
    // Binning of TH2
    int ThetaBin = 15;
    int PhiBin = 30;
    int EnergyBin = 30;
    int EnergyMinBin = -3;
    int EnergyMaxBin = 2;
    
    // Bin sizes
    double dPhi = 2*PI/PhiBin;
    double dTheta = PI/2/ThetaBin;
    
    double GeV = 0.001; // MeV -> GeV
    
    // Momentum vectors
    std::vector<double> StartMomentum (3);
    std::vector<double> HitMomentumTPC (3);
    std::vector<double> DecayElectMomentum(3);
    
    std::string TreeName = "MCTree";
    TFile* File = new TFile (FileName.c_str());
    TTree* Tree = (TTree*)File->Get(TreeName.c_str());
    
    
    int EventsElapsed;
    double TimeElapsed;
    
    std::vector<int> EventNumber;
    std::vector<int>*pEventNumber = &EventNumber;
    
    std::vector<int> PrimaryParticleCharge;
    std::vector<int>*pPrimaryParticleCharge = &PrimaryParticleCharge;
    
    std::vector<std::string> PrimaryParticleName;
    std::vector<std::string>*pPrimaryParticleName = &PrimaryParticleName;
    
    std::vector<double> PrimaryParticleTime;
    std::vector<double>*pPrimaryParticleTime = &PrimaryParticleTime;
    
    std::vector<double> PrimaryParticleEnergy;
    std::vector<double>*pPrimaryParticleEnergy = &PrimaryParticleEnergy;
    
    std::vector<double> PrimaryParticleMomentum_x, PrimaryParticleMomentum_y, PrimaryParticleMomentum_z;
    std::vector<double>*pPrimaryParticleMomentum_x = &PrimaryParticleMomentum_x,*pPrimaryParticleMomentum_y = &PrimaryParticleMomentum_y,*pPrimaryParticleMomentum_z = &PrimaryParticleMomentum_z;
    
    std::vector<double> PrimaryParticlePosition_x, PrimaryParticlePosition_y, PrimaryParticlePosition_z;
    std::vector<double>*pPrimaryParticlePosition_x = &PrimaryParticlePosition_x,*pPrimaryParticlePosition_y = &PrimaryParticlePosition_y,*pPrimaryParticlePosition_z = &PrimaryParticlePosition_z;
    
    Tree->SetBranchAddress("TimeElapsed", &TimeElapsed);
    Tree->SetBranchAddress("EventsElapsed", &EventsElapsed);
    
    Tree->SetBranchAddress("EventNumber",&pEventNumber);
    Tree->SetBranchAddress("PrimaryParticleCharge", &pPrimaryParticleCharge);
    Tree->SetBranchAddress("PrimaryParticleName", &pPrimaryParticleName);
    Tree->SetBranchAddress("PrimaryParticleTime", &pPrimaryParticleTime);
    Tree->SetBranchAddress("PrimaryParticleEnergy", &pPrimaryParticleEnergy);
    
    
    Tree->SetBranchAddress("PrimaryParticleMomentum_x", &pPrimaryParticleMomentum_x);
    Tree->SetBranchAddress("PrimaryParticleMomentum_y", &pPrimaryParticleMomentum_y);
    Tree->SetBranchAddress("PrimaryParticleMomentum_z", &pPrimaryParticleMomentum_z);
    
    
    Tree->SetBranchAddress("PrimaryParticlePosition_x", &pPrimaryParticlePosition_x);
    Tree->SetBranchAddress("PrimaryParticlePosition_y", &pPrimaryParticlePosition_y);
    Tree->SetBranchAddress("PrimaryParticlePosition_z", &pPrimaryParticlePosition_z);
        
    // Get all entries
    Tree->GetEntry(0);
    
    
    EventsElapsed = PrimaryParticleEnergy.size();
    
    std::cout << "Time elapsed in Simulation: " << TimeElapsed << " s" << std::endl;
    std::cout << "Number of Events processed: " << EventsElapsed << std::endl;
    std::cout << "Photon Rate: " << (double)EventsElapsed/TimeElapsed << " /s" << std::endl;
    std::cout << "Photon Flux: " << (double)EventsElapsed/TimeElapsed/SubBoxArea << "/cm^2 /s" << std::endl;
    std::cout << "Readout frames processed: " << (int)(TimeElapsed/0.0016) << std::endl;
    
    // 1D Histograms
    TH1F *IntegratedEnergy = new TH1F ("Integral Energy Spectrum","",EnergyBin,EnergyMinBin,EnergyMaxBin);
    TH1F *StartEnergyDistribution = new TH1F ("Differential Energy Spectrum","",EnergyBin,EnergyMinBin,EnergyMaxBin);
    TH1F *Phi = new TH1F ("Azimuth Angle Intensity","Photon Azimuth Angle Intensity",PhiBin,-180,180);
    TH1F *Theta = new TH1F ("Zenith Angle Intensity","Photon Zenith Angle Intensity",ThetaBin,0,90);
    TH1F *ThetaFlux = new TH1F ("Zenith Angle Differential Flux","Photon Zenith Angle Differential Flux",ThetaBin,0,90);
    
    // 2D Histograms
    TH2D *PhiVsTheta = new TH2D ("Flux of #varphi & #theta","Photon Directional Intensity",PhiBin,-180,180,ThetaBin,0,90);
    
    // Correction Histogram
    TH1F *EnergyCorrection = new TH1F ("Energy Correction","Photon Energy Correction",EnergyBin,EnergyMinBin,EnergyMaxBin);
    
    // Sanity Check histograms
    TF1 *SinF = new TF1("SinTheta","sin(x*4.0*atan(1.0)/180)",0,90);
    TH1F *SinH = new TH1F("SinTheta", "SinTheta", ThetaBin, 0, 90);
//     SinH -> Eval(SinF, "R");
    
    TH2F *Sin2D = new TH2F ("ThetaCorrection","Theta Correction",PhiBin,-180,180,ThetaBin,0,90);
    
    // Energy spectra from mesurements
    TF1 *EnergySpectrumFit1 = new TF1("ESpectrum1","4.3e-3*(x*1000)^(-1.09)*1000",1e-3,1e-1); //  Ryan et al. 1979 pp4 (bib: CosmicGammaSpectrumFit1)
    TF1 *EnergySpectrumFit2 = new TF1("ESpectrum2","1.3e-7*x^(-1.8)*10*exp(-(1000-10)/185)*1000",1e-1,2); // Daniel and Stephens 1974 pp10 (bib: CosmicGammaSpectrumFit2)
    TF1 *EnergySpectrumFit3 = new TF1("ESpectrum3","8.02e-7*x^(-2.6)*14.3*exp(-(1000-14.3)/155)*1000",30,100); // Anand et al. 1973 (bib: CosmicGammaSpectrumFit3)
    
    // Information about the spectrum fits
    std::cout << "Spectrum Fit 1: A = " << 4.3e-3 << " in GeV units: " << 4.3e-3*std::pow(1000,-1.09)*1000 << " DepthFactor " << std::exp(-(1000-1000)/190)  << std::endl;
    std::cout << "Spectrum Fit 2: A = " << 1.3e-7 << " in GeV units: " << 1.3e-7*10000 << " DepthFactor " << std::exp(-(1000-10)/185)  << std::endl;
    std::cout << "Spectrum Fit 3: A = " << 8.02e-7 << " in GeV units: " << 8.02e-7*14.3*1000 << " DepthFactor " << std::exp(-(1000-14.3)/155)  << std::endl;
    
    double ThetaBinSize = 90/(double)ThetaBin/180*PI;
    for(unsigned int i = 1; i <= PhiBin; i++)
    {
        for(unsigned int j = 1; j <= ThetaBin; j++)
        {
            Sin2D -> SetBinContent( i, j, std::sin(ThetaBinSize*((double)j - 0.5)) );
            Sin2D -> SetBinError( i, j, std::sin(ThetaBinSize*((double)j - 0.5)) );
            
            if(i == 1) 
            {
                SinH -> SetBinContent( j, std::sin(ThetaBinSize*((double)j - 0.5)) );
                SinH -> SetBinError( j, 0 );
            }
        }
    }
    
    for (int i=0;i<EventsElapsed;i++)
    {
        // Rotate Vector onto the x-y plane and calculate cos(phi) and sin(phi)
        CosPhi = -PrimaryParticleMomentum_x[i]/sqrt(std::pow(PrimaryParticleMomentum_x[i],2)+std::pow(PrimaryParticleMomentum_z[i],2));
        SinPhi = -PrimaryParticleMomentum_z[i]/sqrt(std::pow(PrimaryParticleMomentum_x[i],2)+std::pow(PrimaryParticleMomentum_z[i],2));
        
        // Use elegant formula for phi extraction (only if -pi <= phi <= pi)
        PhiAngle = std::copysign(std::acos(CosPhi),SinPhi)/PI*180;

        // Fill 1D histograms
        IntegratedEnergy -> Fill(std::log10(PrimaryParticleEnergy[i]*GeV));
        StartEnergyDistribution -> Fill(std::log10(PrimaryParticleEnergy[i]*GeV));
        Phi -> Fill(PhiAngle);
        Theta -> Fill(std::acos(std::abs(PrimaryParticleMomentum_y[i]))/PI*180/*, 1/std::sin(std::acos(-PrimaryParticleMomentum_y[i]))*/); // 1/(sin(thata)) weight
        ThetaFlux -> Fill(acos(std::abs(PrimaryParticleMomentum_y[i]))/PI*180);
        
        // Fill 2D histogram
        PhiVsTheta->Fill(PhiAngle,std::acos(-PrimaryParticleMomentum_y[i])/PI*180/*,1/std::sin(std::acos(-PrimaryParticleMomentum_y[i]))*/);
    }
    
    for(int ebin = 1; ebin <= EnergyBin; ebin++)
    {
        // Integrate bins
        IntegratedEnergy -> SetBinContent(ebin,IntegratedEnergy->Integral(ebin,EnergyBin));
        
        // Calculate error bars correctly
        double BinError = 0.0;
        for(int back_count = EnergyBin; back_count >= ebin; back_count--)
        {
            BinError += std::pow(IntegratedEnergy->GetBinError(back_count),2);
        }
        // Fill errors into bins
        IntegratedEnergy -> SetBinError(ebin,std::sqrt(BinError));
    }
    
    TAxis *axis = StartEnergyDistribution -> GetXaxis();
    int bins = axis -> GetNbins();
    double from = axis -> GetXmin();
    double to = axis -> GetXmax();
    double width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins+1]; //Does not work with just double newbins[bins+1]
    for (int i = 0; i <= bins; i++)
    {
        new_bins[i] = std::pow(10, from + i * width);
        if (i > 0)
        {
            // Fill bin width histogram
            EnergyCorrection -> SetBinContent( i, new_bins[i]-new_bins[i-1] );
        }
    }

//     std::cout << TimeElapsed << std::endl;
    
    // Manually divide histograms, because root has a shitty implementation that messes up the errors
    for(unsigned int i = 1; i < EnergyBin+1; i++)
    {
        StartEnergyDistribution -> SetBinContent( i, StartEnergyDistribution -> GetBinContent(i) / EnergyCorrection -> GetBinContent(i) );
        StartEnergyDistribution -> SetBinError( i, StartEnergyDistribution -> GetBinError(i) / EnergyCorrection -> GetBinContent(i) );
    }
    
    // Scale Histograms, according to thesis cosmic ray theory part
    IntegratedEnergy -> Scale(1/(SubBoxArea*TimeElapsed*2*PI));
    StartEnergyDistribution -> Scale(1/(SubBoxArea*TimeElapsed*2*PI));
    Phi -> Scale(1/(SubBoxArea*TimeElapsed*dPhi));
    Theta -> Divide(SinH);
    Theta ->  Scale(1/(SubBoxArea*TimeElapsed*2*PI*dTheta)); 
    ThetaFlux -> Scale(ThetaBin*2/PI/TimeElapsed/SubBoxArea); // Radian Scale
    PhiVsTheta -> Divide(Sin2D);
    PhiVsTheta -> Scale(1/(SubBoxArea*TimeElapsed*dPhi*dTheta));
    
    // Fit I*cos(theta)^n
    TF1  *FitFunction = new TF1("FitFunction","[0]*cos(x/180*4.0*atan(1.0))^[1]",0,90);
    Theta -> Fit("FitFunction","N");
    
    TLegend* FitLegend = new TLegend(0.50,0.50,0.85,0.85);
    FitLegend->SetNColumns(1);
    FitLegend->SetLineStyle ( 0 );
    FitLegend->SetLineColorAlpha ( 0,0 );
    FitLegend->SetFillStyle ( 0 );
    FitLegend->SetMargin ( 0.2 );
    FitLegend->SetEntrySeparation(0.3);
    FitLegend->SetTextSize(0.03);
    FitLegend->SetTextSize(0.046);
    
    FitLegend -> AddEntry(Theta, " CRY Output","lep");
    FitLegend -> AddEntry(FitFunction, " #chi^{2}-Fit: I(#theta') = I(0) cos^{n_{#gamma}}(#theta') ","L");
    FitLegend -> AddEntry((TObject*)0, " I(0) = (1.083 #pm 0.001) #times 10^{-2}"," ");
    FitLegend -> AddEntry((TObject*)0, " n_{#gamma} = 3.217 #pm 0.004"," ");
    
    // Scaling sanity check by integrating, should give the same for every integral
//     Theta -> Multiply(SinH);
//     PhiVsTheta -> Multiply(Sin2D);
//     std::cout << Phi->Integral()*dPhi << std::endl;
//     std::cout << 2*PI*Theta->Integral()*dTheta << std::endl;
//     std::cout << PhiVsTheta->Integral()*dPhi*dTheta << std::endl;
    
    // Make a legend for fits and data points of differential energy spectrum
    TLegend* Legend = new TLegend(0.20,0.10,0.670,0.20);
    Legend->SetNColumns(1);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.3);
    Legend->SetTextSize(0.03);
    
    Legend -> AddEntry(StartEnergyDistribution, " CRY Output","lep");
    Legend -> AddEntry(EnergySpectrumFit1, " j_{#gamma}(E) #propto E^{-1.09} (Ryan et al. 1979)","L");
    Legend -> AddEntry(EnergySpectrumFit2, " j_{#gamma}(E) #propto E^{-1.8 } (Daniel & Stephens 1974)","L");
    Legend -> AddEntry(EnergySpectrumFit3, " j_{#gamma}(E) #propto E^{-2.6 } (Anand et al. 1973)","L");
    
    TCanvas *c1 = new TCanvas("Integrated Energy","Integrated Energy",900,1600);
    c1 -> SetLogx();
    c1 -> SetLogy();
    c1 -> SetTopMargin(0.05);
    c1 -> SetBottomMargin(0.06);
    c1 -> SetLeftMargin(0.15);
    c1 -> SetRightMargin(0.05);
    TLatex *IntEnergyTitle = new TLatex(3e-3,1.3e-2,"#scale[1.3]{Photon Integral Energy Spectrum}"); // Title
    IntegratedEnergy -> SetLineColor(28);
    IntegratedEnergy -> GetXaxis() -> Set(bins, new_bins);
    IntegratedEnergy -> GetXaxis() -> SetTitle("Photon Kinetic Energy E [GeV]");
    IntegratedEnergy -> GetXaxis() -> SetLabelOffset(-0.022);
    IntegratedEnergy -> GetXaxis() -> SetTitleOffset(0.6);
    IntegratedEnergy -> GetYaxis() -> SetTitle("Photon Integral Energy Spectrum J_{#gamma}(E) [cm^{-2} s^{-1} sr^{-1}]");
    IntegratedEnergy -> GetYaxis() -> SetRangeUser(5e-10,1e-2);
    IntegratedEnergy -> GetYaxis() -> SetTitleOffset(1.4);
    IntegratedEnergy -> Draw();
    IntEnergyTitle -> Draw("SAME");
    c1->SaveAs((StoragePath+"CRYIntegratedEnergy.pdf").c_str());
    
    TCanvas *c2 = new TCanvas("Kinetic Energy","Kinetic Energy",900,1600);
    c2 -> SetLogx();
    c2 -> SetLogy();
    c2 -> SetTopMargin(0.05);
    c2 -> SetBottomMargin(0.06);
    c2 -> SetLeftMargin(0.15);
    c2 -> SetRightMargin(0.05);
    TLatex *StartEnergyTitle = new TLatex(1.5e-3,16,"#scale[1.3]{Photon Differential Energy Spectrum}"); // Title
    StartEnergyDistribution -> SetLineColor(28);
    StartEnergyDistribution -> GetXaxis() -> Set(bins, new_bins);
    StartEnergyDistribution -> GetXaxis() -> SetTitle("Photon Kinetic Energy E [GeV]");
    StartEnergyDistribution -> GetXaxis() -> SetLabelOffset(-0.022);
    StartEnergyDistribution -> GetXaxis() -> SetTitleOffset(0.6);
    StartEnergyDistribution -> GetYaxis() -> SetTitle("Photon Differential Energy Spectrum j_{#gamma}(E) [cm^{-2} s^{-1} sr^{-1} GeV^{-1}]");
    StartEnergyDistribution -> GetYaxis() -> SetRangeUser(3e-11,1e1);
    StartEnergyDistribution -> GetYaxis() -> SetTitleOffset(1.4);
    StartEnergyDistribution -> Draw();
    StartEnergyTitle -> Draw("SAME");
    EnergySpectrumFit1 -> SetLineColor(46);
    EnergySpectrumFit1 -> Draw("SAME");
    EnergySpectrumFit2 -> SetLineColor(9);
    EnergySpectrumFit2 -> Draw("SAME");
    EnergySpectrumFit3 -> SetLineColor(8);
    EnergySpectrumFit3 -> Draw("SAME");
    StartEnergyDistribution -> Draw("SAME");
    gPad->RedrawAxis();
    Legend -> Draw();
    c2->SaveAs((StoragePath+"CRYKineticEnergy.pdf").c_str()); 
    
    TCanvas *c3 = new TCanvas("Phi","Phi",700,500);
//     c3 -> SetLeftMargin(0.12);
    c3 -> SetRightMargin(0.03);
    Phi -> SetLineColor(28);
    Phi -> GetXaxis() -> SetTitle("Azimuth Angle #varphi' [#circ]");
    Phi -> GetYaxis() -> SetTitle("Directional Intensity I_{#gamma}(#varphi') [cm^{-2} s^{-1} sr^{-1}]");
    Phi -> GetYaxis() -> SetTitleOffset(0.8);
    Phi -> GetYaxis() -> SetRangeUser(0,3.2e-3);
    Phi -> GetYaxis() -> SetMaxDigits(3);
    Phi -> Draw();
    c3->SaveAs((StoragePath+"CRYPhi.pdf").c_str()); 
    
    TCanvas *c4 = new TCanvas("Theta Flux","Theta Flux",700,500);
    c4-> SetLeftMargin(0.11);
    c4 -> SetRightMargin(0.09);
    ThetaFlux -> SetLineColor(28);
    ThetaFlux -> GetXaxis() -> SetTitle("Zenith Angle #theta' [#circ]");
    ThetaFlux -> GetYaxis() -> SetTitle("Photon Differential Flux #frac{dJ_{#gamma}(#theta')}{d#theta'} [cm^{-2} s^{-1} rad^{-1}]");
    ThetaFlux -> GetYaxis() -> SetMaxDigits(3);
    ThetaFlux -> GetYaxis() -> SetTitleOffset(0.85);
    ThetaFlux -> Draw();
//     c4->SaveAs((StoragePath+"CRYThetaFlux.pdf").c_str()); 
    
    TCanvas *c5 = new TCanvas("Theta","Theta",700,500);
//     c5 -> SetLeftMargin(0.11);
    c5 -> SetRightMargin(0.03);
    Theta -> SetLineColor(28);
    Theta -> GetXaxis() -> SetTitle("Zenith Angle #theta' [#circ]");
    Theta -> GetYaxis() -> SetTitle("Directional Intensity I_{#gamma}(#theta') [cm^{-2} s^{-1} sr^{-1}]");
    Theta -> GetYaxis() -> SetTitleOffset(0.9);
    Theta -> GetYaxis() -> SetRangeUser(0,12e-3);
    Theta -> GetYaxis() -> SetMaxDigits(3);
    Theta -> Draw();
    FitFunction -> SetLineColor(12);
    FitFunction -> Draw("SAME");
    Theta -> Draw("SAME");
    gPad->RedrawAxis();
    FitLegend -> Draw();
    c5->SaveAs((StoragePath+"CRYTheta.pdf").c_str()); 
    
    TCanvas *c6 = new TCanvas("Gamma Flux","Gamma Flux",700,500);
    c6 -> SetLeftMargin(0.09);
    c6 -> SetRightMargin(0.15);
    PhiVsTheta -> GetXaxis() -> SetTitle("Azimuth Angle #varphi' [#circ]");
    PhiVsTheta -> GetYaxis() -> SetTitle("Zenith Angle #theta' [#circ]");
    PhiVsTheta -> GetZaxis() -> SetTitle("Photon Directional Intensity I_{#gamma} [cm^{-2} s^{-1} sr^{-1}]");
    PhiVsTheta -> GetZaxis() -> SetRangeUser(0,1.2e-2);
    PhiVsTheta -> GetZaxis() -> SetMaxDigits(2);
    PhiVsTheta -> GetYaxis() -> SetTitleOffset(0.9);
    PhiVsTheta -> SetTickLength(0.01,"Z");
    PhiVsTheta -> Draw("colz");
    PhiVsTheta -> Draw("colz");
//     c6->SaveAs((StoragePath+"CRYPhiTheta.pdf").c_str());
}
