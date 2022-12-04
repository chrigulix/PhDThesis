// Particle analysis Program

// C++ Headers
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <random>
#include <utility>
#include <array>

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
#include "TRandom.h"
#include "TLine.h"

// Define globals
void analysis(std::string FileName);
const double PI = 4.0*atan(1.0);
const double e = 1.60217733e-19;
const double epsilon_0 = 8.85418782e-12;
const double epsilon_r = 1.53;

// Unit Converstion
const double GeV = 0.001; // MeV -> GeV
const double SqdegToSr = std::pow(PI/180,2); // Square Degrees -> Steradian

// TPC active colume boundaries
const double ActiveVolumeX = 256.35; // cm
const double ActiveVolumeY = 233; // cm
const double ActiveVolumeZ = 1036.8; // cm

// Fiducial cuts
const double FiducialCutTop = 20; // cm
const double FiducialCutSides = 10; // cm

// Readout window length
const double ReadoutWindow = 0.00225; // s

// Input file name and output path
const std::string FileName = "/home/christoph/Programming/CosmicAnalysis/Cosmic/plots/v14/MC_uboone.root";
const std::string StoragePath = "../images/CosmicGammaBackground/";

void DrawCosmicGammaBackground()
{
    // Cut values
    double EnergyCutValue = 0.2; // GeV
    
    // Binning
    int EnergyBin = 30;
    int EnergyMinBin = -4; // log number i.e. 10^-3 GeV
    int EnergyMaxBin = 1; // log number i.e. 10^2 GeV
    int EdepPathBin = 45;
    int XPosBin = 256/8;
    int YPosBin = 233/8;
    int ZPosBin = 1037/8;
    
    int EPhiBin = 30;
    int EThetaBin = 15;
    
    double dEPhi = 2*PI/EPhiBin;
    double dETheta = PI/2/EThetaBin;
    
    // Max and Min
    double PathMax = 450; // cm
    double DepMax = 1.0;//1040*2.1;

    // Open File
    std::string TreeName = "MCTree";
    TFile* File = new TFile (FileName.c_str());
    
    // Copy root tree from file into memory
    TTree* Tree = (TTree*)File->Get(TreeName.c_str());
 
    // Time and Event counters
    int EventsElapsed;
    double TimeElapsed;
    
    // Initialize vectors and pointers, for tree branch reading
    std::vector<int> HitNumber;
    std::vector<int>*pHitNumber = &HitNumber;
    
    std::vector<int> ComptonHitNumber;
    std::vector<int>*pComptonHitNumber = &ComptonHitNumber;

    std::vector<std::string> ProcessName;
    std::vector<std::string>*pProcessName = &ProcessName;

    std::vector<double> HitEnergy;
    std::vector<double>*pHitEnergy = &HitEnergy;
 
    std::vector<double> EnergyDeposit;
    std::vector<double>*pEnergyDeposit = &EnergyDeposit;
 
    std::vector<double> TrackLength;
    std::vector<double>*pTrackLength = &TrackLength;
 
    std::vector<double> HitMomentum_x, HitMomentum_y, HitMomentum_z;    
    std::vector<double> *pHitMomentum_x = &HitMomentum_x,*pHitMomentum_y = &HitMomentum_y,*pHitMomentum_z = &HitMomentum_z;
 
    std::vector<double> HitPosition_x, HitPosition_y, HitPosition_z;
    std::vector<double>*pHitPosition_x = &HitPosition_x,*pHitPosition_y = &HitPosition_y,*pHitPosition_z = &HitPosition_z;
 
    std::vector<double> ComptonEnergy;
    std::vector<double> *pComptonEnergy = &ComptonEnergy;
    
    std::vector<double> ComptonPosition_x, ComptonPosition_y, ComptonPosition_z;
    std::vector<double>*pComptonPosition_x = &ComptonPosition_x,*pComptonPosition_y = &ComptonPosition_y,*pComptonPosition_z = &ComptonPosition_z;
    
    std::vector<double> ComptonMomentum_x, ComptonMomentum_y, ComptonMomentum_z;
    std::vector<double>*pComptonMomentum_x = &ComptonMomentum_x,*pComptonMomentum_y = &ComptonMomentum_y,*pComptonMomentum_z = &ComptonMomentum_z;

    // Set branch addresses, and read data into memory
    Tree->SetBranchAddress("TimeElapsed", &TimeElapsed);
    Tree->SetBranchAddress("EventsElapsed", &EventsElapsed);
    Tree->SetBranchAddress("HitNumber", &pHitNumber);
    Tree->SetBranchAddress("ComptonHitNumber", &pComptonHitNumber);
    Tree->SetBranchAddress("ProcessName", &pProcessName);
    Tree->SetBranchAddress("EnergyDeposit", &pEnergyDeposit);
    Tree->SetBranchAddress("TrackLength", &pTrackLength);
    
    // Gamma related variables
    Tree->SetBranchAddress("HitEnergy", &pHitEnergy);
    Tree->SetBranchAddress("HitMomentum_x", &pHitMomentum_x);
    Tree->SetBranchAddress("HitMomentum_y", &pHitMomentum_y); 
    Tree->SetBranchAddress("HitMomentum_z", &pHitMomentum_z);
    Tree->SetBranchAddress("HitPosition_x", &pHitPosition_x);
    Tree->SetBranchAddress("HitPosition_y", &pHitPosition_y);
    Tree->SetBranchAddress("HitPosition_z", &pHitPosition_z);

    // Electron related variables
    Tree->SetBranchAddress("ComptonEnergy", &pComptonEnergy);
    Tree->SetBranchAddress("ComptonPosition_x", &pComptonPosition_x);
    Tree->SetBranchAddress("ComptonPosition_y", &pComptonPosition_y);
    Tree->SetBranchAddress("ComptonPosition_z", &pComptonPosition_z);
    Tree->SetBranchAddress("ComptonMomentum_x", &pComptonMomentum_x);
    Tree->SetBranchAddress("ComptonMomentum_y", &pComptonMomentum_y);
    Tree->SetBranchAddress("ComptonMomentum_z", &pComptonMomentum_z);
    
    // Read whole tree
    Tree->GetEntry(0);
 
    int NumberOfHits = HitNumber.size();
    int NumberOfComptons = ComptonHitNumber.size();
    
    // Output summery about the simulation
    std::cout << "--------------------Simulation Details--------------------" << std::endl;
    std::cout << "Number of cosmic gammas generated: " << EventsElapsed << std::endl;
    std::cout << "Number of gamma hits in active volume: " << HitNumber.size() << std::endl;
    std::cout << "Number of compton and photoelectrons: " << ComptonHitNumber.size() << std::endl;
    std::cout << "Time elapsed [s]: " << TimeElapsed << std::endl;
    std::cout << "Gamma hit rate [s^-1]: " << HitNumber.size()/TimeElapsed << std::endl;
    std::cout << "Compton and photoelectron rate [s^-1]: " << ComptonHitNumber.size()/TimeElapsed << std::endl;
    std::cout << "Gamma hit rate per readout (2.25 ms): " << EventsElapsed/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Compton and photoelectron rate per readout (2.25 ms): " << ComptonHitNumber.size()/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Compton rate per readout error " << 1/std::sqrt(ComptonHitNumber.size())*ComptonHitNumber.size()/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Readout frames processed (2.25 ms): " << TimeElapsed/ReadoutWindow << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    
    // Initialize electron histograms
    TH1F *ElectronEnergyFiducial = new TH1F ("DifferentialEnergyRateFiducial","",EnergyBin,EnergyMinBin,EnergyMaxBin);
    TH1F *IntElectronEnergyFiducial = new TH1F ("IntegralEnergyRateFiducial","",EnergyBin,EnergyMinBin,EnergyMaxBin);
    
    // Fiducial
    TH1F *ElectronPhi = new TH1F ("ElectronPhi","Differential Directional Rate in #varphi",EPhiBin,-180,180);
    TH1F *ElectronTheta = new TH1F ("ElectronTheta","Differential Directional Rate in #theta",EThetaBin,0,180);
    TH2F *ElectronPhiTheta = new TH2F ("ElectronPhiTheta","Differential Directional Rate",EPhiBin,-180,180,EThetaBin,0,180);
    
    // Energy Cut 200 MeV
    TH1F *ElectronPhiECut = new TH1F ("ElectronPhiECut","Differential Directional Rate in #varphi (E #geq 200 MeV)",EPhiBin,-180,180);
    TH1F *ElectronThetaECut = new TH1F ("ElectronThetaECut","Differential Directional Rate in #theta (E #geq 200 MeV)",EThetaBin,0,180);
    TH2F *ElectronPhiThetaECut = new TH2F ("ElectronPhiThetaECut","Differential Directional Rate (E #geq 200 MeV)",EPhiBin,-180,180,EThetaBin,0,180);
    
    TH1F *ElectronTrackLength = new TH1F ("Track","Differential Track Length Rate",EdepPathBin,0,PathMax);
    
    // Create sin(theta) scaling histogram for ElectronPhiTheta and sanity checks
    TH1F *ElectronSinH = new TH1F("ESinTheta", "SinTheta", EThetaBin, 0, 180);
    
    // 2D Scaling Histograms
    TH2F *ElectronSin2D = new TH2F ("ThetaCorrection","Theta Correction",EPhiBin,-180,180,EThetaBin,0,180);

    // Fill scaling histograms for theta distributions
    double EThetaBinSize = 180/(double)EThetaBin/180*PI;
    for(unsigned int i = 1; i <= EPhiBin; i++)
    {
        for(unsigned int j = 1; j <= EThetaBin; j++)
        {
            ElectronSin2D -> SetBinContent( i, j, std::sin(EThetaBinSize*((double)j - 0.5)) );
            ElectronSin2D -> SetBinError(i,j,0);
            if(i == 1)
            {
                ElectronSinH -> SetBinContent( j, std::sin(EThetaBinSize*((double)j - 0.5)) );
                ElectronSinH -> SetBinError(j, 0);
            }
        }
    }
    
    // Fiducial cut diagrams
    TH2F *ElectronPosXY = new TH2F ("Upstream View","Upstream View of Detector",XPosBin,0,ActiveVolumeX,YPosBin,-ActiveVolumeY/2,ActiveVolumeY/2);
    TH2F *ElectronPosZY = new TH2F ("Side View","Side View of Detector",ZPosBin,0,ActiveVolumeZ,YPosBin,-ActiveVolumeY/2,ActiveVolumeY/2);
    
    std::vector<TLine*> FiducialBoxXY;
    FiducialBoxXY.push_back(new TLine(FiducialCutSides,ActiveVolumeY/2-FiducialCutTop,ActiveVolumeX-FiducialCutSides,ActiveVolumeY/2-FiducialCutTop));
    FiducialBoxXY.push_back(new TLine(ActiveVolumeX-FiducialCutSides,ActiveVolumeY/2-FiducialCutTop,ActiveVolumeX-FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides));
    FiducialBoxXY.push_back(new TLine(ActiveVolumeX-FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides,FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides));
    FiducialBoxXY.push_back(new TLine(FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides,FiducialCutSides,ActiveVolumeY/2-FiducialCutTop));
    
    std::vector<TLine*> FiducialBoxZY;
    FiducialBoxZY.push_back(new TLine(FiducialCutSides,ActiveVolumeY/2-FiducialCutTop,ActiveVolumeZ-FiducialCutSides,ActiveVolumeY/2-FiducialCutTop));
    FiducialBoxZY.push_back(new TLine(ActiveVolumeZ-FiducialCutSides,ActiveVolumeY/2-FiducialCutTop,ActiveVolumeZ-FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides));
    FiducialBoxZY.push_back(new TLine(ActiveVolumeZ-FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides,FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides));
    FiducialBoxZY.push_back(new TLine(FiducialCutSides,-ActiveVolumeY/2+FiducialCutSides,FiducialCutSides,ActiveVolumeY/2-FiducialCutTop));
    
    // Temporary histogram for log-bin energy scaling
    TH1F *EnergyCorrection = new TH1F ("Energy Correction","Energy Correction",EnergyBin,EnergyMinBin,EnergyMaxBin);
    
    // Random number generator sanity check histogram
    TH1F *Random = new TH1F ("Random","Random",100,0,1);
    TLine* RandomAverage = new TLine (0,ComptonEnergy.size()/100,1,ComptonEnergy.size()/100);
    
    
    // Temporary variables 
    double ElectronCosPhi;
    double ElectronSinPhi;
    double ElectronPhiAngle;
    double ElectronThetaAngle;
    
    // Initialize random generator for x-fiducial cut
    unsigned long int seed = 1337420;
    TRandomRanlux48 RandomGenerator(seed); // Ranlux48 Best but slow
    double RandomNumber = 1;
    
    // Check numbers
    unsigned long int FiducialCut = 0;
    unsigned long int FiducialSurvived = 0;
    unsigned long int EnergyCut = 0;
    unsigned long int EnergySurvived = 0;
    unsigned long int ThetaCut = 0;
    unsigned long int ThetaSurvived = 0;
    
    
    // Initialize random x position to simulate drift
    double XPosition;

    for (int e_index = 0; e_index < ComptonEnergy.size(); e_index++)
    {
        // Generate random number
        RandomNumber = RandomGenerator.Rndm();
        
        // Sanity check, fill random number histogram
        Random -> Fill(RandomNumber);
        
        // Generate x-Position
        XPosition = RandomNumber*ActiveVolumeX;
        
        // Change units, ATTENTION shit is stored in mm and not cm
        ComptonPosition_x[e_index] /= 10.0;
        ComptonPosition_y[e_index] /= 10.0;
        ComptonPosition_z[e_index] /= 10.0;
        
        // Fill 2D Histograms with positions        
        ElectronPosXY -> Fill(XPosition,ComptonPosition_y[e_index]);
        ElectronPosZY -> Fill(ComptonPosition_z[e_index]+ActiveVolumeZ/2,ComptonPosition_y[e_index]);
        
        // Calculate cosine and sine of phi
        ElectronCosPhi = ComptonMomentum_x[e_index]/std::sqrt( std::pow(ComptonMomentum_x[e_index],2) + std::pow(ComptonMomentum_y[e_index],2) );
        ElectronSinPhi = ComptonMomentum_y[e_index]/std::sqrt( std::pow(ComptonMomentum_x[e_index],2) + std::pow(ComptonMomentum_y[e_index],2) );
        
        // Use acos of CosPhi as angle and multiply the sign of SinPhi. This gives a distribution from -180 -> 180 degrees 
        ElectronPhiAngle = std::copysign(std::acos(ElectronCosPhi),ElectronSinPhi)/PI*180;
        
        // Use acos of Momentum_z as theta angle
        ElectronThetaAngle = std::acos(ComptonMomentum_z[e_index])/PI*180;
        
        // Cut count ++
        FiducialCut++;
        EnergyCut++;
        ThetaCut++;
        
        // Apply fiducial cuts
        if( (ComptonPosition_y[e_index] <= (ActiveVolumeY/2-FiducialCutTop)) && 
            (std::abs(ComptonPosition_y[e_index]) <= (ActiveVolumeY/2-FiducialCutSides))  && 
            (std::abs(ComptonPosition_z[e_index]) <= (ActiveVolumeZ/2-FiducialCutSides)) &&
            (XPosition >= FiducialCutSides) &&  (XPosition <= (ActiveVolumeX-FiducialCutSides)) )
        {
            // Fill histograms after events made it through the cut
            ElectronEnergyFiducial -> Fill(std::log10(ComptonEnergy[e_index]*GeV));
            IntElectronEnergyFiducial -> Fill(std::log10(ComptonEnergy[e_index]*GeV));
            ElectronPhi -> Fill(ElectronPhiAngle);
            ElectronTheta -> Fill(ElectronThetaAngle);
            // Fill PhiTheta with a weight of 1/sin(theta) for per steradian scaling (not so nice result)
            ElectronPhiTheta -> Fill( ElectronPhiAngle, ElectronThetaAngle);
            ElectronTrackLength -> Fill(TrackLength[e_index]/10);
            
            // Count events
            FiducialSurvived++;
            FiducialCut--; // 4287746
            
            // Energy Cut
            if(ComptonEnergy[e_index]*GeV >= EnergyCutValue)
            {
                ElectronPhiECut -> Fill(ElectronPhiAngle);
                ElectronThetaECut -> Fill(ElectronThetaAngle/*, 1/std::sin(ElectronThetaAngle/180*PI)*/);
                ElectronPhiThetaECut -> Fill( ElectronPhiAngle, ElectronThetaAngle);
                
                // Count events
                EnergySurvived++;
                EnergyCut--;
                
                // Theta cut
                if(ElectronThetaAngle <= 45)
                {
                    ThetaSurvived++;
                    ThetaCut--;
                }
            }
        }
    }
    
    // Energy rate for plots
    double EnergyRateAtCut = (double)EnergySurvived/TimeElapsed; 
    
    // Output summery about the cuts performed post simulation
    std::cout << "-----------------------Cut Summary-----------------------" << std::endl;
    std::cout << "Events surviving fiducial cut: " << FiducialSurvived << std::endl;
    std::cout << "Events cut by fiducial cut: " << FiducialCut << std::endl;
    std::cout << "Fraction survicing fiducial cut: " << (double)FiducialSurvived/(double)(FiducialSurvived+FiducialCut) << std::endl;
    std::cout << "Fiducial rate [s^-1] : " << (double)FiducialSurvived/TimeElapsed << std::endl;
    std::cout << "Fiducial events per readout window : " << (double)FiducialSurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Error per readout window : " << 1/std::sqrt(FiducialSurvived)*(double)FiducialSurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << "Events surviving energy cut: " << EnergySurvived << std::endl;
    std::cout << "Events cut by energy cut: " << EnergyCut << std::endl;
    std::cout << "Fraction survicing energy cut: " << (double)EnergySurvived/(double)(FiducialSurvived) << std::endl;
    std::cout << "Rate after all cuts [s^-1] : " << (double)EnergySurvived/TimeElapsed << std::endl;
    std::cout << "Events per readout window after cuts : " << (double)EnergySurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Error readout window after cuts : " << 1/std::sqrt(EnergySurvived)*(double)EnergySurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    std::cout << "Events surviving theta cut: " << ThetaSurvived << std::endl;
    std::cout << "Events cut by theta cut: " << ThetaCut << std::endl;
    std::cout << "Fraction survicing theta cut: " << (double)ThetaSurvived/(double)(EnergySurvived) << std::endl;
    std::cout << "Rate after all cuts [s^-1] : " << (double)ThetaSurvived/TimeElapsed << std::endl;
    std::cout << "Events per readout window after cuts : " << (double)ThetaSurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "Events per readout window after cuts : " << 1/std::sqrt(ThetaSurvived)*(double)ThetaSurvived/TimeElapsed*ReadoutWindow << std::endl;
    std::cout << "----------------------------------------------------------" << std::endl;
    
    // Loop over energy bins in order to integrate the energy
    for(int ebin = 1; ebin <= EnergyBin; ebin++)
    {
        // Integrate bins
        IntElectronEnergyFiducial -> SetBinContent(ebin,IntElectronEnergyFiducial->Integral(ebin,EnergyBin));
        
        // Calculate error bars correctly
        double BinError = 0.0;
        for(int back_count = EnergyBin; back_count >= ebin; back_count--)
        {
            BinError += std::pow(IntElectronEnergyFiducial->GetBinError(back_count),2);
        }
        // Fill errors into bins
        IntElectronEnergyFiducial -> SetBinError(ebin,std::sqrt(BinError));
    }
    
    // Temp
    Random -> Sumw2();
    
//     ElectronEnergyFiducial -> Sumw2();
//     IntElectronEnergyFiducial -> Sumw2();
    
//     ElectronPhi -> Sumw2();
//     ElectronTheta -> Sumw2();
    

    // Produce logarithmic x-axis bins
    TAxis *axis = ElectronEnergyFiducial -> GetXaxis();
    int bins = axis -> GetNbins();
    double from = axis -> GetXmin();
    double to = axis -> GetXmax();
    double width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins+1]; //Does not work with just double newbins[bins+1
    // Fill new axis bins
    for (int i = 0; i <= bins; i++)
    {
        new_bins[i] = std::pow(10, from + i * width);
        if (i > 0)
        {
            // Fill bin width histogram
            EnergyCorrection -> SetBinContent( i, new_bins[i]-new_bins[i-1] );
        }
    }
    
    // Scaling factors
    double BinVolumeXY = ActiveVolumeX/(double)XPosBin * ActiveVolumeY/(double)YPosBin * ActiveVolumeZ; // cm^3
    double BinVolumeZY = ActiveVolumeX * ActiveVolumeY/(double)YPosBin * ActiveVolumeZ/(double)ZPosBin; // cm^3
    
    // Manually divide histograms, because root has a shitty implementation
    for(unsigned int i = 1; i <= EnergyBin; i++)
    {
        ElectronEnergyFiducial -> SetBinContent( i, ElectronEnergyFiducial -> GetBinContent(i) / EnergyCorrection -> GetBinContent(i) );
        ElectronEnergyFiducial -> SetBinError( i, ElectronEnergyFiducial -> GetBinError(i) / EnergyCorrection -> GetBinContent(i) );
    }
    
    // Scaling histograms
    ElectronEnergyFiducial -> Scale(1/TimeElapsed);
    IntElectronEnergyFiducial -> Scale(1/TimeElapsed);
    
    ElectronPosXY -> Scale(1/TimeElapsed/BinVolumeXY);
    ElectronPosZY -> Scale(1/TimeElapsed/BinVolumeZY);
    // Get Minimum and maximum
    double PosMin = std::min(ElectronPosXY->GetMinimum(),ElectronPosZY->GetMinimum());
    double PosMax = std::min(ElectronPosXY->GetMaximum(),ElectronPosZY->GetMaximum());
    
    ElectronPhi -> Scale(1/(TimeElapsed*2*dEPhi)); 
    ElectronTheta -> Divide(ElectronSinH);
    ElectronTheta -> Scale(1/(TimeElapsed*2*PI*dETheta));
    ElectronPhiECut -> Scale(1/(TimeElapsed*2*dEPhi));
    ElectronThetaECut -> Divide(ElectronSinH);
    ElectronThetaECut -> Scale(1/(TimeElapsed*2*PI*dETheta));
    
    
    ElectronPhiTheta -> Divide(ElectronSin2D);
    ElectronPhiTheta -> Scale(1/(TimeElapsed*dEPhi*dETheta));
    ElectronPhiThetaECut -> Divide(ElectronSin2D);
    ElectronPhiThetaECut -> Scale(1/(TimeElapsed*dEPhi*dETheta));
    
    ElectronTrackLength -> Scale(EdepPathBin/PathMax/TimeElapsed);
    
    // Scaling sanity check, all three numbers should be the same (or at least close to the same)
//     ElectronTheta -> Multiply(ElectronSinH);
//     ElectronPhiTheta -> Multiply(ElectronSin2D);
//     std::cout << 2*ElectronPhi->Integral()*dEPhi << std::endl;
//     std::cout << 2*PI*ElectronTheta->Integral()*dETheta << std::endl;
//     std::cout << ElectronPhiTheta->Integral()*dEPhi*dETheta << std::endl;

    
    TCanvas *c1 = new TCanvas("Single Electron Rate","Single Electron Rate ",900,1600);
    c1 -> SetLogx();
    c1 -> SetLogy();
    c1 -> SetTopMargin(0.05);
    c1 -> SetBottomMargin(0.06);
    c1 -> SetLeftMargin(0.16);
    c1 -> SetRightMargin(0.05);
    TLatex *ElectronEnergyTitle = new TLatex(1.3e-4,1.5e8,"#scale[1.3]{Electron Differential Energy Spectrum}"); // Title
    TBox *ExclusionBox = new TBox(1e-4,5e-4,EnergyCutValue,1e8);
    ExclusionBox -> SetFillColorAlpha(16,1);
    ElectronEnergyFiducial -> SetLineColor(kRed+1);
    ElectronEnergyFiducial -> GetXaxis() -> Set(bins, new_bins);
    ElectronEnergyFiducial -> GetXaxis() -> SetTitle("Electron Kinetic Energy E [GeV]");
    ElectronEnergyFiducial -> GetXaxis() -> SetLabelOffset(-0.022);
    ElectronEnergyFiducial -> GetXaxis() -> SetTitleOffset(0.6);
    ElectronEnergyFiducial -> GetYaxis() -> SetTitle("Electron Differential Energy Rate #frac{dr_{e}(E)}{dE} [s^{-1} GeV^{-1}]");
    ElectronEnergyFiducial -> GetYaxis() -> SetRangeUser(5e-4,1e8);
    ElectronEnergyFiducial -> GetYaxis() -> SetTitleOffset(1.4);
    ElectronEnergyFiducial -> Draw();
    ExclusionBox -> Draw("SAME");
    ElectronEnergyTitle -> Draw("SAME");
    ElectronEnergyFiducial -> Draw("SAME");
    gPad->RedrawAxis();
    c1->SaveAs((StoragePath+"ElectronEnergy.pdf").c_str());
    
    
    TCanvas *c2 = new TCanvas("Integrated Single Electron Rate","Electron Integrated Single Rate",900,1600);
    c2 -> SetLogx();
    c2 -> SetLogy();
    c2 -> SetTopMargin(0.05);
    c2 -> SetBottomMargin(0.06);
    c2 -> SetLeftMargin(0.15);
    c2 -> SetRightMargin(0.05);
    TLatex *IntEnergyTitle = new TLatex(2.3e-4,1.3e4,"#scale[1.3]{Electron Integral Energy Spectrum}"); // Title
    TLine *RateLine = new TLine(7e-5,EnergyRateAtCut,EnergyCutValue,EnergyRateAtCut);
    RateLine -> SetLineColor(12);
    RateLine -> SetLineStyle(7);
    RateLine -> SetLineWidth(2);
    TLine *CutEnergyLine = new TLine(EnergyCutValue,5e-4,EnergyCutValue,EnergyRateAtCut);
    CutEnergyLine -> SetLineColor(12);
    CutEnergyLine -> SetLineStyle(7);
    CutEnergyLine -> SetLineWidth(2);
    TLatex *RateText = new TLatex(2.5e-5,0.65,"#scale[1.0]{0.75}");
    IntElectronEnergyFiducial -> SetLineColor(kRed+1);
    IntElectronEnergyFiducial -> GetXaxis() -> Set(bins, new_bins);
    IntElectronEnergyFiducial -> GetXaxis() -> SetTitle("Electron Kinetic Energy E [GeV]");
    IntElectronEnergyFiducial -> GetXaxis() -> SetLabelOffset(-0.022);
    IntElectronEnergyFiducial -> GetXaxis() -> SetTitleOffset(0.6);
    IntElectronEnergyFiducial -> GetYaxis() -> SetTitle("Electron Integral Energy Rate R_{e} [s^{-1}]");
    IntElectronEnergyFiducial -> GetYaxis() -> SetRangeUser(5e-4,1e4);
    IntElectronEnergyFiducial -> GetYaxis() -> SetTitleOffset(1.4);
    IntElectronEnergyFiducial -> Draw();
    IntEnergyTitle -> Draw("SAME");
    RateLine -> Draw("SAME");
    CutEnergyLine -> Draw("SAME");
    RateText -> Draw("SAME");
    gPad->RedrawAxis();
    c2->SaveAs((StoragePath+"IntegratedElectronEnergy.pdf").c_str());
    
    
    TCanvas *c3 = new TCanvas("Electron Distribution Upstream","Electron Distribution Upstream",2*256*1.07,2*233);
    c3 -> SetLogz();
    c3 -> SetLeftMargin(0.06);
    c3 -> SetRightMargin(0.21);
    ElectronPosXY -> GetXaxis() -> SetTitle("X Coordinate [cm]");
    ElectronPosXY -> GetZaxis() -> SetTitle("Differential Volume Rate #frac{dr_{e}}{dV} [cm^{-3} s^{-1}]");
    ElectronPosXY -> GetZaxis() -> SetTickLength(0.02);
    ElectronPosXY -> GetZaxis() -> SetMaxDigits(1);
    ElectronPosXY -> GetZaxis() -> SetRangeUser(PosMin,PosMax);
    ElectronPosXY -> GetZaxis() -> SetTitleOffset(1.3);
    ElectronPosXY -> Draw("colz");
    gPad->Update();
    TPaletteAxis *paletteXY = (TPaletteAxis*)ElectronPosXY->GetListOfFunctions()->FindObject("palette");
    paletteXY->SetX1NDC(0.81);
    paletteXY->SetX2NDC(0.86);
    paletteXY->SetY1NDC(0.1);
    paletteXY->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();
    for(auto Line : FiducialBoxXY)
    {
        Line -> SetLineColor(kRed+1);
        Line -> SetLineWidth(3);
        Line -> Draw("SAME");
    }
    c3->SaveAs((StoragePath+"UpstreamView.pdf").c_str());
    
    
    TCanvas *c4 = new TCanvas("Electron Distribution Side","Electron Distribution Side",1037*0.87,2*233);
    c4 -> SetLogz();
    c4 -> SetLeftMargin(0.06);
    c4 -> SetRightMargin(0.01);
    ElectronPosZY -> GetXaxis() -> SetTitle("Z-Coordinate [cm]");
    ElectronPosZY -> GetYaxis() -> SetTitle("Y-Coordinate [cm]");
    ElectronPosZY -> GetYaxis() -> SetTitleOffset(0.6);
    ElectronPosZY -> GetZaxis() -> SetRangeUser(PosMin,PosMax);
    ElectronPosZY -> GetZaxis() -> SetTickLength(0.02);
    ElectronPosZY -> Draw("colz");
    gPad->Update();
    TPaletteAxis *paletteZY = (TPaletteAxis*)ElectronPosZY->GetListOfFunctions()->FindObject("palette");
    paletteZY->SetX1NDC(1);
    paletteZY->SetX2NDC(1);
    paletteZY->SetY1NDC(1);
    paletteZY->SetY2NDC(1);
    gPad->Modified();
    gPad->Update();
    for(auto Line : FiducialBoxZY)
    {
        Line -> SetLineColor(kRed+1);
        Line -> SetLineWidth(3);
        Line -> Draw("SAME");
    }
    c4->SaveAs((StoragePath+"SideView.pdf").c_str());
    
    
    TCanvas *c5 = new TCanvas("Electron Phi Angle","Electron Phi Angle",700,500);
    c5 -> SetLeftMargin(0.11);
    c5 -> SetRightMargin(0.03);
    ElectronPhi -> SetLineColor(kRed+1);
    ElectronPhi -> GetXaxis() -> SetTitle("#varphi Angle [#circ]");
    ElectronPhi -> GetYaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#varphi)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronPhi -> GetYaxis() -> SetTitleOffset(0.9);
    ElectronPhi -> GetYaxis() -> SetMaxDigits(3);
    ElectronPhi -> GetYaxis() -> SetRangeUser(0.6e3,1.8e3);
    ElectronPhi -> Draw();
    c5->SaveAs((StoragePath+"ElectronPhi.pdf").c_str());
    
    
    TCanvas *c6 = new TCanvas("Electron Theta Angle","Electron Theta Angle",700,500);
    c6 -> SetLeftMargin(0.11);
    c6 -> SetRightMargin(0.03);
    ElectronTheta -> SetLineColor(kRed+1);
    ElectronTheta -> GetXaxis() -> SetTitle("#theta Angle [#circ]");
    ElectronTheta -> GetYaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#theta)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronTheta -> GetYaxis() -> SetTitleOffset(0.9);
    ElectronTheta -> GetYaxis() -> SetMaxDigits(3);
    ElectronTheta -> GetYaxis() -> SetRangeUser(1.4e3,2.4e3);
    ElectronTheta -> Draw();
    c6->SaveAs((StoragePath+"ElectronTheta.pdf").c_str());
    
    
    TCanvas *c7 = new TCanvas("Electron Angle Distribution","Electron Angular Distribution",1600,900);
    c7 -> SetLeftMargin(0.07);
    c7 -> SetRightMargin(0.14);
    ElectronPhiTheta -> GetXaxis() -> SetTitle("#varphi Angle [#circ]");
    ElectronPhiTheta -> GetXaxis() -> SetMaxDigits(3);
    ElectronPhiTheta -> GetYaxis() -> SetTitle("#theta Angle [#circ]");
    ElectronPhiTheta -> GetYaxis() -> SetTitleOffset(0.7);
    ElectronPhiTheta -> GetYaxis() -> SetMaxDigits(3);
    ElectronPhiTheta -> GetZaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#varphi,#theta)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronPhiTheta -> GetZaxis() -> SetTitleOffset(0.8);
    ElectronPhiTheta -> Draw("COLZ");
    gPad->Update();
    TPaletteAxis* EPhiThetaPalette = (TPaletteAxis*)ElectronPhiTheta -> GetListOfFunctions() -> FindObject("palette");
    EPhiThetaPalette -> GetAxis() -> SetMaxDigits(2);
    gPad->Modified();
    gPad->Update();
    c7->SaveAs((StoragePath+"ElectronPhiTheta.pdf").c_str());
    
    
    TCanvas *c8 = new TCanvas("Electron Phi Angle E-Cut","Electron Phi Angle (E > 200 MeV)",700,500);
    c8 -> SetLeftMargin(0.11);
    c8 -> SetRightMargin(0.03);
    ElectronPhiECut -> SetLineColor(9);
    ElectronPhiECut -> GetXaxis() -> SetTitle("#varphi Angle [#circ]");
    ElectronPhiECut -> GetXaxis() -> SetMaxDigits(3);
    ElectronPhiECut -> GetYaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#varphi)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronPhiECut -> GetYaxis() -> SetTitleOffset(0.9);
    ElectronPhiECut -> GetYaxis() -> SetMaxDigits(1);
    ElectronPhiECut -> GetYaxis() -> SetRangeUser(0,0.3);
    ElectronPhiECut -> Draw();
    c8->SaveAs((StoragePath+"ElectronPhiECut.pdf").c_str());
    
    
    TCanvas *c9 = new TCanvas("Electron Theta Angle E-Cut","Electron Theta Angle  (E > 200 MeV)",700,500);
    c9 -> SetLeftMargin(0.11);
    c9 -> SetRightMargin(0.03);
    ElectronThetaECut -> SetLineColor(9);
    ElectronThetaECut -> GetXaxis() -> SetMaxDigits(3);
    ElectronThetaECut -> GetXaxis() -> SetTitle("#theta Angle [#circ]");
    ElectronThetaECut -> GetYaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#theta)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronThetaECut -> GetYaxis() -> SetTitleOffset(0.9);
    ElectronThetaECut -> GetYaxis() -> SetMaxDigits(1);
    ElectronThetaECut -> GetYaxis() -> SetRangeUser(0,0.3);
    ElectronThetaECut -> Draw();
    c9->SaveAs((StoragePath+"ElectronThetaECut.pdf").c_str());
    
    
    TCanvas *c10 = new TCanvas("Electron Angle Distribution E-Cut","Electron Angular Distribution (E > 200 MeV)",1600,900);
    c10 -> SetLeftMargin(0.07);
    c10 -> SetRightMargin(0.14);
    ElectronPhiThetaECut -> GetXaxis() -> SetTitle("#varphi Angle [#circ]");
    ElectronPhiThetaECut -> GetXaxis() -> SetMaxDigits(3);
    ElectronPhiThetaECut -> GetYaxis() -> SetTitle("#theta Angle [#circ]");
    ElectronPhiThetaECut -> GetYaxis() -> SetTitleOffset(0.7);
    ElectronPhiThetaECut -> GetYaxis() -> SetMaxDigits(3);
    ElectronPhiThetaECut -> GetZaxis() -> SetTitle("e^{-} Diff. Directional Rate #frac{dr_{e}(#varphi,#theta)}{d#Omega} [sr^{-1} s^{-1}]");
    ElectronPhiThetaECut -> GetZaxis() -> SetTitleOffset(0.8);
    ElectronPhiThetaECut -> Draw("COLZ");
    c10->SaveAs((StoragePath+"ElectronPhiThetaECut.pdf").c_str());

    
    TCanvas *c11 = new TCanvas("Random","Random",1400,500);
    Random -> GetXaxis() -> SetTitle("N");
    Random -> GetYaxis() -> SetTitle("Random Number");
    Random -> Draw();
    RandomAverage -> SetLineColor(1);
    RandomAverage -> Draw("same");

    
    TCanvas *c12 = new TCanvas("Track length","Track length",700,500);
    c12 -> SetLeftMargin(0.11);
    c12 -> SetRightMargin(0.03);
    c12 -> SetLogy();
    ElectronTrackLength -> SetLineColor(46);
    ElectronTrackLength -> GetXaxis() -> SetTitle("Track length l [cm]");
    ElectronTrackLength -> GetXaxis() -> SetMaxDigits(3);
    ElectronTrackLength -> GetYaxis() -> SetTitle("Differential Rate #frac{dr_{e}(l)}{dl} [s^{-1} cm^{-1}]");
    ElectronTrackLength -> GetYaxis() -> SetTitleOffset(0.9);
    ElectronTrackLength -> Draw();
//     c12->SaveAs((StoragePath+"ElectronTrackLength.pdf").c_str());
}

