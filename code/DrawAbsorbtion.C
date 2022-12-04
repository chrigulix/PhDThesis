// This macro draws the absorption coefficient of photons in liquid argon

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <string>

const double rho = 1.396; // density in g/cm^3

const double x_min = 1e-3;
const double x_max = 1e5;
const double y_min = 1e-4;
const double y_max = 1e4;

void DrawAbsorbtion()
{
    // Creat tree and fill with data
    TTree* DataTree = new TTree("Data Tree","Data Tree");
    DataTree -> ReadFile("argon.txt","Energy/F:CohScat/F:IncohScat/F:PE/F:PPn/F:PPe/F:TotCohScat/F:TotIncohScatt/F");
        
    // Initialize data holding structure for changing format
    float Energy[DataTree->GetEntries()];
    float CohScat[DataTree->GetEntries()];
    float IncohScat[DataTree->GetEntries()];
    float PE[DataTree->GetEntries()];
    float PPn[DataTree->GetEntries()];
    float PPe[DataTree->GetEntries()];
    float TotCohScat[DataTree->GetEntries()];
    float TotIncohScatt[DataTree->GetEntries()];
    
    // Temporary variable, because Root is such a fucking mess
    float E, CS, IS, P_E, PP_n, PP_e, TCS, TIS; 
    DataTree->SetBranchAddress("Energy", &E);
    DataTree->SetBranchAddress("CohScat", &CS);
    DataTree->SetBranchAddress("IncohScat", &IS);
    DataTree->SetBranchAddress("PE", &P_E);
    DataTree->SetBranchAddress("PPn", &PP_n);
    DataTree->SetBranchAddress("PPe", &PP_e);
    DataTree->SetBranchAddress("TotCohScat", &TCS);
    DataTree->SetBranchAddress("TotIncohScatt", &TIS);
    
    // Loop over tree entries
    for(int i = 0; i<DataTree->GetEntries(); i++)
    {
        // Get tree entry at index i
        DataTree->GetEntry(i);
        
        // Store tree data into arrays
        Energy[i] = E;
        CohScat[i] = CS*rho;
        IncohScat[i] = IS*rho;
        PE[i] = P_E*rho;
        PPn[i] = PP_n*rho;
        PPe[i] = PP_e*rho;
        TotCohScat[i] = TCS*rho;
        TotIncohScatt[i] = TIS*rho;
    }
    
    // Create and Fill Graphs
    TGraph* CohScatGraph = new TGraph( DataTree->GetEntries(), Energy, CohScat);
    TGraph* IncohScatGraph = new TGraph( DataTree->GetEntries(), Energy, IncohScat);
    TGraph* PEGraph = new TGraph( DataTree->GetEntries(), Energy, PE);
    TGraph* PPnGraph = new TGraph( DataTree->GetEntries(), Energy, PPn);
    TGraph* PPeGraph = new TGraph( DataTree->GetEntries(), Energy, PPe);
    TGraph* TotCohScatGraph = new TGraph( DataTree->GetEntries(), Energy, TotCohScat);
    TGraph* TotIncohScattGraph = new TGraph( DataTree->GetEntries(), Energy, TotIncohScatt);
    
    // Set graph properties SetLineStyle(0);
    CohScatGraph->SetLineColor(42);
    CohScatGraph->SetLineWidth(2);
    CohScatGraph->SetLineStyle(2);
    
    IncohScatGraph->SetLineColor(46);
    IncohScatGraph->SetLineWidth(2);
    IncohScatGraph->SetLineStyle(2);
    
    PEGraph->SetLineColor(8);
    PEGraph->SetLineWidth(2);
    PEGraph->SetLineStyle(2);
    
    PPnGraph->SetLineColor(9);
    PPnGraph->SetLineWidth(2);
    PPnGraph->SetLineStyle(2);
    
    PPeGraph->SetLineColor(38);
    PPeGraph->SetLineWidth(2);
    PPeGraph->SetLineStyle(2);
    
    TotCohScatGraph->SetLineColor(1);
    TotCohScatGraph->SetLineWidth(2);
    TotIncohScattGraph->SetLineColor(1);
    TotIncohScattGraph->SetLineWidth(2);
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.55,0.60,0.85,0.85);
    Legend->SetNColumns(1);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.3);
    Legend->SetTextSize(0.03);
//     Legend->SetHeader("Particles");
    
    // Fill legend
//     Legend -> AddEntry(CohScatGraph, " Rayleigh Scattering","L"); // coherent scattering
    Legend -> AddEntry(PEGraph, " Photoelectric Absorption","L");
    Legend -> AddEntry(IncohScatGraph, " Compton Scattering","L"); // Non-coherent scattering
    Legend -> AddEntry(PPnGraph, " Pair Production in Nuc. Field","L");
    Legend -> AddEntry(PPeGraph, " Pair Production in Elec. Field","L");
//     Legend -> AddEntry(TotCohScatGraph, " Total Attenuation","L");
    Legend -> AddEntry(TotIncohScattGraph, " Total Attenuation","L");
    
    // Creat multi graph
    TMultiGraph* MultiGraph = new TMultiGraph();
    
    // Fill multi graph
//     MultiGraph -> Add(CohScatGraph);
    MultiGraph -> Add(IncohScatGraph);
//     MultiGraph -> Add(PEGraph);
    MultiGraph -> Add(PPnGraph);
    MultiGraph -> Add(PPeGraph);
//     MultiGraph -> Add(TotCohScatGraph);
//     MultiGraph -> Add(TotIncohScattGraph);            
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Bethe-Bloch",0,0,1400,1000);
    C0 -> SetLogx();
    C0 -> SetLogy();
    
    MultiGraph -> SetTitle("Photon Attenuation Coefficient");
    MultiGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    MultiGraph -> GetXaxis() -> SetTitle("Photon Energy E [MeV]");
    MultiGraph -> GetYaxis() -> SetTitle("Attenuation Coefficient #mu [cm^{-1}]");
    MultiGraph -> Draw("AC");
    
    PEGraph -> Draw("SAME L");
    TotIncohScattGraph -> Draw("SAME L");
    
    TLatex* RadLengthTxt = new TLatex(3e3,0.15,"#mu = #frac{7}{9} #frac{#rho}{X_{0}}");
    RadLengthTxt -> Draw("SAME");
    
    TLine* RadLimit = new TLine(1e3,y_min,1e3,0.5);
    RadLimit -> SetLineWidth(2);
    RadLimit -> SetLineColor(1);
    RadLimit -> SetLineStyle(7);
    RadLimit -> Draw("SAME");
    
    Legend -> Draw();
    
    C0 -> SaveAs("../images/Detector/PhotonAbsorption.pdf");
}
