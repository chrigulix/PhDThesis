// Bethe-Bloch formula drawing for different materials.

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

const double C_1 = 5.2146; // delta function value for liquid argon
const double a = 0.1956; // delta function value for liquid argon
const double k = 3.0000; // delta function value for liquid argon
const double x_0 = 0.2000; // delta function value for liquid argon
const double x_1 = 3.0000; // delta function value for liquid argon

const double Z = 18; // of argon
const double A =  39.948; // g/mol of argon
const double I = 188.0e-6; // MeV of argon

const unsigned long int SampleNo = 200;
const double x_min = 1e-1;
const double x_max = 1e5;
const double y_min = 1;
const double y_max = 100;

double BetheBloch(double x, double* Parameter)
{   
    // Particle properties
    double m_part = Parameter[0];
    double z_part = Parameter[1];
    
    // Particle momentum
    double p_part = x/m_part;
    
    // Calculate different values 
    double  ArProp = std::pow(z_part,2)*Z/A;
    
    double K = 4*pi*std::pow(h_bar,2)*N_A*std::pow(alpha,2)/m_e;
    K/=(e*1e2); // From J -> MeV and m^2 -> cm^2: K = 0.397975 MeV mol^-1 cm^2.
        
    double W_max = 2*m_e_MeV*std::pow(p_part,2)/(1+2*m_e_MeV/m_part*std::sqrt(std::pow(p_part,2)+1)+std::pow(m_e_MeV/m_part,2)); // MeV
    double W_cut = 0.05;
    
    double LogTerm = std::log( 2*m_e_MeV*std::pow(p_part,2)*W_max*std::pow(I,-2) );
    
    // Bethe-Bloch without delta
    double FunctionValue = K*ArProp*( 0.5*LogTerm*(1+1/(std::pow(p_part,2))) - 1 ); //0.5*(1+W_cut/W_max)
        
    // Logarithm of particle momentum
    double p_part_log = std::log10(p_part);
    
    // calculate and incorporate delta into function
    if(x_1 <= p_part_log)
    {
        FunctionValue -= K*ArProp*(1 + 1/std::pow(p_part,2)) * 0.5*( 2*std::log(10)*p_part_log - C_1 );
    }
    else if( x_0 <= p_part_log && p_part_log < x_1)
    {
        FunctionValue -= K*ArProp*(1+1/std::pow(p_part,2)) * 0.5*( 2*std::log(10)*p_part_log - C_1 + a*std::pow(x_1-p_part_log,k) );
    }
    
    return FunctionValue*rho;
}

double BremsStrahlung(double p_part, double* Parameter)
{
    // Particle properties
    double m_part = Parameter[0];
    double E_c = Parameter[2];
    
    // Clalculate particle energy using momentum and mass
    double E_part = std::sqrt( std::pow(p_part,2) + std::pow(m_part,2) );
    
    // "Critical momentum", where Brems Strahlung and Bethe-Bloch have the same value
    double p_c = std::sqrt( std::pow(E_c,2) - std::pow(m_part,2) );
    
    // Calculate stopping power of Brems Strahlung at E_c, use Bethe-Bloch
    double StoppingPowerEc = BetheBloch(p_c,Parameter);
    
    // Determine value of Brems Strahlung at E = 0 (looking for b in y = mx + b (m = 1/X_0))
    double ZeroPoint = StoppingPowerEc - E_c/X_0;
    
    // Calculate Brems Strahlung value (y = mx + b)
    double FunctionValue = E_part/X_0;
    
    // Make exception for muon
    if(Parameter[0]>100)
    {
        // Introduce 1/X_0 form factor derived from https://cds.cern.ch/record/224555/files/cer-000137340
        FunctionValue = E_part/X_0*5.66e-5;
    }
    
    if(FunctionValue < 0.0) FunctionValue = 0.0;
    
    return FunctionValue*rho;
}

void DrawBetheBloch() 
{
    // Create vector containing relevant particle information: name, mass, charge, critical energy (e and mu only, pi estimated, rest irrelevant), and colour in plot
    std::vector<std::tuple<std::string,double, double, double,unsigned int>> Particle;
    
    // Fill vector with particle information
    Particle.push_back(std::make_tuple("e^{-}",0.510998918,-1.0,32.84,12)); // electron former colour 1
    Particle.push_back(std::make_tuple("#mu^{-}",105.6583755,-1.0,4.85e5,9)); // muon former colour 2
    Particle.push_back(std::make_tuple("#pi^{-}",139.57061,-1.0,1e10,8)); // pion former colour 3
    Particle.push_back(std::make_tuple("K^{-}",493.677,-1.0,1e10,49)); // kaon former colour 4
    Particle.push_back(std::make_tuple("p^{+}",938.27208816,1.0,1e10,46)); // proton former colour 6
    Particle.push_back(std::make_tuple("#alpha^{2+}",3727.3794066,-2.0,1e10,38)); // alpha former colour 7
    
    // Create vector containing different functions for different particles
    std::vector<TGraph*> BetheBlochGraph;
    
    // Setup Legend
    TLegend* Legend = new TLegend(0.72,0.68,0.89,0.85);
    Legend->SetNColumns(2);
    Legend->SetLineStyle ( 0 );
    Legend->SetLineColorAlpha ( 0,0 );
    Legend->SetFillStyle ( 0 );
    Legend->SetMargin ( 0.2 );
    Legend->SetEntrySeparation(0.2);
    Legend->SetHeader("Particles");
    
    // Prepare Graph fill variables
    double XValue[SampleNo], YValue[SampleNo], YBrems[SampleNo], Parameter[3], RunningX[1];
    
    // Create Log tick length
    double XAxisTic = (std::log10(x_max) - std::log10(x_min))/(double)(SampleNo-1);
    
    // Loop over different particles
    for(unsigned int i_part = 0; i_part < Particle.size(); i_part++)
    {
        // Fill Parameter
        Parameter[0] = std::get<1>(Particle.at(i_part));
        Parameter[1] = std::get<2>(Particle.at(i_part));
        Parameter[2] = std::get<3>(Particle.at(i_part));
        
        // Loop over all Graph samples
        for(unsigned long int i = 0; i < SampleNo; i++)
        {
            // Calculate log spaced x-values
            XValue[i] = std::pow(10, std::log10(x_min) + i*XAxisTic);
            
            // Use Function to calculate y-values
            YValue[i] = BetheBloch(XValue[i],Parameter);
            
            // Add Brems Strahlung to electron and muon
            if(i_part < 1)
            {
                YBrems[i] = BremsStrahlung(XValue[i],Parameter);
            }
            else if (i_part == 1) // for muon add brems strahlung
            {
                YValue[i] += BremsStrahlung(XValue[i],Parameter);
            }
            
            // "Cut" all Beta*gamma values below 0.1
            if(XValue[i]/Parameter[0] < 0.05 )
            {
                // In case of beta*gamma underflow set y-value to 10 times y-axis maximum so they are out of the plot
                YValue[i] = y_max*10;
            }
        }
        
        // If Electron 
        if(i_part < 1)
        {
            // Fill Brems
            BetheBlochGraph.push_back( new TGraph(SampleNo,XValue,YBrems) );
            BetheBlochGraph.back()->SetLineColor(std::get<4>(Particle.at(i_part)));
            BetheBlochGraph.back()->SetLineWidth(2);
            BetheBlochGraph.back()->SetLineStyle(4);
            
            // Fill Bethe
            BetheBlochGraph.push_back( new TGraph(SampleNo,XValue,YValue) );
            BetheBlochGraph.back()->SetLineColor(std::get<4>(Particle.at(i_part)));
            BetheBlochGraph.back()->SetLineWidth(2);
            BetheBlochGraph.back()->SetLineStyle(2);

            // Add Brems and Bethe-Bloch for combined plot
            for(unsigned long int i = 0; i < SampleNo; i++)
            {
                YValue[i] += YBrems[i];
            }
        }
        
        // Fill and configure graphs
        BetheBlochGraph.push_back( new TGraph(SampleNo,XValue,YValue) );
        BetheBlochGraph.back()->SetLineColor(std::get<4>(Particle.at(i_part)));
        BetheBlochGraph.back()->SetLineWidth(3);
        
        // Add legend entry
        Legend -> AddEntry( BetheBlochGraph.back(), (" "+std::get<0>(Particle.at(i_part))).c_str(),"L" );
    }
    
    // Calculate critical momentum for electron
    double ElectronP_c = std::sqrt( std::pow(std::get<3>(Particle.at(0)),2) - std::pow(std::get<1>(Particle.at(0)),2) );
    
    // Create E_c line and text
    TLine* EcLine = new TLine(ElectronP_c,y_min-0.06,ElectronP_c,rho*1.64);
    EcLine -> SetLineWidth(2);
    EcLine -> SetLineColor(1);
    EcLine -> SetLineStyle(0);
    TLatex *EcText = new TLatex(std::get<3>(Particle.at(0))-4,y_min-0.18,"p_{c}");
    
    // Make Make Texts
    TText* Ionisation = new TText(3,rho*1.18,"Ionisation");
    Ionisation -> SetTextAngle(10);
    Ionisation -> SetTextSize(0.03);
    
    TText* Radiative = new TText(std::get<3>(Particle.at(0)),rho*1.9,"Radiative");
    Radiative -> SetTextSize(0.03);
    Radiative -> SetTextAngle(63);
    
    // Create and fill multi graph
    TMultiGraph* MultiGraph = new TMultiGraph();
    
    // Fill the multi graph backwards
    for(int i = BetheBlochGraph.size()-1; i > -1; i--)
    {
        MultiGraph->Add(BetheBlochGraph.at(i));
    }
    
    // Create canvas, set it up, and draw functions
    TCanvas* C0 = new TCanvas("C0","Bethe-Bloch",0,0,1400,1000);
    C0 -> SetLogx();
    C0 -> SetLogy();
    
    MultiGraph -> SetTitle("Charged Particle Energy Dissipation");
    MultiGraph -> GetHistogram() -> GetXaxis() -> SetRangeUser(x_min,x_max);
    MultiGraph -> GetHistogram() -> GetYaxis() -> SetRangeUser(y_min,y_max);
    MultiGraph -> GetXaxis() -> SetTitle("Particle Momentum p [MeV/c]");
    MultiGraph -> GetYaxis() -> SetTitle("Linear Stop. Power  #LT-dE/dx#GT  [MeV cm^{-1}]");
    
    MultiGraph -> Draw("AC");
    EcLine -> Draw("SAME");
    EcText -> Draw("SAME");
    Ionisation -> Draw("SAME");
    Radiative -> Draw("SAME");
    Legend -> Draw("SAME");
    C0 -> SaveAs("../images/Detector/BetheBloch.pdf");
}
