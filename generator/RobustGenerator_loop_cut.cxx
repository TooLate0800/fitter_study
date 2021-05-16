#include <fstream>
#include <iostream>
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

const double Mp = 938.272046;
const double tofm = 25.68189504;

Double_t GetExpectedEnergy( float Eb, float t)
{
    /*Given the beam energy in MeV, and scattering angle theta in deg, calculate the
      expected scattering electron energy E'*/
    Double_t theta = t/180.*TMath::Pi();
    Double_t expectE = 0.;
    
    //this formula neglect the electron mass, but compare to the one that doesn't
    //the value of (1 - ratio between the two) is at the level of ~1e-9
    expectE = Eb*Mp / ( Eb*(1.-cos(theta)) + Mp );
   
    return expectE;
}

Double_t GetQ2(double theta, double EBeam)
{   /*Given the scattering angle theta in deg, beam energy in MeV, calculate the four momentum transfer squared
      in GeV^2, electron mass is not neglected in the formula*/
    return (2.*EBeam*GetExpectedEnergy(EBeam, theta)*(1-cos(theta/180.*TMath::Pi()))/1.e6 /*- 2.*pow(Me/1.e3, 2)*/ );
}

double z_calc(double Q2_in, double Tc_in, double T0_in)
{
  double res=0.;
  double num=sqrt(Tc_in+Q2_in)-sqrt(Tc_in-T0_in);
  double den=sqrt(Tc_in+Q2_in)+sqrt(Tc_in-T0_in);
  res=num/den;

  return res;
}



double GetDipole(double qq)
{
    qq *= tofm;
    double p0 = 12./(0.83*0.83);
    
    return 1./(1.+qq/p0)/(1.+qq/p0);
}

double GetMonopole(double qq)
{
    qq *= tofm;
    double p0=6./(0.83*0.83);
    
    return 1./(1.+qq/p0);
}

double GetGaussian(double qq)
{
    qq *= tofm;
    double p0=6./(0.83*0.83);
    
    return exp(-qq/p0);
}

double GetKellyGE(double x)
{
    x *= tofm;
    double A=1./(tofm*4.*(Mp/1000.)*(Mp/1000.));
    
    double p0 = -0.24*A;
    double p1 = 10.98*A;
    double p2 = 12.82*A*A;
    double p3 = 21.97*A*A*A;
    
    return ((1.+p0*x)/(1.+p1*x+p2*x*x+p3*x*x*x));
}

double GetArrington1(double x)//Venkat_2011
{
    x *= tofm;
    double A=1./(tofm*4.*(Mp/1000.)*(Mp/1000.));
    
    double p0=2.90966*A;
    double p1=-1.11542229*A*A;
    double p2=3.866171e-2*A*A*A;
    double p3=14.5187212*A;
    double p4=40.88333*A*A;
    double p5=99.999998*A*A*A;
    double p6=4.579e-5*A*A*A*A;
    double p7=10.3580447*A*A*A*A*A;
    
    return ((1.+p0*x+p1*x*x+p2*x*x*x)/(1.+p3*x+p4*x*x+p5*x*x*x+p6*x*x*x*x+p7*x*x*x*x*x));
}

double GetArrington2(double x)//Arrington_2004
{
    x *= tofm;
    double A=1./tofm;
    double p0=3.226*A;
    double p1=1.508*A*A;
    double p2=-0.3773*A*A*A;
    double p3=0.611*A*A*A*A;
    double p4=-0.1853*A*A*A*A*A;
    double p5=1.596e-2*A*A*A*A*A*A;
    
    return (1./(1.+p0*x+p1*x*x+p2*x*x*x+p3*x*x*x*x+p4*x*x*x*x*x+p5*x*x*x*x*x*x));
}

double GetArringtonSick(double x)//Arrington_2007
{   
    double f5 = 1. + -0.284*x;
    double f4 = 1. +  1.176*x/f5;
    double f3 = 1. + -1.212*x/f4;
    double f2 = 1. + -0.178*x/f3;
    double f1 = 1. +  3.440*x/f2;
    
    return 1./f1;
}

double GetAlarcon(double x)
{
    x*= tofm;
    double tz=z_calc(x, 2.00252188772, -17.99);
    
    double GE = 0.239448638275 - 1.11264843899*tz + 1.44819766508*pow(tz,2) + 0.514648365159*pow(tz,3) - 2.36672103495*pow(tz,4) + 0.926165750483*pow(tz,5) + 2.05049172945*pow(tz,6) - 4.11073019989*pow(tz,7) + 2.64932410946*pow(tz,8) + 3.51485719222*pow(tz,9) - 7.5760640139*pow(tz,10) + 4.96350589461*pow(tz,11) - 1.14047565701*pow(tz,12);  //0.85 fm
    
    return GE;
}

double GetBernauer2014(double x)
{
    double tz = x;
    double GE =1. - 3.36591660*tz + 1.45487683e+01*pow(tz,2) - 8.87959239e+01*pow(tz,3) + 4.61097705e+02*pow(tz,4) - 1.67562381e+03*pow(tz,5) + 4.07646487e+03*pow(tz,6) - 6.45411460e+03*pow(tz,7) + 6.34035079e+03*pow(tz,8) - 3.49373923e+03*pow(tz,9) + 8.22601568e+02*pow(tz,10);  //0.8868 fm
    return GE;
}

double GetZhihong(double x)
{
    x *= tofm;
    double tz=z_calc(x, 2.00252188772, -17.99);
    double GE = 0.239163298067 - 1.109858574410*tz + 1.444380813060*pow(tz,2) + 0.479569465603*pow(tz,3) - 2.286894741870*pow(tz,4) + 1.126632984980*pow(tz,5) + 1.250619843540*pow(tz,6) - 3.631020471590*pow(tz,7) + 4.082217023790*pow(tz,8) + 0.504097346499*pow(tz,9) - 5.085120460510*pow(tz,10) + 3.967742543950*pow(tz,11) - 0.981529071103*pow(tz,12);  //0.8792 fm
    return GE;
}

double GetBernauerThesis(double x)
{
    double tz = x;
    double GE =1. - 3.3686*tz + 14.5606*pow(tz,2) - 88.1912*pow(tz,3) + 453.6244*pow(tz,4) - 1638.7911*pow(tz,5) + 3980.7174*pow(tz,6) - 6312.6333*pow(tz,7) + 6222.3646*pow(tz,8) - 3443.2251*pow(tz,9) + 814.4112*pow(tz,10);  //0.8872 fm
    return GE;
}


double GetGE(double qq, int FF)
{
    if (FF == 0){
        //dipole form factor -- 0.83
        return GetDipole(qq);
    }
    else if (FF == 1){
        //monopole form factor -- 0.83
        return GetMonopole(qq);
    }
    else if (FF == 2){
        //gaussian form factor -- 0.83
        return GetGaussian(qq);
    }
    else if (FF == 3){
        //Kelly form factor -- 0.863
        return GetKellyGE(qq);
    }
    //else if (FF == 4){
    //    //Arrington form factor 1 -- 0.8779
    //    return GetArrington1(qq);
    //}
    else if (FF == 4){
        //Arrington form factor 2 -- 0.8682
        return GetArrington2(qq);
    }
    else if (FF == 5){
        //Arrington sick -- 0.8965
        return GetArringtonSick(qq);
    }
    else if (FF == 6){
        //Alarcon -- 0.85
        return GetAlarcon(qq);
    }
    else if (FF == 7){
        //Bernauer 2014 -- 0.8868
        return GetBernauer2014(qq);
    }
    else if (FF == 8){
        //Zhihong -- 0.8792
        return GetZhihong(qq);
    }
    //else if (FF == 10){
    //    //Bernauer thesis -- 0.8872
    //    return GetBernauerThesis(qq);
    //}
}

using namespace std;

int main()
//int main(int argc, char *argv[])
{
for (int k=0; k<9;k++){
    int FF = k;//choose different models for data generation
    //if (argc == 1) FF = stoi(argv[1]);
    
    gRandom->SetSeed(0);
   
    ifstream inFile[3];
    inFile[0].open("/var/phy/project/mepg/jz271/fitter_study/data/Carl-norm.dat");//the input file, need to obtain Q^2 points and the stat uncer of GE
    //inFile[0].open("/var/phy/project/mepg/jz271/fitter_study/data/Mainz_RS.dat");//the input file, need to obtain Q^2 points and the stat uncer of GE
    //inFile[0].open("ge_uncertainty_table/1GeV_GE.txt");
    //inFile[1].open("ge_uncertainty_table/2GeV_GE.txt");
//    inFile[2].open("ge_uncertainty_table/3.3GeV_ge_uncertainty_2times.txt");
    
    double EBeam[3][100];
    double theta[3][100];
    double error[3][100];
    double Q2[3][2000];
    double GE[3][2000];
    double dGE[3][2000];
    int count[3];
    int count_cut = 0;   
 
    for (int i=0; i<1; i++){
    //for (int i=0; i<3; i++){
        inFile[i]>>count[i];
        
        for (int j=0; j<count[i]; j++){
            inFile[i]>>Q2[i][j]>>GE[i][j]>>dGE[i][j];//loading the information from the input file
            Q2[i][j] = Q2[i][j]/tofm;//pay attention to the dimension of Q^2 (default:fm^-2)
            GE[i][j] = GetGE(Q2[i][j], FF);//use the selected form factor model to calculate GE value at Q^2 (GeV^2)
            if (Q2[i][j]<0.25){
                count_cut++;
            }
        }
        inFile[i].close();
    }
    
    for (int n=0; n<10000; n++){//repeat for X times/ generate X sets of pseudo-data
    //for (int n=0; n<10000; n++){//repeat for X times/ generate X sets of pseudo-data
        if (n%1000 == 0) cout<<n<<endl;
        for (int i=0; i<1; i++){
            ofstream outFile;
            //outFile.open("dipole_smeared.txt");//define the output file name
            //outFile.open(Form("Mainz_PRad_range_nosmear/Model_%d_table.txt",k+1));//define the output file name
            //outFile.open(Form("PRad_pseudo/Model_%d/%dGeV_table_%d.txt",k+1, i+1, n+1));//define the output file name
            outFile.open(Form("Mainz_pseudo_Cut025GeV/Model_%d_table_%d.txt",k+1, n+1));//define the output file name
            outFile<<count_cut<<endl;
            double norm = gRandom->Gaus(1., 0.002);//overall normalization
            
            for (int j=0; j<count_cut; j++){
                //cout<<Q2[i][j]*tofm<<" "<<GE[i][j]<<endl;
                //double thisGE = GE[i][j]/norm + gRandom->Gaus(0., dGE[i][j]/2);
                double thisGE = GE[i][j]/norm + gRandom->Gaus(0., dGE[i][j]);//smear the stat. uncertainty to the central value of GE
                //double thisGE = GE[i][j];
                double thisdGE = dGE[i][j];
                //outFile<<smearedQ2*tofm<<" "<<thisGE<<" "<<thisdGE<<endl;//write the pseudo-data to the files
                if (Q2[i][j]<0.25){
                    outFile<<Q2[i][j]*tofm<<" "<<thisGE<<" "<<dGE[i][j]<<endl;//write the pseudo-data to the files
                }
            }
            outFile.close();
        }
    }
}
}
