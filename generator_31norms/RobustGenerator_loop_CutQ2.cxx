#include <fstream>
#include <iostream>
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

const double Mp = 938.272046;
const double tofm = 25.68189504;
const double mu = 2.79284736; // Magnetic moment of the proton

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


//Kelly form factor of proton
const double a11_K = -0.24; // Electric form factor
const double b11_K = 10.98;
const double b12_K = 12.82;
const double b13_K = 21.97;

const double a21_K = 0.12; // Magnetic form factor
const double b21_K = 10.97;
const double b22_K = 18.86;
const double b23_K = 6.55;

Double_t GetDipoleGM(double qq)
{
    return mu*1./pow(1. + fabs(qq)/0.71, 2);
}

Double_t GetKellyGM(double qq)
{
    /*return the Kelly magnetic form factor value at Q2 (GeV^2)*/
    double t = fabs(qq) / (4. * (Mp*Mp/1.e6)); // Tau
    return mu*(1. + a21_K*t)/(1. + b21_K*t + b22_K*pow(t, 2) + b23_K*pow(t, 3));
}

Double_t GetArrington2004GM(double qq)
{
    const double p2  = 3.19;
    const double p4  = 1.355;
    const double p6  = 0.151;
    const double p8  = -1.14e-2;
    const double p10 = 5.33e-4;
    const double p12 = -9e-6; 
    return mu*(1./(1. + p2*qq + p4*qq*qq + p6*pow(qq,3) + p8*pow(qq, 4) + p10*pow(qq, 5) + p12*pow(qq, 6) ) );
}

Double_t GetArrington2007GM(double qq)
{
    double b1 = 3.173;
    double b2 = -0.314;
    double b3 = -1.165;
    double b4 = 5.619;
    double b5 = -1.087;
    double f5 = 1. + b5*qq;
    double f4 = 1. + b4*qq/f5;
    double f3 = 1. + b3*qq/f4;
    double f2 = 1. + b2*qq/f3;
    double f1 = 1. + b1*qq/f2;

    return mu/f1;
}

Double_t GetBernauerThesisGM(double qq)
{
    double GM = 1. + -2.5952*qq + 1.0222*pow(qq, 2) + 23.4945*pow(qq, 3)
                   + -93.0372*pow(qq, 4) + 140.7984*pow(qq, 5) + -0.3656*pow(qq, 6)
                   + -305.6759*pow(qq, 7) + 444.6251*pow(qq, 8) + -273.6688*pow(qq, 9)
                   + 64.5811*pow(qq, 10) ;
    return mu*GM;
}

Double_t GetXiaohuiGM(double qq)
{
    const double a21 = -1.43573; // Magnetic form factor
    const double a22 = 1.19052066;
    const double a23 = 2.5455841e-1;
    const double b21 = 9.70703681;
    const double b22 = 3.7357e-4;
    const double b23 = 6.0e-8;
    const double b24 = 9.9527277;
    const double b25 = 12.7977739;
    double t = fabs(qq)/(4.*Mp*Mp/1.e6); // Tau
    return mu*(1. + a21*t + a22*t*t + a23*t*t*t)/(1. + b21*t + b22*t*t + b23*t*t*t + b24*t*t*t*t + b25*t*t*t*t*t);
}

Double_t GetBernauerUnBoundGM2(double qq)
{
    double GM = 1. + -2.46540*qq + -0.72681*pow(qq, 2) + 35.31550*pow(qq, 3)
                   + -136.38607*pow(qq, 4) + 228.82482*pow(qq, 5) + -98.11370*pow(qq, 6)
                   +  -234.09993*pow(qq, 7) + 349.90224*pow(qq, 8) + -122.36794*pow(qq, 9)
                   + -56.48797*pow(qq, 10) + 35.79737*pow(qq, 11);
    return mu*GM;
}

Double_t GetBernauerBoundGM2(double qq)
{
    double GM = 1. + -2.76018*qq + 4.97927*pow(qq, 2) + -5.19625*pow(qq, 3)
                   + 2.19337*pow(qq, 4) + -0.00000*pow(qq, 5) + 0.50348*pow(qq, 6)
                   + -0.53301*pow(qq, 7) + 0.00000*pow(qq, 8) + -0.00001*pow(qq, 9)
                   + 0.00001*pow(qq, 10) +  -0.00000*pow(qq, 11);
    return mu*GM;
}

Double_t GetBernauerUnBoundGM1(double qq)
{
    double GM = 1. + -2.52353*qq + -0.70801*pow(qq, 2) + 40.15560*pow(qq, 3)
                   + -176.65515*pow(qq, 4) + 380.27777*pow(qq, 5) + -392.56122*pow(qq, 6)
                   + 11.52512*pow(qq, 7) + 442.38507*pow(qq, 8) +  -492.08183*pow(qq, 9)
                   + 230.28595*pow(qq, 10) + -40.91969*pow(qq, 11);
    return mu*GM;
}

Double_t GetBernauerBoundGM1(double qq)
{
    double GM = 1. + -2.79955*qq + 5.18800*pow(qq, 2) + -5.74202*pow(qq, 3)
                   + 2.80554*pow(qq, 4) + -0.00001*pow(qq, 5) + 0.01034*pow(qq, 6)
                   +  -0.27663*pow(qq, 7) + 0.00000*pow(qq, 8) + -0.00092*pow(qq, 9)
                   + 0.00127*pow(qq, 10) + -0.00003*pow(qq, 11);
    return mu*GM;
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

Double_t GetGM(double qq, int GMMODEL)
{
    if (GMMODEL == 0){
        return GetKellyGM(qq);
    }
    else if (GMMODEL == 1){
        return GetDipoleGM(qq);
    }
    else if (GMMODEL == 2){
        return GetXiaohuiGM(qq);
    }
    else if (GMMODEL == 3){
        return GetArrington2004GM(qq);
    }
    else if (GMMODEL == 4){
        return GetArrington2007GM(qq);
    }
    else if (GMMODEL == 5){
        return GetBernauerThesisGM(qq);
    }
    else if (GMMODEL == 6){
        return GetBernauerUnBoundGM1(qq);
    }
    else if (GMMODEL == 7){
        return GetBernauerBoundGM1(qq);
    }
    else if (GMMODEL == 8){
        return GetBernauerUnBoundGM2(qq);
    }
    else if (GMMODEL == 9){
        return GetBernauerBoundGM2(qq);
    }
    else return GetKellyGM(qq);
}


using namespace std;

int main()
//int main(int argc, char *argv[])
{
for (int k=0; k<9;k++){
    int FF = k;//choose different models for data generation
    //if (argc == 1) FF = stoi(argv[1]);
    
    gRandom->SetSeed(0);
   
    ifstream inFile[34];
    ifstream inf;

 
    double theta[34][2000];
    double Q2_GeV[34][2000];
    double xs[34][2000];
    double thisxs[34][2000];
    double dxs[34][2000];
    double GE[34][2000];
    double GM[34][2000];
    //double dGE[34][2000];
    int count[34];
    double temp;
    double Q2_fm[2000];
    double d_GE[2000];
    double diff[2000];
    int counts;

    
    //inf.open("../data/Carl-norm.dat"); 
    //inf>>counts;
    //for (int i=0; i<counts; i++){
    //    inf>>Q2_fm[i]>>temp>>d_GE[i];
    //}


    for (int i=0; i<32; i++){//sets
        inFile[i].open(Form("../data/cs_data/Cut_06GeV/%d.txt",i));
        inFile[i]>>count[i];
        
        for (int j=0; j<count[i]; j++){
            inFile[i]>>theta[i][j]>>Q2_GeV[i][j]>>xs[i][j]>>dxs[i][j];//loading the information from the input file
            GE[i][j] = GetGE(Q2_GeV[i][j], FF);//use the selected form factor model to calculate GE value at Q^2 (GeV^2)
//            for (int l=0; l<counts; l++){
//                diff[l] = fabs(Q2_GeV[i][j]*tofm-Q2_fm[l]);//compare Q2 points in two files
//            }
//            dGE[i][j]=d_GE[0];
//            for (int l=0; l<counts-1; l++){
//                if (diff[l+1]<diff[l]){dGE[i][j]=d_GE[l+1];}
//            }//dGE is taken from the smallest diff one
        }
        inFile[i].close();
    }
    
    for (int n=0; n<10000; n++){//repeat for X times/ generate X sets of pseudo-data
        if (n%10 == 0) cout<<n<<endl;
        int GMFF = gRandom->Integer(9);
        for (int i=0; i<32; i++){//sets
            ofstream outFile;
            //outFile.open(Form("xs_pseudo_GMfromRatio_Model8/Model%d_norm%d_table_%d.txt",k+1, i+1, n+1));//define the output file name
            outFile.open(Form("xs_pseudo_GMfromRatio_CutQ2_06GeV_1e4/Model%d_norm%d_table_%d.txt",k+1, i+1, n+1));//define the output file name
            //outFile.open(Form("Mainz_PRad_range/Model_%d/%dGeV_table_%d.txt",k+1, i+1, n+1));//define the output file name
            outFile<<count[i]<<endl;
            double norm = gRandom->Gaus(1., 0.005);//overall normalization
            //double norm = gRandom->Gaus(1., 0.002);//overall normalization
            
            for (int j=0; j<count[i]; j++){
                double Ratio = 1 - Q2_GeV[i][j]/8.;
                double GM = mu*GE[i][j]/Ratio;
                //double GM = GetGM(Q2_GeV[i][j],GMFF);//Random
                //double Q2_fm = Q2_GeV[i][j]*tofm;
                //double GM = GM=mu*(1.-8.46217e-01*Q2_fm/6.+1.07376e-01*Q2_fm)/(1.+1.07376e-01*Q2_fm);;//R11
                //double GM = GetGM(Q2_GeV[i][j],1);//dipole
                double GE_std_dip = 1/pow(1.+Q2_GeV[i][j]/0.71,2);
                double GM_std_dip = mu*GE_std_dip;
                double tau = Q2_GeV[i][j]/4./pow(Mp*1e-3, 2);
                double t = theta[i][j] / 180. * TMath::Pi();
                double epsilon = 1./(1. + 2.*(1.+tau)*pow(tan(t/2.), 2));
                thisxs[i][j] = (epsilon*pow(GE[i][j],2)+tau*pow(GM,2))/(epsilon*pow(GE_std_dip,2)+tau*pow(GM_std_dip,2));
                //double thisGE = GE[i][j]/norm + gRandom->Gaus(0., dGE[i][j]/2);
                double smeared_xs = thisxs[i][j]/norm + gRandom->Gaus(0., dxs[i][j]);//smear the stat. uncertainty to the central value of GE
                //double thisGE = GE[i][j];
                //outFile<<smearedQ2*tofm<<" "<<thisGE<<" "<<thisdGE<<endl;//write the pseudo-data to the files
                outFile<<theta[i][j]<<" "<<Q2_GeV[i][j]<<" "<<smeared_xs<<" "<<dxs[i][j]<<endl;//write the pseudo-data to the files
            }
            outFile.close();
        }
    }
}
}
