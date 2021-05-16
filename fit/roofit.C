#include<iostream>
#include<fstream>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

#include "stdlib.h"
#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TTree.h"

#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEventList.h"

//#include "data_one.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


Double_t Z_expansion(Double_t *x, Double_t *par){
      Float_t xx =x[0];
      double Tc=2.;
      double num=sqrt(Tc+xx)-sqrt(Tc);
      double den=sqrt(Tc+xx)+sqrt(Tc);
      double z=num/den;
      double GE=par[0]*(1. - par[1]*par[1]*z*2.0*2./3.+par[2]*z*z);
      return GE;
}


int roofit(){
    gStyle->SetOptFit(111111);
    //auto c=new TCanvas();
    //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+[2]*x)",0.,1.5);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x+[11]*x*x*x*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[3]*(1+[0]*x+[1]*x*x)/(1+[2]*x)",0.,45);//rational_2_1
    //TF1 *f1 = new TF1 ("f1","[0]*pow((1+x/[1]),-2)",0.,25);//dipole


    //TF1 *f1 = new TF1 ("f1","[0]/pow((1+[1]*x/12),2)",0.,26);//dipole
    //TF1 *f1 = new TF1 ("f1","[0]/(1+[1]*x/6)",0.,26);//monopole
    //TF1 *f1 = new TF1 ("f1","[0]*exp(-[1]*x/6)",0.,26);//gaussian
    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x*x)",0.,26);//poly2
    //TF1 *f1 = new TF1 ("f1","[0]*(1/(1+([1]*x/6)/(1+[2]*x)))",0.,26);//CF2
    //TF1 *f1 = new TF1 ("f1","[0]*(1/(1+([1]*x/6)/(1+([2]*x)/(1+[3]*x))))",0.,26);//CF3
    TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x*x+[3]*x*x*x)",0.,26);//poly3
    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",0.,26);//poly4
    //TF1 *f1 = new TF1 ("f1","[0]*(1/(1+([1]*x/6)/(1+([2]*x)/(1+([3]*x)/(1+[4]*x)))))",0.,26);//CF4


    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x)/(1+[2]*x)",0.,26);//rational_1_1
    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x)/(1+[2]*x+[3]*x*x)",0.,26);//rational_1_2
    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x+[3]*x*x)/(1+[2]*x)",0.,26);//rational_1_2
    //TF1 *f1 = new TF1 ("f1","[0]*(1-[1]*x/6+[2]*x+[3]*x*x)/(1+[2]*x+[4]*x*x)",0.,26);//rational_2_2
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,26);//rational_1_1
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x+[3]*x*x)",0.,26);//rational_1_2
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x+[3]*x*x)/(1+[1]*x)",0.,26);//rational_2_1
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x+[3]*x*x)/(1+[1]*x+[4]*x*x)",0.,26);//rational_2_2


    //TF1 *f1 = new TF1 ("f1","[4]*(1+[0]*x)/(1+[1]*x+[2]*x*x+[3]*x*x*x)",0.,45);//rational_1_3
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,45);//modified_A_B
    //TF1 *f1 = new TF1 ("f1","real([4]*pow((1+[0]*x),[2])/pow((1+[1]*x),[3]))",0.,45);//modified_A_B
    //TF1 *f1 = new TF1 ("f1","real([3]*pow((1+[0]*x),[2])/(1+[1]*x))",0.,45);//modified_A
    //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+0.0227*x*x+0.005*x*x*x)",0.,1.5);
    //TF1 *f1 = new TF1 ("f1",Z_expansion,0.,26,3);//polyZ
    //TF1 *f1= gROOT->GetFunction("myfunc");
    //f1->SetRange(0.0,26.)


    f1-> SetParameter(0,0.997);
    f1-> SetParameter(1,0.8404*0.8404);
    //f1-> SetParameter(2,-0.6083);
    //polyZ
    //f1-> SetParameter(2,0.1);
    //CF3
    //f1-> SetParameter(2,-0.01);
    //f1-> SetParameter(3,1.0);
    //CF4
    //f1-> SetParameter(2,-0.01);
    //f1-> SetParameter(3,1.0);
    //f1-> SetParameter(4,0.1);
    //f1-> SetParLimits(1,0.0,1.5);

    TCanvas* c1;
    c1 = new TCanvas("c1","c1",100,10,800,500);
    c1->SetGrid();
    TGraphErrors graph("../data/Carl-norm.dat","%lg %lg %lg");//input file
    //TGraphErrors graph("../generator/Mainz_PRad_range/Model_1/1GeV_table_2791.txt","%lg %lg %lg");//input file
    //TGraphErrors graph("../generator/test_input/dipole_smeared.txt","%lg %lg %lg");//input file
    //TGraphErrors graph("/home/jz271/PRad/model_generator/0.7GeV_table_onlydat.txt","%lg %lg %lg");
    graph.SetMarkerStyle(kCircle);
    graph.SetFillColor(0);
    graph.Print();
    graph.DrawClone("AP");
    graph.Fit("f1","B");
    f1->Draw("same");
    cout<<f1->GetChisquare()/(1422-3)<<endl;
    //cout<<f1->GetParError(1)<<endl;
    double p0 = f1->GetParameter(0);
    double p1 = f1->GetParameter(1);
    //double Rfit = 6*sqrt(p1*p1-p0*p0);//rational
    //double Rfit=sqrt(-6*p1);//poly
    //double Rfit=sqrt(12/p1);//dipole
    double Rfit=sqrt(p1);//dipole,monopole,gaussian
    cout<<Rfit<<endl;
    return 0;
}
