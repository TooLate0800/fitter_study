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

int roofit(){
    gStyle->SetOptFit(111111);
    //auto c=new TCanvas();
    //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+[2]*x)",0.,1.5);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x+[11]*x*x*x*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x)",0.,26);
    //TF1 *f1 = new TF1 ("f1","[3]*(1+[0]*x+[1]*x*x)/(1+[2]*x)",0.,45);//rational_2_1
    TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,0.5);//rational_1_1
    //TF1 *f1 = new TF1 ("f1","[4]*(1+[0]*x)/(1+[1]*x+[2]*x*x+[3]*x*x*x)",0.,45);//rational_1_3
    //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,45);//modified_A_B
    //TF1 *f1 = new TF1 ("f1","real([4]*pow((1+[0]*x),[2])/pow((1+[1]*x),[3]))",0.,45);//modified_A_B
    //TF1 *f1 = new TF1 ("f1","real([3]*pow((1+[0]*x),[2])/(1+[1]*x))",0.,45);//modified_A
    //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+0.0227*x*x+0.005*x*x*x)",0.,1.5);
    //f1-> SetParameter(0,-0.047);
    //f1-> SetParameter(1,0.6626);
    //f1-> SetParameter(2,1.0);
    //f1-> SetParameter(3,0.80);
    //f1-> SetParameter(4,1.00);
    //f1-> SetParameter(0,-0.065);
    //f1-> SetParameter(1,0.6626);
    //f1-> SetParameter(2,0.0227);
    //f1-> SetParameter(3,0.005);
    TCanvas* c1;
    c1 = new TCanvas("c1","c1",100,10,800,500);
    c1->SetGrid();
    TGraphErrors graph("/home/jz271/DRad/generator/PRad_II_10times_Kelly/2GeV_table_nosmear_onlydat.txt","%lg %lg %lg");//input file
    //TGraphErrors graph("/home/jz271/PRad/model_generator/0.7GeV_table_onlydat.txt","%lg %lg %lg");
    graph.SetMarkerStyle(kCircle);
    graph.SetFillColor(0);
    graph.Print();
    graph.DrawClone("AP");
    graph.Fit("f1","B");
    f1->Draw("same");
    cout<<f1->GetChisquare()<<endl;
    return 0;
}
