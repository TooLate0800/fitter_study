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

#include "data1.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


double GE_fun(const double &Q2, const double &R, const double &fp)
{
  double GE;
  double R2=R*R;
  //GE=fp/(1.+R2*Q2/12)/(1.+R2*Q2/12);//definition of dipole
  //GE=fp/(1.+R2*Q2/6);//definition of monopole
  //GE=fp*exp(-R2*Q2/6);//definition of gaussian
  //GE=fp*(1.-R2*Q2/6.);//poly_1nd
  GE=fp*(1/(1+(R2*Q2/6)));//CF1

  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];//radius, one of the fitting parameters
  double fp1=fitpara[1];//the floating parameter
  
  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
  //for(int i=1;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,fp1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);//def of chi2
  }

  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(1000000);  // for GSL 
  //min->SetTolerance(1);
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  ROOT::Math::Functor f(&chi2,2);
  min->SetFunction(f);


  double step[2],variable[2];
  step[0]=0.0001;step[1]=0.0001;
  variable[0]=0.83;variable[1]=9.99364e-01;//define the initial value of the parameters

  //define the fitting range of each parameter
  min->SetLimitedVariable(0, "Rp", variable[0], step[0],  0.0, 1.5);
  min->SetLimitedVariable(1, "fp1", variable[1], step[1], 0.9, 1.1);

  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  fp=xsfit[1];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  fperr=xsfiterr[1];


  return 0;
}


int main(Int_t argc, char *argv[]){
  for (int j = 1; j <  10; j++){
 
     for(int i = 0; i <  10000; i++){//repetition times
      string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/Mainz_PRad_range/Model_%d/1GeV_table_%d.txt",j,i+1);//input files

      
      double modiFactor = 1.;



      readdata(fname1, modiFactor);


              cout << "minimize" << endl;
              xyminimizer();

              cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
              ofstream outFile1;
              outFile1.open(Form("fit_result_PRad_range/CF1_Model%d_1e4.txt",j), std::ios_base::app);
              outFile1<<Rfit<<endl;
              outFile1.close();
      }
  }
  

  return 0;
}
