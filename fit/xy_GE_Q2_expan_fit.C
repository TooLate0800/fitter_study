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

#include "data3.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

double GE_fun(const double &Q2, const double &R, const double &pd, const double &fp)
{
  double GE;
  double R2=R*R;

  //GE=1.-Q2*R2/6.;
  //GE=fp/(1.+R2*Q2/6.);//monopole
  //GE=fp/((1.+Q2*R2/12.)*(1.+Q2*R2/12.));//dipole
  GE=fp*(1.-R2*Q2/6.+pd*Q2*Q2);

  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];
  double pd=fitpara[1];
  double fp1=fitpara[2];
  double fp2=fitpara[3];
  double fp3=fitpara[4];

  // double fp2=fitpara[2];

  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,fp1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  for(int i=0;i<ndata_2;i++){
    Q2v=Q2_2[i];GE0=GE_2[i];GEd=dGE_2[i];
    GE_cal=GE_fun(Q2v,R,pd,fp2);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  for(int i=0;i<ndata_3;i++){
    Q2v=Q2_3[i];GE0=GE_3[i];GEd=dGE_3[i];
    GE_cal=GE_fun(Q2v,R,pd,fp3);
    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  ROOT::Math::Functor f(&chi2,5);
  min->SetFunction(f);


  double step[5],variable[5];
  step[0]=0.001;step[1]=0.0001;step[2]=0.0001;step[3]=0.0001;step[4]=0.0001;
  variable[0]=0.87;variable[1]=9.19364e-02;variable[2]=1.0;variable[3]=1.0;variable[4]=1.0;

  min->SetLimitedVariable(0, "Rp", variable[0], step[0], 0.5, 2.5);
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -2., 2.);
  min->SetLimitedVariable(2, "fp1", variable[2], step[2], 0.9, 1.1);
  min->SetLimitedVariable(3, "fp2", variable[3], step[3], 0.9, 1.1);
  min->SetLimitedVariable(4, "fp3", variable[4], step[4], 0.9, 1.1);


  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pd_fit=xsfit[1];
  fp=xsfit[2];
  fp_2=xsfit[3];
  fp_3=xsfit[4];


  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pd_fiterr=xsfiterr[1];
  fperr=xsfiterr[2];
  fperr_2=xsfiterr[3];
  fperr_3=xsfiterr[4];


  return 0;
}


int main(Int_t argc, char *argv[]){
  /*string fname1="no_syst_1GeV_full.txt";
  string fname2="no_syst_2GeV_full.txt";
  */

 for(int i = 0; i <  10000; i++){

  string fname1=Form("/home/jz271/PRad/model_generator/robust_table/700_table_%d.txt",i+1);
  string fname2=Form("/home/jz271/PRad/model_generator/robust_table/1400_table_%d.txt",i+1);
  string fname3=Form("/home/jz271/PRad/model_generator/robust_table/2100_table_%d.txt",i+1);


  if (argc == 4){
      fname1 = argv[1];
      fname2 = argv[2];
      fname3 = argv[3];
  }
  
  double modiFactor = 1.;
  if (argc == 5){
      fname1 = argv[1];
      fname2 = argv[2];
      fname3 = argv[3];
      
      istringstream ss( argv[3] );
      double n;
      ss >> n;
      modiFactor += n/33;
      cout << modiFactor << endl;
  }


  readdata(fname1, fname2, fname3, modiFactor);


  double para[5];
  para[0]=1.907962e-00; para[1]=1.03484e-01; para[2]=9.99736e-01; para[3]=9.98915e-01;para[4]=9.98915e-01;
  double chi2_test=chi2(para);
  cout << "chi2 test= " << chi2_test << endl;


  cout << "minimize" << endl;
  xyminimizer();

  cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
  
  ofstream outFile1;
  outFile1.open("fitter_study/Q2_Zhihong.txt", std::ios_base::app);
  outFile1<<Rfit<<endl;
  outFile1.close();
  
}
  return 0;
}
