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

#include "data.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

double GE_fun(const double &Q2, const double &R, const double &pd, const double &fp, const double &pd1)
{
  double GE;
  double R2=R*R;

  // GE=1.-Q2*R2/6.;
  //GE=fp*(1.-R2*Q2/6.+pd*Q2)/(1.+pd*Q2+pd1*Q2*Q2);//R12
  //GE=fp*(1.-R2*Q2/6.+pd*Q2+pd1*Q2*Q2)/(1.+pd*Q2);//R21
  GE=fp*(1/(1+(R2*Q2/6)/(1+pd*Q2/(1+pd1*Q2))));//CF3
  //GE=fp*(1.-R2*Q2/6.+pd*Q2*Q2+pd1*Q2*Q2*Q2);//poly3
  //GE=fp*(1.-R2*Q2/6.+pd*Q2+pd1*Q2*Q2)/(1.+pd*Q2);
  //double dipole
  //double a1 = 0.5*(R2/6-2/pd1)*(1/pd-1/pd1);
  //GE=fp*(a1*pow((1.+Q2/pd),-2)+(1-a1)*pow((1.+Q2/pd1),-2));
  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];
  double pd=fitpara[1];
  double fp1=fitpara[2];
  //double fp2=fitpara[3];
  double pd1=fitpara[3];

  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,fp1,pd1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  //for(int i=0;i<ndata_2;i++){
  //  Q2v=Q2_2[i];GE0=GE_2[i];GEd=dGE_2[i];
  //  GE_cal=GE_fun(Q2v,R,pd,fp2,pd1);

  //  result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  //}

  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000); // for Minuit
  min->SetMaxIterations(100000);  // for GSL 
  min->SetTolerance(1);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  //ROOT::Math::Functor f(&chi2,5);
  ROOT::Math::Functor f(&chi2,4);
  min->SetFunction(f);


  double step[4],variable[4];
  step[0]=0.1;step[1]=0.001;step[2]=0.1;step[3]=0.1;//step[4]=0.00001;
  variable[0]=0.87;variable[1]=9.19364e-02;variable[2]=1.0;variable[3]=10.;//variable[4]=9.19364e-02;

  //min->SetLimitedVariable(0, "Rp", variable[0], step[0], -500.0, 500.0);
  min->SetLimitedVariable(0, "Rp", variable[0], step[0], 0.0, 1.5);
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -500.0, 500.0);
  min->SetLimitedVariable(2, "fp1", variable[2], step[2], 0.9, 1.1);
  //min->SetLimitedVariable(3, "fp2", variable[3], step[3], 0.9, 1.1);
  min->SetLimitedVariable(3, "pd1", variable[3], step[3], -500.0, 500.0);

  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pd_fit=xsfit[1];
  fp=xsfit[2];
  //fp_2=xsfit[3];
  pd1_fit=xsfit[3];
  //pd2_fit=xsfit[5];

  chi2fit=f(xsfit);

  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pd_fiterr=xsfiterr[1];
  fperr=xsfiterr[2];
  //fperr_2=xsfiterr[3];
  pd1_fiterr=xsfiterr[3];
  //pd2_fiterr=xsfiterr[5];


  return 0;
}


int main(Int_t argc, char *argv[]){
 for(int i = 0; i <  1; i++){
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/Model_1/1GeV_table_%d.txt",i+1);//input files
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/Mainz_PRad_range/Model_1/1GeV_table_%d.txt",i+1);//input files
  string fname1="/var/phy/project/mepg/jz271/fitter_study/data/Carl-norm.dat";//input files
  
  if (argc == 2){
      fname1 = argv[1];
  }
  
  double modiFactor = 1.;

  readdata(fname1, modiFactor);


  double para[4];
  para[0]=8.07962e-01; para[1]=1.03484e-01; para[2]=9.99736e-01; para[3]=9.98915e-01; //para[4]=1.03484e-02; para[5]=1.03484e-02; 
  double chi2_test=chi2(para);
  cout << "chi2 test= " << chi2_test << endl;


  cout << "minimize" << endl;
  xyminimizer();

  cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
  
  ofstream outFile1;
  outFile1.open("fit_result/R12check_Model1_1e4.txt", std::ios_base::app);
  //outFile1<<Rfit<<endl;
  outFile1.close();
 

  //double fitR = sqrt(6*(2*Rfit/pd_fit+2*(1-Rfit)/pd1_fit));//double dipole
  //double fitRerr = sqrt(pow(2*(Rfiterr/pd_fit-1/pd1_fit),2)+pow(pd_fiterr*2*Rfit/(pd_fit*pd_fit),2)+pow(pd1_fiterr*2*(Rfit-1)/(pd1_fit*pd1_fit),2))/2;
  //cout<<fitR<<" "<<fitRerr<<endl; 
  }
  return 0;
}
