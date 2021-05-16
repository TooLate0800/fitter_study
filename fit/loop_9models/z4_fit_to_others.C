#include<iostream>
#include<fstream>
#include <iomanip>
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

void z_trans()
{
  double num,den,Tc;
  Tc=2.;
  for(int i=0;i<ndata;i++){
    num=sqrt(Tc+Q2[i])-sqrt(Tc);
    den=sqrt(Tc+Q2[i])+sqrt(Tc);
    z_arr[i]=num/den;
  }

//  for(int i=0;i<ndata_2;i++){
//    num=sqrt(Tc+Q2_2[i])-sqrt(Tc);
//    den=sqrt(Tc+Q2_2[i])+sqrt(Tc);
//    z_arr2[i]=num/den;
//  }

}

double GE_fun(const double &Q2, const double &R, const double &pd, const double& pd1, const double& pd2, const double &fp)
//double GE_fun(const double &z, const double &R, const double &pd, const double& pd1, const double& pd2, const double &fp)
{
  double GE;
  double R2=R*R;

  // GE=fp*(1.-R2*Q2/6.+pd*Q2)/(1.+pd*Q2);
  //GE=fp*(1. + R*z*+pd*z*z + pd1*z*z*z+ pd2*z*z*z*z);
  //GE=fp*(1. - R2*z*2.0*2./3.+pd*z*z + pd1*z*z*z+ pd2*z*z*z*z);
  //GE=fp*(1.-R2*Q2/6.+pd*Q2+pd1*Q2*Q2)/(1.+pd*Q2+pd2*Q2*Q2);//Rational_2_2
  //GE=fp*(1.-R2*Q2/6.+pd*Q2*Q2+pd1*Q2*Q2*Q2+pd2*Q2*Q2*Q2*Q2);//poly_4
  GE=fp*(1/(1+(R2*Q2/6)/(1+pd*Q2/(1+pd1*Q2/(1+pd2*Q2)))));//CF4




  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];
  double pd=fitpara[1];
  double pd1 = fitpara[2];
  double pd2 = fitpara[3];
  double fp1=fitpara[4];
  //double fp2=fitpara[5];
  // double fp2=fitpara[2];

  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,pd1,pd2,fp1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  //for(int i=0;i<ndata_2;i++){
  //  zv=z_arr2[i];GE0=GE_2[i];GEd=dGE_2[i];
  //  GE_cal=GE_fun(zv,R,pd,pd1,pd2,fp2);

  //  result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  //}

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
  step[0]=0.001;step[1]=0.0001;step[2]=0.0001;step[3]=0.0001;step[4]=0.0001;//step[5] = 0.0001;
  variable[0]=0.87;variable[1]=9.19364e-02;variable[2]=2.90499e+00;variable[3]=1.0;variable[4]=1.0;//variable[5] = 1.0;

  //min->SetLimitedVariable(0, "Rp", variable[0], step[0], -500., 500.);
  min->SetLimitedVariable(0, "Rp", variable[0], step[0], 0.5, 2.5);
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -500., 500.);
  min->SetLimitedVariable(2, "pd1", variable[2], step[2], -500., 500.);
  min->SetLimitedVariable(3, "pd2", variable[3], step[3], -500., 500.);
  min->SetLimitedVariable(4, "fp1", variable[4], step[4], 0.9, 1.1);
 // min->SetLimitedVariable(5, "fp2", variable[5], step[5], 0.9, 1.1);

  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pd_fit=xsfit[1];
  fp=xsfit[4];
 // fp_2=xsfit[5];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pd_fiterr=xsfiterr[1];
  fperr=xsfiterr[4];
 // fperr_2=xsfiterr[5];


  return 0;
}


int main(int argc, char *argv[]){
  /*string fname1="no_syst_1GeV_full.txt";
  string fname2="no_syst_2GeV_full.txt";
  */
  for (int j = 1; j <  10; j++){
 
     for(int i = 0; i <  10000; i++){
      string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/Mainz_PRad_range/Model_%d/1GeV_table_%d.txt",j,i+1);//input files

      readdata(fname1, 1);
      //z_trans();


      double para[5];
      para[0]=8.22168e-01; para[1]=-1.45965e+00; para[2] = 2.90499e+00; para[2]= 1.; para[3]= 9.98208e-01;para[4]=1.03484e-02; //para[5]=1.03484e-02;
      double chi2_test=chi2(para);
      cout << "chi2 test= " << chi2_test << endl;


      cout << "minimize" << endl;
      xyminimizer();

      cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
      
      ofstream outFile1;
      outFile1.open(Form("fit_result_PRad_range/z4_Model%d_1e4.txt",j), std::ios_base::app);

      //outFile1<<Rfit<<endl;

      outFile1.close();
      }
   }
  return 0;
}
