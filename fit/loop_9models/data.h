#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

//data
Int_t ndata,ndata_2;
Double_t Q2[2000],GE[2000],dGE[2000];
Double_t Q2_2[2000],GE_2[2000],dGE_2[2000];
Double_t z_arr[2000],z_arr2[2000];
Double_t fp, fp_2;
Double_t fperr, fperr_2;
Double_t pd_fit, pd1_fit, pd2_fit;
Double_t pd_fiterr, pd1_fiterr, pd2_fiterr;
Double_t chi2fit,Rfit;
Double_t Rfiterr;
Int_t fit_type;
Int_t npower;
Int_t dist_type;

int readdata(string fname, double modi_factor)
{
  TString filename=fname;

  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> ndata;

  for(int i=0;i<ndata;i++){
    infile >> Q2[i] >> GE[i] >> dGE[i];
    //if (i >= 19) { dGE[i] *= modi_factor; }
    //if (i == ndata-1) dGE[i] *= modi_factor;
    //dGE[i] *= modi_factor;
  }
  infile.close();

  //filename=fname_2;

  //infile.open(filename);
  //if(!infile.is_open()){
  //  cout << "File not exist." << endl;
  //  return 0;
  //}
  //infile >> ndata_2;

  for(int i=0;i<ndata_2;i++){
    infile >> Q2_2[i] >> GE_2[i] >> dGE_2[i];
    //if (i>=29) { dGE_2[i] *= modi_factor;}
    //if (i<15) dGE_2[i] *= modi_factor;
  }
  infile.close();

  return 0;
}
