#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

//data
Int_t ndata,ndata_2,ndata_3;
Double_t Q2[2000],GE[2000],dGE[2000];
Double_t Q2_2[1000],GE_2[1000],dGE_2[1000];
Double_t Q2_3[1000],GE_3[1000],dGE_3[1000];
Double_t z_arr[1000],z_arr2[1000], z_arr3[1000];
Double_t fp, fp_2, fp_3;
Double_t fperr, fperr_2, fperr_3;
Double_t pd_fit;
Double_t pd_fiterr;
Double_t chi2fit,Rfit;
Double_t Rfiterr;
Int_t fit_type;
Int_t npower;
Int_t dist_type;

int readdata(string fname, string fname_2, string fname_3, double modi_factor)
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
  }
  infile.close();

  filename=fname_2;

  infile.open(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> ndata_2;

  for(int i=0;i<ndata_2;i++){
    infile >> Q2_2[i] >> GE_2[i] >> dGE_2[i];
  }
  infile.close();

  filename=fname_3;

  infile.open(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> ndata_3;

  for(int i=0;i<ndata_3;i++){
    infile >> Q2_3[i] >> GE_3[i] >> dGE_3[i];
  }
  infile.close();

  return 0;
}
