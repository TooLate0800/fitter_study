#include <algorithm>
#include <string>
#include <iostream>

void Read_cut(){

ifstream inf("CrossSections.dat");
int k = 0;
int count[34], end_points[35];
double energy[1422], angle[1422], Q2[1422], cs[1422], dcs[1422];
double temp;
std::string spec[1422], norms[1422], norm_name[34];
int j = 0;
double en,an,QQ,c,dc;
std::string sp, nor;
for (int i=0; i<1422; i++){
    inf>>en>>sp>>an>>QQ>>c>>dc>>temp>>temp>>temp>>temp>>nor>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
    //cout<<nor<<endl;
    if (QQ < 0.7){
        energy[j] = en;
        spec[j] = sp;
        angle[j] = an;
        Q2[j] = QQ;
        cs[j] = c;
        dcs[j] = dc;
        norms[j] = nor;
        j++;
    }
    //cout<<angle[i]<<" "<<spec[i]<<" "<<norms[i]<<endl;
}
inf.close();
int total_cout = j;
cout<<total_cout<<endl;

//count[0] = 102;
end_points[0] = 0;
char mark = ":";
//end_points[34] = total_cout+1;
int k = 0;
for (int i=0; i<total_cout-1; i++){
    if (norms[i+1] != norms[i]){
       k++;
       end_points[k] = i+1;
       //cout<<i<<" "<<norms[i]<<endl;
       std::replace(norms[i].begin(), norms[i].end(), ':', ' ');
       norm_name[k-1] = norms[i];
       //cout<<end_points[k]<<" "<<norm_name[k-1]<<" "<<norm_name[k-1].length()<<endl;  
    }    
}
end_points[k+1] = total_cout;
std::replace(norms[total_cout-1].begin(), norms[i].end(), ':', ' ');
norm_name[k] = norms[total_cout-1];
//cout<<norm_name[k]<<endl;
cout<<k+1<<endl;//number of data set

ofstream outf("Cut_07GeV/norm_table.dat");
outf<<k+1<<endl;

for (int i=0; i<k+1; i++){
    ofstream outFile;
    count[i] = end_points[i+1]-end_points[i];
    if (norm_name[i].length() < 3){
    outf<<i<<" "<<count[i]<<" "<<norm_name[i]<<" "<<"-1"<<endl;}
    else{
    outf<<i<<" "<<count[i]<<" "<<norm_name[i]<<endl;}
    outFile.open(Form("Cut_07GeV/%d.txt",i));
    //outFile.open(Form("%d_%s.txt",i, norm_name[i]));
    outFile<<count[i]<<endl;
    for (int j=end_points[i]; j<end_points[i+1]; j++){
         outFile<<angle[j]<<" "<<Q2[j]<<" "<<cs[j]<<" "<<dcs[j]<<endl;         
    }
    outFile.close();
}
outf.close();



}

