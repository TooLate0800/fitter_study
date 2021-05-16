#include <string>
#include <iostream>

void Read(){

ifstream inf("CrossSections.dat");
int k = 0;
int count[34], end_points[35];
double energy[1422], angle[1422], Q2[1422], cs[1422], dcs[1422];
double temp;
std::string spec[1422], norms[1422], norm_name[34];
for (int i=0; i<1422; i++){
    inf>>energy[i]>>spec[i]>>angle[i]>>Q2[i]>>cs[i]>>dcs[i]>>temp>>temp>>temp>>temp>>norms[i]>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
    //cout<<angle[i]<<" "<<spec[i]<<" "<<norms[i]<<endl;
}
inf.close();

count[0] = 102;
end_points[0] = 0;
end_points[34] = 1422;
int k = 0;
for (int i=0; i<1421; i++){
    if (norms[i+1] != norms[i]){
       k++;
       end_points[k] = i+1;
       //cout<<i<<" "<<norms[i]<<endl;
       norm_name[k-1] = norms[i];
       cout<<end_points[k]<<" "<<norm_name[k-1]<<endl;  
    }    
}
norm_name[33] = norms[1421];
//cout<<k<<endl;

ofstream outf("norm_table.dat");

for (int i=0; i<34; i++){
    ofstream outFile;
    count[i] = end_points[i+1]-end_points[i];
    outf<<i<<" "<<count[i]<<" "<<norm_name[i]<<endl;
    outFile.open(Form("%d.txt",i));
    //outFile.open(Form("%d_%s.txt",i, norm_name[i]));
    outFile<<count[i]<<endl;
    for (int j=end_points[i]; j<end_points[i+1]; j++){
         outFile<<angle[j]<<" "<<Q2[j]<<" "<<cs[j]<<" "<<dcs[j]<<endl;         
    }
    outFile.close();
}
outf.close();



}

