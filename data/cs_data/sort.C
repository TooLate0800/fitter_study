#include <algorithm>
#include <string>
#include <iostream>

void count_norm(){

ifstream inf("Cut_07GeV/norm_table.dat");
Int_t number_of_sets;
Int_t number_of_norms = 0;
Int_t number_of_points[40];
Int_t norm1[40];
Int_t norm2[40];
Int_t norm[80];
double temp;
  inf >> number_of_sets;
  cout<<number_of_sets<<endl;
  for(int i=0;i<number_of_sets;i++){
    inf >> temp >> number_of_points[i] >> norm1[i] >> norm2[i];
    //cout<<norm1[i]<<" "<<norm2[i];
  }
  inf.close();
  number_of_norms = number_of_sets*2;
for(int i=0;i<number_of_sets;i++){
    for(int j=0;j<i;j++){
        if (norm1[i] == norm1[j] || norm1[i] == norm2[j]){number_of_norms--;break;}//neglect repeated numbers
    }
    for(int j=0;j<i;j++){
        if (norm2[i] == norm2[j] || norm2[i] == norm1[j]){number_of_norms--;break;}
    }
    norm[i] = norm1[i];
    norm[i+number_of_sets] = norm2[i];//combine into 1 array
}
number_of_norms--;//neglect -1;
cout<<number_of_norms<<endl;

cout<<"====================="<<endl;
sort(norm,norm+2*number_of_sets);
for(int i=0;i<2*number_of_sets;i++){
    cout<<norm[i]<<endl;
}
cout<<"====================="<<endl;
for (int j=0; j<number_of_sets; j++){
    cout<<norm1[j]<<" "<<norm2[j]<<endl; 
}

cout<<"====================="<<endl;
int new_norm = 0;
int start = 0;
for(int i=1;i<2*number_of_sets;i++){
    if (norm[i-1] ! = -1 && norm[i] > (norm[i-1]+1)){
        new_norm = norm[i-1] + 1;
        start = i+1;
        //cout<<i<<" "<<norm[i]<<" "<<new_norm<<endl;
        for (int j=0; j<number_of_sets; j++){
            if (norm1[j] == norm[i]){norm1[j]=new_norm;}
            if (norm2[j] == norm[i]){norm2[j]=new_norm;}
        }
        break; 
    }
}
cout<<start<<endl;
//cout<<"====================="<<endl;
if (start != 0){
    for(int i=start;i<2*number_of_sets;i++){
        if (norm[i-1] ! = -1 && norm[i] > (norm[i-1])){
            new_norm++;//wrong
            //cout<<i<<" "<<norm[i]<<" "<<new_norm<<endl;
            for (int j=0; j<number_of_sets; j++){
                if (norm1[j] == norm[i]){norm1[j]=new_norm;}
                if (norm2[j] == norm[i]){norm2[j]=new_norm;}
            }
        }
    }
}
cout<<"====================="<<endl;
for (int j=0; j<number_of_sets; j++){
    cout<<norm1[j]<<" "<<norm2[j]<<endl; 
}

}

