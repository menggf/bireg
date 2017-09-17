#include "util.h"
#include<iostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<algorithm> 
#include<cmath>
#include<sstream>
#include<Rcpp.h>
#include<omp.h>

// [[Rcpp::plugins(openmp)]]

using namespace std;

// [[Rcpp::export]]

Rcpp::NumericVector bcacpp(Rcpp::NumericMatrix exp, Rcpp::IntegerVector tfs, Rcpp::NumericVector par ){
    int min_gene = int(par(0)); //mininum number of genes, e.g. 50
    int max_gene = int(par(1)); //mininum number of genes, e.g. 50
    int min_patient = int(par(2)); //mininum number of patients, e.g. 50
	int step = int(par(3)); //mininum number of patients, e.g. 50
    float cutoff = float(par(4)); // the core number for multiple threads
    int cpu = int(par(5)); // the cpu number for multiple threads
    
  	//cout<<par(0)<<" "<<par(2)<<" "<<par(3) <<endl;
  	vector <float> results;
	int num_gene = exp.nrow();
	int num_patient = exp.ncol();
	
	Rcpp::CharacterVector gene_names = rownames(exp);
	Rcpp::CharacterVector patient_names = colnames(exp);
	//Rcpp::IntegerVector tfs = match(tfs_names, gene_names); 
	int num_tf= tfs.length();
	float **input= new float*[num_gene]; // record the expression data
    for(int i=0;i< num_gene;i++)
        input[i]= new float[num_patient];
    for(int i=0;i< num_gene;i++)
        for(int j=0;j< num_patient;j++)
        	input[i][j]=exp(i,j);
	
	
	const int nProcessors=cpu; // set the core number
    omp_set_num_threads(nProcessors);

#pragma omp parallel for
	//for(int xxx=0; xxx< 1; xxx++){
	for(int xxx=0; xxx< num_tf; xxx++){
		float cor[num_gene]; // record the expression data
    	int wh=int(tfs(xxx))-1; //used as core
    	//cout<<patient_names(wh)<<endl;
    	int row[num_gene];
		int col[num_patient];
		init(col, num_patient, 1);
		int left_row = 0;
		int left_col = num_patient;
		
		vector <int> pas;
		vector <int> nges;
		vector <float> ar;
		int mrow=0;
		while(left_col >= min_patient){
			cormatrix(input, cor, wh, num_gene, num_patient, col, left_col);
			int has_row=0;
			init(row, num_gene, 0);
			for(int i=0; i < num_gene; i++)
				if(cor[i] >= cutoff){
					has_row++;
					row[i]=1;
				}
			//cout<<left_col<<" "<<has_row<<endl;
			if(has_row >= max_gene || has_row <  mrow-5 )
				break;
			
			mrow=has_row;
			int ord[num_gene];
			for(int i=0;i<num_gene;i++)
				ord[i]=i;
			myorder(cor,ord,0, num_gene - 1);
			
			int used[num_gene];
			init(used, num_gene, 0);
			for(int j=0;j < min_gene;j++)
				used[ord[j]]=1;
			
			used[wh]=0;
			int removal[num_patient]; // to record wich patient to be removed
			init(removal, num_patient,-1);
			removepatient(input,removal, used,col, wh, num_patient, num_gene, min_patient, min_gene-1, cutoff); // rank the removed patients
			
			int tag=1;
			for(int ee=1; ee <= step; ee++){ // remove at most "step" patients
				init(row, num_gene, 0);
				if(left_col -1 < min_patient){ // reach the min_patients
					tag=0;
					break;
				}
				for(int j=0; j < num_patient; j++){ // remove one patients
					if(col[removal[j]]==0)
						continue;
					col[removal[j]]=0;
					left_col--;
					pas.push_back(removal[j]);
					break;
				}
				int tst=0;
				has_row=0;
				vector <float> rs;
				for(int j=0; j < num_gene; j++){ // check if statisfy min_genes
					float rr=newspearman(input, wh, ord[j], col, num_patient, left_col);
					if(rr > cutoff){
						has_row++;
						row[ord[j]]=1;
						rs.push_back(rr);
					}
					else{
						tst++;
					}
					if(tst > 20){
						break;
					}
				}
				if(rs.size()==0)
					ar.push_back(0);
				else
					ar.push_back(average(rs));
				nges.push_back(has_row);
				if(has_row >= max_gene){
					break;
				}
			}
			if(tag==0)
				break;
		}
		#pragma omp critical(datainfo)
{
		results.push_back(wh+1);
		results.push_back(pas.size());		
		for(int i=0;i<pas.size();i++)
			results.push_back(pas[i]+1);
		for(int i=0;i<nges.size();i++)
			results.push_back(nges[i]);
		for(int i=0;i<nges.size();i++)
			results.push_back(ar[i]);
}
	}
	for(int i=0;i< num_gene;i++)
        delete input[i];
    delete input;
	int nn=results.size();
	Rcpp::NumericVector zz(nn);
	for(int i=0; i < nn; i++)
		zz(i)=results[i];
    return zz;
    //return results;
}



