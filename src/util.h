#include<cmath>
#include<iostream>
#include<vector>
#include<map>
#include<fstream>
#include<sstream>
#include<cstdlib>

using namespace std;

template <class T>
int max_pos(T *x, int n){
	int i=0;
	int wh=0;
	for(i=1;i<n;i++)
		if(x[i]>x[wh])
			wh=i;
	return wh;
}
template <class T>
int min_pos(T *x, int n){
	int i=0;
	int wh=0;
	for(i=1;i<n;i++)
		if(x[i] < x[wh])
			wh=i;
	return wh;
}

string to_string (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}
vector <string> split(string s, const string delim) {
    vector<string> result;
    size_t pos=0;
    
    while(1){
        size_t temp=s.find(delim, pos);
        if(temp==string::npos){
            result.push_back(s.substr(pos,1000));
            return result;
        }
        if(temp==pos)
            result.push_back("");
        else
            result.push_back(s.substr(pos,temp-pos));
        pos=temp+1;
    }
    return result;
}

template <class T>
void init(T *x,int n,T v){
    for(int i=0;i< n;i++)
        x[i]=v;
}
void merge(float *a, int *rank, int r,int m,int l)
{
	int i,j,k;
	int n1=m-r+1,n2=l-m;
	int *b=new int[n1];
	int *c=new int[n2];
	for (i = 0; i<n1; i++) {
		b[i]=rank[i+r];
	}
	for (i = 0; i<n2; i++) {
		c[i]=rank[i+m+1];
	}
	i=j=0;
	k=r;
	while (i<n1&&j<n2)
	{
		if(a[b[i]] < a[c[j]]){
			rank[k]=c[j];
			j++;
		}
		else{
			rank[k]=b[i];
			i++;
		}
		k++;
	}
	while (i<n1){
		rank[k]=b[i];
		i++;
		k++;
	}
	while (j<n2){
		rank[k]=c[j];
		k++;
		j++;
	}
	delete b;
	delete c;
}
void merge(int *a, int *rank, int r,int m,int l)
{
	int i,j,k;
	int n1=m-r+1,n2=l-m;
	int *b=new int[n1];
	int *c=new int[n2];
	for (i = 0; i<n1; i++) {
		b[i]=rank[i+r];
	}
	for (i = 0; i<n2; i++) {
		c[i]=rank[i+m+1];
	}
	i=j=0;
	k=r;
	while (i<n1&&j<n2)
	{
		if(a[b[i]] < a[c[j]]){
			rank[k]=c[j];
			j++;
		}
		else{
			rank[k]=b[i];
			i++;
		}
		k++;
	}
	while (i<n1){
		rank[k]=b[i];
		i++;
		k++;
	}
	while (j<n2){
		rank[k]=c[j];
		k++;
		j++;
	}
	delete b;
	delete c;
}
template <class T>
void myorder(T *a,int *rank, int r,int l)
{
	if (r<l) {
		int m=(r+l)/2;
		myorder(a,rank, r,m);
		myorder(a,rank, m+1,l);
		merge(a,rank, r,m, l);
	}
}

template <class T>
void myorder2(T *cc, int *od,int n){ // get the rank
	int *temp=new int[n];
	int i;
	for(i=0;i<n;i++)
		temp[i]=i;
    myorder(cc, temp, 0, n-1);
    for(i=0;i<n;i++)
    	od[temp[i]]=i;
    delete temp;
}
float pick(float *x,int n,int wh){
    int temp[n],i;
    for(i=0;i<n;i++)
        temp[i]=i;
    myorder(x,temp,0, n-1);
    return x[temp[wh-1]];    
}
float spearman(float **input, int from, int to, int num){ // spearman's correlation
	int order1[num], order2[num],i;
	for(i=0;i<num;i++){
		order1[i]=i;
		order2[i]=i;
	}
	myorder2(input[from], order1, num);
	myorder2(input[to], order2, num);
	float d=0.0;
	for(int i=0;i<num;i++)
		d =d + pow(order1[i]-order2[i],2);
	
	return(1 - 6 * d/(num*(pow(num,2)-1)));
}
float newspearman(float **input, int from, int to, int *col, int num_patient, int num){ // spearman's correlation by allow exclude some patients
	int i,tt=0;
	float *x=new float[num];
	float *y=new float[num];
	for(i=0;i<num_patient;i++){
		if(col[i]==0)
			continue;
		x[tt]=input[from][i];
		y[tt]=input[to][i];
		tt++;
	}
	if(tt!=num){
		std::cout<<"error!\n";
		std::cout<<tt<<"\t"<<num<<"\n";
	}
	int order1[num], order2[num];
	for(i=0; i<num; i++){
		order1[i]=i;
		order2[i]=i;
	}
	myorder2(x, order1, num);
	myorder2(y, order2, num);
	delete x;
	delete y;
	float d=0.0;
	for(int i=0;i<num;i++)
		d=d + pow(order1[i]-order2[i],2);
	return(1- 6*d/(num*(pow(num,2)-1)));
}

int myremoval(int target,int partner,int ncol, float **array, int min_patient,float cutoff,int *result, int info=0){ 
    int order1[ncol], order2[ncol],order3[ncol],i=0;
    for(i=0;i<ncol;i++){
		order1[i]=i;
		order2[i]=i;
		order3[i]=i;
	}
    myorder2(array[target], order1, ncol); // get the rank
    myorder2(array[partner], order2, ncol); // get the rank
    //for(int i=0;i<10;i++)
    //	cout<<"\t"<<i<<"\t"<<order1[i]<<"\t"<<array[target][i]<<endl;
    
    float d[ncol];
    for(i=0;i<ncol;i++)
    	d[i]=pow(order1[i]-order2[i],2);
    myorder(d, order3, 0, ncol-1); // order the rank difference
    
    int remove=0;
    for(i=0;i<ncol;i++){
    	float r=newspearman(array, target,partner,result, ncol, ncol-remove);
    	//cout<<i<<"\t"<<r<<"\t"<< d[order3[i]]<<"\t"<<order3[i]<<endl;
    	if(r > cutoff)
    		break;
    	result[order3[i]]=0;
    	remove++;
    	if(ncol-remove < min_patient)
    		return 0;
    }
    for(i=0;i<ncol;i++){
    	if(i >= remove)
    	 	result[i]=-1; // not remove
    	 else
    	 	result[i]=order3[i]; // here small number should be remove
    }
    return ncol-remove;
}

void vote(short **cc,int *order,int ncol,int n,int window=50){ //use steps to find the, here "n" is the gene numbers
    int i,j,k;
    int has[ncol]; // make sure not to count the used ones;
    init(has,ncol,0);
    int ss[ncol]; //count the occurrence
    int rr=0;
    int temp[ncol];
    for(i=window;i < ncol;i=i+window){
        init(ss, ncol, 0);
        for(j=0;j<i;j++){
            for(k=0;k< n;k++){
                if(cc[k][j] == -1) // make sure not -1
                    continue;
                ss[cc[k][j]]++; // count the rank occurence in the windows
            }
        }
        for(j=0;j<ncol;j++)
            temp[j]=j;
        myorder(ss,temp,0,ncol-1);
        int tt=0;
        for(j=0;j<ncol;j++){
            if(has[temp[j]]==1)
                continue;
            if(ss[temp[j]]<=((float(n)/10) > 2 ? (float(n)/10):2  ))
                break;
            order[rr]=temp[j];
            has[temp[j]]=1;
            tt++;
            rr++;
            if(rr==ncol-1)
                break;
            if(tt==window)
                break;
        }
    }
    for(i=0;i<ncol;i++){
        if(has[temp[i]]==1)
            continue;
        order[rr]=temp[i];
        rr++;
        if(rr==ncol)
            break;
    }
}

int connectivity(short **mx, int *remove, int num, float cutoff, int min_gene){
	init(remove,num,0);
	int have=num;
	int sum[num];
	for(int i=0;i<num; i++){
		sum[i]=0;
		for(int j=0;j<num;j++)
			sum[i]+=mx[i][j];
	}
	while(have >= min_gene){
		int wh=min_pos(sum, num);
		if(sum[wh] >= have*cutoff)
			return have;
		remove[wh]=1;
		have--;
		for(int i=0;i<num;i++)
			if(remove[i]==1)
				sum[i]=num;
			else
				sum[i]-=mx[wh][i];
	}
	return -1;
}
void iscore(float **sim, float *row, int num_gene, int min_gene,  vector <int> &c1, float mycutoff){
	for(int i=0;i < num_gene;i++){
    	if(row[i] > 0.0 ) // not to count the used ones
      		continue;
      	int ccc1=0; //counts 
       	for(int j=0;j< num_gene;j++){
       		if(sim[i][j] < mycutoff || row[j]> 0.0)
       			continue;
           	ccc1++;
        }
        if(ccc1 >= min_gene)
        	c1.push_back(i);
    }
}
vector <int> writeout(int left_row, int *row, int left_col, int *col, int num_gene, int num_patient, int wh){
	vector <int> result;
	result.push_back(wh);
	result.push_back(left_row);
	int pp=0;
	for(int i=0; i< num_gene; i++){
		if(row[i]==0)
			continue;
		pp++;
		result.push_back(i);
	}
	if(pp!=left_row){
		cout<<"Error: row : "<< pp <<" "<< left_row <<endl;
	}
	pp=0;
	result.push_back(left_col);
	for(int i=0; i< num_patient; i++){
		if(col[i]==0)
			continue;
		pp++;
		result.push_back(i);
	}
	if(pp!=left_col){
		cout<<"Error: col: "<< pp <<" "<< left_col <<endl;
	}
	return result;
}
int limitto(int *has_tag,float *has, int num_gene, int limited_num){
	float ct=pick(has,num_gene,limited_num);
	int x=0;
	for(int j=0;j<num_gene;j++)
    	if(has[j] < ct )
        	has_tag[j]=0;
       	else
        	x++;	
    return x; 
}

void removepatient(float **input,int *removal, int *has_tag, int wh, int num_patient, int num_gene, int min_patient,int has_used_num,float cutoff ){
	short **record= new short*[has_used_num]; //a matrix to record the removal rank
	for(int j=0;j < has_used_num;j++)
    	record[j]=new short[num_patient];
    int zz=0;
    for(int j=0;j<num_gene;j++){ // to count or vote, use every gene as partner
    	if(has_tag[j]==0) // not similar
        	continue;
		int *use=new int[num_patient];
		init(use, num_patient,1);
		int out= myremoval(wh,j,num_patient, input, min_patient,cutoff, use);
		if(out!=0){
			for(int k=0;k<num_patient;k++)
		        record[zz][k]=use[k]; 
		}	
		else{
		    for(int k=0;k<num_patient;k++)
		        record[zz][k]=-1;	 
		}
		delete use;
        zz++;
	}
    int step=30;
    vote(record, removal,num_patient,has_used_num,step); // here we use a  widows of step to find which ones are very like to remove and rank them   
    for(int j=0;j<has_used_num;j++)
    	delete record[j];
    delete record;         
}
void removepatient(float **input,int *removal, int *has_tag, int *col, int wh, int num_patient, int num_gene, int min_patient,int has_used_num,float cutoff ){
	short **record= new short*[has_used_num]; //a matrix to record the removal rank
	for(int j=0;j < has_used_num;j++)
    	record[j]=new short[num_patient];
    int zz=0;
    for(int j=0;j<num_gene;j++){ // to count or vote, use every gene as partner
    	if(has_tag[j]==0) // not similar
        	continue;
		int *use=new int[num_patient];
		init(use, num_patient,1);
		int out= myremoval(wh,j,num_patient, input, min_patient,cutoff, use);
		if(out!=0){
			for(int k=0;k<num_patient;k++)
		        if(col[k]!=0)
		        	record[zz][k]=use[k];
		        else
		        	record[zz][k]=-1;
		}	
		else{
		    for(int k=0;k<num_patient;k++)
		        record[zz][k]=-1;	 
		}
		delete use;
        zz++;
	}
    int step=50;
    vote(record, removal,num_patient,has_used_num,step); // here we use a  widows of step to find which ones are very like to remove and rank them   
    for(int j=0;j<has_used_num;j++)
    	delete record[j];
    delete record;         
}

int refresh(float **input, int *removal, int *has, int *recordcol, int wh, int num_patient, int num_gene, int min_patient, int min_gene, float cutoff){
	int remove=0; // the number of removed patient
    int max_have=0; // 
    int recordremove=-1; // the patinent number of removed record
    int include_num=0; // the number of included genes
    int col[num_patient];
            init(col, num_patient, 1);
            for(int k=0; k < num_patient; k++){ //vote to determined which will be remove
            	if( num_patient - remove < min_patient)
                    return recordremove;
                include_num=0;
                if(removal[k] == -1)
                    continue;
                for(int g=0; g<num_gene; g++){ // count how many patients has good similar after remove one patients
                    if(has[g]==0)
                        continue;
                    float r=newspearman(input, wh, g, col, num_patient, num_patient-remove);// the correlation after remove one patient
                    if(r > cutoff)
                        include_num++;
                }
                if(recordremove != -1 && remove > 0.4 * num_patient)
                	return recordremove;
                if(include_num > min_gene && include_num * ( num_patient - remove ) > max_have + 1){ // we will record the information with max coverage
                   max_have=include_num * (num_patient-remove);
                   recordremove = remove;
                   for(int g=0; g < num_patient; g++){
                       if(col[g] == 0){
                            recordcol[g]=0;
                        }
                        else{
                            recordcol[g]=1;
                        }
                   }
                   if(remove > 0.4 * num_patient)
                   		break;
                }
                remove++;
                col[removal[k]]=0;
            }
            return recordremove;
}
int compact(float **sim, int *has, int *recordcol, int wh, int recordremove, int num_gene, int num_patient, int min_gene, float overlap, float cutoff ){
    int include_num=0;
    map <int,int> tt; 
	for(int j=0;j<num_gene;j++){
		if(has[j]==0)
			continue;
		tt[j]=include_num;
		include_num++; // include itself, only for the un-used gene
	}
	short **mx=new short*[include_num]; // a matrix indicate wether the genes have enough similarity
	for(int j=0;j<include_num;j++){
	   	mx[j]=new short[include_num];
	   	for(int k=0;k<include_num;k++)
			mx[j][k]=0;	
	}        
	for(int x=0; x < num_gene;x++){
		if(has[x]==0)
			continue;
		for(int y=x;y<num_gene;y++){
			if(has[y]==0)
		    	continue;
			if(sim[x][y] < cutoff)
				continue;
			mx[tt[x]][tt[y]]=1;	
			mx[tt[y]][tt[x]]=1;	
		}
	}	
	int rm[include_num];
	init(rm, include_num, 0);
	int new_include_num=connectivity(mx, rm, include_num, overlap, min_gene);
	for(int j=0;j<include_num;j++)
		delete mx[j]; 
	delete mx;
	for(int i=0;i<num_gene;i++){
		if(has[i]==0)
			continue;
		if(rm[tt[i]]==1)
			has[i]=0;
	}
	return new_include_num;	           
}
int compact(float **input, float **sim, int *has, int *recordcol, int wh, int recordremove, int num_gene, int num_patient, int min_gene, float overlap, float cutoff ){ 
    map <int,int> tt; 
    int include_num=0;
    for(int k=0;k<num_gene;k++){
    	has[k]=0;
        float rr=newspearman(input, wh, k, recordcol,num_patient, num_patient-recordremove);
        if(rr < cutoff)
        	continue;
        tt[k]=include_num;
        include_num++;
        has[k]=1;
    }
	
	return include_num;	           
}

void cormatrix(float **input, float *cor, int wh, int num_gene,int num_patient, int *col, int left_col){
	for(int j=0; j < num_gene; j++)
		cor[j]=newspearman(input, wh, j, col, num_patient, left_col);
}

float average(vector <float> x){
	float sum=0;
	int n=x.size();
	for(int i=0;i< n; i++)
		sum+=x[i];
	return sum/n;
}





