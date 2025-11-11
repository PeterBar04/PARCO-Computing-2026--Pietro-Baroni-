#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "timer.h"

using namespace std;

//Matrix in COO format read from .mtx
typedef struct Matrix{
	int rows;
	int cols;
	int nnz;
	
	//COO format
	vector<int> row_index; 		//row index array
	vector<int> column_index;   //column index array
	vector<double> data;	
	
	void print(){
		cout << "Matrix COO format\n" <<endl;
		cout << "row_index: ";
		for(int i: row_index){
			cout << i << " ";
		} cout << endl;
		
		cout << "column_index: ";
		for(int i: column_index){
			cout << i << " ";
		} cout << endl;
		
		cout << "data: ";
		for(int i: data){
			cout << i << " ";
		} cout << endl;
	}		
}Matrix;


//---------------------------------------//
//Matrix in CSR format
typedef struct CSR{
	vector<int> pointer; //row  array
	vector<int> index;   //column index array
	vector<double> data;
	
	void print(){
		cout << "POINTER: ";
		for(int i: pointer){
			cout << i << " ";
		} cout << endl;
		
		cout << "INDEX: ";
		for(int i: index){
			cout << i << " ";
		} cout << endl;
		
		cout << "DATA: ";
		for(int i: data){
			cout << i << " ";
		} cout << endl;
	}
};

//---------------------------------------//
//Read matrix from file and get COO format
void read_matrix_from_file(Matrix& matrix){
	std::ifstream file("matrix1.mtx");
	
	// Ignore comments headers
	while (file.peek() == '%') file.ignore(2048, '\n');
	
	// Read number of rows and columns
	file >> matrix.rows >> matrix.cols >> matrix.nnz;
	
	// fill the matrix with data
	for (int l = 0; l < matrix.nnz; l++)
	{
	    double data;
	    int row, col;
	    file >> row >> col >> data;    
	    
	    matrix.row_index.push_back(row-1); 		//mtx indexes start from 1
	    matrix.column_index.push_back(col-1);
	    matrix.data.push_back(data);
	}
	
	file.close();
}

//---------------------------------------//
//Convert Matrix from COO to CSR
void convert_COO_in_CSR(Matrix& matrix, CSR& csr){

	int nrows = matrix.rows;
    int nnz   = matrix.nnz;

    // --- Combina i vettori in un unico array di triple
    struct Entry {
        int row;
        int col;
        double val;
    };

    vector<Entry> entries(nnz);
    for (int i = 0; i < nnz; i++) {
        entries[i] = {matrix.row_index[i], matrix.column_index[i], matrix.data[i]};
    }

    // --- 2 Ordina per riga, poi per colonna
    sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
        if (a.row == b.row)
            return a.col < b.col;
        return a.row < b.row;
    });

    // ---  Copia nei vettori CSR (ordinati)
    csr.data.resize(nnz);
    csr.index.resize(nnz);
    csr.pointer.assign(nrows + 1, 0);

    for (int i = 0; i < nnz; i++) {
        csr.data[i] = entries[i].val;
        csr.index[i] = entries[i].col;
        csr.pointer[entries[i].row + 1]++;  // conta quanti elementi per riga
    }

    // ---  Calcola la somma cumulativa (prefix sum)
    for (int i = 0; i < nrows; i++) {
        csr.pointer[i + 1] += csr.pointer[i];
    }

}

//---------------------------------------//
//Initialize vector for multiplication with random values from 1 to 100
void init_vec(vector<double> &v, int cols){
	for(int i=0; i<cols; i++){
		v[i] = rand() % 100 + 1;
	}
}

//---------------------------------------//
//Multiply matrix and vector
void matrix_vector_mul(CSR csr, vector<double> &v, vector<double> &r,
						const string &schedule_type, int chunk_size){

	int rows = (int)r.size();
	
	omp_sched_t sched_kind;
    if (schedule_type == "static") sched_kind = omp_sched_static;
    else if (schedule_type == "dynamic") sched_kind = omp_sched_dynamic;
    else if (schedule_type == "guided") sched_kind = omp_sched_guided;
    else sched_kind = omp_sched_auto;

    omp_set_schedule(sched_kind, chunk_size);
	
	double start, end, elapsed;
	GET_TIME(start);
	
	#pragma omp parallel for
	for(int ip=0; ip<rows; ip++){		
		double sum = 0.0;
		    for (int k = csr.pointer[ip]; k < csr.pointer[ip+1]; k++) {
		        sum += csr.data[k] * v[ csr.index[k] ];
		    }
		    r[ip] = sum;
	}
	
	GET_TIME(end);
	
	elapsed = (end-start)/1000;	
	printf("Time %e\n", elapsed);
}


int main(int argc, char** argv) {
	
	srand(time(0));
	
	if (argc < 3) {
        cout << "Usage: ./deliverable1 <schedule_type> <chunk_size> [num_runs]\n" << endl;
        return -1;
    }
	
	string schedule_type = argv[1];
    int chunk_size = atoi(argv[2]);
	int num_runs = atoi(argv[3]);
	
	Matrix matrix;
	CSR csr;

	read_matrix_from_file(matrix);
	
	vector<double> vec(matrix.cols, 0); //initialized to 0
	vector<double> result(matrix.cols, 0); 
	
	convert_COO_in_CSR(matrix,csr);
	
	init_vec(vec,matrix.cols);

	matrix_vector_mul(csr, vec, result, schedule_type, chunk_size);
	
	return 0;
}

