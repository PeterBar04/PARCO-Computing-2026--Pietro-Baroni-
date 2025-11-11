#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "timer.h"

using namespace std;

#define NUM_RUNS 10 //Number of runs to test time of parallel code

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
struct CSR{
	vector<int> pointer; //row  array
	vector<int> index;   //column index array
	vector<double> data;
	
	void print(){
		cout << "Matrix CSR format\n" <<endl;
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
	std::ifstream file("1138_bus.mtx");
	
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
/*
	int nrows = matrix.rows;
    int nnz   = matrix.nnz;

	csr.index = matrix.column_index;
	csr.data = matrix.data;
	
	csr.pointer.assign(nrows+1, 0); //initialize at 0
		
	// Conta quanti elementi per riga
	for (int i=0; i < nnz; i++)
	    csr.pointer[matrix.row_index[i]+1]++;
	//now each csr.pointer[i+1] contains the number of elements of row i
	
	// Calcola somma cumulativa, inizio di ciascuna riga
	for (int i=0; i < nrows; i++)
	    csr.pointer[i+1] += csr.pointer[i];
*/

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
//Multiplication of matrix and vector
void matrix_vector_mul(CSR csr, vector<double> &v, vector<double> &r, double& time, bool serial){

	int rows = (int)r.size();
	
	double start=0;
	double end=0;
	
	if(serial){
		GET_TIME(start);
		for(int ip=0; ip<rows; ip++){		
			double sum = 0.0;
		    for (int k = csr.pointer[ip]; k < csr.pointer[ip+1]; k++) {
		        sum += csr.data[k] * v[ csr.index[k] ];
		    }
		    r[ip] = sum;
		}
		GET_TIME(end);
	} else{
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
	}
		
	time = end - start;
}


int main(int argc, char** argv) {
	
	srand(time(0));
	
	Matrix matrix;
	CSR csr;

	double time_serial = 0;
	double time_parallel=0; 
	double total_time_parallel=0;

	read_matrix_from_file(matrix);
	
	vector<double> vec(matrix.cols, 0); //initialized to 0
	vector<double> result(matrix.cols, 0); 
	
	convert_COO_in_CSR(matrix,csr);
	
	init_vec(vec,matrix.cols);
	
	//compute time for serial code
	matrix_vector_mul(csr, vec, result, time_serial, 1);
	printf("Serial Time = %f seconds\n\n", time_serial);
	
	//compute the time of running the matrix vector multiplication with openMP
	for (int num_threads = 1; num_threads <= 64; num_threads *= 2) {
        
		omp_set_num_threads(num_threads);
        total_time_parallel=0;
        
		for(int run=0; run< NUM_RUNS; run++){
			matrix_vector_mul(csr, vec, result, time_parallel, 0);
			total_time_parallel += time_parallel;
		}
		
		// Compute average parallel time, speedup, and efficiency
	    double avg_time_parallel = total_time_parallel / NUM_RUNS;
	    double avg_speedup = time_serial / avg_time_parallel;
	    double avg_efficiency = avg_speedup / num_threads;
		
		// Print results for each thread count
        printf("%11d | %17f | %7.2f | %9.2f%%\n", 
               num_threads, avg_time_parallel, avg_speedup, avg_efficiency * 100);

	}
	
	


	return 0;
}
