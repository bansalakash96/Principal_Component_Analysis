#include <octave/oct.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <chrono>

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count()) ;
const long int r_max= 4294967295;   

double uniform_random(double a,double b)
{
    double res;
  

    
    res=((double)rng()/(r_max))*(b-a)+a;
 
    return res;
  
}

/* REMEMBER : We are using routines of octave and INDEX starts from 0 HERE AS WELL */

// array       : It should be updated with the values in the absolute Descending Order when this array comes out of this function
// array_index : This should give you the indices of array where the elements were stored before sorting them.
void SORT_absolute_Descending( double array[] , int array_index[] ,int n ){

	double abs_array[n];
	double sorted_array[n];

	for(int i=0;i<n;i++){

			abs_array[i] = fabs(array[i]);

                   }
 
	for( int loop=0; loop<n; loop++ ){

			double max_value = ( *( std::max_element( abs_array , abs_array + n ) ) );

					for( int i=0;i<n;i++ ){

							if( abs_array[i] == max_value ){

								sorted_array[loop] = array[i];
 								array_index[loop] = i;
           
           						abs_array[i] = 0;
	     								           }	
     						             }

							     }

	for(int i=0;i<n;i++){

		array[i] = sorted_array[i];

					}		     

}

std::complex <double> NormComplex(ComplexMatrix x,int n){ 

	std::complex <double> ans=0.0;//Important to initialise with 0

	for(int i=0;i<n;i++){

		ans += x(i,0)*x(i,0); 
	}

	return std::sqrt(ans);
}

 std::complex <double> DotComplex(ComplexMatrix x,ComplexMatrix y,int n){

	std::complex <double> ans;//replace by complex

	for(int i=0;i<n;i++){

		ans += x(i,0)*y(i,0); 
	}

	return ans;
}


/* Column_Mum will start from 0 to n-1 */
ComplexMatrix column_vector(ComplexMatrix a,int n,int Column_Num){

	ComplexMatrix result(n,1);

	for(int i=0;i<n;i++){

		result(i,0) = a(i,Column_Num);
	}

	return result;

}


Matrix readusual(const char* filename){
int n;
std::cout << "Enter the order of square matrix "<< std::endl;
std::cin >> n;
Matrix a(n,n,0);
std ::ifstream in1;
in1.open(filename);

for (int i=0;i<n; i++){
for (int j=0;j<n;j++){
in1 >> a(i,j);
}
}
//cout <<"square matrix : "<< a << endl;
return a;
}


int main(){

int n = 30;
char *file_matrix="example_1_covariance_matrix.dat";


Matrix a = Matrix(n,n);

a = readusual(file_matrix);


std::cout << "Matrix A of size  (" << n << "x" << n << ") is : " << std::endl;

/*******************************************************************************************************************************/

 /* ACTUAL EIGENPAIRS */

EIG eig;
eig = EIG(a);

ComplexMatrix eigen_values  = eig.eigenvalues();
ComplexMatrix eigen_vectors = eig.right_eigenvectors();

// std::cout << "eigenvalues before sort are " << std::endl;
// std::cout << eigen_values << std::endl;

// std::cout << "eigenvectors for Matrix : a  : before sort are as follows : " << std::endl;
// std::cout << eigen_vectors << std::endl;
// std::cout << std::endl;
std::cout << std::endl;

/* Finding real part of eigenvalues in order to use a function for absolute values of eigenvalues in the descending order*/
/************************************************************************************/

double eigen_real_Step1[n];
int array_index_Step1[n];
for(int i=0;i<n;i++){
	eigen_real_Step1[i] = eigen_values(i,0).real();
}

SORT_absolute_Descending( eigen_real_Step1 , array_index_Step1 , n);

for(int i=0;i<n;i++){
	eigen_values(i,0) = eigen_real_Step1[i];
}

ComplexMatrix eigen_vectors_after_swap(n,n);

/* Sorting eigenvectors in descending order of absolute magnitude of eigenvalues */
/******************************************************************************************/
	for(int j=0;j<n;j++){
			for(int i=0;i<n;i++){
				eigen_vectors_after_swap(i,j) = eigen_vectors(i,array_index_Step1[j]);
								}
							  }

std::cout << std::endl;
std::cout << "Eigenvalues after sorting are: " << std::endl;
std::cout << eigen_values << std::endl;


std::cout << "Eigenvectors after swap based on descending order of eigenvalues for their absolute values are: " << std::endl;


// std :: cout <<  eigen_vectors_after_swap << std :: endl;


std :: cout <<  column_vector(eigen_vectors_after_swap,n,0) << std :: endl;

std::cout << std::endl;
std::cout << std::endl;
std::cout << std::endl;
std :: cout <<  column_vector(eigen_vectors_after_swap,n,1) << std :: endl;


return 0;

}
