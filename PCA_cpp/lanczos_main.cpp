#include "lanczos_functions.h"

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

int n=30;
char *file_matrix="example_1_covariance_matrix.dat";


int n_start = 1;
int Num_of_vectors_in_updated_B = n_start;//Total Size of orthogonal matrix to be passed to Gram-Schmidt



Matrix a = Matrix(n,n);

a = readusual(file_matrix);
std::cout << a << std::endl;


/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/

ComplexMatrix B(n,Num_of_vectors_in_updated_B);
ComplexMatrix eigen_vectors_after_swap(n,Num_of_vectors_in_updated_B);
ComplexMatrix eigen_values(Num_of_vectors_in_updated_B,1);

int num=n;//Number of vectors in B

ComplexMatrix alpha(num,1);
ComplexMatrix beta(num,1);

ComplexMatrix w(n,1);
ComplexMatrix v(n,1);

for(int i=0;i<n;i++){
	w(i,0) = uniform_random(-1,1);
}
Complex norm_init = NormComplex( w , n );

for(int i=0;i<n;i++){
	w(i,0) = w(i,0)/norm_init;
}

for(int i=0;i<n;i++){
	v(i,0) = 0.0;
}

/* Toeplitz Matrix-Vector Multiplication */
// v = v + toeplitz_symmetric_mat_vec_multiply(a,w,n);
/*************************************************/


/* Conventional Matrix-Vector Multiplication */
v = v + a*w;
/*************************************************/

for(int i=0;i<n;i++){
	B(i,0) = w(i,0);
}

alpha(0,0) = DotComplex( w , v , n );

v = v - alpha(0,0)*w;
beta(0,0) = NormComplex( v , n );

// std::cout << "beta " << beta << std::endl;

std:: ofstream file;
char fileName[250];

// while( beta( Num_of_vectors_in_updated_B -1 ,0).real() > 0.0001 ){
for( int p=1;p<5;p++ ){

Num_of_vectors_in_updated_B +=1;

// toeplitz_lanczos_main( a, n, Num_of_vectors_in_updated_B, B, eigen_vectors_after_swap, eigen_values, alpha, beta, w, v);
conventional_lanczos_main( a, n, Num_of_vectors_in_updated_B, B, eigen_vectors_after_swap, eigen_values, alpha, beta, w, v);
 
// sprintf(fileName,"./eigen_values_lanczos_main/%d.csv",Num_of_vectors_in_updated_B ) ;

sprintf( fileName, "./eigen_values_lanczos_main/test_july_13/%d.csv", Num_of_vectors_in_updated_B ) ;


file.open(fileName) ;

file << Num_of_vectors_in_updated_B   ;

for(int i=0;i<Num_of_vectors_in_updated_B;i++){
	file << " " << eigen_values(i,0).real()  ;
	}

file << std::endl;
file.close();

}


std::cout <<"eigenvectors " << std::endl;
std::cout << eigen_vectors_after_swap << std::endl;




std:: ofstream file_eigen;
char fileName_eigen[250];

sprintf( fileName_eigen, "./top_2_eigen_vectors.csv") ;
file_eigen.open(fileName_eigen) ;


for(int i=0;i<n;i++){
	file_eigen << eigen_vectors_after_swap(i,0).real() << "," << eigen_vectors_after_swap(i,1).real() << std::endl; 
}



file_eigen.close();




return 0;

}







