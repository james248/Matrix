package com.jcdeck.matrix;

import java.util.Arrays;
import java.util.Random;

/**
 * This holds a two dimensional array of doubles
 * 
 * @author James Decker
 *
 */
public class Matrix {
	
    protected int m;             // number of rows
    protected int n;             // number of columns
    protected double[][] data;   // M-by-N array

    /**
     * Constructs an {@code m}-by-{@code n} matrix of 0's. If either dimension is
     * 0 or negative, nothing will happen;
     * 
     * @param m number of rows
     * @param n number of columns
     */
    public Matrix(int m, int n) {
    	if(m>=0 && n>=0){
	        this.m = m;
	        this.n = n;
	        data = new double[m][n];
    	}
    }

    /**
     * Constructs a matrix based on 2d array. If data is null or if has
     * a dimension of 0 in either rows or columns, nothing will happen.
     * 
     * @param data the data that will be put in the matrix
     */
    public Matrix(double[][] data) {
    	if(data != null && data.length>0){
	        if(data[0].length > 0){
		        m = data.length;
		        n = data[0].length;
		        this.data = new double[m][n];
		        for (int i = 0; i < m; i++)
		            for (int j = 0; j < n; j++)
		                    this.data[i][j] = data[i][j];
	        }
    	}
    }

    /**
     * Constructs a copy of {@code other}
     * 
     * @param other a matrix to copy
     */
    protected Matrix(Matrix other) {
    	this(other.data);
    }

    
    /**
     * Constructs and return a random {@code m}-by-{@code n} matrix with values between 0 and 1 with
     * even distribution. If either dimension is less than or equal to zero, null will be returned.
     * 
     * @param m The number of rows in the matrix
     * @param n The number of columns in the matrix
     * @return A new matrix with random values between 0 and 1
     */
    public static Matrix random(int m, int n) {
    	
    	if(m<=0 && n<=0)
        	return null;
    	
        Matrix A = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A.data[i][j] = Math.random();
        return A;
    }

    /**
     * Constructs and returns a random {@code M}-by-{@code N} matrix with Gaussian values
     * with median 0 and variance 1. If either dimension is less than or equal to zero,
     * null is returned
     * 
     * @param m The number of rows in the matrix
     * @param n The number of columns in the matrix
     * @return A new matrix with random Gaussian values
     */
    public static Matrix randomGaussian(int m, int n){
    	if(m<=0 && n<=0)
        	return null;
    	Random rd = new Random();
    	Matrix A = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A.data[i][j] = rd.nextGaussian();
        return A;
    }
    
    /**
     * Constructs and returns the {@code N}-by-{@code N} identity matrix. If {@code N} is
     * less than or equal to zero, null is returned.
     * 
     * @param N The width and height of the matrix
     * @return {@code N}-by-{@code N} identity matrix
     */
    public static Matrix identity(int N) {
    	if(N <= 0)
    		return null;
        Matrix I = new Matrix(N, N);
        for (int i = 0; i < N; i++)
            I.data[i][i] = 1;
        return I;
    }

    /**
     * swap rows {@code i} and {@code j}. If either value is out of bounds, nothing will happen.
     * 
     * @param i first row to swap
     * @param j second row to swap
     */
    private void swap(int i, int j) {
    	if(i<0 || i>=this.m || j<0 || j>=this.m)
    		return;
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    /**
     * Creates and returns the transpose of the invoking matrix
     * 
     * @return The transpose of the invoking matrix
     */
    public Matrix transpose() {
        Matrix A = new Matrix(n, m);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    /**
     * return C = A + B
     * 
     * @param B The matrix to add to the existing matrix
     * @return The sum of two matrices
     */
    public Matrix plus(Matrix B) {
        Matrix A = this;
        if (B.m != A.m || B.n != A.n)
        	return null;
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }

    /**
     * Return C = A - B. Returns null if the matrices are different dimensions
     * 
     * @param B the matrix to subtract from the existing matrix.
     * @return The resulting matrix after subtracting B from A
     */
    public Matrix minus(Matrix B) {
        Matrix A = this;
        if (B.m != A.m || B.n != A.n){
        	//throw new RuntimeException("Illegal matrix dimensions.");
        	return null;
        }
        
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        
        return C;
    }

    /**
     * Returns true if every element in other is equal to this matrix. If either matrix has illegal
     * dimensions or {@code other} is not a matrix, false will be returned.
     * 
     * @param other matrix to check equivalence
     * @return true if matrices are equal.
     */
    @Override
    public boolean equals(Object other) {
    	
    	if(!this.legalMatrix())
    		return false;
    	
    	Matrix B;
    	try{
    		B = (Matrix) other;
    		
    		if(!B.legalMatrix())
        		return false;
    		
    	}catch(Exception e){
    		System.out.println("The Object Passed was not of type Matrix");
    		return false;
    	}
        Matrix A = this;
        
        if (B.m != A.m || B.n != A.n)
        	return false;
        
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        
        return true;
    }

    /**
     * Returns {@code C} = {@code A} * {@code B}. Dot product for 1d matrices. If either
     * matrix has illegal dimension, null is returned.
     * 
     * @param B The matrix to multiply this by
     * @return The product of {@code A} and {@code B}
     */
    public Matrix times(Matrix B) {
    	
    	if(!this.legalMatrix() || !B.legalMatrix())
    		return null;
    	
        Matrix A = this;
        if (A.n != B.m) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.m, B.n);
        for (int i = 0; i < C.m; i++)
            for (int j = 0; j < C.n; j++)
                for (int k = 0; k < A.n; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }

    /**
     * Return a new matrix where every element in this one is multiplied by d
     * 
     * @param d A scaler to multiply with this matrix
     * @return The Product of A and d
     */
    public Matrix multiplyBy(double d){
    	
    	Matrix C = new Matrix(m, n);
    	for(int i = 0; i<data.length; i++)
    		for(int j = 0; j<data[0].length; j++)
    			C.set(data[i][j] * d, i, j);
    	return C;
    	
    }
    
    /**
     * Element wise multiplication. Returns a new matrix where each element
     * in C is that element in A time that element in {@code B}. returns null if
     * matrices are not the same direction.
     * 
     * @param B
     * @return new matrix of the same dimensions where each element has been multiplied by that element in {@code B}
     */
    public Matrix multiplyBy(Matrix B){
    	
    	if(!this.legalMatrix())
    		return null;
    	
    	Matrix C = new Matrix(m, n);
    	for(int i = 0; i<m; i++)
    		for(int j = 0; j<n; j++)
    			C.set(data[i][j]*B.data[i][j], i, j);
    	return C;
    	
    }
    
    /**
     * Element wise exponential function. If the matrix has illegal dimensions, null
     * is returned
     * 
     * @param d The exponent all elements of the matrix will be raised to
     * @return An new matrix of the same dimensions where each element has been raised to the {@code d}
     */
    public Matrix raiseTo(double d){
    	
    	if(!this.legalMatrix())
    		return null;
    	
    	Matrix c = new Matrix(m, n);
    	for(int i = 0; i<m; i++)
    		for(int j = 0; j<n; j++)
    			c.set(Math.pow(data[i][j], d), i, j);
    	return c;
    }
    
    /**
     * Calculates the sum of all the elements in the matrix. If the matrix
     * has illegal dimension, NaN is returned
     * 
     * @return the sum of all elements in the matrix
     */
    public double sum(){
    	
    	if(!this.legalMatrix())
    		return Double.NaN;
    	
    	double sum = 0;
    	for(int i = 0; i<m; i++)
    		for(int j = 0; j<n; j++)
    			sum += data[i][j];
    	return sum;
    }
    
    
    
    // return x = A^-1 b, assuming A is square and has full rank
    public Matrix solve2(Matrix rhs) {
        if (m != n || rhs.m != n || rhs.n != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this);
        Matrix b = new Matrix(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < n; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < n; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < n; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < n; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < n; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrix x = new Matrix(n, 1);
        for (int j = n - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < n; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;
   
    }

    /**
     * Prints the matrix to standard output. If the matrix has illegal dimension,
     * "Illegal matrix dimensions" will be printed
     */
    public void print() {
    	
    	if(!this.legalMatrix()){
    		System.out.println("Illegal matrix dimensions");
    		return;
    	}
    	
        for (int i = 0; i < m; i++)
            System.out.println(Arrays.toString(data[i]));
        
    }
    
    /**
     * Prints the matrix to the console rounding each element to {@code decimals}. If
     * the matrix has illegal dimension it will print "Illegal matrix dimensions".
     * 
     * @param decimals number of decimal digits to print
     * @see print()
     */
    public void print(int decimals){
    	
    	if(!this.legalMatrix()){
    		System.out.println("Illegal matrix dimensions");
    		return;
    	}
    	
    	for (int i = 0; i < m; i++) {
    		System.out.print("[");
    		for(int j = 0; j<n-1; j++){
    			System.out.print(((int) (data[i][j]*Math.pow(10, decimals))) / Math.pow(10, decimals));
    			System.out.print(", ");
    		}
			System.out.print(((int) (data[i][n-1]*Math.pow(10, decimals))) / Math.pow(10, decimals));
    		System.out.println("]");
        }
    }

    /**
     * Returns the value at row {@code i} and column {@code j}. If illegal
     * arguments are given or the matrix has illegal dimension, NaN will be
     * returned.
     * 
     * @param i row to get value from
     * @param j column to get value from
     * @return the value at row {@code i} and column {@code j}
     */
    public double get(int i, int j){

    	if(!this.legalMatrix())
    		return Double.NaN;
    	
    	if(i<0 || i>=this.getM() || j<0 || j>=this.getN())
    		return Double.NaN;
    	
    	return data[i][j];
    }
    public int getM(){
    	return m;
    }
    public int getN(){
    	return n;
    }
    
    /**
     * Sets the value at row {@code i} and column {@code j} to {@code val}.
     * If the matrix has illegal dimensions or will throw an out of bounds
     * exceptions, nothing will be changed.
     * 
     * @param val new value for one element of the matrix
     * @param i row index to change
     * @param j column index to change
     */
    public void set(double val, int i, int j){
    	
    	if(!this.legalMatrix())
    		return;
    	
    	if(i<0 || i>=this.getM() || j<0 || j>=this.getN())
    		return;
    	
    	data[i][j] = val;
    }
    
    /**
     * Returns a matrix with row i and column j removed. This 
     * new matrix will be one smaller in both dimensions.
     * 
     * @param i The row to remove
     * @param j The column to remove
     * @return A matrix with a row and column removed
     */
    public Matrix getWithout(int i, int j){
    	
    	if(i==-1)
    		return new Matrix(this);
    	
    	double[][] smallerMatrix = new double[this.m-1][this.n-1];
    	
    	//fill smaller matrix with the passed matrix without i and j
    	for(int k = 0; k<m-1; k++){
    		for(int h = 0; h<n-1; h++){
    			//h and k are the index of the smaller matrix that will be modified
    			smallerMatrix[k][h] = data[(k>=i)?k+1:k] [(h>=j)?h+1:h];
    		}
    	}
    	
    	
    	Matrix mat = new Matrix(smallerMatrix);
    	return mat;
    	
    }
    
    /**
     * Solves the matrix. If x and y are -1 nothing is removed
     * 
     * @param y removes this row before solving
     * @param x removes this column before solving
     * @return the value of the matrix
     */
    public static double solveDeterminant(Matrix matrix){
    	
    	//make sure it is a square
    	if(matrix.m != matrix.n)
    		throw new RuntimeException("Illegal Matrix Demention");
    	
    	
    	//if it has gotten down to a 2x2 matrix return the solution
    	if (matrix.m == 2){
    		//matrix.print();
    		return matrix.get(0, 0)*matrix.get(1, 1) - matrix.get(1, 0)*matrix.get(0, 1);
    	}
    	
    	final int size = matrix.m;
    	
    	double sum = 0;
    	for(int i = 0; i<size; i++){
    		sum += matrix.get(0, i) * solveDeterminant(matrix.getWithout(0, i)) * ((i%2==1)?-1:1);
    		//System.out.println(sum);
    	}
    	
    	
    	return sum;
    	
    }
    
    /**
     * Solves the minor of {@code a}. If {@code a} has illegal dimensions,
     * null is returned
     * 
     * @param a matrix to solve the minor of
     * @return minor of {@code a}
     */
    public static Matrix solveMinor(Matrix a){

    	if(!a.legalMatrix())
    		return null;
    	
    	Matrix co = new Matrix(a);
    	
    	final int size = co.m;
    	
    	for(int i = 0; i<size; i++){
    		for(int j = 0; j<size; j++){
    			co.set(solveDeterminant(a.getWithout(i, j)), i, j);
    		}
    	}
    	
    	return co;
    }
    
    /**
     * Multiplies every other element in {@code matrix} by -1. If {@code matrix}
     * has illegal dimensions, null is returned
     * 
     * @param a matrix to operate on
     * @return a after operation
     */
    public static Matrix checkerboardNegation(Matrix a){
    	
    	if(!a.legalMatrix())
    		return null;
    	
    	Matrix x = new Matrix(a);
    	
    	for(int i = 0; i<a.m; i++)
    		for(int j = 0; j<a.n; j++){
    			if((i+j)%2 == 1){
    				x.set(x.get(i, j)*-1, i, j);
    			}
    		}
    	
    	
    	return x;
    }
    
    /**
     * Calculates the inverse Matrix of {@code a}. If a has illegal dimensions
     * null is returned.
     * 
     * @param a matrix to invert
     * @return the inverse of {@code a}
     */
    public static Matrix invert(Matrix a) {
    	
    	if(!a.legalMatrix())
    		return null;
    	
    	
    	System.out.println("matrix.Matrix.invert(Matrix)");
        
    	if(a.m != a.n)
    		throw new RuntimeException("Invalid dimetions");
    	
    	System.out.println("    Solving Determinint...");
    	double determinant = solveDeterminant(a);
    	
    	System.out.println("    Solving Minor...");
    	Matrix minor = solveMinor(a);
    	
    	System.out.println("    Solving cofactor...");
    	Matrix cofactor = checkerboardNegation(minor);
    	
    	
    	System.out.println("    Solving Adjoint...");
    	Matrix adjoint = cofactor.transpose();
    	
    	
    	System.out.println("    Multiplying by determinant...");
    	adjoint.multiplyBy(1/determinant);
    	
    	return adjoint;
    	
    }
    
    
    /**
     * Returns true if the matrix has valid dimensions. Returns false
     * if data is null, or if either dimensions is less than or equal to 0.
     * Should be called on each matrix being used before every operation.
     * 
     * @return true is the matrix can be used in operations
     */
    private boolean legalMatrix(){
		if(data==null || m<=0 || n<=0){
			throw new RuntimeException("Illegal Matrix");
			//return false;
    	}
    	return true;
    }
    
}