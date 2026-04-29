#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


#define N 100
#define EPS 1e-6
#define a_value 10.0

void fill(double A[N][N], double b[N])
{
	for (int i = 0; i < N; i++)
	{
		b[i] = i + 1;

		for (int j = 0; j < N; j++)
		{
			if (i == j) A[i][j] = a_value;
			else if (j < i || j == i + 1) A[i][j] = 1.0 / (i + 1);
			else A[i][j] = 0.0;
		}
	}
}

void print_matrix(double A[N][N])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < ((i + 2 < N) ? (i + 2) : N); j++)
		{
			printf("%lf ", A[i][j]);
		}
		printf("\n");
	}
}
//----------------------------------------------------------------------------------------------------------------------
double vector_norm(double x[N])
{
	double norm;

	double sum = 0.0;

	for (int i = 0; i < N; i++)
	{
		sum += x[i] * x[i];
	}

	norm = sqrt(sum);

	return norm;
}

void normalize(double x[N])
{
	double n = vector_norm(x);
	for (int i = 0; i < N; i++)
	{
		x[i] /= n;
	}
}

double scalar_product(double a[N], double b[N])
{
	double scalar_product = 0.0;

	for (int i = 0; i < N; i++)
	{
		scalar_product += a[i] * b[i];
	}

	return scalar_product;
}

bool cont(double x[N], double x_prev[N])
{
	double resid[N];

	for (int i = 0; i < N; i++)
	{
		resid[i] = x[i] - x_prev[i];
	}

	return vector_norm(resid) / vector_norm(x_prev) > EPS;
}

void nullity(double A[N][N], double b[N], double x[N], FILE *file)
{
	double result[N];
	double nullity[N];
	double nullity_norm;
	
	int flag;

	for (int i = 0; i < N; i++)
	{
		result[i] = 0.0;
		nullity[i] = 0.0;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			result[i] += A[i][j] * x[j];
		}
	}

	printf("\tIf you want to see the whole nullity vector, please, enter 1, else - enter 0\n");

	scanf_s("%d", &flag);

	for (int i = 0; i < N; i++)
	{
		nullity[i] = fabs(b[i] - result[i]);
		if (flag == 1) printf("\tnullity[%d] = %.15f\n", i, nullity[i]);
		if (flag == 1) fprintf(file, "\tnullity[%d] = %.15f\n", i, nullity[i]);
	}

	nullity_norm = vector_norm(nullity);

	printf("\t||nullity|| = %.15f\n\n", nullity_norm);
	fprintf(file, "\t||nullity|| = %.15f\n\n", nullity_norm);
}
//----------------------------------------------------------------------------------------------------------------------
void invert_matrix(double A[N][N], double inverse[N][N]) 
{
	double augmented[N][2 * N];

	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			augmented[i][j] = A[i][j];
			augmented[i][j + N] = (i == j) ? 1.0 : 0.0;
		}
	}

	for (int i = 0; i < N; i++) 
	{
		int max_row = i;
		for (int k = i + 1; k < N; k++) 
		{
			if (fabs(augmented[k][i]) > fabs(augmented[max_row][i])) 
			{
				max_row = k;
			}
		}

		if (fabs(augmented[max_row][i]) < EPS) 
		{
			printf("Inverted natrix cannot be found\n");
			break;
		}

		if (i != max_row) 
		{
			for (int j = 0; j < 2 * N; j++) 
			{
				double temp = augmented[i][j];
				augmented[i][j] = augmented[max_row][j];
				augmented[max_row][j] = temp;
			}
		}

		double pivot = augmented[i][i];
		for (int j = 0; j < 2 * N; j++) 
		{
			augmented[i][j] /= pivot;
		}

		for (int k = 0; k < N; k++) 
		{
			if (k != i) 
			{
				double factor = augmented[k][i];
				for (int j = 0; j < 2 * N; j++) 
				{
					augmented[k][j] -= factor * augmented[i][j];
				}
			}
		}
	}

	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			inverse[i][j] = augmented[i][j + N];
		}
	}
}

void multiply_matrix_vector(double A[N][N], double x[N], double y[N])
{
	for (int i = 0; i < N; i++)
	{
		y[i] = 0.0;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			y[i] += A[i][j] * x[j];
		}
	}
}

void matrix_transpose(double A[N][N], double A_T[N][N])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A_T[j][i] = A[i][j];
		}
	}

}

void matrix_multiply(double A[N][N], double B[N][N], double result[N][N])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			result[i][j] = scalar_product(A[i], B[j]);
		}
	}
}
//----------------------------------------------------------------------------------------------------------------------
void seidel(double A[N][N], double b[N], double x[N], FILE *file)
{
	double x_prev[N];

	int flag;

	int iterations = 0;

	for (int i = 0; i < N; i++)
	{
		x[i] = 1.0;
	}

	do
	{
		for (int i = 0; i < N; i++)
		{
			x_prev[i] = x[i];
		}

		for (int i = 0; i < N; i++)
		{
			double sum1 = 0.0, sum2 = 0.0;

			for (int j = 0; j < i; j++)
			{
				sum1 += A[i][j] * x[j];
			}

			for (int j = i + 1; j < N; j++)
			{
				sum2 += A[i][j] * x_prev[j];
			}

			x[i] = (b[i] - sum1 - sum2) / A[i][i];
		}

		iterations++;

	} while (cont(x, x_prev));

	printf("THE SEIDEL METHOD\n");

	printf("\tIf you want to see the calculated root, please, enter 1\n");
	printf("\tIf you want to proceed, please, enter 0\n");

	fprintf(file, "THE SEIDEL METHOD\n");

	scanf_s("%d", &flag);

	if (flag == 1)
	{
		printf("\tSeidel method gave the root:\n");
		fprintf(file, "\tSeidel method gave the root:\n");

		for (int i = 0; i < N; i++)
		{
			printf("\tx[%d] = %.15f\n", i, x[i]);
			fprintf(file, "\tx[%d] = %.15f\n", i, x[i]);
		}
	}

	printf("\tThe vector of nullity is (b - Ax*), where x* is the found solution\n");
	fprintf(file, "\tThe vector of nullity is (b - Ax*), where x* is the found solution\n");

	nullity(A, b, x, file);
}

void gauss(double A[N][N], double b[N], double x[N], FILE *file)
{
	int flag;

	for (int i = 0; i < N; i++)
	{
		if (fabs(A[i][i]) < EPS)
		{
			printf("\tZero pivot element detected. Method fails.\n");
			return;
		}

		for (int k = i + 1; k < N; k++)
		{
			double factor = A[k][i] / A[i][i];

			for (int j = i; j < N; j++)
			{
				A[k][j] -= factor * A[i][j];
			}

			b[k] -= factor * b[i];
		}
	}

	for (int i = N - 1; i >= 0; i--)
	{
		x[i] = b[i];

		for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

		x[i] /= A[i][i];
	}

	printf("THE GAUSS METHOD\n");
	fprintf(file, "THE GAUSS METHOD\n");

	printf("\tIf you want to see the calculated root, please, enter 1\n");
	printf("\tIf you want to proceed, please, enter 0\n");

	scanf_s("%d", &flag);

	if (flag == 1)
	{
		printf("\tGauss method gave the root:\n");
		fprintf(file, "\tGauss method gave the root:\n");

		for (int i = 0; i < N; i++)
		{
			printf("\tx[%d] = %.15f\n", i, x[i]);
			fprintf(file, "\tx[%d] = %.15f\n", i, x[i]);
		}
	}

	printf("\tThe vector of nullity is (b - Ax*), where x* is the found solution\n");
	fprintf(file, "\tThe vector of nullity is (b - Ax*), where x* is the found solution\n");

	nullity(A, b, x, file);
}
//----------------------------------------------------------------------------------------------------------------------
double power_method(double A[N][N])
{
	double x[N], y[N], lambda = 0.0, lambda_prev = 0.0;
	int iter = 0;

	for (int i = 0; i < N; i++) x[i] = 1.0;
	normalize(x);

	do {
		lambda_prev = lambda;

		multiply_matrix_vector(A, x, y);
		lambda = scalar_product(y, x) / scalar_product(x, x);

		for (int i = 0; i < N; i++) x[i] = y[i];
		normalize(x);

		iter++;
	} while (fabs(lambda - lambda_prev) > EPS);

	return lambda;
}

double inverse_power_method(double A_1[N][N]) 
{
	double x[N], y[N], lambda = 0.0, lambda_prev = 0.0;
	int iter = 0;

	for (int i = 0; i < N; i++) x[i] = 1.0;
	normalize(x);

	do {
		lambda_prev = lambda;

		multiply_matrix_vector(A_1, x, y); // y = A^{-1} * x
		lambda = scalar_product(y, x) / scalar_product(x, x);

		for (int i = 0; i < N; i++) x[i] = y[i];
		normalize(x);

		iter++;
	} while (fabs(1.0 / lambda - 1.0 / lambda_prev) > EPS);

	return (1.0 /lambda);
}
//----------------------------------------------------------------------------------------------------------------------
double matrix_norm_1(double A[N][N])
{
	double max_sum = 0.0;
	for (int j = 0; j < N; j++) 
	{
		double col_sum = 0.0;
		for (int i = 0; i < N; i++) 
		{
			col_sum += fabs(A[i][j]);
		}
		if (col_sum > max_sum) 
		{
			max_sum = col_sum;
		}
	}
	return max_sum;
}

double matrix_norm_infinity(double A[N][N])
{
	double max_sum = 0.0;
	for (int i = 0; i < N; i++) 
	{
		double row_sum = 0.0;
		for (int j = 0; j < N; j++) 
		{
			row_sum += fabs(A[i][j]);
		}
		if (row_sum > max_sum) 
		{
			max_sum = row_sum;
		}
	}
	return max_sum;
}

double matrix_norm_2(double A[N][N])
{
	double A_T[N][N];
	double multiplied[N][N];
	double lambda_max = 0.0;

	double result;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A_T[i][j] = 0.0;
		}
	}

	matrix_transpose(A, A_T);

	matrix_multiply(A_T, A, multiplied);

	lambda_max = power_method(multiplied);

	result = sqrt(lambda_max);

	return result;
}

void mu(double A[N][N], double A_1[N][N], FILE *file)
{
	double mu_1, mu_infty, mu_2;

	mu_1 = matrix_norm_1(A) * matrix_norm_1(A_1);
	mu_infty = matrix_norm_infinity(A) * matrix_norm_infinity(A_1);
	mu_2 = matrix_norm_2(A) * matrix_norm_2(A_1);

	printf("\tCondition number for three matrix norms:\n");
	printf("\t\tmu_1 = %f, mu_infinity = %f, mu_2 = %f\n\n", mu_1, mu_infty, mu_2);

	fprintf(file, "\tCondition number for three matrix norms:\n");
	fprintf(file, "\t\tmu_1 = %f, mu_infinity = %f, mu_2 = %f\n\n", mu_1, mu_infty, mu_2);
}
//----------------------------------------------------------------------------------------------------------------------
int main(void)
{
	FILE* file;

	double A[N][N], b[N], x[N], A_copy[N][N], b_copy[N], A_1[N][N], lambda_max = 0.0, lambda_min = 0.0, norm_1, norm_infty, norm_2;
	
	file = fopen("output.txt", "w");

	if (!file)
	{
		printf("Error opening file for writing. \n");
		return 0;
	}

	fill(A, b);
	fill(A_copy, b_copy);

	invert_matrix(A, A_1);

	lambda_max = power_method(A);
	lambda_min = inverse_power_method(A_1);

	printf("INFORMATION ABOUT INITIAL MATRIX A\n\n");
	fprintf(file,"INFORMATION ABOUT INITIAL MATRIX A\n\n");

	printf("\tlambda_max = %f, lambda_min = %f\n\n", lambda_max, lambda_min);
	fprintf(file,"\tlambda_max = %f, lambda_min = %f\n\n", lambda_max, lambda_min);

	norm_1 = matrix_norm_1(A);
	norm_infty = matrix_norm_infinity(A);
	norm_2 = matrix_norm_2(A);

	printf("\t||A||_1 = %f\n", norm_1);
	printf("\t||A||_infinity = %f\n", norm_infty);
	printf("\t||A||_2 = %f\n\n", norm_2);

	fprintf(file, "\t||A||_1 = %f\n", norm_1);
	fprintf(file, "\t||A||_infinity = %f\n", norm_infty);
	fprintf(file, "\t||A||_2 = %f\n\n", norm_2);

	mu(A, A_1, file);

	seidel(A, b, x, file);

	gauss(A_copy, b_copy, x, file);

	return 0;
}


