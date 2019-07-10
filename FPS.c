#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

int dec_prime(int N);
int Segmentation(double **array, int N);
int Fast_Poisson_Solver(double **array, int N, int p, int M);
int intial(double *x_r, double *x_i, double *y_r, double *y_i, int M);
int Initial_dst(double *x, double *y,int N);
int Bit_Reverse_Integer(int N, int p, double *x_r, double *y_r);
int FFT(double *y_r, double *y_i ,int N, int a);
int Transpose(double **array, int N);
int Diagonal(double **array, int N);
int iDST2D(double **array, int N, int p, int M);
int iDST(double *x_r, double *x_i, double *y_r, double *y_i, int N, int p, int M);
int DST2D(double **array, int N, int p, int M);
int DST(double *x_r, double *x_i, double *y_r, double *y_i, int N, int p, int M);
double angle(int p, int n);

int main()
{
	int i, j, N, M, p, M_1;
	double **array;
	clock_t t1, t2;
	printf("請輸入N(必須是2的次方數)：");
	scanf("%d", &N);
	M_1=N-1;
	M=2*N;
	printf("未知的點數M為：%d\n", M_1*M_1);
	printf("測試函數：f(x,y)=(-5.0)*(M_PI*M_PI)*(sin(M_PI*x))*(sin(2.0*M_PI*y))");
	p=dec_prime(M);
	array=(double **)malloc((M_1)*sizeof(double*));
	for (i=0;i<M_1;i++)
	{
		array[i]=(double *)malloc((M)*sizeof(double));
	}
	Segmentation(array,N);
	t1 = clock();
	Fast_Poisson_Solver(array,N,p,M);
	t2 = clock();
	printf("\nFast Poisson Solver: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	
	/* //把算出來的數值解列印出來          
	for(j=0;j<N-1;j++)
	{
		for(i=0;i<N-1;i++)
		{
			printf("\narray[%d][%d] %f", j,i,array[j][i]);
		}		
	}*/
	 
	free(array);
}

int Segmentation(double **array, int N)
{
	int i, j;
	double k, e;
	k=(N-1.0)/N;
	for(j=0;j<N-1;j++)
	{
		for(i=0;i<N-1;i++)
		{
				array[j][i]=(-5.0)*(M_PI*M_PI)*(sin(M_PI*k))*(sin(2.0*M_PI*(i+1.0)/N));
		}
		k=k-(1.0/N);
	}
}

int Fast_Poisson_Solver(double **array, int N, int p, int M)
{
	iDST2D(array,N,p,M);
	Diagonal(array,N);
	DST2D(array,N,p,M);
}


int iDST(double *x_r, double *x_i, double *y_r, double *y_i, int N, int p, int M)
{
	int i;
	Initial_dst(x_r, y_r, M);
	for(i=0;i<M;i++)
	{
		x_r[i]=y_r[i];
	}
	Bit_Reverse_Integer(M,p,x_r,y_r);
	FFT(y_r,y_i, M,p);
}

int DST(double *x_r, double *x_i, double *y_r, double *y_i, int N, int p, int M)
{
	int i;
	Initial_dst(x_r, y_r, M);
	for(i=0;i<M;i++)
	{
		x_r[i]=y_r[i];
	}
	Bit_Reverse_Integer(M,p,x_r,y_r);
	FFT(y_r,y_i, M,p);
}

int DST2D(double **array, int N, int p, int M)
{
	int i, j, k;
	double *x_r, *x_i, *y_r, *y_i;
	x_r = (double *) malloc(M*sizeof(double));
	x_i = (double *) malloc(M*sizeof(double));
	y_r = (double *) malloc(M*sizeof(double));
	y_i = (double *) malloc(M*sizeof(double));
	for(i=0;i<N-1;i++)
	{
		intial(x_r,x_i,y_r,y_i,M);
		for(j=0;j<N-1;j++)
		{
			x_r[j]=array[i][j];
		}
		DST(x_r,x_i,y_r,y_i,N,p,M);
		array[i][0]=(-0.5)*y_i[1];
		for(k=2;k<=(M-2)/2;k++)
		{
			array[i][k-1]=(-0.5)*y_i[k];
		}
	}
	Transpose(array,N);
	for(i=0;i<N-1;i++)
	{
		intial(x_r,x_i,y_r,y_i,M);
		for(j=0;j<N-1;j++)
		{
			x_r[j]=array[i][j];
		}
		DST(x_r,x_i,y_r,y_i,N,p,M);
		array[i][0]=(-0.5)*y_i[1];
		for(k=2;k<=(M-2)/2;k++)
		{
			array[i][k-1]=(-0.5)*y_i[k];
		}
	}
	Transpose(array,N);
	for(j=0;j<N-1;j++)
	{
		for(i=0;i<N-1;i++)
		{
			array[j][i]=(4.0/(N*N))*array[j][i];
		}		
	}
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
}

int iDST2D(double **array, int N, int p, int M)
{
	int i, j, k;
	double *x_r, *x_i, *y_r, *y_i, h;
	x_r = (double *) malloc(M*sizeof(double));
	x_i = (double *) malloc(M*sizeof(double));
	y_r = (double *) malloc(M*sizeof(double));
	y_i = (double *) malloc(M*sizeof(double));
	for(i=0;i<N-1;i++)
	{
		intial(x_r,x_i,y_r,y_i,M);
		for(j=0;j<N-1;j++)
		{
			x_r[j]=array[i][j];
		}
		iDST(x_r,x_i,y_r,y_i,N,p,M);
		array[i][0]=(-0.5)*y_i[1];
		for(k=2;k<=(M-2)/2;k++)
		{
			array[i][k-1]=(-0.5)*y_i[k];
		}
	}
	Transpose(array,N);
	for(i=0;i<N-1;i++)
	{
		intial(x_r,x_i,y_r,y_i,M);
		for(j=0;j<N-1;j++)
		{
			x_r[j]=array[i][j];
		}
		iDST(x_r,x_i,y_r,y_i,N,p,M);
		array[i][0]=(-0.5)*y_i[1];
		for(k=2;k<=(M-2)/2;k++)
		{
			array[i][k-1]=(-0.5)*y_i[k];
		}
	}
	Transpose(array,N);
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
}

int Transpose(double **array, int N)
{
	int i, j;
	double v;
	for(i=0;i<N-1;++i)
	{
		for(j=i+1;j<N-1;++j)
		{
			v = array[i][j];
			array[i][j] = array[j][i];
			array[j][i] = v;
		}
	}
	return 0;
}

int Diagonal(double **array, int N)
{
	double lambda_x,lambda_y;
	int i, j;
	for(i=0;i<N-1;i++)
	{
		for(j=0;j<N-1;j++)
		{
			array[i][j]=-(array[i][j])/(((4.0)*N*N*(sin((i+1)*M_PI/(2.0*N)))*sin((i+1)*M_PI/(2.0*N)))+((4.0)*N*N*(sin((j+1)*M_PI/(2.0*N)))*sin((j+1)*M_PI/(2.0*N))));
		}
	}
}

int intial(double *x_r, double *x_i, double *y_r, double *y_i, int M)
{
	int i;
	for(i=0;i<M;i++)
	{
		x_r[i]=0.0;
		y_r[i]=0.0;
		x_i[i]=0.0;
		y_i[i]=0.0;
	}
}

int Initial_dst(double *x, double *y,int N)
{
	int n;
	for(n=0;n<(N-2)/2;++n)
	{
		y[n+1]=x[n];
	}
	for(n=(N-4)/2;n>=0;--n)
	{
		y[N-n-1]=-x[n];
	}
}

int Bit_Reverse_Integer(int N, int p, double *x_r, double *y_r)
{
	int group=1, i,j,k,n;
	int size=N;	
	while(p>0)
	{

	for (i = 0; i<group; i++)
	{
		for(j=0;j<2;j++)
		{
			for(k=0;k<size/2;k++)
			{
				y_r[i*size+(j*size)/2+k] = x_r[i*size+j+2*k];
		    }
		}
	}
	for(n=0;n<N;n++)
	{
		x_r[n] = y_r[n];
	}

		group=group*2;
		size = size/2;
		p--;
	}
	return 0;
}

int FFT(double *y_r, double *y_i,int N, int a)
{
	int n,p,k,q1,i;
	double w_r, w_i,theta;
	double t_r, t_i;
	n = 2;
	while(a>0)
	{
		for(k=0;k<n/2;k++)
		{
			theta = -2.0*k*M_PI/n;
			w_r = cos(theta);
			w_i = sin(theta);
			for(p=k;p<N;p+=n)
			{
				q1 = p + n/2;
				t_r = w_r*y_r[q1]-w_i*y_i[q1];
				t_i = w_r*y_i[q1]+w_i*y_r[q1];
				y_r[q1] = y_r[p] - t_r;
    			y_i[q1] = y_i[p] - t_i;
    			y_r[p] += t_r;
    			y_i[p] += t_i;
		
			}
		}
	n = n*2;
	a--;
    }
	return 0; 
}

int dec_prime(int N)
{
	int p;
	p = 0;
	while(N%2 == 0)
	{
		N = N/2;
		p += 1;
	}
	return p;
}
