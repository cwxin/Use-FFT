#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int Generate_N(int p,int q, int r);
int Initial(double *x, double *y, int M);
int Initial_dct(double *x, double *y, int N);
int Bit_Reverse_Integer(int N,int p,int q,int r,double *x_r, double *y_r);
int FFT(double *x_r, double *x_i, double *y_r, double *y_i,double *z_r,double *z_i, int N ,int p,int q,int r);
double angle(int p, int n);
int Print_Complex_Vector(double *y_r, double *y_i, double N);
int dec_prime(int N,int *prime);

int main()
{
	int p, N, q, r, n, M, *prime, N1;
	double *y_r, *y_i, *x_r, *x_i,*z_r,*z_i;

	printf("請注意，在2N+2的質因數只能有2、3、5，此程式才能執行。");
	printf("\nPlease input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	M = Generate_N(p, q, r);
	N = 2*M+2;
	printf("N=(2^%d)*(3^%d)*(5^%d) = %d\n",p,q,r,M);

	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));
	z_r = (double *) malloc(N*sizeof(double));
	z_i = (double *) malloc(N*sizeof(double));
	prime = (int *) malloc(3*sizeof(int));
	
	for(n=0;n<N;n++)
	{
		x_r[n]=0.0;
		x_i[n]=0.0;
	}
	
	Initial(x_r, x_i, M);
	
	for(n=0;n<N;n++)
	{
		y_r[n]=0.0;
		y_i[n]=0.0;
	}
	Initial_dct(x_r, y_r, N);
	
	for(n=0;n<N;n++)
	{
		x_r[n]=y_r[n];
		x_i[n]=0.0;
	}
		
	N1=dec_prime(N,prime);
	p = prime[0];
	q = prime[1];
	r = prime[2];

	if(N1 == 1)
	{
	Bit_Reverse_Integer(N,p,q,r,x_r,y_r);
	FFT(x_r, x_i, y_r, y_i,z_r,z_i, N,p,q,r);
	Print_Complex_Vector(y_r, y_i, N);
	}
	else
	{
		printf("請再輸入一次。"); 
	}
	return 0;
}

int Generate_N(int p, int q, int r)
{
	int i, N = 1;
	for(i=0;i<p;++i) N = N * 2;
	for(i=0;i<q;++i) N = N * 3;
	for(i=0;i<r;++i) N = N * 5;	
	return N;
}

int Initial(double *x, double *y, int M)
{
	int n;
	for(n=0;n<M;++n)
	{
		x[n] = n;
		y[n] = 0.0;
	}
}

int Initial_dct(double *x, double *y,int N)
{
	int n;
	for(n=0;n<N/2;++n)
	{
		y[n+1]=x[n];
	}
	for(n=(N-2)/2;n>1;n--)
	{
		y[N-n]=-x[n-1];
	}
}

int dec_prime(int N, int *prime)
{
	int p, q, r, n, w;
	p = 0;
	q = 0;
	r = 0;
	while(N%2 == 0)
	{
		N = N/2;
		p += 1;
	}

	while(N%3 == 0)
	{
		N = N/3;
		q += 1;
	}
	
	while(N%5 == 0)
	{
			N = N/5;
			r += 1;
	}
	prime[0] = p;
	prime[1] = q;
	prime[2] = r;
	return N;
}

int Bit_Reverse_Integer(int N,int p,int q, int r,double *x_r,double *y_r)
{

	int group=1, i,j,k,n;
	int size=N;

	while(r>0)
	{

	for (i = 0; i<group; i++)
	{
		for(j=0;j<5;j++)
		{
			for(k=0;k<size/5;k++)
			{
				y_r[i*size+(j*size)/5+k] = x_r[i*size+j+5*k];
		    }
		}
	}
	for(n=0;n<N;n++)
	{
		x_r[n] = y_r[n];
	}

		group=group*5;
		size = size/5;
		r--;
	}
	
	while(q>0)
	{

	for (i = 0; i<group; i++)
	{

		for(j=0;j<3;j++)
		{
			for(k=0;k<size/3;k++)
			{
				y_r[i*size+(j*size)/3+k] = x_r[i*size+j+3*k];
		    }
		}
	}
	for(n=0;n<N;n++)
	{
		x_r[n] = y_r[n];
	}

		group=group*3;
		size = size/3;
		q--;
	}
		
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

int FFT(double *x_r, double *x_i, double *y_r, double *y_i,double *z_r,double *z_i ,int N, int a,int q,int r)
{
	int h;
	int n,p,k,q1,i,q4,q2,q3;
	double w_r, w_i,w_r1,w_i1,w_r2,w_i2,w_r3,w_i3,w_r4,w_i4;
	double t_r, t_i, t_r1, t_r2, t_i1, t_i2,t_r3,t_i3,t_r4,t_i4;
	n = 2;
	while(a>0)
	{
		for(k=0;k<n/2;k++)
		{
			for(p=k;p<N;p+=n)
			{
				q1 = p + n/2;
				w_r = cos(angle(0*p,n));
				w_i = sin(angle(0*p,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*p,n));
				w_i1 = sin(angle(1*p,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				z_r[p] = t_r+t_r1;
				z_i[p] = t_i+t_i1;

				w_r = cos(angle(0*q1,n));
				w_i = sin(angle(0*q1,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q1,n));
				w_i1 = sin(angle(1*q1,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				z_r[q1] = t_r+t_r1;
				z_i[q1] = t_i+t_i1;
			}
		}
			for (i=0;i<N;i++)
			{
				y_r[i] = z_r[i];
		        y_i[i] = z_i[i];
			}
	n = n*2;
	a--;
    }
    
	n = n/2*3;
	
	while(q>0)
	{
		for(k=0;k<n/3;k++)
		{
			for(p=k;p<N;p+=n)
			{
				q1 = p + n/3;
				q2 = p + 2*n/3;
				w_r = cos(angle(0*p,n));
				w_i = sin(angle(0*p,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*p,n));
				w_i1 = sin(angle(1*p,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*p,n));
				w_i2 = sin(angle(2*p,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				z_r[p] = t_r+t_r1+t_r2;
				z_i[p] = t_i+t_i1+t_i2;

				w_r = cos(angle(0*q1,n));
				w_i = sin(angle(0*q1,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q1,n));
				w_i1 = sin(angle(1*q1,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q1,n));
				w_i2 = sin(angle(2*q1,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				z_r[q1] = t_r+t_r1+t_r2;
				z_i[q1] = t_i+t_i1+t_i2;

				w_r = cos(angle(0*q2,n));
				w_i = sin(angle(0*q2,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q2,n));
				w_i1 = sin(angle(1*q2,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q2,n));
				w_i2 = sin(angle(2*q2,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				z_r[q2] = t_r+t_r1+t_r2;
				z_i[q2] = t_i+t_i1+t_i2; 
			}
		}
			for (i=0;i<N;i++)
			{
				y_r[i] = z_r[i];
		        y_i[i] = z_i[i];
			}
	n = n*3;
	q--;
    }
	n = n/3*5;
	
	while(r>0)
	{
		for(k=0;k<n/5;k++)
		{
			for(p=k;p<N;p+=n)
			{
				q1 = p + n/5;
				q2 = p + 2*n/5;
				q3 = p + 3*n/5;
				q4 = p + 4*n/5;
				w_r = cos(angle(0*p,n));
				w_i = sin(angle(0*p,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*p,n));
				w_i1 = sin(angle(1*p,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*p,n));
				w_i2 = sin(angle(2*p,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				w_r3 = cos(angle(3*p,n));
				w_i3 = sin(angle(3*p,n));
				t_r3 = w_r3*y_r[q3]-w_i3*y_i[q3];
				t_i3 = w_r3*y_i[q3]+w_i3*y_r[q3];
				w_r4 = cos(angle(4*p,n));
				w_i4 = sin(angle(4*p,n));
				t_r4 = w_r4*y_r[q4]-w_i4*y_i[q4];
				t_i4 = w_r4*y_i[q4]+w_i4*y_r[q4];
				z_r[p] = t_r+t_r1+t_r2+t_r3+t_r4;
				z_i[p] = t_i+t_i1+t_i2+t_i3+t_i4;
				
				w_r = cos(angle(0*q1,n));
				w_i = sin(angle(0*q1,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q1,n));
				w_i1 = sin(angle(1*q1,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q1,n));
				w_i2 = sin(angle(2*q1,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				w_r3 = cos(angle(3*q1,n));
				w_i3 = sin(angle(3*q1,n));
				t_r3 = w_r3*y_r[q3]-w_i3*y_i[q3];
				t_i3 = w_r3*y_i[q3]+w_i3*y_r[q3];
				w_r4 = cos(angle(4*q1,n));
				w_i4 = sin(angle(4*q1,n));
				t_r4 = w_r4*y_r[q4]-w_i4*y_i[q4];
				t_i4 = w_r4*y_i[q4]+w_i4*y_r[q4];
				z_r[q1] = t_r+t_r1+t_r2+t_r3+t_r4;
				z_i[q1] = t_i+t_i1+t_i2+t_i3+t_i4;

				w_r = cos(angle(0*q2,n));
				w_i = sin(angle(0*q2,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q2,n));
				w_i1 = sin(angle(1*q2,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q2,n));
				w_i2 = sin(angle(2*q2,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				w_r3 = cos(angle(3*q2,n));
				w_i3 = sin(angle(3*q2,n));
				t_r3 = w_r3*y_r[q3]-w_i3*y_i[q3];
				t_i3 = w_r3*y_i[q3]+w_i3*y_r[q3];
				w_r4 = cos(angle(4*q2,n));
				w_i4 = sin(angle(4*q2,n));
				t_r4 = w_r4*y_r[q4]-w_i4*y_i[q4];
				t_i4 = w_r4*y_i[q4]+w_i4*y_r[q4];
				z_r[q2] = t_r+t_r1+t_r2+t_r3+t_r4;
				z_i[q2] = t_i+t_i1+t_i2+t_i3+t_i4;
				
				w_r = cos(angle(0*q3,n));
				w_i = sin(angle(0*q3,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q3,n));
				w_i1 = sin(angle(1*q3,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q3,n));
				w_i2 = sin(angle(2*q3,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				w_r3 = cos(angle(3*q3,n));
				w_i3 = sin(angle(3*q3,n));
				t_r3 = w_r3*y_r[q3]-w_i3*y_i[q3];
				t_i3 = w_r3*y_i[q3]+w_i3*y_r[q3];
				w_r4 = cos(angle(4*q3,n));
				w_i4 = sin(angle(4*q3,n));
				t_r4 = w_r4*y_r[q4]-w_i4*y_i[q4];
				t_i4 = w_r4*y_i[q4]+w_i4*y_r[q4];
				z_r[q3] = t_r+t_r1+t_r2+t_r3+t_r4;
				z_i[q3] = t_i+t_i1+t_i2+t_i3+t_i4;
				
				w_r = cos(angle(0*q4,n));
				w_i = sin(angle(0*q4,n));
				t_r = w_r*y_r[p]-w_i*y_i[p];
				t_i = w_r*y_i[p]+w_i*y_r[p];
				w_r1 = cos(angle(1*q4,n));
				w_i1 = sin(angle(1*q4,n));
				t_r1 = w_r1*y_r[q1]-w_i1*y_i[q1];
				t_i1 = w_r1*y_i[q1]+w_i1*y_r[q1];
				w_r2 = cos(angle(2*q4,n));
				w_i2 = sin(angle(2*q4,n));
				t_r2 = w_r2*y_r[q2]-w_i2*y_i[q2];
				t_i2 = w_r2*y_i[q2]+w_i2*y_r[q2];
				w_r3 = cos(angle(3*q4,n));
				w_i3 = sin(angle(3*q4,n));
				t_r3 = w_r3*y_r[q3]-w_i3*y_i[q3];
				t_i3 = w_r3*y_i[q3]+w_i3*y_r[q3];
				w_r4 = cos(angle(4*q4,n));
				w_i4 = sin(angle(4*q4,n));
				t_r4 = w_r4*y_r[q4]-w_i4*y_i[q4];
				t_i4 = w_r4*y_i[q4]+w_i4*y_r[q4];
				z_r[q4] = t_r+t_r1+t_r2+t_r3+t_r4;
				z_i[q4] = t_i+t_i1+t_i2+t_i3+t_i4;				
			}
		}
			for (i=0;i<N;i++)
			{
				y_r[i] = z_r[i];
		        y_i[i] = z_i[i];
			}
	n = n*5;
	r--;
	}
	return 0; 
}

double angle(int p, int n)
{
	double a;
	a = -2.0*p*M_PI/n;
	return a;
}

int Print_Complex_Vector(double *y_r, double *y_i, double N)
{
	int n;
	printf("%d : %f\n", 0, (-0.5)*y_i[1]);	
	for(n=1;n<(N-2)/2;n++)
	{
		printf("%d : %f\n", n, (-0.5)*(y_i[n+1]));		
	}
	return 0;
}
