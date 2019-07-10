#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int Generate_N(int p,int q, int r);
int Bit_Reverse_Integer(int N,int p,int q,int r,double *x_r, double *y_r);
int Initial(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int FFT(double *x_r, double *x_i, double *y_r, double *y_i,double *z_r,double *z_i, int N ,int p,int q,int r);
int Print_Complex_Vector(double *y_r, double *y_i, int N);
double angle(int p, int n);

int main()
{
	int p, N, q, r,n;
	double *y_r, *y_i, *x_r, *x_i,*z_r,*z_i;
	clock_t t1, t2;
	
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=(2^%d)*(3^%d)*(5^%d) = %d\n",p,q,r,N);
	
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));
	z_r = (double *) malloc(N*sizeof(double));
	z_i = (double *) malloc(N*sizeof(double));
	
	Initial(x_r, x_i, y_r, y_i, N);

	t1 = clock();
	Bit_Reverse_Integer(N,p,q,r,x_r,y_r);
	FFT(x_r, x_i, y_r, y_i,z_r,z_i, N,p,q,r);
	t2 = clock();
	
	//Print_Complex_Vector(y_r, y_i, N);   /*Print the result by FFT*/
	printf("Fast FT: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	
	free;
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

int Initial(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		x_r[n] = n;
		x_i[n] = 0.0;
		y_r[n] = 0.0;
		y_i[n] = 0.0;	
	}
}

int Bit_Reverse_Integer(int N,int p,int q, int r,double *x_r,double *y_r)
{

	int group=1, i,j,k,n;
	int size=N;

	while(r>0)
	{
		//printf("group=%d,size=%d\n",group,size);

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
		//printf("group=%d,size=%d\n",group,size);

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
		//printf("group=%d,size=%d\n",group,size);

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

int Print_Complex_Vector(double *y_r, double *y_i, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		if (y_i[n] >= 0) printf("%d : %f + %.20f i\n", n, y_r[n], y_i[n]);		
		else printf("%d : %f %.20f i\n", n, y_r[n], y_i[n]);
	}
	return 0;
}

double angle(int p, int n)
{
	double a;
	a = -2.0*p*M_PI/n;
	return a;
}
