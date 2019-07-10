#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int Odd_terms(int N);
int Initial(int *x_r, int *y_r, int *w_r, int *z_r, int N_i);
int Iuput_Initial(int *x_r,int *y_r,int N);
int Prime_Omega(int N, int N_i, int *w_r);
int Check_Prime(int N_p);
int Find_Omega(int N_p, int N_i, int *w_r);
int Check_Omega(int *w_r, int N_i);
int Bit_Reverse_Integer(int N_i, int *x_r, int *z_r);
int FFT(int *x_r, int *w_r, int N_i, int a, int prime);
int IFFT(int *x_r, int *w_r, int N_i, int a, int prime);
int result(int *x_r, int N_i);

int main()
{
	int N, N_i, *x_r, *y_r, *z_r, *w_r, n, P, power;
	printf("想要做N位數*N位數的乘法(請輸入N)：");
	scanf("%d", &N);
	if (N<=0)
	{
		printf("請輸入大於零的正整數");
	}
	N_i=Odd_terms(N);
	x_r = (int *) malloc(N_i*sizeof(int));
	y_r = (int *) malloc(N_i*sizeof(int));
	z_r = (int *) malloc(N_i*sizeof(int));
	w_r = (int *) malloc(N_i*sizeof(int));
	Initial(x_r, y_r, w_r, z_r, N_i);
	Iuput_Initial(x_r, y_r, N);
	w_r[0]=1;
	P=Prime_Omega(N, N_i, w_r);
	power=Bit_Reverse_Integer(N_i, x_r, z_r);
	FFT(x_r, w_r, N_i, power, P);
	power=Bit_Reverse_Integer(N_i, y_r, z_r);
	FFT(y_r, w_r, N_i, power, P);
	for(n=1;n<N_i;n++)
	{
		z_r[n]=(x_r[N_i-n]*y_r[N_i-n])%P;
	} 	
	z_r[0]=(x_r[0]*y_r[0])%P;
	power=Bit_Reverse_Integer(N_i, z_r, x_r);
	IFFT(z_r, w_r, N_i, power, P);
	for(n=1;n<N_i;n++)
	{
		x_r[n]=z_r[N_i-n];
	}
		x_r[0]=z_r[0];

	for(n=0;n<N_i;n++)
	{
		while(x_r[n]%N_i!=0)
		{
			x_r[n]=x_r[n]+P;
		}
		x_r[n]=x_r[n]/N_i;
	}
	printf("P=%d,", P);
	printf("w=%d\n", w_r[1]);
	result(x_r, N_i);
}

int Odd_terms(int N)
{
	int N_R, k, K_R, i;
	N_R=1;
	k=0;
	K_R=1;
	
	if(N==1)
	{
		K_R=4;
	}
	else
	{	
	while(N>N_R)
	{
		N_R=N_R*2;
		k=k+1;
	}
	k=k+1;
	}
	for(i=1;i<k+1;i++)
	{
		K_R=K_R*2;
	}
	return K_R;
}

int Initial(int *x_r, int *y_r, int *w_r, int *z_r, int N_i)
{
	int i;
	for(i=0;i<N_i;i++)
	{
		x_r[i]=0;
		y_r[i]=0;
		w_r[i]=0;
		z_r[i]=0;
	}
}

int Iuput_Initial(int *x_r, int *y_r, int N)
{
	int i, j;
	j=0;
	printf("\n請把想要乘的數字並且由個位數開始輸入:\n");
	printf("\n請務必只輸入一個數字(0~9)，否則會有問題。\n\n"); 
	for(i=0;i<N;i++)
	{
		printf("第一個數字的第%d位數:", i+1);
		scanf("%d", &x_r[i]);
	}
	printf("\n");
	for(i=0;i<N;i++)
	{
		printf("第二個數字的第%d位數:", i+1);
		scanf("%d", &y_r[i]);
	}
	printf("\n");
	printf("我們想要知道");
	for(i=N-1;i>=0;i--)
	{
		if(x_r[i]!=0)
		{
			j=i;
			break;
		}		
	}
	for(i=j;i>=0;i--)
	{
		printf("%d", x_r[i]);
	}
	printf("*");
	for(i=N-1;i>=0;i--)
	{
		if(y_r[i]!=0)
		{
			j=i;
			break;
		}		
	}
	for(i=j;i>=0;i--)
	{
		printf("%d", y_r[i]);
	}
	printf("=？\n");
}

int Prime_Omega(int N, int N_i, int *w_r)
{
	int k, i, r, h, N_p, O_r, N_k;
	k=1;
	r=0;
	i=0;
	h=2;
	N_p=N_i;
	while(81*N>N_p+1)
	{
		k=k+1;
		N_p=N_i*k;
	}
	k=-1;
	N_k=N_p+1;
	while(r==0)
	{
		k=k+1;
		N_p=N_k+N_i*k;
		r=Check_Prime(N_p);
	}
	while(h>1)
	{
		Find_Omega(N_p, N_i, w_r);
		h=Check_Omega(w_r, N_i);
		k=0;
		r=0;
		i=i+1;
		N_k=N_p;
		while(r==0)
		{
			k=k+1;
			N_p=N_k+N_i*k;
			r=Check_Prime(N_p);
		}
	}
	return N_p-N_i*k;
}

int Check_Prime(int N_p)
{
	int i;
	for(i=2;i<=sqrt(N_p);i++)
	{
		if(N_p%i==0)
		{
			return 0;
		}
	}
	return 1;
}

int Find_Omega(int N_p, int N_i, int *w_r)
{
	int i, j, k, l, O_i;
	i=2;
	k=0;
	while((k%N_p)!=1)
	{
		k=1;
		for(j=0;j<N_i;j++)
		{
			k=k*i;
			k=k%N_p;
		}
		i=i+1;
	}
	O_i=i-1;
	w_r[1]=O_i;
	l=O_i;
	for(i=2;i<N_i;i++)
	{
		l=(l*w_r[1])%N_p;
		w_r[i]=l;
	}	
}

int Check_Omega(int *w_r, int N_i)
{
	int i, k;
	k=0;
	for(i=0;i<N_i;i++)
	{
		if(w_r[i]==1)
		{
			k=k+1;
		}
	}
	return k;
}

int Bit_Reverse_Integer(int N_i, int *x_r, int *z_r)
{
	int group, size, i, j, k, n, p, number, power;
	group=1;
	size=N_i;
	number=N_i;
	
	p=0;

	while(number!=1)
	{
		number=number/2;
		p=p+1;
	}
	power=p;
	
	while(p>0)
	{
	for (i = 0; i<group; i++)
	{
		for(j=0;j<2;j++)
		{
			for(k=0;k<size/2;k++)
			{
				z_r[i*size+(j*size)/2+k] = x_r[i*size+j+2*k];
		    }
		}
	}
	for(n=0;n<N_i;n++)
	{
		x_r[n] = z_r[n];
	}

		group=group*2;
		size = size/2;
		p--;
	}
	return power;
}

int FFT(int *x_r, int *w_r, int N_i, int a, int prime)
{
	int n,p,k,q1,i;
	int t_r;
	n = 2;
	while(a>0)
	{
		for(k=0;k<n/2;k++)
		{
			for(p=k;p<N_i;p+=n)
			{
				q1 = p + n/2;
				if(N_i-(k*N_i/n)==N_i)
				{t_r = (w_r[0]*x_r[q1])%prime;}
				else
				{t_r = (w_r[N_i-(k*N_i/n)]*x_r[q1])%prime;}
				x_r[q1] = (x_r[p]+w_r[N_i/2]*t_r)%prime;
				x_r[p] = (x_r[p]+t_r)%prime;
			}
		}
	n = n*2;
	a--;
    }
}

int IFFT(int *x_r, int *w_r, int N_i, int a, int prime)
{
	int n,p,k,q1,i;
	int t_r;
	n = 2;
	while(a>0)
	{
		for(k=0;k<n/2;k++)
		{
			for(p=k;p<N_i;p+=n)
			{
				q1 = p + n/2;
				t_r = (w_r[k*N_i/n]*x_r[q1])%prime;
				x_r[q1] = (x_r[p]+w_r[N_i/2]*t_r)%prime;
				x_r[p] = (x_r[p]+t_r)%prime;
			}
		}
	n = n*2;
	a--;
    }
}

int result(int *x_r, int N_i)
{
	int i, j, h;
	for(i=0;i<N_i;i++)
	{
		j=0;
		h=x_r[i];
		x_r[i]=h%10;
		h=h/10;
		while(h>0)
		{
			j=j+1;
			x_r[i+j]=x_r[i+j]+h%10;
			h=h/10;
		}
	}
	for(i=N_i-1;i>=0;i--)
	{
		if(x_r[i]!=0)
		{
			j=i;
			break;
		}		
	}
	printf("結果為：");
	for(i=j;i>=0;i--)
	{
	printf("%d", x_r[i]);
	}
}
