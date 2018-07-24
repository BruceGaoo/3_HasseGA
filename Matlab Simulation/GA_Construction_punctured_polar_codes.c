//�����ص�����
//ÿ�ζ����¿����ڴ棬��ֹ�����ڴ湲��֮�������
//����GA���������׼�����
//MaxValue�����úͶ��ַ��и��ֳ��������ö���һ������
//���ַ��е���ʼ�㣨y_0,y_1���Ĳ�Ӧ�ô�����ֹ����,����׼
//�����ʺܵ͵�ʱ����Ҫ��һ�����߶��ַ��ľ��ȣ�����׼�����������Ϣλ������ʶ�һ����
#include "mex.h"
#include "matrix.h"
#include "math.h"
#define MaxValue 100.0
//////////////////////////////////////////////////////////////////////////
//	
//	SC Decoder
//	
//
//	N		code length; n = logN/log2;
//	LR_msg	likelihood ratio message of bit-reversed codeword
//	frozen	position of frozen set
//	Lout	likelihood ratio message of decoded bit
//	hd_dec	decoded word
//	
//////////////////////////////////////////////////////////////////////////
double psi(double x)
{ double f;
  if (x>=10)
  f = sqrt(3.1415926/x)*(1-10/(7*x))*exp(-x/4);
  else
  f = exp(-0.4527*pow(x,0.86)+0.0218); 
  return f;
}

double psi_inv(double y)
{ double x_i,y_i,y_1,y_2;
  double x_2=MaxValue;
  double x_1=10;
  double x;

  if (y > 0.0385) //������
      x = pow(((0.0218-log(y))/0.4527),(1/0.86)); 
  else //���ַ����
     if (y<2.4264e-12)//��ֵΪ��x=MaxValueʱ��x��ȡֵ
      x = MaxValue;
     else
    {      
      do
 {   x_i = (x_1+x_2)/2;
     y_i = sqrt(3.1415956/x_i)*(1-10/(7*x_i))*exp(-x_i/4);
     if (y_i > y) 
     {x_1 = x_i;    
     }
     else
     {x_2 = x_i; 
     }
     
 }while(fabs(y_i - y)>1e-15);
  x = x_i;
	 }
  return x;
}


long double node_f_SPA(long double a, long double b)//�����phi����
{
	long double f,psi_a,psi_b;
	long  double ta, tb, y;
if (a == 0 || b == 0)
    f = 0;
else  
{
    
	psi_a = psi(a);
	psi_b = psi(b);
    y = 1 - (1-psi_a)*(1-psi_b);
	f = psi_inv(y);
}  
	if (f>MaxValue)
	{
		f = MaxValue;
	}
	else if (f<-MaxValue)
	{
		f = -MaxValue;
	}

	return f;
}


//////////////////////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////////////////////
int CalPos(int Lambda, int Phi, int Beta)
{
	int result;
	int tmp;

	tmp = (int)pow(2.0, (long double)Lambda);
	result = Beta * tmp + Phi;

	return result;
}

void recusivelyUpdateB(int Lambda, int Phi, int** bb,int n)
{
	int Beta;
	int Psi;
	int index;
	int p1, p2, p3, p4;

	if ((Phi % 2) == 1)
	{
		Psi = (int)(Phi / 2);

		index = (int)pow(2.0, (long double)(n - Lambda));

		for (Beta = 0; Beta<index; Beta++)
		{
			p1 = CalPos((Lambda - 1), Psi, (2 * Beta));
			p2 = CalPos(Lambda, (Phi - 1), Beta);
			p3 = CalPos(Lambda, Phi, Beta);
			p4 = CalPos((Lambda - 1), Psi, (2 * Beta + 1));

			bb[Lambda - 1][p1] = bb[Lambda][p2] ^ bb[Lambda][p3];
			bb[Lambda - 1][p4] = bb[Lambda][p3];
		}


		if (Psi % 2 == 1)
		{
			recusivelyUpdateB((Lambda - 1), Psi, bb,n);
		}
	}
}


//////////////////////////////////////////////////////////////////////////
//	
//	x	decide
//	n	number of
//	y	��y[0]~y[n-1]from low to high
//////////////////////////////////////////////////////////////////////////
void d2b(int x, int* y,int n)
{
	int i, index;

	for (i = n - 1; i >= 0; i--)
	{
		index = (int)pow(2.0, (long double)i);
		y[i] = (int)(x / index);
		x = x - y[i] * index;
	}
}

//void SC_decoding_SPA(long double* LLR_msg, long double* frozen_set, long double *hd_dec,int N, int n, int fr_n)
void SC_decoding_SPA(long double* LLR_msg, long double *hd_dec,int N, int n)

{
	//static short flag;
	//static int** PartialSum;
	//static double** mid_LLR;	//	
	//static int** func_index;	//	
    
//	short flag;
	int** PartialSum;
	long double** mid_LLR;	//	
	int** func_index;	//	
    
	int i, j, k;
	int index;
	
	int p2;
	int s;
	int *p1;
	//int err_fro_counter=0;
	
	
	int cnt_frozen = 0;			//	

	int phi;
	int tmp;


	//////////////////////////////////////////////////////////////////////////
	//	//////////////////////////////////////////////////
	//if (flag == 0)//���ֿ����ڴ�
	//{
		PartialSum = (int**)malloc(sizeof(int*) * (n+1));//initial partialsum matrix (n+1)*N
		for (i=0; i<=n; i++)
		{
			PartialSum[i] = (int *)malloc(sizeof(int) * N);
		}

		mid_LLR = (long double**)malloc(sizeof(long double*) * (n+1));//initial mid_LLR matrix (n+1)*N
		for (i=0; i<=n; i++)
		{
			mid_LLR[i] = (long double*)malloc(sizeof(long double) * N);
		}

		func_index = (int**)malloc(sizeof(int*) * N); //initial func_index matrix N*n
		for (i=0; i<N; i++)
		{
			func_index[i] = (int*)malloc(sizeof(int) * n);
		}

		//flag = 1;
	//}


	//////////////////////////////////////////////////////////////////////////
	//	
	//////////////////////////////////////////////////////////////////////////
	for (i=0; i<N; i++)
	{
		mid_LLR[n][i] = LLR_msg[i];//mid_LLR �ĵ�n�е�ֵΪ�ŵ����յ���LLR
	}

	p1 = (int*)malloc(sizeof(int)* n);
	for (i=0; i<n; i++)
	{
		p1[i] = 0;
	}

	//	
	for (i=0; i<N; i++)				//	˳�������i����Ϣλ
	{
		d2b(i, func_index[i],n);		//	

		for (j=n-1; j>=0; j--)		//	
		{
			if (((i>0)&&(func_index[i][j]^func_index[i-1][j])) || i==0)//������������1��˵���ò�ڵ��Ѿ���������ˣ�����ֱ������
			{
				index = (int) pow(2.0,(long double)j);	//	��j��һ����Ҫ����2^j���ڵ�

				if (func_index[i][j]==1)			//	
				{
					for (k=0; k<index; k++)//Ӧ����ÿ��ʱ϶һ�����Լ�����ٸ��ڵ����˼  //	��j��һ����Ҫ����2^j���ڵ�
					{
						p2 = p1[j] + k - index;

						phi = i / index - 1;

						tmp = CalPos((n-j), phi, k);

						s = PartialSum[n-j][tmp];

						mid_LLR[j][(p1[j])] = mid_LLR[j+1][p2] + mid_LLR[j+1][p2+1];//��ֵ��Ϊ������������
	if (mid_LLR[j][(p1[j])]>MaxValue)
	{
		mid_LLR[j][(p1[j])] = MaxValue;
	}
	else if (mid_LLR[j][(p1[j])]<-MaxValue)
	{
		mid_LLR[j][(p1[j])] = -MaxValue;
	}



						p1[j] ++;
					}
				}
				else if (func_index[i][j]==0)		//	
				{
					for (k=0; k<index; k++)
					{
						p2 = p1[j] + k;             //�ڴ˴�k��ʵʼ�յ���p1[j]

						mid_LLR[j][(p1[j])] = node_f_SPA(mid_LLR[j+1][p2], mid_LLR[j+1][p2+1]);

						p1[j] ++;
					}
				}
			}
		}//	end for loop l

		//	
		//if (i==frozen_set[cnt_frozen])

	

				PartialSum[n][i] = 0;//partialsum matrix ȫ0����

		

		if ((i%2)==1)
		{
			recusivelyUpdateB(n, i, PartialSum,n);
		}
	}//	end for loop i


	for (i=0; i<N; i++)
	{
		hd_dec[i] = mid_LLR[0][i];
	}
	//*outNum=err_fro_counter;
  for (i=0; i<=n; i++)
		{
			free(PartialSum[i]);
            free(mid_LLR[i]);
		}  
    
free(PartialSum);
free(mid_LLR);

  for (i=0; i<N; i++)
		{
			free(func_index[i]);
		}

    free(func_index);
    free(p1);
	

}


//���룺ÿ��λ�õľ�ֵ�����
//����׵�λ�õľ�ֵΪ0������0
//����λ�þ�ֵΪ����2
//����ȫ���д���Ҳ�첻�ˣ���ΪҲ��һ��һ�����㷨
//ֻ��Ҫ�����ֵ�ͺ��ˣ���
void mexFunction(int output_size, mxArray *output[], int input_size, const mxArray *input[])
{
	double *LLR_msg = mxGetPr(input[0]); //	����ľ�ֵ������׵�λ��Ϊ0
//	double *frozen_set = mxGetPr(input[1]);
	double *dN = mxGetPr(input[1]);  //
	double *dn = mxGetPr(input[2]);  //
 //  double *nframe = mxGetPr(input[4]);
 //   double *re_i = mxGetPr(input[5]);
//	double *frozen_n = mxGetPr(input[4]);
	//////////////////////////////////////
	int N = (int)(*dN);
	int n = (int)(*dn);
 //   int n_frame = (int)(*nframe);
 //   int rei = (int)(*re_i);
//	int fr_n = (int)(*frozen_n);
	int i;
	double *outData;
	double *hd_dec;//����ľ�ֵ
	//double *outNum;
	//double *Error_Pos;
	/////////////////////////////////////
   // if (n_frame == 97 && rei == 137)
   //     n_frame = n_frame;
    
	output[0] = mxCreateDoubleMatrix(N, 1, mxREAL);

	outData = mxGetPr(output[0]);
	
	hd_dec = (double*)malloc(sizeof(double)* N);
	//output[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	//outNum = mxGetPr(output[1]);
	//output[2] = mxCreateDoubleMatrix(N, 1, mxREAL);
	//Error_Pos = mxGetPr(output[2]);
	//for (i=0;i<N;i++)
	//{
	//	Error_Pos[i] = 0;
	//}
//	SC_decoding_SPA(LLR_msg, frozen_set, hd_dec,N, n, fr_n);
    SC_decoding_SPA(LLR_msg, hd_dec,N, n);
	for (i = 0; i<N; i++)
	{
		outData[i] = (double)hd_dec[i];
	}
free(hd_dec);	
}