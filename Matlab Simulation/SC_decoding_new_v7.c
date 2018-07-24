//仅返回单参数
//每次都重新开辟内存，防止出现内存共享之类的问题
//针对v1进行精简，去除掉测试用的代码
//考虑了f运算时，数值过小问题,采用直接判定小于某个数的方式
//针对上述判决方法做了改进
#include "mex.h"
#include "matrix.h"
#include "math.h"
#define MaxValue 32.0
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
long double node_f_SPA(long double a, long double b)
{
	long double f;
	long  double ta, tb;

	ta = tanh(a / 2.0);
	tb = tanh(b / 2.0);

	f = log((1.0 + ta*tb) / (1.0 - ta*tb));
	if (fabs(f) < 1e-14 && (ta == 0 || tb == 0))        
        if (f > 0)
        f = 1e-14;    //就相当于记录个符号
      else
        f = -1e-14;
    
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

long double node_g_SPA(long double a, long double b, int s)
{
	long double g;

	if ((s != 0) && (s != 1))
	{
		printf("Error");
	}

	g = b + (1 - 2 * s) * a;

	if (g>MaxValue)
	{
		g = MaxValue;
	}
	else if (g<-MaxValue)
	{
		g = -MaxValue;
	}

	return g;
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
//	y	，y[0]~y[n-1]from low to high
//////////////////////////////////////////////////////////////////////////
void d2b(int x, int* y,int n) 
//实际效果是得到一个翻转后的二进制序列。ig:n=3 x=5 100 001。
{
	int i, index;

	for (i = n - 1; i >= 0; i--)
	{
		index = (int)pow(2.0, (long double)i);
		y[i] = (int)(x / index);
		x = x - y[i] * index;
	}
}

void SC_decoding_SPA(long double* LLR_msg, long double* frozen_set, long double *hd_dec,int N, int n, int fr_n)
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
	//if (flag == 0)//各种开辟内存
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
		mid_LLR[n][i] = LLR_msg[i];//mid_LLR 的第n行的值为信道接收到的LLR
	}

	p1 = (int*)malloc(sizeof(int)* n); 
	for (i=0; i<n; i++)
	{
		p1[i] = 0;
	}

	//	
	for (i=0; i<N; i++)				//	顺序译码第i个信息位
	{
		d2b(i, func_index[i],n);		//	

		for (j=n-1; j>=0; j--)		//	
		{
			if (((i>0)&&(func_index[i][j]^func_index[i-1][j])) || i==0)//如果异或结果不是1，说明该层节点已经被计算过了，所以直接跳过
			{
				index = (int) pow(2.0,(long double)j);	//	第j层一次需要计算2^j个节点

				if (func_index[i][j]==1)			//	
				{
					for (k=0; k<index; k++)//应该是每个时隙一共可以计算多少个节点的意思  //	第j层一次需要计算2^j个节点
					{
						p2 = p1[j] + k - index;

						phi = i / index - 1;

						tmp = CalPos((n-j), phi, k);

						s = PartialSum[n-j][tmp];

						mid_LLR[j][(p1[j])] = node_g_SPA(mid_LLR[j+1][p2], mid_LLR[j+1][p2+1], s);

						p1[j] ++;
					}
				}
				else if (func_index[i][j]==0)		//	
				{
					for (k=0; k<index; k++)
					{
						p2 = p1[j] + k;             //在此处k其实始终等于p1[j]

						mid_LLR[j][(p1[j])] = node_f_SPA(mid_LLR[j+1][p2], mid_LLR[j+1][p2+1]);

						p1[j] ++;
					}
				}
			}
		}//	end for loop l

		//	
		//if (i==frozen_set[cnt_frozen])
		if (frozen_set != NULL && i==frozen_set[cnt_frozen] && cnt_frozen < fr_n)
        {
			PartialSum[n][i] = 0;
			cnt_frozen ++;
			//if (mid_LLR[0][i]<0)
			//{
			//	err_fro_counter++;
			//	Error_Pos[i]=1;//记录冻结集译码错误的位置
			//}

		}
		else 
		{
			if (mid_LLR[0][i]>=0)
			{
				PartialSum[n][i] = 0;//partialsum matrix 矩阵的第n行表示信息位
			}
			else
			{
				PartialSum[n][i] = 1;
			}

		}

		if ((i%2)==1)
		{
			recusivelyUpdateB(n, i, PartialSum,n);
		}
	}//	end for loop i


	for (i=0; i<N; i++)
	{
		hd_dec[i] = PartialSum[n][i];
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





void mexFunction(int output_size, mxArray *output[], int input_size, const mxArray *input[])
{
	double *LLR_msg = mxGetPr(input[0]); //	记录各层的似然比信息
	double *frozen_set = mxGetPr(input[1]);
	double *dN = mxGetPr(input[2]);  //
	double *dn = mxGetPr(input[3]);  //
 //  double *nframe = mxGetPr(input[4]);
 //   double *re_i = mxGetPr(input[5]);
	double *frozen_n = mxGetPr(input[4]);
	//////////////////////////////////////
	int N = (int)(*dN);
	int n = (int)(*dn);
 //   int n_frame = (int)(*nframe);
 //   int rei = (int)(*re_i);
	int fr_n = (int)(*frozen_n);
	int i;
	double *outData;
	double *hd_dec;
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
    SC_decoding_SPA(LLR_msg, frozen_set, hd_dec,N, n, fr_n);
	for (i = 0; i<N; i++)
	{
		outData[i] = (double)hd_dec[i];
	}
free(hd_dec);	
}