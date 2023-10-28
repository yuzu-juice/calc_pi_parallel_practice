/**********
HPC
calc_pi_Gauss_Legendre

a0 = 1
b0 = 1/sq2
t0 = 1/4
p0 = 1

a = (a0+b0)/2
b = sqrt(a0*b0)
t = t0 - p0 * (a0-a)^2
p = 2*p0

pi = (a+b)^2 / (4*t)


コンパイル：gcc calpara.c -lgmp -fopenmp
**********/
#include <stdio.h>
#include <omp.h>
#include <gmp.h>

 //ループ回数
#define N 35

//変数の宣言
mpf_t a0, b0, t0, p0, aa;
mpf_t a, b, t, p;
int i, j;  //ループ変数
char *sf;
mp_exp_t exponent;
long keta = 10000000000; //求めたい桁数
int moji = 100; //1行の出力文字数


int main()
{
  FILE *fp;
  fp = fopen("pi100oku.txt", "w");
  //精度の指定
  mpf_set_default_prec(3.33*keta);

  
  //変数の初期化と代入
  mpf_init_set_ui(a0, 1);
  mpf_init(b0);
  mpf_sqrt_ui(b0, 2);
  mpf_ui_div(b0, 1, b0);
  mpf_init_set_d(t0, 0.25);
  mpf_init_set_ui(p0, 1); 
  mpf_init(a); 
  mpf_init(b);
  mpf_init(t);
  mpf_init(p);
  mpf_init(aa);

  //printf("%d\n", omp_get_max_threads());

  for (i=0; i<N; ++i)
  {
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        mpf_mul(b, a0, b0);
        mpf_sqrt(b, b);
        //#pragma omp critical
	//gmp_fprintf(stdout, "myid=%d, b=%Ff\n", omp_get_thread_num(), b);
      }

      #pragma omp section
      {
        mpf_mul_ui(p, p0, 2);
        //#pragma omp critical
    	//gmp_fprintf(stdout, "myid=%d, p=%Ff\n", omp_get_thread_num(), p);
      }
 
      #pragma omp section
      {
        mpf_add(a, a0, b0);
        mpf_div_ui(a, a, 2);
        mpf_sub(aa, a0, a);
        mpf_pow_ui(aa, aa, 2);
        mpf_mul(p0, p0, aa);
        mpf_sub(t, t0, p0);
        //#pragma omp critical
	//gmp_fprintf(stdout, "myid=%d, t=%Ff\n", omp_get_thread_num(), t);
	//#pragma omp critical
	//gmp_fprintf(stdout, "myid=%d, a=%Ff\n", omp_get_thread_num(), a);
      }
    }
    #pragma omp barrier

    mpf_set(a0, a);
    mpf_set(b0, b);
    mpf_set(t0, t);
    mpf_set(p0, p);
  }

  mpf_add(a, a, b);
  mpf_pow_ui(a, a, 2);
  mpf_mul_ui(t, t, 4);
  mpf_div(a, a, t);

  sf = mpf_get_str(NULL, &exponent, 10, 0, a); 
  
  //結果の出力
  gmp_fprintf(fp, "%s\n", "3.");
  for (i=0; i<keta/moji; ++i)
  {
    for (j=0; j<moji; ++j)
    {
      gmp_fprintf(fp, "%c", sf[i*moji+j+1]);
    }
    gmp_fprintf(fp, "\n");
  }

  fclose(fp);
  
  //終了処理
  mpf_clear(a0);
  mpf_clear(b0);
  mpf_clear(t0);
  mpf_clear(p0);
  mpf_clear(a);
  mpf_clear(aa);
  mpf_clear(b);
  mpf_clear(t);
  mpf_clear(p);
  
  return 0;
}
