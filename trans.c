/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */


/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }

}

void mini_trans(int startX, int startY, int A[32][32], int B[32][32])
{
    int blockSize = 8;
  // printf("doing block: %d, %d\n", startX, startY);
    int k, l, tmp;

    for (k = startX; k < (startX + blockSize); k++) {
        for (l = startY; l < (startY + blockSize); l++) {
            tmp = A[k][l];
            B[l][k] = tmp;
        }
    }

}
 /*
  * block - A blockwise transpose function
  */
 char block_desc[] = "block transpose";
 void blocka(int M, int N, int A[N][M], int B[M][N])
 {
     int rblock, cblock;
     int blockSize = (N == 32? 8:16);
     int blockN = blockSize;
     int blockM = blockSize;
     int a = 0;
     int t = 0;
     // printf("blockN: %d", blockN);
     for (rblock= 0; rblock < N; rblock += blockN) {
       if ((rblock + blockSize) > N) {
         blockN = N - rblock;
       } else {
         blockN = blockSize;
       }
       for (cblock = 0; cblock < M; cblock += blockM) {

         if ((cblock + blockSize) > M) {
           blockM = M - cblock;
         } else {
           blockM = blockSize;
         }
         //mini-block transpose
         //printf("%d x %d block at %d, %d\n", blockN, blockM, i, j);
         int r, c, tmp;
         for (c = cblock; c < (blockM+cblock); c++) {
             for (r = rblock; r < (blockN+rblock); r++) {
	      if( r != c)
	      {
                 tmp = A[c][r];
                 B[r][c] = tmp;
	      }
	     else
	     {
	     t = A[r][c];
	     a = r;
	     }
             }
             if(rblock == cblock)
	     {
	       B[a][a] = t;
	     }
         }
       }
     }
 }
 void blockb(int M, int N, int A[N][M], int B[M][N])
 {
     int rblock, cblock;
     int blockSize = (N == 32? 8:16);
     int blockN = blockSize;
     int blockM = blockSize;
     int a = 0;
     int t = 0;
     // printf("blockN: %d", blockN);
     for (rblock= 0; rblock < N; rblock += blockN) {
       if ((rblock + blockSize) > N) {
         blockN = N - rblock;
       } else {
         blockN = blockSize;
       }
       for (cblock = 0; cblock < M; cblock += blockM) {

         if ((cblock + blockSize) > M) {
           blockM = M - cblock;
         } else {
           blockM = blockSize;
         }
         //mini-block transpose
         //printf("%d x %d block at %d, %d\n", blockN, blockM, i, j);
         int r, c, tmp;
         for (c = cblock; c < (blockM+cblock); c++) {
             for (r = rblock; r < (blockN+rblock); r++) {
	      if( r != c)
	      {
                 tmp = A[r][c];
                 B[c][r] = tmp;
	      }
	     else
	     {
	     t = A[r][c];
	     a = r;
	     }
             }
             if(rblock == cblock)
	     {
	       B[a][a] = t;
	     }
         }
       }
     }
 }


 void block32of64(int startX, int startY, int A[64][64], int B[64][64])
 {
     int i, j;
     // int blockSize = 8;
     int blockN = 4;
     int blockM = 4;
     // printf("blockN: %d", blockN);
     for (i= startX; i < (startX + 32); i += blockN) {
       // if ((i + blockSize) > N) {
         // blockN = N - i;
       // } else {
         // blockN = 8;
       // }
       for (j = startY; j < (startY + 32); j += blockM) {

         // if ((j + blockSize) > M) {
           // blockM = M - j;
         // } else {
           // blockM = 8;
         // }
         //mini-block transpose
         printf("-%d x %d block at %d, %d\n", blockN, blockM, i, j);
         int k, l, tmp;
         for (l = j; l < (j + blockM); l++) {
             for (k = i; k < (i + blockN); k++) {
                 tmp = A[k][l];
                 B[l][k] = tmp;
             }
         }
       }
     }
 }
 void block64(int M, int N, int A[N][M], int B[M][N]) {
   /*int blocksize = 4;
   int rblock;
   int cblock;
   int t0, t1, t2, t3, t4, t5;
   
   
   for(rblock = 0; rblock < N; rblock += blocksize)
   {
     for(cblock = 0; cblock < M; cblock += blocksize)
     {
       for(int r = 0; r < 2; r++)
       {
	 for(int c = 0; c < 2; c++)
	 {
	   B[cblock + c][rblock + r] = A[rblock+r][cblock+c]; 
	 }
       }
       
     }
   }*/
   
   int colRun;
   int rowRun;
   for(colRun=0; colRun<64; colRun+=8 ){
       	for(rowRun=0; rowRun<64; rowRun+=8 ){

	  int a0, a1, a2, a3, a4, a5, a6, a7;
	  int k;
        /* 
        The first loop is for a 4 collumn 8 row of A. In this loop we  use the supporting variables to store all elements of a row. Then we try to transpose with the right positions in B. However, we are doing on a 4row x 8 collumn B to avoid cache miss so not all elements will be transposed correctly. For example A[3][5] cannot be transposed to B[5][3]. Thus, elements which can't be transposed correctly will be stored in another location for later use.
        */ 
        	for(k=0; k<4; k++){ 
        		a0 = A[colRun+k][rowRun+0];
        		a1 = A[colRun+k][rowRun+1];
        		a2 = A[colRun+k][rowRun+2];
        		a3 = A[colRun+k][rowRun+3];
        		a4 = A[colRun+k][rowRun+4];
        		a5 = A[colRun+k][rowRun+5];
        		a6 = A[colRun+k][rowRun+6];
        		a7 = A[colRun+k][rowRun+7];

        		// In the code, I comment "Good job" for the elements that are transposed correctly. "Later use" for later assignment
        		B[rowRun+0][colRun+k+0] = a0;   // good job
        		B[rowRun+0][colRun+k+4] = a5;	// later use
        		B[rowRun+1][colRun+k+0] = a1;	// good job
        		B[rowRun+1][colRun+k+4] = a6;	//later use
        		B[rowRun+2][colRun+k+0] = a2;	// good job
        		B[rowRun+2][colRun+k+4] = a7;	//later use
        		B[rowRun+3][colRun+k+0] = a3;	// good job
        		B[rowRun+3][colRun+k+4] = a4;	// later use
        	}


            /* Part B, moving sub-matrix b to sub-matrix c
             * and moving A->B for sub-matrix b and move matrix d
             */
            /*
            Now that we have dealt with the first 4 col 8 arrow of A. The next job to deal with the "later use" assignment above. The "later use" assignments that we did above have taken a lot of places, so we need to bring these elements to their right positions. 
             */
        	a0 = A[colRun+4][rowRun+4];
        	a1 = A[colRun+5][rowRun+4];
        	a2 = A[colRun+6][rowRun+4];
        	a3 = A[colRun+7][rowRun+4];
        	a4 = A[colRun+4][rowRun+3];
        	a5 = A[colRun+5][rowRun+3];
        	a6 = A[colRun+6][rowRun+3];
        	a7 = A[colRun+7][rowRun+3];


        	B[rowRun+4][colRun+0] = B[rowRun+3][colRun+4];  // B[4][0] = a4 = A[0][4] For example
        	B[rowRun+4][colRun+4] = a0;  // B[4][4] = A[4][4] For example
        	B[rowRun+3][colRun+4] = a4;
        	B[rowRun+4][colRun+1] = B[rowRun+3][colRun+5];
        	B[rowRun+4][colRun+5] = a1;
        	B[rowRun+3][colRun+5] = a5;
        	B[rowRun+4][colRun+2] = B[rowRun+3][colRun+6];
        	B[rowRun+4][colRun+6] = a2;
        	B[rowRun+3][colRun+6] = a6;
        	B[rowRun+4][colRun+3] = B[rowRun+3][colRun+7];
        	B[rowRun+4][colRun+7] = a3;
        	B[rowRun+3][colRun+7] = a7;

        	// this loops deal with the the remaning elements .
        	for(k=0;k<3;k++){


        		a0 = A[colRun+4][rowRun+5+k];
        		a1 = A[colRun+5][rowRun+5+k];
        		a2 = A[colRun+6][rowRun+5+k];
        		a3 = A[colRun+7][rowRun+5+k];
        		a4 = A[colRun+4][rowRun+k];
        		a5 = A[colRun+5][rowRun+k];
        		a6 = A[colRun+6][rowRun+k];
        		a7 = A[colRun+7][rowRun+k];


        		B[rowRun+5+k][colRun+0] = B[rowRun+k][colRun+4];
        		B[rowRun+5+k][colRun+4] = a0;
        		B[rowRun+k][colRun+4] = a4;
        		B[rowRun+5+k][colRun+1] = B[rowRun+k][colRun+5];
        		B[rowRun+5+k][colRun+5] = a1;
        		B[rowRun+k][colRun+5] = a5;
        		B[rowRun+5+k][colRun+2] = B[rowRun+k][colRun+6];
        		B[rowRun+5+k][colRun+6] = a2;
        		B[rowRun+k][colRun+6] = a6;
        		B[rowRun+5+k][colRun+3] = B[rowRun+k][colRun+7];
        		B[rowRun+5+k][colRun+7] = a3;
        		B[rowRun+k][colRun+7] = a7;


        	}


        }
} 
 }


 /*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
 char transpose_submit_desc[] = "Transpose submission";
 void transpose_submit(int M, int N, int A[N][M], int B[M][N])
 {
   if(N==64)
    block64(M, N, A, B);
   else if(N==32)
     blocka(M,N,A,B);
   else
     blockb(M,N,A,B);
 }

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);


    //registerTransFunction(block, block_desc);
    /* Register any additional transpose functions */
    //registerTransFunction(trans, trans_desc);

}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                printf("error at: %d, %d\n", i, j);
                return 0;
            }
        }
    }
    return 1;
}
