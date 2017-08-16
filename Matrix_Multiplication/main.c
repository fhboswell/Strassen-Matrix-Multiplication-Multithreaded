//
//  main.c
//  Martix
//
//  Created by Franklin Henry Boswell on 7/20/17.
//  Copyright © 2017 Franklin Henry Boswell. All rights reserved.
//

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <pthread.h>






void square_MAT_Mul_recursion_multi(int *MATA, int *MATB, int*MATC, int N, int O_width);
void strassen_recursion_multi_plus(int *MATA, int *MATB, int*MAT, int  N, int cutoff);
void strassen_recursion_multi(int *MATA, int *MATB, int*MAT, int  N, int cutoff);

void square_MAT_Mul_recursion(int *MATA, int *MATB, int*MATC, int N, int O_width, int cutoff);
void strassen_recursion(int *MATA, int *MATB, int*MAT, int  N, int cutoff);


void naive_IKJ_Square(int *MAT1, int *MAT2, int *MAT4, int N);


void addition(int *MAT1, int *MAT2, int*MAT, int N);
void subtraction(int *MAT1, int *MAT2, int*MAT, int N);

void addition_SUB_MAT(int *MAT1, int *MAT2, int*MAT, int N);
void subtraction_SUB_MAT(int *MAT1, int *MAT2, int*MAT, int N);

void addition_SUB_MAT_OUT(int *MAT1, int *MAT2, int*MAT, int N);
void subtraction_SUB_MAT_OUT(int *MAT1, int *MAT2, int*MAT, int N);

//tests the acuracy of the implementation against the naive algorithm implementation
void proof (int *MATA, int *MATB, int*MAT, int  N);

void printMAT(int *MAT, int N);
void print_Sub_MAT(int *MAT, int O);

void test();
double get_time();
void reshape(int *MATAA, int *MATBB, int N);


//functions for multi threading

void *addition_ThreadRunner(void *arg);
void *subtraction_ThreadRunner(void *arg);

void *makeC11(void *arg);
void *makeC12(void *arg);
void *makeC21(void *arg);
void *makeC22(void *arg);


void *mulRunner(void *arg);

void *naive(void *arg);

//structs for passing data

struct threadData{
    int *subMAT1;
    int *subMAT2;
    int *outMAT;
    int N;
    int ThreadNum;//for demonstration
};

struct threadData_C{
    
    int *P1 ;
    int *P2;
    int *P3;
    int *P4;
    int *P5;
    int *P6;
    int *P7;
    
    int *A11A;
    int *B22B;
    int *A22A;
    int *B11B;
    
    int *C11;
    int *C12;
    int *C21;
    int *C22;
    
    
    int N;
    int ThreadNum;
    
    
    
};



struct threadData_SM{
    
    
    int *A1;
    int *A2;
    int *B1;
    int *B2;
    
    int *C;
    int O_width;
    
    int N;
    int ThreadNum;
    
    
    
};

struct threadData_P{
    
    
    int *A;
    int *B;
    
    
    int *C;
    
    
    int N;
    int ThreadNum;
    
    
    
};



int main(int argc, const char * argv[]) {
    
    
    test();
    
    
    
}




void strassen_recursion_multi_plus(int *MATA, int *MATB, int*MAT, int N, int cutoff){
    
    
    
    
    int *S1 = (int *)malloc(N * N * sizeof(int));
    int *S2 = (int *)malloc(N * N * sizeof(int));
    int *S3 = (int *)malloc(N * N * sizeof(int));
    int *S4 = (int *)malloc(N * N * sizeof(int));
    int *S5 = (int *)malloc(N * N * sizeof(int));
    int *S6 = (int *)malloc(N * N * sizeof(int));
    int *S7 = (int *)malloc(N * N * sizeof(int));
    int *S8 = (int *)malloc(N * N * sizeof(int));
    int *S9 = (int *)malloc(N * N * sizeof(int));
    int *S10 = (int *)malloc(N * N * sizeof(int));
    
    
    int *B11 = (MATB + 0*N + 0);
    int *B12 = (MATB + 0*N + ((N/2)));
    int *B21 = (MATB + N*((N/2)) + 0);
    int *B22 =  (MATB + N*((N/2))+((N/2)));
    
    
    
    int *A11 = (MATA + 0*N + 0);
    int *A12 = (MATA + 0*N + ((N/2)));
    int *A21 = (MATA + N*((N/2)) + 0);
    int *A22 =  (MATA + N*((N/2))+((N/2)));
    
    int *C11 = (MAT + 0*N + 0);
    int *C12 = (MAT + 0*N + ((N/2)));
    int *C21 = (MAT + N*((N/2)) + 0);
    int *C22 =  (MAT + N*((N/2))+((N/2)));
    
    
    
    N = N/2;
    
    struct threadData adds[6];
    struct threadData subs[4];
    
    adds[0] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A12, .outMAT = S2, .N = N, .ThreadNum = 0 };
    adds[1] =(struct threadData){ .subMAT1 = A21, .subMAT2 = A22, .outMAT = S3, .N = N, .ThreadNum = 1};
    adds[2] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A22, .outMAT = S5, .N = N, .ThreadNum = 2 };
    adds[3] =(struct threadData){ .subMAT1 = B11, .subMAT2 = B22, .outMAT = S6, .N = N,.ThreadNum = 3 };
    adds[4] =(struct threadData){ .subMAT1 = B21, .subMAT2 = B22, .outMAT = S8, .N = N, .ThreadNum = 4 };
    adds[5] =(struct threadData){ .subMAT1 = B11, .subMAT2 = B12, .outMAT = S10, .N = N, .ThreadNum = 5 };
    
    subs[0] =(struct threadData){ .subMAT1 = B12, .subMAT2 = B22, .outMAT = S1, .N = N, .ThreadNum = 6 };
    subs[1] =(struct threadData){ .subMAT1 = B21, .subMAT2 = B11, .outMAT = S4, .N = N, .ThreadNum = 7};
    subs[2] =(struct threadData){ .subMAT1 = A12, .subMAT2 = A22, .outMAT = S7, .N = N, .ThreadNum = 8 };
    subs[3] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A21, .outMAT = S9, .N = N,.ThreadNum = 9};
    
    
    
    
    pthread_t adders[6];
    pthread_t subtractors[4];
    
    for (int i = 0; i < 6; i++){
        
      		pthread_create(&adders[i], NULL, addition_ThreadRunner, &adds[i]);
        
        
    }
    for (int j = 0; j < 4; j++){
        
      		pthread_create(&subtractors[j], NULL, subtraction_ThreadRunner, &subs[j]);
        
        
    }
    
    
    
    
    for (int k = 0; k < 6; k++)
    {
        pthread_join(adders[k], NULL);
        
        
    }
    
    for (int l = 0; l < 4; l++)
    {
      		pthread_join(subtractors[l], NULL);
        
        
    }
    
    
    if(N > 1){
        int *P1 = (int *)calloc(N * N, sizeof(int));
        int *P2 = (int *)calloc(N * N, sizeof(int));
        int *P3 = (int *)calloc(N * N, sizeof(int));
        int *P4 = (int *)calloc(N * N, sizeof(int));
        int *P5 = (int *)calloc(N * N, sizeof(int));
        int *P6 = (int *)calloc(N * N, sizeof(int));
        int *P7 = (int *)calloc(N * N, sizeof(int));
        
        int *A11A = (int *)malloc(N * N * sizeof(int));
        int *B22B = (int *)malloc(N * N * sizeof(int));
        int *A22A = (int *)malloc(N * N * sizeof(int));
        int *B11B = (int *)malloc(N * N * sizeof(int));
        
        reshape(A11, A11A, N);
        reshape(B11, B11B, N);
        reshape(B22, B22B, N);
        reshape(A22, A22A, N);
        
        if(N  <= cutoff){
            //printMAT(P1, N);
            struct threadData_P make_P[7];
            pthread_t make_P_th[7];
            
            make_P[0] =(struct threadData_P){ .A = A11A, .B = S1, .C = P1, .N = N, .ThreadNum = 0 };
            make_P[1] =(struct threadData_P){ .A = S2, .B = B22B, .C = P2, .N = N, .ThreadNum = 1 };
            make_P[2] =(struct threadData_P){ .A = S3, .B = B11B, .C = P3, .N = N, .ThreadNum = 2 };
            make_P[3] =(struct threadData_P){ .A = A22A, .B = S4, .C = P4, .N = N, .ThreadNum = 3 };
            make_P[4] =(struct threadData_P){ .A = S5, .B = S6, .C = P5, .N = N, .ThreadNum = 4 };
            make_P[5] =(struct threadData_P){ .A = S7, .B = S8, .C = P6, .N = N, .ThreadNum = 5 };
            make_P[6] =(struct threadData_P){ .A = S9, .B = S10, .C = P7, .N = N, .ThreadNum = 6 };
            /*
             naive_IKJ_Square(A11A, S1, P1, N);
             naive_IKJ_Square(S2, B22B, P2, N);
             naive_IKJ_Square(S3, B11B, P3, N);
             naive_IKJ_Square(A22A, S4, P4, N);
             naive_IKJ_Square(S5, S6, P5, N);
             naive_IKJ_Square(S7, S8, P6, N);
             naive_IKJ_Square(S9, S10, P7, N);
             */
            
            for (int j = 0; j < 7; j++){
                
                pthread_create(&make_P_th[j], NULL, naive, &make_P[j]);
                
                
            }
            
            
            
            
            for (int k = 0; k < 7; k++)
            {
                pthread_join(make_P_th[k], NULL);
                
                
            }
            //printMAT(P1,N);
            
            
            
            
        }else{
            
            strassen_recursion_multi_plus(A11A, S1, P1, N, cutoff);
            strassen_recursion_multi_plus(S2, B22B, P2, N, cutoff);
            strassen_recursion_multi_plus(S3, B11B, P3, N, cutoff);
            strassen_recursion_multi_plus(A22A, S4, P4, N, cutoff);
            strassen_recursion_multi_plus(S5, S6, P5, N, cutoff);
            strassen_recursion_multi_plus(S7, S8, P6, N, cutoff);
            strassen_recursion_multi_plus(S9, S10, P7, N, cutoff);
        }
        
        
        struct threadData_C make_C =(struct threadData_C){ .P1 = P1, .P2 = P2, .P3 = P3,.P4 = P4, .P5 = P5, .P6 = P6, .P7 = P7, .N = N, .A11A = A11A, .B22B = B22B, .A22A= A22A, .B11B =B11B, .C11=C11, .C12 = C12, .C21 = C21, .C22 = C22, .ThreadNum = 0 };
        
        pthread_t Make_C_SUBs[4];
        
        
        pthread_create(&Make_C_SUBs[0], NULL, makeC11, &make_C);
        pthread_create(&Make_C_SUBs[1], NULL, makeC12, &make_C);
        pthread_create(&Make_C_SUBs[2], NULL, makeC21, &make_C);
        pthread_create(&Make_C_SUBs[3], NULL, makeC22, &make_C);
        
        pthread_join(Make_C_SUBs[0], NULL);
        pthread_join(Make_C_SUBs[1], NULL);
        pthread_join(Make_C_SUBs[2], NULL);
        pthread_join(Make_C_SUBs[3], NULL);
        
        
        
        
        free(P1);
        free(P2);
        free(P3);
        free(P4);
        free(P5);
        free(P6);
        free(P7);
        
        free(A11A);
        free(B22B);
        free(A22A);
        free(B11B);
        
        
        
        
        
        
    }else{
        
        int P1 = (*A11) * (*S1);
        int P2 = (*S2) * (*B22);
        int P3 = (*S3) * (*B11);
        int P4 = (*A22) * (*S4);
        int P5 = (*S5) * (*S6);
        int P6 = (*S7) * (*S8);
        int P7 = (*S9) * (*S10);
        
        *(MAT) = P5 + P4 - P2 + P6;
        *(MAT + 1) = P1 + P2;
        *(MAT + 2) = P3 + P4;
        *(MAT + 3) = P5 + P1 - P3 - P7;
        
        
    }
    
    free(S1);
    free(S2);
    free(S3);
    free(S4);
    free(S5);
    free(S6);
    free(S7);
    free(S8);
    free(S9);
    free(S10);
}
void strassen_recursion_multi(int *MATA, int *MATB, int*MAT, int N,  int cutoff){
    
    
    
    
    int *S1 = (int *)malloc(N * N * sizeof(int));
    int *S2 = (int *)malloc(N * N * sizeof(int));
    int *S3 = (int *)malloc(N * N * sizeof(int));
    int *S4 = (int *)malloc(N * N * sizeof(int));
    int *S5 = (int *)malloc(N * N * sizeof(int));
    int *S6 = (int *)malloc(N * N * sizeof(int));
    int *S7 = (int *)malloc(N * N * sizeof(int));
    int *S8 = (int *)malloc(N * N * sizeof(int));
    int *S9 = (int *)malloc(N * N * sizeof(int));
    int *S10 = (int *)malloc(N * N * sizeof(int));
    
    
    int *B11 = (MATB + 0*N + 0);
    int *B12 = (MATB + 0*N + ((N/2)));
    int *B21 = (MATB + N*((N/2)) + 0);
    int *B22 =  (MATB + N*((N/2))+((N/2)));
    
    
    
    int *A11 = (MATA + 0*N + 0);
    int *A12 = (MATA + 0*N + ((N/2)));
    int *A21 = (MATA + N*((N/2)) + 0);
    int *A22 =  (MATA + N*((N/2))+((N/2)));
    
    int *C11 = (MAT + 0*N + 0);
    int *C12 = (MAT + 0*N + ((N/2)));
    int *C21 = (MAT + N*((N/2)) + 0);
    int *C22 =  (MAT + N*((N/2))+((N/2)));
    
    
    
    N = N/2;
    
    struct threadData adds[6];
    struct threadData subs[4];
    
    adds[0] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A12, .outMAT = S2, .N = N, .ThreadNum = 0 };
    adds[1] =(struct threadData){ .subMAT1 = A21, .subMAT2 = A22, .outMAT = S3, .N = N, .ThreadNum = 1};
    adds[2] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A22, .outMAT = S5, .N = N, .ThreadNum = 2 };
    adds[3] =(struct threadData){ .subMAT1 = B11, .subMAT2 = B22, .outMAT = S6, .N = N,.ThreadNum = 3 };
    adds[4] =(struct threadData){ .subMAT1 = B21, .subMAT2 = B22, .outMAT = S8, .N = N, .ThreadNum = 4 };
    adds[5] =(struct threadData){ .subMAT1 = B11, .subMAT2 = B12, .outMAT = S10, .N = N, .ThreadNum = 5 };
    
    subs[0] =(struct threadData){ .subMAT1 = B12, .subMAT2 = B22, .outMAT = S1, .N = N, .ThreadNum = 6 };
    subs[1] =(struct threadData){ .subMAT1 = B21, .subMAT2 = B11, .outMAT = S4, .N = N, .ThreadNum = 7};
    subs[2] =(struct threadData){ .subMAT1 = A12, .subMAT2 = A22, .outMAT = S7, .N = N, .ThreadNum = 8 };
    subs[3] =(struct threadData){ .subMAT1 = A11, .subMAT2 = A21, .outMAT = S9, .N = N,.ThreadNum = 9};
    
    
    
    
    pthread_t adders[6];
    pthread_t subtractors[4];
    
    for (int i = 0; i < 6; i++){
        
      		pthread_create(&adders[i], NULL, addition_ThreadRunner, &adds[i]);
        
        
    }
    for (int j = 0; j < 4; j++){
        
      		pthread_create(&subtractors[j], NULL, subtraction_ThreadRunner, &subs[j]);
        
        
    }
    
    
    
    
    for (int k = 0; k < 6; k++)
    {
        pthread_join(adders[k], NULL);
        
        
    }
    
    for (int l = 0; l < 4; l++)
    {
      		pthread_join(subtractors[l], NULL);
        
        
    }
    
    
    if(N > 1){
        int *P1 = (int *)calloc(N * N, sizeof(int));
        int *P2 = (int *)calloc(N * N, sizeof(int));
        int *P3 = (int *)calloc(N * N, sizeof(int));
        int *P4 = (int *)calloc(N * N, sizeof(int));
        int *P5 = (int *)calloc(N * N, sizeof(int));
        int *P6 = (int *)calloc(N * N, sizeof(int));
        int *P7 = (int *)calloc(N * N, sizeof(int));
        
        int *A11A = (int *)malloc(N * N * sizeof(int));
        int *B22B = (int *)malloc(N * N * sizeof(int));
        int *A22A = (int *)malloc(N * N * sizeof(int));
        int *B11B = (int *)malloc(N * N * sizeof(int));
        
        reshape(A11, A11A, N);
        reshape(B11, B11B, N);
        reshape(B22, B22B, N);
        reshape(A22, A22A, N);
        
        if(N  <= cutoff){
            
            naive_IKJ_Square(A11A, S1, P1, N);
            naive_IKJ_Square(S2, B22B, P2, N);
            naive_IKJ_Square(S3, B11B, P3, N);
            naive_IKJ_Square(A22A, S4, P4, N);
            naive_IKJ_Square(S5, S6, P5, N);
            naive_IKJ_Square(S7, S8, P6, N);
            naive_IKJ_Square(S9, S10, P7, N);
            
            
            
        }else{
            
            strassen_recursion_multi(A11A, S1, P1, N, cutoff);
            strassen_recursion_multi(S2, B22B, P2, N, cutoff);
            strassen_recursion_multi(S3, B11B, P3, N, cutoff);
            strassen_recursion_multi(A22A, S4, P4, N, cutoff);
            strassen_recursion_multi(S5, S6, P5, N, cutoff);
            strassen_recursion_multi(S7, S8, P6, N, cutoff);
            strassen_recursion_multi(S9, S10, P7, N, cutoff);
        }
        
        
        struct threadData_C make_C =(struct threadData_C){ .P1 = P1, .P2 = P2, .P3 = P3,.P4 = P4, .P5 = P5, .P6 = P6, .P7 = P7, .N = N, .A11A = A11A, .B22B = B22B, .A22A= A22A, .B11B =B11B, .C11=C11, .C12 = C12, .C21 = C21, .C22 = C22, .ThreadNum = 0 };
        
        pthread_t Make_C_SUBs[4];
        
        
        pthread_create(&Make_C_SUBs[0], NULL, makeC11, &make_C);
        pthread_create(&Make_C_SUBs[1], NULL, makeC12, &make_C);
        pthread_create(&Make_C_SUBs[2], NULL, makeC21, &make_C);
        pthread_create(&Make_C_SUBs[3], NULL, makeC22, &make_C);
        
        pthread_join(Make_C_SUBs[0], NULL);
        pthread_join(Make_C_SUBs[1], NULL);
        pthread_join(Make_C_SUBs[2], NULL);
        pthread_join(Make_C_SUBs[3], NULL);
        
        
        
        
        free(P1);
        free(P2);
        free(P3);
        free(P4);
        free(P5);
        free(P6);
        free(P7);
        
        free(A11A);
        free(B22B);
        free(A22A);
        free(B11B);
        
        
        
        
        
        
    }else{
        
        int P1 = (*A11) * (*S1);
        int P2 = (*S2) * (*B22);
        int P3 = (*S3) * (*B11);
        int P4 = (*A22) * (*S4);
        int P5 = (*S5) * (*S6);
        int P6 = (*S7) * (*S8);
        int P7 = (*S9) * (*S10);
        
        *(MAT) = P5 + P4 - P2 + P6;
        *(MAT + 1) = P1 + P2;
        *(MAT + 2) = P3 + P4;
        *(MAT + 3) = P5 + P1 - P3 - P7;
        
        
    }
    
    free(S1);
    free(S2);
    free(S3);
    free(S4);
    free(S5);
    free(S6);
    free(S7);
    free(S8);
    free(S9);
    free(S10);
}



void strassen_recursion(int *MATA, int *MATB, int*MAT, int N, int cutoff){
    
    
    
    
    int *S1 = (int *)malloc(N * N * sizeof(int));
    int *S2 = (int *)malloc(N * N * sizeof(int));
    int *S3 = (int *)malloc(N * N * sizeof(int));
    int *S4 = (int *)malloc(N * N * sizeof(int));
    int *S5 = (int *)malloc(N * N * sizeof(int));
    int *S6 = (int *)malloc(N * N * sizeof(int));
    int *S7 = (int *)malloc(N * N * sizeof(int));
    int *S8 = (int *)malloc(N * N * sizeof(int));
    int *S9 = (int *)malloc(N * N * sizeof(int));
    int *S10 = (int *)malloc(N * N * sizeof(int));
    
    
    int *B11 = (MATB + 0*N + 0);
    int *B12 = (MATB + 0*N + ((N/2)));
    int *B21 = (MATB + N*((N/2)) + 0);
    int *B22 =  (MATB + N*((N/2))+((N/2)));
    
    
    
    int *A11 = (MATA + 0*N + 0);
    int *A12 = (MATA + 0*N + ((N/2)));
    int *A21 = (MATA + N*((N/2)) + 0);
    int *A22 =  (MATA + N*((N/2))+((N/2)));
    
    int *C11 = (MAT + 0*N + 0);
    int *C12 = (MAT + 0*N + ((N/2)));
    int *C21 = (MAT + N*((N/2)) + 0);
    int *C22 =  (MAT + N*((N/2))+((N/2)));
    
    
    
    N = N/2;
    
    
    subtraction_SUB_MAT(B12, B22, S1, N);
    addition_SUB_MAT(A11, A12, S2, N);
    addition_SUB_MAT(A21, A22, S3, N);
    subtraction_SUB_MAT(B21, B11, S4, N);
    addition_SUB_MAT(A11, A22, S5, N);
    addition_SUB_MAT(B11, B22, S6, N);
    subtraction_SUB_MAT(A12,A22,S7,N);
    addition_SUB_MAT(B21, B22, S8, N);
    subtraction_SUB_MAT(A11, A21,S9, N);
    addition_SUB_MAT(B11, B12, S10, N);
    
    if(N > 1){
        int *P1 = (int *)calloc(N * N, sizeof(int));
        int *P2 = (int *)calloc(N * N, sizeof(int));
        int *P3 = (int *)calloc(N * N, sizeof(int));
        int *P4 = (int *)calloc(N * N, sizeof(int));
        int *P5 = (int *)calloc(N * N, sizeof(int));
        int *P6 = (int *)calloc(N * N, sizeof(int));
        int *P7 = (int *)calloc(N * N, sizeof(int));
        
        int *A11A = (int *)malloc(N * N * sizeof(int));
        int *B22B = (int *)malloc(N * N * sizeof(int));
        int *A22A = (int *)malloc(N * N * sizeof(int));
        int *B11B = (int *)malloc(N * N * sizeof(int));
        
        reshape(A11, A11A, N);
        reshape(B11, B11B, N);
        reshape(B22, B22B, N);
        reshape(A22, A22A, N);
        
        //printMAT(S4, N);
        
        
        
        if(N  <= cutoff){
            //printMAT(P1, N);
            naive_IKJ_Square(A11A, S1, P1, N);
            naive_IKJ_Square(S2, B22B, P2, N);
            naive_IKJ_Square(S3, B11B, P3, N);
            naive_IKJ_Square(A22A, S4, P4, N);
            naive_IKJ_Square(S5, S6, P5, N);
            naive_IKJ_Square(S7, S8, P6, N);
            naive_IKJ_Square(S9, S10, P7, N);
            
            
            
        }else{
            
            strassen_recursion(A11A, S1, P1, N, cutoff);
            strassen_recursion(S2, B22B, P2, N, cutoff);
            strassen_recursion(S3, B11B, P3, N, cutoff);
            strassen_recursion(A22A, S4, P4, N, cutoff);
            strassen_recursion(S5, S6, P5, N, cutoff);
            strassen_recursion(S7, S8, P6, N, cutoff);
            strassen_recursion(S9, S10, P7, N, cutoff);
        }
        
        
        
        //c11
        int *C11_Result1 = (int *)malloc(N * N * sizeof(int));
        int *C11_Result2 = (int *)malloc(N * N * sizeof(int));
        addition(P5, P4, C11_Result1, N);
        subtraction(C11_Result1, P2, C11_Result2, N);
        addition_SUB_MAT_OUT(C11_Result2, P6, C11, N);
        
        addition_SUB_MAT_OUT(P1, P2, C12, N);
        
        addition_SUB_MAT_OUT(P3, P4, C21, N);
        
        int *C22_Result1 = (int *)malloc(N * N * sizeof(int));
        int *C22_Result2 = (int *)malloc(N * N * sizeof(int));
        addition(P5, P1, C22_Result1, N);
        subtraction(C22_Result1, P3, C22_Result2, N);
        subtraction_SUB_MAT_OUT(C22_Result2, P7, C22, N);
        
        
        free(P1);
        free(P2);
        free(P3);
        free(P4);
        free(P5);
        free(P6);
        free(P7);
        
        free(A11A);
        free(B22B);
        free(A22A);
        free(B11B);
        
        
        
        
    }else{
        
        int P1 = (*A11) * (*S1);
        int P2 = (*S2) * (*B22);
        int P3 = (*S3) * (*B11);
        int P4 = (*A22) * (*S4);
        int P5 = (*S5) * (*S6);
        int P6 = (*S7) * (*S8);
        int P7 = (*S9) * (*S10);
        
        *(MAT) = P5 + P4 - P2 + P6;
        *(MAT + 1) = P1 + P2;
        *(MAT + 2) = P3 + P4;
        *(MAT + 3) = P5 + P1 - P3 - P7;
        
        
    }
    
    free(S1);
    free(S2);
    free(S3);
    free(S4);
    free(S5);
    free(S6);
    free(S7);
    free(S8);
    free(S9);
    free(S10);
}


void square_MAT_Mul_recursion_multi(int *MATA, int *MATB, int*MATC, int N, int O_width){
    if(N == 1){
        (*MATC) = (*MATA) * (*MATB);
        
    }else{
        int *B11 = (MATB + 0*N + 0);
        int *B12 = (MATB + 0*N + ((N/2)));
        int *B21 = (MATB + O_width*((N/2)) + 0);
        int *B22 =  (MATB + O_width*((N/2))+((N/2)));
        
        
        
        int *A11 = (MATA + 0*N + 0);
        int *A12 = (MATA + 0*N + ((N/2)));
        int *A21 = (MATA + O_width*((N/2)) + 0);
        int *A22 =  (MATA + O_width*((N/2))+((N/2)));
        
        int *C11 = (MATC + 0*N + 0);
        int *C12 = (MATC + 0*N + ((N/2)));
        int *C21 = (MATC + N*((N/2)) + 0);
        int *C22 = (MATC + N*((N/2))+((N/2)));
        
        
        N = N/2;
        if(N*4 < O_width){
            int *C11C = (int *)malloc(N * N * sizeof(int));
            int *C12C = (int *)malloc(N * N * sizeof(int));
            int *C21C = (int *)malloc(N * N * sizeof(int));
            int *C22C = (int *)malloc(N * N * sizeof(int));
            
            int *C11CC = (int *)malloc(N * N * sizeof(int));
            int *C12CC = (int *)malloc(N * N * sizeof(int));
            int *C21CC = (int *)malloc(N * N * sizeof(int));
            int *C22CC = (int *)malloc(N * N * sizeof(int));
            
            square_MAT_Mul_recursion_multi(A11, B11, C11C, N, O_width);
            square_MAT_Mul_recursion_multi(A12, B21, C11CC, N,  O_width);
            addition_SUB_MAT_OUT(C11C, C11CC, C11, N);
            
            square_MAT_Mul_recursion_multi(A11, B12, C12C, N,  O_width);
            square_MAT_Mul_recursion_multi(A12, B22, C12CC, N, O_width);
            addition_SUB_MAT_OUT(C12C, C12CC, C12, N);
            
            square_MAT_Mul_recursion_multi(A21, B11, C21C, N, O_width);
            square_MAT_Mul_recursion_multi(A22, B21, C21CC, N, O_width);
            addition_SUB_MAT_OUT(C21C, C21CC, C21, N);
            
            square_MAT_Mul_recursion_multi(A21, B12, C22C, N, O_width);
            square_MAT_Mul_recursion_multi(A22, B22, C22CC, N, O_width);
            addition_SUB_MAT_OUT(C22C, C22CC, C22, N);
            
            free(C11C);
            free(C12C);
            free(C21C);
            free(C22C);
            free(C11CC);
            free(C12CC);
            free(C21CC);
            free(C22CC);
            
        }else{
            
            struct threadData_SM mulThreads[4];
            pthread_t muls[4];
            
            mulThreads[0] =(struct threadData_SM){ .A1 = A11, .A2 = A12, .B1 = B11, .B2 = B21, .C = C11, .O_width = O_width, .N = N, .ThreadNum = 0 };
            mulThreads[1] =(struct threadData_SM){ .A1 = A11, .A2 = A12, .B1 = B12, .B2 = B22, .C = C12, .O_width = O_width, .N = N, .ThreadNum = 0 };
            mulThreads[2] =(struct threadData_SM){ .A1 = A21, .A2 = A22, .B1 = B11, .B2 = B21, .C = C21, .O_width = O_width, .N = N, .ThreadNum = 0 };
            mulThreads[3] =(struct threadData_SM){ .A1 = A21, .A2 = A22, .B1 = B12, .B2 = B22, .C = C22, .O_width = O_width, .N = N, .ThreadNum = 0 };
            
            
            for (int j = 0; j < 4; j++){
                
                pthread_create(&muls[j], NULL, mulRunner, &mulThreads[j]);
                
                
            }
            
            
            for (int k = 0; k < 4; k++)
            {
                pthread_join(muls[k], NULL);
                
                
            }
        }
        
    }
}
void square_MAT_Mul_recursion(int *MATA, int *MATB, int*MATC, int N, int O_width, int cutoff){
    if(N == 1){
        (*MATC) = (*MATA) * (*MATB);
    }else{
        int *B11 = (MATB + 0*N + 0);
        int *B12 = (MATB + 0*N + ((N/2)));
        int *B21 = (MATB + O_width*((N/2)) + 0);
        int *B22 =  (MATB + O_width*((N/2))+((N/2)));
        
        
        
        int *A11 = (MATA + 0*N + 0);
        int *A12 = (MATA + 0*N + ((N/2)));
        int *A21 = (MATA + O_width*((N/2)) + 0);
        int *A22 =  (MATA + O_width*((N/2))+((N/2)));
        
        int *C11 = (MATC + 0*N + 0);
        int *C12 = (MATC + 0*N + ((N/2)));
        int *C21 = (MATC + N*((N/2)) + 0);
        int *C22 = (MATC + N*((N/2))+((N/2)));
        
        
        N = N/2;
        
        
        
        
        int *C11C = (int *)malloc(N * N * sizeof(int));
        int *C12C = (int *)malloc(N * N * sizeof(int));
        int *C21C = (int *)malloc(N * N * sizeof(int));
        int *C22C = (int *)malloc(N * N * sizeof(int));
        
        int *C11CC = (int *)malloc(N * N * sizeof(int));
        int *C12CC = (int *)malloc(N * N * sizeof(int));
        int *C21CC = (int *)malloc(N * N * sizeof(int));
        int *C22CC = (int *)malloc(N * N * sizeof(int));
        
        square_MAT_Mul_recursion(A11, B11, C11C, N, O_width, cutoff);
        square_MAT_Mul_recursion(A12, B21, C11CC, N,  O_width, cutoff);
        addition_SUB_MAT_OUT(C11C, C11CC, C11, N);
        
        square_MAT_Mul_recursion(A11, B12, C12C, N,  O_width, cutoff);
        square_MAT_Mul_recursion(A12, B22, C12CC, N, O_width, cutoff);
        addition_SUB_MAT_OUT(C12C, C12CC, C12, N);
        
        square_MAT_Mul_recursion(A21, B11, C21C, N, O_width, cutoff);
        square_MAT_Mul_recursion(A22, B21, C21CC, N, O_width, cutoff);
        addition_SUB_MAT_OUT(C21C, C21CC, C21, N);
        
        square_MAT_Mul_recursion(A21, B12, C22C, N, O_width, cutoff);
        square_MAT_Mul_recursion(A22, B22, C22CC, N, O_width, cutoff);
        addition_SUB_MAT_OUT(C22C, C22CC, C22, N);
        
        free(C11C);
        free(C12C);
        free(C21C);
        free(C22C);
        free(C11CC);
        free(C12CC);
        free(C21CC);
        free(C22CC);
        
        
    }
}





void naive_IKJ_Square(int *MAT1, int *MAT2, int*MAT4, int N){
    
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
                *(MAT4 + i*N + j) += *(MAT1 + i*N + k) * *(MAT2 + k*N + j);
            }
        }
    }
    
    
}



void addition(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N + j) = *(MAT1 + i*N + j) + *(MAT2 + i*N + j);
        }
        
    }
}
void subtraction(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N + j) = *(MAT1 + i*N + j) - *(MAT2 + i*N + j);
            // printf(" %d ",*(MAT2 + i*N + j));
        }
        
    }
}

void addition_SUB_MAT(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N + j) = *(MAT1 + i*N*2 + j) + *(MAT2 + i*N*2 + j);
        }
        
    }
}
void subtraction_SUB_MAT(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N + j) = *(MAT1 + i*N*2 + j) - *(MAT2 + i*N*2 + j);
        }
        
    }
}
void addition_SUB_MAT_OUT(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N*2 + j) = *(MAT1 + i*N + j) + *(MAT2 + i*N + j);
        }
        
    }
}
void subtraction_SUB_MAT_OUT(int *MAT1, int *MAT2, int*MAT, int N){
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(MAT + i*N*2 + j) = *(MAT1 + i*N + j) - *(MAT2 + i*N + j);
        }
        
    }
}



void *addition_ThreadRunner(void *arg){
    struct threadData *data = arg;
    for (int i = 0; i < data->N; i++) {
        for (int j = 0; j < data->N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(data->outMAT + i*data->N + j) = *(data->subMAT1 + i*data->N*2 + j) + *(data->subMAT2 + i*data->N*2 + j);
        }
        
    }
    
    
    // printf("THread = %d \n", data->ThreadNum);
    
    pthread_exit(0);
}

void *subtraction_ThreadRunner(void *arg){
    struct threadData *data = arg;
    
    for (int i = 0; i < data->N; i++) {
        for (int j = 0; j < data->N; j++) {
            // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
            *(data->outMAT + i*data->N + j) = *(data->subMAT1 + i*data->N*2 + j) - *(data->subMAT2 + i*data->N*2 + j);
        }
        
    }
    
    
    // printf("THread = %d \n", data->ThreadNum);
    
    pthread_exit(0);
}

void *makeC11(void *arg){
    struct threadData_C *data = arg;
    int *C11_Result1 = (int *)malloc(data->N * data->N * sizeof(int));
    int *C11_Result2 = (int *)malloc(data->N * data->N * sizeof(int));
    addition(data->P5, data->P4, C11_Result1, data->N);
    subtraction(C11_Result1, data->P2, C11_Result2, data->N);
    addition_SUB_MAT_OUT(C11_Result2, data->P6, data->C11, data->N);
    free(C11_Result1);
    free(C11_Result2);
    
    pthread_exit(0);
}
void *makeC12(void *arg){
    struct threadData_C *data = arg;
    addition_SUB_MAT_OUT(data->P1, data->P2, data->C12, data->N);
    pthread_exit(0);
}
void *makeC21(void *arg){
    struct threadData_C *data = arg;
    addition_SUB_MAT_OUT(data->P3, data->P4, data->C21, data->N);
    
    
    pthread_exit(0);
}
void *makeC22(void *arg){
    struct threadData_C *data = arg;
    int *C22_Result1 = (int *)malloc(data->N * data->N * sizeof(int));
    int *C22_Result2 = (int *)malloc(data->N * data->N * sizeof(int));
    addition(data->P5, data->P1, C22_Result1, data->N);
    subtraction(C22_Result1, data->P3, C22_Result2, data->N);
    subtraction_SUB_MAT_OUT(C22_Result2, data->P7, data->C22, data->N);
    free(C22_Result1);
    free(C22_Result2);
    
    pthread_exit(0);
}

void *mulRunner(void *arg){
    
    
    struct threadData_SM *data = arg;
    
    int *CC1 = (int *)calloc(data->N * data->N, sizeof(int));
    int *CC2 = (int *)calloc(data->N * data->N,  sizeof(int));
    
    square_MAT_Mul_recursion_multi(data->A1, data->B1, CC1, data->N, data->O_width);
    square_MAT_Mul_recursion_multi(data->A2, data->B2, CC2, data->N,  data->O_width);
    addition_SUB_MAT_OUT(CC1, CC2, data->C, data->N);
    
    free(CC1);
    free(CC2);
    
    pthread_exit(0);
    
}
void * naive(void *arg){
    struct threadData_P*data = arg;
    for (int i = 0; i < data->N; i++) {
        for (int k = 0; k < data->N; k++) {
            for (int j = 0; j < data->N; j++) {
                // MAT3[i][j] += MAT1[i][k] * MAT2[k][j];
                *(data->C + i*data->N + j) += *(data->A + i*data->N + k) * *(data->B + k*data->N + j);
            }
        }
    }
    pthread_exit(0);
    
    
}



void reshape(int *MATAA, int *MATBB, int N){
    for (int i = 0; i <  N; i++){
        for (int j = 0; j < N; j++){
            *(MATBB + i*N + j) = *(MATAA + i*N*2 + j);
            
        }
    }
}


void printMAT(int *MAT,int N){
    
    printf("N = %d \n", N);
    for (int i = 0; i <  N; i++){
        printf("\n");
        for (int j = 0; j < N; j++){
            printf(" %d ",*(MAT + i*N + j));
            
        }
    }
}


double get_time() //https://stackoverflow.com/questions/2349776/how-can-i-benchmark-c-code-easily
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}
void proof (int *MAT1, int *MAT2, int*MAT, int  N){
    
    int *MAT4 = (int *)calloc(N * N, sizeof(int));
    int *MAT5 = (int *)calloc(N * N, sizeof(int));
    
    int cutoff;
    if(N < 4000){
        cutoff = 128;
        
    }else{
        cutoff = 512;
    }
    
    
    
    
    strassen_recursion_multi_plus(MAT1, MAT2, MAT4, N, cutoff);
    square_MAT_Mul_recursion_multi(MAT1, MAT2, MAT5, N, N);
    
    
    
    
    
    naive_IKJ_Square(MAT1, MAT2, MAT, N);
    
    
    
    //printMAT(MAT4, N);
    //printMAT(MAT, N);
    int strassen = 0;
    int recursion = 0;
    
    for (int i = 0; i <  N; i++){
        for (int j = 0; j < N; j++){
            if(*(MAT4 + i*N + j) !=  *(MAT + i*N + j)){
                printf("fail");
                strassen = 1;
            }
            if(*(MAT5 + i*N + j) !=  *(MAT + i*N + j)){
                printf("fail");
                recursion = 1;
            }
        }
    }
    if(strassen == 0){
        printf("\nStrassen Algorithm output is acurate");
    }else{
        printf("\nStrassen Algorithm output fails");
    }
    if(recursion == 0){
        printf("\nRecursion Algorithmoutput is acurate");
    }else{
        printf("\nRecursion output fails");
    }
    printf("\nDone Proofing");
    
    
    
    
    
}
void test(){
    
    int N = 1024;
    
    
    
    
    //printMAT(MAT2, N);
    //printMAT(MAT1, N);
    //square_MAT_Mul_recursion_multi(MAT1, MAT2, MAT4, N, N);
    //naive_IKJ_Square(MAT1, MAT2, MAT4, N);
    
    //proof(MAT1, MAT2, MAT4, N);
    
    printf("\n");
    
    double bench = 0;
    double start =  0;
    
    
    int test_sizes[] ={ 32, 64, 128, 256, 512, 1024, 2048, 4096,};
    
    printf("naive_IKJ_Square\n");
    for(int i = 0; i < 8; i++){
        N = test_sizes[i];
        
        int *MAT1 = (int *)malloc(N * N * sizeof(int));
        int *MAT2 = (int *)malloc(N * N * sizeof(int));
        int *MAT4= (int *)calloc(N * N, sizeof(int));
        
        
        int cutoff;
        if(N < 4100){
            cutoff = 128;
            
        }else{
            cutoff = 256;
        }
        srand(time(NULL));
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT1 + i*N + j) = rand() % 31 + 10;
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT2 + i*N + j) = rand() % 45;
        
        
        
        start =get_time();
        naive_IKJ_Square(MAT1, MAT2, MAT4, N);
        bench = (get_time() - start);
        printf("Size = %d X %d\t", test_sizes[i],test_sizes[i]);
        if(i == 0 || i ==1){
            printf("\t");
        }
        printf("Time = %f \n", bench);
        
        free(MAT1);
        free(MAT2);
        free(MAT4);
    }
    
    
    printf("strassen_recursion\n");
    for(int i = 0; i < 8; i++){
        N = test_sizes[i];
        
        int *MAT1 = (int *)malloc(N * N * sizeof(int));
        int *MAT2 = (int *)malloc(N * N * sizeof(int));
        int *MAT4= (int *)calloc(N * N, sizeof(int));
        
        
        int cutoff;
        if(N < 4100){
            cutoff = 128;
            
        }else{
            cutoff = 256;
        }
        srand(time(NULL));
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT1 + i*N + j) = rand() % 31 + 10;
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT2 + i*N + j) = rand() % 45;
        
        
        
        start =get_time();
        strassen_recursion(MAT1, MAT2, MAT4, N, cutoff);
        bench = (get_time() - start);
        printf("Size = %d X %d\t", test_sizes[i],test_sizes[i]);
        if(i == 0 || i ==1){
            printf("\t");
        }
        printf("Time = %f \n", bench);
        
        free(MAT1);
        free(MAT2);
        free(MAT4);
    }
    
    
    
    
    printf("strassen_recursion_multi\n");
    for(int i = 0; i < 8; i++){
        N = test_sizes[i];
        
        int *MAT1 = (int *)malloc(N * N * sizeof(int));
        int *MAT2 = (int *)malloc(N * N * sizeof(int));
        int *MAT4= (int *)calloc(N * N, sizeof(int));
        
        
        int cutoff;
        if(N < 4100){
            cutoff = 128;
            
        }else{
            cutoff = 256;
        }
        srand(time(NULL));
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT1 + i*N + j) = rand() % 31 + 10;
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT2 + i*N + j) = rand() % 45;
        
        
        
        start =get_time();
        strassen_recursion_multi(MAT1, MAT2, MAT4, N, cutoff);
        bench = (get_time() - start);
        printf("Size = %d X %d\t", test_sizes[i],test_sizes[i]);
        if(i == 0 || i ==1){
            printf("\t");
        }
        printf("Time = %f \n", bench);
        
        free(MAT1);
        free(MAT2);
        free(MAT4);
    }
    
    printf("strassen_recursion_multi_plus\n");
    for(int i = 0; i < 8; i++){
        N = test_sizes[i];
        
        int *MAT1 = (int *)malloc(N * N * sizeof(int));
        int *MAT2 = (int *)malloc(N * N * sizeof(int));
        int *MAT4= (int *)calloc(N * N, sizeof(int));
        
        
        int cutoff;
        if(N < 4100){
            cutoff = 128;
            
        }else{
            cutoff = 256;
        }
        srand(time(NULL));
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT1 + i*N + j) = rand() % 31 + 10;
        
        
        for (int i = 0; i <  N; i++)
            for (int j = 0; j < N; j++)
                *(MAT2 + i*N + j) = rand() % 45;
        
        
        
        start =get_time();
        strassen_recursion_multi_plus(MAT1, MAT2, MAT4, N, cutoff);
        bench = (get_time() - start);
        printf("Size = %d X %d\t", test_sizes[i],test_sizes[i]);
        if(i == 0 || i ==1){
            printf("\t");
        }
        printf("Time = %f \n", bench);
        
        free(MAT1);
        free(MAT2);
        free(MAT4);
    }
    
    
}



