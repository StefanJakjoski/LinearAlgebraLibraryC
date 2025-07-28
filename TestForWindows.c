#include "LinearAlgebra.h"
#include <stdbool.h>
#include <windows.h>

#define ANSI_GREEN "\e[0;32m"
#define ANSI_RED "\e[0;31m"
#define ANSI_DEF "\x1b[0m"

void Eval(bool res){
    if(res){ 
        printf("\x1b[32m Clear\x1b[0m\n");
    }else{
        printf("\x1b[31m Error\x1b[0m\n");
    }
}

bool CreateVectorTest(){
    Vector* v;
    CreateVector(&v, 2);

    if(v->vectorDimension != 2)
        return false;
    
    return !v->vectorArray[0] && !v->vectorArray[1];
}

bool CompareVectorsTest(){
    Vector* v;
    CreateVector(&v, 2);
    
    Vector* d;
    CreateVector(&d, 2);
    
    return CompareVectors(v, d);
}

bool ArrayToVectorTest(){
    double* arr = (double *) calloc(2, sizeof(double));
    arr[0] = 1.0; arr[1] = 2.0;
    Vector* v;
    ArrayToVector(&v, arr, 2);

    Vector* d;
    CreateVector(&d, 2);
    d->vectorArray[0] = 1; d->vectorArray[1] = 2;

    return CompareVectors(v, d);
}

bool DuplicateVectorTest(){
    Vector* d;
    CreateVector(&d, 2);
    d->vectorArray[0] = 1; d->vectorArray[1] = 2;

    Vector* v;
    DuplicateVector(&v, d);
    
    return CompareVectors(v, d);
}

bool SumVectorsTest(){
    Vector* d;
    CreateVector(&d, 2);
    d->vectorArray[0] = 1; d->vectorArray[1] = 2;
    
    Vector* v;
    CreateVector(&v, 2);
    v->vectorArray[0] = 3; v->vectorArray[1] = 3;
    
    Vector* vd;
    CreateVector(&vd, 2);
    vd->vectorArray[0] = 4; vd->vectorArray[1] = 5;

    return CompareVectors(vd, SumVectors2(v, d));
}

bool SubtractVectorsTest(){
    Vector* d;
    CreateVector(&d, 2);
    d->vectorArray[0] = 1; d->vectorArray[1] = 2;
    
    Vector* v;
    CreateVector(&v, 2);
    v->vectorArray[0] = 3; v->vectorArray[1] = 3;
    
    Vector* vd;
    CreateVector(&vd, 2);
    vd->vectorArray[0] = 2; vd->vectorArray[1] = 1;

    return CompareVectors(vd, SubtractVectors2(v, d));
}

bool MultiplyVectorTest(){
    Vector* d;
    CreateVector(&d, 2);
    d->vectorArray[0] = 1; d->vectorArray[1] = 2;
    
    Vector* vc;
    CreateVector(&vc, 2);
    vc->vectorArray[0] = 2; vc->vectorArray[1] = 4;

    return CompareVectors(vc, MultiplyVector2(d, 2));
}

bool CrossProductTest(){
    Vector* a;
    double* aA = (double *) calloc(3, sizeof(double));
    aA[0] = 1; aA[1] = 3; aA[2] = 2;
    ArrayToVector(&a, aA, 3);
    Vector* b;
    double* bA = (double *) calloc(3, sizeof(double));
    bA[0] = 1.5; bA[1] = 4; bA[2] = 13;
    ArrayToVector(&b, bA, 3);
    Vector* c;
    CrossProduct(&c, a, b);

    Vector* comp;
    CreateVector(&comp, 3);
    comp->vectorArray[0] = 31; comp->vectorArray[1] = 10; comp->vectorArray[2] = -0.5;

    return CompareVectors(comp, c);
}

bool VectorLengthTest(){
    Vector* comp;
    CreateVector(&comp, 3);
    comp->vectorArray[0] = 31; comp->vectorArray[1] = 10; comp->vectorArray[2] = -0.5;

    return VectorLength(comp) == sqrt(31*31 + 100 + 0.25);
}

bool UnitVectorTest(){
    Vector* comp;
    CreateVector(&comp, 3);
    comp->vectorArray[0] = 31; comp->vectorArray[1] = 10; comp->vectorArray[2] = -0.5;

    return VectorLength(UnitVector2(comp)) == 1.0;
}

bool DotProductTest(){
    Vector* a;
    double* aA = (double *) calloc(3, sizeof(double));
    aA[0] = 1; aA[1] = 3; aA[2] = 2;
    ArrayToVector(&a, aA, 3);
    Vector* b;
    double* bA = (double *) calloc(3, sizeof(double));
    bA[0] = 1.5; bA[1] = 4; bA[2] = 13;
    ArrayToVector(&b, bA, 3);

    double dot;
    if(!DotProduct(&dot, a, b))
        return false;

    return dot == 39.5;
}


int main(){
    //Enable colored text
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD dwMode = 0;

    // Get the current console mode
    GetConsoleMode(hConsole, &dwMode);

    // Enable virtual terminal processing
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    SetConsoleMode(hConsole, dwMode);
    //END Enable colored text

    printf("Testing...\n");
    printf("CreateVector() "); Eval(CreateVectorTest());
    printf("CompareVectors() "); Eval(CompareVectorsTest());
    printf("ArrayToVector() "); Eval(ArrayToVectorTest());
    printf("DuplicateVector() "); Eval(DuplicateVectorTest());
    printf("SumVectors() "); Eval(SumVectorsTest());
    printf("SubtractVectors() "); Eval(SubtractVectorsTest());
    printf("MultiplyVector() "); Eval(MultiplyVectorTest());
    printf("CrossProduct() "); Eval(CrossProductTest());
    printf("VectorLength() "); Eval(VectorLengthTest());
    printf("UnitVector() "); Eval(UnitVectorTest());
    printf("DotProduct() "); Eval(DotProductTest());
}