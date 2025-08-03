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
   
    
    bool check = !vGet(v, 0) && !vGet(v, 1);
    FreeVector(v);

    return check;
}

bool CompareVectorsTest(){
    Vector* v;
    CreateVector(&v, 2);
    
    Vector* d;
    CreateVector(&d, 2);
    
    bool check = CompareVectors(v, d);
    FreeVector(v); FreeVector(d);

    return check;
}

bool ArrayToVectorTest(){
    double* arr = (double *) calloc(2, sizeof(double));
    arr[0] = 1.0; arr[1] = 2.0;
    Vector* v;
    ArrayToVector(&v, arr, 2);

    Vector* d;
    CreateVector(&d, 2);
    vSet(d, 0, 1.0); vSet(d, 1, 2.0);

    bool check = CompareVectors(v, d);
    FreeVector(v); FreeVector(d);

    return check;
}

bool DuplicateVectorTest(){
    Vector* d;
    CreateVector(&d, 2);
    vSet(d, 0, 1.0); vSet(d, 1, 2.0);


    Vector* v;
    DuplicateVector(&v, d);
    
    bool check = CompareVectors(v, d);
    FreeVector(v); FreeVector(d);

    return check;
}

bool SumVectorsTest(){
    Vector* d;
    CreateVector(&d, 2);
    vSet(d, 0, 1.0); vSet(d, 1, 2.0);
    
    Vector* v;
    CreateVector(&v, 2);
    vSet(v, 0, 3.0); vSet(v, 1, 3.0);
    
    Vector* vd;
    CreateVector(&vd, 2);
    vSet(vd, 0, 4.0); vSet(vd, 1, 5.0);

    bool check = CompareVectors(vd, SumVectors2(v, d));
    FreeVector(d); FreeVector(v); FreeVector(vd);

    return check;
}

bool SubtractVectorsTest(){
    Vector* d;
    CreateVector(&d, 2);
    vSet(d, 0, 1.0); vSet(d, 1, 2.0);
    
    Vector* v;
    CreateVector(&v, 2);
    vSet(v, 0, 3.0); vSet(v, 1, 3.0);
    
    Vector* vd;
    CreateVector(&vd, 2);
    vSet(vd, 0, 2.0); vSet(vd, 1, 1.0);

    bool check = CompareVectors(vd, SubtractVectors2(v, d));
    FreeVector(d); FreeVector(v); FreeVector(vd);

    return check;
}

bool MultiplyVectorTest(){
    Vector* d;
    CreateVector(&d, 2);
    vSet(d, 0, 1.0); vSet(d, 1, 2.0);
    
    Vector* vc;
    CreateVector(&vc, 2);
    vSet(vc, 0, 2.0); vSet(vc, 1, 4.0);

    bool check = CompareVectors(vc, MultiplyVector2(d, 2));
    FreeVector(d); FreeVector(vc);

    return check;
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
    vSet(comp, 0, 31.0); vSet(comp, 1, 10.0); vSet(comp, 2, -0.5);

    bool check = CompareVectors(comp, c);
    FreeVector(a); FreeVector(b); FreeVector(c); FreeVector(comp);

    return check;
}

bool VectorLengthTest(){
    Vector* comp;
    CreateVector(&comp, 3);
    vSet(comp, 0, 31.0); vSet(comp, 1, 10.0); vSet(comp, 2, -0.5);

    bool check = VectorLength(comp) == sqrt(31*31 + 100 + 0.25);
    FreeVector(comp);

    return check;
}

bool UnitVectorTest(){
    Vector* comp;
    CreateVector(&comp, 3);
    vSet(comp, 0, 31.0); vSet(comp, 1, 10.0); vSet(comp, 2, -0.5);

    bool check = VectorLength(UnitVector2(comp)) == 1.0;
    FreeVector(comp);

    return check;
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
    
    FreeVector(a); FreeVector(b);

    return dot == 39.5;
}

bool ProjectVectorTest(){
    Vector *a = NULL; Vector *b = NULL;
    if(!CreateVector(&a, 2) || !CreateVector(&b, 2))
        return false;

    vSet(a, 0, 1.0); vSet(a, 1, 1.0);
    vSet(b, 0, 2.0); vSet(b, 1, 0.0);

    Vector* projba = NULL;
    if(!ProjectVector(&projba, b, a))
        return 0;
    
    bool check = vGet(projba, 0) && !vGet(projba, 1);
    FreeVector(a); FreeVector(b); FreeVector(projba);

    return check;
}

bool CreateMatrixTest(){
    Matrix* m = NULL;
    CreateMatrix(&m, 2, 2);

    if(m->matrixRows != 2 || m->matrixColumns != 2)
        return false;
    
    Matrix* i = NULL;
    if(!CreateIdentityMatrix(&i, 2))
        return false;
    
    bool basicCheck = !mGet(m,0,0) && !mGet(m,0,1) && !mGet(m,1,0) && !mGet(m,1,1);
    bool identityCheck = mGet(i,0,0) && !mGet(i,0,1) && !mGet(i,1,0) && mGet(i,1,1);
    FreeMatrix(m); FreeMatrix(i);

    return basicCheck && identityCheck;
}

bool CompareMatricesTest(){
    Matrix* a = NULL; Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&a, 2, 1) || !CreateMatrix(&b, 2, 1) || !CreateMatrix(&c, 2, 1))
        return false;

    //a->matrixArray[0][0] = 1.0; a->matrixArray[1][0] = 2.0;
    //b->matrixArray[0][0] = 1.0; b->matrixArray[1][0] = 2.0;
    //c->matrixArray[0][0] = 1.0; c->matrixArray[1][0] = 23.4;

    mSet(a, 0, 0, 1.0); mSet(a, 0, 0, 2.0);
    mSet(b, 0, 0, 1.0); mSet(b, 0, 0, 2.0);
    mSet(c, 0, 0, 1.0); mSet(c, 0, 0, 23.4);

    bool check = CompareMatrices(a, b);
    bool check2 = !CompareMatrices(a, c);
    FreeMatrix(a); FreeMatrix(b); FreeMatrix(c);

    return check && check2;
}

bool TransposeMatrixTest(){
    Matrix* a = NULL; Matrix* b = NULL;
    if(!CreateMatrix(&a, 2, 3) || !CreateMatrix(&b, 3, 2))
        return false;
    
    /*
    a->matrixArray[0][0] = 1; a->matrixArray[0][1] = 2; a->matrixArray[0][2] = 3;
    a->matrixArray[1][0] = 4; a->matrixArray[1][1] = 5; a->matrixArray[1][2] = 6;

    b->matrixArray[0][0] = 1; b->matrixArray[0][1] = 4;
    b->matrixArray[1][0] = 2; b->matrixArray[1][1] = 5;
    b->matrixArray[2][0] = 3; b->matrixArray[2][1] = 6;
    */

    mSet(a, 0, 0, 1); mSet(a, 0, 1, 2.0); mSet(a, 0, 2, 3.0);
    mSet(a, 1, 0, 4.0); mSet(a, 1, 1, 5.0); mSet(a, 1, 2, 6.0);

    
    mSet(b, 0, 0, 1.0); mSet(b, 0, 1, 4.0);
    mSet(b, 1, 0, 2.0); mSet(b, 1, 1, 5.0);
    mSet(b, 2, 0, 3.0); mSet(b, 2, 1, 6.0);

    Matrix* aT = NULL;
    if(!TransposeMatrix(&aT, a))
        return false;

    bool check = CompareMatrices(b, aT);
    FreeMatrix(a); FreeMatrix(b); FreeMatrix(aT);

    return check;
}

bool ArrayToMatrixTest(){
    Matrix* c = NULL;
    if(!CreateMatrix(&c, 2, 2))
        return false;
    
    //c->matrixArray[0][0] = 1; c->matrixArray[0][1] = 23.4;
    //c->matrixArray[1][0] = -321; c->matrixArray[1][1] = 323.4;
    mSet(c, 0, 0, 1.0); mSet(c, 0, 1, 23.4);
    mSet(c, 1, 0, -321.0); mSet(c, 1, 1, 323.4);


    double** array = (double **) calloc(2, sizeof(double *));
    array[0] = (double *) calloc(2, sizeof(double));
    array[1] = (double *) calloc(2, sizeof(double));
    array[0][0] = 1.0; array[0][1] = 23.4;
    array[1][0] = -321; array[1][1] = 323.4;
    Matrix* comp = NULL;
    ArrayToMatrix(&comp, array, 2, 2);

    bool check = CompareMatrices(comp, c);
    FreeMatrix(c); FreeMatrix(comp);

    return check;
}

bool DuplicateMatrixTest(){
    Matrix* c = NULL;
    if(!CreateMatrix(&c, 2, 2))
        return false;
    
    //c->matrixArray[0][0] = 1; c->matrixArray[0][1] = 23.4;
    //c->matrixArray[1][0] = -321; c->matrixArray[1][1] = 323.4;
    mSet(c, 0, 0, 1.0); mSet(c, 0, 1, 23.4);
    mSet(c, 1, 0, -321.0); mSet(c, 1, 1, 323.4);

    Matrix* cD = DuplicateMatrix(c);
    if(cD == NULL)
        return false;

    bool check = CompareMatrices(c, cD);
    FreeMatrix(c); FreeMatrix(cD);

    return check;
}

bool SumMatricesTest(){
    Matrix* a = NULL; Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&a, 2, 1) || !CreateMatrix(&b, 2, 1) || !CreateMatrix(&c, 2, 1))
        return false;
    
    //a->matrixArray[0][0] = 1.0; a->matrixArray[1][0] = 2.0;
    //b->matrixArray[0][0] = 1.0; b->matrixArray[1][0] = 2.0;
    //c->matrixArray[0][0] = 2.0; c->matrixArray[1][0] = 4.0;
    mSet(a, 0, 0, 1.0); mSet(a, 1, 0, 2.0);
    mSet(b, 0, 0, 1.0); mSet(b, 1, 0, 2.0);
    mSet(c, 0, 0, 2.0); mSet(c, 1, 0, 4.0);
    
    Matrix* ab = NULL;
    if(!SumMatrices(&ab, a, b))
        return false;

    bool check = CompareMatrices(c, ab);
    FreeMatrix(a); FreeMatrix(b); FreeMatrix(c); FreeMatrix(ab);

    return check;
}

bool MultiplyMatricesTest(){
    Matrix* a = NULL; Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&a, 2, 3) || !CreateMatrix(&b, 3, 2) || !CreateMatrix(&c, 2, 2))
        return false;

    /*
    a->matrixArray[0][0] = 1.0; a->matrixArray[0][1] = 2.0; a->matrixArray[0][2] = 3.0;
    a->matrixArray[1][0] = 4.0; a->matrixArray[1][1] = 5.0; a->matrixArray[1][2] = 6.0;

    b->matrixArray[0][0] = 3.0; b->matrixArray[0][1] = 4.0;
    b->matrixArray[1][0] = 5.0; b->matrixArray[1][1] = 6.0;
    b->matrixArray[2][0] = 1.0; b->matrixArray[2][1] = 2.0;

    c->matrixArray[0][0] = 16.0; c->matrixArray[0][1] = 22.0;
    c->matrixArray[1][0] = 43.0; c->matrixArray[1][1] = 58.0;
    */
    mSet(a, 0, 0, 1); mSet(a, 0, 1, 2.0); mSet(a, 0, 2, 3.0);
    mSet(a, 1, 0, 4.0); mSet(a, 1, 1, 5.0); mSet(a, 1, 2, 6.0);

    mSet(b, 0, 0, 3.0); mSet(b, 0, 1, 4.0);
    mSet(b, 1, 0, 5.0); mSet(b, 1, 1, 6.0);
    mSet(b, 2, 0, 1.0); mSet(b, 2, 1, 2.0);

    mSet(c, 0, 0, 16.0); mSet(c, 0, 1, 22.0);
    mSet(c, 1, 0, 43.0); mSet(c, 1, 1, 58.0);

    Matrix* ab = NULL;
    if(!MultiplyMatrices(&ab, a, b))
        return false;
    bool check = CompareMatrices(c, ab);
    FreeMatrix(a); FreeMatrix(b); FreeMatrix(c); FreeMatrix(ab);

    return check;
}

bool MatrixScalarMultiplicationTest(){
    Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&b, 3, 2) || !CreateMatrix(&c, 3, 2))
        return false;
    /*
    b->matrixArray[0][0] = 3; b->matrixArray[0][1] = 4;
    b->matrixArray[1][0] = 5; b->matrixArray[1][1] = 6;
    b->matrixArray[2][0] = 1; b->matrixArray[2][1] = 2;
    
    c->matrixArray[0][0] = 9; c->matrixArray[0][1] = 12;
    c->matrixArray[1][0] = 15; c->matrixArray[1][1] = 18;
    c->matrixArray[2][0] = 3; c->matrixArray[2][1] = 6;
    */
    mSet(b, 0, 0, 3.0); mSet(b, 0, 1, 4.0);
    mSet(b, 1, 0, 5.0); mSet(b, 1, 1, 6.0);
    mSet(b, 2, 0, 1.0); mSet(b, 2, 1, 2.0);
    
    mSet(c, 0, 0, 9.0); mSet(c, 0, 1, 12.0);
    mSet(c, 1, 0, 15.0); mSet(c, 1, 1, 18.0);
    mSet(c, 2, 0, 3.0); mSet(c, 2, 1, 6.0);

    Matrix* b2 = NULL;
    if(!MatrixScalarMultiplication(&b2, b, 3.0))
        return false;
    
    bool check = CompareMatrices(c, b2);
    FreeMatrix(b); FreeMatrix(c); FreeMatrix(b2);

    return check;
}

bool SwapMatrixRowsTest(){
    Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&b, 3, 2) || !CreateMatrix(&c, 3, 2))
        return false;

    /*
    b->matrixArray[0][0] = 3; b->matrixArray[0][1] = 4;
    b->matrixArray[1][0] = 5; b->matrixArray[1][1] = 6;
    b->matrixArray[2][0] = 1; b->matrixArray[2][1] = 2;
    
    c->matrixArray[0][0] = 3; c->matrixArray[0][1] = 4;
    c->matrixArray[1][0] = 1; c->matrixArray[1][1] = 2;
    c->matrixArray[2][0] = 5; c->matrixArray[2][1] = 6;
    */
    mSet(b, 0, 0, 3.0); mSet(b, 0, 1, 4.0);
    mSet(b, 1, 0, 5.0); mSet(b, 1, 1, 6.0);
    mSet(b, 2, 0, 1.0); mSet(b, 2, 1, 2.0);
    
    mSet(c, 0, 0, 3.0); mSet(c, 0, 1, 4.0);
    mSet(c, 1, 0, 1.0); mSet(c, 1, 1, 2.0);
    mSet(c, 2, 0, 5.0); mSet(c, 2, 1, 6.0);

    if(!SwapMatrixRows(b, 1, 2))
        return false;
    
    bool check = CompareMatrices(c, b);
    FreeMatrix(b); FreeMatrix(c);

    return check;
}

bool RowManipulationTest(){
    Matrix* b = NULL; Matrix* c = NULL;
    if(!CreateMatrix(&b, 2, 2) || !CreateMatrix(&c, 2, 2))
        return false;

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            //b->matrixArray[i][j] = 1;
            //c->matrixArray[i][j] = 2;
            mSet(b, i, j, 1.0); mSet(c, i, j, 2.0);
        }
    }

    AddRow(b, 1, 0);
    MultiplyRow(b, 0, 2);

    bool check = CompareMatrices(c, b);
    FreeMatrix(b); FreeMatrix(c);

    return check;
}

bool DeterminantRecursiveTest(){
    Matrix* m = NULL;
    if(!CreateIdentityMatrix(&m, 4))
        return false;
    
    Matrix* m2 = NULL;
    if(!MatrixScalarMultiplication(&m2, m, 3.0))
        return false;

    bool check = DeterminantRecursive(m) == 1 && DeterminantRecursive(m2) == 81;
    FreeMatrix(m); FreeMatrix(m2);

    return check;
}

bool InvertSquareMatrixTest(){
    Matrix* c = NULL;
    if(!CreateMatrix(&c, 2, 2))
        return false;
    //c->matrixArray[0][0] = 4.0; c->matrixArray[0][1] = 0.0;
    //c->matrixArray[1][0] = -3.0; c->matrixArray[1][1] = 323.4;

    mSet(c, 0, 0, 4.0); mSet(c, 0, 1, 0.0);
    mSet(c, 1, 0, -3.0); mSet(c, 1, 1, 323.4);

    Matrix* inv = NULL;
    if(!InvertSquareMatrix(&inv, c)){
        return false;
    }
    

    Matrix* check = MultiplyMatrices2(inv, c);
    Matrix* id = CreateIdentityMatrix2(2);

    bool check2 = CompareMatrices(check, id);
    FreeMatrix(c); FreeMatrix(inv); FreeMatrix(check); FreeMatrix(id);

    return check2;
}

bool MultiplyMatrixVectorTest(){
    Matrix* c = CreateMatrix2(2, 2);
    //c->matrixArray[0][0] = 4.5; c->matrixArray[0][1] = 0.0;
    //c->matrixArray[1][0] = -3.0; c->matrixArray[1][1] = 323.4;
    
    mSet(c, 0, 0, 4.5); mSet(c, 0, 1, 0.0);
    mSet(c, 1, 0, -3.0); mSet(c, 1, 1, 323.4);

    Vector* v = CreateVector2(2);
    vSet(v, 0, 2.0); vSet(v, 1, 0.0);

    Vector* cv = CreateVector2(2);
    vSet(cv, 0, 9.0); vSet(cv, 1, -6.0);

    bool check = CompareVectors(cv, MultiplyMatrixVector2(c, v));
    FreeMatrix(c); FreeVector(v); FreeVector(cv);

    return check;
}

bool MPInverseTest(){
    Matrix* c = CreateMatrix2(2, 2);
    //c->matrixArray[0][0] = 4.5; c->matrixArray[0][1] = 0;
    //c->matrixArray[1][0] = -3.0; c->matrixArray[1][1] = 323.4;
    
    mSet(c, 0, 0, 4.5); mSet(c, 0, 1, 0.0);
    mSet(c, 1, 0, -3.0); mSet(c, 1, 1, 323.4);

    Matrix* mp, *inv;
    if(!MoorePenroseInverse(&mp, c))
        return false;
    if(!InvertSquareMatrix(&inv, c))
        return false;
    
    bool check = CompareMatrices(mp, inv);
    FreeMatrix(c); FreeMatrix(mp); FreeMatrix(inv);

    return check;
}

bool QRDecompositionTest(){
    Matrix* c = CreateMatrix2(3,3);
    //c->matrixArray[0][0] = 12.0; c->matrixArray[0][1] = -51.0; c->matrixArray[0][2] = 4.0;
    //c->matrixArray[1][0] = 6.0; c->matrixArray[1][1] = 167.0; c->matrixArray[1][2] = -68.0;
    //c->matrixArray[2][0] = -4.0; c->matrixArray[2][1] = 24.0; c->matrixArray[2][2] = -41.0;
    
    mSet(c, 0, 0, 12.0); mSet(c, 0, 1, -51.0); mSet(c, 0, 2, 4.0);
    mSet(c, 1, 0, 6.0); mSet(c, 1, 1, 167.0); mSet(c, 1, 2, -68.0);
    mSet(c, 2, 0, -4.0); mSet(c, 2, 1, 24.0); mSet(c, 2, 2, -41.0);

    Matrix *Q = NULL; Matrix* R = NULL;
    QRDecomposition(&Q, &R, c);

    Matrix* A = NULL;
    if(!MultiplyMatrices(&A, Q, R))
        return false;

    bool check = CompareMatrices(A, c);
    FreeMatrix(c); FreeMatrix(Q); FreeMatrix(R); FreeMatrix(A);

    return check;   
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
    printf("Size of Vector = %llu\n", sizeof(Vector));
    
    printf("Testing Vectors...\n");
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
    printf("ProjectVector() "); Eval(ProjectVectorTest());
    
    printf("\nTesting Matrices...\n");
    printf("CreateMatrix() "); Eval(CreateMatrixTest());
    printf("CompareMatrices() "); Eval(CompareMatricesTest());
    printf("TransposeMatrix() "); Eval(TransposeMatrixTest());
    printf("ArrayToMatrix() "); Eval(ArrayToMatrixTest());
    printf("DuplicateMatrix() "); Eval(DuplicateMatrixTest());
    printf("MultiplyMatrices() "); Eval(MultiplyMatricesTest());
    printf("SumMatrices() "); Eval(SumMatricesTest());
    printf("MatrixScalarMultiplication() "); Eval(MatrixScalarMultiplicationTest());
    printf("SwapMatrixRows() "); Eval(SwapMatrixRowsTest());
    printf("Add/MultiplyRow() "); Eval(RowManipulationTest());
    printf("DeterminantRecursive() "); Eval(DeterminantRecursiveTest());
    printf("InvertSquareMatrix() "); Eval(InvertSquareMatrixTest());
    printf("MultiplyMatrixVector() "); Eval(MultiplyMatrixVectorTest());
    printf("MoorePenroseInverse() "); Eval(MPInverseTest());
    printf("QRDecomposition() "); Eval(QRDecompositionTest());
}