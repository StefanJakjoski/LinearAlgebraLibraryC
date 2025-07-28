#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double* vectorArray;
    int vectorDimension;
} Vector;


// CreateVector FUNCTION
// DEFINE VECTOR ARRAY WITH POSITIVE DIMENSION AT SPECIFIED POINTER
//
// Vector* newVector;
// if(!CreateVector(&newVector, (int) vectorDimension))
//   errorHandler(); 
//
int CreateVector(Vector** res, int vectorSize){
    if(vectorSize <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return 0;
    }

    Vector* buffer = (Vector *) calloc(1, sizeof(Vector));
    buffer->vectorDimension = vectorSize;
    buffer->vectorArray = (double *) calloc(vectorSize, sizeof(double));
    *res = buffer;
    return 1;
}

// CreateVector2 FUNCTION
// RETURNS POINTER TO VECTOR ARRAY WITH SPECIFIED DIMENSION
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* newVector = CreateVector2((int) vectorDimension);
//
Vector* CreateVector2(int vectorSize){
    Vector* v;
    if(!CreateVector(&v, vectorSize))
        return NULL;

    return v;
}

// CompareVectors FUNCTION
// COMPARES DIMENSIONS AND INDIVIDUAL ELEMENTS OF 2 SPECIFIED VECTORS
// RETURNS 0 IF VECTORS ARE NOT IDENTICAL, RETURNS 1 OTHERWISE
//
// Vector* newVector = CreateVector2((int) vectorDimension);
//
int CompareVectors(Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension){
        fprintf(stderr, "Improper vector dimensions.\n");
        return 0;
    }

    for(int i = 0; i < a->vectorDimension; i++){
        if(a->vectorArray[i] != b->vectorArray[i]){
            return 0;
        }
    }

    return 1;
}

// VectorToString FUNCTION
// PRINT CONTENTS OF VECTOR a TO STANDARD OUTPUT
//
// if(!VectorToString((Vector *) a))
//   errorHandler(); 
//
int VectorToString(Vector* a){
    if(a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return 0;
    }

    printf("{");
    for(int i = 0; i < a->vectorDimension; i++){
        if(i >= 1){ printf(", "); }
        printf("%f", a->vectorArray[i]);
    }
    printf("}\n");

    return 1;
}

// ArrayToVector FUNCTION
// DEFINE VECTOR ARRAY WITH POSITIVE DIMENSION AT SPECIFIED POINTER
// FROM ARRAY OF TYPE DOUBLE
//
// Vector* newVector;
// if(!ArrayToVector(&newVector, (double *) array, (int) dimensionSize))
//   errorHandler(); 
//
int ArrayToVector(Vector** res, double* array, int arrayDimension){
    if(!CreateVector(res, arrayDimension)){
        fprintf(stderr, "Error during new vector creation.\n");
        return 0;
    }

    Vector* buffer = *res;
    memcpy(buffer->vectorArray, array, arrayDimension*sizeof(double));
    buffer->vectorDimension = arrayDimension;

    return 1;
}

// DuplicateVector FUNCTION
// DUPLICATE VECTOR ARRAY TO SPECIFIED POINTER
//
// Vector* newVector;
// if(!DuplicateVector(&newVector, (Vector *) original))
//   errorHandler(); 
//
int DuplicateVector(Vector** res, Vector* original){
    if(original->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return 0;
    }

    Vector* buffer = CreateVector2(original->vectorDimension);
    memcpy(buffer->vectorArray, original->vectorArray, original->vectorDimension * sizeof(double));
    *res = buffer;

    return 1;
}

// DuplicateVector2 FUNCTION
// RETURN DUPLICATE OF SPECIFIED VECTOR ARRAY
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* newVector = DuplicateVector2((Vector *) original);
//
Vector* DuplicateVector2(Vector* original){
    Vector* v;
    if(!DuplicateVector(&v, original))
        return NULL;
    
    return v;
}

// SumVectors FUNCTION 
// MATHEMATICAL SUMMATION OF VECTORS a AND b WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR a
//
// (Effectively: a += b)
// if(!SumVectors((Vector *) a, (Vector *) b))
//   errorHandler(); 
//
int SumVectors(Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension || a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimensions.\n");
        return 0;
    }

    for(int i = 0; i < a->vectorDimension; i++){
        a->vectorArray[i] += b->vectorArray[i];
    }

    return 1;
}

// SumVectors2 FUNCTION 
// RETURNS MATHEMATICAL SUMMATION OF VECTORS a AND b AS NEW VECTOR
// RETURNS NULL IF DIMENSIONS <= 0
//
// Vector* c = SumVectors2((Vector *) a, (Vector *) b); 
//
Vector* SumVectors2(Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension || a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimensions.\n");
        return NULL;
    }

    Vector* v = DuplicateVector2(a);
    if(!SumVectors(v, b))
        return NULL;
    
    return v;
}

// SubtractVectors FUNCTION 
// MATHEMATICAL SUBTRACTION OF VECTORS a AND b WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR a
//
// if(!SubtractVectors((Vector *) a, (Vector *) b))
//   errorHandler(); 
//
int SubtractVectors(Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension || a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimensions.\n");
        return 0;
    }

    for(int i = 0; i < a->vectorDimension; i++){
        a->vectorArray[i] -= b->vectorArray[i];
    }

    return 1;
}

// SubtractVectors2 FUNCTION 
// RETURNS MATHEMATICAL SUBTRACTION OF VECTORS a AND b AS NEW VECTOR
// RETURNS NULL IF DIMENSIONS <= 0
//
// Vector* c = SubtractVectors2((Vector *) a, (Vector *) b); 
//
Vector* SubtractVectors2(Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension || a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimensions.\n");
        return NULL;
    }

    Vector* v = DuplicateVector2(a);
    if(!SubtractVectors(v, b))
        return NULL;
    
    return v;
}

// MultiplyVector FUNCTION 
// MATHEMATICAL MULTIPLICATION OF VECTOR v AND SCALAR COEFFICIENT c WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR v
//
// if(!MultiplyVector((Vector *) v, (int) b))
//   errorHandler(); 
//
int MultiplyVector(Vector* a, const double coefficient){
    if(a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return 0;
    }

    for(int i = 0; i < a->vectorDimension; i++){
        a->vectorArray[i] *= coefficient;
    }

    return 1;
}

// MultiplyVector2 FUNCTION 
// OUTPUTS MATHEMATICAL MULTIPLICATION OF VECTOR v AND SCALAR COEFFICIENT c
// PRESERVES ORIGINAL VECTOR v
// 
// Vector* product = MultiplyVector2((Vector *) v, (double) c);
//
Vector* MultiplyVector2(Vector* a, const double coefficient){
    if(a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return NULL;
    }

    Vector* v = DuplicateVector2(a);
    if(!MultiplyVector(v, coefficient))
        return NULL;
    
    return v;
}

// CrossProduct FUNCTION
// CALCULATE CROSS PRODUCT OF 2 VECTORS AND SAVE TO POINTER
//
// Vector* crossProductVector;
// if(!CrossProduct(&crossProductVector, (Vector *) a, (Vector *) b))
//   errorHandler(); 
//
int CrossProduct(Vector** newVector, Vector* a, Vector* b){
    int properVectorDimension = 3;
    if(a->vectorDimension != properVectorDimension || b->vectorDimension != properVectorDimension){
        fprintf(stderr, "Improper vector dimensions.\n");
        return 0;
    }

    if(!CreateVector(newVector, properVectorDimension)){
        fprintf(stderr, "Error during new vector creation.\n");
        return 0;
    }

    Vector *buffer = *newVector;
    buffer->vectorArray[0] = a->vectorArray[1]*b->vectorArray[2] - a->vectorArray[2]*b->vectorArray[1];
    buffer->vectorArray[1] = a->vectorArray[0]*b->vectorArray[2] - a->vectorArray[2]*b->vectorArray[0];
    buffer->vectorArray[2] = a->vectorArray[0]*b->vectorArray[1] - a->vectorArray[1]*b->vectorArray[0];

    return 1;
}

// DotProduct FUNCTION
// CALCULATE DOT PRODUCT OF 2 VECTORS TO SPECIFIED POINTER.
//
// double dotProduct;
// if(!DotProduct(&dotProduct, (Vector *) a, (Vector *) b))
//   errorHandler(); 
//
int DotProduct(double* res, Vector* a, Vector* b){
    if(a->vectorDimension != b->vectorDimension || a->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimensions.\n");
        return 0;
    }

    double buffer = 0;
    for(int i = 0; i < a->vectorDimension; i++){
        buffer += a->vectorArray[i]*b->vectorArray[i];
    }
    *res = buffer;

    return 1;
}

// DotProduct2 FUNCTION
// RETURN DOT PRODUCT OF 2 VECTORS AS DOUBLE
//
// double dotProduct = DotProduct2((Vector *) a, (Vector *) b);
// if(isnan(dotProduct))
//   errorHandler();
//
double DotProduct2(Vector* a, Vector* b){
    double dotProduct;
    if(!DotProduct(&dotProduct, a, b))
        return NAN;
    
    return dotProduct;
}

// VectorLength FUNCTION
// CALCULATE LENGTH (OR NORM) OF SPECIFIED VECTOR AND RETURN VALUE.
// RETURNS -1 IF VECTOR DIMENSION <= 0
//
// double norm;
// if(norm = VectorLength((Vector *) v) < 0)
//   errorHandler(); 
//
double VectorLength(Vector* v){
    if(v->vectorDimension <= 0){
        fprintf(stderr, "Improper vector dimension.\n");
        return -1;
    }

    double res = 0;
    for(int i = 0; i < v->vectorDimension; i++){
        res += pow(v->vectorArray[i], 2);
    }

    return sqrt(res);
}

// UnitVector FUNCTION
// CALCULATE UNIT VECTOR FROM SPECIFIED VECTOR AND OVERWRITE ORIGINAL
//
// Vector* v;
// if(!UnitVector((Vector *) v))
//   errorHandler(); 
//
int UnitVector(Vector* original){
    double norm = VectorLength(original);
    if(norm < 0)
        return 0;
    
    if(!MultiplyVector(original, 1.0/norm))
        return 0;

    return 1;
}

// UnitVector2 FUNCTION
// CALCULATE UNIT VECTOR FROM SPECIFIED VECTOR RETURN VECTOR OUTPUT
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* unitVector = unitVector2((Vector *) v);
//
Vector* UnitVector2(Vector* original){
    Vector* v;
    if(!DuplicateVector(&v, original))
        return NULL;

    if(!UnitVector(v))
        return NULL;

    return v;
}

#endif // LINEAR_ALGEBRA_H