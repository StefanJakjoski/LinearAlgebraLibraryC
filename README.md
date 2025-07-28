# LinearAlgebraLibraryC
A library written in C compatible with linear algebra. 
Contains basic functionalities used in Linear Algebra. For the time being, the library is limited to vector support.

LIST OF FUNCTIONS:
// CreateVector FUNCTION
// DEFINE VECTOR ARRAY WITH POSITIVE DIMENSION AT SPECIFIED POINTER
//
// Vector* newVector;
// if(!CreateVector(&newVector, (int) vectorDimension))
//   errorHandler(); 
//

// CreateVector2 FUNCTION
// RETURNS POINTER TO VECTOR ARRAY WITH SPECIFIED DIMENSION
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* newVector = CreateVector2((int) vectorDimension);
//

// CompareVectors FUNCTION
// COMPARES DIMENSIONS AND INDIVIDUAL ELEMENTS OF 2 SPECIFIED VECTORS
// RETURNS 0 IF VECTORS ARE NOT IDENTICAL, RETURNS 1 OTHERWISE
//
// Vector* newVector = CreateVector2((int) vectorDimension);
//

// VectorToString FUNCTION
// PRINT CONTENTS OF VECTOR a TO STANDARD OUTPUT
//
// if(!VectorToString((Vector *) a))
//   errorHandler(); 
//

// ArrayToVector FUNCTION
// DEFINE VECTOR ARRAY WITH POSITIVE DIMENSION AT SPECIFIED POINTER
// FROM ARRAY OF TYPE DOUBLE
//
// Vector* newVector;
// if(!ArrayToVector(&newVector, (double *) array, (int) dimensionSize))
//   errorHandler(); 
//

// DuplicateVector FUNCTION
// DUPLICATE VECTOR ARRAY TO SPECIFIED POINTER
//
// Vector* newVector;
// if(!DuplicateVector(&newVector, (Vector *) original))
//   errorHandler(); 
//

// DuplicateVector2 FUNCTION
// RETURN DUPLICATE OF SPECIFIED VECTOR ARRAY
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* newVector = DuplicateVector2((Vector *) original);
//

// SumVectors FUNCTION 
// MATHEMATICAL SUMMATION OF VECTORS a AND b WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR a
//
// if(!SumVectors((Vector *) a, (Vector *) b))
//   errorHandler(); 
//

// SumVectors2 FUNCTION 
// RETURNS MATHEMATICAL SUMMATION OF VECTORS a AND b AS NEW VECTOR
// RETURNS NULL IF DIMENSIONS <= 0
//
// Vector* c = SumVectors2((Vector *) a, (Vector *) b); 
//

// SubtractVectors FUNCTION 
// MATHEMATICAL SUBTRACTION OF VECTORS a AND b WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR a
//
// if(!SubtractVectors((Vector *) a, (Vector *) b))
//   errorHandler(); 
//

// SubtractVectors2 FUNCTION 
// RETURNS MATHEMATICAL SUBTRACTION OF VECTORS a AND b AS NEW VECTOR
// RETURNS NULL IF DIMENSIONS <= 0
//
// Vector* c = SubtractVectors2((Vector *) a, (Vector *) b); 
//

// MultiplyVector FUNCTION 
// MATHEMATICAL MULTIPLICATION OF VECTOR v AND SCALAR COEFFICIENT c WRITTEN OVER VECTOR a
// DOES NOT PRESERVE ORIGINAL VECTOR v
//
// if(!MultiplyVector((Vector *) v, (int) b))
//   errorHandler(); 
//

// MultiplyVector2 FUNCTION 
// OUTPUTS MATHEMATICAL MULTIPLICATION OF VECTOR v AND SCALAR COEFFICIENT c
// PRESERVES ORIGINAL VECTOR v
// 
// Vector* product = MultiplyVector2((Vector *) v, (double) c);
//

// CrossProduct FUNCTION
// CALCULATE CROSS PRODUCT OF 2 VECTORS AND SAVE TO POINTER
//
// Vector* crossProductVector;
// if(!CrossProduct(&crossProductVector, (Vector *) a, (Vector *) b))
//   errorHandler(); 
//

// DotProduct FUNCTION
// CALCULATE DOT PRODUCT OF 2 VECTORS TO SPECIFIED POINTER.
//
// double dotProduct;
// if(!DotProduct(&dotProduct, (Vector *) a, (Vector *) b))
//   errorHandler(); 
//

// DotProduct2 FUNCTION
// RETURN DOT PRODUCT OF 2 VECTORS AS DOUBLE
//
// double dotProduct = DotProduct2((Vector *) a, (Vector *) b);
// if(isnan(dotProduct))
//   errorHandler();
//

// VectorLength FUNCTION
// CALCULATE LENGTH (OR NORM) OF SPECIFIED VECTOR AND RETURN VALUE.
// RETURNS -1 IF VECTOR DIMENSION <= 0
//
// double norm;
// if(norm = VectorLength((Vector *) v) < 0)
//   errorHandler(); 
//

// UnitVector FUNCTION
// CALCULATE UNIT VECTOR FROM SPECIFIED VECTOR AND OVERWRITE ORIGINAL
//
// Vector* v;
// if(!UnitVector((Vector *) v))
//   errorHandler(); 
//

// UnitVector2 FUNCTION
// CALCULATE UNIT VECTOR FROM SPECIFIED VECTOR RETURN VECTOR OUTPUT
// RETURNS NULL IF DIMENSION <= 0
//
// Vector* unitVector = unitVector2((Vector *) v);
//