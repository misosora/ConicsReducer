#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Sets values of matrix according to coeficients of original equation
void setMatrix(double **matrix, double a, double b, double c, double d, double e, double f){
    matrix[0][0] = a;
    matrix[1][1] = c;
    matrix[2][2] = f;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            if(i != j){
                if(i+j == 1){
                    matrix[i][j] = b/2.0;
                }
                else if(i+j == 2){
                    matrix[i][j] = d/2.0;
                }
                else{
                    matrix[i][j] = e/2.0;
                }
            }
        }
    }
}

// Calculates determinant of the first four values of matrix
double calcDet(double **matrix){
    double positive = matrix[0][0] * matrix[1][1];
    double negative = - (matrix[0][1] * matrix[1][0]);

    double det = positive + negative;

    return det;
}

// Eliminates linear terms (x and/or y)
void translation(double **matrix, double *a, double *b, double *c, double *d, double *e, double *f){
    double det = calcDet(matrix);
    if(det == 0){
        printf("Translation is impossible.\n\n");
        return;
    }

    double detEq = (matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]);
    double h = ((matrix[0][1] * matrix[1][2]) - (matrix[1][1] * matrix[0][2])) / detEq;
    double k = ((matrix[0][2] * matrix[1][0]) - (matrix[1][2] * matrix[0][0])) / detEq;

    (*f) = (matrix[2][0] * h) + (matrix[2][1] * k) + matrix[2][2];
    (*d) = 0;
    (*e) = 0;
}

// Eliminates mixed term (xy)
void rotation(double *a, double *b, double *c, double *d, double *e, double *f){
    double cotg2 = ((*a) - (*c)) / (*b);
    double inSin2 = sqrt(1+pow(cotg2, 2));

    double firstEq = (*a) + (*c);
    double secEq = (*b) * inSin2;

    (*a) = (firstEq + secEq)/2;
    (*c) = firstEq - (*a);
    (*b) = 0;

    if((*d) != 0 || (*e) != 0){
        double sin2 = 1/inSin2;
        double cos2 = cotg2*sin2;

        double cos = sqrt((1+cos2)/2);
        double sin = sqrt(1-pow(cos, 2));

        double aux = (*d);
        (*d) = ((*d) * cos) + ((*e) * sin);
        (*e) = (aux * (-sin)) + ((*e) * cos);
    }
}

void operations(double **matrix, double *a, double *b, double *c, double *d, double *e, double *f){
    if((*d) != 0 || (*e) != 0){
        translation(matrix, a, b, c, d, e, f);
    }

    if((*b) != 0){
        rotation(a, b, c, d, e, f);
    }
}

void freeMatrix(double **matrix){
    for(int i = 0; i < 3; i++){
        free(matrix[i]);
    }
    free(matrix);
}

int main(int argc, char *argv[]){
    double a, b, c, d, e, f;
    scanf("%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f);
    
    double **matrix = (double **) malloc(3*sizeof(double *));
    for(int i = 0; i < 3; i++){
        matrix[i] = (double *) malloc(3*sizeof(double));
    }
    setMatrix(matrix, a, b, c, d, e, f);

    operations(matrix, &a, &b, &c, &d, &e, &f);
    
    freeMatrix(matrix);

    return 0;
}
