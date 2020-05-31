#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void print(float **matrix, int n, int m);
float determinant(float **matrix, int n);
void cofactor(float **matrix, float **aux, int p, int q, int n);
void adj(float **matrix, float **adj,int n);
int inverse(float **matrix, float **inverse,int n);
void product(float **result, float **a, float **b, int na, int ma, int nb, int mb);
void identity(float **matrix, int n);
void matrix_copy(float **A, float **B, int n, int m);
void clear_matrix(float **M, int n, int m);

void v_m_product(float *result, float *v, float **m, int nv, int nn, int nm);
void m_v_product(float *result, float **matrix, float*vector, int n, int m, int nv);
float dot_product(float *lv, float *rv, int n);
void print_vector(float *v, int n);

void copyColunm(float* vector, float **matrix, int n, int j);
void setColumn(float **matrix, float *vector, int n, int j);

int findRowOfMin(float** matrix, float* b, int n, int j);
float findMax(float *c, int n);
int indexOf(float value,float *c, int n);



int main(){
  /*
    Formato de entrada:
    
    n (Cantidad de coeficientes)
    m (Cantidad de restricciones)
    c1 c2 .. .. .. .. cn (Coeficientes de la funcion objetivo z) (vector)
    a1,1 a1,2........ a1,n (Coeficientes de la matriz de restricciones) (matriz nxn)
    a1,1 a1,2........ a1,n
    .. .. .. .. .. .. .. 
    .. .. .. .. .. .. .. 
    .. .. .. .. .. .. .. 
    .. .. .. .. .. .. .. 
    am,1 am,2 .. .. .. am,n
    r1 r2 .. .. .. .. rm (Restricciones >= , <= , =) (char) 
    b1 b2 .. .. .. .. bm (Coeficientes del vector a la derecha de las restricciones)(vector)
    m / M (Indica si se debe maximizar (M) o minimizar (m)) (char)
  */

  int n,m;
  float *C,**A,**B,**I,*b, **AI;
  char *r;
  char action, aux_char;

  scanf("%d",&n);  
  scanf("%d",&m);

  //matrices
  A = (float**)malloc(m*sizeof(float*)); //matriz de valores de las restricciones
  B = (float**)malloc(m*sizeof(float*)); //matriz de la base
  AI = (float**)malloc(m*sizeof(float*)); //matriz de A e I concatenadas
  I = (float**)malloc(m*sizeof(float*));
  float** RESULT = (float**)malloc(m*sizeof(float*)); //matriz temporal para guardar resultados intermedios
  //vectores
  float *c_b = (float*)malloc(m*sizeof(float)); //vector de coeficientes en la base
  float *x_b = (float*)malloc(m*sizeof(float)); //vector de de variables basicas
  float *c = (float*)malloc(n*sizeof(float)); //vector de coeficientes de la funcion objetivo
  b = (float*)malloc(m*sizeof(float)); //vector de valores de las restricciones
  
  
  r = (char*)malloc(m*sizeof(char)); //vector de tipo de restriccion

  // Inicializacion de la matriz

  for(int i=0;i<m;i++){
    A[i] = (float*)malloc(n*sizeof(float));
    I[i] = (float*)malloc(n*sizeof(float));
    B[i] = (float*)malloc(m*sizeof(float));
    AI[i] = (float*)malloc((m+n)*sizeof(float));
    RESULT[i] = (float*)malloc((m+n) * sizeof(float));
  }

  // Se leen los coeficientes de  la funcion objetivo

  for(int i=0;i<n;i++) scanf("%f ",&c[i]);

  // Se leen los coeficientes de la matriz de restricciones

  for(int i=0; i<m; i++){
    for(int j=0; j<n ;j++){
      scanf("%f ",&A[i][j]);
    }
  }

  // Se leen lost tipos de restricciones (mayor, menor, igual o estricto) 

  for(int i=0; i<m; i++) scanf("%c ", &r[i]);

  // Se leen los coeficientes a la derecha de las restricciones

  for(int i=0; i<m; i++) scanf("%f ", &b[i]);

  // Lee si el problema es de maximizacion o minimizacion

  scanf("%c ",&action);


  float **inv;
  inv = (float**)malloc(m*sizeof(float*));
  for(int i=0;i<m;i++) inv[i] = (float*)malloc(m*sizeof(float));
  float *column = malloc(m*sizeof(float));
  float **AP = (float**)malloc(m*sizeof(float*));
  for(int i = 0; i < m; i++) AP[i] = (float*)malloc(n*sizeof(float));
  
  int iter_number = 1;
  //Base inicial y coeficientes de la base
  for(int i = 0; i < m; i++){
    c_b[i] = 0;
    x_b[i] = b[i];
  }
    // Se crea la matriz de identidad
  identity(B, m);
  inverse(B, inv, n);

    while(1){
    float max_coef_z = findMax(c, n); //buscamos el valor del coeficiente Z de mayor magnitud de crecimiento
    if(abs(max_coef_z) == 0.0) {
      break;
    }
    int coef_z_index = indexOf(max_coef_z, c, n); // buscamos el indice de dicho valor
      printf("iteracion %d: ---------------\n", iter_number++);
      printf("MAX COEF Z: %f\n", max_coef_z);

      product(RESULT, inv, A, m, m, m, n); //multiplicamos B inversa * A 

    int rowmin = findRowOfMin(RESULT, x_b, m, coef_z_index); // ahora en la columna de dicho coeficiente, buscamos la restricción mas cercana (la de menor magnitud b[i]/M[i, coef_z_index])
    if(rowmin == -1){
      printf("Hay un coeficiente que maximiza Z pero no hay ninguna variable que pueda salir (todas son negativas o cero).");
      break;
    }
    copyColunm(column, A, m, coef_z_index); //copiamos la columna
    setColumn(B, column, m, rowmin); //y la sustituimos en la columna de la variable que sale de la base

    printf("max c[i]: %.2f\n", max_coef_z);
    printf("variable que entra: c[i=%d]\n", coef_z_index);
    printf("variable que sale: x_b[j=%d]\n", rowmin);
    printf("A:\n");
    print(A, m, n);
    printf("\n");
    printf("B:\n");
    print(B, m, m);
    printf("\n");

    inverse(B,inv,m); //hallamos la inversa y la almacenamos en inv
    printf("inversa:\n");
    print(inv,m,m);
    printf("\n");

    //sacamos y metemos de la base
    float swap = c_b[rowmin];
    c_b[rowmin] = c[coef_z_index];
    c[coef_z_index] = swap;

    m_v_product(x_b, inv, b, m,m,m);
    print_vector(x_b, m);
    print_vector(c_b, m);
    printf("Z: %.2f\n", dot_product(c_b, x_b, m));      
    }


  //Deallocation phase
  for(int i = 0; i < m; i++){
    free(A[i]);
    free(I[i]);
    free(B[i]);
    free(AI[i]);
    free(RESULT[i]);
    free(inv[i]);
    free(AP[i]);
  }

  free(c); free(A); free(B); free(b); free(I); free(r); free(inv); free(AI); free(RESULT); free(column); free(x_b); free(c_b); free(AP);

}

void clear_matrix(float **M, int n, int m){
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      M[i][j] = 0;
    }
  }
}

void matrix_copy(float **A, float **B, int n, int m){
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < m; ++j)
    {
      A[i][j] = B[i][j];
    }
  }
}

void identity(float ** matrix, int n){
    for(int i = 0 ;i < n; i++){
      for(int j = 0; j < n; j++){
        if(i==j) matrix[i][j]= 1;
        else matrix[i][j] = 0;
      }
    }
}

void v_m_product(float *result, float *v, float **m, int nv, int nn, int nm){
  if(nv != nn){
    printf("Error in v_m_product: row vector v length must be equal to number of rows of a matrix's column");
    exit(1);
  }
  for(int i = 0; i < nn; i++){
    result[i] = 0;
    for(int j = 0; j < nm; j++){
      result[i] += v[j] * m[j][i];
    }
  }
}

void m_v_product(float *result, float **matrix, float*vector, int n, int m, int nv){
  for(int i = 0; i < n; i++){
    result[i] = 0;
    for(int j = 0; j < m; j++){
      result[i]+= matrix[i][j] * vector[j];
    }
  }
}

float dot_product(float *lv, float *rv, int n){
  float result = 0;
  for(int i = 0; i < n; i++){
    result += lv[i]*rv[i];
  }
  return result;
}

void product(float **result, float **a, float **b, int na, int ma, int nb, int mb){
  int i, j, k;
  float acc;
  if(ma != nb){
    printf("error: en el producto de matrices, el numero de columnas de la primera matriz debe ser igual al numero de filas de la segunda matriz");
    exit(1);
  }
  for(i = 0; i < ma; i++){
    for(j = 0; j < nb; j++){
      for(k = 0; k < ma; k++){
        result[i][j]+= a[i][k]*b[k][j];
      }
    }
  }
}

void print_vector(float *v, int n){
  for(int i = 0; i < n; i++){
    printf("%.2f ", v[i]);
  }
  printf("\n");
}


float findMax(float *c, int n){
  float max = FLT_MIN;
  for (int i = 0; i < n; i++){
    if (c[i] > max) max = c[i];
  }
  return max;
}

int findRowOfMin(float** matrix, float* b, int n, int j){
  float min = FLT_MAX;
  float current;
  int index = -1;
  for(int i = 0; i < n; i++){
    if(matrix[i][j] > 0){
      current = b[i] / matrix[i][j];
      if(current < min){
        min = current;
        index = i;
      }
    }
  }
  return index;
}

int indexOf(float value,float *c, int n){
  int i;
  for(i = 0; i < n; i++){
    if(c[i] == value) return i;
  }
  return -1;
}

void copyColunm(float* vector, float **matrix, int n, int j){
  int i;
  for(i = 0; i < n; i++){
    vector[i] = matrix[i][j];
  }
}

void setColumn(float **matrix, float *vector, int n, int j){
  for(int i = 0; i < n; i++){
    matrix[i][j] = vector[i];
  }
}


void print(float **matrix, int n, int m){
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      printf("%f ",matrix[i][j]);
    }
    printf("\n");
  }
}


void cofactor(float **matrix, float **aux, int p, int q, int n){
  int i=0, j=0;
  for(int k=0;k<n;k++){
    for(int l=0;l<n;l++){
      if(k != p && l != q){
        aux[i][j++] = matrix[k][l];
        if(j == n-1){
          j = 0;
          i++;
        }
      }
    }
  }
}

float determinant(float **matrix,int n){
  float det = 0.0;
  if(n==1) return matrix[0][0];
  float **temp;
  temp = (float**)malloc(n*sizeof(float*));
  for(int i=0;i<n;i++)  temp[i] = (float*)malloc(n*sizeof(float));
  int sign = 1;
  for(int k=0;k<n;k++){
    cofactor(matrix,temp,0,k,n);
    det += sign * matrix[0][k] * determinant(temp,n-1);
    sign *= -1;
  }
  for(int i=0;i<n;i++) free(temp[i]);
  free(temp);
  return det;
}

/*
float determinant(float **matrix, int n){
  float acc = 0, positiveProduct, negativeProduct;
  for(int i = 0; i < n; i++){
    positiveProduct = 1;
    negativeProduct = -1;
    for(int j = 0; j < n; j++){
      positiveProduct *= matrix[j][(i+j) % n];
      negativeProduct *= matrix[n-j-1][(i+j) % n];
    }
    acc += positiveProduct + negativeProduct;
  }
  return acc;
}
*/

void adj(float **matrix, float **adj,int n){
  if(n==1){
    adj[0][0] = 1;
    return;
  }
  int sign = 1;
  float **temp;
  temp = (float**)malloc(n*sizeof(float*));
  for(int i=0;i<n;i++)  temp[i] = (float*)malloc(n*sizeof(float));

  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      cofactor(matrix,temp,i,j,n);
      sign = ((i+j)%2==0) ? 1:-1;
      adj[j][i] = sign * determinant(temp,n-1);
    }
  }
  for(int i = 0; i < n; i++) free(temp[i]);
  free(temp);
}

int inverse(float **matrix, float **inverse,int n){
  float det = determinant(matrix,n);
  if(det==0) return 0;
  float **adju;
  adju = (float**)malloc(n*sizeof(float*));
  for(int i=0;i<n;i++) adju[i] = (float*)malloc(n*sizeof(float));
  
  adj(matrix,adju,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      inverse[i][j] = adju[i][j]/det;
    }
  }
  for(int i=0;i<n;i++) free(adju[i]);
  free(adju);
  return 1;
}