// dllmain.cpp : DLL 애플리케이션의 진입점을 정의합니다.
#include "pch.h"


/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under BSD 3-Clause license,
  * the "License"; You may not use this file except in compliance with the
  * License. You may obtain a copy of the License at:
  *                        opensource.org/licenses/BSD-3-Clause
  *
  ******************************************************************************
  */
  /* USER CODE END Header */
  /* Includes ------------------------------------------------------------------*/


//#include "TinyEKF.h"

/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */




  /* Cholesky-decomposition matrix-inversion code, adapated from
     http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static int choldc1(double* a, double* p, int n) {
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i * n + j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i * n + k] * a[j * n + k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j * n + i] = sum / p[i];
            }
        }
    }

    return 0; /* success */
}

static int choldcsl(double* A, double* a, double* p, int n)
{
    int i, j, k; double sum;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i * n + j] = A[i * n + j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i * n + i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j * n + k] * a[k * n + i];
            }
            a[j * n + i] = sum / p[j];
        }
    }

    return 0; /* success */
}


static int cholsl(double* A, double* a, double* p, int n)
{
    int i, j, k;
    if (choldcsl(A, a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i * n + j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i * n + i] *= a[i * n + i];
        for (k = i + 1; k < n; k++) {
            a[i * n + i] += a[k * n + i] * a[k * n + i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i * n + j] += a[k * n + i] * a[k * n + j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i * n + j] = a[j * n + i];
        }
    }

    return 0; /* success */
}

static void zeros(double* a, int m, int n)
{
    int j;
    for (j = 0; j < m * n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double* a, int m, int n, const char* fmt)
{
    int i, j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j)
            printf(f, a[i * n + j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double* a, double* b, double* c, int arows, int acols, int bcols)
{
    int i, j, l;

    for (i = 0; i < arows; ++i)
        for (j = 0; j < bcols; ++j) {
            c[i * bcols + j] = 0;
            for (l = 0; l < acols; ++l)
                c[i * bcols + j] += a[i * acols + l] * b[l * bcols + j];
        }
}

static void mulvec(double* a, double* x, double* y, int m, int n)
{
    int i, j;

    for (i = 0; i < m; ++i) {
        y[i] = 0;
        for (j = 0; j < n; ++j)
            y[i] += x[j] * a[i * n + j];
    }
}

static void transpose(double* a, double* at, int m, int n)
{
    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j) {
            at[j * m + i] = a[i * n + j];
        }
}

/* A <- A + B */
static void accum(double* a, double* b, int m, int n)
{
    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            a[i * n + j] += b[i * n + j];
}

/* C <- A + B */
static void add(double* a, double* b, double* c, int n)
{
    int j;

    for (j = 0; j < n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double* a, double* b, double* c, int n)
{
    int j;

    for (j = 0; j < n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(double* a, int m, int n)
{
    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            a[i * n + j] = -a[i * n + j];
}

static void mat_addeye(double* a, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        a[i * n + i] += 1;
}

/* TinyEKF code ------------------------------------------------------------------- */

//#include "tiny_ekf.h"

typedef struct {

    double* x;    /* state vector */

    double* P;  /* prediction error covariance */
    double* Q;  /* process noise covariance */
    double* R;  /* measurement error covariance */

    double* G;  /* Kalman gain; a.k.a. K */

    double* F;  /* Jacobian of process model */
    double* H;  /* Jacobian of measurement model */

    double* Ht; /* transpose of measurement Jacobian */
    double* Ft; /* transpose of process Jacobian */
    double* Pp; /* P, post-prediction, pre-update */

    double* fx;  /* output of user defined f() state-transition function */
    double* hx;  /* output of user defined h() measurement function */

    /* temporary storage */
    double* tmp0;
    double* tmp1;
    double* tmp2;
    double* tmp3;
    double* tmp4;
    double* tmp5;

} ekf_t2;

static void unpack(void* v, ekf_t2* ekf, int n, int m)
{
    /* skip over n, m in data structure */
    char* cptr = (char*)v;
    cptr += 2 * sizeof(int);

    double* dptr = (double*)cptr;
    ekf->x = dptr;
    dptr += n;
    ekf->P = dptr;
    dptr += n * n;
    ekf->Q = dptr;
    dptr += n * n;
    ekf->R = dptr;
    dptr += m * m;
    ekf->G = dptr;
    dptr += n * m;
    ekf->F = dptr;
    dptr += n * n;
    ekf->H = dptr;
    dptr += m * n;
    ekf->Ht = dptr;
    dptr += n * m;
    ekf->Ft = dptr;
    dptr += n * n;
    ekf->Pp = dptr;
    dptr += n * n;
    ekf->fx = dptr;
    dptr += n;
    ekf->hx = dptr;
    dptr += m;
    ekf->tmp0 = dptr;
    dptr += n * n;
    ekf->tmp1 = dptr;
    dptr += n * m;
    ekf->tmp2 = dptr;
    dptr += m * n;
    ekf->tmp3 = dptr;
    dptr += m * m;
    ekf->tmp4 = dptr;
    dptr += m * m;
    ekf->tmp5 = dptr;
}

void ekf_init(void* v, int n, int m)
{
    /* retrieve n, m and set them in incoming data structure */
    int* ptr = (int*)v;
    *ptr = n;
    ptr++;
    *ptr = m;

    /* unpack rest of incoming structure for initlization */
    ekf_t2 ekf;
    unpack(v, &ekf, n, m);

    /* zero-out matrices */
    zeros(ekf.P, n, n);
    zeros(ekf.Q, n, n);
    zeros(ekf.R, m, m);
    zeros(ekf.G, n, m);
    zeros(ekf.F, n, n);
    zeros(ekf.H, m, n);
}

int ekf_step1(void* v, double* z)
{
    /* unpack incoming structure */

    int* ptr = (int*)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t2 ekf;
    unpack(v, &ekf, n, m);

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);
    transpose(ekf.F, ekf.Ft, n, n);
    mulmat(ekf.tmp0, ekf.Ft, ekf.Pp, n, n, n);
    accum(ekf.Pp, ekf.Q, n, n);

    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(ekf.H, ekf.Ht, m, n);
    mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);
    mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);
    mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);
    accum(ekf.tmp3, ekf.R, m, m);
    if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
    mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, ekf.hx, ekf.tmp5, m);
    mulvec(ekf.G, ekf.tmp5, ekf.tmp2, n, m);
    ekf.tmp2[3] = 0;
    add(ekf.fx, ekf.tmp2, ekf.x, n);

    /* P_k = (I - G_k H_k) P_k */
    mulmat(ekf.G, ekf.H, ekf.tmp0, n, m, n);
    negate(ekf.tmp0, n, n);
    mat_addeye(ekf.tmp0, n);
    mulmat(ekf.tmp0, ekf.Pp, ekf.P, n, n, n);

    /* success */
    return 0;
}

int ekf_step2(void* v, double* z)
{
    /* unpack incoming structure */

    int* ptr = (int*)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t2 ekf;
    unpack(v, &ekf, n, m);

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);
    transpose(ekf.F, ekf.Ft, n, n);
    mulmat(ekf.tmp0, ekf.Ft, ekf.Pp, n, n, n);
    accum(ekf.Pp, ekf.Q, n, n);

    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(ekf.H, ekf.Ht, m, n);
    mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);
    mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);
    mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);
    accum(ekf.tmp3, ekf.R, m, m);
    if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
    mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, ekf.hx, ekf.tmp5, m);
    mulvec(ekf.G, ekf.tmp5, ekf.tmp2, n, m);
    ekf.tmp2[1] = 0;
    ekf.tmp2[2] = 0;
    add(ekf.fx, ekf.tmp2, ekf.x, n);

    /* P_k = (I - G_k H_k) P_k */
    mulmat(ekf.G, ekf.H, ekf.tmp0, n, m, n);
    negate(ekf.tmp0, n, n);
    mat_addeye(ekf.tmp0, n);
    mulmat(ekf.tmp0, ekf.Pp, ekf.P, n, n, n);

    /* success */
    return 0;
}

#define Nsta 4

/* observables */
#define Mobs 3

typedef struct {

    int n;          /* number of state values */
    int m;          /* number of observables */

    double x[Nsta];    /* state vector */

    double P[Nsta][Nsta];  /* prediction error covariance */
    double Q[Nsta][Nsta];  /* process noise covariance */
    double R[Mobs][Mobs];  /* measurement error covariance */

    double G[Nsta][Mobs];  /* Kalman gain; a.k.a. K */

    double F[Nsta][Nsta];  /* Jacobian of process model */
    double H[Mobs][Nsta];  /* Jacobian of measurement model */

    double Ht[Nsta][Mobs]; /* transpose of measurement Jacobian */
    double Ft[Nsta][Nsta]; /* transpose of process Jacobian */
    double Pp[Nsta][Nsta]; /* P, post-prediction, pre-update */

    double fx[Nsta];   /* output of user defined f() state-transition function */
    double hx[Mobs];   /* output of user defined h() measurement function */

    /* Personal Value */
//    double *xp[Nsta];
    /* temporary storage */
    double tmp0[Nsta][Nsta];
    double tmp1[Nsta][Mobs];
    double tmp2[Mobs][Nsta];
    double tmp3[Mobs][Mobs];
    double tmp4[Mobs][Mobs];
    double tmp5[Mobs];

} ekf_t;
/**
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constant N and
 * and then #include <TinyEKF.h>  You will also need to implement a model() method for your application.
 */
class TinyEKF {

private:

    ekf_t ekf;

protected:

    /**
      * The current state.
      */
    double* x;

    /**
     * Initializes a TinyEKF object.
     */
    TinyEKF() {
        ekf_init(&this->ekf, Nsta, Mobs);
        this->x = this->ekf.x;
    }

    /**
     * Deallocates memory for a TinyEKF object.
     */
    ~TinyEKF() { }

    /**
     * Implement this function for your EKF model.
     * @param fx gets output of state-transition function <i>f(x<sub>0 .. n-1</sub>)</i>
     * @param F gets <i>n &times; n</i> Jacobian of <i>f(x)</i>
     * @param hx gets output of observation function <i>h(x<sub>0 .. n-1</sub>)</i>
     * @param H gets <i>m &times; n</i> Jacobian of <i>h(x)</i>
     */
    virtual void model(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) = 0;
    virtual void model2(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) = 0;
    /**
     * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setP(int i, int j, double value)
    {
        this->ekf.P[i][j] = value;
    }

    /**
     * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setQ(int i, int j, double value)
    {
        this->ekf.Q[i][j] = value;
    }

    /**
     * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setR(int i, int j, double value)
    {
        this->ekf.R[i][j] = value;
    }

public:

    /**
     * Returns the state element at a given index.
     * @param i the index (at least 0 and less than <i>n</i>
     * @return state value at index
     */
    double getX(int i)
    {
        return this->ekf.x[i];
    }

    /**
     * Sets the state element at a given index.
     * @param i the index (at least 0 and less than <i>n</i>
     * @param value value to set
     */
    void setX(int i, double value)
    {
        this->ekf.x[i] = value;
    }

    /**
      Performs one step of the prediction and update.
     * @param z observation vector, length <i>m</i>
     * @return true on success, false on failure caused by non-positive-definite matrix.
     */
    bool step(double* z)
    {
        this->model(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H);

        return ekf_step1(&this->ekf, z) ? false : true;
    }
    bool step2(double* z)
    {
        this->model2(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H);

        return ekf_step2(&this->ekf, z) ? false : true;
    }
};

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/


/* USER CODE BEGIN PFP */

/* USER Global Value */
double a[3], w[3], h[3];
double pi, theta, psi = 0.0;
double t = 0.01;
double g = 9.81;
double S_gps_x, S_gps_y, S_gps_z = 0.0;// initial location
double V_x, V_y, V_z = 0.0; // initial velocity
bool returns; // KF's return (T/F)
double return_x[4]; // KF's return value (quarternion)
double xp[4]; // model2's initial value (model1's result)

/* Experiment TXT */
double S[100][3] = { 0.0, };
double A[100][3] = { 0.0, };
double W[100][3] = { 0.0, };
double H[100][3] = { 0.0, };


/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
class IMUEKF : public TinyEKF {

public:

    IMUEKF() {
        // state vector init = Quaternion
        this->setX(0, 0.5);
        this->setX(1, 0.5);
        this->setX(2, 0.5);
        this->setX(3, 0.5);

        // error cov
        this->setP(0, 0, 0.0055);
        this->setP(1, 1, 0.0055);
        this->setP(2, 2, 0.0055);
        this->setP(3, 3, 0.0055);

        // process noise
        this->setQ(0, 0, pow(0.4, 2));
        this->setQ(1, 1, pow(0.4, 2));
        this->setQ(2, 2, pow(0.4, 2));
        this->setQ(3, 3, pow(0.4, 2));

        // measurement noise
        this->setR(0, 0, 0.01); // Range Noise
        this->setR(1, 1, 100); // Angle Noise
        this->setR(2, 2, pow(0.4, 2));
        this->setR(3, 3, pow(0.4, 2));
    }

protected:


    void model(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) {

        // Process model hat[X_predict]
        fx[0] = this->x[0] + (this->x[1] * (-w[0] * t / 2)) + (this->x[2] * (-w[1] * t / 2)) + (this->x[3] * (-w[2] * t / 2));
        fx[1] = this->x[1] + (this->x[0] * (w[0] * t / 2)) + (this->x[2] * (w[2] * t / 2)) + (this->x[3] * (-w[1] * t / 2));
        fx[2] = this->x[2] + (this->x[0] * (-w[1] * t / 2)) + (this->x[1] * (-w[2] * t / 2)) + (this->x[3] * (w[0] * t / 2));
        fx[3] = this->x[3] + (this->x[0] * (w[2] * t / 2)) + (this->x[1] * (w[1] * t / 2)) + (this->x[2] * (-w[0] * t / 2));

        //state transition matrix
        F[0][0] = 1;
        F[0][1] = -w[0] * t / 2;
        F[0][2] = -w[1] * t / 2;
        F[0][3] = -w[2] * t / 2;

        F[1][0] = w[0] * t / 2;
        F[1][1] = 1;
        F[1][2] = w[2] * t / 2;
        F[1][3] = -w[1] * t / 2;

        F[2][0] = -w[1] * t / 2;
        F[2][1] = -w[2] * t / 2;
        F[2][2] = 1;
        F[2][3] = w[0] * t / 2;

        F[3][0] = w[2] * t / 2;
        F[3][1] = w[1] * t / 2;
        F[3][2] = -w[0] * t / 2;
        F[3][3] = 1;


        // Measurement function simplifies the relationship between state and sensor readings for convenience.
        // A more realistic measurement function would distinguish between state value and measured value; e.g.:

        // step 1 hx
        hx[0] = (2 * fx[1] * fx[3] - 2 * fx[0] * fx[2]) * 9.8;
        hx[1] = (2 * fx[0] * fx[1] + 2 * fx[2] * fx[3]) * 9.8;
        hx[2] = (pow(fx[0], 2) - pow(fx[1], 2) - pow(fx[2], 2) + pow(fx[3], 2)) * 9.8;


        // step 1 Jacobian of measurement function
        H[0][0] = -2 * fx[2];
        H[0][1] = 2 * fx[3];
        H[0][2] = -2 * fx[0];
        H[0][3] = 2 * fx[1];

        H[1][0] = 2 * fx[1];
        H[1][1] = 2 * fx[0];
        H[1][2] = 2 * fx[3];
        H[1][3] = 2 * fx[2];

        H[2][0] = 2 * fx[0];
        H[2][1] = -2 * fx[1];
        H[2][2] = -2 * fx[2];
        H[2][3] = 2 * fx[3];

    }

    void model2(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) {}

};


class IMUEKF2 : public TinyEKF {

public:

    IMUEKF2() {
        // state vector init = Quaternion
        this->setX(0, xp[0]);
        this->setX(1, xp[1]);
        this->setX(2, xp[2]);
        this->setX(3, xp[3]);

        // error cov
        this->setP(0, 0, 0.0055);
        this->setP(1, 1, 0.0055);
        this->setP(2, 2, 0.0055);
        this->setP(3, 3, 0.0055);

        // process noise
        this->setQ(0, 0, pow(0.4, 2));
        this->setQ(1, 1, pow(0.4, 2));
        this->setQ(2, 2, pow(0.4, 2));
        this->setQ(3, 3, pow(0.4, 2));

        // measurement noise
        this->setR(1, 1, 100); // Angle Noise
        this->setR(2, 2, pow(0.4, 2));
        this->setR(3, 3, pow(0.4, 2));
    }

protected:

    void model(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) {}

    void model2(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta]) {
        // Process model hat[X_predict]
        fx[0] = this->x[0] + (this->x[1] * (-w[0] * t / 2)) + (this->x[2] * (-w[1] * t / 2)) + (this->x[3] * (-w[2] * t / 2));
        fx[1] = this->x[0] + (this->x[0] * (w[0] * t / 2)) + (this->x[2] * (w[2] * t / 2)) + (this->x[3] * (-w[1] * t / 2));
        fx[2] = this->x[0] + (this->x[0] * (-w[1] * t / 2)) + (this->x[1] * (-w[2] * t / 2)) + (this->x[3] * (w[0] * t / 2));
        fx[3] = this->x[0] + (this->x[0] * (w[2] * t / 2)) + (this->x[1] * (w[1] * t / 2)) + (this->x[2] * (-w[0] * t / 2));

        //state transition matrix
        F[0][0] = 1;
        F[0][1] = -w[0] * t / 2;
        F[0][2] = -w[1] * t / 2;
        F[0][3] = -w[2] * t / 2;

        F[1][0] = w[0] * t / 2;
        F[1][1] = 1;
        F[1][2] = w[2] * t / 2;
        F[1][3] = -w[1] * t / 2;

        F[2][0] = -w[1] * t / 2;
        F[2][1] = -w[2] * t / 2;
        F[2][2] = 1;
        F[2][3] = w[0] * t / 2;

        F[3][0] = w[2] * t / 2;
        F[3][1] = w[1] * t / 2;
        F[3][2] = -w[0] * t / 2;
        F[3][3] = 1;

        // Measurement function simplifies the relationship between state and sensor readings for convenience.
        // A more realistic measurement function would distinguish between state value and measured value; e.g.:

        // step 2 hx
        hx[0] = 2 * fx[1] * fx[2] + 2 * fx[0] * fx[3];
        hx[1] = pow(fx[0], 2) - pow(fx[1], 2) - pow(fx[2], 2) - pow(fx[3], 2);
        hx[2] = 2 * fx[2] * fx[3] - 2 * fx[0] * fx[1];


        // step 1 Jacobian of measurement function
        H[0][0] = 2 * fx[3];
        H[0][1] = 2 * fx[2];
        H[0][2] = 2 * fx[1];
        H[0][3] = 2 * fx[0];

        H[1][0] = 2 * fx[0];
        H[1][1] = -2 * fx[1];
        H[1][2] = -2 * fx[2];
        H[1][3] = 2 * fx[3];

        H[2][0] = -2 * fx[1];
        H[2][1] = -2 * fx[0];
        H[2][2] = 2 * fx[3];
        H[2][3] = 2 * fx[2];

    }
};

void quaternion_to_euler(double q0, double q1, double q2, double q3) {

    pi = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (pow(q1, 2.0) + pow(q2, 2.0)));
    theta = asin(2 * (q0 * q2 - q1 * q3) / (pow(q0, 2.0) + pow(q1, 2.0) + pow(q2, 2.0) + pow(q3, 2.0)));
    psi = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (pow(q2, 2.0) + pow(q3, 2.0)));

}


typedef double(*arrPointer)[3];

arrPointer make_string(arrPointer AA, arrPointer WW, arrPointer HH)
{

    IMUEKF ekf;
    IMUEKF2 ekf2;

    for (int i = 0; i < 100; i++) {

        a[0] = AA[i][0];
        a[1] = AA[i][1];
        a[2] = AA[i][2];

        w[0] = WW[i][0];
        w[1] = WW[i][1];
        w[2] = WW[i][2];

        h[0] = HH[i][0];
        h[1] = HH[i][1];
        h[2] = HH[i][2];



        double FilterData[3] = { a[0], a[1], a[2] };
        returns = ekf.step(FilterData);

        ekf2.setX(0, ekf.getX(0));
        ekf2.setX(1, ekf.getX(1));
        ekf2.setX(2, ekf.getX(2));
        ekf2.setX(3, ekf.getX(3));

        double FilterData2[3] = { h[0], h[1], h[2] };
        returns = ekf2.step2(FilterData2);


        return_x[0] = ekf2.getX(0);
        return_x[1] = ekf2.getX(1);
        return_x[2] = ekf2.getX(2);
        return_x[3] = ekf2.getX(3);

        // Update Roll, Pitch, Yaw by return x (quaternion)
        quaternion_to_euler(return_x[0], return_x[1], return_x[2], return_x[3]);

        double acc_x = a[0] + w[1] * (V_z)-w[2] * (V_y)-g * sin(theta);
        double acc_y = a[1] + w[0] * (V_z)+w[2] * (V_x)+g * sin(pi) * cos(theta);
        double acc_z = a[2] + w[0] * (V_y)-w[1] * (V_x)+g * cos(pi) * cos(theta);

        double acc_x_world = cos(theta) * cos(pi) * acc_x + cos(theta) * sin(pi) * acc_y - sin(theta) * acc_z;
        double acc_y_world = (sin(psi) * sin(theta) * cos(pi) - cos(psi) * sin(pi)) * acc_x + (sin(pi) * sin(theta) * sin(psi) + cos(pi) * cos(psi)) * acc_y + sin(psi) * cos(theta) * acc_z;
        double acc_z_world = (cos(pi) * sin(theta) * cos(psi) + sin(pi) * sin(psi)) * acc_x + (cos(psi) * sin(theta) * sin(pi) - sin(psi) * cos(pi)) * acc_y + cos(psi) * cos(theta) * acc_z;

        V_x = V_x + acc_x_world * t;
        V_y = V_y + acc_y_world * t;
        V_z = V_z + acc_z_world * t;

        S_gps_x = S_gps_x + V_x * t;
        S_gps_y = S_gps_y + V_y * t;
        S_gps_z = S_gps_z + V_z * t;

        S[i][0] = S_gps_x;
        S[i][1] = S_gps_y;
        S[i][2] = S_gps_z;
     }

    //std::string str;
    //str = "[ ";
    //for (int idx = 0; idx < 100; idx++) {
    //    str.append(std::to_string(S[idx][0]));
    //    str.append(", ");
    //    str.append(std::to_string(S[idx][1]));
    //    str.append(", ");
    //    str.append(std::to_string(S[idx][2]));
    //    if (idx < 99) {
    //        str.append(", ");
    //    }
    //}
    //str.append("]");
    //return str;
    return S;
}



//scp -P 7022 -r Dll1.dll byungkeun@147.47.208.20:/home/byungkeun/PROGRAMS




////      HAL_UART_Transmit(&huart2, imubufftx, 45, 100);  // send received data
//   }
////      memset(imubufftx, 0, sizeof(imubufftx));
//  }

//   if (RcvStat == HAL_OK) {  // receive check
//      HAL_UART_Transmit(&huart2, UsartData, 1, 100) ;  // send received data
//   }

  /* USER CODE END 3 */


/**
  * @brief System Clock Configuration
  * @retval None
  */


  /************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
