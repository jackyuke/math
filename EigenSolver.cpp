#include <float.h>                           // required for DBL_EPSILON
#include <math.h>                            // required for fabs()
#include <stdio.h>
#include "EigenSolver.h"

using namespace RenderEngine;

int jcbi(double a[], int n, double v[], double eps, int jt)
{
	int i,j,p = 0,q = 0,u,w,t,s,l;
	double fm,cn,sn,omega,x,y,d;
	l=1;
	for (i=0; i<=n-1; i++)
	{
		v[i*n+i]=1.0;
		for (j=0; j<=n-1; j++)
			if (i!=j) v[i*n+j]=0.0;
	}
	while (1==1)
	{
		fm=0.0;
		for (i=1; i<=n-1; i++)
			for (j=0; j<=i-1; j++)
			{
				d=fabs(a[i*n+j]);
				if ((i!=j)&&(d>fm))
				{
					fm=d;
					p=i;
					q=j;
				}
			}
		if (fm<eps)  return(1);
		if (l>jt)  return(-1);
		l=l+1;
		u=p*n+q;
		w=p*n+p;
		t=q*n+p;
		s=q*n+q;
		x=-a[u];
		y=(a[s]-a[w])/2.0;
		omega=x/sqrt(x*x+y*y);
		if (y<0.0) omega=-omega;
		sn=1.0+sqrt(1.0-omega*omega);
		sn=omega/sqrt(2.0*sn);
		cn=sqrt(1.0-sn*sn);
		fm=a[w];
		a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
		a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
		a[u]=0.0;
		a[t]=0.0;
		for (j=0; j<=n-1; j++)
			if ((j!=p)&&(j!=q))
			{
				u=p*n+j;
				w=q*n+j;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		for (i=0; i<=n-1; i++)
			if ((i!=p)&&(i!=q))
			{
				u=i*n+p;
				w=i*n+q;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		for (i=0; i<=n-1; i++)
		{
			u=i*n+p;
			w=i*n+q;
			fm=v[u];
			v[u]=fm*cn+v[w]*sn;
			v[w]=-fm*sn+v[w]*cn;
		}
	}
	return(1);
}

void EigenSolver::Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors,
                                       double *A, int n)
{
	int i, j, k, m;
	double *pAk, *pAm, *p_r, *p_e;
	double threshold_norm;
	double threshold;
	double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
	double sin_2phi, cos_2phi, cot_2phi;
	double dum1;
	double dum2;
	double dum3;
	double max;

	// Take care of trivial cases

	if ( n < 1) return;
	if ( n == 1)
	{
		eigenvalues[0] = *A;
		*eigenvectors = 1.0;
		return;
	}

	// Initialize the eigenvalues to the identity matrix.

	for (p_e = eigenvectors, i = 0; i < n; i++)
		for (j = 0; j < n; p_e++, j++)
			if (i == j) *p_e = 1.0;
			else *p_e = 0.0;

	// Calculate the threshold and threshold_norm.

	for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++)
		for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
	threshold = sqrt(threshold + threshold);
	threshold_norm = threshold * DBL_EPSILON;
	max = threshold + 1.0;
	while (threshold > threshold_norm)
	{
		threshold /= 10.0;
		if (max < threshold) continue;
		max = 0.0;
		for (pAk = A, k = 0; k < (n-1); pAk += n, k++)
		{
			for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++)
			{
				if ( fabs(*(pAk + m)) < threshold ) continue;

				// Calculate the sin and cos of the rotation angle which
				// annihilates A[k][m].

				cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
				dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
				if (cot_2phi < 0.0) dum1 = -dum1;
				tan_phi = -cot_2phi + dum1;
				tan2_phi = tan_phi * tan_phi;
				sin2_phi = tan2_phi / (1.0 + tan2_phi);
				cos2_phi = 1.0 - sin2_phi;
				sin_phi = sqrt(sin2_phi);
				if (tan_phi < 0.0) sin_phi = - sin_phi;
				cos_phi = sqrt(cos2_phi);
				sin_2phi = 2.0 * sin_phi * cos_phi;
				cos_2phi = cos2_phi - sin2_phi;

				// Rotate columns k and m for both the matrix A
				//     and the matrix of eigenvectors.

				p_r = A;
				dum1 = *(pAk + k);
				dum2 = *(pAm + m);
				dum3 = *(pAk + m);
				*(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
				*(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
				*(pAk + m) = 0.0;
				*(pAm + k) = 0.0;
				for (i = 0; i < n; p_r += n, i++)
				{
					if ( (i == k) || (i == m) ) continue;
					if ( i < k ) dum1 = *(p_r + k);
					else dum1 = *(pAk + i);
					if ( i < m ) dum2 = *(p_r + m);
					else dum2 = *(pAm + i);
					dum3 = dum1 * cos_phi + dum2 * sin_phi;
					if ( i < k ) *(p_r + k) = dum3;
					else *(pAk + i) = dum3;
					dum3 = - dum1 * sin_phi + dum2 * cos_phi;
					if ( i < m ) *(p_r + m) = dum3;
					else *(pAm + i) = dum3;
				}
				for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++)
				{
					dum1 = *(p_e + k);
					dum2 = *(p_e + m);
					*(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
					*(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
				}
			}
			for (i = 0; i < n; i++)
				if ( i == k ) continue;
				else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
		}
	}
	for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k);
}

#define SQ(x) ((x) * (x))
void EigenSolver::Jacobi3(double eigenvalues[], double *eigenvectors, Vector3* pts, int n)
{
	// calculate mean value
	Vector3 m = Vector3::Zero;
	for (int i = 0; i < n; i++)
	{
		m += pts[i];
	}

	m[0] = m[0] / n;
	m[1] = m[1] / n;
	m[2] = m[2] / n;

	//printf("%f, %f, %f\n", m[0], m[1], m[2]);

	// construct covariance matrix
	double a[9] = {0};

	for (int i = 0; i < n; i++)
	{
		a[0] += SQ((pts[i])[0] - m[0]);
		a[1] += ((pts[i])[0] - m[0]) * ((pts[i])[1] - m[1]);
		a[3] = a[1];
		a[2] += ((pts[i])[0] - m[0]) * ((pts[i])[2] - m[2]);
		a[6] = a[2];
		a[4] += SQ((pts[i])[1] - m[1]);
		a[8] += SQ((pts[i])[2] - m[2]);
		a[7] += ((pts[i])[2] - m[2]) * ((pts[i])[1] - m[1]);
		a[5] = a[7];

	}

	for (int i = 0; i < 9; i++)
		a[i] /= n;

	for (int i = 0; i < 9; i++)
		printf("%f \n", a[i]);
	//Jacobi_Cyclic_Method(eigenvalues, eigenvectors, a, n);
	jcbi(a, 3, eigenvectors, 0.000001, 100);

	/*for (int i = 0; i < 9; i++)
		printf("%f \n", a[i]);*/
}
