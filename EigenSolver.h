#ifndef EIGENSOLVER_H_
#define EIGENSOLVER_H_

#include "MathLib.h"

namespace RenderEngine
{
	class EigenSolver
	{
	public:
		static void Jacobi3(double eigenvalues[], double *eigenvectors, Vector3* pts, int n);

		static void Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors,
		                                 double *A, int n);
	};
}

#endif
