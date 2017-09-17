#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MathLib.h"

namespace RenderEngine
{
	class Quaternion
	{
	public:
		inline Quaternion ();
		inline Quaternion (double fW, double fX, double fY, double fZ);
		inline Quaternion (const Quaternion& rkQ);

		inline Quaternion (const Vector3& rkAxis, double fAngle);

		Quaternion (const Vector3 akRotColumn[3]);

		operator const double* () const;
		operator double* ();
		double operator[] (int i) const;
		double& operator[] (int i);
		double W () const;
		double& W ();
		double X () const;
		double& X ();
		double Y () const;
		double& Y ();
		double Z () const;
		double& Z ();

		Quaternion& operator= (const Quaternion& rkQ);

		bool operator== (const Quaternion& rkQ) const;
		bool operator!= (const Quaternion& rkQ) const;
		bool operator<  (const Quaternion& rkQ) const;
		bool operator<= (const Quaternion& rkQ) const;
		bool operator>  (const Quaternion& rkQ) const;
		bool operator>= (const Quaternion& rkQ) const;
		Quaternion operator+ (const Quaternion& rkQ) const;
		Quaternion operator- (const Quaternion& rkQ) const;
		Quaternion operator* (const Quaternion& rkQ) const;
		Quaternion operator* (double fScalar) const;
		Quaternion operator/ (double fScalar) const;
		Quaternion operator- () const;
		Quaternion& operator+= (const Quaternion& rkQ);
		Quaternion& operator-= (const Quaternion& rkQ);
		Quaternion& operator*= (double fScalar);
		Quaternion& operator/= (double fScalar);
		inline Quaternion& FromRotationMatrix (const Vector3 akRotColumn[3]);
		void ToRotationMatrix ( Vector3 akRotColumn[3] ) const;
		void ToRotationMatrix ( double rkRot[9] ) const;
		inline Quaternion& FromAxisAngle ( const Vector3& rkAxis, double fAngle);
		void ToAxisAngle ( Vector3& rkAxis, double& rfAngle) const;
		double Length () const;
		double SquaredLength () const;
		double Dot (const Quaternion& rkQ) const;
		double Normalize ();
		Quaternion Inverse () const;
		Quaternion Conjugate () const;
		Quaternion Exp () const;
		Quaternion Log () const;
		Vector3 Rotate (const Vector3& rkVector) const;
		Quaternion& Align (const Vector3& rkV1, const Vector3& rkV2);
		Quaternion& Slerp (double fT, const Quaternion& rkP, const Quaternion& rkQ);
		Quaternion& SlerpExtraSpins (double fT, const Quaternion& rkP, const Quaternion& rkQ, int iExtraSpins);

		static const Quaternion IDENTITY;
		static const Quaternion ZERO;

	private:

		int CompareArrays (const Quaternion& rkQ) const;

		static int ms_iNext[3];

		double m_afTuple[4];
	};

	inline Quaternion operator* (double fScalar, const Quaternion& rkQ);

	//----------------------------------------------------------------------------
	inline Quaternion::Quaternion ()
	{
		// uninitialized for performance in array construction
	}
	//----------------------------------------------------------------------------

	inline Quaternion::Quaternion (double fW, double fX, double fY, double fZ)
	{
		m_afTuple[0] = fW;
		m_afTuple[1] = fX;
		m_afTuple[2] = fY;
		m_afTuple[3] = fZ;
	}
	//----------------------------------------------------------------------------

	inline Quaternion::Quaternion (const Quaternion& rkQ)
	{
		size_t uiSize = 4 * sizeof(double);
		memcpy(m_afTuple, rkQ.m_afTuple, uiSize);
	}
	//----------------------------------------------------------------------------

	inline Quaternion::Quaternion (const Vector3& rkAxis, double fAngle)
	{
		FromAxisAngle(rkAxis,fAngle);
	}
	//----------------------------------------------------------------------------

	inline Quaternion::Quaternion (const Vector3 akRotColumn[3])
	{
		FromRotationMatrix(akRotColumn);
	}

	//----------------------------------------------------------------------------

	inline Quaternion::operator const double* () const
	{
		return m_afTuple;
	}
	//----------------------------------------------------------------------------

	inline Quaternion::operator double* ()
	{
		return m_afTuple;
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::operator[] (int i) const
	{
		return m_afTuple[i];
	}
	//----------------------------------------------------------------------------

	inline double& Quaternion::operator[] (int i)
	{
		return m_afTuple[i];
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::W () const
	{
		return m_afTuple[0];
	}
	//----------------------------------------------------------------------------

	inline double& Quaternion::W ()
	{
		return m_afTuple[0];
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::X () const
	{
		return m_afTuple[1];
	}
	//----------------------------------------------------------------------------

	inline double& Quaternion::X ()
	{
		return m_afTuple[1];
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::Y () const
	{
		return m_afTuple[2];
	}
	//----------------------------------------------------------------------------

	inline double& Quaternion::Y ()
	{
		return m_afTuple[2];
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::Z () const
	{
		return m_afTuple[3];
	}
	//----------------------------------------------------------------------------

	inline double& Quaternion::Z ()
	{
		return m_afTuple[3];
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::operator= (const Quaternion& rkQ)
	{
		size_t uiSize = 4*sizeof(double);
		memcpy(m_afTuple, rkQ.m_afTuple, uiSize);
		return *this;
	}
	//----------------------------------------------------------------------------

	inline int Quaternion::CompareArrays (const Quaternion& rkQ) const
	{
		return memcmp(m_afTuple,rkQ.m_afTuple,4*sizeof(double));
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator== (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) == 0;
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator!= (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) != 0;
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator< (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) < 0;
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator<= (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) <= 0;
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator> (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) > 0;
	}
	//----------------------------------------------------------------------------

	inline bool Quaternion::operator>= (const Quaternion& rkQ) const
	{
		return CompareArrays(rkQ) >= 0;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator+ (const Quaternion& rkQ) const
	{
		Quaternion kSum;
		for (int i = 0; i < 4; i++)
		{
			kSum.m_afTuple[i] = m_afTuple[i] + rkQ.m_afTuple[i];
		}
		return kSum;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator- (const Quaternion& rkQ) const
	{
		Quaternion kdiff;
		for (int i = 0; i < 4; i++)
		{
			kdiff.m_afTuple[i] = m_afTuple[i] - rkQ.m_afTuple[i];
		}
		return kdiff;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator* (const Quaternion& rkQ) const
	{
		// NOTE:  Multiplication is not generally commutative, so in most
		// cases p*q != q*p.

		Quaternion kProd;

		kProd.m_afTuple[0] =
		    m_afTuple[0]*rkQ.m_afTuple[0] -
		    m_afTuple[1]*rkQ.m_afTuple[1] -
		    m_afTuple[2]*rkQ.m_afTuple[2] -
		    m_afTuple[3]*rkQ.m_afTuple[3];

		kProd.m_afTuple[1] =
		    m_afTuple[0]*rkQ.m_afTuple[1] +
		    m_afTuple[1]*rkQ.m_afTuple[0] +
		    m_afTuple[2]*rkQ.m_afTuple[3] -
		    m_afTuple[3]*rkQ.m_afTuple[2];

		kProd.m_afTuple[2] =
		    m_afTuple[0]*rkQ.m_afTuple[2] +
		    m_afTuple[2]*rkQ.m_afTuple[0] +
		    m_afTuple[3]*rkQ.m_afTuple[1] -
		    m_afTuple[1]*rkQ.m_afTuple[3];

		kProd.m_afTuple[3] =
		    m_afTuple[0]*rkQ.m_afTuple[3] +
		    m_afTuple[3]*rkQ.m_afTuple[0] +
		    m_afTuple[1]*rkQ.m_afTuple[2] -
		    m_afTuple[2]*rkQ.m_afTuple[1];

		return kProd;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator* (double fScalar) const
	{
		Quaternion kProd;
		for (int i = 0; i < 4; i++)
		{
			kProd.m_afTuple[i] = fScalar*m_afTuple[i];
		}
		return kProd;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator/ (double fScalar) const
	{
		Quaternion kQuot;
		int i;

		if (fScalar != (double)0.0)
		{
			double fInvScalar = ((double)1.0)/fScalar;
			for (i = 0; i < 4; i++)
			{
				kQuot.m_afTuple[i] = fInvScalar*m_afTuple[i];
			}
		}
		else
		{
			for (i = 0; i < 4; i++)
			{
				kQuot.m_afTuple[i] = (double)99999999.0f;
			}
		}

		return kQuot;
	}
	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::operator- () const
	{
		Quaternion kNeg;
		for (int i = 0; i < 4; i++)
		{
			kNeg.m_afTuple[i] = -m_afTuple[i];
		}
		return kNeg;
	}
	//----------------------------------------------------------------------------

	inline Quaternion operator* (double fScalar, const Quaternion& rkQ)
	{
		Quaternion kProd;
		for (int i = 0; i < 4; i++)
		{
			kProd[i] = fScalar*rkQ[i];
		}
		return kProd;
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::operator+= (const Quaternion& rkQ)
	{
		for (int i = 0; i < 4; i++)
		{
			m_afTuple[i] += rkQ.m_afTuple[i];
		}
		return *this;
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::operator-= (const Quaternion& rkQ)
	{
		for (int i = 0; i < 4; i++)
		{
			m_afTuple[i] -= rkQ.m_afTuple[i];
		}
		return *this;
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::operator*= (double fScalar)
	{
		for (int i = 0; i < 4; i++)
		{
			m_afTuple[i] *= fScalar;
		}
		return *this;
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::operator/= (double fScalar)
	{
		int i;

		if (fScalar != (double)0.0)
		{
			double fInvScalar = ((double)1.0)/fScalar;
			for (i = 0; i < 4; i++)
			{
				m_afTuple[i] *= fInvScalar;
			}
		}
		else
		{
			for (i = 0; i < 4; i++)
			{
				m_afTuple[i] = (double)99999999.0f;
			}
		}

		return *this;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	inline void Quaternion::ToRotationMatrix ( Vector3 rkRot[3] ) const
	{
		double fTx  = ((double)2.0)*m_afTuple[1];
		double fTy  = ((double)2.0)*m_afTuple[2];
		double fTz  = ((double)2.0)*m_afTuple[3];
		double fTwx = fTx*m_afTuple[0];
		double fTwy = fTy*m_afTuple[0];
		double fTwz = fTz*m_afTuple[0];
		double fTxx = fTx*m_afTuple[1];
		double fTxy = fTy*m_afTuple[1];
		double fTxz = fTz*m_afTuple[1];
		double fTyy = fTy*m_afTuple[2];
		double fTyz = fTz*m_afTuple[2];
		double fTzz = fTz*m_afTuple[3];

		rkRot[0][0] = (double)1.0-(fTyy+fTzz);
		rkRot[0][1] = fTxy-fTwz;
		rkRot[0][2] = fTxz+fTwy;

		rkRot[1][0] = fTxy+fTwz;
		rkRot[1][1] = (double)1.0-(fTxx+fTzz);
		rkRot[1][2] = fTyz-fTwx;

		rkRot[2][0] = fTxz-fTwy;
		rkRot[2][1] = fTyz+fTwx;
		rkRot[2][2] = (double)1.0-(fTxx+fTyy);
	}

	//----------------------------------------------------------------------------
	inline void Quaternion::ToRotationMatrix ( double rkRot[9] ) const
	{
		double fTx  = ((double)2.0)*m_afTuple[1];
		double fTy  = ((double)2.0)*m_afTuple[2];
		double fTz  = ((double)2.0)*m_afTuple[3];
		double fTwx = fTx*m_afTuple[0];
		double fTwy = fTy*m_afTuple[0];
		double fTwz = fTz*m_afTuple[0];
		double fTxx = fTx*m_afTuple[1];
		double fTxy = fTy*m_afTuple[1];
		double fTxz = fTz*m_afTuple[1];
		double fTyy = fTy*m_afTuple[2];
		double fTyz = fTz*m_afTuple[2];
		double fTzz = fTz*m_afTuple[3];

		rkRot[0] = (double)1.0-(fTyy+fTzz);
		rkRot[1] = fTxy-fTwz;
		rkRot[2] = fTxz+fTwy;

		rkRot[3] = fTxy+fTwz;
		rkRot[4] = (double)1.0-(fTxx+fTzz);
		rkRot[5] = fTyz-fTwx;

		rkRot[6] = fTxz-fTwy;
		rkRot[7] = fTyz+fTwx;
		rkRot[8] = (double)1.0-(fTxx+fTyy);
	}

	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::FromRotationMatrix (const Vector3 rkRot[3])
	{
		// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
		// article "Quaternion Calculus and Fast Animation".

		double fTrace = rkRot[0][0] + rkRot[1][1] + rkRot[2][2];
		double fRoot;

		if (fTrace > 0.0)
		{
			// |w| > 1/2, may as well choose w > 1/2
			fRoot = sqrt(fTrace + 1.0);  // 2w
			m_afTuple[0] = 0.5*fRoot;
			fRoot = 0.5/fRoot;  // 1/(4w)
			m_afTuple[1] = (rkRot[2][1]-rkRot[1][2])*fRoot;
			m_afTuple[2] = (rkRot[0][2]-rkRot[2][0])*fRoot;
			m_afTuple[3] = (rkRot[1][0]-rkRot[0][1])*fRoot;
		}
		else
		{
			// |w| <= 1/2
			int i = 0;
			if (rkRot[1][1] > rkRot[0][0])
			{
				i = 1;
			}
			if (rkRot[2][2] > rkRot[i][i])
			{
				i = 2;
			}
			int j = ms_iNext[i];
			int k = ms_iNext[j];

			fRoot = sqrt(rkRot[i][i]-rkRot[j][j]-rkRot[k][k]+1.0);
			double* apfQuat[3] = { &m_afTuple[1], &m_afTuple[2], &m_afTuple[3] };
			*apfQuat[i] = 0.5*fRoot;
			fRoot = 0.5/fRoot;
			m_afTuple[0] = (rkRot[k][j]-rkRot[j][k])*fRoot;
			*apfQuat[j] = (rkRot[j][i]+rkRot[i][j])*fRoot;
			*apfQuat[k] = (rkRot[k][i]+rkRot[i][k])*fRoot;
		}

		return *this;
	}

	/*
	//----------------------------------------------------------------------------

	inline void Quaternion::ToRotationMatrix (Vector3 akRotColumn[3]) const
	{
		hsl::matrix4<double> kRot;
		ToRotationMatrix(kRot);
		for (int iCol = 0; iCol < 3; iCol++)
		{
			akRotColumn[iCol][0] = kRot(0,iCol);
			akRotColumn[iCol][1] = kRot(1,iCol);
			akRotColumn[iCol][2] = kRot(2,iCol);
		}
	}
	*/

	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::FromAxisAngle (const Vector3 & rkAxis, double fAngle)
	{
		// assert:  axis[] is unit length
		//
		// The quaternion representing the rotation is
		//   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

		double fHalfAngle = ((double)0.5)*fAngle;
		double fSin = (double)sin(fHalfAngle);
		m_afTuple[0] = (double)cos(fHalfAngle);
		m_afTuple[1] = fSin*rkAxis.x;
		m_afTuple[2] = fSin*rkAxis.y;
		m_afTuple[3] = fSin*rkAxis.z;

		return *this;
	}
	//----------------------------------------------------------------------------

	inline void Quaternion::ToAxisAngle (Vector3& rkAxis, double& rfAngle) const
	{
		// The quaternion representing the rotation is
		//   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

		double fSqrLength = m_afTuple[1]*m_afTuple[1] + m_afTuple[2]*m_afTuple[2]
		                    + m_afTuple[3]*m_afTuple[3];
		if (fSqrLength > (double)0.00001)
		{
			rfAngle = ((double)2.0)*acos(Math::Clamp(m_afTuple[0], -1.0, 1.0));
			double fInvLength = 1.0 / sqrt( fSqrLength);
			rkAxis[0] = m_afTuple[1]*fInvLength;
			rkAxis[1] = m_afTuple[2]*fInvLength;
			rkAxis[2] = m_afTuple[3]*fInvLength;
		}
		else
		{
			// angle is 0 (mod 2*pi), so any axis will do
			rfAngle = (double)0.0;
			rkAxis[0] = (double)1.0;
			rkAxis[1] = (double)0.0;
			rkAxis[2] = (double)0.0;
		}
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::Length () const
	{
		return (double)sqrt(
		           m_afTuple[0]*m_afTuple[0] +
		           m_afTuple[1]*m_afTuple[1] +
		           m_afTuple[2]*m_afTuple[2] +
		           m_afTuple[3]*m_afTuple[3]);
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::SquaredLength () const
	{
		return
		    m_afTuple[0]*m_afTuple[0] +
		    m_afTuple[1]*m_afTuple[1] +
		    m_afTuple[2]*m_afTuple[2] +
		    m_afTuple[3]*m_afTuple[3];
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::Dot (const Quaternion& rkQ) const
	{
		double fDot = (double)0.0;
		for (int i = 0; i < 4; i++)
		{
			fDot += m_afTuple[i]*rkQ.m_afTuple[i];
		}
		return fDot;
	}
	//----------------------------------------------------------------------------

	inline double Quaternion::Normalize ()
	{
		double fLength = Length();
		int i;

		if (fLength > (double)0.00001)
		{
			double fInvLength = ((double)1.0)/fLength;
			for (i = 0; i < 4; i++)
			{
				m_afTuple[i] *= fInvLength;
			}
		}
		else
		{
			fLength = (double)0.0;
			for (i = 0; i < 4; i++)
			{
				m_afTuple[i] = (double)0.0;
			}
		}

		return fLength;
	}

	//----------------------------------------------------------------------------

	inline Quaternion Quaternion::Inverse () const
	{
		Quaternion kInverse;

		double fNorm = 0.0;
		int i;
		for (i = 0; i < 4; i++)
		{
			fNorm += m_afTuple[i]*m_afTuple[i];
		}

		if (fNorm > 0.0)
		{
			double fInvNorm = 1.0/fNorm;
			kInverse.m_afTuple[0] = m_afTuple[0]*fInvNorm;
			kInverse.m_afTuple[1] = -m_afTuple[1]*fInvNorm;
			kInverse.m_afTuple[2] = -m_afTuple[2]*fInvNorm;
			kInverse.m_afTuple[3] = -m_afTuple[3]*fInvNorm;
		}
		else
		{
			// return an invalid result to flag the error
			for (i = 0; i < 4; i++)
			{
				kInverse.m_afTuple[i] = 0.0;
			}
		}

		return kInverse;
	}

	//----------------------------------------------------------------------------

	inline Vector3 Quaternion::Rotate (const Vector3& rkVector) const
	{
		// Given a vector u = (x0,y0,z0) and a unit length quaternion
		// q = <w,x,y,z>, the vector v = (x1,y1,z1) which represents the
		// rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
		// multiplication and where u is treated as the quaternion <0,x0,y0,z0>.
		// Note that q^{-1} = <w,-x,-y,-z>, so no double work is required to
		// invert q.  Now
		//
		//   q*u*q^{-1} = q*<0,x0,y0,z0>*q^{-1}
		//     = q*(x0*i+y0*j+z0*k)*q^{-1}
		//     = x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
		//
		// As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and 2*k*q^{-1} are the columns
		// of the rotation matrix computed in Quaternion::ToRotationMatrix.
		// The vector v is obtained as the product of that rotation matrix with
		// vector u.  As such, the quaternion representation of a rotation
		// matrix requires less space than the matrix and more time to compute
		// the rotated vector.  Typical space-time tradeoff...

		double kRot[9];
		ToRotationMatrix( kRot );
		Vector3 result;

		result.x = rkVector.Dot( Vector3( kRot[0], kRot[3], kRot[6] ) );
		result.y = rkVector.Dot( Vector3( kRot[1], kRot[4], kRot[7] ) );
		result.z = rkVector.Dot( Vector3( kRot[2], kRot[5], kRot[8] ) );

		return result;
	}
	//----------------------------------------------------------------------------

	inline Quaternion& Quaternion::Align (const Vector3& rkV1, const Vector3& rkV2)
	{
		// If V1 and V2 are not parallel, the axis of rotation is the unit-length
		// vector U = Cross(V1,V2)/Length(Cross(V1,V2)).  The angle of rotation,
		// A, is the angle between V1 and V2.  The quaternion for the rotation is
		// q = cos(A/2) + sin(A/2)*(ux*i+uy*j+uz*k) where U = (ux,uy,uz).
		//
		// (1) Rather than extract A = acos(Dot(V1,V2)), multiply by 1/2, then
		//     compute sin(A/2) and cos(A/2), we reduce the computational costs by
		//     computing the bisector B = (V1+V2)/Length(V1+V2), so cos(A/2) =
		//     Dot(V1,B).
		//
		// (2) The rotation axis is U = Cross(V1,B)/Length(Cross(V1,B)), but
		//     Length(Cross(V1,B)) = Length(V1)*Length(B)*sin(A/2) = sin(A/2), in
		//     which case sin(A/2)*(ux*i+uy*j+uz*k) = (cx*i+cy*j+cz*k) where
		//     C = Cross(V1,B).
		//
		// If V1 = V2, then B = V1, cos(A/2) = 1, and U = (0,0,0).  If V1 = -V2,
		// then B = 0.  This can happen even if V1 is approximately -V2 using
		// floating point arithmetic, since Vector3::Normalize checks for
		// closeness to zero and returns the zero vector accordingly.  The test
		// for exactly zero is usually not recommend for floating point
		// arithmetic, but the implementation of Vector3::Normalize guarantees
		// the comparison is robust.  In this case, the A = pi and any axis
		// perpendicular to V1 may be used as the rotation axis.

		Vector3 kBisector = rkV1 + rkV2;
		kBisector.Normalize();

		double fCosHalfAngle = rkV1.Dot(kBisector);
		Vector3 kCross;

		m_afTuple[0] = fCosHalfAngle;

		if (fCosHalfAngle != 0.0)
		{
			kCross = kBisector.Cross(rkV1);//rkV1.Cross(kBisector);
			m_afTuple[1] = kCross.x;
			m_afTuple[2] = kCross.y;
			m_afTuple[3] = kCross.z;
		}
		else
		{
			double fInvLength;
			if (fabs(rkV1.x) >= fabs(rkV1.y))
			{
				// V1.x or V1.z is the largest magnitude component
				fInvLength = 1.0 / sqrt(rkV1.x*rkV1.x + rkV1.z*rkV1.z);
				m_afTuple[1] = -rkV1.z*fInvLength;
				m_afTuple[2] = 0.0;
				m_afTuple[3] = +rkV1.x*fInvLength;
			}
			else
			{
				// V1.y or V1.z is the largest magnitude component
				fInvLength = 1.0 / sqrt(rkV1.y*rkV1.y + rkV1.z*rkV1.z);
				m_afTuple[1] = 0.0;
				m_afTuple[2] = +rkV1.z*fInvLength;
				m_afTuple[3] = -rkV1.y*fInvLength;
			}
		}

		return *this;
	}

	inline Quaternion& Quaternion::Slerp (double fT, const Quaternion& rkP, const Quaternion& rkQ)
	{
		double fCos = Math::Clamp( rkP.Dot(rkQ), -1.0, 1.0 );

		Quaternion B = rkQ;

		if ( fCos < 0.0 )
		{
			fCos = -fCos;
			B = -B;
		}

		double fAngle = (double)acos(fCos);

		double fCoeff0;
		double fCoeff1;
		if ( (double)fabs(fAngle) >= (double)0.00001 )
		{
			double fSin = (double)sin(fAngle);
			double fInvSin = 1.0/fSin;
			fCoeff0 = (double)sin((1.0-fT)*fAngle)*fInvSin;
			fCoeff1 = (double)sin(fT*fAngle)*fInvSin;
		}
		else
		{
			fCoeff0 = 1.0 - fT;
			fCoeff1 = fT;
		}

		*this = rkP * fCoeff0 + B * fCoeff1;
		return *this;
	}

	inline Quaternion& Quaternion::SlerpExtraSpins(double fT,
	        const Quaternion& rkP, const Quaternion& rkQ, int iExtraSpins)
	{
		double fCos = Math::Clamp( rkP.Dot(rkQ), -1.0, 1.0 );
		double fAngle = acos(fCos);

		if (fabs(fAngle) >= (double)0.00001)
		{
			double fSin = sin(fAngle);
			double fPhase = Math::PI*iExtraSpins*fT;
			double fInvSin = 1.0/fSin;
			double fCoeff0 = sin((1.0-fT)*fAngle-fPhase)*fInvSin;
			double fCoeff1 = sin(fT*fAngle + fPhase)*fInvSin;
			*this = fCoeff0*rkP + fCoeff1*rkQ;
		}
		else
		{
			*this = rkP;
		}

		return *this;
	}
}

#endif
