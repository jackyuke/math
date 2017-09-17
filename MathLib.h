#ifndef _MATHLIB_H_
#define _MATHLIB_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>

#define EPISILON 1e-5

#define EPSILON_VALUE FLT_EPSILON

// ==
#define FLOAT_EQUAL(x, y, tol) (fabs((x) - (y)) < (tol))
// >
#define FLOAT_GT(x, y, tol) ((x) - (y) > (tol))
// >=
#define FLOAT_GE(x, y, tol) (FLOAT_GT(x, y, tol) || FLOAT_EQUAL(x, y, tol))
// <
#define FLOAT_LT(x, y, tol) ((x) < (y) - (tol))
// <=
#define FLOAT_LE(x, y, tol) (FLOAT_LT(x, y, tol) || FLOAT_EQUAL(x, y, tol))

class Math
{
public:
    static double GaussianSample ( double x, double mean, double std )
    {
        static const double sqrt2PI = sqrt ( 2.0 * PI );

        double dist = x - mean;
        return exp ( - ( dist * dist ) / ( 2.0 * std * std ) ) / ( std *  sqrt2PI );
    }

    static double GaussianSampleUnit ( double x, double mean, double std )
    {
        double dist = x - mean;
        return exp ( - ( dist * dist ) / ( 2.0 * std * std ) );
    }

    static double Clamp ( double x, const double a, const double b )
    {
        if ( x <= a )
        {
            return a;
        }
        if ( x >= b )
        {
            return b;
        }

        return x;
    }

    static double LinearClamp ( double x, const double a, const double b )
    {
        if ( x <= a )
        {
            return 0.0;
        }
        if ( x >= b )
        {
            return 1.0;
        }

        return ( x - a ) / ( b - a );
    }

    static double EaseClamp ( double x, const double a, const double b )
    {
        if ( x <= a )
        {
            return 0.0;
        }
        if ( x >= b )
        {
            return 1.0;
        }

        double t = ( x - a ) / ( b - a );
        return 0.5 * ( sin ( ( t - 0.5 ) * PI ) + 1.0 );
    }

    static double LinearInterp ( double t, const double a, const double b )
    {
        return a * ( 1.0 - t ) + b * t;
    }

    static double Min ( const double a, const double b )
    {
        return ( a < b ) ? a : b;
    }

    static double Max ( const double a, const double b )
    {
        return ( a > b ) ? a : b;
    }

    static const double PI;
    static const double DEG_TO_RAD;
};

class Box2D
{
public:
    Box2D();
    Box2D ( double width, double height );
    ~Box2D();

    inline double width()
    {
        return _right - _left;
    }

    inline double height()
    {
        return _bottom - _top;
    }

    inline double left()
    {
        return _left;
    }

    inline double top()
    {
        return _top;
    }

    inline double bottom()
    {
        return _bottom;
    }

    inline double right()
    {
        return _right;
    }

    inline void setWidth ( const double newWidth )
    {
        _right = _left + newWidth;
    }

    inline void setHeight ( const double newHeight )
    {
        _bottom = _top + newHeight;
    }

private:
    double _left, _top, _bottom, _right;
};

class Vector3
{
public:
    Vector3()
    {
        x = y = z = 0.0;
    }
    Vector3 ( double _c )
    {
        x = y = z = _c;
    }
    Vector3 ( double _x, double _y, double _z )
    {
        x = _x;
        y = _y;
        z = _z;
    }
    Vector3 ( const double xyz[3] )
    {
        x = xyz[0];
        y = xyz[1];
        z = xyz[2];
    }
    Vector3 ( const float xyz[3] )
    {
        x = ( double ) xyz[0];
        y = ( double ) xyz[1];
        z = ( double ) xyz[2];
    }
    Vector3 ( const Vector3& _rhs )
    {
        x = _rhs.x;
        y = _rhs.y;
        z = _rhs.z;
    }
    ~Vector3() {}

    const Vector3& operator= ( const Vector3& _rhs )
        {
            x = _rhs.x;
            y = _rhs.y;
            z = _rhs.z;
            return *this;
        }
    inline Vector3 operator-() const
    {
        Vector3 _v;
        _v.x = -x;
        _v.y = -y;
        _v.z = -z;
        return _v;
    }
    inline Vector3 operator- ( const Vector3& _rhs ) const
    {
        Vector3 _v;
        _v.x = x - _rhs.x;
        _v.y = y - _rhs.y;
        _v.z = z - _rhs.z;
        return _v;
    }
    inline Vector3 operator+ ( const Vector3& _rhs ) const
    {
        return Vector3 ( x + _rhs.x, y + _rhs.y, z + _rhs.z );
    }
    inline const Vector3& operator+= ( const Vector3& _rhs )
    {
        x += _rhs.x;
        y += _rhs.y;
        z += _rhs.z;
        return *this;
    }
    inline Vector3 operator* ( const double _rhs ) const
    {
        Vector3 _v;
        _v.x = x * _rhs;
        _v.y = y * _rhs;
        _v.z = z * _rhs;
        return _v;
    }
    inline Vector3 operator/ ( const double _rhs ) const
    {
        Vector3 _v;
        _v.x = x / _rhs;
        _v.y = y / _rhs;
        _v.z = z / _rhs;
        return _v;
    }
    inline double& operator[] ( const int i )
    {
        return ( ( double* ) this ) [i];
    }
    inline const double& operator[] ( const int i ) const
    {
        return ( ( double* ) this ) [i];
    }
    inline void ToArray ( double array[3] ) const
    {
        array[0] = x;
        array[1] = y;
        array[2] = z;
    }
    inline void ToArray ( float array[3] ) const
    {
        array[0] = ( float ) x;
        array[1] = ( float ) y;
        array[2] = ( float ) z;
    }
    inline double Length() const
    {
        return sqrt ( x*x + y*y + z*z );
    }
    inline double LengthSquared() const
    {
        return ( x*x + y*y + z*z );
    }
    inline const Vector3& Normalize()
    {
        double len = Length();
        if ( len > EPSILON_VALUE )
        {
            x = x / len;
            y = y / len;
            z = z / len;
        }
        else
        {
            x = y = z = 0.0;
        }
        return *this;
    }

    inline Vector3 Cross ( const Vector3& _b )	const
    {
        return Vector3 ( y*_b.z - z*_b.y, z*_b.x - x*_b.z, x*_b.y - y*_b.x );
    }
    inline double Dot ( const Vector3& _b ) const
    {
        return x*_b.x + y*_b.y + z*_b.z;
    }
    inline bool operator== ( const Vector3& _rhs ) const
    {
        return FLOAT_EQUAL(x, _rhs.x, EPISILON) && FLOAT_EQUAL(y, _rhs.y, EPISILON) && FLOAT_EQUAL(z, _rhs.z, EPISILON);
    }
    inline bool operator!= ( const Vector3& _rhs ) const
    {
        return !FLOAT_EQUAL(x, _rhs.x, EPISILON) || !FLOAT_EQUAL(y, _rhs.y, EPISILON) || !FLOAT_EQUAL(z, _rhs.z, EPISILON);
    }

    static const Vector3 Zero;

public:
    double x, y, z;
};

typedef class Vector3 Point3;
typedef std::vector<Point3> PointList;

typedef enum
{
    eInFront,
    eCoinciding,
    eInBack
}PointPosition;

class Plane
{
public:
    Plane()
        : m_normal()
        , m_position()
    {}
    Plane(const Vector3 &normal, const Vector3 &position)
        : m_normal(normal)
        , m_position(position)
    {}

    Plane(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3)
    {}

    Vector3 getNormal() const {return m_normal;}
    Vector3 getPosition() const {return m_position;}
    void setNormal(const Vector3 &normal) {m_normal = normal;}
    void setPosition(const Vector3 &position) {m_position = position;}

    double distance(const Vector3 &pos) const {return fabs((m_position - pos).Dot(m_normal));}
private:
    Vector3 m_normal, m_position;
};

class Line
{
public:
    Line (const Vector3 &p1, const Vector3 &p2 );
    ~Line ( void );

    bool calIntersectPointByPlane ( Plane* plane, Vector3& p );

private:
    Vector3 _p1, _p2;
};

/* Axis aligned bounding box in 3d space */
class AABB3D
{
public:
    AABB3D();
    inline AABB3D ( const Vector3& minPoint, const Vector3& maxPoint )
    {
        _min = minPoint;
        _max = maxPoint;
    }
    inline AABB3D ( double* minPoint, double* maxPoint )
    {
        _min = Vector3 ( minPoint );
        _max = Vector3 ( maxPoint );
    }
    AABB3D ( const AABB3D& rhs );
    ~AABB3D();

    inline
        void setMin ( Vector3 v )
    {
        _min = v;
    }

    inline
        void setMax ( Vector3 v )
    {
        _max = v;
    }

    inline
        Vector3 getMin() const
    {
        return _min;
    }

    inline
        Vector3 getMax() const
    {
        return _max;
    }

    AABB3D& operator= ( const AABB3D& rhs );

    Plane getSplitPlane();

    void combineAABB( const AABB3D& other );

private:
    Vector3 _min;
    Vector3 _max;
};

// represent a coordinate system
// x1 x2 x3
// y1 y2 y3
// z1 z2 z3
class Matrix3
{
public:
    Matrix3();
    Matrix3 ( double m[] );
    Matrix3 ( double m1, double m2, double m3,
              double m4, double m5, double m6,
              double m7, double m8, double m9 );
    Matrix3 ( const Matrix3& rhs )
    {
        copy ( rhs );
    }
    inline Matrix3& operator= ( const Matrix3& rhs )
        {
            if ( this != &rhs )
                copy ( rhs );
            return *this;
        }

    void setX ( Vector3& x );
    void setY ( Vector3& y );
    void setZ ( Vector3& z );
    Vector3 getX() const;
    Vector3 getY() const;
    Vector3 getZ() const;

    static const Matrix3 NORMALXYZ;

private:
    inline void copy ( const Matrix3& rhs )
    {
        _m[0] = rhs._m[0];
        _m[1] = rhs._m[1];
        _m[2] = rhs._m[2];
        _m[3] = rhs._m[3];
        _m[4] = rhs._m[4];
        _m[5] = rhs._m[5];
        _m[6] = rhs._m[6];
        _m[7] = rhs._m[7];
        _m[8] = rhs._m[8];
    }

private:
    double _m[9];
};

/*4 * 4 column matster  matrix*/
class Matrix4
{
public:
    Matrix4();
    Matrix4 ( double m1, double m2, double m3, double m4,
              double m5, double m6, double m7, double m8,
              double m9, double m10, double m11, double m12,
              double m13, double m14, double m15, double m16 )
    {
        _m[0] = m1;
        _m[1] = m2;
        _m[2] = m3;
        _m[3] = m4;
        _m[4] = m5;
        _m[5] = m6;
        _m[6] = m7;
        _m[7] = m8;
        _m[8] = m9;
        _m[9] = m10;
        _m[10] = m11;
        _m[11] = m12;
        _m[12] = m13;
        _m[13] = m14;
        _m[14] = m15;
        _m[15] = m16;
    }

    Matrix4 ( double* m )
    {
        for ( int i = 0; i < 16; i++ )
            _m[i] = m[i];
    }

    Matrix4 ( const Matrix4& rhs )
    {
        _copy ( rhs );
    }

    ~Matrix4();

    inline
        void getData ( float* m ) const
    {
        for ( int i = 0; i < 16; i++ )
            m[i] = _m[i];
    }

    inline
        const double* data()
    {
        return _m;
    }

    inline
        Matrix4 operator* ( const Matrix4& mr )
    {
        return
            Matrix4 (
                _m[0] * mr._m[0] + _m[1] * mr._m[4] + _m[2] * mr._m[8] + _m[3] * mr._m[12],
                _m[0] * mr._m[1] + _m[1] * mr._m[5] + _m[2] * mr._m[9] + _m[3] * mr._m[13],
                _m[0] * mr._m[2] + _m[1] * mr._m[6] + _m[2] * mr._m[10] + _m[3] * mr._m[14],
                _m[0] * mr._m[3] + _m[1] * mr._m[7] + _m[2] * mr._m[11] + _m[3] * mr._m[15],
                _m[4] * mr._m[0] + _m[5] * mr._m[4] + _m[6] * mr._m[8] + _m[7] * mr._m[12],
                _m[4] * mr._m[1] + _m[5] * mr._m[5] + _m[6] * mr._m[9] + _m[7] * mr._m[13],
                _m[4] * mr._m[2] + _m[5] * mr._m[6] + _m[6] * mr._m[10] + _m[7] * mr._m[14],
                _m[4] * mr._m[3] + _m[5] * mr._m[7] + _m[6] * mr._m[11] + _m[7] * mr._m[15],
                _m[8] * mr._m[0] + _m[9] * mr._m[4] + _m[10] * mr._m[8] + _m[11] * mr._m[12],
                _m[8] * mr._m[1] + _m[9] * mr._m[5] + _m[10] * mr._m[9] + _m[11] * mr._m[13],
                _m[8] * mr._m[2] + _m[9] * mr._m[6] + _m[10] * mr._m[10] + _m[11] * mr._m[14],
                _m[8] * mr._m[3] + _m[9] * mr._m[7] + _m[10] * mr._m[11] + _m[11] * mr._m[15],
                _m[12] * mr._m[0] + _m[13] * mr._m[4] + _m[14] * mr._m[8] + _m[15] * mr._m[12],
                _m[12] * mr._m[1] + _m[13] * mr._m[5] + _m[14] * mr._m[9] + _m[15] * mr._m[13],
                _m[12] * mr._m[2] + _m[13] * mr._m[6] + _m[14] * mr._m[10] + _m[15] * mr._m[14],
                _m[12] * mr._m[3] + _m[13] * mr._m[7] + _m[14] * mr._m[11] + _m[15] * mr._m[15] );
    }

    Vector3 rMul ( const Vector3& v ) const;
    Matrix3 rMul ( const Matrix3& m ) const;

    Matrix4& operator= ( const Matrix4& rhs )
        {
            _copy ( rhs );
            return *this;
        }

    inline
        bool operator== ( const Matrix4& rhs )
    {
        return FLOAT_EQUAL ( _m[0], rhs._m[0], EPISILON )
        && FLOAT_EQUAL ( _m[1], rhs._m[1], EPISILON )
        && FLOAT_EQUAL ( _m[2], rhs._m[2], EPISILON )
        && FLOAT_EQUAL ( _m[3], rhs._m[3], EPISILON )
        && FLOAT_EQUAL ( _m[4], rhs._m[4], EPISILON )
        && FLOAT_EQUAL ( _m[5], rhs._m[5], EPISILON )
        && FLOAT_EQUAL ( _m[6], rhs._m[6], EPISILON )
        && FLOAT_EQUAL ( _m[7], rhs._m[7], EPISILON )
        && FLOAT_EQUAL ( _m[8], rhs._m[8], EPISILON )
        && FLOAT_EQUAL ( _m[9], rhs._m[9], EPISILON )
        && FLOAT_EQUAL ( _m[10], rhs._m[10], EPISILON )
        && FLOAT_EQUAL ( _m[11], rhs._m[11], EPISILON )
        && FLOAT_EQUAL ( _m[12], rhs._m[12], EPISILON )
        && FLOAT_EQUAL ( _m[13], rhs._m[13], EPISILON )
        && FLOAT_EQUAL ( _m[14], rhs._m[14], EPISILON )
        && FLOAT_EQUAL ( _m[15], rhs._m[15], EPISILON );
    }

    inline
        Matrix4 transpose()
    {
        return Matrix4 ( _m[0], _m[4], _m[8], _m[12],
                         _m[1], _m[5], _m[9], _m[13],
                         _m[2], _m[6], _m[10], _m[14],
                         _m[3], _m[7], _m[11], _m[15] );
    }

    Matrix4 inverse();
    inline double operator[] ( int i )
    {
        return _m[i];
    }

    static const Matrix4 ZERO;
    static const Matrix4 IDENTITY;

private:
    inline void _copy ( const Matrix4& rhs )
    {
        for ( int i = 0; i < 16; i++ )
            _m[i] = rhs._m[i];
    }

private:
    double _m[16];
};

inline Vector3 operator* ( const Vector3& v, const Matrix4& m )
{
    return m.rMul ( v );
}

inline Matrix3 operator* ( const Matrix3& m, const Matrix4& t )
{
    return t.rMul ( m );
}

class Triangle
{
public:
    Triangle();
    ~Triangle();

    const Point3& operator[] ( int index ) const
    {
        return m_vertices[index];
    }

    Point3& operator[] ( int index )
    {
        return m_vertices[index];
    }

    const Point3& getNormal() const
    {
        return m_normal;
    }

private:
    Point3 m_vertices[3];
    Vector3 m_normal;
};

typedef class Triangle TriangleMesh;

typedef struct _HitTestOption
{
    bool HitPos;
    bool HitNeg;
}HitTestOption;

class Ray
{
public:
    Ray();
    ~Ray();

    bool hitTest ( const Triangle& mesh, HitTestOption* opt, double& dist );

private:
    Point3 m_start;
    Vector3 m_direction;
};

#endif
