#include "MathLib.h"

const double Math::PI = 3.1415926;
const double Math::DEG_TO_RAD = PI / 180.0;

/*
 * AABB3D
 */
AABB3D::AABB3D()
        : _min ( DBL_MAX, DBL_MAX, DBL_MAX )
        , _max ( -DBL_MAX, -DBL_MAX, -DBL_MAX )
{
}

AABB3D::~AABB3D()
{
}

AABB3D::AABB3D ( const AABB3D& rhs )
{
    operator= ( rhs );
}

AABB3D& AABB3D::operator= ( const AABB3D & rhs )
{
    if ( &rhs == this )
        return *this;

    _min = rhs._min;
    _max = rhs._max;
    return *this;
}

Plane AABB3D::getSplitPlane()
{
    //1.find the longest edge
    int max_index = 0;
    if ( ( _max[max_index] - _min[max_index] ) < ( _max[1] - _min[1] ) )
        max_index = 1;
    if ( ( _max[max_index] - _min[max_index] ) < ( _max[2] - _min[2] ) )
        max_index = 2;
    //2.split this edge
    double p1[3], p2[3];
    _min.ToArray ( p1 );
    _max.ToArray ( p2 );
    p1[max_index] = ( p1[max_index] + p2[max_index] ) / 2.0;
    p2[max_index] = p1[max_index];
    return Plane ( Vector3 ( p1 ),
                   Vector3 ( p1[0], p1[1], p2[2] ),
                   Vector3 ( p2 ) );
}

void AABB3D::combineAABB ( const AABB3D& other )
{
    for ( int i = 0; i < 3; i++ )
    {
        if ( FLOAT_GT ( _min[i], other._min[i], EPSILON_VALUE ) )
            _min[i] = other._min[i];
        if ( FLOAT_GT ( other._max[i], _max[i], EPSILON_VALUE ) )
            _max[i] = other._max[i];
    }
}

/*
 * Box2D
 */
Box2D::Box2D()
        : _left ( 0.0 )
        , _top ( 0.0 )
        , _bottom ( 0.0 )
        , _right ( 0.0 )
{
}

Box2D::Box2D ( double width, double height )
        : _left ( 0.0 )
        , _top ( 0.0 )
        , _bottom ( height )
        , _right ( width )
{
}

Box2D::~Box2D()
{
}

/*
 * Line
 */
Line::Line (const Vector3 &p1, const Vector3 &p2 )
    : _p1(p1)
    , _p2(p2)
{
}

Line::~Line ( void )
{
}

bool Line::calIntersectPointByPlane ( Plane* plane, Vector3& p )
{
    Vector3    u = _p2 - _p1;
    Vector3    w = _p1 - ( plane->getPosition() );
    Vector3 pNormalofPlane = plane->getNormal();

    double     D = pNormalofPlane.Dot ( u );
    double     N = -pNormalofPlane.Dot ( w );

    if ( fabs ( D ) < 0.000001 )
    {
        return false;
    }

    double t = N / D;

    p = _p1 + ( _p2 - _p1 ) * t;

    return true;
}

/*
 * Matrix3
 */
Matrix3::Matrix3()
{
    for ( int i = 0; i < 9; i++ )
        _m[i] = 0.0;
}

Matrix3::Matrix3 ( double m[] )
{
    for ( int i = 0; i < 9; i++ )
        _m[i] = m[i];
}

Matrix3::Matrix3 ( double m1, double m2, double m3,
                   double m4, double m5, double m6,
                   double m7, double m8, double m9 )
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
}

void Matrix3::setX ( Vector3& x )
{
    x.ToArray ( _m );
}

void Matrix3::setY ( Vector3& y )
{
    y.ToArray ( &_m[3] );
}

void Matrix3::setZ ( Vector3& z )
{
    z.ToArray ( &_m[6] );
}

Vector3 Matrix3::getX() const
{
    return Vector3 ( _m );
}

Vector3 Matrix3::getY() const
{
    return Vector3 ( &_m[3] );
}

Vector3 Matrix3::getZ() const
{
    return Vector3 ( &_m[6] );
}

const Matrix3 Matrix3::NORMALXYZ ( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 );

/*
 * Matrix4
 */
const Matrix4 Matrix4::ZERO (
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f );

const Matrix4 Matrix4::IDENTITY (
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f );



Matrix4::Matrix4()
{
    for ( int i = 0; i < 16; i++ )
        _m[i] = 0.0f;
}

Matrix4::~Matrix4()
{
}

Vector3 Matrix4::rMul ( const Vector3& v ) const
{
    return Vector3 ( _m[0] * v[0] + _m[1] * v[1] + _m[2] * v[2] + _m[3],
                     _m[4] * v[0] + _m[5] * v[1] + _m[6] * v[2] + _m[7],
                     _m[8] * v[0] + _m[9] * v[1] + _m[10] * v[2] + _m[11] );
}

Matrix3 Matrix4::rMul ( const Matrix3& m ) const
{
    Vector3 x = m.getX();
    x = x * ( *this );
    Vector3 y = m.getY();
    y = y * ( *this );
    Vector3 z = m.getZ();
    z = z * ( *this );
    Matrix3 r;
    r.setX ( x );
    r.setY ( y );
    r.setZ ( z );
    return r;
}

Matrix4 Matrix4::inverse()
{
    float m00 = _m[0], m01 = _m[1], m02 = _m[2], m03 = _m[3];
    float m10 = _m[4], m11 = _m[5], m12 = _m[6], m13 = _m[7];
    float m20 = _m[8], m21 = _m[9], m22 = _m[10], m23 = _m[11];
    float m30 = _m[12], m31 = _m[13], m32 = _m[14], m33 = _m[15];

    float v0 = m20 * m31 - m21 * m30;
    float v1 = m20 * m32 - m22 * m30;
    float v2 = m20 * m33 - m23 * m30;
    float v3 = m21 * m32 - m22 * m31;
    float v4 = m21 * m33 - m23 * m31;
    float v5 = m22 * m33 - m23 * m32;

    float t00 = + ( v5 * m11 - v4 * m12 + v3 * m13 );
    float t10 = - ( v5 * m10 - v2 * m12 + v1 * m13 );
    float t20 = + ( v4 * m10 - v2 * m11 + v0 * m13 );
    float t30 = - ( v3 * m10 - v1 * m11 + v0 * m12 );

    float invDet = 1 / ( t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03 );

    float d00 = t00 * invDet;
    float d10 = t10 * invDet;
    float d20 = t20 * invDet;
    float d30 = t30 * invDet;

    float d01 = - ( v5 * m01 - v4 * m02 + v3 * m03 ) * invDet;
    float d11 = + ( v5 * m00 - v2 * m02 + v1 * m03 ) * invDet;
    float d21 = - ( v4 * m00 - v2 * m01 + v0 * m03 ) * invDet;
    float d31 = + ( v3 * m00 - v1 * m01 + v0 * m02 ) * invDet;

    v0 = m10 * m31 - m11 * m30;
    v1 = m10 * m32 - m12 * m30;
    v2 = m10 * m33 - m13 * m30;
    v3 = m11 * m32 - m12 * m31;
    v4 = m11 * m33 - m13 * m31;
    v5 = m12 * m33 - m13 * m32;

    float d02 = + ( v5 * m01 - v4 * m02 + v3 * m03 ) * invDet;
    float d12 = - ( v5 * m00 - v2 * m02 + v1 * m03 ) * invDet;
    float d22 = + ( v4 * m00 - v2 * m01 + v0 * m03 ) * invDet;
    float d32 = - ( v3 * m00 - v1 * m01 + v0 * m02 ) * invDet;

    v0 = m21 * m10 - m20 * m11;
    v1 = m22 * m10 - m20 * m12;
    v2 = m23 * m10 - m20 * m13;
    v3 = m22 * m11 - m21 * m12;
    v4 = m23 * m11 - m21 * m13;
    v5 = m23 * m12 - m22 * m13;

    float d03 = - ( v5 * m01 - v4 * m02 + v3 * m03 ) * invDet;
    float d13 = + ( v5 * m00 - v2 * m02 + v1 * m03 ) * invDet;
    float d23 = - ( v4 * m00 - v2 * m01 + v0 * m03 ) * invDet;
    float d33 = + ( v3 * m00 - v1 * m01 + v0 * m02 ) * invDet;

    return Matrix4 ( d00, d01, d02, d03,
                     d10, d11, d12, d13,
                     d20, d21, d22, d23,
                     d30, d31, d32, d33 );
}

Ray::Ray()
{
}

Ray::~Ray()
{
}

const Vector3 Vector3::Zero ( 0, 0, 0 );

const double EPSILON = 1e-4f;
inline double fabs ( double a )
{
    return a < 0.0 ? -a : a;
}

bool Ray::hitTest ( const Triangle& mesh, HitTestOption* opt, double& dist )
{
    double t;
    {
        double denom = mesh.getNormal().Dot ( m_direction );
        //检查相交平面
        if ( denom > EPSILON )
        {
            if ( !opt->HitNeg )
            {
                dist = 0;
                return false;
            }
        }
        else if ( denom < -EPSILON )
        {
            if ( !opt->HitPos )
            {
                dist = 0;
                return false;
            }
        }
        else
        {
            dist = 0;
            return false;
        }

        t = mesh.getNormal().Dot ( mesh[0] - m_start );
        if ( t < 0 )
        {
            dist = 0;
            return false;
        }
    }

    size_t i0, i1;
    {
        double n0 = fabs ( mesh.getNormal() [0] );
        double n1 = fabs ( mesh.getNormal() [1] );
        double n2 = fabs ( mesh.getNormal() [2] );
        i0 = 1;
        i1 = 2;
        if ( n1 > n2 )
        {
            if ( n1 > n0 ) i0 = 0;
        }
        else
        {
            if ( n2 > n0 ) i1 = 0;
        }
    }

    //检查交点是否在三角形内
    double u1 = mesh[1][i0] - mesh[0][i0];
    double v1 = mesh[1][i1] - mesh[0][i1];
    double u2 = mesh[2][i0] - mesh[0][i0];
    double v2 = mesh[2][i1] - mesh[0][i1];
    double u0 = t * m_direction[i0] + m_start[i0] - mesh[0][i0];
    double v0 = t * m_direction[i1] + m_start[i1] - mesh[0][i1];
    double alpha = u0 * v2 - u2 * v0;
    double beta = u1 * v0 - u0 * v1;
    double area = u1 * v2 - u2 * v1;
    double tolerance = -EPSILON * area;

    if ( area > 0 )
    {
        if ( alpha < tolerance || beta < tolerance || alpha + beta > area - tolerance )
        {
            dist = 0;
            return false;
        }
    }
    else
    {
        if ( alpha > tolerance || beta > tolerance || alpha + beta < area - tolerance )
        {
            dist = 0;
            return false;
        }
    }

    dist = t;
    return true;
}

Triangle::Triangle()
{

}

Triangle::~Triangle()
{

}
