#include "vecmath.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#define _USE_MATH_DEFINES
#define M_PI       3.14159265358979323846   // pi

Matrix2f::Matrix2f(float fill)
{
    for (int i = 0; i < 4; ++i)
    {
        m_elements[i] = fill;
    }
}

Matrix2f::Matrix2f(float m00, float m01,
                   float m10, float m11)
{
    m_elements[0] = m00;
    m_elements[1] = m10;

    m_elements[2] = m01;
    m_elements[3] = m11;
}

Matrix2f::Matrix2f(const Vector2f &v0, const Vector2f &v1, bool setColumns)
{
    if (setColumns)
    {
        setCol(0, v0);
        setCol(1, v1);
    }
    else
    {
        setRow(0, v0);
        setRow(1, v1);
    }
}

Matrix2f::Matrix2f(const Matrix2f &rm)
{
    memcpy(m_elements, rm.m_elements, 2 * sizeof(float));
}

Matrix2f &Matrix2f::operator=(const Matrix2f &rm)
{
    if (this != &rm)
    {
        memcpy(m_elements, rm.m_elements, 2 * sizeof(float));
    }
    return *this;
}

const float &Matrix2f::operator()(int i, int j) const
{
    return m_elements[j * 2 + i];
}

float &Matrix2f::operator()(int i, int j)
{
    return m_elements[j * 2 + i];
}

Vector2f Matrix2f::getRow(int i) const
{
    return Vector2f(
        m_elements[i],
        m_elements[i + 2]);
}

void Matrix2f::setRow(int i, const Vector2f &v)
{
    m_elements[i] = v.x();
    m_elements[i + 2] = v.y();
}

Vector2f Matrix2f::getCol(int j) const
{
    int colStart = 2 * j;

    return Vector2f(
        m_elements[colStart],
        m_elements[colStart + 1]);
}

void Matrix2f::setCol(int j, const Vector2f &v)
{
    int colStart = 2 * j;

    m_elements[colStart] = v.x();
    m_elements[colStart + 1] = v.y();
}

float Matrix2f::determinant()
{
    return Matrix2f::determinant2x2(
        m_elements[0], m_elements[2],
        m_elements[1], m_elements[3]);
}

Matrix2f Matrix2f::inverse(bool *pbIsSingular, float epsilon)
{
    float determinant = m_elements[0] * m_elements[3] - m_elements[2] * m_elements[1];

    bool isSingular = (fabs(determinant) < epsilon);
    if (isSingular)
    {
        if (pbIsSingular != NULL)
        {
            *pbIsSingular = true;
        }
        return Matrix2f();
    }
    else
    {
        if (pbIsSingular != NULL)
        {
            *pbIsSingular = false;
        }

        float reciprocalDeterminant = 1.0f / determinant;

        return Matrix2f(
            m_elements[3] * reciprocalDeterminant, -m_elements[2] * reciprocalDeterminant,
            -m_elements[1] * reciprocalDeterminant, m_elements[0] * reciprocalDeterminant);
    }
}

void Matrix2f::transpose()
{
    float m01 = (*this)(0, 1);
    float m10 = (*this)(1, 0);

    (*this)(0, 1) = m10;
    (*this)(1, 0) = m01;
}

Matrix2f Matrix2f::transposed() const
{
    return Matrix2f(
        (*this)(0, 0), (*this)(1, 0),
        (*this)(0, 1), (*this)(1, 1));
}

Matrix2f::operator float *()
{
    return m_elements;
}

void Matrix2f::print()
{
    printf("[ %.4f %.4f ]\n[ %.4f %.4f ]\n",
           m_elements[0], m_elements[2],
           m_elements[1], m_elements[3]);
}

// static
float Matrix2f::determinant2x2(float m00, float m01,
                               float m10, float m11)
{
    return (m00 * m11 - m01 * m10);
}

// static
Matrix2f Matrix2f::ones()
{
    Matrix2f m;
    for (int i = 0; i < 4; ++i)
    {
        m.m_elements[i] = 1;
    }

    return m;
}

// static
Matrix2f Matrix2f::identity()
{
    Matrix2f m;

    m(0, 0) = 1;
    m(1, 1) = 1;

    return m;
}

// static
Matrix2f Matrix2f::rotation(float degrees)
{
    float c = cos(degrees);
    float s = sin(degrees);

    return Matrix2f(
        c, -s,
        s, c);
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Matrix2f operator*(float f, const Matrix2f &m)
{
    Matrix2f output;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            output(i, j) = f * m(i, j);
        }
    }

    return output;
}

Matrix2f operator*(const Matrix2f &m, float f)
{
    return f * m;
}

Vector2f operator*(const Matrix2f &m, const Vector2f &v)
{
    Vector2f output(0, 0);

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            output[i] += m(i, j) * v[j];
        }
    }

    return output;
}

Matrix2f operator*(const Matrix2f &x, const Matrix2f &y)
{
    Matrix2f product; // zeroes

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                product(i, k) += x(i, j) * y(j, k);
            }
        }
    }

    return product;
}

Matrix3f::Matrix3f(float fill)
{
    for (int i = 0; i < 9; ++i)
    {
        m_elements[i] = fill;
    }
}

Matrix3f::Matrix3f(float m00, float m01, float m02,
                   float m10, float m11, float m12,
                   float m20, float m21, float m22)
{
    m_elements[0] = m00;
    m_elements[1] = m10;
    m_elements[2] = m20;

    m_elements[3] = m01;
    m_elements[4] = m11;
    m_elements[5] = m21;

    m_elements[6] = m02;
    m_elements[7] = m12;
    m_elements[8] = m22;
}

Matrix3f::Matrix3f(const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, bool setColumns)
{
    if (setColumns)
    {
        setCol(0, v0);
        setCol(1, v1);
        setCol(2, v2);
    }
    else
    {
        setRow(0, v0);
        setRow(1, v1);
        setRow(2, v2);
    }
}

Matrix3f::Matrix3f(const Matrix3f &rm)
{
    memcpy(m_elements, rm.m_elements, 9 * sizeof(float));
}

Matrix3f &Matrix3f::operator=(const Matrix3f &rm)
{
    if (this != &rm)
    {
        memcpy(m_elements, rm.m_elements, 9 * sizeof(float));
    }
    return *this;
}

const float &Matrix3f::operator()(int i, int j) const
{
    return m_elements[j * 3 + i];
}

float &Matrix3f::operator()(int i, int j)
{
    return m_elements[j * 3 + i];
}

Vector3f Matrix3f::getRow(int i) const
{
    return Vector3f(
        m_elements[i],
        m_elements[i + 3],
        m_elements[i + 6]);
}

void Matrix3f::setRow(int i, const Vector3f &v)
{
    m_elements[i] = v.x();
    m_elements[i + 3] = v.y();
    m_elements[i + 6] = v.z();
}

Vector3f Matrix3f::getCol(int j) const
{
    int colStart = 3 * j;

    return Vector3f(
        m_elements[colStart],
        m_elements[colStart + 1],
        m_elements[colStart + 2]);
}

void Matrix3f::setCol(int j, const Vector3f &v)
{
    int colStart = 3 * j;

    m_elements[colStart] = v.x();
    m_elements[colStart + 1] = v.y();
    m_elements[colStart + 2] = v.z();
}

Matrix2f Matrix3f::getSubmatrix2x2(int i0, int j0) const
{
    Matrix2f out;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            out(i, j) = (*this)(i + i0, j + j0);
        }
    }

    return out;
}

void Matrix3f::setSubmatrix2x2(int i0, int j0, const Matrix2f &m)
{
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            (*this)(i + i0, j + j0) = m(i, j);
        }
    }
}

float Matrix3f::determinant() const
{
    return Matrix3f::determinant3x3(
        m_elements[0], m_elements[3], m_elements[6],
        m_elements[1], m_elements[4], m_elements[7],
        m_elements[2], m_elements[5], m_elements[8]);
}

Matrix3f Matrix3f::inverse(bool *pbIsSingular, float epsilon) const
{
    float m00 = m_elements[0];
    float m10 = m_elements[1];
    float m20 = m_elements[2];

    float m01 = m_elements[3];
    float m11 = m_elements[4];
    float m21 = m_elements[5];

    float m02 = m_elements[6];
    float m12 = m_elements[7];
    float m22 = m_elements[8];

    float cofactor00 = Matrix2f::determinant2x2(m11, m12, m21, m22);
    float cofactor01 = -Matrix2f::determinant2x2(m10, m12, m20, m22);
    float cofactor02 = Matrix2f::determinant2x2(m10, m11, m20, m21);

    float cofactor10 = -Matrix2f::determinant2x2(m01, m02, m21, m22);
    float cofactor11 = Matrix2f::determinant2x2(m00, m02, m20, m22);
    float cofactor12 = -Matrix2f::determinant2x2(m00, m01, m20, m21);

    float cofactor20 = Matrix2f::determinant2x2(m01, m02, m11, m12);
    float cofactor21 = -Matrix2f::determinant2x2(m00, m02, m10, m12);
    float cofactor22 = Matrix2f::determinant2x2(m00, m01, m10, m11);

    float determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02;

    bool isSingular = (fabs(determinant) < epsilon);
    if (isSingular)
    {
        if (pbIsSingular != NULL)
        {
            *pbIsSingular = true;
        }
        return Matrix3f();
    }
    else
    {
        if (pbIsSingular != NULL)
        {
            *pbIsSingular = false;
        }

        float reciprocalDeterminant = 1.0f / determinant;

        return Matrix3f(
            cofactor00 * reciprocalDeterminant, cofactor10 * reciprocalDeterminant, cofactor20 * reciprocalDeterminant,
            cofactor01 * reciprocalDeterminant, cofactor11 * reciprocalDeterminant, cofactor21 * reciprocalDeterminant,
            cofactor02 * reciprocalDeterminant, cofactor12 * reciprocalDeterminant, cofactor22 * reciprocalDeterminant);
    }
}

void Matrix3f::transpose()
{
    float temp;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = i + 1; j < 3; ++j)
        {
            temp = (*this)(i, j);
            (*this)(i, j) = (*this)(j, i);
            (*this)(j, i) = temp;
        }
    }
}

Matrix3f Matrix3f::transposed() const
{
    Matrix3f out;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            out(j, i) = (*this)(i, j);
        }
    }

    return out;
}

Matrix3f::operator float *()
{
    return m_elements;
}

void Matrix3f::print()
{
    printf("[ %.4f %.4f %.4f ]\n[ %.4f %.4f %.4f ]\n[ %.4f %.4f %.4f ]\n",
           m_elements[0], m_elements[3], m_elements[6],
           m_elements[1], m_elements[4], m_elements[7],
           m_elements[2], m_elements[5], m_elements[8]);
}

// static
float Matrix3f::determinant3x3(float m00, float m01, float m02,
                               float m10, float m11, float m12,
                               float m20, float m21, float m22)
{
    return (
        m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20));
}

// static
Matrix3f Matrix3f::ones()
{
    Matrix3f m;
    for (int i = 0; i < 9; ++i)
    {
        m.m_elements[i] = 1;
    }

    return m;
}

// static
Matrix3f Matrix3f::identity()
{
    Matrix3f m;

    m(0, 0) = 1;
    m(1, 1) = 1;
    m(2, 2) = 1;

    return m;
}

// static
Matrix3f Matrix3f::rotateX(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(
        1, 0, 0,
        0, c, -s,
        0, s, c);
}

// static
Matrix3f Matrix3f::rotateY(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(
        c, 0, s,
        0, 1, 0,
        -s, 0, c);
}

// static
Matrix3f Matrix3f::rotateZ(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(
        c, -s, 0,
        s, c, 0,
        0, 0, 1);
}

// static
Matrix3f Matrix3f::scaling(float sx, float sy, float sz)
{
    return Matrix3f(
        sx, 0, 0,
        0, sy, 0,
        0, 0, sz);
}

// static
Matrix3f Matrix3f::uniformScaling(float s)
{
    return Matrix3f(
        s, 0, 0,
        0, s, 0,
        0, 0, s);
}

// static
Matrix3f Matrix3f::rotation(const Vector3f &rDirection, float radians)
{
    Vector3f normalizedDirection = rDirection.normalized();

    float cosTheta = cos(radians);
    float sinTheta = sin(radians);

    float x = normalizedDirection.x();
    float y = normalizedDirection.y();
    float z = normalizedDirection.z();

    return Matrix3f(
        x * x * (1.0f - cosTheta) + cosTheta, y * x * (1.0f - cosTheta) - z * sinTheta, z * x * (1.0f - cosTheta) + y * sinTheta,
        x * y * (1.0f - cosTheta) + z * sinTheta, y * y * (1.0f - cosTheta) + cosTheta, z * y * (1.0f - cosTheta) - x * sinTheta,
        x * z * (1.0f - cosTheta) - y * sinTheta, y * z * (1.0f - cosTheta) + x * sinTheta, z * z * (1.0f - cosTheta) + cosTheta);
}

// static
Matrix3f Matrix3f::rotation(const Quat4f &rq)
{
    Quat4f q = rq.normalized();

    float xx = q.x() * q.x();
    float yy = q.y() * q.y();
    float zz = q.z() * q.z();

    float xy = q.x() * q.y();
    float zw = q.z() * q.w();

    float xz = q.x() * q.z();
    float yw = q.y() * q.w();

    float yz = q.y() * q.z();
    float xw = q.x() * q.w();

    return Matrix3f(
        1.0f - 2.0f * (yy + zz), 2.0f * (xy - zw), 2.0f * (xz + yw),
        2.0f * (xy + zw), 1.0f - 2.0f * (xx + zz), 2.0f * (yz - xw),
        2.0f * (xz - yw), 2.0f * (yz + xw), 1.0f - 2.0f * (xx + yy));
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Vector3f operator*(const Matrix3f &m, const Vector3f &v)
{
    Vector3f output(0, 0, 0);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            output[i] += m(i, j) * v[j];
        }
    }

    return output;
}

Matrix3f operator*(const Matrix3f &x, const Matrix3f &y)
{
    Matrix3f product; // zeroes

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                product(i, k) += x(i, j) * y(j, k);
            }
        }
    }

    return product;
}

Matrix4f::Matrix4f( float fill )
{
	for( int i = 0; i < 16; ++i )
	{
		m_elements[ i ] = fill;
	}
}

Matrix4f::Matrix4f( float m00, float m01, float m02, float m03,
				   float m10, float m11, float m12, float m13,
				   float m20, float m21, float m22, float m23,
				   float m30, float m31, float m32, float m33 )
{
	m_elements[ 0 ] = m00;
	m_elements[ 1 ] = m10;
	m_elements[ 2 ] = m20;
	m_elements[ 3 ] = m30;
	
	m_elements[ 4 ] = m01;
	m_elements[ 5 ] = m11;
	m_elements[ 6 ] = m21;
	m_elements[ 7 ] = m31;

	m_elements[ 8 ] = m02;
	m_elements[ 9 ] = m12;
	m_elements[ 10 ] = m22;
	m_elements[ 11 ] = m32;

	m_elements[ 12 ] = m03;
	m_elements[ 13 ] = m13;
	m_elements[ 14 ] = m23;
	m_elements[ 15 ] = m33;
}

Matrix4f& Matrix4f::operator/=(float d)
{
	for(int ii=0;ii<16;ii++){
		m_elements[ii]/=d;
	}
	return *this;
}

Matrix4f::Matrix4f( const Vector4f& v0, const Vector4f& v1, const Vector4f& v2, const Vector4f& v3, bool setColumns )
{
	if( setColumns )
	{
		setCol( 0, v0 );
		setCol( 1, v1 );
		setCol( 2, v2 );
		setCol( 3, v3 );
	}
	else
	{
		setRow( 0, v0 );
		setRow( 1, v1 );
		setRow( 2, v2 );
		setRow( 3, v3 );
	}
}

Matrix4f::Matrix4f( const Matrix4f& rm )
{
	memcpy( m_elements, rm.m_elements, 16 * sizeof( float ) );
}

Matrix4f& Matrix4f::operator = ( const Matrix4f& rm )
{
	if( this != &rm )
	{
		memcpy( m_elements, rm.m_elements, 16 * sizeof( float ) );
	}
	return *this;
}

const float& Matrix4f::operator () ( int i, int j ) const
{
	return m_elements[ j * 4 + i ];
}

float& Matrix4f::operator () ( int i, int j )
{
	return m_elements[ j * 4 + i ];
}

Vector4f Matrix4f::getRow( int i ) const
{
	return Vector4f
	(
		m_elements[ i ],
		m_elements[ i + 4 ],
		m_elements[ i + 8 ],
		m_elements[ i + 12 ]
	);
}

void Matrix4f::setRow( int i, const Vector4f& v )
{
	m_elements[ i ] = v.x();
	m_elements[ i + 4 ] = v.y();
	m_elements[ i + 8 ] = v.z();
	m_elements[ i + 12 ] = v.w();
}

Vector4f Matrix4f::getCol( int j ) const
{
	int colStart = 4 * j;

	return Vector4f
	(
		m_elements[ colStart ],
		m_elements[ colStart + 1 ],
		m_elements[ colStart + 2 ],
		m_elements[ colStart + 3 ]
	);
}

void Matrix4f::setCol( int j, const Vector4f& v )
{
	int colStart = 4 * j;

	m_elements[ colStart ] = v.x();
	m_elements[ colStart + 1 ] = v.y();
	m_elements[ colStart + 2 ] = v.z();
	m_elements[ colStart + 3 ] = v.w();
}

Matrix2f Matrix4f::getSubmatrix2x2( int i0, int j0 ) const
{
	Matrix2f out;

	for( int i = 0; i < 2; ++i )
	{
		for( int j = 0; j < 2; ++j )
		{
			out( i, j ) = ( *this )( i + i0, j + j0 );
		}
	}

	return out;
}

Matrix3f Matrix4f::getSubmatrix3x3( int i0, int j0 ) const
{
	Matrix3f out;

	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
			out( i, j ) = ( *this )( i + i0, j + j0 );
		}
	}

	return out;
}

void Matrix4f::setSubmatrix2x2( int i0, int j0, const Matrix2f& m )
{
	for( int i = 0; i < 2; ++i )
	{
		for( int j = 0; j < 2; ++j )
		{
			( *this )( i + i0, j + j0 ) = m( i, j );
		}
	}
}

void Matrix4f::setSubmatrix3x3( int i0, int j0, const Matrix3f& m )
{
	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
			( *this )( i + i0, j + j0 ) = m( i, j );
		}
	}
}

float Matrix4f::determinant() const
{
	float m00 = m_elements[ 0 ];
	float m10 = m_elements[ 1 ];
	float m20 = m_elements[ 2 ];
	float m30 = m_elements[ 3 ];

	float m01 = m_elements[ 4 ];
	float m11 = m_elements[ 5 ];
	float m21 = m_elements[ 6 ];
	float m31 = m_elements[ 7 ];

	float m02 = m_elements[ 8 ];
	float m12 = m_elements[ 9 ];
	float m22 = m_elements[ 10 ];
	float m32 = m_elements[ 11 ];

	float m03 = m_elements[ 12 ];
	float m13 = m_elements[ 13 ];
	float m23 = m_elements[ 14 ];
	float m33 = m_elements[ 15 ];

	float cofactor00 =  Matrix3f::determinant3x3( m11, m12, m13, m21, m22, m23, m31, m32, m33 );
	float cofactor01 = -Matrix3f::determinant3x3( m12, m13, m10, m22, m23, m20, m32, m33, m30 );
	float cofactor02 =  Matrix3f::determinant3x3( m13, m10, m11, m23, m20, m21, m33, m30, m31 );
	float cofactor03 = -Matrix3f::determinant3x3( m10, m11, m12, m20, m21, m22, m30, m31, m32 );

	return( m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02 + m03 * cofactor03 );
}

Matrix4f Matrix4f::inverse( bool* pbIsSingular, float epsilon ) const
{
	float m00 = m_elements[ 0 ];
	float m10 = m_elements[ 1 ];
	float m20 = m_elements[ 2 ];
	float m30 = m_elements[ 3 ];

	float m01 = m_elements[ 4 ];
	float m11 = m_elements[ 5 ];
	float m21 = m_elements[ 6 ];
	float m31 = m_elements[ 7 ];

	float m02 = m_elements[ 8 ];
	float m12 = m_elements[ 9 ];
	float m22 = m_elements[ 10 ];
	float m32 = m_elements[ 11 ];

	float m03 = m_elements[ 12 ];
	float m13 = m_elements[ 13 ];
	float m23 = m_elements[ 14 ];
	float m33 = m_elements[ 15 ];

    float cofactor00 =  Matrix3f::determinant3x3( m11, m12, m13, m21, m22, m23, m31, m32, m33 );
    float cofactor01 = -Matrix3f::determinant3x3( m12, m13, m10, m22, m23, m20, m32, m33, m30 );
    float cofactor02 =  Matrix3f::determinant3x3( m13, m10, m11, m23, m20, m21, m33, m30, m31 );
    float cofactor03 = -Matrix3f::determinant3x3( m10, m11, m12, m20, m21, m22, m30, m31, m32 );
    
    float cofactor10 = -Matrix3f::determinant3x3( m21, m22, m23, m31, m32, m33, m01, m02, m03 );
    float cofactor11 =  Matrix3f::determinant3x3( m22, m23, m20, m32, m33, m30, m02, m03, m00 );
    float cofactor12 = -Matrix3f::determinant3x3( m23, m20, m21, m33, m30, m31, m03, m00, m01 );
    float cofactor13 =  Matrix3f::determinant3x3( m20, m21, m22, m30, m31, m32, m00, m01, m02 );
    
    float cofactor20 =  Matrix3f::determinant3x3( m31, m32, m33, m01, m02, m03, m11, m12, m13 );
    float cofactor21 = -Matrix3f::determinant3x3( m32, m33, m30, m02, m03, m00, m12, m13, m10 );
    float cofactor22 =  Matrix3f::determinant3x3( m33, m30, m31, m03, m00, m01, m13, m10, m11 );
    float cofactor23 = -Matrix3f::determinant3x3( m30, m31, m32, m00, m01, m02, m10, m11, m12 );
    
    float cofactor30 = -Matrix3f::determinant3x3( m01, m02, m03, m11, m12, m13, m21, m22, m23 );
    float cofactor31 =  Matrix3f::determinant3x3( m02, m03, m00, m12, m13, m10, m22, m23, m20 );
    float cofactor32 = -Matrix3f::determinant3x3( m03, m00, m01, m13, m10, m11, m23, m20, m21 );
    float cofactor33 =  Matrix3f::determinant3x3( m00, m01, m02, m10, m11, m12, m20, m21, m22 );

	float determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02 + m03 * cofactor03;

	bool isSingular = ( fabs( determinant ) < epsilon );
	if( isSingular )
	{
		if( pbIsSingular != NULL )
		{
			*pbIsSingular = true;
		}
		return Matrix4f();
	}
	else
	{
		if( pbIsSingular != NULL )
		{
			*pbIsSingular = false;
		}

		float reciprocalDeterminant = 1.0f / determinant;

		return Matrix4f
			(
				cofactor00 * reciprocalDeterminant, cofactor10 * reciprocalDeterminant, cofactor20 * reciprocalDeterminant, cofactor30 * reciprocalDeterminant,
				cofactor01 * reciprocalDeterminant, cofactor11 * reciprocalDeterminant, cofactor21 * reciprocalDeterminant, cofactor31 * reciprocalDeterminant,
				cofactor02 * reciprocalDeterminant, cofactor12 * reciprocalDeterminant, cofactor22 * reciprocalDeterminant, cofactor32 * reciprocalDeterminant,
				cofactor03 * reciprocalDeterminant, cofactor13 * reciprocalDeterminant, cofactor23 * reciprocalDeterminant, cofactor33 * reciprocalDeterminant
			);
	}
}

void Matrix4f::transpose()
{
	float temp;

	for( int i = 0; i < 3; ++i )
	{
		for( int j = i + 1; j < 4; ++j )
		{
			temp = ( *this )( i, j );
			( * this )( i, j ) = ( *this )( j, i );
			( *this )( j, i ) = temp;
		}
	}
}

Matrix4f Matrix4f::transposed() const
{
	Matrix4f out;
	for( int i = 0; i < 4; ++i )
	{
		for( int j = 0; j < 4; ++j )
		{
			out( j, i ) = ( *this )( i, j );
		}
	}

	return out;
}

Matrix4f::operator float* ()
{
	return m_elements;
}

Matrix4f::operator const float* ()const
{
	return m_elements;
}


void Matrix4f::print()
{
	printf( "[ %.4f %.4f %.4f %.4f ]\n[ %.4f %.4f %.4f %.4f ]\n[ %.4f %.4f %.4f %.4f ]\n[ %.4f %.4f %.4f %.4f ]\n",
		m_elements[ 0 ], m_elements[ 4 ], m_elements[ 8 ], m_elements[ 12 ],
		m_elements[ 1 ], m_elements[ 5 ], m_elements[ 9 ], m_elements[ 13 ],
		m_elements[ 2 ], m_elements[ 6 ], m_elements[ 10], m_elements[ 14 ],
		m_elements[ 3 ], m_elements[ 7 ], m_elements[ 11], m_elements[ 15 ] );
}

// static
Matrix4f Matrix4f::ones()
{
	Matrix4f m;
	for( int i = 0; i < 16; ++i )
	{
		m.m_elements[ i ] = 1;
	}

	return m;
}

// static
Matrix4f Matrix4f::identity()
{
	Matrix4f m;
	
	m( 0, 0 ) = 1;
	m( 1, 1 ) = 1;
	m( 2, 2 ) = 1;
	m( 3, 3 ) = 1;

	return m;
}

// static
Matrix4f Matrix4f::translation( float x, float y, float z )
{
	return Matrix4f
	(
		1, 0, 0, x,
		0, 1, 0, y,
		0, 0, 1, z,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::translation( const Vector3f& rTranslation )
{
	return Matrix4f
	(
		1, 0, 0, rTranslation.x(),
		0, 1, 0, rTranslation.y(),
		0, 0, 1, rTranslation.z(),
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::rotateX( float radians )
{
	float c = cos( radians );
	float s = sin( radians );

	return Matrix4f
	(
		1, 0, 0, 0,
		0, c, -s, 0,
		0, s, c, 0,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::rotateY( float radians )
{
	float c = cos( radians );
	float s = sin( radians );

	return Matrix4f
	(
		c, 0, s, 0,
		0, 1, 0, 0,
		-s, 0, c, 0,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::rotateZ( float radians )
{
	float c = cos( radians );
	float s = sin( radians );

	return Matrix4f
	(
		c, -s, 0, 0,
		s, c, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::rotation( const Vector3f& rDirection, float radians )
{
	Vector3f normalizedDirection = rDirection.normalized();
	
	float cosTheta = cos( radians );
	float sinTheta = sin( radians );

	float x = normalizedDirection.x();
	float y = normalizedDirection.y();
	float z = normalizedDirection.z();

	return Matrix4f
	(
		x * x * ( 1.0f - cosTheta ) + cosTheta,			y * x * ( 1.0f - cosTheta ) - z * sinTheta,		z * x * ( 1.0f - cosTheta ) + y * sinTheta,		0.0f,
		x * y * ( 1.0f - cosTheta ) + z * sinTheta,		y * y * ( 1.0f - cosTheta ) + cosTheta,			z * y * ( 1.0f - cosTheta ) - x * sinTheta,		0.0f,
		x * z * ( 1.0f - cosTheta ) - y * sinTheta,		y * z * ( 1.0f - cosTheta ) + x * sinTheta,		z * z * ( 1.0f - cosTheta ) + cosTheta,			0.0f,
		0.0f,											0.0f,											0.0f,											1.0f
	);
}

// static
Matrix4f Matrix4f::rotation( const Quat4f& q )
{
	Quat4f qq = q.normalized();

	float xx = qq.x() * qq.x();
	float yy = qq.y() * qq.y();
	float zz = qq.z() * qq.z();

	float xy = qq.x() * qq.y();
	float zw = qq.z() * qq.w();

	float xz = qq.x() * qq.z();
	float yw = qq.y() * qq.w();

	float yz = qq.y() * qq.z();
	float xw = qq.x() * qq.w();

	return Matrix4f
	(
		1.0f - 2.0f * ( yy + zz ),		2.0f * ( xy - zw ),				2.0f * ( xz + yw ),				0.0f,
		2.0f * ( xy + zw ),				1.0f - 2.0f * ( xx + zz ),		2.0f * ( yz - xw ),				0.0f,
		2.0f * ( xz - yw ),				2.0f * ( yz + xw ),				1.0f - 2.0f * ( xx + yy ),		0.0f,
		0.0f,							0.0f,							0.0f,							1.0f
	);
}

// static
Matrix4f Matrix4f::scaling( float sx, float sy, float sz )
{
	return Matrix4f
	(
		sx, 0, 0, 0,
		0, sy, 0, 0,
		0, 0, sz, 0,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::uniformScaling( float s )
{
	return Matrix4f
	(
		s, 0, 0, 0,
		0, s, 0, 0,
		0, 0, s, 0,
		0, 0, 0, 1
	);
}

// static
Matrix4f Matrix4f::randomRotation( float u0, float u1, float u2 )
{
	return Matrix4f::rotation( Quat4f::randomRotation( u0, u1, u2 ) );
}

// static
Matrix4f Matrix4f::lookAt( const Vector3f& eye, const Vector3f& center, const Vector3f& up )
{
	// z is negative forward
	Vector3f z = ( eye - center ).normalized();
	Vector3f y = up;
	Vector3f x = Vector3f::cross( y, z );

	// the x, y, and z vectors define the orthonormal coordinate system
	// the affine part defines the overall translation
	Matrix4f view;

	view.setRow( 0, Vector4f( x, -Vector3f::dot( x, eye ) ) );
	view.setRow( 1, Vector4f( y, -Vector3f::dot( y, eye ) ) );
	view.setRow( 2, Vector4f( z, -Vector3f::dot( z, eye ) ) );
	view.setRow( 3, Vector4f( 0, 0, 0, 1 ) );

	return view;
}

// static
Matrix4f Matrix4f::orthographicProjection( float width, float height, float zNear, float zFar, bool directX )
{
	Matrix4f m;

	m( 0, 0 ) = 2.0f / width;
	m( 1, 1 ) = 2.0f / height;
	m( 3, 3 ) = 1.0f;

	m( 0, 3 ) = -1;
	m( 1, 3 ) = -1;

	if( directX )
	{
		m( 2, 2 ) = 1.0f / ( zNear - zFar );
		m( 2, 3 ) = zNear / ( zNear - zFar );
	}
	else
	{
		m( 2, 2 ) = 2.0f / ( zNear - zFar );
		m( 2, 3 ) = ( zNear + zFar ) / ( zNear - zFar );
	}

	return m;
}

// static
Matrix4f Matrix4f::orthographicProjection( float left, float right, float bottom, float top, float zNear, float zFar, bool directX )
{
	Matrix4f m;

	m( 0, 0 ) = 2.0f / ( right - left );
	m( 1, 1 ) = 2.0f / ( top - bottom );
	m( 3, 3 ) = 1.0f;

	m( 0, 3 ) = ( left + right ) / ( left - right );
	m( 1, 3 ) = ( top + bottom ) / ( bottom - top );

	if( directX )
	{
		m( 2, 2 ) = 1.0f / ( zNear - zFar );
		m( 2, 3 ) = zNear / ( zNear - zFar );
	}
	else
	{
		m( 2, 2 ) = 2.0f / ( zNear - zFar );
		m( 2, 3 ) = ( zNear + zFar ) / ( zNear - zFar );
	}

	return m;
}

// static
Matrix4f Matrix4f::perspectiveProjection( float fLeft, float fRight,
										 float fBottom, float fTop,
										 float fZNear, float fZFar,
										 bool directX )
{
	Matrix4f projection; // zero matrix

	projection( 0, 0 ) = ( 2.0f * fZNear ) / ( fRight - fLeft );
	projection( 1, 1 ) = ( 2.0f * fZNear ) / ( fTop - fBottom );
	projection( 0, 2 ) = ( fRight + fLeft ) / ( fRight - fLeft );
	projection( 1, 2 ) = ( fTop + fBottom ) / ( fTop - fBottom );
	projection( 3, 2 ) = -1;

	if( directX )
	{
		projection( 2, 2 ) = fZFar / ( fZNear - fZFar);
		projection( 2, 3 ) = ( fZNear * fZFar ) / ( fZNear - fZFar );
	}
	else
	{
		projection( 2, 2 ) = ( fZNear + fZFar ) / ( fZNear - fZFar );
		projection( 2, 3 ) = ( 2.0f * fZNear * fZFar ) / ( fZNear - fZFar );
	}

	return projection;
}

// static
Matrix4f Matrix4f::perspectiveProjection( float fovYRadians, float aspect, float zNear, float zFar, bool directX )
{
	Matrix4f m; // zero matrix

	float yScale = 1.f / tanf( 0.5f * fovYRadians );
	float xScale = yScale / aspect;

	m( 0, 0 ) = xScale;
	m( 1, 1 ) = yScale;
	m( 3, 2 ) = -1;

	if( directX )
	{
		m( 2, 2 ) = zFar / ( zNear - zFar );
		m( 2, 3 ) = zNear * zFar / ( zNear - zFar );
	}
	else
	{
		m( 2, 2 ) = ( zFar + zNear ) / ( zNear - zFar );
		m( 2, 3 ) = 2.f * zFar * zNear / ( zNear - zFar );
	}

	return m;
}

// static
Matrix4f Matrix4f::infinitePerspectiveProjection( float fLeft, float fRight,
												 float fBottom, float fTop,
												 float fZNear, bool directX )
{
	Matrix4f projection;

	projection( 0, 0 ) = ( 2.0f * fZNear ) / ( fRight - fLeft );
	projection( 1, 1 ) = ( 2.0f * fZNear ) / ( fTop - fBottom );
	projection( 0, 2 ) = ( fRight + fLeft ) / ( fRight - fLeft );
	projection( 1, 2 ) = ( fTop + fBottom ) / ( fTop - fBottom );
	projection( 3, 2 ) = -1;

	// infinite view frustum
	// just take the limit as far --> inf of the regular frustum
	if( directX )
	{
		projection( 2, 2 ) = -1.0f;
		projection( 2, 3 ) = -fZNear;		
	}
	else
	{
		projection( 2, 2 ) = -1.0f;
		projection( 2, 3 ) = -2.0f * fZNear;
	}

	return projection;
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Vector4f operator * ( const Matrix4f& m, const Vector4f& v )
{
	Vector4f output( 0, 0, 0, 0 );

	for( int i = 0; i < 4; ++i )
	{
		for( int j = 0; j < 4; ++j )
		{
			output[ i ] += m( i, j ) * v[ j ];
		}
	}

	return output;
}

Matrix4f operator * ( const Matrix4f& x, const Matrix4f& y )
{
	Matrix4f product; // zeroes

	for( int i = 0; i < 4; ++i )
	{
		for( int j = 0; j < 4; ++j )
		{
			for( int k = 0; k < 4; ++k )
			{
				product( i, k ) += x( i, j ) * y( j, k );
			}
		}
	}

	return product;
}

// static
const Quat4f Quat4f::ZERO = Quat4f( 0, 0, 0, 0 );

// static
const Quat4f Quat4f::IDENTITY = Quat4f( 1, 0, 0, 0 );

Quat4f::Quat4f()
{
	m_elements[ 0 ] = 0;
	m_elements[ 1 ] = 0;
	m_elements[ 2 ] = 0;
	m_elements[ 3 ] = 0;
}

Quat4f::Quat4f( float w, float x, float y, float z )
{
	m_elements[ 0 ] = w;
	m_elements[ 1 ] = x;
	m_elements[ 2 ] = y;
	m_elements[ 3 ] = z;
}

Quat4f::Quat4f( const Quat4f& rq )
{
	m_elements[ 0 ] = rq.m_elements[ 0 ];
	m_elements[ 1 ] = rq.m_elements[ 1 ];
	m_elements[ 2 ] = rq.m_elements[ 2 ];
	m_elements[ 3 ] = rq.m_elements[ 3 ];
}

Quat4f& Quat4f::operator = ( const Quat4f& rq )
{
	if( this != ( &rq ) )
	{
		m_elements[ 0 ] = rq.m_elements[ 0 ];
		m_elements[ 1 ] = rq.m_elements[ 1 ];
		m_elements[ 2 ] = rq.m_elements[ 2 ];
		m_elements[ 3 ] = rq.m_elements[ 3 ];
	}
    return( *this );
}

Quat4f::Quat4f( const Vector3f& v )
{
	m_elements[ 0 ] = 0;
	m_elements[ 1 ] = v[ 0 ];
	m_elements[ 2 ] = v[ 1 ];
	m_elements[ 3 ] = v[ 2 ];
}

Quat4f::Quat4f( const Vector4f& v )
{
	m_elements[ 0 ] = v[ 0 ];
	m_elements[ 1 ] = v[ 1 ];
	m_elements[ 2 ] = v[ 2 ];
	m_elements[ 3 ] = v[ 3 ];
}

const float& Quat4f::operator [] ( int i ) const
{
	return m_elements[ i ];
}

float& Quat4f::operator [] ( int i )
{
	return m_elements[ i ];
}

float Quat4f::w() const
{
	return m_elements[ 0 ];
}

float Quat4f::x() const
{
	return m_elements[ 1 ];
}

float Quat4f::y() const
{
	return m_elements[ 2 ];
}

float Quat4f::z() const
{
	return m_elements[ 3 ];
}

Vector3f Quat4f::xyz() const
{
	return Vector3f
	(
		m_elements[ 1 ],
		m_elements[ 2 ],
		m_elements[ 3 ]
	);
}

Vector4f Quat4f::wxyz() const
{
	return Vector4f
	(
		m_elements[ 0 ],
		m_elements[ 1 ],
		m_elements[ 2 ],
		m_elements[ 3 ]
	);
}

float Quat4f::abs() const
{
	return sqrt( absSquared() );	
}

float Quat4f::absSquared() const
{
	return
	(
		m_elements[ 0 ] * m_elements[ 0 ] +
		m_elements[ 1 ] * m_elements[ 1 ] +
		m_elements[ 2 ] * m_elements[ 2 ] +
		m_elements[ 3 ] * m_elements[ 3 ]
	);
}

void Quat4f::normalize()
{
	float reciprocalAbs = 1.f / abs();

	m_elements[ 0 ] *= reciprocalAbs;
	m_elements[ 1 ] *= reciprocalAbs;
	m_elements[ 2 ] *= reciprocalAbs;
	m_elements[ 3 ] *= reciprocalAbs;
}

Quat4f Quat4f::normalized() const
{
	Quat4f q( *this );
	q.normalize();
	return q;
}

void Quat4f::conjugate()
{
	m_elements[ 1 ] = -m_elements[ 1 ];
	m_elements[ 2 ] = -m_elements[ 2 ];
	m_elements[ 3 ] = -m_elements[ 3 ];
}

Quat4f Quat4f::conjugated() const
{
	return Quat4f
	(
		 m_elements[ 0 ],
		-m_elements[ 1 ],
		-m_elements[ 2 ],
		-m_elements[ 3 ]
	);
}

void Quat4f::invert()
{
	Quat4f inverse = conjugated() * ( 1.0f / absSquared() );

	m_elements[ 0 ] = inverse.m_elements[ 0 ];
	m_elements[ 1 ] = inverse.m_elements[ 1 ];
	m_elements[ 2 ] = inverse.m_elements[ 2 ];
	m_elements[ 3 ] = inverse.m_elements[ 3 ];
}

Quat4f Quat4f::inverse() const
{
	return conjugated() * ( 1.0f / absSquared() );
}


Quat4f Quat4f::log() const
{
	float len =
		sqrt
		(
			m_elements[ 1 ] * m_elements[ 1 ] +
			m_elements[ 2 ] * m_elements[ 2 ] +
			m_elements[ 3 ] * m_elements[ 3 ]
		);

	if( len < 1e-6 )
	{
		return Quat4f( 0, m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ] );
	}
	else
	{
		float coeff = acos( m_elements[ 0 ] ) / len;
		return Quat4f( 0, m_elements[ 1 ] * coeff, m_elements[ 2 ] * coeff, m_elements[ 3 ] * coeff );
	}
}

Quat4f Quat4f::exp() const
{
	float theta =
		sqrt
		(
			m_elements[ 1 ] * m_elements[ 1 ] +
			m_elements[ 2 ] * m_elements[ 2 ] +
			m_elements[ 3 ] * m_elements[ 3 ]
		);

	if( theta < 1e-6 )
	{
		return Quat4f( cos( theta ), m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ] );
	}
	else
	{
		float coeff = sin( theta ) / theta;
		return Quat4f( cos( theta ), m_elements[ 1 ] * coeff, m_elements[ 2 ] * coeff, m_elements[ 3 ] * coeff );		
	}
}

Vector3f Quat4f::getAxisAngle( float* radiansOut )
{
	float theta = acos( w() ) * 2;
	float vectorNorm = sqrt( x() * x() + y() * y() + z() * z() );
	float reciprocalVectorNorm = 1.f / vectorNorm;

	*radiansOut = theta;
	return Vector3f
	(
		x() * reciprocalVectorNorm,
		y() * reciprocalVectorNorm,
		z() * reciprocalVectorNorm
	);
}

void Quat4f::setAxisAngle( float radians, const Vector3f& axis )
{
	m_elements[ 0 ] = cos( radians / 2 );

	float sinHalfTheta = sin( radians / 2 );
	float vectorNorm = axis.length();
	float reciprocalVectorNorm = 1.f / vectorNorm;

	m_elements[ 1 ] = axis.x() * sinHalfTheta * reciprocalVectorNorm;
	m_elements[ 2 ] = axis.y() * sinHalfTheta * reciprocalVectorNorm;
	m_elements[ 3 ] = axis.z() * sinHalfTheta * reciprocalVectorNorm;
}

void Quat4f::print()
{
	printf( "< %.4f + %.4f i + %.4f j + %.4f k >\n",
		m_elements[ 0 ], m_elements[ 1 ], m_elements[ 2 ], m_elements[ 3 ] );
}

// static
float Quat4f::dot( const Quat4f& q0, const Quat4f& q1 )
{
	return
	(
		q0.w() * q1.w() +
		q0.x() * q1.x() +
		q0.y() * q1.y() +
		q0.z() * q1.z()
	);
}

// static
Quat4f Quat4f::lerp( const Quat4f& q0, const Quat4f& q1, float alpha )
{
	return( ( q0 + alpha * ( q1 - q0 ) ).normalized() );
}

// static
Quat4f Quat4f::slerp( const Quat4f& a, const Quat4f& b, float t, bool allowFlip )
{
	float cosAngle = Quat4f::dot( a, b );

	float c1;
	float c2;

	// Linear interpolation for close orientations
	if( ( 1.0f - fabs( cosAngle ) ) < 0.01f )
	{
		c1 = 1.0f - t;
		c2 = t;
	}
	else
	{
		// Spherical interpolation
		float angle = acos( fabs( cosAngle ) );
		float sinAngle = sin( angle );
		c1 = sin( angle * ( 1.0f - t ) ) / sinAngle;
		c2 = sin( angle * t ) / sinAngle;
	}

	// Use the shortest path
	if( allowFlip && ( cosAngle < 0.0f ) )
	{
		c1 = -c1;
	}

	return Quat4f( c1 * a[ 0 ] + c2 * b[ 0 ], c1 * a[ 1 ] + c2 * b[ 1 ], c1 * a[ 2 ] + c2 * b[ 2 ], c1 * a[ 3 ] + c2 * b[ 3 ] );
}

// static
Quat4f Quat4f::squad( const Quat4f& a, const Quat4f& tanA, const Quat4f& tanB, const Quat4f& b, float t )
{
	Quat4f ab = Quat4f::slerp( a, b, t );
	Quat4f tangent = Quat4f::slerp( tanA, tanB, t, false );
	return Quat4f::slerp( ab, tangent, 2.0f * t * ( 1.0f - t ), false );
}

// static
Quat4f Quat4f::cubicInterpolate( const Quat4f& q0, const Quat4f& q1, const Quat4f& q2, const Quat4f& q3, float t )
{
	// geometric construction:
	//            t
	//   (t+1)/2     t/2
	// t+1        t	        t-1

	// bottom level
	Quat4f q0q1 = Quat4f::slerp( q0, q1, t + 1 );
	Quat4f q1q2 = Quat4f::slerp( q1, q2, t );
	Quat4f q2q3 = Quat4f::slerp( q2, q3, t - 1 );

	// middle level
	Quat4f q0q1_q1q2 = Quat4f::slerp( q0q1, q1q2, 0.5f * ( t + 1 ) );
	Quat4f q1q2_q2q3 = Quat4f::slerp( q1q2, q2q3, 0.5f * t );

	// top level
	return Quat4f::slerp( q0q1_q1q2, q1q2_q2q3, t );
}

// static
Quat4f Quat4f::logDifference( const Quat4f& a, const Quat4f& b )
{
	Quat4f diff = a.inverse() * b;
	diff.normalize();
	return diff.log();
}

// static
Quat4f Quat4f::squadTangent( const Quat4f& before, const Quat4f& center, const Quat4f& after )
{
	Quat4f l1 = Quat4f::logDifference( center, before );
	Quat4f l2 = Quat4f::logDifference( center, after );
	
	Quat4f e;
	for( int i = 0; i < 4; ++i )
	{
		e[ i ] = -0.25f * ( l1[ i ] + l2[ i ] );
	}
	e = center * ( e.exp() );

	return e;
}

// static
Quat4f Quat4f::fromRotationMatrix( const Matrix3f& m )
{
	float x;
	float y;
	float z;
	float w;

	// Compute one plus the trace of the matrix
	float onePlusTrace = 1.0f + m( 0, 0 ) + m( 1, 1 ) + m( 2, 2 );

	if( onePlusTrace > 1e-5 )
	{
		// Direct computation
		float s = sqrt( onePlusTrace ) * 2.0f;
		x = ( m( 2, 1 ) - m( 1, 2 ) ) / s;
		y = ( m( 0, 2 ) - m( 2, 0 ) ) / s;
		z = ( m( 1, 0 ) - m( 0, 1 ) ) / s;
		w = 0.25f * s;
	}
	else
	{
		// Computation depends on major diagonal term
		if( ( m( 0, 0 ) > m( 1, 1 ) ) & ( m( 0, 0 ) > m( 2, 2 ) ) )
		{
			float s = sqrt( 1.0f + m( 0, 0 ) - m( 1, 1 ) - m( 2, 2 ) ) * 2.0f;
			x = 0.25f * s;
			y = ( m( 0, 1 ) + m( 1, 0 ) ) / s;
			z = ( m( 0, 2 ) + m( 2, 0 ) ) / s;
			w = ( m( 1, 2 ) - m( 2, 1 ) ) / s;
		}
		else if( m( 1, 1 ) > m( 2, 2 ) )
		{
			float s = sqrt( 1.0f + m( 1, 1 ) - m( 0, 0 ) - m( 2, 2 ) ) * 2.0f;
			x = ( m( 0, 1 ) + m( 1, 0 ) ) / s;
			y = 0.25f * s;
			z = ( m( 1, 2 ) + m( 2, 1 ) ) / s;
			w = ( m( 0, 2 ) - m( 2, 0 ) ) / s;
		}
		else
		{
			float s = sqrt( 1.0f + m( 2, 2 ) - m( 0, 0 ) - m( 1, 1 ) ) * 2.0f;
			x = ( m( 0, 2 ) + m( 2, 0 ) ) / s;
			y = ( m( 1, 2 ) + m( 2, 1 ) ) / s;
			z = 0.25f * s;
			w = ( m( 0, 1 ) - m( 1, 0 ) ) / s;
		}
	}

	Quat4f q( w, x, y, z );
	return q.normalized();
}

// static
Quat4f Quat4f::fromRotatedBasis( const Vector3f& x, const Vector3f& y, const Vector3f& z )
{
	return fromRotationMatrix( Matrix3f( x, y, z ) );
}

// static
Quat4f Quat4f::randomRotation( float u0, float u1, float u2 )
{
	float z = u0;
	float theta = static_cast< float >( 2.f * M_PI * u1 );
	float r = sqrt( 1.f - z * z );
	float w = static_cast< float >( M_PI * u2 );

	return Quat4f
	(
		cos( w ),
		sin( w ) * cos( theta ) * r,
		sin( w ) * sin( theta ) * r,
		sin( w ) * z
	);
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Quat4f operator + ( const Quat4f& q0, const Quat4f& q1 )
{
	return Quat4f
	(
		q0.w() + q1.w(),
		q0.x() + q1.x(),
		q0.y() + q1.y(),
		q0.z() + q1.z()
	);
}

Quat4f operator - ( const Quat4f& q0, const Quat4f& q1 )
{
	return Quat4f
	(
		q0.w() - q1.w(),
		q0.x() - q1.x(),
		q0.y() - q1.y(),
		q0.z() - q1.z()
	);
}

Quat4f operator * ( const Quat4f& q0, const Quat4f& q1 )
{
	return Quat4f
	(
		q0.w() * q1.w() - q0.x() * q1.x() - q0.y() * q1.y() - q0.z() * q1.z(),
		q0.w() * q1.x() + q0.x() * q1.w() + q0.y() * q1.z() - q0.z() * q1.y(),
		q0.w() * q1.y() - q0.x() * q1.z() + q0.y() * q1.w() + q0.z() * q1.x(),
		q0.w() * q1.z() + q0.x() * q1.y() - q0.y() * q1.x() + q0.z() * q1.w()
	);
}

Quat4f operator * ( float f, const Quat4f& q )
{
	return Quat4f
	(
		f * q.w(),
		f * q.x(),
		f * q.y(),
		f * q.z()
	);
}

Quat4f operator * ( const Quat4f& q, float f )
{
	return Quat4f
	(
		f * q.w(),
		f * q.x(),
		f * q.y(),
		f * q.z()
	);
}


// static
const Vector2f Vector2f::ZERO = Vector2f( 0, 0 );

// static
const Vector2f Vector2f::UP = Vector2f( 0, 1 );

// static
const Vector2f Vector2f::RIGHT = Vector2f( 1, 0 );

Vector2f::Vector2f( float f )
{
    m_elements[0] = f;
    m_elements[1] = f;
}

Vector2f::Vector2f( float x, float y )
{
    m_elements[0] = x;
    m_elements[1] = y;
}

Vector2f::Vector2f( const Vector2f& rv )
{
    m_elements[0] = rv[0];
    m_elements[1] = rv[1];
}

Vector2f& Vector2f::operator = ( const Vector2f& rv )
{
 	if( this != &rv )
	{
        m_elements[0] = rv[0];
        m_elements[1] = rv[1];
    }
    return *this;
}

const float& Vector2f::operator [] ( int i ) const
{
    return m_elements[i];
}

float& Vector2f::operator [] ( int i )
{
    return m_elements[i];
}

float& Vector2f::x()
{
    return m_elements[0];
}

float& Vector2f::y()
{
    return m_elements[1];
}

float Vector2f::x() const
{
    return m_elements[0];
}	

float Vector2f::y() const
{
    return m_elements[1];
}

Vector2f Vector2f::xy() const
{
    return *this;
}

Vector2f Vector2f::yx() const
{
    return Vector2f( m_elements[1], m_elements[0] );
}

Vector2f Vector2f::xx() const
{
    return Vector2f( m_elements[0], m_elements[0] );
}

Vector2f Vector2f::yy() const
{
    return Vector2f( m_elements[1], m_elements[1] );
}

Vector2f Vector2f::normal() const
{
    return Vector2f( -m_elements[1], m_elements[0] );
}

float Vector2f::abs() const
{
    return sqrt(absSquared());
}

float Vector2f::absSquared() const
{
    return m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1];
}

void Vector2f::normalize()
{
    float norm = abs();
    m_elements[0] /= norm;
    m_elements[1] /= norm;
}

Vector2f Vector2f::normalized() const
{
    float norm = abs();
    return Vector2f( m_elements[0] / norm, m_elements[1] / norm );
}

void Vector2f::negate()
{
    m_elements[0] = -m_elements[0];
    m_elements[1] = -m_elements[1];
}

Vector2f::operator const float* () const
{
    return m_elements;
}

Vector2f::operator float* ()
{
    return m_elements;
}

void Vector2f::print() const
{
	printf( "< %.4f, %.4f >\n",
		m_elements[0], m_elements[1] );
}

Vector2f& Vector2f::operator += ( const Vector2f& v )
{
	m_elements[ 0 ] += v.m_elements[ 0 ];
	m_elements[ 1 ] += v.m_elements[ 1 ];
	return *this;
}

Vector2f& Vector2f::operator -= ( const Vector2f& v )
{
	m_elements[ 0 ] -= v.m_elements[ 0 ];
	m_elements[ 1 ] -= v.m_elements[ 1 ];
	return *this;
}

Vector2f& Vector2f::operator *= ( float f )
{
	m_elements[ 0 ] *= f;
	m_elements[ 1 ] *= f;
	return *this;
}

// static
float Vector2f::dot( const Vector2f& v0, const Vector2f& v1 )
{
    return v0[0] * v1[0] + v0[1] * v1[1];
}

// static
Vector3f Vector2f::cross( const Vector2f& v0, const Vector2f& v1 )
{
	return Vector3f
		(
			0,
			0,
			v0.x() * v1.y() - v0.y() * v1.x()
		);
}

// static
Vector2f Vector2f::lerp( const Vector2f& v0, const Vector2f& v1, float alpha )
{
	return alpha * ( v1 - v0 ) + v0;
}

//////////////////////////////////////////////////////////////////////////
// Operator overloading
//////////////////////////////////////////////////////////////////////////

Vector2f operator + ( const Vector2f& v0, const Vector2f& v1 )
{
    return Vector2f( v0.x() + v1.x(), v0.y() + v1.y() );
}

Vector2f operator - ( const Vector2f& v0, const Vector2f& v1 )
{
    return Vector2f( v0.x() - v1.x(), v0.y() - v1.y() );
}

Vector2f operator * ( const Vector2f& v0, const Vector2f& v1 )
{
    return Vector2f( v0.x() * v1.x(), v0.y() * v1.y() );
}

Vector2f operator / ( const Vector2f& v0, const Vector2f& v1 )
{
    return Vector2f( v0.x() / v1.x(), v0.y() / v1.y() );
}

Vector2f operator - ( const Vector2f& v )
{
    return Vector2f( -v.x(), -v.y() );
}

Vector2f operator * ( float f, const Vector2f& v )
{
    return Vector2f( f * v.x(), f * v.y() );
}

Vector2f operator * ( const Vector2f& v, float f )
{
    return Vector2f( f * v.x(), f * v.y() );
}

Vector2f operator / ( const Vector2f& v, float f )
{
    return Vector2f( v.x() / f, v.y() / f );
}

bool operator == ( const Vector2f& v0, const Vector2f& v1 )
{
    return( v0.x() == v1.x() && v0.y() == v1.y() );
}

bool operator != ( const Vector2f& v0, const Vector2f& v1 )
{
    return !( v0 == v1 );
}


// static
const Vector3f Vector3f::ZERO = Vector3f( 0, 0, 0 );

// static
const Vector3f Vector3f::UP = Vector3f( 0, 1, 0 );

// static
const Vector3f Vector3f::RIGHT = Vector3f( 1, 0, 0 );

// static
const Vector3f Vector3f::FORWARD = Vector3f( 0, 0, -1 );

Vector3f::Vector3f( float f )
{
    m_elements[0] = f;
    m_elements[1] = f;
    m_elements[2] = f;
}

Vector3f::Vector3f( float x, float y, float z )
{
    m_elements[0] = x;
    m_elements[1] = y;
    m_elements[2] = z;
}

Vector3f::Vector3f( const Vector2f& xy, float z )
{
	m_elements[0] = xy.x();
	m_elements[1] = xy.y();
	m_elements[2] = z;
}

Vector3f::Vector3f( float x, const Vector2f& yz )
{
	m_elements[0] = x;
	m_elements[1] = yz.x();
	m_elements[2] = yz.y();
}

Vector3f::Vector3f( const Vector3f& rv )
{
    m_elements[0] = rv[0];
    m_elements[1] = rv[1];
    m_elements[2] = rv[2];
}

Vector3f& Vector3f::operator = ( const Vector3f& rv )
{
    if( this != &rv )
    {
        m_elements[0] = rv[0];
        m_elements[1] = rv[1];
        m_elements[2] = rv[2];
    }
    return *this;
}

const float& Vector3f::operator [] ( int i ) const
{
    return m_elements[i];
}

float& Vector3f::operator [] ( int i )
{
    return m_elements[i];
}

float& Vector3f::x()
{
    return m_elements[0];
}

float& Vector3f::y()
{
    return m_elements[1];
}

float& Vector3f::z()
{
    return m_elements[2];
}

float Vector3f::x() const
{
    return m_elements[0];
}

float Vector3f::y() const
{
    return m_elements[1];
}

float Vector3f::z() const
{
    return m_elements[2];
}

Vector2f Vector3f::xy() const
{
	return Vector2f( m_elements[0], m_elements[1] );
}

Vector2f Vector3f::xz() const
{
	return Vector2f( m_elements[0], m_elements[2] );
}

Vector2f Vector3f::yz() const
{
	return Vector2f( m_elements[1], m_elements[2] );
}

Vector3f Vector3f::xyz() const
{
	return Vector3f( m_elements[0], m_elements[1], m_elements[2] );
}

Vector3f Vector3f::yzx() const
{
	return Vector3f( m_elements[1], m_elements[2], m_elements[0] );
}

Vector3f Vector3f::zxy() const
{
	return Vector3f( m_elements[2], m_elements[0], m_elements[1] );
}

float Vector3f::length() const
{
	return sqrt( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] );
}

float Vector3f::squaredLength() const
{
    return
        (
            m_elements[0] * m_elements[0] +
            m_elements[1] * m_elements[1] +
            m_elements[2] * m_elements[2]
        );
}

void Vector3f::normalize()
{
	float norm = length();
	m_elements[0] /= norm;
	m_elements[1] /= norm;
	m_elements[2] /= norm;
}

Vector3f Vector3f::normalized() const
{
	float norm = length();
	return Vector3f
		(
			m_elements[0] / norm,
			m_elements[1] / norm,
			m_elements[2] / norm
		);
}

Vector2f Vector3f::homogenized() const
{
	return Vector2f
		(
			m_elements[ 0 ] / m_elements[ 2 ],
			m_elements[ 1 ] / m_elements[ 2 ]
		);
}

void Vector3f::negate()
{
	m_elements[0] = -m_elements[0];
	m_elements[1] = -m_elements[1];
	m_elements[2] = -m_elements[2];
}

Vector3f::operator const float* () const
{
    return m_elements;
}

Vector3f::operator float* ()
{
    return m_elements;
}

void Vector3f::print() const
{
	printf( "< %.4f, %.4f, %.4f >\n",
		m_elements[0], m_elements[1], m_elements[2] );
}

Vector3f& Vector3f::operator += ( const Vector3f& v )
{
	m_elements[ 0 ] += v.m_elements[ 0 ];
	m_elements[ 1 ] += v.m_elements[ 1 ];
	m_elements[ 2 ] += v.m_elements[ 2 ];
	return *this;
}

Vector3f& Vector3f::operator -= ( const Vector3f& v )
{
	m_elements[ 0 ] -= v.m_elements[ 0 ];
	m_elements[ 1 ] -= v.m_elements[ 1 ];
	m_elements[ 2 ] -= v.m_elements[ 2 ];
	return *this;
}

Vector3f& Vector3f::operator *= ( float f )
{
	m_elements[ 0 ] *= f;
	m_elements[ 1 ] *= f;
	m_elements[ 2 ] *= f;
	return *this;
}

Vector3f Vector3f::operator%(Vector3f &b)
{
	return Vector3f(m_elements[1] * b.z() - m_elements[2] * b.y(), m_elements[2] * b.x() - m_elements[0] * b.z(), m_elements[0] * b.y() - m_elements[1] * b.x());
}

// static
float Vector3f::dot( const Vector3f& v0, const Vector3f& v1 )
{
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

// static
Vector3f Vector3f::cross( const Vector3f& v0, const Vector3f& v1 )
{
    return Vector3f
        (
            v0.y() * v1.z() - v0.z() * v1.y(),
            v0.z() * v1.x() - v0.x() * v1.z(),
            v0.x() * v1.y() - v0.y() * v1.x()
        );
}

// static
Vector3f Vector3f::lerp( const Vector3f& v0, const Vector3f& v1, float alpha )
{
	return alpha * ( v1 - v0 ) + v0;
}

// static
Vector3f Vector3f::cubicInterpolate( const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const Vector3f& p3, float t )
{
	// geometric construction:
	//            t
	//   (t+1)/2     t/2
	// t+1        t	        t-1

	// bottom level
	Vector3f p0p1 = Vector3f::lerp( p0, p1, t + 1 );
	Vector3f p1p2 = Vector3f::lerp( p1, p2, t );
	Vector3f p2p3 = Vector3f::lerp( p2, p3, t - 1 );

	// middle level
	Vector3f p0p1_p1p2 = Vector3f::lerp( p0p1, p1p2, 0.5f * ( t + 1 ) );
	Vector3f p1p2_p2p3 = Vector3f::lerp( p1p2, p2p3, 0.5f * t );

	// top level
	return Vector3f::lerp( p0p1_p1p2, p1p2_p2p3, t );
}

Vector3f operator + ( const Vector3f& v0, const Vector3f& v1 )
{
    return Vector3f( v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2] );
}

Vector3f operator - ( const Vector3f& v0, const Vector3f& v1 )
{
    return Vector3f( v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2] );
}

Vector3f operator * ( const Vector3f& v0, const Vector3f& v1 )
{
    return Vector3f( v0[0] * v1[0], v0[1] * v1[1], v0[2] * v1[2] );
}

Vector3f operator / ( const Vector3f& v0, const Vector3f& v1 )
{
    return Vector3f( v0[0] / v1[0], v0[1] / v1[1], v0[2] / v1[2] );
}

Vector3f operator - ( const Vector3f& v )
{
    return Vector3f( -v[0], -v[1], -v[2] );
}

Vector3f operator * ( float f, const Vector3f& v )
{
    return Vector3f( v[0] * f, v[1] * f, v[2] * f );
}

Vector3f operator * ( const Vector3f& v, float f )
{
    return Vector3f( v[0] * f, v[1] * f, v[2] * f );
}

Vector3f operator / ( const Vector3f& v, float f )
{
    return Vector3f( v[0] / f, v[1] / f, v[2] / f );
}

bool operator == ( const Vector3f& v0, const Vector3f& v1 )
{
    return( v0.x() == v1.x() && v0.y() == v1.y() && v0.z() == v1.z() );
}

bool operator != ( const Vector3f& v0, const Vector3f& v1 )
{
    return !( v0 == v1 );
}

Vector4f::Vector4f( float f )
{
	m_elements[ 0 ] = f;
	m_elements[ 1 ] = f;
	m_elements[ 2 ] = f;
	m_elements[ 3 ] = f;
}

Vector4f::Vector4f( float fx, float fy, float fz, float fw )
{
	m_elements[0] = fx;
	m_elements[1] = fy;
	m_elements[2] = fz;
	m_elements[3] = fw;
}

Vector4f::Vector4f( float buffer[ 4 ] )
{
	m_elements[ 0 ] = buffer[ 0 ];
	m_elements[ 1 ] = buffer[ 1 ];
	m_elements[ 2 ] = buffer[ 2 ];
	m_elements[ 3 ] = buffer[ 3 ];
}

Vector4f::Vector4f( const Vector2f& xy, float z, float w )
{
	m_elements[0] = xy.x();
	m_elements[1] = xy.y();
	m_elements[2] = z;
	m_elements[3] = w;
}

Vector4f::Vector4f( float x, const Vector2f& yz, float w )
{
	m_elements[0] = x;
	m_elements[1] = yz.x();
	m_elements[2] = yz.y();
	m_elements[3] = w;
}

Vector4f::Vector4f( float x, float y, const Vector2f& zw )
{
	m_elements[0] = x;
	m_elements[1] = y;
	m_elements[2] = zw.x();
	m_elements[3] = zw.y();
}

Vector4f::Vector4f( const Vector2f& xy, const Vector2f& zw )
{
	m_elements[0] = xy.x();
	m_elements[1] = xy.y();
	m_elements[2] = zw.x();
	m_elements[3] = zw.y();
}

Vector4f::Vector4f( const Vector3f& xyz, float w )
{
	m_elements[0] = xyz.x();
	m_elements[1] = xyz.y();
	m_elements[2] = xyz.z();
	m_elements[3] = w;
}

Vector4f::Vector4f( float x, const Vector3f& yzw )
{
	m_elements[0] = x;
	m_elements[1] = yzw.x();
	m_elements[2] = yzw.y();
	m_elements[3] = yzw.z();
}

Vector4f::Vector4f( const Vector4f& rv )
{
	m_elements[0] = rv.m_elements[0];
	m_elements[1] = rv.m_elements[1];
	m_elements[2] = rv.m_elements[2];
	m_elements[3] = rv.m_elements[3];
}

Vector4f& Vector4f::operator = ( const Vector4f& rv )
{
	if( this != &rv )
	{
		m_elements[0] = rv.m_elements[0];
		m_elements[1] = rv.m_elements[1];
		m_elements[2] = rv.m_elements[2];
		m_elements[3] = rv.m_elements[3];
	}
	return *this;
}

const float& Vector4f::operator [] ( int i ) const
{
	return m_elements[ i ];
}

float& Vector4f::operator [] ( int i )
{
	return m_elements[ i ];
}

float& Vector4f::x()
{
	return m_elements[ 0 ];
}

float& Vector4f::y()
{
	return m_elements[ 1 ];
}

float& Vector4f::z()
{
	return m_elements[ 2 ];
}

float& Vector4f::w()
{
	return m_elements[ 3 ];
}

float Vector4f::x() const
{
	return m_elements[0];
}

float Vector4f::y() const
{
	return m_elements[1];
}

float Vector4f::z() const
{
	return m_elements[2];
}

float Vector4f::w() const
{
	return m_elements[3];
}

Vector2f Vector4f::xy() const
{
	return Vector2f( m_elements[0], m_elements[1] );
}

Vector2f Vector4f::yz() const
{
	return Vector2f( m_elements[1], m_elements[2] );
}

Vector2f Vector4f::zw() const
{
	return Vector2f( m_elements[2], m_elements[3] );
}

Vector2f Vector4f::wx() const
{
	return Vector2f( m_elements[3], m_elements[0] );
}

Vector3f Vector4f::xyz() const
{
	return Vector3f( m_elements[0], m_elements[1], m_elements[2] );
}

Vector3f Vector4f::yzw() const
{
	return Vector3f( m_elements[1], m_elements[2], m_elements[3] );
}

Vector3f Vector4f::zwx() const
{
	return Vector3f( m_elements[2], m_elements[3], m_elements[0] );
}

Vector3f Vector4f::wxy() const
{
	return Vector3f( m_elements[3], m_elements[0], m_elements[1] );
}

Vector3f Vector4f::xyw() const
{
	return Vector3f( m_elements[0], m_elements[1], m_elements[3] );
}

Vector3f Vector4f::yzx() const
{
	return Vector3f( m_elements[1], m_elements[2], m_elements[0] );
}

Vector3f Vector4f::zwy() const
{
	return Vector3f( m_elements[2], m_elements[3], m_elements[1] );
}

Vector3f Vector4f::wxz() const
{
	return Vector3f( m_elements[3], m_elements[0], m_elements[2] );
}

float Vector4f::abs() const
{
	return sqrt( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] + m_elements[3] * m_elements[3] );
}

float Vector4f::absSquared() const
{
	return( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] + m_elements[3] * m_elements[3] );
}

void Vector4f::normalize()
{
	float norm = sqrt( m_elements[0] * m_elements[0] + m_elements[1] * m_elements[1] + m_elements[2] * m_elements[2] + m_elements[3] * m_elements[3] );
	m_elements[0] = m_elements[0] / norm;
	m_elements[1] = m_elements[1] / norm;
	m_elements[2] = m_elements[2] / norm;
	m_elements[3] = m_elements[3] / norm;
}

Vector4f Vector4f::normalized() const
{
	float length = abs();
	return Vector4f
		(
			m_elements[0] / length,
			m_elements[1] / length,
			m_elements[2] / length,
			m_elements[3] / length
		);
}

void Vector4f::homogenize()
{
	if( m_elements[3] != 0 )
	{
		m_elements[0] /= m_elements[3];
		m_elements[1] /= m_elements[3];
		m_elements[2] /= m_elements[3];
		m_elements[3] = 1;
	}
}

Vector4f Vector4f::homogenized() const
{
	if( m_elements[3] != 0 )
	{
		return Vector4f
			(
				m_elements[0] / m_elements[3],
				m_elements[1] / m_elements[3],
				m_elements[2] / m_elements[3],
				1
			);
	}
	else
	{
		return Vector4f
			(
				m_elements[0],
				m_elements[1],
				m_elements[2],
				m_elements[3]
			);
	}
}

void Vector4f::negate()
{
	m_elements[0] = -m_elements[0];
	m_elements[1] = -m_elements[1];
	m_elements[2] = -m_elements[2];
	m_elements[3] = -m_elements[3];
}

Vector4f::operator const float* () const
{
	return m_elements;
}

Vector4f::operator float* ()
{
	return m_elements;
}

void Vector4f::print() const
{
	printf( "< %.4f, %.4f, %.4f, %.4f >\n",
		m_elements[0], m_elements[1], m_elements[2], m_elements[3] );
}

// static
float Vector4f::dot( const Vector4f& v0, const Vector4f& v1 )
{
	return v0.x() * v1.x() + v0.y() * v1.y() + v0.z() * v1.z() + v0.w() * v1.w();
}

// static
Vector4f Vector4f::lerp( const Vector4f& v0, const Vector4f& v1, float alpha )
{
	return alpha * ( v1 - v0 ) + v0;
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Vector4f operator + ( const Vector4f& v0, const Vector4f& v1 )
{
	return Vector4f( v0.x() + v1.x(), v0.y() + v1.y(), v0.z() + v1.z(), v0.w() + v1.w() );
}

Vector4f operator - ( const Vector4f& v0, const Vector4f& v1 )
{
	return Vector4f( v0.x() - v1.x(), v0.y() - v1.y(), v0.z() - v1.z(), v0.w() - v1.w() );
}

Vector4f operator * ( const Vector4f& v0, const Vector4f& v1 )
{
	return Vector4f( v0.x() * v1.x(), v0.y() * v1.y(), v0.z() * v1.z(), v0.w() * v1.w() );
}

Vector4f operator / ( const Vector4f& v0, const Vector4f& v1 )
{
	return Vector4f( v0.x() / v1.x(), v0.y() / v1.y(), v0.z() / v1.z(), v0.w() / v1.w() );
}

Vector4f operator - ( const Vector4f& v )
{
	return Vector4f( -v.x(), -v.y(), -v.z(), -v.w() );
}

Vector4f operator * ( float f, const Vector4f& v )
{
	return Vector4f( f * v.x(), f * v.y(), f * v.z(), f * v.w() );
}

Vector4f operator * ( const Vector4f& v, float f )
{
	return Vector4f( f * v.x(), f * v.y(), f * v.z(), f * v.w() );
}

Vector4f operator / ( const Vector4f& v, float f )
{
    return Vector4f( v[0] / f, v[1] / f, v[2] / f, v[3] / f );
}

bool operator == ( const Vector4f& v0, const Vector4f& v1 )
{
    return( v0.x() == v1.x() && v0.y() == v1.y() && v0.z() == v1.z() && v0.w() == v1.w() );
}

bool operator != ( const Vector4f& v0, const Vector4f& v1 )
{
    return !( v0 == v1 );
}
