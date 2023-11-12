#ifndef VECMATH
#define VECMATH

#include <cstdio>

class Vector2f;

// 2x2 Matrix, stored in column major order (OpenGL style)
class Matrix2f
{
public:
    // Fill a 2x2 matrix with "fill", default to 0.
    Matrix2f(float fill = 0.f);
    Matrix2f(float m00, float m01,
             float m10, float m11);

    // setColumns = true ==> sets the columns of the matrix to be [v0 v1]
    // otherwise, sets the rows
    Matrix2f(const Vector2f &v0, const Vector2f &v1, bool setColumns = true);

    Matrix2f(const Matrix2f &rm);            // copy constructor
    Matrix2f &operator=(const Matrix2f &rm); // assignment operator
    // no destructor necessary

    const float &operator()(int i, int j) const;
    float &operator()(int i, int j);

    Vector2f getRow(int i) const;
    void setRow(int i, const Vector2f &v);

    Vector2f getCol(int j) const;
    void setCol(int j, const Vector2f &v);

    float determinant();
    Matrix2f inverse(bool *pbIsSingular = NULL, float epsilon = 0.f);

    void transpose();
    Matrix2f transposed() const;

    // ---- Utility ----
    operator float *(); // automatic type conversion for GL
    void print();

    static float determinant2x2(float m00, float m01,
                                float m10, float m11);

    static Matrix2f ones();
    static Matrix2f identity();
    static Matrix2f rotation(float degrees);

private:
    float m_elements[4];
};

// Scalar-Matrix multiplication
Matrix2f operator*(float f, const Matrix2f &m);
Matrix2f operator*(const Matrix2f &m, float f);

// Matrix-Vector multiplication
// 2x2 * 2x1 ==> 2x1
Vector2f operator*(const Matrix2f &m, const Vector2f &v);

// Matrix-Matrix multiplication
Matrix2f operator*(const Matrix2f &x, const Matrix2f &y);

class Matrix2f;
class Quat4f;
class Vector3f;

// 3x3 Matrix, stored in column major order (OpenGL style)
class Matrix3f
{
public:
    // Fill a 3x3 matrix with "fill", default to 0.
    Matrix3f(float fill = 0.f);
    Matrix3f(float m00, float m01, float m02,
             float m10, float m11, float m12,
             float m20, float m21, float m22);

    // setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2]
    // otherwise, sets the rows
    Matrix3f(const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, bool setColumns = true);

    Matrix3f(const Matrix3f &rm);            // copy constructor
    Matrix3f &operator=(const Matrix3f &rm); // assignment operator
    // no destructor necessary

    const float &operator()(int i, int j) const;
    float &operator()(int i, int j);

    Vector3f getRow(int i) const;
    void setRow(int i, const Vector3f &v);

    Vector3f getCol(int j) const;
    void setCol(int j, const Vector3f &v);

    // gets the 2x2 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    Matrix2f getSubmatrix2x2(int i0, int j0) const;

    // sets a 2x2 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    void setSubmatrix2x2(int i0, int j0, const Matrix2f &m);

    float determinant() const;
    Matrix3f inverse(bool *pbIsSingular = NULL, float epsilon = 0.f) const; // TODO: invert in place as well

    void transpose();
    Matrix3f transposed() const;

    // ---- Utility ----
    operator float *(); // automatic type conversion for GL
    void print();

    static float determinant3x3(float m00, float m01, float m02,
                                float m10, float m11, float m12,
                                float m20, float m21, float m22);

    static Matrix3f ones();
    static Matrix3f identity();
    static Matrix3f rotateX(float radians);
    static Matrix3f rotateY(float radians);
    static Matrix3f rotateZ(float radians);
    static Matrix3f scaling(float sx, float sy, float sz);
    static Matrix3f uniformScaling(float s);
    static Matrix3f rotation(const Vector3f &rDirection, float radians);

    // Returns the rotation matrix represented by a unit quaternion
    // if q is not normalized, it it normalized first
    static Matrix3f rotation(const Quat4f &rq);

private:
    float m_elements[9];
};

// Matrix-Vector multiplication
// 3x3 * 3x1 ==> 3x1
Vector3f operator*(const Matrix3f &m, const Vector3f &v);

// Matrix-Matrix multiplication
Matrix3f operator*(const Matrix3f &x, const Matrix3f &y);

class Matrix2f;
class Matrix3f;
class Quat4f;
class Vector3f;
class Vector4f;

// 4x4 Matrix, stored in column major order (OpenGL style)
class Matrix4f
{
public:
    // Fill a 4x4 matrix with "fill".  Default to 0.
    Matrix4f(float fill = 0.f);
    Matrix4f(float m00, float m01, float m02, float m03,
             float m10, float m11, float m12, float m13,
             float m20, float m21, float m22, float m23,
             float m30, float m31, float m32, float m33);

    // setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2 v3]
    // otherwise, sets the rows
    Matrix4f(const Vector4f &v0, const Vector4f &v1, const Vector4f &v2, const Vector4f &v3, bool setColumns = true);

    Matrix4f(const Matrix4f &rm);            // copy constructor
    Matrix4f &operator=(const Matrix4f &rm); // assignment operator
    Matrix4f &operator/=(float d);
    // no destructor necessary

    const float &operator()(int i, int j) const;
    float &operator()(int i, int j);

    Vector4f getRow(int i) const;
    void setRow(int i, const Vector4f &v);

    // get column j (mod 4)
    Vector4f getCol(int j) const;
    void setCol(int j, const Vector4f &v);

    // gets the 2x2 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    Matrix2f getSubmatrix2x2(int i0, int j0) const;

    // gets the 3x3 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    Matrix3f getSubmatrix3x3(int i0, int j0) const;

    // sets a 2x2 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    void setSubmatrix2x2(int i0, int j0, const Matrix2f &m);

    // sets a 3x3 submatrix of this matrix to m
    // starting with upper left corner at (i0, j0)
    void setSubmatrix3x3(int i0, int j0, const Matrix3f &m);

    float determinant() const;
    Matrix4f inverse(bool *pbIsSingular = NULL, float epsilon = 0.f) const;

    void transpose();
    Matrix4f transposed() const;

    // ---- Utility ----
    operator float *();             // automatic type conversion for GL
    operator const float *() const; // automatic type conversion for GL

    void print();

    static Matrix4f ones();
    static Matrix4f identity();
    static Matrix4f translation(float x, float y, float z);
    static Matrix4f translation(const Vector3f &rTranslation);
    static Matrix4f rotateX(float radians);
    static Matrix4f rotateY(float radians);
    static Matrix4f rotateZ(float radians);
    static Matrix4f rotation(const Vector3f &rDirection, float radians);
    static Matrix4f scaling(float sx, float sy, float sz);
    static Matrix4f uniformScaling(float s);
    static Matrix4f lookAt(const Vector3f &eye, const Vector3f &center, const Vector3f &up);
    static Matrix4f orthographicProjection(float width, float height, float zNear, float zFar, bool directX);
    static Matrix4f orthographicProjection(float left, float right, float bottom, float top, float zNear, float zFar, bool directX);
    static Matrix4f perspectiveProjection(float fLeft, float fRight, float fBottom, float fTop, float fZNear, float fZFar, bool directX);
    static Matrix4f perspectiveProjection(float fovYRadians, float aspect, float zNear, float zFar, bool directX);
    static Matrix4f infinitePerspectiveProjection(float fLeft, float fRight, float fBottom, float fTop, float fZNear, bool directX);

    // Returns the rotation matrix represented by a quaternion
    // uses a normalized version of q
    static Matrix4f rotation(const Quat4f &q);

    // returns an orthogonal matrix that's a uniformly distributed rotation
    // given u[i] is a uniformly distributed random number in [0,1]
    static Matrix4f randomRotation(float u0, float u1, float u2);

private:
    float m_elements[16];
};

// Matrix-Vector multiplication
// 4x4 * 4x1 ==> 4x1
Vector4f operator*(const Matrix4f &m, const Vector4f &v);

// Matrix-Matrix multiplication
Matrix4f operator*(const Matrix4f &x, const Matrix4f &y);

class Vector3f;
class Vector4f;

class Quat4f
{
public:
    static const Quat4f ZERO;
    static const Quat4f IDENTITY;

    Quat4f();

    // q = w + x * i + y * j + z * k
    Quat4f(float w, float x, float y, float z);

    Quat4f(const Quat4f &rq);            // copy constructor
    Quat4f &operator=(const Quat4f &rq); // assignment operator
    // no destructor necessary

    // returns a quaternion with 0 real part
    Quat4f(const Vector3f &v);

    // copies the components of a Vector4f directly into this quaternion
    Quat4f(const Vector4f &v);

    // returns the ith element
    const float &operator[](int i) const;
    float &operator[](int i);

    float w() const;
    float x() const;
    float y() const;
    float z() const;
    Vector3f xyz() const;
    Vector4f wxyz() const;

    float abs() const;
    float absSquared() const;
    void normalize();
    Quat4f normalized() const;

    void conjugate();
    Quat4f conjugated() const;

    void invert();
    Quat4f inverse() const;

    // log and exponential maps
    Quat4f log() const;
    Quat4f exp() const;

    // returns unit vector for rotation and radians about the unit vector
    Vector3f getAxisAngle(float *radiansOut);

    // sets this quaternion to be a rotation of fRadians about v = < fx, fy, fz >, v need not necessarily be unit length
    void setAxisAngle(float radians, const Vector3f &axis);

    // ---- Utility ----
    void print();

    // quaternion dot product (a la vector)
    static float dot(const Quat4f &q0, const Quat4f &q1);

    // linear (stupid) interpolation
    static Quat4f lerp(const Quat4f &q0, const Quat4f &q1, float alpha);

    // spherical linear interpolation
    static Quat4f slerp(const Quat4f &a, const Quat4f &b, float t, bool allowFlip = true);

    // spherical quadratic interoplation between a and b at point t
    // given quaternion tangents tanA and tanB (can be computed using squadTangent)
    static Quat4f squad(const Quat4f &a, const Quat4f &tanA, const Quat4f &tanB, const Quat4f &b, float t);

    static Quat4f cubicInterpolate(const Quat4f &q0, const Quat4f &q1, const Quat4f &q2, const Quat4f &q3, float t);

    // Log-difference between a and b, used for squadTangent
    // returns log( a^-1 b )
    static Quat4f logDifference(const Quat4f &a, const Quat4f &b);

    // Computes a tangent at center, defined by the before and after quaternions
    // Useful for squad()
    static Quat4f squadTangent(const Quat4f &before, const Quat4f &center, const Quat4f &after);

    static Quat4f fromRotationMatrix(const Matrix3f &m);

    static Quat4f fromRotatedBasis(const Vector3f &x, const Vector3f &y, const Vector3f &z);

    // returns a unit quaternion that's a uniformly distributed rotation
    // given u[i] is a uniformly distributed random number in [0,1]
    // taken from Graphics Gems II
    static Quat4f randomRotation(float u0, float u1, float u2);

private:
    float m_elements[4];
};

Quat4f operator+(const Quat4f &q0, const Quat4f &q1);
Quat4f operator-(const Quat4f &q0, const Quat4f &q1);
Quat4f operator*(const Quat4f &q0, const Quat4f &q1);
Quat4f operator*(float f, const Quat4f &q);
Quat4f operator*(const Quat4f &q, float f);

class Vector3f;

class Vector2f
{
public:
    static const Vector2f ZERO;
    static const Vector2f UP;
    static const Vector2f RIGHT;

    Vector2f(float f = 0.f);
    Vector2f(float x, float y);

    // copy constructors
    Vector2f(const Vector2f &rv);

    // assignment operators
    Vector2f &operator=(const Vector2f &rv);

    // no destructor necessary

    // returns the ith element
    const float &operator[](int i) const;
    float &operator[](int i);

    float &x();
    float &y();

    float x() const;
    float y() const;

    Vector2f xy() const;
    Vector2f yx() const;
    Vector2f xx() const;
    Vector2f yy() const;

    // returns ( -y, x )
    Vector2f normal() const;

    float abs() const;
    float absSquared() const;
    void normalize();
    Vector2f normalized() const;

    void negate();

    // ---- Utility ----
    operator const float *() const; // automatic type conversion for OpenGL
    operator float *();             // automatic type conversion for OpenGL
    void print() const;

    Vector2f &operator+=(const Vector2f &v);
    Vector2f &operator-=(const Vector2f &v);
    Vector2f &operator*=(float f);

    static float dot(const Vector2f &v0, const Vector2f &v1);

    static Vector3f cross(const Vector2f &v0, const Vector2f &v1);

    // returns v0 * ( 1 - alpha ) * v1 * alpha
    static Vector2f lerp(const Vector2f &v0, const Vector2f &v1, float alpha);

private:
    float m_elements[2];
};

// component-wise operators
Vector2f operator+(const Vector2f &v0, const Vector2f &v1);
Vector2f operator-(const Vector2f &v0, const Vector2f &v1);
Vector2f operator*(const Vector2f &v0, const Vector2f &v1);
Vector2f operator/(const Vector2f &v0, const Vector2f &v1);

// unary negation
Vector2f operator-(const Vector2f &v);

// multiply and divide by scalar
Vector2f operator*(float f, const Vector2f &v);
Vector2f operator*(const Vector2f &v, float f);
Vector2f operator/(const Vector2f &v, float f);

bool operator==(const Vector2f &v0, const Vector2f &v1);
bool operator!=(const Vector2f &v0, const Vector2f &v1);

class Vector2f;

class Vector3f
{
public:
    static const Vector3f ZERO;
    static const Vector3f UP;
    static const Vector3f RIGHT;
    static const Vector3f FORWARD;

    Vector3f(float f = 0.f);
    Vector3f(float x, float y, float z);

    Vector3f(const Vector2f &xy, float z);
    Vector3f(float x, const Vector2f &yz);

    // copy constructors
    Vector3f(const Vector3f &rv);

    // assignment operators
    Vector3f &operator=(const Vector3f &rv);

    // no destructor necessary

    // returns the ith element
    const float &operator[](int i) const;
    float &operator[](int i);

    float &x();
    float &y();
    float &z();

    float x() const;
    float y() const;
    float z() const;

    Vector2f xy() const;
    Vector2f xz() const;
    Vector2f yz() const;

    Vector3f xyz() const;
    Vector3f yzx() const;
    Vector3f zxy() const;

    float length() const;
    float squaredLength() const;

    void normalize();
    Vector3f normalized() const;

    Vector2f homogenized() const;

    void negate();

    // ---- Utility ----
    operator const float *() const; // automatic type conversion for OpenGL
    operator float *();             // automatic type conversion for OpenGL
    void print() const;

    Vector3f &operator+=(const Vector3f &v);
    Vector3f &operator-=(const Vector3f &v);
    Vector3f &operator*=(float f);
    Vector3f operator%(Vector3f &b);

    static float dot(const Vector3f &v0, const Vector3f &v1);
    static Vector3f cross(const Vector3f &v0, const Vector3f &v1);

    // computes the linear interpolation between v0 and v1 by alpha \in [0,1]
    // returns v0 * ( 1 - alpha ) * v1 * alpha
    static Vector3f lerp(const Vector3f &v0, const Vector3f &v1, float alpha);

    // computes the cubic catmull-rom interpolation between p0, p1, p2, p3
    // by t \in [0,1].  Guarantees that at t = 0, the result is p0 and
    // at p1, the result is p2.
    static Vector3f cubicInterpolate(const Vector3f &p0, const Vector3f &p1, const Vector3f &p2, const Vector3f &p3, float t);

    float max() const
    {
        return m_elements[0] > m_elements[1] && m_elements[0] > m_elements[2]
                   ? m_elements[0]
               : m_elements[1] > m_elements[2] ? m_elements[1]
                                               : m_elements[2];
    }

private:
    float m_elements[3];
};

// component-wise operators
Vector3f operator+(const Vector3f &v0, const Vector3f &v1);
Vector3f operator-(const Vector3f &v0, const Vector3f &v1);
Vector3f operator*(const Vector3f &v0, const Vector3f &v1);
Vector3f operator/(const Vector3f &v0, const Vector3f &v1);

// unary negation
Vector3f operator-(const Vector3f &v);

// multiply and divide by scalar
Vector3f operator*(float f, const Vector3f &v);
Vector3f operator*(const Vector3f &v, float f);
Vector3f operator/(const Vector3f &v, float f);

bool operator==(const Vector3f &v0, const Vector3f &v1);
bool operator!=(const Vector3f &v0, const Vector3f &v1);

class Vector4f
{
public:
    Vector4f(float f = 0.f);
    Vector4f(float fx, float fy, float fz, float fw);
    Vector4f(float buffer[4]);

    Vector4f(const Vector2f &xy, float z, float w);
    Vector4f(float x, const Vector2f &yz, float w);
    Vector4f(float x, float y, const Vector2f &zw);
    Vector4f(const Vector2f &xy, const Vector2f &zw);

    Vector4f(const Vector3f &xyz, float w);
    Vector4f(float x, const Vector3f &yzw);

    // copy constructors
    Vector4f(const Vector4f &rv);

    // assignment operators
    Vector4f &operator=(const Vector4f &rv);

    // no destructor necessary

    // returns the ith element
    const float &operator[](int i) const;
    float &operator[](int i);

    float &x();
    float &y();
    float &z();
    float &w();

    float x() const;
    float y() const;
    float z() const;
    float w() const;

    Vector2f xy() const;
    Vector2f yz() const;
    Vector2f zw() const;
    Vector2f wx() const;

    Vector3f xyz() const;
    Vector3f yzw() const;
    Vector3f zwx() const;
    Vector3f wxy() const;

    Vector3f xyw() const;
    Vector3f yzx() const;
    Vector3f zwy() const;
    Vector3f wxz() const;

    float abs() const;
    float absSquared() const;
    void normalize();
    Vector4f normalized() const;

    // if v.z != 0, v = v / v.w
    void homogenize();
    Vector4f homogenized() const;

    void negate();

    // ---- Utility ----
    operator const float *() const; // automatic type conversion for OpenGL
    operator float *();             // automatic type conversion for OpenG
    void print() const;

    static float dot(const Vector4f &v0, const Vector4f &v1);
    static Vector4f lerp(const Vector4f &v0, const Vector4f &v1, float alpha);

private:
    float m_elements[4];
};

// component-wise operators
Vector4f operator+(const Vector4f &v0, const Vector4f &v1);
Vector4f operator-(const Vector4f &v0, const Vector4f &v1);
Vector4f operator*(const Vector4f &v0, const Vector4f &v1);
Vector4f operator/(const Vector4f &v0, const Vector4f &v1);

// unary negation
Vector4f operator-(const Vector4f &v);

// multiply and divide by scalar
Vector4f operator*(float f, const Vector4f &v);
Vector4f operator*(const Vector4f &v, float f);
Vector4f operator/(const Vector4f &v, float f);

bool operator==(const Vector4f &v0, const Vector4f &v1);
bool operator!=(const Vector4f &v0, const Vector4f &v1);

#endif // VECMATH
