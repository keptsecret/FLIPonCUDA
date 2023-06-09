#ifndef FOC_VECMATH_H
#define FOC_VECMATH_H

#include "foc.h"

namespace foc {

/* -------------------------------------------------------------------------------
 * Vectors
 */

template <typename T>
class Vector2 {
public:
	Vector2() :
			x(0), y(0) {}

	FOC_CPU_GPU
	Vector2(T _x, T _y) :
			x(_x), y(_y) {
	}

	FOC_CPU_GPU explicit Vector2(const Point2<T>& p);

	FOC_CPU_GPU explicit Vector2(const Point3<T>& p);

	FOC_CPU_GPU
	Vector2<T> operator+(const Vector2<T>& v) const {
		return Vector2(x + v.x, y + v.y);
	}

	FOC_CPU_GPU
	Vector2<T>& operator+=(const Vector2<T>& v) const {
		x += v.x;
		y += v.y;
		return *this;
	}

	FOC_CPU_GPU
	Vector2<T> operator-(const Vector2<T>& v) const {
		return Vector2(x - v.x, y - v.y);
	}

	FOC_CPU_GPU
	Vector2<T> operator-() const {
		return Vector2(-x, -y);
	}

	FOC_CPU_GPU
	Vector2<T>& operator-=(const Vector2<T>& v) const {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	FOC_CPU_GPU
	Vector2<T> operator*(T s) const {
		return Vector2(s * x, s * y);
	}

	FOC_CPU_GPU
	Vector2<T>& operator*=(T s) const {
		x *= s;
		y *= s;
		return *this;
	}

	FOC_CPU_GPU
	Vector2<T> operator/(T s) const {
		float inv = 1.0f / s;
		return Vector2(x * inv, y * inv);
	}

	FOC_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		float inv = 1.0f / s;
		x *= inv;
		y *= inv;
		return *this;
	}

	FOC_CPU_GPU
	bool operator==(const Vector2<T>& v) const {
		return x == v.x && y == v.y;
	}

	FOC_CPU_GPU
	bool operator!=(const Vector2<T>& v) const {
		return x != v.x || y != v.y;
	}

	FOC_CPU_GPU
	T operator[](int i) const {
		if (i == 0) {
			return x;
		}
		return y;
	}

	FOC_CPU_GPU
	T& operator[](int i) {
		if (i == 0) {
			return x;
		}
		return y;
	}

	FOC_CPU_GPU
	float lengthSquared() const { return x * x + y * y; }

	FOC_CPU_GPU
	float length() const { return std::sqrt(lengthSquared()); }

	FOC_CPU_GPU
	friend std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
		os << "[ " << v.x << ", " << v.y << " ]";
		return os;
	}

	T x, y;
};

template <typename T>
class Vector3 {
public:
	Vector3() :
			x(0), y(0), z(0) {}

	FOC_CPU_GPU
	Vector3(T _x, T _y, T _z) :
			x(_x), y(_y), z(_z) {
	}

	FOC_CPU_GPU explicit Vector3(const Point3<T>& p);

	FOC_CPU_GPU
	Vector3<T> operator+(const Vector3<T>& v) const {
		return Vector3(x + v.x, y + v.y, z + v.z);
	}

	FOC_CPU_GPU
	Vector3<T>& operator+=(const Vector3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	FOC_CPU_GPU
	Vector3<T> operator-(const Vector3<T>& v) const {
		return Vector3(x - v.x, y - v.y, z - v.z);
	}

	FOC_CPU_GPU
	Vector3<T> operator-() const {
		return Vector3(-x, -y, -z);
	}

	FOC_CPU_GPU
	Vector3<T>& operator-=(const Vector3<T>& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	FOC_CPU_GPU
	Vector3<T> operator*(T s) const {
		return Vector3(s * x, s * y, s * z);
	}

	FOC_CPU_GPU
	Vector3<T>& operator*=(T s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	FOC_CPU_GPU
	Vector3<T> operator/(T s) const {
		float inv = 1.0f / s;
		return Vector3(x * inv, y * inv, z * inv);
	}

	FOC_CPU_GPU
	Vector3<T>& operator/=(T s) {
		float inv = 1.0f / s;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	FOC_CPU_GPU
	bool operator==(const Vector3<T>& v) const {
		return x == v.x && y == v.y && z == v.z;
	}

	FOC_CPU_GPU
	bool operator!=(const Vector3<T>& v) const {
		return x != v.x || y != v.y || z != v.z;
	}

	FOC_CPU_GPU
	T operator[](int i) const {
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	FOC_CPU_GPU
	T& operator[](int i) {
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	FOC_CPU_GPU
	float lengthSquared() const { return x * x + y * y + z * z; }

	FOC_CPU_GPU
	float length() const { return std::sqrt(lengthSquared()); }

	T x, y, z;
};

/* -------------------------------------------------------------------------------
 * Point
 */

template <typename T>
class Point3 {
public:
	Point3<T>() :
			x(0), y(0), z(0) {}

	FOC_CPU_GPU
	Point3<T>(T _x, T _y, T _z) :
			x(_x), y(_y), z(_z) {
	}

	template <typename U>
	explicit Point3(const Point3<U>& p) :
			x((T)p.x), y((T)p.y), z((T)p.z) {
	}

	template <typename U>
	FOC_CPU_GPU explicit operator Vector3<U>() const {
		return Vector3<U>(x, y, z);
	}

	FOC_CPU_GPU
	Point3<T> operator+(const Vector3<T>& v) const {
		return Point3<T>(x + v.x, y + v.y, z + v.z);
	}

	FOC_CPU_GPU
	Point3<T>& operator+=(const Vector3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	FOC_CPU_GPU
	Point3<T> operator+(const Point3<T>& v) const {
		return Point3<T>(x + v.x, y + v.y, z + v.z);
	}

	FOC_CPU_GPU
	Point3<T>& operator+=(const Point3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	FOC_CPU_GPU
	Vector3<T> operator-(const Point3<T>& p) const {
		return Vector3<T>(x - p.x, y - p.y, z - p.z);
	}

	FOC_CPU_GPU
	Point3<T> operator-(const Vector3<T>& v) const {
		return Point3<T>(x - v.x, y - v.y, z - v.z);
	}

	FOC_CPU_GPU
	Point3<T>& operator-=(const Vector3<T>& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	template <typename U>
	FOC_CPU_GPU Point3<T> operator*(U f) const {
		return Point3<T>(f * x, f * y, f * z);
	}

	template <typename U>
	FOC_CPU_GPU Point3<T>& operator*=(U f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	template <typename U>
	FOC_CPU_GPU Point3<T> operator/(U f) const {
		float inv = 1.0f / f;
		return Point3<T>(inv * x, inv * y, inv * z);
	}

	template <typename U>
	FOC_CPU_GPU Point3<T>& operator/=(U f) {
		float inv = 1.0f / f;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	FOC_CPU_GPU
	T operator[](int i) const {
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	FOC_CPU_GPU
	T& operator[](int i) {
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	FOC_CPU_GPU
	bool operator==(const Point3<T>& p) const {
		return x == p.x && y == p.y && z == p.z;
	}

	FOC_CPU_GPU
	bool operator!=(const Point3<T>& p) const {
		return x != p.x || y != p.y || z != p.z;
	}

	T x, y, z;
};

template <typename T>
class Point2 {
public:
	Point2<T>() :
			x(0), y(0) {}

	FOC_CPU_GPU
	Point2<T>(T _x, T _y) :
			x(_x), y(_y) {
	}

	FOC_CPU_GPU
	explicit Point2(const Point3<T>& p) :
			x(p.x), y(p.y) {
	}

	template <typename U>
	explicit Point2(const Point2<U>& p) {
		x = (T)p.x;
		y = (T)p.y;
	}

	template <typename U>
	explicit Point2(const Vector2<U>& p) {
		x = (T)p.x;
		y = (T)p.y;
	}

	template <typename U>
	explicit operator Vector2<U>() const {
		return Vector2<U>(x, y);
	}

	FOC_CPU_GPU
	Point2<T> operator+(const Vector2<T>& v) const {
		return Point2<T>(x + v.x, y + v.y);
	}

	FOC_CPU_GPU
	Point2<T>& operator+=(const Vector2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	FOC_CPU_GPU
	Point2<T> operator+(const Point2<T>& v) const {
		return Point2<T>(x + v.x, y + v.y);
	}

	FOC_CPU_GPU
	Point2<T>& operator+=(const Point2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	FOC_CPU_GPU
	Vector2<T> operator-(const Point2<T>& p) const {
		return Vector2<T>(x - p.x, y - p.y);
	}

	FOC_CPU_GPU
	Point2<T> operator-(const Vector2<T>& v) const {
		return Point2<T>(x - v.x, y - v.y);
	}

	FOC_CPU_GPU
	Point2<T>& operator-=(const Vector2<T>& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	template <typename U>
	FOC_CPU_GPU Point2<T> operator*(U f) const {
		return Point2<T>(f * x, f * y);
	}

	template <typename U>
	FOC_CPU_GPU Point2<T>& operator*=(U f) {
		x *= f;
		y *= f;
		return *this;
	}

	FOC_CPU_GPU
	Vector2<T> operator/(T s) const {
		float inv = 1 / s;
		return Vector2(x * inv, y * inv);
	}

	FOC_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		float inv = 1 / s;
		x *= inv;
		y *= inv;
		return *this;
	}

	FOC_CPU_GPU
	bool operator==(const Point2<T>& v) const {
		return x == v.x && y == v.y;
	}

	FOC_CPU_GPU
	bool operator!=(const Point2<T>& v) const {
		return x != v.x || y != v.y;
	}

	FOC_CPU_GPU
	T operator[](int i) const {
		if (i == 0) {
			return x;
		}
		return y;
	}

	FOC_CPU_GPU
	T& operator[](int i) {
		if (i == 0) {
			return x;
		}
		return y;
	}

	T x, y;
};

/* -------------------------------------------------------------------------------
 * Bounds (AABBs)
 */

template <typename T>
class Bounds3 {
public:
	Bounds3() {
		T min_num = std::numeric_limits<T>::lowest();
		T max_num = std::numeric_limits<T>::max();
		p_min = Point3<T>(max_num, max_num, max_num);
		p_max = Point3<T>(min_num, min_num, min_num);
	}

	FOC_CPU_GPU
	explicit Bounds3(const Point3<T>& p) :
			p_min(p), p_max(p) {}

	FOC_CPU_GPU
	Bounds3(const Point3<T>& p1, const Point3<T>& p2) :
			p_min(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)), p_max(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)) {}

	FOC_CPU_GPU
	inline const Point3<T>& operator[](int i) const {
		return (i == 0) ? p_min : p_max;
	}

	FOC_CPU_GPU
	inline Point3<T>& operator[](int i) {
		return (i == 0) ? p_min : p_max;
	}

	bool operator==(const Bounds3<T>& b) const {
		return b.p_min == p_min && b.p_max == p_max;
	}

	bool operator!=(const Bounds3<T>& b) const {
		return b.p_min != p_min || b.p_max != p_max;
	}

	FOC_CPU_GPU
	Point3<T> corner(int n_corner) const {
		return Point3<T>((*this)[(n_corner & 1)].x,
				(*this)[(n_corner & 2) ? 1 : 0].y,
				(*this)[(n_corner & 4) ? 1 : 0].z);
	}

	FOC_CPU_GPU
	Vector3<T> diagonal() const { return p_max - p_min; }

	FOC_CPU_GPU
	T surfaceArea() const {
		Vector3<T> d = diagonal();
		return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
	}

	FOC_CPU_GPU
	T volume() const {
		Vector3<T> d = diagonal();
		return d.x * d.y * d.z;
	}

	FOC_CPU_GPU
	int maximumAxis() const {
		Vector3<T> d = diagonal();
		if (d.x > d.y && d.x > d.z) {
			return 0;
		} else if (d.y > d.z) {
			return 1;
		} else {
			return 2;
		}
	}

	FOC_CPU_GPU
	Vector3<T> offset(const Point3<T>& p) const {
		Vector3<T> o = p - p_min;
		if (p_max.x > p_min.x) {
			o.x /= p_max.x - p_min.x;
		}
		if (p_max.y > p_min.y) {
			o.y /= p_max.y - p_min.y;
		}
		if (p_max.z > p_min.z) {
			o.z /= p_max.z - p_min.z;
		}
		return o;
	}

	FOC_CPU_GPU
	void boundingSphere(Point3<T>* center, float* radius) const {
		*center = (p_min + p_max) / 2;
		*radius = bInside(*center, *this) ? distance(*center, p_max) : 0;
	}

	FOC_CPU_GPU
	bool intersectP(Point3f o, Vector3f d, float t_max = Infinity, float* hitt0 = nullptr, float* hitt1 = nullptr) const;

	FOC_CPU_GPU
	bool intersectP(Point3f o, Vector3f d, float t_max, Vector3f& inv_dir, const int dir_is_neg[3]) const;

public:
	Point3<T> p_min, p_max;
};

/*---------------------------------------------------------------------------------*/
/*
 * Inline geometry functions
 */

// Vector2
template <typename T>
FOC_CPU_GPU Vector2<T>::Vector2(const Point2<T>& p) :
		x(p.x), y(p.y) {
}

template <typename T>
FOC_CPU_GPU Vector2<T>::Vector2(const Point3<T>& p) :
		x(p.x), y(p.y) {
}

template <typename T>
FOC_CPU_GPU inline Vector2<T> operator*(T s, const Vector2<T>& v) {
	return v * s;
}

template <typename T>
FOC_CPU_GPU inline Vector2<T> abs(const Vector2<T>& v) {
	return Vector2<T>(std::abs(v.x), std::abs(v.y));
}

template <typename T>
FOC_CPU_GPU inline T dot(const Vector2<T>& v1, const Vector2<T>& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
FOC_CPU_GPU inline T absDot(const Vector2<T>& v1, const Vector2<T>& v2) {
	return std::abs(dot(v1, v2));
}

template <typename T>
FOC_CPU_GPU inline Vector2<T> normalize(const Vector2<T>& v) {
	return v / v.length();
}

template <typename T>
FOC_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

// Vector3
template <typename T>
FOC_CPU_GPU Vector3<T>::Vector3(const Point3<T>& p) :
		x(p.x), y(p.y), z(p.z) {
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> operator*(T s, const Vector3<T>& v) {
	return v * s;
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> abs(const Vector3<T>& v) {
	return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
FOC_CPU_GPU inline T dot(const Vector3<T>& v1, const Vector3<T>& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
FOC_CPU_GPU inline T absDot(const Vector3<T>& v1, const Vector3<T>& v2) {
	return std::abs(dot(v1, v2));
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> cross(const Vector3<T>& v1, const Vector3<T>& v2) {
	double v1x = v1.x, v1y = v1.y, v1z = v1.z;
	double v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector3<T>(v1y * v2z - v1z * v2y,
			v1z * v2x - v1x * v2z,
			v1x * v2y - v1y * v2z);
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> normalize(const Vector3<T>& v) {
	return v / v.length();
}

template <typename T>
FOC_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

template <typename T>
FOC_CPU_GPU inline T minComponent(const Vector3<T>& v) {
	return std::min(v.x, std::min(v.y, v.z));
}

template <typename T>
FOC_CPU_GPU inline T maxComponent(const Vector3<T>& v) {
	return std::max(v.x, std::max(v.y, v.z));
}

template <typename T>
FOC_CPU_GPU inline int maxDimension(const Vector3<T>& v) {
	return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> min(const Vector3<T>& v1, const Vector3<T>& v2) {
	return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> max(const Vector3<T>& v1, const Vector3<T>& v2) {
	return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

template <typename T>
FOC_CPU_GPU inline Vector3<T> permute(const Vector3<T>& v, int x, int y, int z) {
	return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
FOC_CPU_GPU inline void coordinateSystem(const Vector3<T>& v1, Vector3<T>* v2, Vector3<T>* v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		*v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	} else {
		*v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	}
	*v3 = cross(v1, *v2);
}

// Point2
template <typename T>
FOC_CPU_GPU inline float distance(const Point2<T>& p1, const Point2<T>& p2) {
	return (p1 - p2).length();
}

template <typename T>
FOC_CPU_GPU inline float distanceSquared(const Point2<T>& p1, const Point2<T>& p2) {
	return (p1 - p2).lengthSquared();
}

template <typename T, typename U>
FOC_CPU_GPU inline Point2<T> operator*(U f, const Point2<T>& p) {
	return p * f;
}

template <typename T>
FOC_CPU_GPU inline Point2<T> lerp(float t, const Point2<T>& p1, const Point2<T>& p2) {
	return (1 - t) * p1 + t * p2;
}

template <typename T>
FOC_CPU_GPU inline Point2<T> min(const Point2<T>& v1, const Point2<T>& v2) {
	return Point2<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y));
}

template <typename T>
FOC_CPU_GPU inline Point2<T> max(const Point2<T>& v1, const Point2<T>& v2) {
	return Point2<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y));
}

template <typename T>
FOC_CPU_GPU inline Point2<T> floor(const Point2<T>& p) {
	return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
FOC_CPU_GPU inline Point2<T> ceil(const Point2<T>& p) {
	return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
FOC_CPU_GPU inline Point2<T> abs(const Point2<T>& p) {
	return Point2<T>(std::abs(p.x), std::abs(p.y));
}

template <typename T>
FOC_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point2<T>& v) {
	os << "[ " << v.x << ", " << v.y << " ]";
	return os;
}

// Point3
template <typename T>
FOC_CPU_GPU inline float distance(const Point3<T>& p1, const Point3<T>& p2) {
	return (p1 - p2).length();
}

template <typename T>
FOC_CPU_GPU inline float distanceSquared(const Point3<T>& p1, const Point3<T>& p2) {
	return (p1 - p2).lengthSquared();
}

template <typename T, typename U>
FOC_CPU_GPU inline Point3<T> operator*(U f, const Point3<T>& p) {
	return p * f;
}

template <typename T>
FOC_CPU_GPU inline Point3<T> lerp(float t, const Point3<T>& p1, const Point3<T>& p2) {
	return (1 - t) * p1 + t * p2;
}

template <typename T>
FOC_CPU_GPU inline Point3<T> min(const Point3<T>& v1, const Point3<T>& v2) {
	return Point3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

template <typename T>
FOC_CPU_GPU inline Point3<T> max(const Point3<T>& v1, const Point3<T>& v2) {
	return Point3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

template <typename T>
FOC_CPU_GPU inline Point3<T> floor(const Point3<T>& p) {
	return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
FOC_CPU_GPU inline Point3<T> ceil(const Point3<T>& p) {
	return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
FOC_CPU_GPU inline Point3<T> abs(const Point3<T>& p) {
	return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
FOC_CPU_GPU inline Point3<T> permute(const Point3<T>& p, int x, int y, int z) {
	return Point3<T>(p[x], p[y], p[z]);
}

template <typename T>
FOC_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

// Bounds3
template <typename T>
FOC_CPU_GPU inline Bounds3<T> bUnion(const Bounds3<T>& b, const Point3<T>& p) {
	return Bounds3<T>(min(b.p_min, p), max(b.p_max, p));
}

template <typename T>
FOC_CPU_GPU inline Bounds3<T> bUnion(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	return Bounds3<T>(min(b1.p_min, b2.p_min), max(b1.p_max, b2.p_max));
}

template <typename T>
FOC_CPU_GPU inline Bounds3<T> bIntersect(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	Bounds3<T> b;
	b.p_min = max(b1.p_min, b2.p_min);
	b.p_max = min(b1.p_max, b2.p_max);
	return b;
}

template <typename T>
FOC_CPU_GPU inline bool bOverlaps(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	bool x = (b1.p_max.x >= b2.p_min.x) && (b1.p_min.x <= b2.p_max.x);
	bool y = (b1.p_max.y >= b2.p_min.y) && (b1.p_min.y <= b2.p_max.y);
	bool z = (b1.p_max.z >= b2.p_min.z) && (b1.p_min.z <= b2.p_max.z);
	return x && y && z;
}

template <typename T>
FOC_CPU_GPU inline bool bInside(const Point3<T>& p, const Bounds3<T>& b) {
	return p.x >= b.p_min.x && p.x <= b.p_max.x &&
			p.y >= b.p_min.y && p.y <= b.p_max.y &&
			p.z >= b.p_min.z && p.z <= b.p_max.z;
}

template <typename T>
FOC_CPU_GPU inline bool bInsideExclusive(const Point3<T>& p, const Bounds3<T>& b) {
	return p.x >= b.p_min.x && p.x < b.p_max.x &&
			p.y >= b.p_min.y && p.y < b.p_max.y &&
			p.z >= b.p_min.z && p.z < b.p_max.z;
}

template <typename T, typename U>
FOC_CPU_GPU inline bool bExpand(const Bounds3<T>& b, U delta) {
	return Bounds3<T>(b.p_min - Vector3<T>(delta, delta, delta), b.p_max + Vector3<T>(delta, delta, delta));
}

template <typename T, typename U>
FOC_CPU_GPU inline Point3<T> bGetNearestPointInsideBounds(const Point3<T>& p, U eps, const Bounds3<T>& b) {
	Point3<T> pos = p;

	if (pos.x >= b.p_min.x + eps && pos.y >= b.p_min.y + eps && pos.z >= b.p_min.z + eps &&
			pos.x <= b.p_max.x + eps && pos.y <= b.p_max.y + eps && pos.z <= b.p_max.z + eps) {
		return pos;
	}

	pos.x = fmax(pos.x, b.p_min.x + eps);
	pos.y = fmax(pos.y, b.p_min.y + eps);
	pos.z = fmax(pos.z, b.p_min.z + eps);

	pos.x = fmin(pos.x, b.p_max.x - eps);
	pos.y = fmin(pos.y, b.p_max.y - eps);
	pos.z = fmin(pos.z, b.p_max.z - eps);

	return pos;
}

template <typename T>
FOC_CPU_GPU inline bool Bounds3<T>::intersectP(Point3f o, Vector3f d, float t_max, float* hitt0, float* hitt1) const {
	float t0 = 0, t1 = t_max;
	for (int i = 0; i < 3; i++) {
		float inv_dir = 1 / d[i];
		float t_near = (p_min[i] - o[i]) * inv_dir;
		float t_far = (p_max[i] - o[i]) * inv_dir;

		if (t_near > t_far) {
			float tmp = t_near;
			t_near = t_far;
			t_far = tmp;
		}
		t_far *= 1 + 2 * gamma(3);
		t0 = t_near > t0 ? t_near : t0;
		t1 = t_far < t1 ? t_far : t1;
		if (t0 > t1) {
			return false;
		}
	}
	if (hitt0) {
		*hitt0 = t0;
	}
	if (hitt1) {
		*hitt1 = t1;
	}
	return true;
}

template <typename T>
FOC_CPU_GPU inline bool Bounds3<T>::intersectP(Point3f o, Vector3f d, float ray_t_max, Vector3f& inv_dir, const int dir_is_neg[3]) const {
	const Bounds3f& bounds = *this;
	float t_min = (bounds[dir_is_neg[0]].x - o.x) * inv_dir.x;
	float t_max = (bounds[1 - dir_is_neg[0]].x - o.x) * inv_dir.x;
	float ty_min = (bounds[dir_is_neg[1]].y - o.y) * inv_dir.y;
	float ty_max = (bounds[1 - dir_is_neg[1]].y - o.y) * inv_dir.y;

	t_max *= 1 + 2 * gamma(3);
	ty_max *= 1 + 2 * gamma(3);

	if (t_min > ty_max || ty_min > t_max) {
		return false;
	}
	if (ty_min > t_min) {
		t_min = ty_min;
	}
	if (ty_max < t_max) {
		t_max = ty_min;
	}

	float tz_min = (bounds[dir_is_neg[2]].z - o.z) * inv_dir.z;
	float tz_max = (bounds[1 - dir_is_neg[2]].z - o.z) * inv_dir.z;

	tz_max *= 1 + 2 * gamma(3);

	if (t_min > tz_max || tz_min > t_max) {
		return false;
	}
	if (tz_min > t_min) {
		t_min = tz_min;
	}
	if (tz_max < t_max) {
		t_max = tz_min;
	}

	return (t_min < ray_t_max) && (t_max > 0);
}

} // namespace foc

#endif // FOC_VECMATH_H