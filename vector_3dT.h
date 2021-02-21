//
//  vector_3dT.h
//  linear_algebra
//

#ifndef __vector3d_T_H__
#define __vector3d_T_H__

#include <iostream>
#include <cstring>
#include <cmath>
#include <initializer_list>

template <typename T> class vector3d;
template <typename T> std::ostream& operator<<(std::ostream& os, const vector3d<T>& v);

typedef vector3d<double> vector3D;
typedef vector3d<float>  vector3F;
typedef vector3d<int>    vector3I;
typedef vector3d<long>   vector3L;
typedef vector3D vec3;


template <typename T>
class vector3d {  // class that serves as both vectors and points
public:
  vector3d();
  vector3d(const std::string& name, int dims);
  vector3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//---------------------------------------------------------------------
  T  operator[](int i) const;
  T& operator[](int i);
//---------------------------------------------------------------------
  void name(const std::string& name);
  std::string name() const;
//---------------------------------------------------------------------
  vector3d<T>& operator+=(const vector3d<T>& v);
  vector3d<T>& operator-=(const vector3d<T>& v);
//---------------------------------------------------------------------
  vector3d<T>& operator+=(T k);
  vector3d<T>& operator-=(T k);
  vector3d<T>& operator*=(T k);
  vector3d<T>& operator/=(T k);
//---------------------------------------------------------------------
  vector3d<T> operator-();
//---------------------------------------------------------------------
  friend vector3d<T> operator+(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return vector3d<T>(u.name_ + "+" + v.name_, u.dims_,
                    { u[0] + v[0], u[1] + v[1], u[2] + v[2], 0} );
  }
  friend vector3d<T> operator-(const vector3d<T>& u, const vector3d<T>& v) {
    // implement code here     NOTE: FRIEND TEMPLATES ** MUST BE INLINED or use <> NOTATION
		//                               as in operator<< <>() case below
		//                               TAKE MY ADVICE -- GO WITH INLINING with friend templates
    check_equal_dims(u,v);
    return vector3d<T>(u.name_ + "-" + v.name_, u.dims_,
                    { u[0] - v[0], u[1] - v[1], u[2] - v[2], 0} );
  }
//---------------------------------------------------------------------
  friend vector3d<T> operator+(T k, const vector3d<T>& v) {
    return vector3d<T>(std::to_string(k) + "+" + v.name_, v.dims_, { k + v[0], k + v[1], k + v[2], 0 });
  }
  friend vector3d<T> operator+(const vector3d<T>& v, T k) { return k + v; }
//---------------------------------------------------------------------
  friend vector3d<T> operator-(T k, const vector3d<T>& v) {
    // implement code here
    return vector3d(std::to_string(k) + "-" + v.name_, v.dims_, {k-v[0], k-v[1], k-v[2]});
  }
  friend vector3d<T> operator-(const vector3d<T>& v, T k) {
		// implement code here
    //return vector3d(v.name_ + "-" + std::to_string(k), v.dims_, {k-v[0], k-v[1], k-v[2]};
    return -k+v;
 	}
  //---------------------------------------------------------------------
  friend vector3d<T> operator*(T k, const vector3d<T>& v) {
    // implement code here
    return v*k;
    //return vector3d(std::to_string(k) + "*" + v.name_, v.dims_, {k*v[0], k*v[1], k*v[2]});
  }
  friend vector3d<T> operator*(const vector3d<T>& v, T k) {
    // implement code here
    //return v * k;
    return vector3d(std::to_string(k) + "*" + v.name_, v.dims_, {k*v[0], k*v[1], k*v[2]});
	}
  //---------------------------------------------------------------------
  friend vector3d<T> operator/(const vector3d<T>& v, T k) {
    // implement code here
    if(k==0)
    {
      throw new std::invalid_argument("dividing by zero not possible");

    }
    return vector3d(v.name_ + "/" + std::to_string(k), v.dims_, {v[0]/k, v[1]/k, v[2]/k});
  }
//---------------------------------------------------------------------
  friend bool operator<(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return u.magnitude() < v.magnitude();
  }
  friend bool operator==(const vector3d<T>& u, const vector3d<T>& v) {
    // implement code here
    check_equal_dims(u,v);
    return u[0] == v[0] && u[1]==v[1] && u[2]==v[2];

  }
  friend bool operator!=(const vector3d<T>& u, const vector3d<T>& v) {
    // implement code here
    return !(u==v);
	}
//---------------------------------------------------------------------
  double dot(const vector3d<T>& other) const;
  double magnitude() const;
  double angle(const vector3d<T>& other) const;
  vector3d<T> cross(const vector3d<T>& other) const;
//---------------------------------------------------------------------
  static vector3d<T> zero();
//---------------------------------------------------------------------
	// uses <> notation for outlining AND has to prototype the class and method above the class
	// in general, avoid this and use inline friend templates
  friend std::ostream& operator<< <>(std::ostream& os, const vector3d<T>& v);

private:
  friend void check_equal_dims(const vector3d& u, const vector3d& v) {
    if (u.dims_ != v.dims_) { throw new std::invalid_argument("vector3d dims mismatch"); }
  }
  void check_bounds(int i) const;

private:
  std::string name_;
  int dims_;
  double data_[4];
};


template <typename T> vector3d<T>::vector3d() : vector3d("no_name", 3) {}  // 3d default dims
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims)
  : name_(name), dims_(dims) {
    memset(data_, 0, dims_ * sizeof(double));
    data_[3] = 0;  // vectors have 0 at end, pts have 1
  }
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims,
                                            const std::initializer_list<T>& li)
  : vector3d(name, dims) {
    int i = 0;
    for (T value : li) {
      if (i > dims_) { break; }
      data_[i++] = value;
    }
    data_[3] = 0;
  }

//-----------------------------------------------------------------------
template <typename T>
T  vector3d<T>::operator[](int i) const {  check_bounds(i);  return data_[i]; }  // read-only index operator

template <typename T>
T& vector3d<T>::operator[](int i) { /* read-write index operator */
	// implement code here
  check_bounds(i);
  return data_[i];
}
//-----------------------------------------------------------------------
template <typename T>
std::string vector3d<T>::name() const { return name_; }

template <typename T>
void vector3d<T>::name(const std::string& name) { name_ = name; }
//-----------------------------------------------------------------------
template <typename T>  vector3d<T>& vector3d<T>::operator+=(const vector3d<T>& v)  {
  vector3d<T>& u = *this;
  for (int i = 0; i < 3; ++i) { u[i] += v[i]; }
  return *this;
}
template <typename T>  vector3d<T>& vector3d<T>::operator-=(const vector3d<T>& v)  {
   // implement code here
   return *this = -v;

}
//-----------------------------------------------------------------------
template <typename T>  vector3d<T>& vector3d<T>::operator+=(T k)  {
   // implement code here
   vector3d<T> & u = * this;
   //u[0] += k;
   //u[1] += k;
   //u[2] += k;
   for(int i = 0; i<3; i++)
   {
     u[i] += k;
   }
  return *this;
}
template <typename T>  vector3d<T>& vector3d<T>::operator-=(T k)  {
   // implement code here
   *this = -k;
   return *this;
}
template <typename T>  vector3d<T>& vector3d<T>::operator*=(T k)  {
   // implement code here
   vector3d<T> &u = *this;
   for(int i = 0; i<3; i++)
   {
     u[i] *= k;
   }
   return *this;
}
template <typename T>  vector3d<T>& vector3d<T>::operator/=(T k)  {
   // implement code here
   double div = 1.0/k;
   *this *= div;
   return *this;
};

//-----------------------
template <typename T>
vector3d<T> vector3d<T>::operator-() {
	return vector3d("-" + name_, dims_, { -data_[0], -data_[1], -data_[2], 0 });
}
//-----------------------
template <typename T>
double vector3d<T>::dot(const vector3d<T>& v) const {
   // implement code here
   const vector3d<T> &u = *this;
   check_equal_dims(*this, v);
   return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}
//-----------------------
template <typename T>
double vector3d<T>::magnitude() const {  return sqrt(dot(*this));  }

template <typename T>
double vector3d<T>::angle(const vector3d<T>& v) const {
   // implement code here
   double dot = this->dot(v);
   double mag = this->magnitude();
   double other = v.magnitude();
   return acos(dot/ (mag * other));
}
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::cross(const vector3d<T>& v) const {
  check_equal_dims(*this, v);
  if (v.dims_ != 3) {
    throw new std::invalid_argument("cross_product only implemented for vector3d's");
  }
  return vector3d<T>(name_ + " x " + v.name_, dims_,
                    { data_[1]*v[2] - data_[2]*v[1], -(data_[0]*v[2] - data_[2]*v[0]),
                      data_[0]*v[1] - data_[1]*v[0], 0 });
}
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::zero() { return vector3d<T>("zero", 3, {0, 0, 0, 0}); }
//-----------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const vector3d<T>& v) {
  os << "<" << v.name_ << ", ";
  if (v.dims_ == 0) { os << "empty>"; }
  else {
    for (int i = 0; i < v.dims_ + 1; ++i) {
      os << v[i];
      if (i < v.dims_) { os << " "; }
    }
    os << ">";
  }
  return os;
}

template <typename T>
void vector3d<T>::check_bounds(int i) const {  // 1 extra dimension for pts/vectors
   // implement code here
   if(i>dims_)
   {
     throw new std::invalid_argument("out of bounds");
   }
}


#endif
