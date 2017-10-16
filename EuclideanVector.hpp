#ifndef EUCLIDEANVECTOR_HPP
#define EUCLIDEANVECTOR_HPP

#include <iostream>
#include <new>
#include <vector>
#include <list>
#include <initializer_list>
#include <cmath>

namespace evec {

class EuclideanVector {
    // Nonmember functions.
    friend bool operator==(const EuclideanVector& v, const EuclideanVector& w);
    friend bool operator!=(const EuclideanVector& v, const EuclideanVector& w);
    friend EuclideanVector operator+(const EuclideanVector& v, const EuclideanVector& w);
    friend EuclideanVector operator-(const EuclideanVector& v, const EuclideanVector& w);
    friend double operator*(const EuclideanVector& v, const EuclideanVector& w);
    friend EuclideanVector operator*(const EuclideanVector& v, const double& x);
    friend EuclideanVector operator*(double x, const EuclideanVector& v);
    friend EuclideanVector operator/(const EuclideanVector& v, const double& x);
    friend std::ostream& operator<<(std::ostream& os, const EuclideanVector& v);
private:
    // Data members.
    double* start;
    double* finish;
    
    // Allocator functions. Consider creating a base class to handle memory.
    double* allocate(unsigned int dim, double*) {
        double* buffer = (double*)(::operator new((size_t)(dim * sizeof (double))));
        if (buffer == 0) { std::cerr << "Stack overflow." << std::endl; exit(1); }
        return buffer;
    }
    void construct(double* p, const double& value) {
        ::new(static_cast<void*>(&*p)) double(value);
    }
    void deallocate(double* buffer) {
        ::operator delete (buffer);
    }
    double* copy_n(double* start, unsigned int dim, const double& mag) {
        while (dim--) construct(&*start++, mag); return start; // start is one past the last element.
    }
    double* copy_list(double* first, double* last, double* start) {
        while (first != last) construct(&*start++, *first++); return start; // start is one past the last element.
    }
public:
    // Constructors and destructor.
    EuclideanVector(unsigned int dim = 1) {
        start = allocate(dim, (double*)0);
        finish = copy_n(start, dim, 0);
    }
    EuclideanVector(unsigned int dim, double mag) {
        start = allocate(dim, (double*)0);
        finish = copy_n(start, dim, mag);
    }
    EuclideanVector(std::vector<double>::iterator first, std::vector<double>::iterator end) {
        unsigned int dim = (unsigned int)(end - first);
        start = allocate(dim, (double*)0);
        finish = start;
        for (auto it = first; it != end; ++it) construct(&*finish++, *it); // Some difficulties here due to const iterator.
    }
    EuclideanVector(std::list<double>::iterator first, std::list<double>::iterator end) {
        unsigned int dim = 0;
        for (auto it = first; it != end; ++it) dim++;
        start = allocate(dim, (double*)0);
        finish = start;
        for (auto it = first; it != end; ++it) construct(&*finish++, *it); // Some difficulties here due to const iterator.
    }
    EuclideanVector(std::initializer_list<double> init) {
        start = allocate((unsigned int)(init.size()), (double*)0);
        finish = start;
        for (auto it = init.begin(); it != init.end(); ++it) construct(&*finish++, *it); // Some difficulties here due to const iterator.
    }
    EuclideanVector(const EuclideanVector& v) {
        unsigned int dim = (unsigned int)(v.finish - v.start);
        start = allocate(dim, (double*)0);
        finish = copy_list(v.start, v.finish, start); }
    EuclideanVector(EuclideanVector&& v) : start(v.start), finish(v.finish) {
        v.start = v.finish = nullptr;
    }
    ~EuclideanVector() = default;
    
    // Copy assignment.
    EuclideanVector& operator=(const EuclideanVector& v) { // Object already exists so don't allocate memory, keep old object.
        if (this != &v) start = v.start; finish = v.finish;
        return *this;
    }
    // Move assignment.
    EuclideanVector& operator=(EuclideanVector&& v) { // Object already exists so don't allocate memory, discard old object.
        start = finish = nullptr;
        std::swap(*this, v);
        return *this;
    }
    
    // Member functions.
    unsigned int getNumDimensions() { return (unsigned int)(finish - start); }
    double get(unsigned int dim_i) { return *(start + dim_i); }
    double getEuclideanNorm() {
        double sum = 0;
        double* curr = start;
        while (curr != finish) { sum += pow(*curr++, 2); }
        return std::sqrt(sum);
    }
    EuclideanVector createUnitVector() {
        EuclideanVector v{*this}; // call the copy constructor
        double v_norm = v.getEuclideanNorm();
        double* curr = v.start;
        while (curr != v.finish) *curr++ /= v_norm;
        return v;
    }
    
    // Operator overloading.
    const double operator[](unsigned int dim_i) const { return *(start + dim_i); }
    double& operator[](unsigned int dim_i) { return *(start + dim_i); }
    void operator+=(const EuclideanVector& v) {
        double* curr = start;
        double* durr = v.start;
        while (curr != finish) { *curr++ += *durr++; }
    }
    void operator-=(const EuclideanVector& v) {
        double* curr = start;
        double* durr = v.start;
        while (curr != finish) { *curr++ -= *durr++; }
    }
    void operator*=(double x) {
        double* curr = start;
        while (curr != finish) { *curr++ *= x; }
    }
    void operator/=(double x) {
        double* curr = start;
        while (curr != finish) { *curr++ /= x; }
    }
    // Typecast overloading.
    operator std::vector<double>() {
        std::vector<double> res(start, finish);
        return res;
    }
    operator std::list<double>() {
        std::list<double> res(start, finish);
        return res;
    }
};

// Nonmember functions.
bool operator==(const EuclideanVector& v, const EuclideanVector& w) {
    for (int i = 0; i < v.finish - v.start; ++i) {
        if (v[i] != w[i]) return false;
    }
    return true;
}
bool operator!=(const EuclideanVector& v, const EuclideanVector& w) {
    return !(v == w);
}
EuclideanVector operator+(const EuclideanVector& v, const EuclideanVector& w) {
    if ( v.finish - v.start != w.finish - w.start ) { throw std::runtime_error("Unequal Euclidean vector sizes."); }
    EuclideanVector res{v};
    std::ptrdiff_t dim = res.finish - res.start;
    for (int i = 0; i < dim; ++i) res[i] += w[i];
    return res;
}
EuclideanVector operator-(const EuclideanVector& v, const EuclideanVector& w) {
    if ( v.finish - v.start != w.finish - w.start ) { throw std::runtime_error("Unequal Euclidean vector sizes."); }
    EuclideanVector res{v};
    std::ptrdiff_t dim = res.finish - res.start;
    for (int i = 0; i < dim; ++i) res[i] -= w[i];
    return res;
}
double operator*(const EuclideanVector& v, const EuclideanVector& w) {
    double sum = 0;
    unsigned int dim = (unsigned int)(v.finish - v.start);
    for (int i = 0; i < dim; ++i) sum += v[i] * w[i];
    return sum;
}
EuclideanVector operator*(const EuclideanVector& v, const double& x) {
    EuclideanVector res{v};
    std::ptrdiff_t dim = res.finish - res.start;
    for (int i = 0; i < dim; ++i) res[i] *=  x;
    return res;
}
EuclideanVector operator*(double x, const EuclideanVector& v) {
    EuclideanVector res{v};
    std::ptrdiff_t dim = res.finish - res.start;
    for (int i = 0; i < dim; ++i) res[i] *=  x;
    return res;
}
EuclideanVector operator/(const EuclideanVector& v, const double& x) {
    EuclideanVector res{v};
    std::ptrdiff_t dim = res.finish - res.start;
    for (int i = 0; i < dim; ++i) res[i] /=  x;
    return res;
}
std::ostream& operator<<(std::ostream& os, const EuclideanVector& v) {
    os << "[";
    double* curr = v.start;
    double* last = v.finish;
    last -= 1;
    while (curr != v.finish)
        (curr == last) ? os << *curr++ : os << *curr++ << " ";
    os << "]";
    return os;
}

} /* namespace evec */

#endif /* EUCLIDEANVECTOR_HPP */
