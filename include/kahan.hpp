#ifndef KAHAN_SUMMATION_LINEAR_ALGEBRA_HPP 
#define KAHAN_SUMMATION_LINEAR_ALGEBRA_HPP 

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/*
 * Implementations of linear algebraic operations with Eigen matrices/arrays 
 * that employ Kahan's compensated summation algorithm. 
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School 
 * Last updated:
 *     11/25/2021
 */
namespace Eigen {

namespace KahanSum {

template <typename T>
T vectorSum(const std::vector<T>& v)
{
    /*
     * Return the sum of the entries in the given vector, computed using 
     * Kahan's compensated summation algorithm. 
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    T sum, err, delta, newsum; 
    sum = 0; 
    err = 0; 
    for (const T x : v)
    {
        delta = x - err; 
        newsum = sum + delta; 
        err = (newsum - sum) - delta; 
        sum = newsum; 
    }

    return sum; 
}

template <typename Derived, typename T = typename Derived::Scalar>
T trace(const Eigen::MatrixBase<Derived>& A)
{
    /*
     * Return the trace of A, computed using Kahan's compensated summation
     * algorithm.
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization.  
     */
    T sum, err, delta, newsum; 
    sum = 0; 
    err = 0;
    int dim = std::min(A.rows(), A.cols()); 
    for (unsigned i = 0; i < dim; ++i)
    {
        delta = A(i, i) - err; 
        newsum = sum + delta;  
        err = (newsum - sum) - delta; 
        sum = newsum; 
    } 

    return sum; 
}

template <typename Derived, typename T = typename Derived::Scalar>
Eigen::Matrix<T, Derived::RowsAtCompileTime, 1> rowSum(const Eigen::MatrixBase<Derived>& A)
{
    /*
     * Return the column vector containing the row sums of A, computed using
     * Kahan's compensated summation algorithm.
     *
     * This function takes *dense* matrices as input and returns *dense* vectors.  
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    Eigen::Matrix<T, Derived::RowsAtCompileTime, 1> row_sums(A.rows());
    T sum, err, delta, newsum;  
    for (unsigned i = 0; i < A.rows(); ++i)
    {
        sum = 0;
        err = 0; 
        for (unsigned j = 0; j < A.cols(); ++j)
        {
            delta = A(i, j) - err; 
            newsum = sum + delta; 
            err = (newsum - sum) - delta; 
            sum = newsum; 
        }
        row_sums(i) = sum; 
    }

    return row_sums; 
}

template <typename Derived, typename T = typename Derived::Scalar>
T rowSum(const Eigen::MatrixBase<Derived>& A, const int i)
{
    /*
     * Return the i-th row sum of A, computed using Kahan's compensated
     * summation algorithm. 
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    T sum, err, delta, newsum; 
    sum = 0;
    err = 0; 
    for (unsigned j = 0; j < A.cols(); ++j)
    {
        delta = A(i, j) - err; 
        newsum = sum + delta; 
        err = (newsum - sum) - delta; 
        sum = newsum; 
    }

    return sum; 
}

template <typename Derived, typename T = typename Derived::Scalar>
Eigen::Matrix<T, 1, Derived::ColsAtCompileTime> colSum(const Eigen::MatrixBase<Derived>& A)
{
    /*
     * Return the row vector containing the column sums of A, computed using
     * Kahan's compensated summation algorithm.
     *
     * This function takes *dense* matrices as input and returns *dense* vectors.  
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    Eigen::Matrix<T, 1, Derived::ColsAtCompileTime> col_sums(A.cols());
    T sum, err, delta, newsum;  
    for (unsigned i = 0; i < A.cols(); ++i)
    {
        sum = 0;
        err = 0; 
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            delta = A(j, i) - err; 
            newsum = sum + delta; 
            err = (newsum - sum) - delta; 
            sum = newsum; 
        }
        col_sums(i) = sum; 
    }

    return col_sums; 
}

template <typename Derived, typename T = typename Derived::Scalar>
T colSum(const Eigen::MatrixBase<Derived>& A, const int i)
{
    /*
     * Return the i-th column sum of A, computed using Kahan's compensated
     * summation algorithm. 
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    T sum, err, delta, newsum; 
    sum = 0;
    err = 0; 
    for (unsigned j = 0; j < A.rows(); ++j)
    {
        delta = A(j, i) - err; 
        newsum = sum + delta; 
        err = (newsum - sum) - delta; 
        sum = newsum; 
    }

    return sum; 
}

template <typename Derived, typename T = typename Derived::Scalar>
T sum(const Eigen::MatrixBase<Derived>& A)
{
    /*
     * Return the sum of the entries in a matrix, computed using Kahan's
     * compensated summation algorithm. 
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization. 
     */
    T sum, err, delta, newsum; 
    sum = 0; 
    err = 0;
    for (unsigned i = 0; i < A.rows(); ++i)
    {
        for (unsigned j = 0; j < A.cols(); ++j)
        {
            delta = A(i, j) - err; 
            newsum = sum + delta;
            err = (newsum - sum) - delta; 
            sum = newsum;
        }
    }
    
    return sum; 
}

template <typename DerivedA, typename DerivedB, typename T = typename DerivedA::Scalar>
T dotProduct(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedB>& B,
             const int i, const int j)
{
    /*
     * Return the dot product of the i-th row of A and the j-th column of B,
     * computed using Kahan's compensated summation algorithm.
     *
     * The two matrices are assumed to have matching dimensions (i.e.,
     * A.cols() == B.rows()).
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization.  
     */
    T sum, err, prod, delta, newsum; 
    sum = 0;
    err = 0; 
    for (unsigned k = 0; k < A.cols(); ++k)
    {
        prod = A(i, k) * B(k, j);
        delta = prod - err; 
        newsum = sum + delta;  
        err = (newsum - sum) - delta; 
        sum = newsum; 
    } 

    return sum; 
}

template <typename Derived, typename T = typename Derived::Scalar>
T dotProduct(const Eigen::SparseMatrix<T, Eigen::RowMajor>& A, const Eigen::MatrixBase<Derived>& B,
             const int i, const int j)
{
    /*
     * Return the dot product of the i-th row of A and the j-th column of B,
     * computed using Kahan's compensated summation algorithm, where A is 
     * specified as a *compressed sparse row-major* matrix. 
     *
     * The two matrices are assumed to have matching dimensions (i.e.,
     * A.cols() == B.rows()).
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization.  
     */
    T sum, err, prod, delta, newsum; 
    sum = 0;
    err = 0;

    // Pick out only the nonzero entries in the i-th row of A
    for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it) 
    {
        prod = it.value() * B(it.col(), j);
        delta = prod - err; 
        newsum = sum + delta;  
        err = (newsum - sum) - delta; 
        sum = newsum; 
    } 

    return sum; 
}

template <typename DerivedA, typename DerivedB, typename T = typename DerivedA::Scalar>
Eigen::Matrix<T, DerivedA::RowsAtCompileTime, DerivedB::ColsAtCompileTime> multiply(const Eigen::MatrixBase<DerivedA>& A,
                                                                                    const Eigen::MatrixBase<DerivedB>& B) 
{
    /*
     * Return the matrix product A * B, computed using Kahan's compensated 
     * summation algorithm.
     *
     * The two matrices are assumed to have matching dimensions (i.e.,
     * A.cols() == B.rows()).
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization.  
     */ 
    Eigen::Matrix<T, DerivedA::RowsAtCompileTime, DerivedB::ColsAtCompileTime> prod(A.rows(), B.cols());
    for (unsigned i = 0; i < A.rows(); ++i)
    {
        for (unsigned j = 0; j < B.cols(); ++j)
        {
            prod(i, j) = dotProduct(A, B, i, j); 
        }
    }
    
    return prod; 
}

template <typename T, typename Derived>
Eigen::Matrix<T, Eigen::Dynamic, Derived::ColsAtCompileTime> multiply(const Eigen::SparseMatrix<T, Eigen::RowMajor>& A,
                                                                      const Eigen::MatrixBase<Derived>& B)
{
    /*
     * Return the matrix product A * B, computed using Kahan's compensated 
     * summation algorithm, where A is specified as a *compressed sparse 
     * row-major* matrix
     *
     * The two matrices are assumed to have matching dimensions (i.e.,
     * A.cols() == B.rows()).
     *
     * Note that this algorithm may not be effective at preserving floating-
     * point accuracy if the compiler is overly aggressive at optimization.  
     */
    Eigen::Matrix<T, Eigen::Dynamic, Derived::ColsAtCompileTime> prod(A.rows(), B.cols());
    for (unsigned i = 0; i < A.rows(); ++i)
    {
        for (unsigned j = 0; j < B.cols(); ++j)
        {
            prod(i, j) = dotProduct(A, B, i, j); 
        }
    }
    
    return prod; 
}

}   // namespace KahanSum

}   // namespace Eigen

#endif
