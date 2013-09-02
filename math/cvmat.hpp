#ifndef math_cvmat_hpp_included_
#define math_cvmat_hpp_included_

#include <opencv2/core/core.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace math {

//! Converts an ublas matrix to an OpenCV matrix
template<typename T, typename L, typename C>
cv::Mat cvMat(const boost::numeric::ublas::matrix<T,L,C>& input)
{
    cv::Mat mat(input.size1(), input.size2(), CV_64F);
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++)
            mat.at<double>(i, j) = input(i, j);
    }
    return mat;
}

} // namespace math

#endif // math_cvmat_hpp_included_
