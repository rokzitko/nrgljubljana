// Rok Zitko, Nov 2020

#ifndef H5EASY_BITS_UBLAS_MATRIX_HPP
#define H5EASY_BITS_UBLAS_MATRIX_HPP

#include <boost/numeric/ublas/matrix.hpp>

#include "../H5Easy.hpp"
#include "H5Easy_misc.hpp"
#include "H5Easy_scalar.hpp"

namespace H5Easy {

namespace detail {

using namespace boost::numeric;

template <class T>
struct is_ublas_matrix : std::false_type {};
template <class T>
struct is_ublas_matrix<ublas::matrix<T>> : std::true_type {};

template <typename T>
struct io_impl<T, typename std::enable_if<is_ublas_matrix<T>::value>::type> {

    inline static DataSet dump(File& file,
                               const std::string& path,
                               const T& data,
                               const DumpOptions& options) {
        using value_type = typename T::value_type;
        std::vector<size_t> dims(2);
        dims[0] = data.size1();
        dims[1] = data.size2();
        DataSet dataset = initDataset<value_type>(file, path, dims, options);
        dataset.write_raw(&data(0,0));
        if (options.flush()) {
            file.flush();
        }
        return dataset;
    }

};

}  // namespace detail
}  // namespace H5Easy

#endif  // H5EASY_BITS_UBLAS_MATRIX_HPP
