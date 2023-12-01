#ifndef L1Trigger_L1TMuonEndCapPhase2_NdArrayDesc_h
#define L1Trigger_L1TMuonEndCapPhase2_NdArrayDesc_h

#include <array>
#include <initializer_list>
#include <memory>

#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"

namespace emtf::phase2 {

    // Adapted from:
    // https://github.com/tensorflow/tensorflow/blob/master/tensorflow/lite/kernels/internal/common.h
    class NdArrayDesc {
        public:
            NdArrayDesc(std::initializer_list<int>&& extents);

            ~NdArrayDesc();

            const unsigned int& num_dimensions() const;

            const unsigned int& num_elements() const;

            // Compute 1-D index from n-D index: This will begin the recursion
            template<typename... ArrIndices>
                typename std::enable_if<
                (std::is_integral<ArrIndices>::value && ...), unsigned int>::
                type get_index(const ArrIndices&... indices) const {
                    return get_index_next_dim_contrib(0, std::forward<const ArrIndices>(indices)...);
                }

        private:
            unsigned int num_dimensions_;
            unsigned int num_elements_;

            std::unique_ptr<unsigned int[]> extents_;  // the extent of each dimension
            std::unique_ptr<unsigned int[]> strides_;  // the number of elements between consecutive dimensions

            // This will recursively call itself eachtime with one argument less
            template<typename Index, typename ...ArrIndices>
                typename std::enable_if<
                std::is_integral<Index>::value
                && (std::is_integral<ArrIndices>::value && ...), unsigned int>::
                type get_index_next_dim_contrib(
                        const unsigned int& i_dim, const Index& index, const ArrIndices&... indices) const {
                    emtf_assert(index < extents_[i_dim]); // Requested index beyond the dimension's extension.
                    emtf_assert(i_dim < num_dimensions_); // Requested more dimensions than those available

                    unsigned int contrib = index * strides_[i_dim];
                    unsigned int next_contrib = get_index_next_dim_contrib(
                            i_dim + 1,
                            std::forward<const ArrIndices>(indices)...);

                    return contrib + next_contrib;
                }

            // When all arguments run out, this function will be called
            unsigned int get_index_next_dim_contrib(const unsigned int&) const;

    };


}  // namespace emtf::phase2

#endif  // L1Trigger_L1TMuonEndCapPhase2_NdArrayDesc_h not defined
