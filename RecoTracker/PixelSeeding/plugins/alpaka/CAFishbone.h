#ifndef RecoTracker_PixelSeeding_plugins_alpaka_CAFishbone_h
#define RecoTracker_PixelSeeding_plugins_alpaka_CAFishbone_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Math/interface/approx_atan2.h"
#include "HeterogeneousCore/AlpakaInterface/interface/VecArray.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "CACell.h"
#include "CAStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::caPixelDoublets {

  template <typename TrackerTraits>
  using CellNeighbors = caStructures::CellNeighborsT<TrackerTraits>;
  template <typename TrackerTraits>
  using CellTracks = caStructures::CellTracksT<TrackerTraits>;
  template <typename TrackerTraits>
  using CellNeighborsVector = caStructures::CellNeighborsVectorT<TrackerTraits>;
  template <typename TrackerTraits>
  using CellTracksVector = caStructures::CellTracksVectorT<TrackerTraits>;
  template <typename TrackerTraits>
  using OuterHitOfCell = caStructures::OuterHitOfCellT<TrackerTraits>;
  template <typename TrackerTraits>
  using HitsConstView = typename CACellT<TrackerTraits>::HitsConstView;

  template <typename TrackerTraits>
  class CAFishbone {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  HitsConstView<TrackerTraits> hh,
                                  CACellT<TrackerTraits>* cells,
                                  uint32_t const* __restrict__ nCells,
                                  OuterHitOfCell<TrackerTraits> const* isOuterHitOfCellWrap,
                                  int32_t nHits,
                                  bool checkTrack) const {
      constexpr auto maxCellsPerHit = CACellT<TrackerTraits>::maxCellsPerHit;

      int32_t layer2Offset = isOuterHitOfCellWrap->offset;
      // if there are no hits outside of the BPIX1, there is nothing to do
      if (nHits <= layer2Offset)
        return;

      auto const isOuterHitOfCell = isOuterHitOfCellWrap->container;

      float x[maxCellsPerHit], y[maxCellsPerHit], z[maxCellsPerHit], n[maxCellsPerHit];
      uint32_t cc[maxCellsPerHit];
      uint16_t d[maxCellsPerHit];
      uint8_t l[maxCellsPerHit];

      // index Dim<TAcc>::value - 1 runs faster...
      constexpr uint32_t dimIndexX = 1u;
      const uint32_t blockDimensionX(alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[dimIndexX]);
      const auto [firstElementIdxNoStrideX, endElementIdxNoStrideX] =
          cms::alpakatools::element_index_range_in_block(acc, 0u, dimIndexX);

      // Outermost loop on Y
      constexpr uint32_t dimIndexY = 0u;
      const uint32_t gridDimensionY(alpaka::getWorkDiv<alpaka::Grid, alpaka::Elems>(acc)[dimIndexY]);
      const auto [firstElementIdxNoStrideY, endElementIdxNoStrideY] =
          cms::alpakatools::element_index_range_in_grid(acc, 0u, dimIndexY);
      uint32_t firstElementIdxY = firstElementIdxNoStrideY;
      uint32_t endElementIdxY = endElementIdxNoStrideY;

      for (uint32_t idy = firstElementIdxY, nt = nHits - layer2Offset; idy < nt; ++idy) {
        if (not cms::alpakatools::next_valid_element_index_strided(
                idy, firstElementIdxY, endElementIdxY, gridDimensionY, nt))
          break;

        auto const& vc = isOuterHitOfCell[idy];
        auto size = vc.size();
        if (size < 2)
          continue;
        // if alligned kill one of the two.
        // in principle one could try to relax the cut (only in r-z?) for jumping-doublets
        auto const& c0 = cells[vc[0]];
        auto xo = c0.outer_x(hh);
        auto yo = c0.outer_y(hh);
        auto zo = c0.outer_z(hh);
        auto sg = 0;
        for (int32_t ic = 0; ic < size; ++ic) {
          auto& ci = cells[vc[ic]];
          if (ci.unused())
            continue;  // for triplets equivalent to next
          if (checkTrack && ci.tracks().empty())
            continue;
          cc[sg] = vc[ic];
          l[sg] = ci.layerPairId();
          d[sg] = ci.inner_detIndex(hh);
          x[sg] = ci.inner_x(hh) - xo;
          y[sg] = ci.inner_y(hh) - yo;
          z[sg] = ci.inner_z(hh) - zo;
          n[sg] = x[sg] * x[sg] + y[sg] * y[sg] + z[sg] * z[sg];
          ++sg;
        }
        if (sg < 2)
          continue;

        // here we parallelize in X
        uint32_t firstElementIdxX = firstElementIdxNoStrideX;
        uint32_t endElementIdxX = endElementIdxNoStrideX;
        for (uint32_t ic = firstElementIdxX; (int)ic < sg - 1; ++ic) {
          if (not cms::alpakatools::next_valid_element_index_strided(
                  ic, firstElementIdxX, endElementIdxX, blockDimensionX, sg - 1))
            break;

          auto& ci = cells[cc[ic]];
          for (auto jc = ic + 1; (int)jc < sg; ++jc) {
            auto& cj = cells[cc[jc]];
            // must be different detectors (in the same layer)
            // if (d[ic]==d[jc]) continue;
            auto cos12 = x[ic] * x[jc] + y[ic] * y[jc] + z[ic] * z[jc];

            if (d[ic] != d[jc] && cos12 * cos12 >= 0.99999f * (n[ic] * n[jc])) {
              // alligned:  kill farthest (prefer consecutive layers)
              // if same layer prefer farthest (longer level arm) and make space for intermediate hit
              bool sameLayer = l[ic] == l[jc];
              if (n[ic] > n[jc]) {
                if (sameLayer) {
                  cj.kill();  // closest
                  ci.setFishbone(acc, cj.inner_hit_id(), cj.inner_z(hh), hh);
                } else {
                  ci.kill();  // farthest
                  // break;  // removed to improve reproducibility, keep it for reference and tests
                }
              } else {
                if (!sameLayer) {
                  cj.kill();  // farthest
                } else {
                  ci.kill();  // closest
                  cj.setFishbone(acc, ci.inner_hit_id(), ci.inner_z(hh), hh);
                  // break;  // removed to improve reproducibility, keep it for reference and tests
                }
              }
            }
          }  // cj
        }    // ci
      }      // hits
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::caPixelDoublets

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_CAFishbone_h
