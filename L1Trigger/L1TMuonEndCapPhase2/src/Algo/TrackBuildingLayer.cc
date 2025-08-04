#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFContext.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DataUtils.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/Algo/TrackBuildingLayer.h"

using namespace emtf::phase2;
using namespace emtf::phase2::algo;

// Static
seg_theta_t TrackBuildingLayer::calc_theta_median(std::vector<seg_theta_t> thetas) {
  auto i_last = thetas.size() - 1;

  // Sort Thetas
  // This will sort ascending order (lower-value means lower-index)
  data::mergesort(&thetas[0], thetas.size(), [](seg_theta_t lower_index_value, seg_theta_t larger_index_value) -> int {
    return lower_index_value > larger_index_value ? 1 : 0;
  });

  // Check if any model_thm_site is null
  // Since the theta array has been sorted, it's enough
  // to check the last index, because the invalid value will be the max
  seg_theta_t invalid_theta = -1;  // This maps to 255 since it underflows

  bool any_invalid = thetas[i_last] == invalid_theta;

  // Calculate median
  if (any_invalid) {
    // Use the min value as the median if there are any invalid thetas
    return thetas[0];
  } else {
    // Calculate the median if all thetas are valid
    return data::getMedianOfSorted(&thetas[0], thetas.size());
  }
}

// Members
TrackBuildingLayer::TrackBuildingLayer(const EMTFContext& context) : context_(context) {}

void TrackBuildingLayer::apply(const segment_collection_t& segments,
                               const std::vector<road_t>& roads,
                               const algo_id_t& algo,
                               std::vector<track_t>& tracks) const {
  // Apply
  for (unsigned int i_road = 0; i_road < roads.size(); ++i_road) {
    // Get road and track
    const auto& road = roads[i_road];
    auto& track = tracks.emplace_back();

    // Initialize track
    track.phi = 0;
    track.theta = 0;
    track.valid = 0;

    for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
      track.site_segs[site_id] = 0;
      track.site_mask[site_id] = 0;
      track.site_rm_mask[site_id] = 0;
    }

    for (unsigned int i_feature = 0; i_feature < v3::kNumTrackFeatures; ++i_feature) {
      track.features[i_feature] = 0;
    }

    // Short-Circuit: If the road has quality-0 skip it
    if (road.quality == 0) {
      continue;
    }

    // Debug Info
    if (this->context_.config_.verbosity_ > 1) {
      if (i_road == 0) {
        edm::LogInfo("L1TEMTFpp") << std::endl;
        edm::LogInfo("L1TEMTFpp") << "==========================================================================="
                                  << std::endl;
        edm::LogInfo("L1TEMTFpp") << "BEGIN TRACK BUILDING" << std::endl;
        edm::LogInfo("L1TEMTFpp") << "---------------------------------------------------------------------------"
                                  << std::endl;
      }

      edm::LogInfo("L1TEMTFpp") << "***************************************************************************"
                                << std::endl;
      edm::LogInfo("L1TEMTFpp") << "Begin building track " << i_road << std::endl;
    }

    // Build Track
    attachSegmentsByPhi(segments, road, algo, track);

    if (algo == algo_id_t::kBeamHalo) {
      attachSegmentsByRadius(segments, algo, track);
    } else {
      attachSegmentsByTheta(segments, algo, track);
    }

    // Debug Info
    if (this->context_.config_.verbosity_ > 1) {
      edm::LogInfo("L1TEMTFpp") << "End building track " << i_road << std::endl;

      if (i_road == (roads.size() - 1)) {
        edm::LogInfo("L1TEMTFpp") << "---------------------------------------------------------------------------"
                                  << std::endl;
        edm::LogInfo("L1TEMTFpp") << "END TRACK BUILDING" << std::endl;
        edm::LogInfo("L1TEMTFpp") << "==========================================================================="
                                  << std::endl;
      }
    }
  }
}

void TrackBuildingLayer::attachSegmentsByPhi(const segment_collection_t& segments,
                                             const road_t& road,
                                             const algo_id_t& algo,
                                             track_t& track) const {
  // ===========================================================================
  // Convert column center to emtf_phi units
  // ---------------------------------------------------------------------------
  // trk_col: Recall that the hitmap is 288 cols wide, and the full chamber hitmap is 315 cols;
  // the chamber hitmap doesn't fit in the hitmap, so we skipped the first 27 cols.
  // In order to calculate the full hitmap col, we need to add back the 27 cols that we skipped.
  // sector_col: The sector's column is the center col of the phi map adding back the 27 skipped cols.
  const trk_col_t trk_col = road.col + v3::kHitmapCropColStart;
  const trk_col_t sector_col = static_cast<trk_col_t>(v3::kHitmapNCols / 2) + v3::kHitmapCropColStart;

  // Each column is emtf_phi=1<<n wide, therefore half of this would be 1<<(n-1)
  // since shifting to right m, is the same as dividing by 2^m.
  seg_phi_t trk_abs_phi =
      (static_cast<seg_phi_t>(trk_col) << v3::kHitmapColFactorLog2) + (1 << (v3::kHitmapColFactorLog2 - 1));
  seg_phi_t sector_abs_phi =
      (static_cast<seg_phi_t>(sector_col) << v3::kHitmapColFactorLog2) + (1 << (v3::kHitmapColFactorLog2 - 1));

  // ===========================================================================
  // Store Data
  // ---------------------------------------------------------------------------
  track.zone = road.zone;
  track.col = trk_col;
  track.pattern = road.pattern;
  track.quality = road.quality;
  track.phi = trk_abs_phi;

  // Calculate track phi with respect to the sector center
  trk_feature_t trk_rel_phi = static_cast<trk_feature_t>(trk_abs_phi) - static_cast<trk_feature_t>(sector_abs_phi);

  // Collect track features
  track.features[36] = track.quality > 0 ? trk_rel_phi : decltype(trk_rel_phi)(0);
  track.features[38] = track.quality;

  // ===========================================================================
  // Unpack model
  // ---------------------------------------------------------------------------
  const auto& model = context_.model_;
  const model::zones::hitmap_t* model_hm;
  const model::zones::pattern_t* model_pat;

  if (algo == algo_id_t::kPrompt) {
    model_hm = &(model.zones_[track.zone].hitmap);
    model_pat = &(model.zones_[track.zone].prompt_patterns[track.pattern]);
  } else if (algo == algo_id_t::kDisplaced) {
    model_hm = &(model.zones_[track.zone].hitmap);
    model_pat = &(model.zones_[track.zone].disp_patterns[track.pattern]);
  } else if (algo == algo_id_t::kBeamHalo) {
    model_hm = &(model.bh_zones_[track.zone].hitmap);
    model_pat = &(model.bh_zones_[track.zone].patterns[track.pattern]);
  } else {
    model_hm = nullptr;
    model_pat = nullptr;
  }

  // ===========================================================================
  // Get pattern info for each row
  // ---------------------------------------------------------------------------
  std::array<trk_col_t, v3::kHitmapNRows> trk_pat_begin;
  std::array<trk_col_t, v3::kHitmapNRows> trk_pat_center;
  std::array<trk_col_t, v3::kHitmapNRows> trk_pat_end;
  std::array<seg_phi_t, v3::kHitmapNRows> trk_pat_phi;

  for (unsigned int i_row = 0; i_row < v3::kHitmapNRows; ++i_row) {
    // Get the model pattern
    const auto& model_pat_row = (*model_pat)[i_row];

    // Offset the pattern's begin, center, and end by the track column
    trk_pat_begin[i_row] = trk_col + model_pat_row.begin;
    trk_pat_center[i_row] = trk_col + model_pat_row.center;
    trk_pat_end[i_row] = trk_col + model_pat_row.end;
    trk_pat_phi[i_row] = 0;

    // Short-Circuit: If the pattern's center is less than the padding used
    // when matching the pattern to the hitmap then the pattern center is 0.
    // This is because at that point, the center is out-of-bounds.
    if (trk_pat_center[i_row] <= v3::kPatternMatchingPadding)
      continue;

    // When the center is beyond the padding, then the pattern
    // is in-bound, therefore we subtract the padding offset.
    // To get the center in terms of the non-padded row BW we need to remove padding
    // since col-padding + 1 should map to 1 in the non-padded hitmap
    const auto& temp_trk_pat_center = trk_pat_center[i_row] - v3::kPatternMatchingPadding;

    // Convert the pattern center to emtf_phi units
    trk_pat_phi[i_row] = (static_cast<seg_phi_t>(temp_trk_pat_center) << v3::kHitmapColFactorLog2) +
                         (1 << (v3::kHitmapColFactorLog2 - 1));
  }

  // ===========================================================================
  // Select segments using phi only
  // ---------------------------------------------------------------------------

  // Init segment phi differences
  std::array<seg_phi_t, v3::kNumTrackSites> trk_seg_phi_diff;

  for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
    trk_seg_phi_diff[site_id] = 0;
  }

  // clang-format off
    std::vector<std::vector<unsigned int>> site_chambers = {
        {  0,   1,   2,   9,  10,  11,  45}, // ME1/1
        {  3,   4,   5,  12,  13,  14,  46}, // ME1/2
        { 18,  19,  20,  48,  21,  22,  23,  24,  25,  26,  49}, // ME2/1 + ME2/2
        { 27,  28,  29,  50,  30,  31,  32,  33,  34,  35,  51}, // ME3/1 + ME3/2
        { 36,  37,  38,  52,  39,  40,  41,  42,  43,  44,  53}, // ME4/1 + ME4/2
        { 57,  58,  59,  66,  67,  68, 100}, // RE1/2
        { 75,  76,  77,  78,  79,  80, 103}, // RE2/2
        { 81,  82,  83, 104,  84,  85,  86,  87,  88,  89, 105}, // RE3/1 + RE3/2
        { 90,  91,  92, 106,  93,  94,  95,  96,  97,  98, 107}, // RE4/1 + RE4/2
        { 54,  55,  56,  63,  64,  65,  99}, // GE1/1 
        { 72,  73,  74, 102}, // GE2/1
        {108, 109, 110, 111, 112, 113, 114} // ME0
    };

    std::vector<unsigned int> site_chamber_orders = {
        0, 0, 2, 2, 2, 0, 0, 2, 2, 0, 1, 0
    };

    std::vector<std::vector<int>> chamber_orders = {
        {-1, -1,  6, -1,  0,  1, -1,  2,  3, -1,  4,  5},
        { 3, -1, -1,  0, -1, -1,  1, -1, -1,  2, -1, -1},
        { 3, -1, 10,  0,  4,  5,  1,  6,  7,  2,  8,  9}
    };
  // clang-format on

  // Select segments
  int zone_mask = (1u << track.zone);

  for (unsigned int i_row = 0; i_row < model_hm->size(); ++i_row) {  // Begin loop rows

    const auto& model_hm_row = (*model_hm)[i_row];

    const auto& trk_pat_row_begin = trk_pat_begin[i_row];
    const auto& trk_pat_row_end = trk_pat_end[i_row];
    const auto& trk_pat_row_phi = trk_pat_phi[i_row];

    if (this->context_.config_.verbosity_ > 2) {
      edm::LogInfo("L1TEMTFpp") << "Pattern Row:"
                                << " row " << i_row << " begin " << trk_pat_row_begin << " end " << trk_pat_row_end
                                << " phi " << trk_pat_row_phi << std::endl;
    }

    for (const auto& model_hm_site : model_hm_row) {  // Begin loop sites in row

      const int site_id = static_cast<int>(model_hm_site.id);

      auto& site_seg_id = track.site_segs[site_id];
      auto& site_bit = track.site_mask[site_id];
      auto& site_min_phi_diff = trk_seg_phi_diff[site_id];

      const auto& s_chambers = site_chambers[site_id];
      const auto& s_chamber_order_id = site_chamber_orders[site_id];
      const auto& s_chamber_order = chamber_orders[s_chamber_order_id];

      for (const auto& chamber_idx : s_chamber_order) {  // Begin loop chambers in site

        if (chamber_idx == -1)
          continue;

        int chamber_id = s_chambers[chamber_idx];

        for (unsigned int i_ch_seg = 0; i_ch_seg < v3::kChamberSegments; ++i_ch_seg) {  // Begin loop segments

          const int seg_id = chamber_id * v3::kChamberSegments + i_ch_seg;
          const auto& seg = segments[seg_id];

          // Short-Circuit: If the segment is invalid move on
          if (!seg.valid) {
            continue;
          }

          // Short-Circuit: If the segment is not in the zone move on
          if (algo != algo_id_t::kBeamHalo && (seg.zones & zone_mask) != zone_mask) {
            continue;
          }

          // Short-Circuit: If the segment is outside of the pattern move on
          const trk_col_t seg_col = (seg.phi >> 4) + v3::kPatternMatchingPadding;

          if (!(trk_pat_row_begin <= seg_col && seg_col <= trk_pat_row_end)) {
            continue;
          }

          // Calculate abs diff between the pattern's row phi and the segment's phi
          seg_phi_t diff;

          if (trk_pat_row_phi > seg.phi) {
            diff = trk_pat_row_phi - seg.phi;
          } else {
            diff = seg.phi - trk_pat_row_phi;
          }

          if (this->context_.config_.verbosity_ > 2) {
            edm::LogInfo("L1TEMTFpp") << "Site candidate:"
                                      << " site_id " << site_id << " seg_id " << seg_id << " seg_phi " << seg.phi
                                      << " seg_theta1 " << seg.theta1 << " seg_theta2 " << seg.theta2 << " seg_bend "
                                      << seg.bend << std::endl;
          }

          // Short-Circuit: If the difference is larger than the min diff move on
          if (site_bit == 1 && site_min_phi_diff <= diff)
            continue;

          // Select better segment
          site_seg_id = seg_id;
          site_bit = 1;
          site_min_phi_diff = diff;
        }  // End loop segments

      }  // End loop chambers in site

      // Debug Info
      if (this->context_.config_.verbosity_ > 2 && site_bit == 1) {
        edm::LogInfo("L1TEMTFpp") << "Segment attached:"
                                  << " site_id " << site_id << " seg_id " << site_seg_id << " seg_phi "
                                  << segments[site_seg_id].phi << " seg_theta1 " << segments[site_seg_id].theta1
                                  << " seg_theta2 " << segments[site_seg_id].theta2 << " seg_bend "
                                  << segments[site_seg_id].bend << std::endl;
      }
    }  // End loop sites in row
  }  // End loop rows
}

void TrackBuildingLayer::attachSegmentsByTheta(const segment_collection_t& segments,
                                               const algo_id_t& algo,
                                               track_t& track) const {
  // ===========================================================================
  // Constants
  // ---------------------------------------------------------------------------
  seg_theta_t invalid_theta = -1;  // This will map to 255 since it underflows

  // ===========================================================================
  // Calculate theta medians
  // ---------------------------------------------------------------------------
  const auto& model = context_.model_;
  const auto& model_thmc = model.theta_medians_;

  std::vector<seg_theta_t> theta_medians;

  for (const auto& method : model_thmc) {  // Begin loop methods

    std::vector<seg_theta_t> medians;

    for (const auto& group : method) {  // Begin loop groups

      std::vector<seg_theta_t> thetas;

      for (const auto& model_thm_site : group) {  // Begin loop sites
        int site_id = static_cast<int>(model_thm_site.id);

        const auto& site_bit = track.site_mask[site_id];

        // Initialize as invalid theta
        auto& theta = thetas.emplace_back(invalid_theta);

        // Short-Circuit: If no segment was selected, move on.
        if (site_bit == 0)
          continue;

        // Get selected segment's theta value
        const auto& site_seg_id = track.site_segs[site_id];
        const auto& site_seg = segments[site_seg_id];

        if (model_thm_site.theta_id == theta_id_t::kTheta1) {
          theta = site_seg.theta1;
        } else if (model_thm_site.theta_id == theta_id_t::kTheta2) {
          theta = site_seg.theta2;
        }

        // If the segment theta is 0 this is invalid theta value
        if (theta == 0) {
          theta = invalid_theta;
        }
      }  // End loop sites

      // Calculate theta median
      if (this->context_.config_.verbosity_ > 2) {
        for (const auto& theta : thetas) {
          edm::LogInfo("L1TEMTFpp") << "theta " << theta << std::endl;
        }
      }

      auto median = calc_theta_median(thetas);
      medians.push_back(median);

      if (this->context_.config_.verbosity_ > 2) {
        edm::LogInfo("L1TEMTFpp") << "group_median " << median << std::endl;
      }
    }  // End loop groups

    // Calculate method median
    auto median = calc_theta_median(medians);
    theta_medians.push_back(median);

    if (this->context_.config_.verbosity_ > 2) {
      edm::LogInfo("L1TEMTFpp") << "theta_median " << median << std::endl;
    }
  }  // End loop methods

  // ===========================================================================
  // Select track theta
  // ---------------------------------------------------------------------------
  seg_theta_t trk_abs_theta;

  if (track.zone != 2) {
    trk_abs_theta = theta_medians[0];
  } else {
    trk_abs_theta = theta_medians[1];
  }

  // If median is invalid, try station 1 median
  if (trk_abs_theta == invalid_theta) {
    trk_abs_theta = theta_medians[2];
  }

  // If all medians are invalid use 0 (0 is an invalid theta)
  if (trk_abs_theta == invalid_theta) {
    trk_abs_theta = 0;
  }

  // ===========================================================================
  // Store Data
  // ---------------------------------------------------------------------------
  track.theta = trk_abs_theta;
  track.valid = 1;

  // Collect track features
  track.features[37] = track.quality > 0 ? trk_abs_theta : decltype(trk_abs_theta)(0);

  // ===========================================================================
  // Compare segment theta to track theta
  // ---------------------------------------------------------------------------

  // Init segment site theta
  std::array<seg_theta_t, v3::kNumTrackSites> trk_seg_theta;

  for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
    trk_seg_theta[site_id] = 0;
  }

  // Get site theta windows
  std::vector<std::vector<seg_theta_t>> site_theta_window;

  if (algo == algo_id_t::kPrompt) {
    // clang-format off
        site_theta_window = {
            {5, 0, 2, 2, 2, 34, 0, 3, 3, 5, 6, 5},
            {5, 9, 5, 4, 5, 14, 7, 7, 7, 7, 7, 4},
            {11, 6, 5, 6, 6, 10, 8, 8, 9, 8, 0, 0}
        };
    // clang-format on
  } else {
    // clang-format off
        site_theta_window = {
            {14, 40, 4, 3, 3, 45, 0, 4, 4, 15, 8, 13},
            {16, 18, 7, 5, 5, 22, 7, 7, 8, 17, 9, 14},
            {26, 15, 8, 9, 9, 17, 11, 9, 10, 26, 21, 0}
        };
    // clang-format on
  }

  // Detach segments if theta_window < diff, it is invalid
  for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
    auto& site_bit = track.site_mask[site_id];
    auto& site_rm_bit = track.site_rm_mask[site_id];

    // Get Theta Window
    const auto& theta_window = site_theta_window[track.zone][site_id];

    // Short-Circuit: If no segment was selected, move on.
    if (site_bit == 0)
      continue;

    const auto& site_seg_id = track.site_segs[site_id];
    const auto& site_seg = segments[site_seg_id];

    // Init differences with out-of-bounds values
    seg_theta_t diff_1 = theta_window + 1;
    seg_theta_t diff_2 = theta_window + 1;

    // Calculate abs theta 1 diff
    if (site_seg.theta1 != 0) {
      if (site_seg.theta1 < trk_abs_theta) {
        diff_1 = trk_abs_theta - site_seg.theta1;
      } else {
        diff_1 = site_seg.theta1 - trk_abs_theta;
      }
    }

    // Calculate abs theta 2 diff
    if (site_seg.theta2 != 0) {
      if (site_seg.theta2 < trk_abs_theta) {
        diff_2 = trk_abs_theta - site_seg.theta2;
      } else {
        diff_2 = site_seg.theta2 - trk_abs_theta;
      }
    }

    // Select the theta with the smallest difference
    if (diff_1 <= diff_2 && diff_1 < theta_window) {
      // Select theta 1 as the correct theta value
      trk_seg_theta[site_id] = site_seg.theta1;
    } else if (diff_2 < theta_window) {
      // Select theta 2 as the correct theta value
      trk_seg_theta[site_id] = site_seg.theta2;
    } else {
      // Invalidate site if both differences are outside of the theta window
      site_bit = 0;
      site_rm_bit = 1;

      // Debug Info
      if (this->context_.config_.verbosity_ > 4) {
        edm::LogInfo("L1TEMTFpp") << "Segment outside of theta window; detatched:"
                                  << " site_id " << site_id << " seg_id " << site_seg_id << " seg_phi " << site_seg.phi
                                  << " seg_theta1 " << site_seg.theta1 << " seg_theta2 " << site_seg.theta2
                                  << std::endl;
      }
    }
  }

  // ===========================================================================
  // Collect Segment Features
  // ---------------------------------------------------------------------------
  const auto& model_ftc = model.features_;

  int i_feature = 0;

  for (auto& model_ft : model_ftc) {
    for (auto& model_ft_site : model_ft.sites) {
      int site_id = static_cast<int>(model_ft_site);

      const auto& site_seg_id = track.site_segs[site_id];
      const auto& site_bit = track.site_mask[site_id];
      const auto& site_seg = segments[site_seg_id];

      auto& trk_feature = track.features[i_feature++];

      // Short-Circuit: No segment attached
      if (site_bit == 0) {
        continue;
      }

      // Fill features
      if (model_ft.id == feature_id_t::kPhi) {
        // Note: This is the segment's phi with respect to the track's abs phi
        trk_feature = static_cast<trk_feature_t>(site_seg.phi) - static_cast<trk_feature_t>(track.phi);
      } else if (model_ft.id == feature_id_t::kTheta) {
        // Note: This is the segment's theta with respect to the track's abs theta
        trk_feature = static_cast<trk_feature_t>(trk_seg_theta[site_id]) - static_cast<trk_feature_t>(track.theta);
      } else if (model_ft.id == feature_id_t::kBend) {
        trk_feature = site_seg.bend;
      } else if (model_ft.id == feature_id_t::kQuality) {
        trk_feature = site_seg.qual1;
      }
    }
  }

  // Debug Info
  if (this->context_.config_.verbosity_ > 1) {
    edm::LogInfo("L1TEMTFpp") << "Track"
                              << " zone " << track.zone << " col " << track.col << " pat " << track.pattern << " qual "
                              << track.quality << " abs_phi " << track.phi << " abs_theta " << track.theta
                              << " features " << std::endl;

    for (unsigned int i = 0; i < v3::kNumTrackFeatures; ++i) {
      if (i > 0) {
        edm::LogInfo("L1TEMTFpp") << " ";
      }

      edm::LogInfo("L1TEMTFpp") << track.features[i];
    }

    edm::LogInfo("L1TEMTFpp") << std::endl;
  }
}

void TrackBuildingLayer::attachSegmentsByRadius(const segment_collection_t& segments,
                                                const algo_id_t& algo,
                                                track_t& track) const {
  // ===========================================================================
  // Constants
  // ---------------------------------------------------------------------------
  seg_rho_t invalid_rho = -1;  // This will map to 255 since it underflows

  std::vector<std::vector<seg_rho_t>> site_rho_lut = {
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       1,   2,   4,   5,   7,   8,   10,  11,  12,  14,  15,  17,  18,  20,  21,  23,  24,  26,  27,  29,  30,  32,
       33,  35,  37,  38,  40,  41,  43,  44,  46,  48,  49,  51,  52,  54,  56,  57,  59,  60,  62,  64,  65,  67,
       69,  71,  72,  74,  76,  77,  79,  81,  83,  85,  86,  88,  90,  92,  94,  95,  97,  99,  101, 103, 105, 107,
       109, 111, 112, 114, 116, 118, 120, 122, 124, 126, 129, 131, 133, 135, 137, 139, 141, 143, 146, 148, 150, 152,
       155, 157, 159, 161, 164, 166, 169, 171, 173, 176, 178, 181, 183, 186, 188, 191, 194, 196},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   3,   5,   6,   8,   10,  11,
       13,  15,  16,  18,  20,  22,  23,  25,  27,  28,  30,  32,  34,  35,  37,  39,  41,  42,  44,  46,  48,  49,
       51,  53,  55,  57,  59,  60,  62,  64,  66,  68,  70,  72,  73,  75,  77,  79,  81,  83,  85,  87,  89,  91,
       93,  95,  97,  99,  101, 103, 105, 107, 109, 111, 113, 115, 117, 120, 122, 124, 126, 128, 130, 133, 135, 137,
       139, 142, 144, 146, 149, 151, 153, 156, 158, 160, 163, 165, 168, 170, 173, 175, 178, 180, 183, 185, 188, 191,
       193, 196, 199, 201, 204, 207, 210, 212, 215, 218, 221, 224, 227, 230, 233, 236, 239, 242},
      {0,   0,   0,   0,   0,   0,   0,   0,   1,   3,   5,   7,   9,   11,  13,  15,  17,  19,  21,  23,  25,  27,
       29,  31,  33,  35,  37,  39,  41,  43,  45,  47,  49,  51,  53,  55,  58,  60,  62,  64,  66,  68,  70,  72,
       74,  77,  79,  81,  83,  85,  87,  90,  92,  94,  96,  99,  101, 103, 105, 108, 110, 112, 115, 117, 119, 122,
       124, 126, 129, 131, 134, 136, 138, 141, 143, 146, 148, 151, 153, 156, 158, 161, 164, 166, 169, 171, 174, 177,
       179, 182, 185, 188, 190, 193, 196, 199, 202, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 237, 240,
       244, 247, 250, 253, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   0,   0,   0,   2,   4,   6,   8,   10,  12,  15,  17,  19,  21,  23,  26,  28,  30,  32,  35,  37,  39,
       41,  43,  46,  48,  50,  53,  55,  57,  59,  62,  64,  66,  69,  71,  73,  76,  78,  80,  83,  85,  88,  90,
       92,  95,  97,  100, 102, 104, 107, 109, 112, 114, 117, 119, 122, 124, 127, 130, 132, 135, 137, 140, 143, 145,
       148, 151, 153, 156, 159, 161, 164, 167, 170, 172, 175, 178, 181, 184, 187, 189, 192, 195, 198, 201, 204, 207,
       210, 213, 216, 219, 222, 226, 229, 232, 235, 238, 242, 245, 248, 251, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   2,   4,   6,   9,   11,  13,  16,  18,  21,  23,  25,  28,  30,  33,  35,  38,  40,  42,  45,  47,  50,
       52,  55,  57,  60,  62,  65,  67,  70,  72,  75,  77,  80,  82,  85,  87,  90,  93,  95,  98,  100, 103, 106,
       108, 111, 114, 116, 119, 122, 124, 127, 130, 133, 135, 138, 141, 144, 147, 149, 152, 155, 158, 161, 164, 167,
       169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 203, 206, 209, 212, 215, 218, 222, 225, 228, 231, 235,
       238, 241, 245, 248, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   3,   4,   6,   8,   9,   11,  13,
       14,  16,  18,  19,  21,  23,  25,  26,  28,  30,  32,  33,  35,  37,  39,  40,  42,  44,  46,  48,  49,  51,
       53,  55,  57,  59,  60,  62,  64,  66,  68,  70,  72,  74,  76,  77,  79,  81,  83,  85,  87,  89,  91,  93,
       95,  97,  99,  101, 103, 105, 108, 110, 112, 114, 116, 118, 120, 122, 125, 127, 129, 131, 133, 136, 138, 140,
       142, 145, 147, 149, 152, 154, 157, 159, 161, 164, 166, 169, 171, 174, 176, 179, 181, 184, 186, 189, 192, 194,
       197, 200, 203, 205, 208, 211, 214, 217, 220, 222, 225, 228, 231, 234, 237, 240, 244, 247},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   3,   5,   7,   9,   11,  13,  15,  17,  19,  21,  23,  25,
       27,  28,  30,  32,  34,  36,  38,  40,  42,  44,  46,  48,  50,  52,  54,  56,  58,  61,  63,  65,  67,  69,
       71,  73,  75,  77,  79,  81,  84,  86,  88,  90,  92,  94,  97,  99,  101, 103, 106, 108, 110, 112, 115, 117,
       119, 121, 124, 126, 129, 131, 133, 136, 138, 140, 143, 145, 148, 150, 153, 155, 158, 160, 163, 165, 168, 171,
       173, 176, 179, 181, 184, 187, 189, 192, 195, 198, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233,
       236, 239, 242, 245, 248, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   0,   2,   4,   6,   9,   11,  13,  16,  18,  20,  23,  25,  27,  30,  32,  34,  37,  39,  41,  44,  46,
       49,  51,  53,  56,  58,  61,  63,  65,  68,  70,  73,  75,  78,  80,  83,  85,  88,  90,  93,  95,  98,  100,
       103, 106, 108, 111, 113, 116, 119, 121, 124, 127, 129, 132, 135, 137, 140, 143, 146, 148, 151, 154, 157, 159,
       162, 165, 168, 171, 174, 177, 180, 182, 185, 188, 191, 194, 197, 200, 204, 207, 210, 213, 216, 219, 222, 225,
       229, 232, 235, 238, 242, 245, 248, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {3,   6,   8,   11,  13,  16,  18,  21,  24,  26,  29,  31,  34,  36,  39,  41,  44,  47,  49,  52,  54,  57,
       60,  62,  65,  67,  70,  73,  75,  78,  81,  83,  86,  89,  91,  94,  97,  100, 102, 105, 108, 111, 113, 116,
       119, 122, 125, 128, 130, 133, 136, 139, 142, 145, 148, 151, 154, 157, 160, 163, 166, 169, 172, 175, 178, 181,
       184, 187, 190, 193, 196, 200, 203, 206, 209, 212, 216, 219, 222, 226, 229, 232, 236, 239, 242, 246, 249, 253,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   0,   1,   3,   4,   5,   7,   8,   10,  11,  12,  14,  15,  17,  18,  20,  21,  23,  24,  25,  27,  28,
       30,  31,  33,  34,  36,  37,  39,  40,  42,  43,  45,  47,  48,  50,  51,  53,  54,  56,  58,  59,  61,  62,
       64,  66,  67,  69,  71,  72,  74,  76,  77,  79,  81,  83,  84,  86,  88,  90,  91,  93,  95,  97,  99,  101,
       102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 145,
       147, 149, 151, 153, 156, 158, 160, 163, 165, 167, 170, 172, 175, 177, 179, 182, 184, 187},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   3,   5,   7,   9,   11,  13,  15,  17,  19,  21,  23,  25,
       27,  28,  30,  32,  34,  36,  38,  40,  42,  44,  46,  48,  50,  52,  54,  56,  58,  61,  63,  65,  67,  69,
       71,  73,  75,  77,  79,  81,  84,  86,  88,  90,  92,  94,  97,  99,  101, 103, 106, 108, 110, 112, 115, 117,
       119, 121, 124, 126, 129, 131, 133, 136, 138, 140, 143, 145, 148, 150, 153, 155, 158, 160, 163, 165, 168, 171,
       173, 176, 179, 181, 184, 187, 189, 192, 195, 198, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233,
       236, 239, 242, 245, 248, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255},
      {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   0,   0,   0,   0,   1,   3,   4,   5,   7,   8,   9,   11,  12,  13,  15,  16,  18,  19,  20,  22,  23,
       25,  26,  27,  29,  30,  32,  33,  34,  36,  37,  39,  40,  42,  43,  45,  46,  48,  49,  51,  52,  54,  55,
       57,  58,  60,  62,  63,  65,  66,  68,  70,  71,  73,  74,  76,  78,  79,  81,  83,  85,  86,  88,  90,  91,
       93,  95,  97,  99,  100, 102, 104, 106, 108, 110, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133,
       135, 137, 139, 141, 144, 146, 148, 150, 152, 155, 157, 159, 161, 164, 166, 168, 171, 173},
  };

  // ===========================================================================
  // Calculate rho medians
  // ---------------------------------------------------------------------------
  const auto& model = context_.model_;
  const auto& model_thmc = model.theta_medians_;

  std::vector<seg_rho_t> rho_medians;

  for (const auto& method : model_thmc) {  // Begin loop methods

    std::vector<seg_rho_t> medians;

    for (const auto& group : method) {  // Begin loop groups

      std::vector<seg_rho_t> rhos;

      for (const auto& model_thm_site : group) {  // Begin loop sites
        int site_id = static_cast<int>(model_thm_site.id);

        const auto& site_bit = track.site_mask[site_id];

        // Initialize as invalid rho
        auto& rho = rhos.emplace_back(invalid_rho);

        // Short-Circuit: If no segment was selected, move on.
        if (site_bit == 0)
          continue;

        // Get selected segment's rho value
        const auto& site_seg_id = track.site_segs[site_id];
        const auto& site_seg = segments[site_seg_id];

        if (model_thm_site.theta_id == theta_id_t::kTheta1) {
          rho = site_rho_lut[site_id][site_seg.theta1];
        } else if (model_thm_site.theta_id == theta_id_t::kTheta2) {
          rho = site_rho_lut[site_id][site_seg.theta2];
        }

        // If the segment rho is 0 this is invalid rho value
        if (rho == 0) {
          rho = invalid_rho;
        }
      }  // End loop sites

      // Calculate rho median
      if (this->context_.config_.verbosity_ > 2) {
        for (const auto& rho : rhos) {
          edm::LogInfo("L1TEMTFpp") << "rho " << rho << std::endl;
        }
      }

      auto median = calc_theta_median(rhos);
      medians.push_back(median);

      if (this->context_.config_.verbosity_ > 2) {
        edm::LogInfo("L1TEMTFpp") << "group_median " << median << std::endl;
      }
    }  // End loop groups

    // Calculate method median
    auto median = calc_theta_median(medians);
    rho_medians.push_back(median);

    if (this->context_.config_.verbosity_ > 2) {
      edm::LogInfo("L1TEMTFpp") << "rho_median " << median << std::endl;
    }
  }  // End loop methods

  // ===========================================================================
  // Select track rho
  // ---------------------------------------------------------------------------
  seg_rho_t trk_abs_rho;

  if (algo == algo_id_t::kBeamHalo || track.zone != 2) {
    trk_abs_rho = rho_medians[0];
  } else {
    trk_abs_rho = rho_medians[1];
  }

  // If median is invalid, try station 1 median
  if (trk_abs_rho == invalid_rho) {
    trk_abs_rho = rho_medians[2];
  }

  // If all medians are invalid use 0 (0 is an invalid rho)
  if (trk_abs_rho == invalid_rho) {
    trk_abs_rho = 0;
  }

  // ===========================================================================
  // Store Data
  // ---------------------------------------------------------------------------
  track.rho = trk_abs_rho;
  track.valid = 1;

  // Collect track features
  track.features[37] = track.quality > 0 ? trk_abs_rho : decltype(trk_abs_rho)(0);

  // ===========================================================================
  // Compare segment rho to track rho
  // ---------------------------------------------------------------------------

  // Init segment site rho
  std::array<seg_rho_t, v3::kNumTrackSites> trk_seg_rho;

  for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
    trk_seg_rho[site_id] = 0;
  }

  // Get site rho windows
  std::vector<seg_rho_t> site_rho_window = {0, 24, 22, 14, 11, 28, 31, 30, 22, 0, 7, 0};

  // Detach segments if rho_window < diff, it is invalid
  for (unsigned int site_id = 0; site_id < v3::kNumTrackSites; ++site_id) {
    auto& site_bit = track.site_mask[site_id];
    auto& site_rm_bit = track.site_rm_mask[site_id];

    // Get Rho Window
    const auto& rho_window = site_rho_window[site_id];

    // Short-Circuit: If no segment was selected, move on.
    if (site_bit == 0)
      continue;

    const auto& site_seg_id = track.site_segs[site_id];
    const auto& site_seg = segments[site_seg_id];

    // Init differences with out-of-bounds values
    seg_rho_t diff_1 = rho_window + 1;
    seg_rho_t diff_2 = rho_window + 1;

    // Calculate abs rho 1 diff
    if (site_seg.theta1 != 0) {
      seg_rho_t seg_rho = site_rho_lut[site_id][site_seg.theta1];

      if (seg_rho < trk_abs_rho) {
        diff_1 = trk_abs_rho - seg_rho;
      } else {
        diff_1 = seg_rho - trk_abs_rho;
      }
    }

    // Calculate abs theta 2 diff
    if (site_seg.theta2 != 0) {
      seg_rho_t seg_rho = site_rho_lut[site_id][site_seg.theta2];

      if (seg_rho < trk_abs_rho) {
        diff_2 = trk_abs_rho - seg_rho;
      } else {
        diff_2 = seg_rho - trk_abs_rho;
      }
    }

    // Select the rho with the smallest difference
    if (diff_1 <= diff_2 && diff_1 < rho_window) {
      // Select rho 1 as the correct rho value
      trk_seg_rho[site_id] = site_rho_lut[site_id][site_seg.theta1];
    } else if (diff_2 < rho_window) {
      // Select rho 2 as the correct rho value
      trk_seg_rho[site_id] = site_rho_lut[site_id][site_seg.theta2];
    } else {
      // Invalidate site if both differences are outside of the rho window
      site_bit = 0;
      site_rm_bit = 1;

      // Debug Info
      if (this->context_.config_.verbosity_ > 4) {
        edm::LogInfo("L1TEMTFpp") << "Segment outside of rho window; detatched:"
                                  << " site_id " << site_id << " seg_id " << site_seg_id << " seg_phi " << site_seg.phi
                                  << " diff_1 " << diff_1 << " diff_2 " << diff_2 << std::endl;
      }
    }
  }

  // ===========================================================================
  // Calculate Direction
  // ---------------------------------------------------------------------------
  ap_uint<8> st1_word = 0;
  if (track.site_mask[1] == 1) {
    st1_word = segments[track.site_segs[1]].bx;
    st1_word += 1;
  } else if (track.site_mask[5] == 1) {
    st1_word = segments[track.site_segs[5]].bx;
    st1_word += 1;
  } 

  ap_uint<8> st2_word = 0;
  if (track.site_mask[2] == 1) {
    st2_word = segments[track.site_segs[2]].bx;
    st2_word += 1;
  } else if (track.site_mask[6] == 1) {
    st2_word = segments[track.site_segs[6]].bx;
    st2_word += 1;
  } else if (track.site_mask[10] == 1) {
    st2_word = segments[track.site_segs[10]].bx;
    st2_word += 1;
  } 

  ap_uint<8> st3_word = 0;
  if (track.site_mask[3] == 1) {
    st3_word = segments[track.site_segs[3]].bx;
    st3_word += 1;
  } else if (track.site_mask[7] == 1) {
    st3_word = segments[track.site_segs[7]].bx;
    st3_word += 1;
  } 

  ap_uint<8> st4_word = 0;
  if (track.site_mask[4] == 1) {
    st4_word = segments[track.site_segs[4]].bx;
    st4_word += 1;
  } else if (track.site_mask[8] == 1) {
    st4_word = segments[track.site_segs[8]].bx;
    st4_word += 1;
  } 

  ap_uint<8> bx_word = st1_word;
  bx_word += st2_word << 2;
  bx_word += st3_word << 4;
  bx_word += st4_word << 6;

  std::vector<trk_direction_t> direction_lut = {
      0,0,0,0,0,1,2,2,0,3,1,2,0,3,3,1,0,1,2,2,1,1,2,2,2,0,2,2,2,0,0,2,
      0,3,1,2,3,3,0,0,1,3,1,2,2,0,0,2,0,3,3,1,3,3,0,0,3,3,3,0,1,3,3,1,
      0,1,2,2,1,1,2,2,2,0,2,2,2,0,0,2,1,1,2,2,1,1,2,2,2,0,2,2,2,0,0,2,
      2,0,2,2,0,0,0,0,2,0,2,2,2,0,0,2,2,0,0,2,0,0,0,0,0,0,0,0,2,0,0,2,
      0,3,1,2,3,3,0,0,1,3,1,2,2,0,0,2,3,3,0,0,3,3,0,0,0,0,0,0,0,0,0,0,
      1,3,1,2,3,3,0,0,1,3,1,2,2,0,0,2,2,0,0,2,0,0,0,0,0,0,0,0,2,0,0,2,
      0,3,3,1,3,3,0,0,3,3,3,0,1,3,3,1,3,3,0,0,3,3,0,0,0,0,0,0,0,0,0,0,
      3,3,3,0,3,3,0,0,3,3,3,0,0,0,0,0,1,3,3,1,3,3,0,0,3,3,3,0,1,3,3,1
  };

  track.direction = direction_lut[bx_word];

  // ===========================================================================
  // Collect Segment Features
  // ---------------------------------------------------------------------------
  const auto& model_ftc = model.features_;

  int i_feature = 0;

  for (auto& model_ft : model_ftc) {
    for (auto& model_ft_site : model_ft.sites) {
      int site_id = static_cast<int>(model_ft_site);

      const auto& site_seg_id = track.site_segs[site_id];
      const auto& site_bit = track.site_mask[site_id];
      const auto& site_seg = segments[site_seg_id];

      auto& trk_feature = track.features[i_feature++];

      // Short-Circuit: No segment attached
      if (site_bit == 0) {
        continue;
      }

      // Fill features
      if (model_ft.id == feature_id_t::kPhi) {
        // Note: This is the segment's phi with respect to the track's abs phi
        trk_feature = static_cast<trk_feature_t>(site_seg.phi) - static_cast<trk_feature_t>(track.phi);
      } else if (model_ft.id == feature_id_t::kTheta) {
        // Note: This is the segment's rho with respect to the track's abs rho
        trk_feature = static_cast<trk_feature_t>(trk_seg_rho[site_id]) - static_cast<trk_feature_t>(track.rho);
      } else if (model_ft.id == feature_id_t::kBend) {
        trk_feature = site_seg.bend;
      } else if (model_ft.id == feature_id_t::kQuality) {
        trk_feature = site_seg.qual1;
      }
    }
  }

  // Debug Info
  if (this->context_.config_.verbosity_ > 1) {
    edm::LogInfo("L1TEMTFpp") << "Track"
                              << " zone " << track.zone << " col " << track.col << " pat " << track.pattern << " qual "
                              << track.quality << " abs_phi " << track.phi << " abs_rho " << trk_abs_rho
                              << " features:" << std::endl;

    for (unsigned int i = 0; i < v3::kNumTrackFeatures; ++i) {
      if (i > 0) {
        edm::LogInfo("L1TEMTFpp") << " ";
      }

      edm::LogInfo("L1TEMTFpp") << track.features[i];
    }

    edm::LogInfo("L1TEMTFpp") << std::endl;
  }
}
