#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFConstants.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFContext.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFTypes.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Algo/HitmapLayer.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/CSCTPConverter.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/CSCTPSelector.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/RPCTPConverter.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/RPCTPSelector.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/GEMTPConverter.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/GEMTPSelector.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/ME0TPConverter.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/ME0TPSelector.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/GE0TPConverter.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/GE0TPSelector.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/TPConverters.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/TPSelectors.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/DAQ/TPrimitives.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/CSCUtils.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/SectorProcessor.h"

using namespace emtf::phase2;

SectorProcessor::SectorProcessor(const EMTFContext& context, const int& endcap, const int& sector)
    : context_(context), endcap_(endcap), sector_(sector), event_(nullptr), proc_bx_(-999) {
  // ===========================================================================
  // Register Selectors/Converters
  // ===========================================================================
  if (this->context_.config_.csc_en_) {
    tp_selectors_[L1TMuon::kCSC] = std::make_unique<CSCTPSelector>(context_, endcap_, sector_);
    tp_converters_[L1TMuon::kCSC] = std::make_unique<CSCTPConverter>(context_, endcap_, sector_);
  }

  if (this->context_.config_.rpc_en_) {
    tp_selectors_[L1TMuon::kRPC] = std::make_unique<RPCTPSelector>(context_, endcap_, sector_);
    tp_converters_[L1TMuon::kRPC] = std::make_unique<RPCTPConverter>(context_, endcap_, sector_);
  }

  if (this->context_.config_.gem_en_) {
    tp_selectors_[L1TMuon::kGEM] = std::make_unique<GEMTPSelector>(context_, endcap_, sector_);
    tp_converters_[L1TMuon::kGEM] = std::make_unique<GEMTPConverter>(context_, endcap_, sector_);
  }

  if (this->context_.config_.me0_en_) {
    tp_selectors_[L1TMuon::kME0] = std::make_unique<ME0TPSelector>(context_, endcap_, sector_);
    tp_converters_[L1TMuon::kME0] = std::make_unique<ME0TPConverter>(context_, endcap_, sector_);
  }

  if (this->context_.config_.ge0_en_) {
    tp_selectors_[L1TMuon::kME0] = std::make_unique<GE0TPSelector>(context_, endcap_, sector_);
    tp_converters_[L1TMuon::kME0] = std::make_unique<GE0TPConverter>(context_, endcap_, sector_);
  }
}

SectorProcessor::~SectorProcessor() {
  // Do Nothing
}

void SectorProcessor::configureEvent(const edm::Event& event) {
  // Event
  event_ = &event;
  proc_bx_ = -999;

  // Reset Window Hits
  hcol_buffer_.clear();
  subsystem_ilink_tpc_maps_.clear();
}

void SectorProcessor::configureBx(const int& input_bx) {
  // Reset BX Maps
  subsystem_ilink_tpc_maps_.clear();

  // Remove BX TPCollections that aren't in the bx window
  // Note that the first entry in hcol_buffer_ is the earliest BX
  const auto min_bx = this->context_.config_.min_bx_;
  const auto delay_bx = this->context_.config_.bx_window_ - 1;
  const auto start_bx = min_bx + delay_bx;

  if (start_bx < input_bx) {
    hcol_buffer_.erase(hcol_buffer_.begin());
  }

  // BX
  proc_bx_ = input_bx - delay_bx;
}

void SectorProcessor::select(const TriggerPrimitive& tp, const TPInfo& tp_info) {
  // Get TPSelector
  auto tp_subsystem = tp.subsystem();
  auto tp_selectors_it = tp_selectors_.find(tp_subsystem);

  // Short-Circuit: Operation not supported
  if (tp_selectors_it == tp_selectors_.end()) {
    edm::LogWarning("L1TEMTFpp") << "TPCollector has been implemented, "
                                 << "but there is no TPSelector for " << tp_subsystem;
    return;
  }

  // Select TP that belongs to this Sector Processor
  auto& subsystem_ilink_tpc_map = subsystem_ilink_tpc_maps_[tp_subsystem];  // reference to subsystem trigger primitive collection
  tp_selectors_it->second->select(tp, tp_info, subsystem_ilink_tpc_map);
}

void SectorProcessor::process(EMTFHitCollection& out_hits,
                              EMTFTrackCollection& out_tracks,
                              EMTFInputCollection& out_inputs) {
  // ===========================================================================
  // Merge subsystem selections
  // ===========================================================================
  ILinkTPCMap bx_ilink_tpc_map;

  for (auto& [subsystem, ilink_tpc_map] : subsystem_ilink_tpc_maps_) {
    copyTP(ilink_tpc_map, bx_ilink_tpc_map);
  }

  // Free memory
  subsystem_ilink_tpc_maps_.clear();

  // ===========================================================================
  // Convert trigger primitives to EMTF Hits
  // ===========================================================================

  // Convert tp into hits
  EMTFHitCollection bx_hcol;

  convertTP(out_hits.size(), bx_ilink_tpc_map, bx_hcol);

  // Append to bx window hits
  hcol_buffer_.push_back(bx_hcol);

  // Free memory
  bx_ilink_tpc_map.clear();

  // Record hits
  out_hits.insert(out_hits.end(), bx_hcol.begin(), bx_hcol.end());

  // ===========================================================================
  // Process tracks between MIN_BX and MAX_BX
  // ===========================================================================
  const auto min_bx = this->context_.config_.min_bx_;
  const auto max_bx = this->context_.config_.max_bx_;

  if (!((min_bx <= proc_bx_) && (proc_bx_ <= max_bx))) {
    return;
  }

  // ===========================================================================
  // Convert EMTF Hits to Segments
  // ===========================================================================

  // Init Segment to Hit Map
  std::map<int, int> seg_to_hit;

  // Convert bx window hits into segments
  segment_collection_t segments;

  populateSegments(hcol_buffer_, seg_to_hit, segments);

  // Build Tracks
  buildTracks(seg_to_hit, segments, false, out_tracks);  // With prompt setup
  buildTracks(seg_to_hit, segments, true, out_tracks);   // With displaced setup

  // ===========================================================================
  // Record segments/hits used in track building
  // ===========================================================================
  if (!seg_to_hit.empty()) {
    EMTFInput::hits_t hit_id_col;
    EMTFInput::segs_t seg_id_col;

    for (const auto& [seg_id, hit_id] : seg_to_hit) {
      seg_id_col.push_back(seg_id);
      hit_id_col.push_back(hit_id);
    }

    EMTFInput emtf_input;

    const int endcap_pm = (endcap_ == 2) ? -1 : endcap_;  // 1: +endcap, -1: -endcap

    emtf_input.setEndcap(endcap_pm);
    emtf_input.setSector(sector_);
    emtf_input.setBx(proc_bx_);
    emtf_input.setHits(hit_id_col);
    emtf_input.setSegs(seg_id_col);

    out_inputs.push_back(emtf_input);
  }
}

void SectorProcessor::copyTP(const ILinkTPCMap& source, ILinkTPCMap& target) const {
  typedef typename ILinkTPCMap::iterator Iterator_t;
  typedef typename ILinkTPCMap::mapped_type Collection_t;

  for (auto& source_kv : source) {
    std::pair<Iterator_t, bool> ins_res = target.insert(source_kv);

    // Short-Circuit: Insertion succeeded, move on
    if (ins_res.second) {
      continue;
    }

    // Merge into target collection
    const Collection_t& source_col = source_kv.second;
    Collection_t& target_col = ins_res.first->second;

    target_col.insert(target_col.end(), source_col.begin(), source_col.end());
  }
}

void SectorProcessor::convertTP(const int& initial_hit_id, const ILinkTPCMap& ilink_tpc_map, EMTFHitCollection& hits) {
  EMTFHitCollection substitutes;

  for (const auto& [ilink, ilink_tpc] : ilink_tpc_map) {  // loop input link trigger primitive collections

    unsigned int cnt_segments = 0;  // Enumerates each segment in the same chamber

    for (const auto& tp_entry : ilink_tpc) {  // loop trigger primitives

      // Unpack Entry
      const auto& tp = tp_entry.tp_;  // const reference
      auto tp_info = tp_entry.info_;  // copy info

      // Unpack trigger primitive
      auto tp_subsystem = tp.subsystem();

      // Get Converter
      auto tp_converters_it = tp_converters_.find(tp_subsystem);

      // Short-Circuit: Operation not supported
      if (tp_converters_it == tp_converters_.end()) {
        edm::LogWarning("L1TEMTFpp") << "TPCollector & TPSelector have been implemented, "
                                     << "but there is no TPConverter for " << tp_subsystem;
        continue;
      }

      // Set Segment Id
      tp_info.segment_id = cnt_segments++;  // Save and increase segment count

      // Convert
      EMTFHit hit;

      tp_converters_it->second->convert(tp, tp_info, hit);

      // Append to hit collections
      if (tp_info.flag_substitute) {
        substitutes.push_back(hit);
      } else {
        hits.push_back(hit);
      }
    }
  }

  // Substitutes are placed at the end of the hit collection
  hits.insert(hits.end(), std::make_move_iterator(substitutes.begin()), std::make_move_iterator(substitutes.end()));

  // Assign Hit Ids
  unsigned int cnt_hits = initial_hit_id;

  for (auto& hit : hits) {
    hit.setId(cnt_hits++);
  }
}

void SectorProcessor::populateSegments(const std::vector<EMTFHitCollection>& bx_window_hcol_buffer,
                                       std::map<int, int>& seg_to_hit,
                                       segment_collection_t& segments) {
  // Initialize
  for (unsigned int seg_id = 0; seg_id < v3::kNumSegments; ++seg_id) {
    segments[seg_id].phi = 0;
    segments[seg_id].bend = 0;
    segments[seg_id].theta1 = 0;
    segments[seg_id].theta2 = 0;
    segments[seg_id].qual1 = 0;
    segments[seg_id].qual2 = 0;
    segments[seg_id].time = 0;
    segments[seg_id].zones = 0;
    segments[seg_id].tzones = 0;
    segments[seg_id].cscfr = 0;
    segments[seg_id].layer = 0;
    segments[seg_id].bx = 0;
    segments[seg_id].valid = 0;
  }

  // Populate
  // Loop hit collections from earliest to latest BX
  std::map<int, unsigned int> bx_window_ch_seg;

  for (const auto& bx_hcol : bx_window_hcol_buffer) {
    std::map<int, unsigned int> bx_ch_seg;

    for (const auto& hit : bx_hcol) {  // Begin loop hits in BX
      // Unpack Hit
      const auto& hit_chamber = hit.emtfChamber();
      const auto& hit_segment = hit.emtfSegment();
      const auto& hit_valid = hit.flagValid();

      emtf_assert(hit_valid);  // segment must be valid

      // Get Channel Segment Count
      unsigned int ch_seg = bx_window_ch_seg[hit_chamber] + hit_segment;

      // Update bx chamber last segment
      bx_ch_seg[hit_chamber] = ch_seg;

      // Short-Circuit: Accept at most 2 segments
      if (!(ch_seg < v3::kChamberSegments)) {
        continue;
      }

      // Calculate Host
      const auto& hit_host = hit.emtfHost();

      // Calculate Timezone
      // Note: Uses Hit BX relative to Sector Processor BX
      const auto& hit_bx = hit.bx();
      const int hit_rel_bx = (hit_bx - proc_bx_);
      const auto hit_timezones = context_.timezone_lut_.getTimezones(hit_host, hit_rel_bx);

      // Calculate algo seg
      const unsigned int seg_id = hit_chamber * v3::kChamberSegments + ch_seg;

      emtf_assert(seg_id < v3::kNumSegments);

      seg_to_hit[seg_id] = hit.id();

      // Populate segment
      segments[seg_id].phi = hit.emtfPhi();
      segments[seg_id].bend = hit.emtfBend();
      segments[seg_id].theta1 = hit.emtfTheta1();
      segments[seg_id].theta2 = hit.emtfTheta2();
      segments[seg_id].qual1 = hit.emtfQual1();
      segments[seg_id].qual2 = hit.emtfQual2();
      segments[seg_id].time = hit.emtfTime();
      segments[seg_id].zones = hit.emtfZones();
      segments[seg_id].tzones = hit_timezones;
      segments[seg_id].cscfr = hit.cscFR();
      segments[seg_id].layer = hit.layer();
      segments[seg_id].bx = hit.bx();
      segments[seg_id].valid = hit.flagValid();

      // Debug Info
      if (this->context_.config_.verbosity_ > 1) {
        edm::LogInfo("L1TEMTFpp") << std::endl
                                  << "Event: " << event_->id() << " Endcap: " << endcap_ << " Sector: " << sector_
                                  << " BX: " << (proc_bx_) << " Hit iLink: " << hit_chamber << " Hit iSeg: " << ch_seg
                                  << " Hit Host " << hit_host << " Hit Rel BX " << (hit_bx - proc_bx_)
                                  << " Hit Timezones " << hit_timezones << std::endl;

        edm::LogInfo("L1TEMTFpp") << " id " << seg_id << " phi " << segments[seg_id].phi << " bend "
                                  << segments[seg_id].bend << " theta1 " << segments[seg_id].theta1 << " theta2 "
                                  << segments[seg_id].theta2 << " qual1 " << segments[seg_id].qual1 << " qual2 "
                                  << segments[seg_id].qual2 << " time " << segments[seg_id].time << " zones "
                                  << segments[seg_id].zones << " timezones " << segments[seg_id].tzones << " cscfr "
                                  << segments[seg_id].cscfr << " layer " << segments[seg_id].layer << " bx "
                                  << segments[seg_id].bx << " valid " << segments[seg_id].valid << std::endl;
      }
    }  // End loop hits from BX

    for (auto& [chamber, ch_seg] : bx_ch_seg) {
      bx_window_ch_seg[chamber] = ch_seg + 1;
    }
  }  // End loop from earliest BX Collection to latest BX Hit Collection
}

void SectorProcessor::buildTracks(const std::map<int, int>& seg_to_hit,
                                  const segment_collection_t& segments,
                                  const bool& displaced_en,
                                  EMTFTrackCollection& out_tracks) {
  // Apply Hitmap Building Layer: Convert segments into hitmaps
  std::vector<hitmap_t> zone_hitmaps;

  context_.hitmap_building_layer_.apply(segments, zone_hitmaps);

  // Apply Pattern Matching Layer: Match patterns to hitmaps to create roads
  std::vector<road_collection_t> zone_roads;

  context_.pattern_matching_layer_.apply(zone_hitmaps, displaced_en, zone_roads);

  // Apply Road Sorting Layer: Find the best roads
  std::vector<road_t> best_roads;

  context_.road_sorting_layer_.apply(v3::kNumTracks, zone_roads, best_roads);

  // Apply Track Building Layer: Match segments to the best roads to create tracks
  std::vector<track_t> tracks;

  context_.track_building_layer_.apply(segments, best_roads, displaced_en, tracks);

  // Apply Duplicate Removal Layer: Removes tracks that share a segment, keeping the one that has the highest quality
  context_.duplicate_removal_layer_.apply(tracks);

  // Apply Parameter Assigment Layer: Run NN on tracks
  context_.parameter_assignment_layer_.apply(displaced_en, tracks);

  // Apply Output Layer
  EMTFTrackCollection bx_tracks;

  context_.output_layer_.apply(endcap_, sector_, proc_bx_, seg_to_hit, tracks, displaced_en, bx_tracks);

  // Record tracks
  out_tracks.insert(out_tracks.end(), bx_tracks.begin(), bx_tracks.end());
}
