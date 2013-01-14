// $Id: $
#ifndef PATPIXELTRACKING_H
#define PATPIXELTRACKING_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "PatKernel/IPatDebugTool.h"
#include "GaudiAlg/ISequencerTimerTool.h"
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/date_time.hpp>
#include <algorithm>

#include "PatPixelHitManager.h"
#include "PatPixelTrack.h"

/** @class PatPixelTracking PatPixelTracking.h
 *  This is the main tracking fo rthe Velo Pixel upgrade
 *
 *  @author Olivier Callot
 *  @date   2011-12-16
 */

 struct shared_conditions {
	  shared_conditions()
          :  package_in(false)
          {}
      boost::interprocess::interprocess_mutex      mutex;
      boost::interprocess::interprocess_condition  package_full;
      boost::interprocess::interprocess_condition  package_empty;
      bool package_in;
 };
typedef boost::interprocess::allocator<char, boost::interprocess::managed_shared_memory::segment_manager>  ShmemAllocator;
typedef boost::interprocess::vector<char, ShmemAllocator> SharedVec;

class PatPixelTracking : public GaudiAlgorithm {
public:
  /// Standard constructor
  PatPixelTracking( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~PatPixelTracking( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

  void searchByPair();

  void makeLHCbTracks( LHCb::Tracks* tracks );

  //== Debugging methods
  bool matchKey( const PatPixelHit* hit ) {
    if ( m_debugTool ) {
      LHCb::LHCbID id = hit->id();
      return m_debugTool->matchKey( id, m_wantedKey );
    }
    return false;
  }

  void removeWorstHit ( double maxChi2 );
  bool addHitsOnSensor( PatPixelSensor* sensor, double xTol, double maxChi2 );
  void printHit( const PatPixelHit* hit, std::string title="" );
  void printTrack( PatPixelTrack& track );
  void printHitOnTrack ( PatPixelHit* hit, bool ifMatch=true );
  template<class T> inline void AddToVector(T element, SharedVec & vec);
  void CombineVectors(std::vector<char> & output);
  template<class T> void combine(SharedVec & output, std::vector<std::vector<T> > &input, unsigned int length);
  void combineAllVectors(SharedVec &outputVector);
  void launchGPU(unsigned int segmentSize);
  void closeGPUmemory();
  void openGPUMemory();

private:
  std::string m_outputLocation;
  PatPixelHitManager* m_hitManager;

  double m_maxXSlope;
  double m_maxYSlope;
  double m_maxZForRBeamCut;
  double m_maxR2Beam;
  double m_extraTol;
  double m_maxChi2ToAdd;
  double m_maxChi2SameSensor;
  int    m_maxMissed;
  double m_maxChi2PerHit;
  double m_maxChi2Short;

  PatPixelTracks m_tracks;
  PatPixelTrack  m_track;

  //== Debug control
  std::string      m_debugToolName;
  int              m_wantedKey;
  IPatDebugTool*   m_debugTool;
  bool             m_isDebug;
  bool             m_debug;

  //== Timing measurement control
  bool             m_doTiming;
  ISequencerTimerTool* m_timerTool;
  int   m_timeTotal;
  int   m_timePrepare;
  int   m_timePairs;
  int   m_timeFinal;
  shared_conditions * cond;

  SharedVec *m_output_vec;

};
#endif // PATPIXELTRACKING_H
