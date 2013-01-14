// $Id: $
#ifndef PATPIXELHITMANAGER_H
#define PATPIXELHITMANAGER_H 1

// Include files
// from Gaudi
#include "PatPixelHit.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "VeloPixDet/DeVeloPix.h"
#include "Event/VeloPixLiteCluster.h"

#include "PatPixelSensor.h"

static const InterfaceID IID_PatPixelHitManager ( "PatPixelHitManager", 1, 0 );
typedef std::vector<std::vector<int> > Int2DimVector;
typedef std::vector<std::vector<double> > Double2DimVector;
/** @class PatPixelHitManager PatPixelHitManager.h
 *  Tool to handle the Pixel velo geometry and hits, from FastVelo
 *
 *  @author Olivier Callot
 *  @date   2012-01-05
 */
class PatPixelHitManager : public GaudiTool, public IIncidentListener {
public:

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_PatPixelHitManager; }

  /// Standard constructor
  PatPixelHitManager( const std::string& type,
                      const std::string& name,
                      const IInterface* parent);

  virtual ~PatPixelHitManager( ); ///< Destructor

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  void buildHits();

  void clearHits();

  void handle ( const Incident& incident );

  PatPixelHits& hits( unsigned int sensor ) { return m_sensors[sensor]->hits(); }

  PatPixelSensor* sensor( unsigned int n ) { return m_sensors[n]; }

  unsigned int firstSensor()   const { return m_firstSensor;   }
  unsigned int lastSensor()    const { return m_lastSensor;    }

  int nbHits()  const { return m_nextInPool - m_pool.begin(); }
  int maxSize() const { return m_maxSize; }

  StatusCode rebuildGeometry();  ///< Recompute the geometry in case of change

  void sortByX();

  std::vector<int>&		get_hitsPerEvent()	{ return m_hitsPerEvent; }
  Int2DimVector& 		get_hit_id()		{ return m_hits_Ids; }
  Int2DimVector& 		get_hit_trackId()	{ return m_hits_trackIds; }
  Int2DimVector& 		get_hit_sensorNum()	{ return m_hits_sensNum; }
  Int2DimVector& 		get_hit_isUsed()	{ return m_hits_isUsed; }
  Double2DimVector& 	get_hit_X()			{ return m_hits_X; }
  Double2DimVector& 	get_hit_Y()			{ return m_hits_Y; }
  Double2DimVector& 	get_hit_Z()			{ return m_hits_Z; }
  Double2DimVector& 	get_hit_W()			{ return m_hits_W; }
  Int2DimVector&		get_sensor_hitStartPos()	{return m_sensor_hitStartPos; }
  Int2DimVector&		get_sensor_hitsNum()		{return m_sensor_hitsNum; }
  Double2DimVector&		get_sensor_Z()				{return m_sensor_Z; }

  void clear_vectors();
  int	get_max_hits()  {return m_max_hits; }

protected:

private:
  //== Containers for hits info

  std::vector<int> m_hitsPerEvent;
  Int2DimVector m_hits_Ids;
  Int2DimVector m_hits_trackIds;
  Int2DimVector m_hits_isUsed;
  Int2DimVector m_hits_sensNum;
  Double2DimVector m_hits_X;
  Double2DimVector m_hits_Y;
  Double2DimVector m_hits_Z;
  Double2DimVector m_hits_W;
  int m_event_number;
  int m_max_hits;
  Int2DimVector m_sensor_hitStartPos;
  Int2DimVector m_sensor_hitsNum;
  Double2DimVector m_sensor_Z;

  DeVeloPix* m_veloPix;

  std::vector<PatPixelHit>  m_pool;
  std::vector<PatPixelHit>::iterator m_nextInPool;
  std::vector<PatPixelSensor*> m_sensors;
  unsigned int m_firstSensor;
  unsigned int m_lastSensor;
  int m_maxSize;
  bool m_eventReady;
};
#endif // PATPIXELHITMANAGER_H
