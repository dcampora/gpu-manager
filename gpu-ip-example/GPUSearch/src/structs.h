#ifndef STRUCTS_H_
#define STRUCTS_H_

struct hitData  { //49 B -  aligned to 56
	int id;
	int trackId;
	int sensor;
	int inputHit;
	double x;
	double y;
	double z;
	double wx;
	bool isUsed;
};

struct sensorInfo{
	int startPosition;
	int hitsNum;
	double z;
};

struct track {
	double m_x0;
	double m_tx;
	double m_y0;
	double m_ty;

	double m_s0;
	double m_sx;
	double m_sz;
	double m_sxz;
	double m_sz2;

	double m_u0;
	double m_uy;
	double m_uz;
	double m_uyz;
	double m_uz2;

	int    internalId;
	int    trackHitsNum; 	//number of all hits for this track
	int	   firstHit;		// location of first hit for this track
	bool   m_backward;
};


#endif /* STRUCTS_H_ */
