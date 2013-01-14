#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "SearchKernel.cuh"
using namespace std;

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }


typedef struct __align__(16)  {
	int startPosition;
	int hitsNum;
	double z;
} sensorInfo;

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

	int internalId;
	int trackHitsNum; // hmm
	int firstHit; 	 //	hmm

	bool m_backward;
};

__device__ __constant__ double 	m_maxXSlope			= 0.400;
__device__ __constant__ double 	m_maxYSlope			= 0.300;
__device__ __constant__ double 	m_maxZForRBeamCut	= 200.0;
__device__ __constant__ double 	m_maxR2Beam			= 1.0 ;
__device__ __constant__ int 	m_maxMissed			= 4;
__device__ __constant__ double 	m_extraTol			= 0.150;
__device__ __constant__ double 	m_maxChi2ToAdd		= 100.0;
__device__ __constant__ double 	m_maxChi2SameSensor	= 16.0;
__device__ __constant__ double  m_maxChi2Short		= 6.0 ;
__device__ __constant__ double  m_maxChi2PerHit		= 16.0;
__device__ __constant__ int 	m_sensNum			= 48;


__device__ double zBeam(track *tr) {
	return -( tr->m_x0 * tr->m_tx + tr->m_y0 * tr->m_ty ) / ( tr->m_tx * tr->m_tx + tr->m_ty * tr->m_ty );
}

__device__ double r2AtZ( double z , track *tr) {
    double xx = tr->m_x0 + z * tr->m_tx;
    double yy = tr->m_y0 + z * tr->m_ty;
    return xx*xx + yy * yy;
 }

__device__ void solve (track *tr) {
	double den = ( tr->m_sz2 * tr->m_s0 - tr->m_sz * tr->m_sz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_tx     = ( tr->m_sxz * tr->m_s0  - tr->m_sx  * tr->m_sz ) / den;
	tr->m_x0     = ( tr->m_sx  * tr->m_sz2 - tr->m_sxz * tr->m_sz ) / den;

	den = ( tr->m_uz2 * tr->m_u0 - tr->m_uz * tr->m_uz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_ty     = ( tr->m_uyz * tr->m_u0  - tr->m_uy  * tr->m_uz ) / den;
	tr->m_y0     = ( tr->m_uy  * tr->m_uz2 - tr->m_uyz * tr->m_uz ) / den;
}

__device__ void addHit ( track *tr, int offset, int *trackIds, double *inX, double *inY, double *inZ, double *inW) {

	trackIds[offset] = tr->internalId;
	tr->trackHitsNum++;

	double z = inZ[offset];
	double x = inX[offset];
	double w = inW[offset];

	tr->m_s0  += w;
	tr->m_sx  += w * x;
	tr->m_sz  += w * z;
	tr->m_sxz += w * x * z;
	tr->m_sz2 += w * z * z;

	double y = inY[offset];

	tr->m_u0  += w;
	tr->m_uy  += w * y;
	tr->m_uz  += w * z;
	tr->m_uyz += w * y * z;
	tr->m_uz2 += w * z * z;

	if( tr->trackHitsNum > 1 ) solve(tr);
}

__device__ void setTrack(track *tr, int hit0offset, int hit1offset, int *trackIds, double *inX, double *inY, double *inZ, double *inW){
	tr->m_backward = false;
	tr->trackHitsNum = 0;

	trackIds[hit0offset] = tr->internalId;
	tr->trackHitsNum++;

	double z = inZ[hit0offset];
	double x = inX[hit0offset];
	double w = inW[hit0offset];

	tr->m_s0  = w;
	tr->m_sx  = w * x;
	tr->m_sz  = w * z;
	tr->m_sxz = w * x * z;
	tr->m_sz2 = w * z * z;

	double y = inY[hit0offset];

	tr->m_u0  = w;
	tr->m_uy  = w * y;
	tr->m_uz  = w * z;
	tr->m_uyz = w * y * z;
	tr->m_uz2 = w * z * z;

	addHit (tr, hit1offset,  trackIds, inX, inY, inZ, inW);
}

__device__ inline double chi2Hit( double x, double y, double hitX, double hitY, double hitW){
	double dx = x - hitX;
	double dy = y - hitY;
	return dx * dx * (hitW) + dy * dy * (hitW);
}
__device__ inline double xAtHit(track *tr, double z )
{
	return tr->m_x0 + tr->m_tx * z;
}

__device__ inline double yAtHit( track *tr, double z  )
{
	return tr->m_y0 + tr->m_ty * z;
}

__device__ inline double chi2Track(track *tr, double hitX, double hitY, double hitZ, double hitW)
{
	return chi2Hit( xAtHit( tr, hitZ ), yAtHit(tr, hitZ ), hitX, hitY, hitW);
}

__device__ inline double chi2(track *tr, double *hitX, double *hitY, double *hitZ, double *hitW, int usedHits[])
{
	double ch = 0.0;
	int nDoF  = -4;
	for (int i =0 ; i<m_sensNum; i++){ // or just... while(i!=m_sensNum && usedHits[i] != -1) ,
		int hitNumber = usedHits[i];
		if(hitNumber== -1){
			break;
		}
		ch += chi2Track(tr, hitX[hitNumber], hitY[hitNumber], hitZ[hitNumber], hitW[hitNumber]);
		nDoF += 2;
	}
	return ch/nDoF;
}

__device__ inline bool addHitsOnSensor( double sensZ, int sensStartPost, int sensHitsNum, double xTol, double maxChi2, track *tr,
								 	    int *tracksIds, double *inX, double *inY, double *inZ, double *inW, int threadId,
								 	    int threadsNumber,int *usedHit ) {
	if (sensHitsNum == 0 ) return false;
	double xGuess = (tr->m_x0 + tr->m_tx * sensZ) - xTol - 1;
	int lastHit = sensStartPost + sensHitsNum -1;
	if((inX[(lastHit*threadsNumber)+threadId]) < xGuess) return false;
	int startHit = sensStartPost; // może lepiej bezpośrednio się do tego odnosić
	unsigned int step = sensHitsNum;
	while ( step > 2 ) {
		step = step/2;
		if (inX[((startHit+step)*threadsNumber)+threadId] < xGuess) startHit +=step;
	}
	bool added = false;
	int tmpOffset = 0;
	for(int iH = startHit; iH<=lastHit ;++iH){ ////
		tmpOffset = (iH*threadsNumber)+threadId;
		double xPred = (tr->m_x0 + tr->m_tx * inZ[tmpOffset]);
		if ( inX[tmpOffset] + xTol < xPred ) continue;
		if ( inX[tmpOffset] - xTol > xPred ) break;
		if ( chi2Track(tr, inX[tmpOffset], inY[tmpOffset], inZ[tmpOffset], inW[tmpOffset]) < maxChi2 ) {
			addHit(tr, tmpOffset,tracksIds, inX, inY, inZ, inW );
			*usedHit = tmpOffset;
			added = true;
		}
	}
	return added;
}

__device__ inline void removeHit(track *tr, int worstHitOffset, double *inX, double *inY, double *inZ, double *inW ){
	tr->trackHitsNum--;

	double z = inZ[worstHitOffset];
	double w = inW[worstHitOffset];
	double x = inX[worstHitOffset];

	tr->m_s0  -= w;
	tr->m_sx  -= w * x;
	tr->m_sz  -= w * z;
	tr->m_sxz -= w * x * z;
	tr->m_sz2 -= w * z * z;

	double y = inY[worstHitOffset];

	tr->m_u0  -= w;
	tr->m_uy  -= w * y;
	tr->m_uz  -= w * z;
	tr->m_uyz -= w * y * z;
	tr->m_uz2 -= w * z * z;

	if( tr->trackHitsNum > 1 ) solve(tr);
}

__device__ inline void removeWorstHit(track* tr, double maxChi2, int usedHits[], int *tracksIds, double *inX, double *inY, double *inZ,
									  double *inW, int *isUsed, int threadsNumber, int threadId)
{
	double topChi2 = 1.e9;
	while( topChi2 > maxChi2 ) {
	    topChi2 = 0.0;
	    int offset = 0;
	    for (int i =0 ; i<m_sensNum; i++){
	    	offset = usedHits[i];
	    	if(offset == -1){
	    		break;
	    	}
	    	double myChi2 = chi2Track(tr, inX[offset], inY[offset], inZ[offset], inW[offset] );
	    	if ( myChi2 > topChi2 ) {
	    		topChi2 = myChi2;
	    		offset = (usedHits[i]*threadsNumber)+threadId;
	    	}
	    }
	    if ( topChi2 > maxChi2 ) {
	      tracksIds[offset] = -1;
	      isUsed[offset] = 0;
	      removeHit(tr, offset, inX, inY, inZ, inW);
	    }
	}
}

__device__ inline bool all3SensorsAreDifferent(int usedHits[], int *sensorsIds) {
    if ( sensorsIds[usedHits[0]] == sensorsIds[usedHits[1]]) return false;
    if ( sensorsIds[usedHits[0]] == sensorsIds[usedHits[2]]) return false;
    if ( sensorsIds[usedHits[1]] == sensorsIds[usedHits[2]]) return false;
    return true;
}

__device__ inline int nbUnused(int usedHits[],int *isUsed) {
	int nn = 0;
    for (int i=0; i<m_sensNum; i++){
    	int hitNumber = usedHits[i];
    	if(hitNumber== -1){
    		break;
    	}
    	if ( ! isUsed[hitNumber] ) ++nn;
    }
    return nn;
}

/*
 * topHitsNum - number of hits in the event with biggest number of hits. Other events are filled with dummy hits to have the same number of hits.
 */
__global__ void searchByPair(void *data, void* resultsKernel, int topHitsNum){

	int threadId = blockDim.x * blockIdx.x + threadIdx.x;
	track m_track = {0};			// move to shared memory
	int trackId = 0;
	int nextHitNb = 0;
	int usedHits[48] = {-1}; //numbers of used hits

	__shared__ int 		*inHitsNum;
	__shared__ int 		*inSensInfoStartPos;
	__shared__ int 		*inSensInfoHitsNum;
	__shared__ double 	*inSensInfoZ;
	__shared__ int 		*hitsIds;
	__shared__ int 		*tracksIds;
	__shared__ int 		*sensorsIds;
	__shared__ int 		*isUsed;
	__shared__ double 	*inX;
	__shared__ double 	*inY;
	__shared__ double 	*inZ;
	__shared__ double 	*inW;

	int eventsNumber = gridDim.x * blockDim.x;			//number of all the events passed to GPU
	int dataOffset =  eventsNumber * topHitsNum; 		// between different hits parameters

	if(threadId % blockDim.x == 0){	//only first thread from the block does that - maybe not the best idea ;)
		inHitsNum 				= 		(int*)data;
		inSensInfoStartPos 		= 		inHitsNum + eventsNumber;
		inSensInfoHitsNum		= 		inSensInfoStartPos + (eventsNumber * m_sensNum);
		inSensInfoZ 			= 		(double*)(inSensInfoHitsNum + (eventsNumber * m_sensNum));
		hitsIds 				= 		(int*)(inSensInfoZ + (eventsNumber * m_sensNum));
		tracksIds 				= 		hitsIds + dataOffset;
		sensorsIds				=		tracksIds + dataOffset;
		isUsed					=		sensorsIds + dataOffset;
		inX 					= 		(double*) (isUsed + dataOffset);
		inY 					= 		inX + dataOffset;
		inZ						=		inY + dataOffset;
		inW						=		inZ + dataOffset;
	}
	__syncthreads();


//	int hitsNum = inHitsNum[threadId];					//real number of hits for this event


	track *resTracks = (track*)(resultsKernel);			//tracks are stored as array of structs - should be stored as arrays of ints (coalesced access)
	sensorInfo sensor0;	//not the best solution
	sensorInfo sensor1;

	int lastSensor = m_sensNum-1;
	int firstSensor = 2;
	for ( int sens0 = lastSensor; firstSensor <= sens0; sens0 -= 1 ) {
		int sens1 = sens0 - 2;
		sensor0.startPosition = inSensInfoStartPos[threadId+(sens0*eventsNumber)];
		sensor0.hitsNum	= inSensInfoHitsNum[threadId+(sens0*eventsNumber)];
		sensor0.z = inSensInfoZ[threadId+(sens0*eventsNumber)];

		sensor1.startPosition = inSensInfoStartPos[threadId+(sens1*eventsNumber)];
		sensor1.hitsNum	= inSensInfoHitsNum[threadId+(sens1*eventsNumber)];
		sensor1.z = inSensInfoZ[threadId+(sens1*eventsNumber)];

		int hit0offset;
		int hit1offset;

		double dxMax = m_maxXSlope * fabs( sensor1.z - sensor0.z );
		double dyMax = m_maxYSlope * fabs( sensor1.z - sensor0.z );

		int first1 = sensor1.startPosition;
		for(int i0=sensor0.startPosition; i0<(sensor0.startPosition + sensor0.hitsNum );++i0){
			hit0offset = (i0*eventsNumber)+threadId;
			if(isUsed[hit0offset]){
				continue;
			}
			double x0 = inX[hit0offset];  //TODO: check if it goes to register of stays in local (global) memory
			double y0 = inY[hit0offset];
			double xMin = x0 - dxMax;
			double xMax = x0 + dxMax;
			for(int i1 = first1; i1<(sensor1.startPosition + sensor1.hitsNum ); i1++){
				memset(usedHits, -1, m_sensNum*sizeof(int));	//czy to na pewno tutaj ma być?
				hit1offset = (i1*eventsNumber)+threadId;
				double x1 = inX[hit1offset];
				if(x1<xMin){
					first1 = i1+1;
					continue;
				}
				if (x1 > xMax) break;
				if (isUsed[hit1offset]) continue;

				double y1  = inY[hit1offset];
				if ( fabs( y1 - y0 ) > dyMax ) continue;

				m_track.internalId = trackId;
				m_track.firstHit = nextHitNb;

				setTrack(&m_track, hit0offset, hit1offset, tracksIds , inX, inY, inZ, inW);
				usedHits[0] = hit0offset;   //!!!
				usedHits[1] = hit1offset;	//!!!

		        if ( sensor0.z < m_maxZForRBeamCut ) {
		          double z_beam  = zBeam(&m_track);
		          if ( z_beam > sensor0.z ) {
		            double r2Beam = r2AtZ( z_beam, &m_track );
		            if ( r2Beam > m_maxR2Beam )  continue;
		          }
		        }

		        int extraStep = 2;
		        int extraSens = sens1-extraStep;
		        int nbMissed = 0;
		        int extraHitId = 2;
		        double lastZ = sensor1.z;
		        //double lastZ = inSensInfoZ[threadId + (sens1*eventsNumber)];
		        while ( extraSens >= 0 ) {
		            double tol     =  m_extraTol;
		            double maxChi2 =  m_maxChi2ToAdd;
		            if ( inSensInfoZ[threadId + (extraSens*eventsNumber)] < lastZ - 100.0 ) {
		            	tol     = 2 * tol;
		                maxChi2 = 2 * maxChi2;
		            }
		            bool added = addHitsOnSensor(inSensInfoZ[threadId+extraSens],inSensInfoStartPos[threadId+extraSens],inSensInfoHitsNum[threadId+extraSens], tol, maxChi2,
		            							&m_track, tracksIds, inX, inY, inZ,inW, threadId, eventsNumber, &usedHits[extraHitId]);
		            if ( added ) {
		            	extraHitId++;
		            	nbMissed = 0;
		                lastZ = inSensInfoZ[threadId + extraSens];
		            } else {
		            	nbMissed += extraStep;
		                extraStep = 1;
		            }
		            if ( m_maxMissed < nbMissed ) break;
		            	extraSens -= extraStep;
		         }

		        //== Try upstream if almost forward tracks
		        if ( sensor0.z > m_maxZForRBeamCut ) {
		          extraStep = 1;
		          extraSens = sens0 + 3;
		          nbMissed = 2;
		          while ( extraSens <= lastSensor ) {
		        	  //extra = &sensors[extraSens];
		        	  int sensOffset = threadId+(extraSens*eventsNumber);
		        	  bool added = addHitsOnSensor(inSensInfoZ[sensOffset],inSensInfoStartPos[sensOffset],inSensInfoHitsNum[sensOffset], m_extraTol, m_maxChi2ToAdd,
		        	  		            		   &m_track, tracksIds, inX, inY, inZ,inW, threadId, eventsNumber, &usedHits[extraHitId]);
		            if ( added ) {
		              nbMissed = 0;
		              extraHitId++;
		            } else {
		              nbMissed += extraStep;
		            }
		            if ( m_maxMissed < nbMissed ) break;
		            extraSens += extraStep;
		          }
		        }
		        removeWorstHit(&m_track,m_maxChi2PerHit,usedHits, tracksIds, inX, inY, inZ, inW, isUsed,eventsNumber,threadId );
		        if ( m_track.trackHitsNum < 3 ) continue;
		        //== Add compatible hits in sens0 and sens1.
		        int tmpOffset;
		        if(i0 != (sensor0.startPosition + sensor0.hitsNum -1)){
		        	tmpOffset = ((i0+1)*eventsNumber)+threadId ;
		        	if ( chi2Track(&m_track, inX[tmpOffset], inY[tmpOffset],inZ[tmpOffset],inW[tmpOffset]) < m_maxChi2SameSensor) {
		        		++i0;
		                addHit(&m_track, hit0offset, tracksIds, inX, inY, inZ, inW );
		            }
		        }
		        if (i1 != (sensor1.startPosition + sensor0.hitsNum -1) ) {
		        	tmpOffset = ((i1+1)*eventsNumber)+threadId ;
		        	if ( chi2Track(&m_track, inX[tmpOffset], inY[tmpOffset],inZ[tmpOffset],inW[tmpOffset]) < m_maxChi2SameSensor) {
		        		++i1;
		        		addHit(&m_track, hit1offset, tracksIds, inX, inY, inZ, inW );
		            }
		        }
		        //== Final check: if only 3 hits, all should be unused and chi2 good.
		        if ( m_track.trackHitsNum == 3 ) {
		        	if ( !all3SensorsAreDifferent(usedHits, sensorsIds) ) {
		                continue;
		            }
		        	if(nbUnused(usedHits,isUsed) != 3){
		        		continue;
		            }
		            if(chi2(&m_track, inX, inY, inZ, inW, usedHits) > m_maxChi2Short){
		            	continue;
		            }
		        } else {
		        	if ( nbUnused(usedHits,isUsed) < .6 * m_track.trackHitsNum ) {
		        		continue;
		        	}
		        }
		        resTracks[trackId] = m_track; ///--- numer ścieżki

				trackId++;
				nextHitNb += m_track.trackHitsNum;

				if ( m_track.trackHitsNum > 3 ) {
					for (int i =0 ; i<48; i++){
						int hitNumber = usedHits[i];
						if(hitNumber== -1){
							break;
					 	}
						isUsed[hitNumber] = 1;
				    }
					break;
				}
			} //i1
		} //i0
	} //sensor0
//	resStats[0] = trackId;
}

inline int _ConvertSMVer2Cores(int major, int minor)
{
    // defines for GPU Architecture types
    typedef struct
    {
        int SM; // 0xMm, M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
    };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }
    return nGpuArchCoresPerSM[7].Cores;
}

inline int getBestDevice()
{
    int current_device     = 0, sm_per_multiproc  = 0;
    int max_compute_perf   = 0, max_perf_device   = 0;
    int device_count       = 0, best_SM_arch      = 0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceCount(&device_count);

    // Find the best major SM Architecture GPU device
    while (current_device < device_count)
    {
        cudaGetDeviceProperties(&deviceProp, current_device);
        if (deviceProp.computeMode != cudaComputeModeProhibited)
        {
            if (deviceProp.major > 0 && deviceProp.major < 9999)
            {
                best_SM_arch = MAX(best_SM_arch, deviceProp.major);
            }
        }
        current_device++;
    }
    // Find the best CUDA capable GPU device
    current_device = 0;
    while (current_device < device_count)
    {
        cudaGetDeviceProperties(&deviceProp, current_device);
        // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
        if (deviceProp.computeMode != cudaComputeModeProhibited)
        {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
            {
                sm_per_multiproc = 1;
            }
            else
            {
                sm_per_multiproc = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
            }

            int compute_perf  = deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;

            if (compute_perf  > max_compute_perf)
            {
                // If we find GPU with SM major > 2, search only these
                if (best_SM_arch > 2)
                {
                    // If our device==dest_SM_arch, choose this, or else pass
                    if (deviceProp.major == best_SM_arch)
                    {
                        max_compute_perf  = compute_perf;
                        max_perf_device   = current_device;
                    }
                }
                else
                {
                    max_compute_perf  = compute_perf;
                    max_perf_device   = current_device;
                }
            }
        }
        ++current_device;
    }
    return max_perf_device;
}

void findCudaDevice()
{
    cudaDeviceProp deviceProp;
    int devID = 0;
    devID = getBestDevice();
    CUDA_CHECK_RETURN(cudaSetDevice(devID));
    CUDA_CHECK_RETURN(cudaGetDeviceProperties(&deviceProp, devID));
    printf("> Using CUDA device [%d]: %s\n", devID, deviceProp.name);
}

/*
 * Function that starts GPU Kernel.
 */
void launchKernel(vector<char> &inputVector, int blocksPerGrid, int threadsPerBlock, unsigned int allHits, int topHitsNum, char **results){
	findCudaDevice(); //choosing the fastest (GFLOPS) card on the system
	void *vecKernel = NULL;
	void *resultsKernel = NULL;

	CUDA_CHECK_RETURN(cudaDeviceReset());
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1); //TODO: we can increase L1 cache if we don't use shared mem
	int resultsInBytes = (allHits)*(sizeof(track) + 2*sizeof(int));   //Allocating memory for as many tracks as hits...we can allocate 1/3 of it but then we must know all the offeset between events
	CUDA_CHECK_RETURN(cudaMalloc((void**) &vecKernel, inputVector.size()));
	CUDA_CHECK_RETURN(cudaMemcpy(vecKernel, &inputVector[0], inputVector.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &resultsKernel, resultsInBytes));

	searchByPair<<<blocksPerGrid, threadsPerBlock>>>(vecKernel,resultsKernel,topHitsNum);
//	searchByPair<<<1,4>>>(vecKernel,resultsKernel, allHits); //testing

	CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());

	*results = (char*)malloc(resultsInBytes);  //TODO: release 'results' memory
	CUDA_CHECK_RETURN(cudaMemcpy(*results, resultsKernel, resultsInBytes, cudaMemcpyDeviceToHost));
	CUDA_CHECK_RETURN(cudaFree((void*) vecKernel));
	CUDA_CHECK_RETURN(cudaFree((void*) resultsKernel));
	CUDA_CHECK_RETURN(cudaDeviceReset());

}
