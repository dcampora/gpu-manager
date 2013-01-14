
#include "gpu-manager.h"

using namespace boost::interprocess;

int main(){
	typedef allocator<char, managed_shared_memory::segment_manager> ShmemAllocator;
	typedef vector<char, ShmemAllocator> Gaudi_vector;
	Gaudi_vector *myvector;
	shared_memory_object::remove("gpu_memory");

	//Create a shared memory object.
	shared_memory_object shm (create_only,"gpu_memory" ,read_write );
	try{
	  shm.truncate(sizeof(shared_conditions));
	  mapped_region region
		 (shm
		 ,read_write
		 );
	  void * addr       = region.get_address();
	  shared_conditions * cond = new (addr) shared_conditions;
	  while(true)
	  {
		  scoped_lock<interprocess_mutex> lock(cond->mutex);
		  if(!cond->package_in){
			  cond->package_empty.wait(lock);
		  }
		  printf("Received package from Gaudi \n");
		  managed_shared_memory segment (open_only ,"gaudi_memory");  //segment name
		  myvector = segment.find<Gaudi_vector>("gaudi_vector").first;
		  cond->package_in = false;
		  //START GPU ALGORITHM!
		  cond->package_full.notify_one();
	  }
	}
	catch(...){
	  shared_memory_object::remove("gpu_memory");
	  throw;
	}
	shared_memory_object::remove("gpu_memory");
	return 0;
}
