#ifndef MAIN_H_
#define MAIN_H_

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/named_condition.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <algorithm>
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <time.h>
#include "tools.h"
#include "SearchKernel.cuh"
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <algorithm>
using namespace boost::interprocess;
struct shared_conditions {
	shared_conditions()
        :  package_in(false)
        {}
    boost::interprocess::interprocess_mutex      mutex;
    boost::interprocess::interprocess_condition  package_full;
    boost::interprocess::interprocess_condition  package_empty;
    bool package_in;
};

#endif /* MAIN_H_ */
