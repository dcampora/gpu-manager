################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/main.cpp \
../src/tools.cpp 

CU_SRCS += \
../src/SearchKernel.cu 

CU_DEPS += \
./src/SearchKernel.d 

OBJS += \
./src/SearchKernel.o \
./src/main.o \
./src/tools.o 

CPP_DEPS += \
./src/main.d \
./src/tools.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -I-lrt -G -g -O0 -gencode arch=compute_30,code=sm_30 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	nvcc --compile -G -I-lrt -O0 -g -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -I-lrt -G -g -O0 -gencode arch=compute_30,code=sm_30 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	nvcc -I-lrt -G -g -O0 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


