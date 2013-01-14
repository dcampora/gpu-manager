################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/gpu-manager.cu 

CU_DEPS += \
./src/gpu-manager.d 

OBJS += \
./src/gpu-manager.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -O3 -gencode arch=compute_10,code=sm_10 -gencode arch=compute_30,code=sm_30 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	nvcc --compile -O3 -gencode arch=compute_10,code=compute_10 -gencode arch=compute_10,code=sm_10 -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


