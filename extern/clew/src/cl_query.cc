/*
 * Copyright 2011-2015 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cstdlib>
#include <cstdio>
#include <vector>

#include "clew.h"

using std::vector;
#define foreach(x, y) for(x : y)

namespace {

// Get all platform IDs.
bool getPlatformIDs(vector<cl_platform_id>* platform_ids) {
  cl_uint num_platforms = 0;
  if (clGetPlatformIDs(0, NULL, &num_platforms) != CL_SUCCESS) {
    return false;
  }
  platform_ids->resize(num_platforms);
  if (clGetPlatformIDs(
      num_platforms, &(*platform_ids)[0], NULL) != CL_SUCCESS) {
    return false;
  }
  return true;
}

// Get device IDs for the given platform.
bool getPlatformDeviceIDs(cl_platform_id platform_id,
                          vector<cl_device_id>* device_ids) {
  cl_uint num_devices = 0;
  if (clGetDeviceIDs(
      platform_id, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices) != CL_SUCCESS) {
    return false;
  }
  device_ids->resize(num_devices);
  if (num_devices == 0) {
    return true;
  }
  if (clGetDeviceIDs(platform_id,
                     CL_DEVICE_TYPE_ALL,
                     num_devices,
                     &(*device_ids)[0],
                     NULL) != CL_SUCCESS) {
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  int result = clewInit();
  if (result != CLEW_SUCCESS) {
    return EXIT_FAILURE;
  }
  // Get all platform IDs.
  vector<cl_platform_id> platform_ids;
  if (!getPlatformIDs(&platform_ids)) {
    return EXIT_FAILURE;
  }
  foreach (cl_platform_id platform_id, platform_ids) {
    vector<cl_device_id> device_ids;
    if (!getPlatformDeviceIDs(platform_id, &device_ids)) {
      continue;
    }
    foreach (cl_device_id device_id, device_ids) {
      char name[1024] = "\0";
      if (clGetDeviceInfo(
          device_id, CL_DEVICE_NAME, sizeof(name), name, NULL) != CL_SUCCESS) {
        continue;
      }
      cl_int max_compute_units = 0;
      if (clGetDeviceInfo(device_id,
                          CL_DEVICE_MAX_COMPUTE_UNITS,
                          sizeof(max_compute_units),
                          &max_compute_units,
                          NULL) != CL_SUCCESS) {
        continue;
      }
      printf("%s:%d\n", name, max_compute_units);
    }
  }
  return EXIT_SUCCESS;
}
