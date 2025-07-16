# This file is distributed under the MIT license.
# See the LICENSE file for details.

add_executable(hpcnrw)
target_sources(hpcnrw PRIVATE main.cpp)
target_link_libraries(hpcnrw PUBLIC visionaray visionaray_common)
