include(FetchContent)
FetchContent_Declare(
  tinyad
  GIT_REPOSITORY git@github.com:the13fools/gpgpt.git
  GIT_TAG 770dc5011371787ef784019f62794ec560623f79
)
FetchContent_GetProperties(tinyad)
if(NOT tinyad_POPULATED)
  # Fetch the content using previously declared details
  FetchContent_Populate(tinyad)
  message(STATUS "tinyad_SOURCE_DIR: ${tinyad_SOURCE_DIR}")
  message(STATUS "tinyad_BINARY_DIR: ${tinyad_BINARY_DIR}")
  add_subdirectory(${tinyad_SOURCE_DIR} ${tinyad_BINARY_DIR})
endif()
FetchContent_MakeAvailable(tinyad)