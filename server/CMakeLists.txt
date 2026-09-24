include(FetchContent)

FetchContent_Declare(
    httplib
    GIT_REPOSITORY https://github.com/yhirose/cpp-httplib.git
    GIT_TAG v0.18.3
)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(httplib nlohmann_json)

add_executable(riskengine_server server_main.cpp)
target_link_libraries(riskengine_server PRIVATE riskengine httplib nlohmann_json::nlohmann_json)
target_compile_definitions(riskengine_server PRIVATE
    RISKENGINE_STATIC_DIR="${CMAKE_CURRENT_SOURCE_DIR}/static")
