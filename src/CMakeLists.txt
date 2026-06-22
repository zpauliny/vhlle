# Library: all sources except main.cpp
set(LIB_SRC
    cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoCMF.cpp eoCMFe.cpp
    eoHadron.cpp eoAZH.cpp eoSmash.cpp trancoeff.cpp fld.cpp hdo.cpp
    s95p.cpp icurqmd.cpp ic.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp
    icDynFlu.cpp rmn.cpp cornelius.cpp icGlauber.cpp icGubser.cpp
    icGlissando.cpp icTrento.cpp icTrento3d.cpp icSuperMC.cpp vtk.cpp
    icTest.cpp particle.cpp simulation.cpp
)

# Library and namespaced alias
add_library(vhlle STATIC ${LIB_SRC})
add_library(vhlle::vhlle ALIAS vhlle)

# Include dirs still needed for consumers
target_include_directories(vhlle
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>
)

# Compile options
target_compile_options(vhlle PRIVATE -Wall -fPIC -O3)

target_link_libraries(vhlle PUBLIC GSL::gsl GSL::gslcblas)

add_executable(hlle_visc main.cpp)
target_link_libraries(hlle_visc PRIVATE vhlle::vhlle)
set_target_properties(hlle_visc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Install with standard defaults
install(TARGETS vhlle hlle_visc
    EXPORT vhlle-targets
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include)

# Parameter tables necessary for runtime are included 
if(VHLLE_PARAMS_DIR)
    add_custom_command(TARGET hlle_visc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${VHLLE_PARAMS_DIR}/eos $<TARGET_FILE_DIR:hlle_visc>/eos
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${VHLLE_PARAMS_DIR}/ic $<TARGET_FILE_DIR:hlle_visc>/ic
        COMMENT "Creating symlinks to parameter tables in ${VHLLE_PARAMS_DIR}"
    )

    install(DIRECTORY ${VHLLE_PARAMS_DIR}/eos DESTINATION share/vhlle)
    install(DIRECTORY ${VHLLE_PARAMS_DIR}/ic  DESTINATION share/vhlle)
endif()