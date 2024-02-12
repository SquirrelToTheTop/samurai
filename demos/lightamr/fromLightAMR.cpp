// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <CLI/CLI.hpp>

// C++ include 
#include <iostream>
#include <filesystem>
#include <vector>
#include <cassert>
#include <cstdint>

// LightAMR
#include "./lightAMR_walks.h"

// Samurai includes 
#include <samurai/cell_array.hpp>
#include <samurai/cell_list.hpp>
#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/samurai.hpp>

namespace fs = std::filesystem;
constexpr int dim = 2;

// Build 2D Mesh as in Samurai's doc
samurai::CellList<dim> getMesh_demo() {
    
    samurai::CellList<dim> cl;

    cl[0][{0}].add_interval({0, 4});
    cl[0][{1}].add_interval({0, 1});
    cl[0][{1}].add_interval({3, 4});
    cl[0][{2}].add_interval({0, 1});
    cl[0][{2}].add_interval({3, 4});
    cl[0][{3}].add_interval({0, 3});

    cl[1][{2}].add_interval({2, 6});
    cl[1][{3}].add_interval({2, 6});
    cl[1][{4}].add_interval({2, 4});
    cl[1][{4}].add_interval({5, 6});
    cl[1][{5}].add_interval({2, 6});
    cl[1][{6}].add_interval({6, 8});
    cl[1][{7}].add_interval({6, 7});

    cl[2][{8}].add_interval({8, 10});
    cl[2][{9}].add_interval({8, 10});
    cl[2][{14}].add_interval({14, 16});
    cl[2][{15}].add_interval({14, 16});

    return cl;

}

samurai::CellList<dim> getMesh_lightamr( const LightAMR & lamr ) {
    
    samurai::CellList<dim> cl;

    auto cells = LightAMR_Utils::WA_uc_ijk_based( lamr, lamr.getNumberOfLevels() );
    std::cout << "\t>[getMesh_lightamr]:: Number of cells : " << cells.size() << std::endl;

    // fix offset ?
    int baseLevel = 2;

    for( size_t ic=0; ic<cells.size(); ++ic ){
        const auto & cell = cells[ ic ];
        cl[ cell.level - baseLevel ][ { static_cast<int>( cell.coord.j ) } ].add_point( static_cast<int>( cell.coord.i ) );
    }

    return cl;

}

int main(int argc, char* argv[])
{
    samurai::initialize(argc, argv);

    // Output parameters
    fs::path path        = fs::current_path();

    std::string fn_addInter = "2d_mesh_addInterval";
    std::string fn_fromLigh = "2d_mesh_fromLightamr";

    CLI::App app{"Create mesh from CellList and save it"};
    app.add_option("--path", path, "Output path")->capture_default_str()->group("Ouput");
    app.add_option("--filenameSamurai", fn_addInter, "File name prefix")->capture_default_str()->group("Ouput");
    app.add_option("--filenameLightAMR", fn_fromLigh, "File name prefix")->capture_default_str()->group("Ouput");
    CLI11_PARSE(app, argc, argv);

    if (!fs::exists(path)){ fs::create_directory(path); }

    { // build AMR mesh based on intervall add function

        // build cell list
        auto cl = getMesh_demo();

        // build cell array from cell list
        const samurai::CellArray<dim> ca{cl};

        // add field 'level'
        auto lf = samurai::make_field<std::size_t, 1>("level", ca);

        // "compute" level from cells
        samurai::for_each_cell(ca, [&](auto& cell){
            lf[cell] = cell.level;
        });

        std::cout << ca << std::endl;
        samurai::save(path, fn_addInter, ca, lf);

    }

    { // build AMR mesh based on lightAMR (cell logical coordinates)

        LightAMR lamr ( 2 ); // dimension
        auto & tree = lamr.getReferenceAMRtree();
        auto & mask = lamr.getReferenceMask();
        auto & cijk = lamr.getReferenceCoarseIJK();
        auto & ncpl = lamr.getReferenceNCellsPerLevel();
        
        tree = { 1,                                            // level-0
                 1, 1, 1, 1,                                   // level-1
                 0,0,0,1, 0,0,1,0, 0,1,0,0, 1,0,0,1,           // level-2
                 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,1,  // level-3
                 0,0,0,0, 0,0,0,0                              // level-4
               };
        mask.resize( tree.size(), 0 );
        cijk = { 0, 0 };
        ncpl = { 1, 4, 16, 20, 8 };        

        // build cell list
        auto cl = getMesh_lightamr( lamr );

        // build cell array from cell list
        const samurai::CellArray<dim> ca{cl};

        // add field 'level'
        auto lf = samurai::make_field<std::size_t, 1>("level", ca);

        // "compute" level from cells
        samurai::for_each_cell(ca, [&](auto& cell){
            lf[cell] = cell.level;
        });

        std::cout << ca << std::endl;
        samurai::save(path, fn_fromLigh, ca, lf);

    }

    samurai::finalize();
    return 0;
}
