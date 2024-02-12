#pragma once 

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <map>
#include <algorithm>
#include <utility>
#include <cassert>

class LightAMR {

    private:
        
        // amr data
        int _nlevels;
        int _levelminDesc; // levelmin of lightAMR desc !! ( can be != simulation levelmin  (>=) )
        int _ndomains;
        int _current_domain;

        int _ndim;
        int _refinement_factor;
        int _nleaf_th;
        uint64_t _ncells;

        // if AMR was not compressed 
        std::vector<uint8_t> _amr_tree;
        std::vector<uint8_t> _mask;

        // contains number of cells by level
        std::vector<uint64_t> _ncells_per_level;

        // contains logical coordinate of coarse cells (level 0)
        std::vector<uint32_t> _coarse_ijk;

    public:

        // constructor
        LightAMR( int d ){
            _ndim               = d;
            _refinement_factor  = 0;
            _nleaf_th           = 0;
            _ncells             = 0;
            _current_domain     = 0;
        }

        // forbid copy construct C++11
        LightAMR(const LightAMR &) = delete;
        LightAMR(LightAMR &&) = default;

        // assignement operator
        LightAMR& operator=(const LightAMR & lamr){ return *this; }

        // destructor
        ~LightAMR() = default;

        // Simple setter ----------------------------------------------------------
        void setDomainID(int id)             { _current_domain = id; }
        void setNumberOfLevel(int)           { assert(false); }
        void setTotalNumberOfCell(int64_t)   { assert(false); }
        void setRefinementFactor( int8_t r ) { _refinement_factor = r; }
        void setDim( int d )                 { _ndim = d; }

        std::tuple<bool, std::string>  isValide() const { assert(false); return std::make_tuple( false, ""); }

        // Getter reference on internal vectors
        std::vector<uint8_t>        & getReferenceAMRtree() { return _amr_tree; }
        std::vector<uint8_t>        & getReferenceMask()    { return _mask; }

        const std::vector<uint8_t>  & getConstReferenceAMRtree() const { return _amr_tree; }
        const std::vector<uint8_t>  & getConstReferenceMask() const { return _mask; }

        std::vector<uint64_t>       & getReferenceNCellsPerLevel() { return _ncells_per_level; }
        const std::vector<uint64_t> & getConstReferenceNCellsPerLevel() const { return _ncells_per_level; }

        std::vector<uint32_t>       & getReferenceCoarseIJK(){ return _coarse_ijk; }
        const std::vector<uint32_t> & getConstReferenceCoarseIJK() const { return _coarse_ijk; }

        // Mesh info
        int    getNumberOfLevels()    const { return static_cast<int>( _ncells_per_level.size() ); }
        int    getDim()               const { return _ndim;  }
        int    getDomainID()          const { return _current_domain; }

        // this include nodes & masked cells
        size_t getNumberOfCells()     const { return _amr_tree.size(); }
        int    getRefinementFactor()  const { return _refinement_factor; }

        // Memory info for LightAMR data structure
        double getMemoryConsumption() const { assert(false); return 0.; }
        double getMemoryAMRDesc()     const { assert(false); return 0.; }
        void   printMemoryInfo()      const { assert(false); }

        // conversion to other data model
        // void export_description_VTK(const std::string&);

        /* -------------- Debug and info -------------- */ 
        void print() const { assert(false); }

};