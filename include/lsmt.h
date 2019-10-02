/*
 * ============================================================================
 *
 *         Author:  Jamshed Khan, jamshed@cs.umd.edu
 *   Organization:  University of Maryland
 *
 * ============================================================================
 */

#ifndef _LSMT_H_
#define _LSMT_H_

#include "coloreddbg.h"
#include "BooPHF.h"



template <typename qf_obj, typename key_obj>
class LSMT
{
    public:
        std::string dir;
        // std::string pendingSamplesList;
        uint scalingFactor;
        uint64_t kmerThreshold;
        uint64_t sampleThreshold;



        LSMT(nlohmann::json &params);
};



template<typename qf_obj, typename key_obj>
LSMT<qf_obj, key_obj>::
    LSMT(nlohmann::json &params)
{
}

#endif