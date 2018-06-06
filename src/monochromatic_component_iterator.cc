//
// Created by Fatemeh Almodaresi on 6/4/18.
//

#include "monochromatic_component_iterator.h"

uint64_t start_time;

monochromatic_component_iterator::monochromatic_component_iterator(const CQF<KeyObject> *g)
        : cqf(g), it(g->begin(0)) {
    // initialize cqf iterator
    k = cqf->keybits()/2;
    std::cerr << " k : " << k << "\n";
    //it = &(cqf->begin(0));
    sdsl::util::assign(visited, sdsl::bit_vector(cqf->slots(), 0));
    std::cerr << "slots: " << cqf->slots() << "\n";
    // initialize the first item and get to done if it's empty
    ++(*this);

}

monochromatic_component_iterator::work_item
monochromatic_component_iterator::front(std::queue <monochromatic_component_iterator::work_item> &w) {
    return w.front();
}

bool monochromatic_component_iterator::done() {return it.done();}
void monochromatic_component_iterator::operator++(void) {
    std::cerr << "++ ";
    while (!it.done()) {
        std::cerr << " , " << it.iter.current << " " << (bool)(visited[it.iter.current]) << " ";
        if ((bool)(visited[it.iter.current])) {
            ++it;
            //std::cerr << " ++it ";
        }
        else {
            break;
            std::cerr << " break\n";
        }
    }
    std::cerr << "\n";
    if (it.done()) return;

    KeyObject keyobj = *it;
    node root(k, HashUtil::hash_64i(keyobj.key, BITMASK(cqf->keybits())));
    monochromatic_component_iterator::work_item neww = {root,0, 0};
    work.push(neww);
}

Mc_stats monochromatic_component_iterator::operator*(void) {
    if (work.empty()) {
        std::cerr << "Throw Exception. We're out of bound for CQF.\n";
        std::exit(1);
    }
    //return getMC();
    Mc_stats res;
    while (!work.empty()) {
        work_item w = front(work);
        work.pop();

        std::cerr << "for w " << w.idx << " : ";
        for (auto &neighbor : neighbors(w)) {
            std::cerr << "n " << neighbor.idx << " ";
            if (neighbor.curr != w.curr) {
                /*if (neighbor.colorid != w.colorid) {
                    res.min_dist = std::min(res.min_dist, manhattanDist(neighbor.colorid, w.colorid));
                }*/
                if (visited[neighbor.idx] == 0) {
                    if (neighbor.colorid == w.colorid) {
                        res.nodeCnt++;
                        work.push(neighbor);
                    }
                }
            }
            else {
                std::cerr << " same || ";
            }
        }
        visited[w.idx] = 1; // set the corresponding bit
        std::cerr << " idx " << w.idx << " visited " << visited[w.idx] << "\n";

    }
    return res;
}

std::set <monochromatic_component_iterator::work_item>
        monochromatic_component_iterator::neighbors(monochromatic_component_iterator::work_item n)  {
    std::set <work_item> result;
    for (const auto b : dna::bases) {
        uint64_t eqid, idx;
        if (exists(b + n.curr, idx, eqid))
            result.insert(work_item(b >> n.curr, idx, eqid));
        if (exists(n.curr + b, idx, eqid))
            result.insert(work_item(n.curr << b, idx, eqid));
    }
    return result;
}

bool monochromatic_component_iterator::exists(edge e, uint64_t& idx, uint64_t& eqid) {
    // no correction for exact CQF. right?
    uint64_t tmp = e.val;
    KeyObject key(HashUtil::hash_64(tmp, BITMASK(cqf->keybits())), 0, 0);
    auto eq_idx = cqf->queryValAndIdx(key);
    if (eq_idx.first) {
        eqid = eq_idx.first;
        idx = eq_idx.second;
        return true;
    }
    return false;
}


/*
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:
 * ===========================================================================
 */

int main ( int argc, char *argv[] ) {

    std::string cqf_file = argv[1];
    std::cout << cqf_file << "\n";
    std::string eq_file = argv[2];
    std::cout << eq_file << "\n";
    CQF<KeyObject> cqf(cqf_file, false);
    monochromatic_component_iterator mci(&cqf);
    //while(!mci.done()) {
    for (auto i = 0; i < 10; i++) {
        std::cout << (*mci).nodeCnt << "\n";
        ++mci;
    }
}
