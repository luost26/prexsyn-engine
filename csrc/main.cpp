#include <iostream>
#include <memory>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

int main() {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    std::cout << "Molecule: " << mol->getNumAtoms() << " atoms\n";
    return 0;
}
