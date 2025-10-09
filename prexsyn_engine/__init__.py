# isort: skip_file
import warnings

# Importing RDKit modules to trigger some initialization
# otherwise, random segfaults can occur
import rdkit.Chem.rdchem
import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdMolChemicalFeatures
import rdkit.Chem.Descriptors

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="^.*already registered.*$", category=RuntimeWarning)
    from . import types
from . import building_block_list
from . import reaction_list
from . import synthesis
from . import indexer
from . import chemspace
from . import featurizer
from . import detokenizer

__version__ = "0.1.0"
