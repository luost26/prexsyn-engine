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
    from . import types as types
from . import building_block_list as building_block_list
from . import reaction_list as reaction_list
from . import synthesis as synthesis
from . import indexer as indexer
from . import chemspace as chemspace
from . import featurizer as featurizer
from . import detokenizer as detokenizer

__version__ = "0.1.0"
