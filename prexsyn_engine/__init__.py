# isort: skip_file
import warnings

# Importing RDKit modules to trigger initialization
# otherwise, segfaults can occur
import rdkit.Chem.rdchem  # pyright: ignore [reportUnusedImport]
import rdkit.Chem.rdChemReactions  # pyright: ignore [reportUnusedImport]
import rdkit.Chem.rdMolChemicalFeatures  # pyright: ignore [reportUnusedImport]
import rdkit.Chem.Descriptors  # pyright: ignore [reportUnusedImport]

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

__version__ = "0.1.1.dev4"
