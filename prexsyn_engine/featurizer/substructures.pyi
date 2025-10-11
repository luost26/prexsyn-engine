from .base import Featurizer

class BRICSFragmentsFeaturizer(Featurizer):
    def __init__(self, name: str, fp_type: str, max_num_fragments: int = 8) -> None: ...
