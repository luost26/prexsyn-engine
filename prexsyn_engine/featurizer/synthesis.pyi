from .base import Featurizer

class PostfixNotationTokenDef:
    PAD: int
    END: int
    START: int
    BB: int
    RXN: int
    num_token_types: int

    def __init__(self, pad: int = ..., end: int = ..., start: int = ..., bb: int = ..., rxn: int = ...) -> None: ...

class PostfixNotationFeaturizer(Featurizer):
    def __init__(self, max_length: int = 16, token_def: PostfixNotationTokenDef = ...) -> None: ...
