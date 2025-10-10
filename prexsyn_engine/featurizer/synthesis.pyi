from .base import Featurizer

class PostfixNotationTokenDef:
    PAD: int
    END: int
    START: int
    BB: int
    RXN: int

class PostfixNotationFeaturizer(Featurizer):
    def __init__(self, max_length: int = 16, token_def: PostfixNotationTokenDef = ...) -> None: ...
