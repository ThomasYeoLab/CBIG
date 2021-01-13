from .embed import has_sklearn

if has_sklearn:
    from .embed import DiffusionMapEmbedding