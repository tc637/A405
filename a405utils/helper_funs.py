def test_scalar(*args):
    """
      return true if every argument is a scalar
    """
    isscalar=True
    for item in args:
        isscalar = isscalar & np.isscalar(item)
    return isscalar
