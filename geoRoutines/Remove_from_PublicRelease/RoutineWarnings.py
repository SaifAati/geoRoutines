class NotGeoreferencedError(Exception):
    """
    Custom rpcm Exception.
    """
    pass


class NoDEMWarning(Warning):
    """
    Custom rpcm warning raised when no SRTM altitude is available.
    """
    pass

class MaxLocalizationIterationsError(Exception):
    """
    Custom rpcm Exception.
    """
    pass