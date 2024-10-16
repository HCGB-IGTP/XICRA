from pkg_resources import get_distribution
try:
	__version__ = get_distribution('XICRA').version
except:
	try:
		__version__ = pkg_resources.resorce_filename('XICRA', VERSION)
	except:
		__version__ = 'local'

__all__ = [
	'modules',
	'scripts',
	'config',
	'other_tools'
]

from XICRA import *

