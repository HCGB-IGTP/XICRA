#from pkg_resources import get_distribution
#try:
#    __version__ = get_distribution('XICRA').version
#except:
#    __version__ = 'local'
### to include when distribution available
## version will be retrieve from setup.py


__all__ = [
	'modules',
	'scripts',
	'config'
]

from XICRA import *

