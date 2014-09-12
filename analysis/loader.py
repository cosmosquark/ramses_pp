from ramses_pp import config

pymses_loaded = config.pymses_enabled
pynbody_loaded = config.pynbody_enabled
yt_loaded = config.yt_enabled

if config.quick_import:
	if config.pymses_enabled:
		from ..modules.pymses import Pymses
	if config.pynbody_enabled:
		from ..modules.pynbody import Pynbody
	if config.yt_enabled:
		from ..modules.yt import YT
else:

	config.list_modules()

	if pymses_loaded:
		try:
			from .pymses import Pymses
		except ImportError as e:
			print 'Unable to import pymses'
			pymses_loaded = False
			print e
	if pynbody_loaded:
		try:
			from .pynbody import Pynbody
		except ImportError as e:
			print 'Unable to import pynbody'
			pynbody_loaded = False
			print e
	if yt_loaded:
		try:
			from .yt import YT
		except ImportError as e:
			print 'Unable to import yt'
			yt_loaded = False
			print e
		if (pymses_loaded or pynbody_loaded or yt_loaded) is False:
			raise RuntimeError("Could not import any modules!")
			
import numpy as np
import json, os, glob, uuid
