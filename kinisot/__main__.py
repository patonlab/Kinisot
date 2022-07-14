# Copied from __main__.py in pip
from __future__ import absolute_import

import os
import sys

if __package__ == '':
    path = os.path.dirname(os.path.dirname(__file__))
    sys.path.insert(0, path)

from kinisot import Kinisot

if __name__ == '__main__':
    sys.exit(Kinisot.main())
