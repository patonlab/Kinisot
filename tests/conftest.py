#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

# Example data lives at the repository root (not inside the installed
# package), so tests run against the checkout regardless of how the
# package itself was installed.
HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.normpath(os.path.join(HERE, '..', 'examples'))


def datapath(path):
    return os.path.join(EXAMPLES, path)
