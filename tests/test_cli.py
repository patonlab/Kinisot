#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""End-to-end tests of the command-line interface (python -m kinisot)."""

import os
import subprocess
import sys

import pytest
from conftest import datapath

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))


def run_cli(args, cwd):
    env = dict(os.environ)
    env['PYTHONPATH'] = REPO_ROOT + os.pathsep + env.get('PYTHONPATH', '')
    return subprocess.run([sys.executable, '-m', 'kinisot'] + args,
                          capture_output=True, text=True, cwd=str(cwd), env=env)


def test_cli_kie_gaussian(tmp_path):
    result = run_cli(['--rct', datapath('gaussian/claisen_gs.out'),
                      '--ts', datapath('gaussian/claisen_ts.out'),
                      '--iso', '5', '-t', '393', '-s', '0.961'], tmp_path)
    assert result.returncode == 0, result.stderr
    dat = (tmp_path / 'Kinisot_output.dat').read_text()
    # corr-KIE golden value for this case, printed at 6 decimals
    assert '1.001940' in dat


def test_cli_eqe_orca(tmp_path):
    result = run_cli(['--rct', datapath('orca/pentane_TT.out'),
                      '--prd', datapath('orca/pentane_GG.out'),
                      '--iso', '6', '-s', '1.0'], tmp_path)
    assert result.returncode == 0, result.stderr
    assert '1.007728' in (tmp_path / 'Kinisot_output.dat').read_text()


def test_cli_version(tmp_path):
    result = run_cli(['--version'], tmp_path)
    assert result.returncode == 0
    assert 'Kinisot v' in result.stdout


@pytest.mark.parametrize("args, fragment", [
    # neither --ts nor --prd
    (['--rct', 'x.out', '--iso', '1'], 'requires either a TS'),
    # both --ts and --prd
    (['--rct', 'x.out', '--ts', 'y.out', '--prd', 'z.out', '--iso', '1'], 'not both'),
    # 3 files but only 2 label sets
    (['--rct', 'a.out', '--rct', 'b.out', '--ts', 'c.out', '--iso', '1', '--iso', '2'], 'labels for each'),
    # unknown flag is rejected (parse_args, not parse_known_args)
    (['--rct', 'x.out', '--ts', 'y.out', '--iso', '1', '--bogus'], 'unrecognized'),
], ids=['no-ts-no-prd', 'ts-and-prd', 'label-count', 'unknown-flag'])
def test_cli_argument_errors(tmp_path, args, fragment):
    result = run_cli(args, tmp_path)
    assert result.returncode == 2  # argparse error exit code
    assert fragment in result.stderr
    # bad arguments must not leave an output file behind
    assert not (tmp_path / 'Kinisot_output.dat').exists()


def test_cli_no_imaginary_ts_exits_cleanly(tmp_path):
    # A ground state passed as --ts: clean error message, exit 1, log closed
    result = run_cli(['--rct', datapath('gaussian/dienophile.out'),
                      '--ts', datapath('gaussian/diene.out'),
                      '--iso', '0', '-s', '1.0'], tmp_path)
    assert result.returncode == 1
    assert 'imaginary frequency' in result.stdout
