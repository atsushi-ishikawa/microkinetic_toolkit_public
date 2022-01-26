#!/usr/bin/env python3


import os
from subprocess import run


NEEDED_FILES = ["test.py", "example.py", "test.csv", "ase.json",
                "molecules.json"]


def run_cmd(cmd):
    _ = run(cmd, shell=True,
            envs=os.environ.copy(), cwd=os.getcwd(),
            capture_output=True)
    print(_.stdout)
    print(_.stderr)
    assert _.check_returncode == 0


def _initialize():
    cwd = os.getcwd()
    for fnm in NEEDED_FILES:
        fp = os.path.join(cwd, fnm)
        assert os.path.exists(fp)


def _delte_files_except():
    abs_fps = []
    cwd = os.getcwd()
    for fnm in NEEDED_FILES:
        fp = os.path.join(cwd, fnm)
        abs_fps.append(fp)
    for candidate in os.listdir(cwd):
        cand_fp = os.path.join(cwd, candidate)
        if cand_fp in abs_fps:
            continue
        else:
            os.remove(cand_fp)


def test_run_all():
    self._initialize()
    self._delte_files_except()
    cmd = "python sample.py"
    run_cmd(cmd)


def test_run_from_pickle():
    self._initialize()
    cmd = "python sample.py"
    run_cmd(cmd)
