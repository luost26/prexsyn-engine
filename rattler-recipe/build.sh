mkdir builddir

$PYTHON -m build -w -n -x -Cbuilddir=builddir

$PYTHON -m pip install dist/*.whl
