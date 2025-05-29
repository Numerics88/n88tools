set -x

case $(uname | tr '[:upper:]' '[:lower:]') in
  linux*)
		# Nothing to do
    ;;
  darwin*)
		export CFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} ${CFLAGS}"
		export CXXFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} ${CXXFLAGS}"
    ;;
  *)
esac

${PYTHON} -m pip install . --no-deps --ignore-installed -vv