set -x

case $(uname | tr '[:upper:]' '[:lower:]') in
  linux*)
		# Nothing to do
    ;;
  darwin*)
		export CFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=14.1 ${CFLAGS}"
		export CXXFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=14.1 ${CXXFLAGS}"
    ;;
  *)
esac

${PYTHON} -m pip install . --no-deps --ignore-installed -vv