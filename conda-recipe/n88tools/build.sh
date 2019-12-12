set -x

case $(uname | tr '[:upper:]' '[:lower:]') in
  linux*)
		# Nothing to do
    ;;
  darwin*)
		export CFLAGS="-isysroot /opt/MacOSX10.9.sdk -mmacosx-version-min=10.9 ${CFLAGS}"
		export CXXFLAGS="-isysroot /opt/MacOSX10.9.sdk -mmacosx-version-min=10.9 ${CXXFLAGS}"
    ;;
  *)
esac

${PYTHON} -m pip install . --no-deps --ignore-installed -vv

