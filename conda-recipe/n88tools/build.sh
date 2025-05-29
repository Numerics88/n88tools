set -x

case $(uname | tr '[:upper:]' '[:lower:]') in
  linux*)
		# Nothing to do
    ;;
  darwin*)
		export CFLAGS="-isysroot /opt/MacOSX11.3.sdk -mmacosx-version-min=11.3 ${CFLAGS}"
		export CXXFLAGS="-isysroot /opt/MacOSX11.3.sdk -mmacosx-version-min=11.3 ${CXXFLAGS}"
    ;;
  *)
esac

${PYTHON} -m pip install . --no-deps --ignore-installed -vv