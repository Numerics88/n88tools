set -x

case $(uname | tr '[:upper:]' '[:lower:]') in
  linux*)
		# Nothing to do
    ;;
  darwin*)
		export CFLAGS="-isysroot /Applications/Xcode_14.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -mmacosx-version-min=14.1 ${CFLAGS}"
		export CXXFLAGS="-isysroot /Applications/Xcode_14.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -mmacosx-version-min=14.1 ${CXXFLAGS}"
    ;;
  *)
esac

${PYTHON} -m pip install . --no-deps --ignore-installed -vv