#!/bin/bash
# Double-click this file (Finder) to build Reticula for macOS on this Mac. Instructions: BUILD_ON_MAC.txt
cd "$(dirname "$0")"
xattr -dr com.apple.quarantine . 2>/dev/null || true
chmod +x pkg/*.sh build_mac.command 2>/dev/null || true
bash pkg/build_mac.sh
status=$?
echo
if [[ $status -eq 0 ]]; then
  echo "Finished. You can close this window."
else
  echo "The build stopped with an error (code $status). Scroll up, or open build_mac.log and send it to the developer."
fi
echo
read -r -p "Press Return to close this window. "
exit $status
