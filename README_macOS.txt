Reticula 2.11.0 for macOS
=========================
Publication figures for framework materials
Built by Stijn Germonpré

This build is for Apple Silicon Macs (M1, M2, M3, M4 — arm64).
macOS 13 Ventura or later is recommended.

Install
-------
1. Open the disk image.
2. Drag Reticula.app to the Applications folder (the shortcut is in the window).
3. Eject the disk image.

First launch (the app is unsigned)
---------------------------------
This release is not signed with an Apple Developer ID, so macOS Gatekeeper
will warn that it "cannot be opened because the developer cannot be verified".
That is expected. Choose one of these:

  • In Finder, Control-click (or right-click) Reticula.app → Open → Open.

  • Or, in Terminal, remove the quarantine flag and then open the app:

      xattr -dr com.apple.quarantine /Applications/Reticula.app
      open /Applications/Reticula.app

After you have opened it once this way, later launches work as usual.

What happens when it starts
---------------------------
A small control window opens and your default browser shows a loading screen,
then the app. Keep the control window open while you work; closing it stops
the server. "Open in browser" brings the app back if you closed the tab.

Project files
-------------
Saving a project writes a .reticula file. After the first time you open one
with Reticula (right-click → Open With), you can double-click .reticula files
to reopen them.

Where your work is stored
-------------------------
  Settings, page state, sample store, logs:
    ~/Library/Application Support/Reticula/

  Figures saved with "Save to folder":
    ~/Library/Application Support/Reticula/exports/

  The log file (useful if something goes wrong) is server_log.txt in that folder.

Requirements and limits
-----------------------
  • This build is unsigned (see above).
  • It was produced on a GitHub-hosted macOS runner and has not been tested
    on a personal Mac. If something fails to start, check the log file and
    get in touch.
  • An Intel (x86_64) build may be added later; this file is the arm64 one.
  • The Windows installer is on the same GitHub release page.

Download
--------
https://github.com/Spru1tje/reticula-releases/releases/tag/v2.11.0
