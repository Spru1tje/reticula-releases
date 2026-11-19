Reticula for Windows (64-bit Windows 10/11)
From raw data to figures, the easy way.
============================================

INSTALL
1. Download Reticula_Setup_<version>.exe from the Releases page and double-click it.
   No administrator rights are needed; it installs for your user only
   (in %LOCALAPPDATA%\Programs\Reticula) and adds a Start-menu entry and, if you tick it,
   a desktop shortcut. You can also tick "Open .reticula project files with Reticula".
2. Windows SmartScreen will probably say "Windows protected your PC", because the program is
   not code-signed. Click "More info" -> "Run anyway". The same may happen the first time you
   start the app itself.

START
- Double-click "Reticula" (Start menu or desktop).
- A small window appears and your default web browser opens at http://127.0.0.1:8501
  (8502, 8503 ... if that port is busy).
- The first start can take 30-60 s (Windows Defender scans the files once). Later starts take a few seconds.
- Everything runs on your own computer; nothing is uploaded. If the Windows firewall asks, you may
  click Cancel: the app only listens on localhost.

STOP
- Click "Quit" in the small window, or close that window.

YOUR FILES
- Documents\Reticula\  holds exports, saved work, the sample database and your settings.
- Uninstalling the program does not delete that folder.

PROBLEMS
- The small window has "Open log" and "Reload app". Please attach the log to a bug report.
- Run without installing: see "Run from source" in the README on GitHub.
