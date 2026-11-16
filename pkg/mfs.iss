#define AppName "Reticula"
; Build:  iscc /DAppVer=2.13.0 /DSourceDir=..\dist\Reticula pkg\mfs.iss     (the release workflow does this)
#ifndef AppVer
  #define AppVer "0.0.0"
#endif
#ifndef SourceDir
  #define SourceDir "..\dist\Reticula"
#endif
; v2.11: single-file installer (one Reticula_Setup_<ver>.exe, no .bin parts) - the update-check asset on GitHub
[Setup]
AppId={{6B1E4C2A-9D3F-4E7B-A1C5-6F2A9B3C1D40}
AppName={#AppName}
AppVersion={#AppVer}
AppPublisher=Stijn Germonpré
AppComments=Publication figures for framework materials
VersionInfoDescription=Reticula - Publication figures for framework materials
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir={#SourcePath}\..\dist_installer
OutputBaseFilename=Reticula_Setup_{#AppVer}
SetupIconFile={#SourcePath}\..\icon\reticula.ico
UninstallDisplayIcon={app}\Reticula.exe
Compression=lzma2/ultra64
SolidCompression=yes
DiskSpanning=no
WizardStyle=modern
ChangesAssociations=yes
DisableProgramGroupPage=yes
[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: checkedonce
Name: "assoc"; Description: "Open .reticula project files with Reticula"; Flags: checkedonce
[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#SourcePath}\README_Windows.txt"; DestDir: "{app}"; Flags: ignoreversion
[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\Reticula.exe"; IconFilename: "{app}\Reticula.exe"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\Reticula.exe"; Tasks: desktopicon
[Registry]
; v2.11: .reticula project files (per user, no admin rights needed)
Root: HKA; Subkey: "Software\Classes\.reticula"; ValueType: string; ValueName: ""; ValueData: "Reticula.Project"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Reticula.Project"; ValueType: string; ValueName: ""; ValueData: "Reticula project"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Reticula.Project\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\Reticula.exe,0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Reticula.Project\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\Reticula.exe"" ""%1"""; Tasks: assoc
[Run]
Filename: "{app}\Reticula.exe"; Description: "Start Reticula now"; Flags: nowait postinstall skipifsilent
[Messages]
WelcomeLabel2=This will install [name/ver] on your computer.%n%nBuilt by Stijn Germonpré%n%nIt is recommended that you close all other applications before continuing.
