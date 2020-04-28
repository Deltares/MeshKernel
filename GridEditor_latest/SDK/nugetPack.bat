SETLOCAL ENABLEDELAYEDEXPANSION

Set packageVersion="1.1.0.0"

rem build Grid Editor
NuGet pack GridEditorDeltaShellPlugin.nuspec -Properties "Configuration=Debug" -Version %packageVersion%
