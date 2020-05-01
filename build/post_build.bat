set libVersion=Debug
rem set guiVersion=Release
set guiVersion=%libVersion%
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "..\GridEditor_latest\src\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.pdb" "..\GridEditor_latest\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.pdb"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "..\GridEditor_latest\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "C:\engines\DeltaShell\GridEditor_latest\src\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.pdb" "C:\engines\DeltaShell\GridEditor_latest\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.pdb"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "C:\engines\DeltaShell\GridEditor_latest\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"


