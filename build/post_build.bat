set libVersion=Debug
set guiVersion=Debug
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "..\GridEditor\src\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.pdb" "..\GridEditor\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.pdb"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "..\GridEditor\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "C:\engines\DeltaShell\GridEditor_latest\src\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"
copy /y ".\%libVersion%\gridgeomStateful_dll.pdb" "C:\engines\DeltaShell\GridEditor_latest\bin\Debug\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.pdb"
copy /y ".\%libVersion%\gridgeomStateful_dll.dll" "C:\engines\DeltaShell\GridEditor_latest\bin\Debug\plugins\DeltaShell.Plugins.GridEditor\Lib\gridgeomStateful_dll.dll"


