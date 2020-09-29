set libVersion=Release
set guiVersion=%libVersion%
copy /y ".\%libVersion%\lib\MeshKernel.dll" "C:\engines\Grid_Editor_plugin\src\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
copy /y ".\%libVersion%\lib\MeshKernel.pdb" "C:\engines\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.pdb"
copy /y ".\%libVersion%\lib\MeshKernel.dll" "C:\engines\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
