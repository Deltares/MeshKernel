set libVersion=Debug
set guiVersion=%libVersion%
copy /y ".\%libVersion%\lib\MeshKernel.dll" "D:\ENGINES\DeltaShell\Grid_Editor_plugin\packages\Deltares.MeshKernel.1.0.0.36\native\Lib"
copy /y ".\%libVersion%\lib\MeshKernel.pdb" "D:\ENGINES\DeltaShell\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.pdb"
copy /y ".\%libVersion%\lib\MeshKernel.dll" "D:\ENGINES\DeltaShell\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
