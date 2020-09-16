set libVersion=Debug
set guiVersion=%libVersion%
copy /y ".\%libVersion%\lib\MeshKernel.dll" "C:\engines\DeltaShell\Grid_editor_front_end\Grid_editor_plugin_v03\Grid_Editor_plugin\src\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
copy /y ".\%libVersion%\lib\MeshKernel.pdb" "C:\engines\DeltaShell\Grid_editor_front_end\Grid_editor_plugin_v03\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.pdb"
copy /y ".\%libVersion%\lib\MeshKernel.dll" "C:\engines\DeltaShell\Grid_editor_front_end\Grid_editor_plugin_v03\Grid_Editor_plugin\bin\%guiVersion%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
