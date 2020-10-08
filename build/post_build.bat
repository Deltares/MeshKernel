set gridEditorPlugin=D:\ENGINES\DeltaShell\Grid_Editor_plugin
set lib=Release
set gui=Release
copy /y ".\%lib%\lib\MeshKernel.dll" "%gridEditorPlugin%\packages\Deltares.MeshKernel.1.0.0.37\native\Lib"
copy /y ".\%lib%\lib\MeshKernel.pdb" "%gridEditorPlugin%\bin\%gui%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.pdb"
copy /y ".\%lib%\lib\MeshKernel.dll" "%gridEditorPlugin%\bin\%gui%\plugins\DeltaShell.Plugins.GridEditor\Lib\MeshKernel.dll"
