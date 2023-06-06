@ echo off

set name=meshkernel

del /f %name%.wixpdb > del.log 2>&1
del /f %name%.wixobj > del.log 2>&1
del /f %name%.msi > del.log 2>&1
del /f del.log 


set wixDir=c:\Program Files (x86)\WiX Toolset v3.11\bin


echo .
echo .
echo Candle %name%.wxs ...
"%wixDir%\candle.exe" %name%.wxs -ext WixUIExtension -out %name%.wixobj


echo .
echo .
echo Light %name%.wixobj ...
"%wixDir%\light.exe" %name%.wixobj -ext WixUIExtension -cultures:en-us -loc %name%.wxl -out %name%.msi

