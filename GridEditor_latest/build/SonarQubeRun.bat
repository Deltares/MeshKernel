set svnDir=%1
set msBuildDir=%2

"%svnDir%\build\SonarQube\SonarScanner.MSBuild.exe" begin /k:"Grid_Editor_trunk" /n:"Grid Editor (trunk)" /d:sonar.scm.provider=svn /d:sonar.host.url="https://sonarqube.deltares.nl" /d:sonar.login="193864686517364e4fd68bf4573f24f98990f0dd" /d:sonar.cs.dotcover.reportsPaths="%svnDir%\bin\Debug\dotCover.html"

%msBuildDir%\MsBuild.exe GridEditor.sln /t:Rebuild /t:clean

%svnDir%\build\dotCover_2019.1.2\dotCover.exe cover /TargetExecutable="%svnDir%\packages\NUnit.Runners.2.7.0\tools\nunit-console.exe" /Output="%svnDir%\bin\Debug\DotCover.html" /ReportType="HTML" /TargetArguments="%svnDir%\bin\Debug\DeltaShell.Plugins.GridEditor.Gui.Tests.dll %svnDir%\bin\Debug\DeltaShell.Plugins.GridEditor.Tests.dll /noshadow /exclude:Build.WorkInProgress

"%svnDir%\build\SonarQube\SonarScanner.MSBuild.exe" end /d:sonar.login="193864686517364e4fd68bf4573f24f98990f0dd"

