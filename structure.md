# LiMDAの構成
## LiMDAのファイル構成
limda <br>
|--src <br> 
| &emsp; |--limda <br> 
| &emsp; | &emsp; |--SimulationFrame.py<br>
| &emsp; | &emsp; |--SimulationFrames.py<br>
| &emsp; | &emsp; |--analisys_frame.py<br>
| &emsp; | &emsp; |--calculate.py<br>
| &emsp; | &emsp; |--const.py<br>
| &emsp; | &emsp; |--export_frame.py<br>
| &emsp; | &emsp; |--export_frames.py]<br>
| &emsp; | &emsp; |--import_frame.py<br>
| &emsp; | &emsp; |--import_frames.py<br>
| &emsp; | &emsp; |--neighbor.cpp <br>
| &emsp; | &emsp; |--neighbor.pyx<br>
| &emsp; | &emsp; |--setup.py<br>
| &emsp; | &emsp; |--testcase_neighbor.py<br>
| &emsp; | <br>
| &emsp; |--scripts <br>
| &emsp; | &emsp; |--change_file_type.py <br>

## SimulationFrameのクラス構成
SimulationFrame <br>
&emsp;├ ImportFrame <br>
&emsp;├ ExportFrame <br>
&emsp;├ Calculate <br>
&emsp;├ AnalisysFrame <br>

## SimulationFramesのクラス構成
SimulationFrames <br>
&emsp;├ ImportFrames <br>
&emsp;├ ExportFrames <br>
