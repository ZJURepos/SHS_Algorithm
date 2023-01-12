# SHS Algorithm

This project contains a third party library and custom library. For the convenience of users, the main interface of the programs is described here. All file paths should be in English. It should be noted that if you want to run the command program, you must finish the development environment configuration firstly. If not, you can download the library files from github repository "https://github.com/ZJURepos/SHS_Algorithm.git". All the code will be open source after the paper is published. 

## Skeletonization

The "main.cpp" is the interface program of Skeletonization. For user convenience, we provide a command-line program "Skeletonization.exe", it reads file's path and output path from configuration file "ske_id.ini".  The users just only change the input file and output path. The skeleton file will be named as "xxx-ske.pcd" or "xxx-ske.ply" .

## DBSCAN 

We provide c++ and python versions of DBSCAN. "DBSCAN.exe" is the c++ version, "DBscan.py" is the python version. If use our algorithm, we suggest use c++ version, it was applied for our algorithm. The python version also can be used, but the codes we provide just remove the content. If users want to recode the algorithm, the python version may be helpful. In order to use c++ DBSCAN, you need to use "File2Txt.py" to change file's format. Because the DBSCAN code just reads and writes ".txt" file. "DBSCAN.exe" will get file path and out put path in "DBscan.ini". 

## Skeleton optimization 

Users can use "skeConnect&Opti&length.py" for skeleton connection and silique length calculation. And users should provide 3 files, source point cloud file, skeleton file and skeleton-DBscan file. For example, if you have cultivar ''xxx'', you should prepare files "xxx.ply", "xxx-ske.ply", "xxx-ske-DB.txt". The program will get file automatically. And users just need to change the paths in lines 1578-1579 in "skeConnect&Opti&length.py". We will get output the classified file named "xxx-classified.txt". You also can get length information in python command line. 

 