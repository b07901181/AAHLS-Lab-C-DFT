# AAHLS-Lab-C-DFT
Implement 1024-size DFT by FFT
## Problem
The problem comes from [Customized Computing Challenge](https://xupsh.github.io/ccc/introduction.html), and this project implements DFT in first part. The detail of problem can be found in the [GitHub link](https://github.com/xupsh/ccc/tree/main/problems/DFT).<br>
The problem is to implement a 1024-size DFT, and we need to pass Csim, Csyn and Cosim. The score of the problem is estimated frequency divided by cycle count.
## Introduction to the Overall System
In this system we perform Radix-2 FFT, and we use a special form for speedup our process. Please refer to **report.pdf** and the [FFT link](https://zh.wikipedia.org/zh-tw/%E5%BA%93%E5%88%A9%EF%BC%8D%E5%9B%BE%E5%9F%BA%E5%BF%AB%E9%80%9F%E5%82%85%E9%87%8C%E5%8F%B6%E5%8F%98%E6%8D%A2%E7%AE%97%E6%B3%95#%E5%96%AE%E4%B8%80%E5%9F%BA%E5%BA%95) for more details.
## Building
Run the following code in command line (or Vitis HLS Command Prompt):
```cmd
vitis_hls -p tcl_script.tcl
```
Or you can build the project following [Lab1](https://github.com/bol-edu/course-lab_1/blob/2022.1/2022.1-Workbook-Lab1.pdf).
## Result
Use 7ns-clock for Csynthesis and Cosimulation, you will get score 6063.81.
