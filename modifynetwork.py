# -*- coding: utf-8 -*-
"""
Created on Sun May 25 22:47:04 2025

@author: ziyan
"""

with open('.\\saves\\nets\\jazz.txt', 'r') as f_in, open('.\\saves\\nets\\jazz-net.txt', 'w') as f_out:
    for line in f_in:
        parts = line.split()  # 默认按空白字符分割
        f_out.write(f"{parts[0]} {parts[1]}\n")  # 只取前两个数字