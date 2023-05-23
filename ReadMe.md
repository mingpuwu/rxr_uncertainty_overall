reactionnum中填写本次要处理的反应序列

Data中命名规范
文件夹名称：反应序列_反应式(反应式中=/<=>要替换成_)
文件夹中文件名称：
    web复制过来的正向反应的excel名称：反应序列_反应式.xlsx
    web复制过来的反向反应的excel名称：反应序列_反应式_reverse.xlsx

运行完1_datacorver.ipynb 后会在相应目录生成：反应序列_反应式.txt
运行完2_reverdataconvert.ipynb 后会在相应目录生成：反应序列_反应式_reverse.txt
运行完3_recal_reverse_rete.ipynb 后会更改反应序列_反应式_reverse.txt文件中的内容


4_reFit.ipynb  是用来对正向反应序列重新fit计算(忽略unused/taotao/henry/Baulch数据)
之后会生成新的反应序列_反应式.txt