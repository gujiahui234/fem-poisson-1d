# 一维 Poisson 方程的 P1 有限元实验

下载 Release 的 ZIP，解压后打开终端。需要 Python 3.12；不需要原作者的目录、Anaconda、Word 或课程讲义。

Windows CMD：

```cmd
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.venv\Scripts\python.exe reproduce.py
```

macOS / Linux：

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python reproduce.py
```

看到 `REPRODUCTION PASSED` 表示实验检查、测试和参考数据比较都通过。首次安装依赖需要联网；计算本身可离线完成。

程序实际组装并求解两个算例的有限元线性系统。默认单元数为 10、20、40、80、160，输出最大误差、L2 误差、能量半范数误差及收敛阶。`results/` 保存 CSV、图表、节点数组和检查记录。

一张图比较两单元解与精确解，另一张图只画最大误差随单元数的变化。L2 和能量误差保留在 CSV 中。

默认复现包含 118 项数值检查和 23 项测试。参考数据保存在 `examples/expected_errors.csv`；跨环境采用数值容差比较，不要求图像文件逐字节相同。当前本地验证环境为 Windows、Python 3.12；Linux 和 Windows 的远程验证结果见 GitHub Actions。

其他网格和算例的命令、误差计算方法及适用范围见 [英文说明](README.md)。这份代码限于两个光滑的一维算例；完整理论证明和二维研究不在此版本中。
