# Mandelbrot Set
北航软件学院研究生-2020春并行程序设计课程（邵兵）-大作业
## 并行求解曼德博集合介绍

[曼德博集合(或曼德博复数集合)](https://baike.baidu.com/item/%E6%9B%BC%E5%BE%B7%E5%8D%9A%E9%9B%86%E5%90%88/5831843?fr=aladdin)是一种在复平面上组成分形的点的集合。

需要提供以下四个版本的优化程序：
- 纯MPI静态调度(static scheduling)，命名为mpi_static
- 纯MPI动态调度(dynamic scheduling)，命名为mpi_dynamic
- 纯OpenMP优化，命名为omp
- MPI+OpenMP混合优化，命名为hybrid

## 运行代码

要求：安装好MPI与OpenMP

```sh
chmod +x run.sh
```

运行串行版本的代码：
```sh
make sequential
time ./run.sh sequential
```

运行所有并行版本的代码：
```sh
make all
time ./run.sh mpi_static
time ./run.sh mpi_dynamic
time ./run.sh omp
time ./run.sh hybrid
```

debug模式，调试结果是否正确
```sh
make all_debug
time ./run.sh mpi_static_debug
time ./run.sh mpi_dynamic_debug
time ./run.sh omp_debug
time ./run.sh hybrid_debug
```

## 注意事项
不保证代码结果最优