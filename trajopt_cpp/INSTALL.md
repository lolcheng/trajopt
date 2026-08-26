# 安装和编译指南

## 1. 安装Ipopt

### 方法1：使用apt安装（推荐，简单快速）

```bash
sudo apt update
sudo apt install coinor-libipopt-dev
```

### 方法2：从源码编译（如果需要最新版本或自定义配置）

如果apt安装的版本不满足需求，可以从源码编译：

```bash
# 安装依赖
sudo apt install gfortran liblapack-dev libblas-dev libmumps-dev

# 下载Ipopt源码
cd ~
wget https://github.com/coin-or/Ipopt/releases/download/releases%2F3.14.12/Ipopt-3.14.12.tgz
tar -xzf Ipopt-3.14.12.tgz
cd Ipopt-3.14.12

# 编译（这可能需要较长时间）
./configure --prefix=/usr/local
make
sudo make install

# 更新库路径
sudo ldconfig
```

### 验证安装

```bash
# 检查Ipopt库
pkg-config --modversion ipopt

# 或者检查头文件
ls /usr/include/coin-or/Ipopt*.hpp
```

## 2. 编译项目

### 步骤1：进入项目目录

```bash
cd /home/yizhe/trajopt/trajopt_cpp
```

### 步骤2：创建构建目录

```bash
mkdir -p build
cd build
```

### 步骤3：配置CMake

```bash
cmake ..
```

如果CMake找不到Ipopt，可能需要指定路径：

```bash
# 如果Ipopt安装在标准位置
cmake ..

# 如果Ipopt安装在自定义位置（例如 /usr/local）
cmake -DIPOPT_DIR=/usr/local ..
```

### 步骤4：编译

```bash
make
```

### 步骤5：运行

```bash
./trajopt_cpp
```

## 3. 常见问题

### 问题1：CMake找不到Ipopt

**错误信息**：
```
CMake Error: Could not find a package configuration file provided by "Ipopt"
```

**解决方法**：

1. 检查Ipopt是否安装：
```bash
pkg-config --modversion ipopt
```

2. 如果已安装但CMake找不到，手动指定路径：
```bash
# 查找Ipopt安装位置
find /usr -name "Ipopt*.hpp" 2>/dev/null

# 在CMakeLists.txt中手动设置路径（见下方）
```

### 问题2：链接错误

**错误信息**：
```
undefined reference to `Ipopt::...`
```

**解决方法**：

1. 检查库文件是否存在：
```bash
find /usr -name "libipopt*" 2>/dev/null
```

2. 更新CMakeLists.txt，手动指定库路径（见下方）

### 问题3：头文件找不到

**错误信息**：
```
fatal error: IpTNLP.hpp: No such file or directory
```

**解决方法**：

1. 查找头文件位置：
```bash
find /usr -name "IpTNLP.hpp" 2>/dev/null
```

2. 更新CMakeLists.txt中的包含路径

## 4. 修改CMakeLists.txt（如果需要）

如果CMake自动检测失败，可以手动指定路径。编辑 `CMakeLists.txt`：

```cmake
# 手动设置Ipopt路径（根据实际安装位置调整）
set(IPOPT_INC_DIR "/usr/include/coin-or")
set(IPOPT_LIB_DIR "/usr/lib")

# 或者使用pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(IPOPT REQUIRED ipopt)

include_directories(${IPOPT_INCLUDE_DIRS})
link_directories(${IPOPT_LIBRARY_DIRS})
```

## 5. 快速测试脚本

创建一个测试脚本 `test_build.sh`：

```bash
#!/bin/bash
set -e

echo "=== 检查依赖 ==="
which cmake || { echo "CMake not found!"; exit 1; }
which g++ || { echo "g++ not found!"; exit 1; }

echo "=== 检查Ipopt ==="
if pkg-config --modversion ipopt > /dev/null 2>&1; then
    echo "Ipopt version: $(pkg-config --modversion ipopt)"
else
    echo "Warning: Ipopt not found via pkg-config"
    echo "Trying to find manually..."
    find /usr -name "Ipopt*.hpp" 2>/dev/null | head -1 || echo "Ipopt headers not found!"
fi

echo "=== 编译项目 ==="
mkdir -p build
cd build
cmake ..
make

echo "=== 编译完成 ==="
echo "运行: ./trajopt_cpp"
```

运行：
```bash
chmod +x test_build.sh
./test_build.sh
```

