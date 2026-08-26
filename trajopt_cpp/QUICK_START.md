# 快速开始指南

## 方法1：使用快速启动脚本（推荐）

```bash
cd /home/yizhe/trajopt/trajopt_cpp
./quick_start.sh
```

脚本会自动：
1. 检查依赖
2. 检查/安装Ipopt
3. 编译项目

## 方法2：手动安装和编译

### 步骤1：安装Ipopt

```bash
sudo apt update
sudo apt install coinor-libipopt-dev
```

### 步骤2：编译项目

```bash
cd /home/yizhe/trajopt/trajopt_cpp
mkdir -p build
cd build
cmake ..
make
```

### 步骤3：运行

```bash
./trajopt_cpp
```

## 验证安装

检查Ipopt是否安装成功：

```bash
pkg-config --modversion ipopt
```

应该输出版本号，例如：`3.14.x`

## 常见问题

### 问题：CMake找不到Ipopt

**解决方法1**：确保已安装
```bash
sudo apt install coinor-libipopt-dev
```

**解决方法2**：检查安装位置
```bash
find /usr -name "Ipopt*.hpp" 2>/dev/null
find /usr -name "libipopt*" 2>/dev/null
```

### 问题：编译错误

如果遇到链接错误，尝试：
```bash
# 更新库路径
sudo ldconfig

# 重新编译
cd build
rm -rf *
cmake ..
make
```

## 详细文档

更多信息请查看：
- `INSTALL.md` - 详细安装指南
- `README.md` - 项目说明
- `ARCHITECTURE.md` - 架构设计

