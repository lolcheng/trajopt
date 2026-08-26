#!/bin/bash
# 快速开始脚本

set -e  # 遇到错误立即退出

echo "=========================================="
echo "TrajOpt C++ 快速开始"
echo "=========================================="

# 检查依赖
echo ""
echo "1. 检查依赖..."
if ! command -v cmake &> /dev/null; then
    echo "错误: CMake 未安装"
    echo "安装: sudo apt install cmake"
    exit 1
fi
echo "  ✓ CMake 已安装"

if ! command -v g++ &> /dev/null; then
    echo "错误: g++ 未安装"
    echo "安装: sudo apt install g++"
    exit 1
fi
echo "  ✓ g++ 已安装"

# 检查Ipopt
echo ""
echo "2. 检查Ipopt..."
if pkg-config --modversion ipopt &> /dev/null; then
    VERSION=$(pkg-config --modversion ipopt)
    echo "  ✓ Ipopt 已安装 (版本: $VERSION)"
else
    echo "  ✗ Ipopt 未找到"
    echo ""
    echo "请先安装Ipopt:"
    echo "  sudo apt update"
    echo "  sudo apt install coinor-libipopt-dev"
    echo ""
    read -p "是否现在安装? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        sudo apt update
        sudo apt install -y coinor-libipopt-dev
    else
        echo "请先安装Ipopt后再运行此脚本"
        exit 1
    fi
fi

# 编译
echo ""
echo "3. 编译项目..."
mkdir -p build
cd build

if [ ! -f CMakeCache.txt ]; then
    echo "  运行 cmake..."
    cmake ..
fi

echo "  运行 make..."
make

echo ""
echo "=========================================="
echo "编译完成！"
echo "=========================================="
echo ""
echo "运行程序:"
echo "  cd build"
echo "  ./trajopt_cpp"
echo ""

