#!/usr/bin/env bash
# 运行分段轨迹优化并可视化结果
# 用法: ./run_and_visualize.sh [地形数据目录] [输出目录]
# 示例: ./run_and_visualize.sh ../../test/terrain_res seg_out

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
TERRAIN_DIR="${1:-../../test/terrain_res}"
OUT_DIR="${2:-seg_out}"

cd "$BUILD_DIR" || exit 1
if [[ ! -x ./seg_trajopt_cpp ]]; then
    echo "错误: 未找到可执行文件 $BUILD_DIR/seg_trajopt_cpp，请先编译 (cd build && cmake .. && make)"
    exit 1
fi

echo "===== 运行轨迹优化 ====="
echo "地形数据: $TERRAIN_DIR"
echo "输出目录: $OUT_DIR"
echo ""
./seg_trajopt_cpp "$TERRAIN_DIR" "$OUT_DIR"
EXIT_OPT=$?

# 0=成功, 1=达到最大迭代但有解，均继续做可视化
if [[ $EXIT_OPT -ne 0 ]] && [[ $EXIT_OPT -ne 1 ]]; then
    echo "优化异常退出 (code $EXIT_OPT)，跳过可视化"
    exit $EXIT_OPT
fi

if [[ ! -f "${OUT_DIR}/final_trajectory.txt" ]] && [[ ! -f "${OUT_DIR}/initial_trajectory.txt" ]]; then
    echo "未找到轨迹文件，跳过可视化"
    exit $EXIT_OPT
fi

echo ""
echo "===== 生成可视化 ====="
if ! python3 "${SCRIPT_DIR}/scripts/plot_seg_trajectory.py" "$OUT_DIR" --no-show; then
    echo "可视化失败"
    exit 2
fi
echo ""
echo "完成。图片: ${BUILD_DIR}/${OUT_DIR}/seg_trajectory_plot.png"
echo "      运动学: ${BUILD_DIR}/${OUT_DIR}/seg_trajectory_kinematics.png"
echo "如需弹窗查看 3D 视角，可单独运行: python3 ${SCRIPT_DIR}/scripts/plot_seg_trajectory.py $OUT_DIR"
exit 0
