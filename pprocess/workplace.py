from extract_for_plot import load_runs, list_observables
from extract_for_plot import plot_field, plot_scalar_bar, plot_across_runs,plot_g2_overlay
import matplotlib.pyplot as plt
data = load_runs("/home/zhangxiangyu/project/all_in_one_git/results/data/out_files")          # 目录里含 reading_guide.nml & b*.nml

# 单个 run：画场
plot_field(data, run="b1", name="corf_G1")
plt.show()
# 单个 run：所有标量柱状图
plot_scalar_bar(data, run="b1")
plt.show()
# 跨 b1..b8：横轴自定义（你的 1×8 数组），标量
x_custom = [0.1,0.2,0.3,0.4]
plot_across_runs(data, name="BF_KE", x=x_custom, xtick_labels=x_custom, xlabel="lambda")
plt.show()
# 跨 b1..b8：场的第 1 个 index（1-based）
plot_across_runs(data, name="corf_G1", field_index=1, x=x_custom, xtick_labels=x_custom, xlabel="lambda")
plt.show()
# 跨 b1..b8：场的第 1 个 index（1-based），G2
plot_g2_overlay(data) 
plot_g2_overlay(data, loglog=True)
