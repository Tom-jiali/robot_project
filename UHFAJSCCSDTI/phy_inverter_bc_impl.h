#ifndef INCLUDED_PHY_INVERTER_BC_IMPL_H
#define INCLUDED_PHY_INVERTER_BC_IMPL_H

#include <gnuradio/block.h>
#include <bitset>
#include <vector>
#include <cstdint>

// 替换为你自己的命名空间
namespace gr {
namespace my_module {

    class phy_inverter_bc_impl : public gr::block
    {
    private:
        static constexpr int BETA = 240;  // 5/6码率下的未编码比特数
        static constexpr int ALPHA = 288; // 编码后的总比特数
        
        // 增广矩阵的一行：240位系数 + 1位常数项
        using RowType = std::bitset<BETA + 1>;

        // ----- 预分配的复用缓存 (彻底消灭 new/delete) -----
        std::vector<uint8_t> d_x_solved;             // 求解结果缓存
        std::vector<RowType> d_augmented_matrix;     // 增广矩阵缓存
        std::vector<int> d_dropped_k_buffer;         // 逆交织索引缓存
        
        int d_bad_carriers_len;

        // 原始 C 矩阵直接存储为 RowType 向量，最高位(BETA位)默认留空/置零
        // 注意：你需要在模块初始化时为它赋值
        std::vector<RowType> d_C_original; 

        // ----- 私有核心算法 -----
        // 接收裸指针和长度，直接操作 d_dropped_k_buffer
        void calculate_dropped_rows(const int32_t* bad_carriers, int len);

        // 高效高斯消元求解器
        bool solve_gf2_linear_system(std::vector<RowType>& augmented_matrix, 
                                     std::vector<uint8_t>& x_out);

    public:
        phy_inverter_bc_impl();
        ~phy_inverter_bc_impl();

        // GNU Radio 调度器必须重写的预测函数
        void forecast(int noutput_items, gr_vector_int &ninput_items_required) override;

        // DSP 主循环
        int general_work(int noutput_items,
                         gr_vector_int &ninput_items,
                         gr_vector_const_void_star &input_items,
                         gr_vector_void_star &output_items) override;
    };

} // namespace my_module
} // namespace gr

#endif /* INCLUDED_PHY_INVERTER_BC_IMPL_H */