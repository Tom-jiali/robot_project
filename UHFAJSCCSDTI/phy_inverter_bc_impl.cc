#include "phy_inverter_bc_impl.h"
#include <gnuradio/io_signature.h>
#include <cstring>   
#include <algorithm> 
#include <initializer_list>

namespace gr {
namespace my_module {

// ---------------------------------------------------------------- //
// 构造函数：完成所有内存预分配与调度配置
// ---------------------------------------------------------------- //
phy_inverter_bc_impl::phy_inverter_bc_impl()
  : gr::block("phy_inverter_bc",
              gr::io_signature::make2(2, 2, sizeof(uint8_t), sizeof(int32_t)), // IN: port0(期望波形比特), port1(恶劣载波索引)
              gr::io_signature::make(1, 1, sizeof(uint8_t)))                   // OUT: port0(求解出的发送比特)
{
    d_bad_carriers_len = 8; // 最多丢弃 8 个载波

    // 预分配内存，杜绝运行时扩容
    d_x_solved.assign(BETA, 0);
    d_augmented_matrix.resize(BETA);
    d_dropped_k_buffer.reserve(48); // 8 个载波 * 6 bit

    // 强制调度器每次索要数据时，必须按 240 的整数倍来要
    set_output_multiple(BETA);

    // [TODO]: 在这里初始化你的 d_C_original 矩阵
}

phy_inverter_bc_impl::~phy_inverter_bc_impl()
{
}

// ---------------------------------------------------------------- //
// 调度预测函数：解决 288 进 240 出的速率失配死锁
// ---------------------------------------------------------------- //
void phy_inverter_bc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
{
    int n_symbols = (noutput_items + BETA - 1) / BETA; 
    ninput_items_required[0] = n_symbols * ALPHA;              
    ninput_items_required[1] = n_symbols * d_bad_carriers_len; 
}

// ---------------------------------------------------------------- //
// 零分配的逆交织计算
// ---------------------------------------------------------------- //
void phy_inverter_bc_impl::calculate_dropped_rows(const int32_t* bad_carriers, int len) 
{
    d_dropped_k_buffer.clear(); // O(1) 操作，复用已有内存

    const int n_cbps = ALPHA; // 288
    const int n_bpsc = 6;     // 64-QAM
    const int s = 3;          

    for (int idx = 0; idx < len; ++idx) {
        int c = bad_carriers[idx];
        if (c < 0 || c > 47) continue; // 脏数据防线

        for (int b = 0; b < n_bpsc; ++b) {
            int j = c * n_bpsc + b;
            int i = s * (j / s) + (j + (16 * j) / n_cbps) % s;
            int k = 16 * i - (n_cbps - 1) * ((16 * i) / n_cbps);
            d_dropped_k_buffer.push_back(k);
        }
    }

    std::sort(d_dropped_k_buffer.begin(), d_dropped_k_buffer.end());
    auto last = std::unique(d_dropped_k_buffer.begin(), d_dropped_k_buffer.end());
    d_dropped_k_buffer.erase(last, d_dropped_k_buffer.end());
}

// ---------------------------------------------------------------- //
// 极速版 GF(2) 高斯消元求解器
// ---------------------------------------------------------------- //
bool phy_inverter_bc_impl::solve_gf2_linear_system(std::vector<RowType>& augmented_matrix, 
                                                   std::vector<uint8_t>& x_out) 
{
    // 1. 高斯消元化为上三角矩阵
    for (int col = 0; col < BETA; ++col) {
        int pivot_row = col;
        // 寻找主元
        while (pivot_row < BETA && !augmented_matrix[pivot_row].test(col)) {
            pivot_row++;
        }

        // 数学陷阱防线：找不到主元，发生秩亏！立刻打断。
        if (pivot_row == BETA) {
            return false; 
        }

        if (pivot_row != col) {
            std::swap(augmented_matrix[col], augmented_matrix[pivot_row]);
        }

        for (int row = col + 1; row < BETA; ++row) {
            if (augmented_matrix[row].test(col)) {
                augmented_matrix[row] ^= augmented_matrix[col]; 
            }
        }
    }

    // 2. 回代求解
    for (int row = BETA - 1; row >= 0; --row) {
        uint8_t sum = augmented_matrix[row][BETA]; 
        for (int col = row + 1; col < BETA; ++col) {
            if (augmented_matrix[row].test(col) && x_out[col]) {
                sum ^= 1; 
            }
        }
        x_out[row] = sum;
    }
    return true;
}

// ---------------------------------------------------------------- //
// DSP 核心工作主循环
// ---------------------------------------------------------------- //
int phy_inverter_bc_impl::general_work(int noutput_items,
                                       gr_vector_int &ninput_items,
                                       gr_vector_const_void_star &input_items,
                                       gr_vector_void_star &output_items)
{
    // 指针类型对齐
    const uint8_t *y_target = (const uint8_t *) input_items[0];     
    const int32_t *bad_carriers = (const int32_t *) input_items[1]; 
    uint8_t *x_out = (uint8_t *) output_items[0];                   

    int available_in_y = ninput_items[0] / ALPHA;
    int available_in_bad = ninput_items[1] / d_bad_carriers_len;
    int available_out = noutput_items / BETA;
    
    // 取短板作为本次批量处理的符号数
    int n_symbols = std::min({available_in_y, available_in_bad, available_out});

    if (n_symbols == 0) return 0; 

    for (int s = 0; s < n_symbols; ++s) {
        int in_offset_y = s * ALPHA;
        int in_offset_bad = s * d_bad_carriers_len;
        int out_offset_x = s * BETA;

        // 直接传裸指针，O(1) 更新 d_dropped_k_buffer
        calculate_dropped_rows(bad_carriers + in_offset_bad, d_bad_carriers_len);
        
        bool drop_mask[ALPHA] = {false};
        for (int row_idx : d_dropped_k_buffer) {
            if (row_idx >= 0 && row_idx < ALPHA) drop_mask[row_idx] = true;
        }

        int matrix_row = 0;
        for(int i = 0; i < ALPHA; ++i) {
            if(!drop_mask[i]) {
                // O(1) 极速拷贝整行基底（避免按 bit 赋值导致的 Cache Miss）
                d_augmented_matrix[matrix_row] = d_C_original[i]; 
                
                // 单独设定增广列 y'
                if (y_target[in_offset_y + i]) {
                    d_augmented_matrix[matrix_row].set(BETA);
                } else {
                    d_augmented_matrix[matrix_row].reset(BETA); 
                }
                
                matrix_row++;
                if (matrix_row == BETA) break; // 拿满退出
            }
        }

        // 鲁棒性盲区防线 1：有效行数提取不足
        if (matrix_row < BETA) {
            std::memset(x_out + out_offset_x, 0, BETA);
            continue; 
        }

        // 极速求解 (内部包含鲁棒性防线 2：秩亏检测)
        if (solve_gf2_linear_system(d_augmented_matrix, d_x_solved)) {
            for(int i = 0; i < BETA; ++i){
                x_out[out_offset_x + i] = d_x_solved[i];
            }
        } else {
            // 发生秩亏，发送全 0，交给后续 JSCC 作为噪声吸收
            std::memset(x_out + out_offset_x, 0, BETA);
        }
    }

    // 报告消耗的输入量
    consume(0, n_symbols * ALPHA); 
    consume(1, n_symbols * d_bad_carriers_len);
    
    return n_symbols * BETA;
}

} // namespace my_module
} // namespace gr