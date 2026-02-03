# Aurelia Accelerator Constraints
# Target Frequency: 250 MHz (Period 4.0ns)

create_clock -period 4.000 -name clk [get_ports clk]

# Input Delays (Assuming somewhat tight coupling)
set_input_delay -clock clk -max 1.000 [get_ports {s_axis_data_tvalid s_axis_data_tdata[*] s_axis_data_tlast}]
set_input_delay -clock clk -min 0.000 [get_ports {s_axis_data_tvalid s_axis_data_tdata[*] s_axis_data_tlast}]

# Output Delays
set_output_delay -clock clk -max 1.000 [get_ports {m_axis_data_tvalid m_axis_data_tdata[*] m_axis_data_tlast}]
set_output_delay -clock clk -min 0.000 [get_ports {m_axis_data_tvalid m_axis_data_tdata[*] m_axis_data_tlast}]

# Relax timing on static config inputs if necessary
set_false_path -from [get_ports {grid_dim[*]}]
