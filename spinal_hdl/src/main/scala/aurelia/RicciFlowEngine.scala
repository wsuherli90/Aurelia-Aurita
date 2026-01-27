package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class RicciFlowEngine extends Component {
  val io = new Bundle {
    val clk = in Bool()
    val rst = in Bool()
    val mem_data_in   = slave Stream(SymmetricMatrix3())
    val mem_data_out  = master Stream(SymmetricMatrix3())
    val y_vec_in = slave Stream(Vector3())
    val voxel_count_out = out(UInt(32 bits))
    val grid_dim = in(Vec(UInt(32 bits), 3))
  }
  
  val area = new ClockingArea(ClockDomain(io.clk, io.rst)) {
      
    val raw_input = io.mem_data_in.translateWith(unpackSymmetric(io.mem_data_in.payload))
      
    val uWindow = new WindowBuffer3D()
    uWindow.io.input << raw_input
    uWindow.io.grid_dim := io.grid_dim
    
    case class WindowData() extends Bundle {
       val w = Window3D()
       val center = Matrix3()
    }
    
    val win_stream = Stream(WindowData())
    win_stream.valid := uWindow.io.output.valid
    win_stream.payload.w := uWindow.io.output.payload
    win_stream.payload.center := uWindow.io.output.payload(1)(1)(1)
    
    uWindow.io.output.ready := win_stream.ready

    val (st_inv, st_bypass) = StreamFork2(win_stream)
    
    val uInv = new MatrixInverse()
    uInv.io.input.valid := st_inv.valid
    uInv.io.input.payload := st_inv.payload.center
    st_inv.ready := uInv.io.input.ready
    
    val joined_inv_win = StreamJoin.arg(uInv.io.output, st_bypass)
    val (s_chern, s_lap, s_ric_bypass) = StreamFork3(joined_inv_win)
    
    val uChern = new ChernConnectionEngine()
    
    uChern.io.input_valid := s_chern.valid
    s_chern.ready := uChern.io.output_ready
    
    def mapAxis(w_axis: Vec[Matrix3]): Vec[Matrix3] = {
       val v = Vec(Matrix3(), 5)
       v(0) := w_axis(0)
       v(1) := w_axis(0)
       v(2) := w_axis(1)
       v(3) := w_axis(2)
       v(4) := w_axis(2)
       v
    }
    
    val w = s_chern.payload._2.w
    uChern.io.win_x := mapAxis(w(1)(1))
    uChern.io.win_y := mapAxis(Vec(w(1)(0)(1), w(1)(1)(1), w(1)(2)(1)))
    uChern.io.win_z := mapAxis(Vec(w(0)(1)(1), w(1)(1)(1), w(2)(1)(1)))
    
    uChern.io.g_inv_in := s_chern.payload._1
    uChern.io.y_vec := io.y_vec_in.payload
    io.y_vec_in.ready := s_chern.fire 
    
    val uLap = new LaplacianOrder4Unit()
    uLap.io.input_valid := s_lap.valid
    s_lap.ready := uLap.io.output_ready
    
    uLap.io.window_3d := s_lap.payload._2.w
    uLap.io.g_inv     := s_lap.payload._1
    
    val join_c_l = StreamJoin.arg(uChern.io.gamma_final, uLap.io.output_laplacian)
    val join_all_ricci = StreamJoin.arg(join_c_l, s_ric_bypass)
    
    val uRicciT = new RicciTensorCore()
    uRicciT.io.input_valid := join_all_ricci.valid
    join_all_ricci.ready := uRicciT.io.output_ready
    
    uRicciT.io.gamma       := join_all_ricci.payload._1._1
    uRicciT.io.laplacian_g := join_all_ricci.payload._1._2
    uRicciT.io.g_inv       := join_all_ricci.payload._2._1
    
    val (s_ricci_calc, s_scalar_bypass) = StreamFork2(join_all_ricci)
    
    uRicciT.io.input_valid := s_ricci_calc.valid
    s_ricci_calc.ready     := uRicciT.io.output_ready
    uRicciT.io.gamma       := s_ricci_calc.payload._1._1
    uRicciT.io.laplacian_g := s_ricci_calc.payload._1._2
    uRicciT.io.g_inv       := s_ricci_calc.payload._2._1
    
    val join_scalar = StreamJoin.arg(uRicciT.io.output_ricci, s_scalar_bypass)
    val uScalar = new RicciScalarUnit()
    
    val (s_to_scalar, s_bypass_metric) = StreamFork2(join_scalar)
    
    uScalar.io.input_valid := s_to_scalar.valid
    s_to_scalar.ready := uScalar.io.output_ready
    uScalar.io.ricci_tensor := s_to_scalar.payload._1
    uScalar.io.g_inv        := s_to_scalar.payload._2._2._1
    
    val join_chem = StreamJoin.arg(uScalar.io.output, s_bypass_metric)
    val (s_to_chem, s_bypass_stress) = StreamFork2(join_chem)
    
    val uChem = new ChemicalPotential()
    uChem.io.ricci_scalar << s_to_chem.translateWith(s_to_chem.payload._1)
    
    val join_stress = StreamJoin.arg(uChem.io.mu_out, s_bypass_stress)
    
    val uStress = new ActiveStressUnit()
    uStress.io.input_valid := join_stress.valid
    join_stress.ready := uStress.io.output_ready
    
    uStress.io.mu_in    := join_stress.payload._1
    uStress.io.g_metric := join_stress.payload._2._2._2._2.center 
    
    val join_final = StreamJoin.arg(uStress.io.output, s_bypass_stress)
    
    val s_update = join_final.translateWith {
        val t_stress = join_final.payload._1
        
        val bypass_lvl1 = join_final.payload._2
        val ricci_tensor_val = bypass_lvl1._2._1
        
        val g_old = bypass_lvl1._2._2._2._2.center
        
        val g_new = Matrix3()
        for(r <- 0 until 3; c <- 0 until 3) {
            val ricci_val = ricci_tensor_val(r)(c)
            val stress_val = t_stress(r)(c)
            val term = (RealConst(-2.0) * ricci_val) + stress_val
            g_new(r)(c) := g_old(r)(c) + (RealConst(TIME_DT_C) * term)
        }
        g_new
    }.s2mPipe()
    
    io.mem_data_out << s_update.translateWith(packSymmetric(s_update.payload))
    
    io.voxel_count_out := io.mem_data_out.fire.asUInt.resize(32 bits)
  }
}
