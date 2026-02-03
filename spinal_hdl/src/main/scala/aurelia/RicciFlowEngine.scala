package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class RicciFlowEngine(config: AureliaConfig = AureliaConfig()) extends Component {
  val io = new Bundle {
    val clk = in Bool()
    val rst = in Bool()
    val mem_data_in   = slave Stream(SymmetricMatrix3())
    val mem_data_out  = master Stream(SymmetricMatrix3())
    val y_vec_in      = slave Stream(Vector3())
    val voxel_count_out = out(UInt(32 bits))
    val grid_dim      = in(Vec(UInt(32 bits), 3))
    
    val axion_out     = master Stream(RealType())
  }
  
  val area = new ClockingArea(ClockDomain(io.clk, io.rst)) {
      
    val raw_input = io.mem_data_in.translateWith(unpackSymmetric(io.mem_data_in.payload))
      
    // Dual Buffer System
    
    val uWindowG = new WindowBuffer(Matrix3(), config)
    uWindowG.io.input << raw_input
    uWindowG.io.grid_dim := io.grid_dim
    
    val uWindowY = new WindowBuffer(Vector3(), config)
    uWindowY.io.input << io.y_vec_in
    uWindowY.io.grid_dim := io.grid_dim
    
    case class WindowData() extends Bundle {
       val w = Window5D()
       val center_g = Matrix3()
       val center_y = Vector3()
    }
    
    val win_stream = Stream(WindowData())
    
    // Join window G and window Y.
    // Assuming identical latency (guaranteed by identical logic and grid dims)
    // We trust both are valid at same cycle if inputs are sync.
    // If inputs drift, we need StreamJoin. Here we assume System DMA feeds both sync.
    
    win_stream.valid := uWindowG.io.output.valid && uWindowY.io.output.valid 
    uWindowG.io.output.ready := win_stream.ready
    uWindowY.io.output.ready := win_stream.ready
    
    win_stream.payload.w := uWindowG.io.output.payload
    win_stream.payload.center_g := uWindowG.io.output.payload(2)(2)(2)
    win_stream.payload.center_y := uWindowY.io.output.payload(2)(2)(2)

    val (st_to_zss, st_main) = StreamFork2(win_stream)
    
    def slice3x3(w5: Window5D): Window3D = {
       val w3 = Window3D()
       for(z <- 0 until 3; y <- 0 until 3; x <- 0 until 3) {
           w3(z)(y)(x) := w5(z+1)(y+1)(x+1)
       }
       w3
    }
    
    val uZSS = new ZeroSkippingScheduler()
    uZSS.io.input_window << st_to_zss.translateWith(slice3x3(st_to_zss.payload.w))
    uZSS.io.threshold := RealConst(1.0e-5) 
    
    val (st_inv, st_bypass) = StreamFork2(st_main)
    
    val uInv = new MatrixInverse()
    uInv.io.input.valid := st_inv.valid
    uInv.io.input.payload := st_inv.payload.center_g
    st_inv.ready := uInv.io.input.ready
    
    // JOIN & DATA GATHERING
    // We define bundles to carry context cleanly to avoid tuple hell.
    
    // Join Inverse Output with Original Window Data
    
    case class InvAndWindow() extends Bundle {
       val inv_g   = Matrix3()
       val det     = RealType()
       val windata = WindowData()
    }
    
    val joined_inv_win = StreamJoin.arg(uInv.io.output, st_bypass).translateWith {
        val res = InvAndWindow()
        res.inv_g   := uInv.io.output.payload._1
        res.det     := uInv.io.output.payload._2
        res.windata := st_bypass.payload
        res
    }
    
    val (s_chern, s_lap, s_ric_bypass) = StreamFork3(joined_inv_win)
    
    val uChern = new ChernConnectionEngine(config)
    
    uChern.io.input_valid := s_chern.valid
    s_chern.ready := uChern.io.output_ready
    
    val w_g = s_chern.payload.windata.w
    uChern.io.win_x := w_g(2)(2) 
    uChern.io.win_y := Vec(w_g(2)(0)(2), w_g(2)(1)(2), w_g(2)(2)(2), w_g(2)(3)(2), w_g(2)(4)(2))
    uChern.io.win_z := Vec(w_g(0)(2)(2), w_g(1)(2)(2), w_g(2)(2)(2), w_g(3)(2)(2), w_g(4)(2)(2))
    
    uChern.io.g_inv_in := s_chern.payload.inv_g
    uChern.io.y_vec := s_chern.payload.windata.center_y 
    
    val uLap = new LaplacianOrder4Unit(config)
    uLap.io.input_valid := s_lap.valid
    s_lap.ready := uLap.io.output_ready
    
    uLap.io.window_3d := s_lap.payload.windata.w
    uLap.io.g_inv     := s_lap.payload.inv_g
    
    // Join Chern and Laplacian results
    val join_c_l = StreamJoin.arg(uChern.io.gamma_final, uLap.io.output_laplacian)
    
    // Join All Data for Ricci Calculation
    // We need: Gamma, Laplacian, Inverse, Determinant, SkipCmd
    
    case class RicciCalcContext() extends Bundle {
       val gamma   = GammaTensor() // Assuming type exists or is Matrix3
       val lap     = Matrix3()
       val context = InvAndWindow()
       val skip    = Bool()
    }
    
    // Step 1: Join (Gamma,Lap) [join_c_l] with Context [s_ric_bypass, which is InvAndWindow]
    // Add FIFO on bypass path to prevent head-of-line blocking from the latency-heavy Laplacian/Chern Units
    val join_calc_ctx = StreamJoin.arg(join_c_l, s_ric_bypass.queue(16))
    
    // Step 2: Join with SkipCmd
    // SkipCmd comes from ZSS which is fast, but we should act defensively.
    val join_with_skip = StreamJoin.arg(join_calc_ctx, uZSS.io.skip_cmd.queue(16)).translateWith {
        val res = RicciCalcContext()
        res.gamma   := join_calc_ctx.payload._1._1
        res.lap     := join_calc_ctx.payload._1._2
        res.context := join_calc_ctx.payload._2
        res.skip    := uZSS.io.skip_cmd.payload
        res
    }
    
    val uRicciT = new RicciTensorCore()
    
    val (s_to_ricci, s_context_bypass) = StreamFork2(join_with_skip)
    
    uRicciT.io.input_valid := s_to_ricci.valid
    s_to_ricci.ready := uRicciT.io.output_ready 
    
    // Clean, readable assignments
    uRicciT.io.gamma       := s_to_ricci.payload.gamma
    uRicciT.io.laplacian_g := s_to_ricci.payload.lap
    uRicciT.io.g_inv       := s_to_ricci.payload.context.inv_g
    uRicciT.io.det_in      := s_to_ricci.payload.context.det
    uRicciT.io.skip_calc   := s_to_ricci.payload.skip
    
    // Output side
    case class RicciOutContext() extends Bundle {
       val ricci   = Matrix3()
       val det_out = RealType()
       val context = InvAndWindow() // Original G, Y, etc
    }
    
    val join_scalar_context = StreamJoin.arg(uRicciT.io.output_ricci, s_context_bypass.queue(16))
    val join_det_out = StreamJoin.arg(join_scalar_context, uRicciT.io.det_out).translateWith {
        val res = RicciOutContext()
        res.ricci   := join_scalar_context.payload._1
        res.context := join_scalar_context.payload._2.context
        res.det_out := uRicciT.io.det_out.payload
        res
    }
    
    val uScalar = new RicciScalarUnit()
    val (s_to_scalar, s_bypass_metric) = StreamFork2(join_det_out)
    
    uScalar.io.input_valid := s_to_scalar.valid
    s_to_scalar.ready := uScalar.io.output_ready
    uScalar.io.ricci_tensor := s_to_scalar.payload.ricci
    uScalar.io.g_inv        := s_to_scalar.payload.context.inv_g
    
    // AXION UNIT
    val uAxion = new AxionFieldUnit()
    val (s_scalar_to_chem, s_scalar_to_axion) = StreamFork2(uScalar.io.output)
    
    uAxion.io.ricci_scalar << s_scalar_to_axion
    
    // Det is in s_bypass_metric payload (RicciOutContext)
    val (s_bm_stress, s_bm_axion) = StreamFork2(s_bypass_metric)
    
    uAxion.io.det_g << s_bm_axion.translateWith(s_bm_axion.payload.det_out)
    io.axion_out    << uAxion.io.axion_out
    
    // Chemical Potential
    // Join Scalar Output with Context
    
    case class ChemContext() extends Bundle {
       val scalar = RealType()
       val ctx    = RicciOutContext()
    }
    
    val join_chem = StreamJoin.arg(s_scalar_to_chem, s_bm_stress).translateWith {
        val res = ChemContext()
        res.scalar := s_scalar_to_chem.payload
        res.ctx    := s_bm_stress.payload
        res
    }
    
    val (s_to_chem, s_bypass_final) = StreamFork2(join_chem)
    
    val uChem = new ChemicalPotential()
    uChem.io.ricci_scalar << s_to_chem.translateWith(s_to_chem.payload.scalar)
    
    val join_stress = StreamJoin.arg(uChem.io.mu_out, s_bypass_final)
    
    val uStress = new ActiveStressUnit()
    uStress.io.input_valid := join_stress.valid
    join_stress.ready := uStress.io.output_ready
    
    uStress.io.mu_in    := join_stress.payload._1
    // Easy access now!
    uStress.io.g_metric := join_stress.payload._2.ctx.context.windata.center_g
    uStress.io.y_vec    := join_stress.payload._2.ctx.context.windata.center_y
    
    val join_final = StreamJoin.arg(uStress.io.output, s_bypass_final)
    
    val s_update = join_final.translateWith {
        val t_stress = join_final.payload._1
        
        val ctx_bundle = join_final.payload._2 // ChemContext
        val ricci_tensor_val = ctx_bundle.ctx.ricci
        val g_old = ctx_bundle.ctx.context.windata.center_g
        
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
