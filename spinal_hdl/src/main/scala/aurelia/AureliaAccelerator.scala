package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class AxionFieldCore extends Component {
    val io = new Bundle {
        val n_center    = in(Vector3())
        val curl_n      = in(Vector3())
        val torsion_mag = in(RealType())
        
        val input_valid = in Bool()
        val output_ready = in Bool()
        val axion_val   = master Stream(RealType())
    }
    
    val COUPLING = RealConst(1.2e-4) 
    val HALF_VAL = RealConst(0.5)
    
    case class AxionInput() extends Bundle {
       val n = Vector3()
       val c = Vector3()
       val t = RealType()
    }
    val in_stream = Stream(AxionInput())
    in_stream.valid := io.input_valid
    in_stream.payload.n := io.n_center
    in_stream.payload.c := io.curl_n
    in_stream.payload.t := io.torsion_mag
    
    
    val s1 = in_stream.translateWith {
        val i = in_stream.payload
        val nx = i.n(0); val ny = i.n(1); val nz = i.n(2)
        val cx = i.c(0); val cy = i.c(1); val cz = i.c(2)
        
        val dot_prod = (nx * cx) + (ny * cy) + (nz * cz)
        (dot_prod, i.t)
    }.s2mPipe()
    
    val s2 = s1.translateWith {
        val (dot, torsion) = s1.payload
        val helicity = dot.abs
        val term = helicity + (HALF_VAL * torsion)
        term * COUPLING
    }.s2mPipe()

    io.axion_val << s2
}

import spinal.lib.bus.amba4.axis._

// The internal Core Logic (renamed from AureliaAccelerator)
class AureliaCore extends Component {
  val io = new Bundle {
    val clk = in Bool() 
    val rst = in Bool()
    
    val host_data_in  = slave Stream(SymmetricMatrix3())
    val host_data_out = master Stream(SymmetricMatrix3())
    
    val y_vec_in      = slave Stream(Vector3())

    val voxel_count_out = out(UInt(32 bits))
    
    val axion_valid = out Bool()
    val axion_val   = out(RealType())
    
    val grid_config = in(Vec(UInt(32 bits), 3))
  }

  val area = new ClockingArea(ClockDomain(io.clk, io.rst)) {
      val uEngine = new RicciFlowEngine()
      
      uEngine.io.mem_data_in << io.host_data_in
      io.host_data_out       << uEngine.io.mem_data_out
      uEngine.io.y_vec_in    << io.y_vec_in
      uEngine.io.grid_dim    := io.grid_config
      
      io.voxel_count_out := uEngine.io.voxel_count_out
      
      val uAxion = new AxionFieldCore()
      uAxion.io.n_center.foreach(_ := RealConst(0.0)) // Fixed .map to .foreach for Unit return
      uAxion.io.curl_n.foreach(_ := RealConst(0.0))
      uAxion.io.torsion_mag := RealConst(1.0)
      uAxion.io.input_valid := True
      uAxion.io.output_ready := True 
      
      io.axion_valid := uAxion.io.axion_val.valid
      io.axion_val   := uAxion.io.axion_val.payload
  }
}

// SOTA Top Level with AXI4-Stream Interface
class AureliaAccelerator extends Component {
  val axiConfig = Axi4StreamConfig(
     dataWidth = 512, // High Bandwidth (HBM/DDR4 friendly)
     useLast   = true,
     useKeep   = true,
     useId     = false,
     useDest   = false
  )

  val io = new Bundle {
     // AXI4-Stream Interfaces (DMA Connected)
     val s_axis_data = slave(Axi4Stream(axiConfig))
     val m_axis_data = master(Axi4Stream(axiConfig))
     
     // Control / Config (AXI-Lite usually, but keeping simple for now)
     val grid_dim    = in(Vec(UInt(32 bits), 3))
     val status_axion = out(RealType())
  }
  
  // Adapter from 512-bit AXI to Internal Objects
  // SymmetricMatrix3 = 6 * 32 bit = 192 bits.
  // We can pack 2 per 512-bit word (384 bits used) or 1.
  // For simplicity, we assume 1 per beat, aligned to LSB.
  // TODO: Implement full width packing for max bandwidth.
  
  val core = new AureliaCore()
  core.io.clk := clockDomain.readClockWire
  core.io.rst := clockDomain.readResetWire
  core.io.grid_config := io.grid_dim
  
  // Map Status
  io.status_axion := core.io.axion_val
  
  // =========================================================
  // Stream Unpacking (AXI -> Core)
  // =========================================================
  val unpackedStream = Stream(SymmetricMatrix3())
  
  unpackedStream.valid := io.s_axis_data.valid
  io.s_axis_data.ready := unpackedStream.ready
  
  // Deserialization Logic
  val data_in_bits = io.s_axis_data.payload.data
  // Mapping first 192 bits to symmetric matrix
  val vec_sym = Vec(SFix(16 exp, -16 exp), 6) // Temporary fixed mapping
  for(i <- 0 until 6) {
     vec_sym(i).raw := data_in_bits(i*32 + 31 downto i*32).asSInt
  }
  unpackedStream.payload := unpackSymmetric(vec_sym)
  
  core.io.host_data_in << unpackedStream
  
  // Stub unused inputs for now (Y-Vector)
  core.io.y_vec_in.valid := False
  core.io.y_vec_in.payload.map(_ := RealConst(0.0))

  // =========================================================
  // Stream Packing (Core -> AXI)
  // =========================================================
  val core_out = core.io.host_data_out
  
  io.m_axis_data.valid := core_out.valid
  core_out.ready       := io.m_axis_data.ready
  
  // Serialization Logic
  val packed_sym = packSymmetric(core_out.payload)
  val out_bits = Bits(512 bits)
  out_bits := 0
  
  for(i <- 0 until 6) {
     out_bits(i*32 + 31 downto i*32) := packed_sym(i).raw.asBits
  }
  
  io.m_axis_data.payload.data := out_bits
  
  // =========================================================
  // TLAST Generation (Critical for DMA)
  // =========================================================
  val tlast_counter = Reg(UInt(32 bits)) init(0)
  val total_voxels  = io.grid_dim(0) * io.grid_dim(1) * io.grid_dim(2)
  
  // Note: io.grid_dim is dynamic. If it changes mid-stream, behavior is undefined.
  // We assume configuration is static during a run.
  
  val is_last_beat = tlast_counter === (total_voxels - 1)
  
  when(io.m_axis_data.fire) {
     when(is_last_beat) {
        tlast_counter := 0
     } otherwise {
        tlast_counter := tlast_counter + 1
     }
  }
  
  io.m_axis_data.payload.last := is_last_beat
  io.m_axis_data.payload.keep := B(0xFFFFFFFFFFFFFFFFL, 64 bits) // Keep all
}

object AureliaAcceleratorVerilog {
  def main(args: Array[String]): Unit = {
    SpinalVerilog(new AureliaAccelerator())
  }
}
