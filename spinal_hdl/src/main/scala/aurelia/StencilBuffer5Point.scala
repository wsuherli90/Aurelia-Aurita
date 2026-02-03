package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class StencilBuffer5Point extends Component {
  val io = new Bundle {
    val input  = slave Stream(Matrix3())
    
    val x_out_idx = out(UInt(32 bits))
    val y_out_idx = out(UInt(32 bits))
    val z_out_idx = out(UInt(32 bits))
    
    val win_x  = out(StencilWindow())
    val win_y  = out(StencilWindow())
    val win_z  = out(StencilWindow())
    
    val output = master Stream(Bool()) 
  }
  
  case class StencilPayload() extends Bundle {
     val win_x = StencilWindow()
     val win_y = StencilWindow()
     val win_z = StencilWindow()
  }
  
  val ROW_SIZE   = AureliaTypes.NX_C
  val PLANE_SIZE = AureliaTypes.NX_C * AureliaTypes.NY_C
  
  val startup_delay = (2 * PLANE_SIZE) + (2 * ROW_SIZE) + 2
  val startup_counter = Counter(Int.MaxValue)
  val is_startup = startup_counter.value < startup_delay
  
  val enable = Bool()
  
  when(is_startup) {
      enable := io.input.valid
      io.input.ready := True 
  } otherwise {
      enable := io.input.valid && io.output.ready
      io.input.ready := io.output.ready 
  }
  
  when(enable && is_startup) {
      startup_counter.increment()
  }

  def makeDelayLine(data_in: Matrix3, en: Bool, depth: Int): Matrix3 = {
     val mem = Mem(Matrix3(), depth)
     val ptr = Counter(depth, en)
     
     mem.write(ptr.value, data_in, enable = en)
     val readVal = mem.readSync(ptr.value, enable = en) 
     readVal
  }
  
  val z_taps = Vec(Matrix3(), 5)
  z_taps(0) := io.input.payload
  
  z_taps(1) := makeDelayLine(z_taps(0), enable, PLANE_SIZE)
  z_taps(2) := makeDelayLine(z_taps(1), enable, PLANE_SIZE)
  z_taps(3) := makeDelayLine(z_taps(2), enable, PLANE_SIZE)
  z_taps(4) := makeDelayLine(z_taps(3), enable, PLANE_SIZE)

  val y_taps = Vec(Matrix3(), 5)
  y_taps(0) := z_taps(2)
  
  y_taps(1) := makeDelayLine(y_taps(0), enable, ROW_SIZE)
  y_taps(2) := makeDelayLine(y_taps(1), enable, ROW_SIZE)
  y_taps(3) := makeDelayLine(y_taps(2), enable, ROW_SIZE)
  y_taps(4) := makeDelayLine(y_taps(3), enable, ROW_SIZE)
  
  val x_taps = Vec(Reg(Matrix3()), 5)
  
  when(enable) {
    x_taps(0) := y_taps(2)
    x_taps(1) := x_taps(0)
    x_taps(2) := x_taps(1) 
    x_taps(3) := x_taps(2)
    x_taps(4) := x_taps(3)
  }
  
  val SYNC_DELAY = (2 * ROW_SIZE) + 2
  
  val z_aligned = Vec(Matrix3(), 5)
  z_aligned(0) := makeDelayLine(z_taps(0), enable, SYNC_DELAY)
  z_aligned(1) := makeDelayLine(z_taps(1), enable, SYNC_DELAY)
  z_aligned(2) := x_taps(2) 
  z_aligned(3) := makeDelayLine(z_taps(3), enable, SYNC_DELAY)
  z_aligned(4) := makeDelayLine(z_taps(4), enable, SYNC_DELAY)

  val y_aligned = Vec(Matrix3(), 5)
  for(i <- 0 until 5) {
     val r1 = RegNextWhen(y_taps(i), enable)
     val r2 = RegNextWhen(r1, enable)
     y_aligned(i) := r2
  }
  
  val final_x = StencilWindow()
  val final_y = StencilWindow()
  val final_z = StencilWindow()
  
  final_x(0) := x_taps(4) 
  final_x(1) := x_taps(3)
  final_x(2) := x_taps(2)
  final_x(3) := x_taps(1)
  final_x(4) := x_taps(0) 
  
  final_y(0) := y_aligned(4)
  final_y(1) := y_aligned(3)
  final_y(2) := y_aligned(2)
  final_y(3) := y_aligned(1)
  final_y(4) := y_aligned(0)
  
  final_z(0) := z_aligned(4)
  final_z(1) := z_aligned(3)
  final_z(2) := z_aligned(2)
  final_z(3) := z_aligned(1)
  final_z(4) := z_aligned(0)
  
  io.win_x := final_x
  io.win_y := final_y
  io.win_z := final_z
  
  io.output.valid := !is_startup && io.input.valid
  io.output.payload := True 
  
  val c_x = Counter(AureliaTypes.NX_C)
  val c_y = Counter(AureliaTypes.NY_C)
  val c_z = Counter(AureliaTypes.NZ_C)
  
  when(enable && !is_startup) {
      c_x.increment()
      when(c_x.willOverflow) {
          c_y.increment()
          when(c_y.willOverflow) {
              c_z.increment()
          }
      }
  }
  
  io.x_out_idx := c_x.value
  io.y_out_idx := c_y.value
  io.z_out_idx := c_z.value
}
