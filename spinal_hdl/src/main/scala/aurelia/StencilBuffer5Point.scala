package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

/**
 * Stencil Buffer (5-Point)
 * Buffered Window Generator for 3D Grid.
 * Input: Stream[Matrix3]
 * Output: Stream[StencilWindow for X, Y, Z]
 * 
 * Flow Control: 
 * - Output is valid only after buffer fills.
 * - Backpressure propagates to input (stops filling buffer).
 */
class StencilBuffer5Point extends Component {
  val io = new Bundle {
    val input  = slave Stream(Matrix3())
    
    // Output Bundle
    val x_out_idx = out(UInt(32 bits))
    val y_out_idx = out(UInt(32 bits))
    val z_out_idx = out(UInt(32 bits))
    
    val win_x  = out(StencilWindow())
    val win_y  = out(StencilWindow())
    val win_z  = out(StencilWindow())
    
    val output = master Stream(Bool()) // Carrying dummy payload, actual data in signals above
    // Or better: Bundle the windows into the Stream Payload
  }
  
  // Define Output Payload Bundle
  case class StencilPayload() extends Bundle {
     val win_x = StencilWindow()
     val win_y = StencilWindow()
     val win_z = StencilWindow()
  }
  // We will drive io.output with this, and also wire individual ports for backward compat/debug
  
  val ROW_SIZE   = AureliaTypes.NX_C
  val PLANE_SIZE = AureliaTypes.NX_C * AureliaTypes.NY_C
  
  // Flow Control Handling
  // We process only when Input is Valid AND Downstream is Ready.
  // Exception: During startup (filling), we don't produce valid output, but we must accept input.
  // So 'ready' to input depends on buffer state?
  // Strategy:
  // - While Filling: Input.Ready = True (Always accept). Output.Valid = False.
  // - Filled: Input.Ready = Output.Ready (Pipeline). Output.Valid = Input.Valid.
  
  val startup_delay = (2 * PLANE_SIZE) + (2 * ROW_SIZE) + 2
  val startup_counter = Counter(Int.MaxValue)
  val is_startup = startup_counter.value < startup_delay
  
  // Enable signal for internal state updates
  // If Startup: enable = input.valid
  // If Running: enable = input.valid && output.ready
  val enable = Bool()
  
  when(is_startup) {
      enable := io.input.valid
      io.input.ready := True // Always consume during fill
  } otherwise {
      enable := io.input.valid && io.output.ready
      io.input.ready := io.output.ready // Backpressure
  }
  
  when(enable && is_startup) {
      startup_counter.increment()
  }

  // =============================================================
  // Memory Delay Lines
  // =============================================================
  
  def makeDelayLine(data_in: Matrix3, en: Bool, depth: Int): Matrix3 = {
     val mem = Mem(Matrix3(), depth)
     val ptr = Counter(depth, en)
     
     mem.write(ptr.value, data_in, enable = en)
     val readVal = mem.readSync(ptr.value, enable = en) // Read old value (N cycles ago)
     // Read-First or Write-First?
     // Spinal Mem readSync with enable: 
     // If enable is high, read occurs. Address is ptr.
     // Write also occurs at ptr.
     // Standard BlockRAM in Write-First mode returns New Data? Read-First returns Old?
     // Spinal default is often Read-First (Old data). We assume Old Data (Delay behavior).
     readVal
  }
  
  // Z Taps
  val z_taps = Vec(Matrix3(), 5)
  z_taps(0) := io.input.payload
  
  z_taps(1) := makeDelayLine(z_taps(0), enable, PLANE_SIZE)
  z_taps(2) := makeDelayLine(z_taps(1), enable, PLANE_SIZE)
  z_taps(3) := makeDelayLine(z_taps(2), enable, PLANE_SIZE)
  z_taps(4) := makeDelayLine(z_taps(3), enable, PLANE_SIZE)

  // Y Taps (Feed from Z2 - Center Plane)
  val y_taps = Vec(Matrix3(), 5)
  y_taps(0) := z_taps(2)
  
  y_taps(1) := makeDelayLine(y_taps(0), enable, ROW_SIZE)
  y_taps(2) := makeDelayLine(y_taps(1), enable, ROW_SIZE)
  y_taps(3) := makeDelayLine(y_taps(2), enable, ROW_SIZE)
  y_taps(4) := makeDelayLine(y_taps(3), enable, ROW_SIZE)
  
  // X Taps (Feed from Y2 - Center Row)
  val x_taps = Vec(Reg(Matrix3()), 5)
  
  when(enable) {
    x_taps(0) := y_taps(2)
    x_taps(1) := x_taps(0)
    x_taps(2) := x_taps(1) // Center Pixel
    x_taps(3) := x_taps(2)
    x_taps(4) := x_taps(3)
  }
  
  // =============================================================
  // Z-Neighbor Alignment
  // =============================================================
  val SYNC_DELAY = (2 * ROW_SIZE) + 2
  
  val z_aligned = Vec(Matrix3(), 5)
  // Use 'enable' to gate these delays too
  z_aligned(0) := makeDelayLine(z_taps(0), enable, SYNC_DELAY)
  z_aligned(1) := makeDelayLine(z_taps(1), enable, SYNC_DELAY)
  z_aligned(2) := x_taps(2) // Center
  z_aligned(3) := makeDelayLine(z_taps(3), enable, SYNC_DELAY)
  z_aligned(4) := makeDelayLine(z_taps(4), enable, SYNC_DELAY)

  // Y-Neighbor Alignment (2 cycles delay to match X shift)
  // We use regs as simple delay
  val y_aligned = Vec(Matrix3(), 5)
  for(i <- 0 until 5) {
     // Delay(val, 2)
     val r1 = RegNextWhen(y_taps(i), enable)
     val r2 = RegNextWhen(r1, enable)
     y_aligned(i) := r2
  }
  
  // =============================================================
  // Outputs
  // =============================================================
  val final_x = StencilWindow()
  val final_y = StencilWindow()
  val final_z = StencilWindow()
  
  // Reverse taps for standard order -2 to +2?
  // x_taps(0) is oldest? No, x_taps(0) is Newest input.
  // Standard Stencil: -2, -1, 0, 1, 2.
  // If flow is right-to-left (Newest at 0), then:
  // index 0 (-2) should be Oldest (taps(4)).
  // index 4 (+2) should be Newest (taps(0)).
  
  final_x(0) := x_taps(4) // Oldest (X-2)
  final_x(1) := x_taps(3)
  final_x(2) := x_taps(2)
  final_x(3) := x_taps(1)
  final_x(4) := x_taps(0) // Newest (X+2)
  
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
  io.output.payload := True // Dummy
  
  // Counters
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
