package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class WindowBuffer[T <: Data](dataType: HardType[T], config: AureliaConfig) extends Component {
  val io = new Bundle {
    val input  = slave Stream(dataType())
    val output = master Stream(Vec(Vec(Vec(dataType(), 5), 5), 5))
    
    val grid_dim = in(Vec(UInt(32 bits), 3)) 
  }

  val NX = io.grid_dim(0)
  val NY = io.grid_dim(1)
  val NZ = io.grid_dim(2)
  
  val MAX_PLANE = config.maxWidth * config.maxWidth 
  
  // Dynamic Memories based on max config
  val plane0 = Mem(dataType(), MAX_PLANE)
  val plane1 = Mem(dataType(), MAX_PLANE)
  val plane2 = Mem(dataType(), MAX_PLANE)
  val plane3 = Mem(dataType(), MAX_PLANE)

  val cx = Counter(32 bits)
  val cy = Counter(32 bits)
  val cz = Counter(32 bits)
  
  // Event-based counter for startup latency
  // 5x5x5 window needs 2 planes + 2 lines + 2 pixels latency?
  // Actually, we need to center the window.
  // Center is at z=0, y=0, x=0 relative to window?
  // Center (2,2,2).
  // We need to wait until we have written enough pixels to populate the window centered at (0,0,0) (conceptually shifted).
  // Total Latency = 2 * (NX * NY) + 2 * NX + 2 + Pipeline?
  // Let's rely on an explicit "primed" flag derived from counters.
  
  val total_pixels_processed = Reg(UInt(32 bits)) init(0)
  
  val write_enable = io.input.fire
  
  when(write_enable) {
     total_pixels_processed := total_pixels_processed + 1
     cx.increment()
     when(cx.value === NX - 1) {
       cx.clear()
       cy.increment()
       when(cy.value === NY - 1) {
         cy.clear()
         cz.increment()
         when(cz.value === NZ - 1) {
           cz.clear()
           total_pixels_processed := 0 // Reset for next frame?
         }
       }
     }
  }
  
  // Calculate dynamic latency threshold
  // We produce valid window when the 'center' pixel is being processed?
  // No, a line buffer typically outputs the data corresponding to the center pixel *when* the bottom-right pixel of the window acts as input?
  // Window (2,2,2) is center.
  // We need data at z+2. The input IS z+2 relative to z.
  // So when we input pixel P, we can form a window centered at P - 2 planes - 2 lines - 2 pixels.
  // So valid output appears when count >= (2*Area + 2*Width + 2).
  
  val latency_threshold = (NX * NY * 2) + (NX * 2) + 2
  
  val is_primed = RegInit(False)
  when(write_enable && (total_pixels_processed === latency_threshold)) {
      is_primed := True
  }
  
  // Note: Resetting is_primed at frame end?
  // For simplicity, we assume continuous stream or we reset.
  // If frame ends, we might need to flush.
  // The user review demanded "Robust AXI". 
  // Robust AXI means we shouldn't drop valid beats.
  // We output if input is valid AND we are primed.
  
  def makeDelay(input: T, max_depth: Int, current_depth: UInt, en: Bool): T = {
     val mem = Mem(dataType(), max_depth)
     val ptr = Counter(32 bits)
     when(en) {
         when(ptr.value === current_depth - 1) {
             ptr.clear()
         } otherwise {
             ptr.increment()
         }
     }
     
     // Robust Read/Write: Read at ptr, Write at ptr.
     // This acts as delay of 'current_depth'.
     val old_data = mem.readSync(ptr.value, enable = en)
     mem.write(ptr.value, input, enable = en)
     old_data
  }
  
  val plane_size = (NX * NY).resize(32 bits) 
  
  val z_plus2  = io.input.payload
  val z_plus1  = makeDelay(z_plus2,  MAX_PLANE, plane_size, io.input.fire)
  val z_center = makeDelay(z_plus1,  MAX_PLANE, plane_size, io.input.fire)
  val z_minus1 = makeDelay(z_center, MAX_PLANE, plane_size, io.input.fire)
  val z_minus2 = makeDelay(z_minus1, MAX_PLANE, plane_size, io.input.fire)
  
  // 5x5 generator also adds latency?
  // makedelay is RAM based. Window5x5Generator uses registers.
  
  case class Window5x5Generator(config_width: Int, source: T, valid: Bool, line_width: UInt) extends Area {
      val row_plus2  = source
      val row_plus1  = makeDelay(row_plus2,  config_width, line_width, valid)
      val row_center = makeDelay(row_plus1,  config_width, line_width, valid)
      val row_minus1 = makeDelay(row_center, config_width, line_width, valid)
      val row_minus2 = makeDelay(row_minus1, config_width, line_width, valid)
      
      val taps = Vec(api_makeWindowRow(row_plus2,  valid),
                     api_makeWindowRow(row_plus1,  valid),
                     api_makeWindowRow(row_center, valid),
                     api_makeWindowRow(row_minus1, valid),
                     api_makeWindowRow(row_minus2, valid))
      
      def result = Vec(taps(0), taps(1), taps(2), taps(3), taps(4))
  }
  
  // 2 registers delay to align center column in row? 
  // api_makeWindowRow returns 5 taps. Center is index 2.
  // Tap 0 (r0) is newest?
  // Row: [New, old, older, center, oldest...]
  // We need to clarify specific alignment, but "valid" logic handles the timing.
  // The pipeline delay inside row generation is 5 regs (shift). Center is at 2 delays.
  // So +2 cycles latency on top of planes.
  
  def api_makeWindowRow(in_px: T, en: Bool): Vec[T] = {
      val r0 = RegNextWhen(in_px, en)
      val r1 = RegNextWhen(r0, en)
      val r2 = RegNextWhen(r1, en)
      val r3 = RegNextWhen(r2, en)
      val r4 = RegNextWhen(r3, en)
      Vec(r0, r1, r2, r3, r4) // r0 is newest (x+2), r2 is center (x), r4 is x-2
  }
  
  // Apply generator
  val w_z_plus2  = Window5x5Generator(config.maxWidth, z_plus2,  io.input.fire, NX)
  val w_z_plus1  = Window5x5Generator(config.maxWidth, z_plus1,  io.input.fire, NX)
  val w_z_center = Window5x5Generator(config.maxWidth, z_center, io.input.fire, NX)
  val w_z_minus1 = Window5x5Generator(config.maxWidth, z_minus1, io.input.fire, NX)
  val w_z_minus2 = Window5x5Generator(config.maxWidth, z_minus2, io.input.fire, NX)
  
  val raw_win = Vec(w_z_plus2.result, w_z_plus1.result, w_z_center.result, w_z_minus1.result, w_z_minus2.result)
  
  val safe_win = Vec(Vec(Vec(dataType(), 5), 5), 5)
  for(z <- 0 until 5; y <- 0 until 5; x <- 0 until 5) {
      safe_win(z)(y)(x) := raw_win(z)(y)(x)
  }
  
  io.output.payload := safe_win
  
  // ROBUST VALIDITY:
  // Valid if input is valid AND we have processed enough pixels to form a window.
  // This signal halts if input halts, perfectly effectively.
  io.output.valid   := io.input.valid && is_primed
  
  // Backpressure: 
  // Since we construct a pipeline that consumes 1 pixel and produces 1 window (after latency),
  // we can accept input if downstream accepts output (or if we are filling up/not primed yet).
  // If not primed, we discard output (valid=0), so downstream ready is irrelevant?
  // No, if downstream is not ready, we cannot Advance Step.
  // So typically: ready = (output.ready || !is_primed)
  // BUT: if we advance when !is_primed, we are "consuming" pixels into the buffer.
  // We should generally always accept pixels to fill the buffer.
  // BUT what if we fill the buffer and then downstream says NOT READY?
  // Then we must stall the shift register.
  // So pipeline enable = (output.ready || !is_primed) && input.valid
  // This logic is implicitly handled by `write_enable` if we bind `io.input.ready` correctly.
  
  io.input.ready := io.output.ready || !is_primed 
}
class WindowBuffer3D extends WindowBuffer(Matrix3(), AureliaConfig()) // Default fallback wrapperbility if needed, but we will update usages.
