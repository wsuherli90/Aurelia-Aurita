package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class WindowBuffer3D extends Component {
  val io = new Bundle {
    val input  = slave Stream(Matrix3())
    val output = master Stream(Window3D())
    
    val grid_dim = in(Vec(UInt(32 bits), 3)) 
  }

  val NX = io.grid_dim(0)
  val NY = io.grid_dim(1)
  val NZ = io.grid_dim(2)
  
  val MAX_NX = 256
  val MAX_NY = 256
  val MAX_PLANE = MAX_NX * MAX_NY
  
  val plane0 = Mem(Matrix3(), MAX_PLANE)
  val plane1 = Mem(Matrix3(), MAX_PLANE)
  
  val cx = Counter(32 bits)
  val cy = Counter(32 bits)
  val cz = Counter(32 bits)
  
  val write_enable = io.input.fire
  val write_addr   = lin_addr.value
  val lin_addr = Counter(32 bits)
  
  when(write_enable) {
     cx.increment()
     lin_addr.increment()
     when(cx.value === NX - 1) {
       cx.clear()
       cy.increment()
       when(cy.value === NY - 1) {
         cy.clear()
         cz.increment()
         lin_addr.clear()
         when(cz.value === NZ - 1) {
           cz.clear()
         }
       }
     }
  }
  
  def makeDelay(input: Matrix3, max_depth: Int, current_depth: UInt, en: Bool): Matrix3 = {
     val mem = Mem(Matrix3(), max_depth)
     val ptr = Counter(32 bits)
     when(en) {
         when(ptr.value === current_depth - 1) {
             ptr.clear()
         } otherwise {
             ptr.increment()
         }
     }
     
     val old_data = mem.readSync(ptr.value, enable = en)
     mem.write(ptr.value, input, enable = en)
     
     old_data
  }
  
  val plane_size = (NX * NY).resize(32 bits) 
  
  val z_plus1  = io.input.payload
  val z_center = makeDelay(z_plus1, MAX_PLANE, plane_size, io.input.fire)
  val z_minus1 = makeDelay(z_center, MAX_PLANE, plane_size, io.input.fire)
  
  case class Window3x3Generator(source: Matrix3, valid: Bool, line_width: UInt) extends Area {
      val row_plus1 = source
      val row_center = makeDelay(row_plus1, MAX_NX, line_width, valid)
      val row_minus1 = makeDelay(row_center, MAX_NX, line_width, valid)
      
      val taps = Vec(api_makeWindowRow(row_plus1, valid),
                     api_makeWindowRow(row_center, valid),
                     api_makeWindowRow(row_minus1, valid))
      
      def result = Vec(taps(0), taps(1), taps(2))
  }
  
  def api_makeWindowRow(in_px: Matrix3, en: Bool): Vec[Matrix3] = {
      val r0 = RegNextWhen(in_px, en)
      val r1 = RegNextWhen(r0, en)
      val r2 = RegNextWhen(r1, en) 
      Vec(r0, r1, r2)
  }
  
  val w_z_plus1 = Window3x3Generator(z_plus1, io.input.fire, NX)
  val w_z_center = Window3x3Generator(z_center, io.input.fire, NX)
  val w_z_minus1 = Window3x3Generator(z_minus1, io.input.fire, NX)
  
  val PIPELINE_DELAY = 4 
  val cx_d = Delay(cx.value, PIPELINE_DELAY)
  val cy_d = Delay(cy.value, PIPELINE_DELAY)
  val cz_d = Delay(cz.value, PIPELINE_DELAY)
  
  val raw_win = Vec(w_z_plus1.result, w_z_center.result, w_z_minus1.result)
  
  val safe_win = Window3D()
  
  for(z <- 0 until 3; y <- 0 until 3; x <- 0 until 3) {
      val is_x_low  = (x == 2) && (cx_d === 0)
      val is_x_high = (x == 0) && (cx_d === NX - 1)
      val is_y_low  = (y == 2) && (cy_d === 0)
      val is_y_high = (y == 0) && (cy_d === NY - 1)
      val is_z_low  = (z == 2) && (cz_d === 0)
      val is_z_high = (z == 0) && (cz_d === NZ - 1)
      
      val clamp_x = is_x_low || is_x_high
      val clamp_y = is_y_low || is_y_high
      val clamp_z = is_z_low || is_z_high
      
      val px = raw_win(z)(y)(x)
      
      val center_px = raw_win(1)(1)(1)
      
      when(clamp_x || clamp_y || clamp_z) {
          safe_win(z)(y)(x) := center_px 
      } otherwise {
          safe_win(z)(y)(x) := px
      }
  }
  
  io.output.payload := safe_win
  io.output.valid   := RegNext(io.input.valid, PIPELINE_DELAY) 
  
  io.input.ready := io.output.ready 
}
