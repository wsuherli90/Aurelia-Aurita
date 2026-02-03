package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class ActiveStressUnit extends Component {
  val io = new Bundle {
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    val mu_in     = in(RealType())
    val g_metric  = in(Matrix3())
    val y_vec     = in(Vector3()) 
    
    val output = master Stream(Matrix3())
  }
  
  case class StressInput() extends Bundle {
     val mu = RealType()
     val g  = Matrix3()
     val y  = Vector3()
  }
  
  val input_stream = Stream(StressInput())
  input_stream.valid := io.input_valid
  input_stream.payload.mu := io.mu_in
  input_stream.payload.g  := io.g_metric
  input_stream.payload.y  := io.y_vec
  
  val uWLC = new WLCStressCore()
  uWLC.io.input_valid := input_stream.valid
  uWLC.io.mu_in := input_stream.payload.mu
  uWLC.io.g_metric := input_stream.payload.g
  uWLC.io.y_vec := input_stream.payload.y
  
  uWLC.io.output_ready := io.output.ready 
  
  io.output << uWLC.io.output
}
